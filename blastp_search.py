#!/usr/bin/env python3
"""
BLASTP search program with per-taxonomy search strategy.

Takes a list of taxonomy names and a FASTA file of query protein sequences,
runs one BLASTP search per (query, taxonomy) pair against NCBI, filters by
E-value, query coverage, and percent identity, and writes a random selection
of matching sequences to per-query FASTA output files.

Running one BLAST per taxonomy ensures that each taxon gets its own hit-list
pool, preventing heavily-annotated species from crowding out rarer ones.

Requirements:
    pip install biopython tqdm
"""

import argparse
import csv
import random
import sys
import time
from pathlib import Path

VERSION = "1.1.0"

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from tqdm import tqdm


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=f"BLASTP search with per-taxonomy filtering (v{VERSION})",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-q", "--query", required=True, metavar="FASTA",
        help="Input FASTA file containing query protein sequences",
    )
    p.add_argument(
        "-t", "--taxonomies", nargs="+", metavar="NAME",
        help="One or more taxonomy names, e.g. 'Bacillus subtilis'",
    )
    p.add_argument(
        "-T", "--taxonomy-file", metavar="FILE",
        help="File with one taxonomy name per line (can be combined with -t)",
    )
    p.add_argument(
        "--email", required=True,
        help="Your e-mail address (required by NCBI)",
    )
    p.add_argument(
        "-e", "--evalue", type=float, default=1e-5,
        help="Maximum E-value cutoff",
    )
    p.add_argument(
        "-c", "--coverage", type=float, default=50.0,
        help="Minimum query coverage in percent (computed from merged multi-HSP intervals)",
    )
    p.add_argument(
        "-i", "--identity", type=float, default=0.0,
        help="Minimum percent identity cutoff (0 = no filter)",
    )
    p.add_argument(
        "-n", "--max-seqs", type=int, default=10,
        help="Maximum number of sequences to return per taxonomy per query (random draw)",
    )
    p.add_argument(
        "-o", "--outdir", default="blastp_results",
        help="Directory where output FASTA files will be written",
    )
    p.add_argument(
        "--db", default="nr",
        help="NCBI BLAST database",
    )
    p.add_argument(
        "--hitlist-size", type=int, default=500,
        help="Number of BLAST hits to fetch per taxonomy (pool for random draw)",
    )
    p.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducible sequence selection",
    )
    p.add_argument(
        "--max-retries", type=int, default=3,
        help="Maximum retry attempts for NCBI calls on transient failures",
    )
    p.add_argument(
        "--retry-delay", type=float, default=5.0,
        help="Base delay in seconds between retries (doubles each attempt)",
    )
    p.add_argument(
        "--resume", action="store_true", default=False,
        help="Resume a previous run, skipping (query, taxonomy) pairs in checkpoint.txt",
    )
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(0)
    args = p.parse_args()

    # Merge taxonomies from -t and -T, deduplicate, preserve order
    # taxonomy_max_seqs holds per-taxonomy -n overrides from the -T file
    taxonomy_max_seqs: dict[str, int] = {}
    taxonomies = list(args.taxonomies) if args.taxonomies else []
    if args.taxonomy_file:
        with open(args.taxonomy_file) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                name = parts[0].strip()
                if len(parts) >= 2:
                    try:
                        taxonomy_max_seqs[name] = int(parts[1].strip())
                    except ValueError:
                        pass
                taxonomies.append(name)
    args.taxonomy_max_seqs = taxonomy_max_seqs
    seen: set[str] = set()
    args.taxonomies = [t for t in taxonomies if not (t in seen or seen.add(t))]
    if not args.taxonomies:
        p.error("Provide at least one taxonomy via -t/--taxonomies or -T/--taxonomy-file")
    return args


# ---------------------------------------------------------------------------
# Retry
# ---------------------------------------------------------------------------

def retry_call(fn, *args, max_attempts: int = 3, base_delay: float = 5.0,
               label: str = "NCBI call", **kwargs):
    """Call fn(*args, **kwargs) with exponential-backoff retries on failure."""
    for attempt in range(max_attempts):
        try:
            return fn(*args, **kwargs)
        except Exception as exc:
            if attempt == max_attempts - 1:
                raise
            delay = base_delay * (2 ** attempt)
            tqdm.write(
                f"  [warning] {label} failed (attempt {attempt + 1}/{max_attempts}): "
                f"{exc}. Retrying in {delay:.0f}s...",
                file=sys.stderr,
            )
            time.sleep(delay)


# ---------------------------------------------------------------------------
# Checkpoint
# ---------------------------------------------------------------------------

def load_checkpoint(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def save_checkpoint(path: Path, key: str) -> None:
    with open(path, "a") as fh:
        fh.write(key + "\n")


# ---------------------------------------------------------------------------
# Taxonomy resolution
# ---------------------------------------------------------------------------

def resolve_taxonomies(taxonomies: list[str], max_retries: int,
                       retry_delay: float) -> dict[str, dict]:
    """
    Resolve taxonomy names to NCBI TaxIDs and canonical scientific names.
    Returns {user_name: {"taxid": str|None, "scientific_name": str}}.
    """
    resolved: dict[str, dict] = {}
    print("Resolving taxonomy names via NCBI...")
    for name in tqdm(taxonomies, desc="Resolving", unit="taxon"):
        try:
            handle = retry_call(
                Entrez.esearch, db="taxonomy", term=name,
                max_attempts=max_retries, base_delay=retry_delay,
                label=f"esearch '{name}'",
            )
            record = Entrez.read(handle)
            handle.close()
            ids = record["IdList"]
        except Exception as exc:
            tqdm.write(f"  [warning] TaxID lookup failed for '{name}': {exc}", file=sys.stderr)
            resolved[name] = {"taxid": None, "scientific_name": name}
            continue

        if not ids:
            tqdm.write(f"  [warning] No TaxID found for '{name}', using as-is", file=sys.stderr)
            resolved[name] = {"taxid": None, "scientific_name": name}
            continue

        if len(ids) > 1:
            tqdm.write(
                f"  [warning] '{name}' matched {len(ids)} taxa, using first (TaxID {ids[0]})",
                file=sys.stderr,
            )

        taxid = ids[0]
        try:
            handle = retry_call(
                Entrez.efetch, db="taxonomy", id=taxid, retmode="xml",
                max_attempts=max_retries, base_delay=retry_delay,
                label=f"efetch taxid {taxid}",
            )
            tax_records = Entrez.read(handle)
            handle.close()
            sci_name = tax_records[0]["ScientificName"]
        except Exception as exc:
            tqdm.write(f"  [warning] Could not fetch canonical name for TaxID {taxid}: {exc}",
                       file=sys.stderr)
            sci_name = name

        resolved[name] = {"taxid": taxid, "scientific_name": sci_name}
        tqdm.write(f"  {name!r} -> TaxID {taxid} ({sci_name})")
        time.sleep(0.34)

    print()
    return resolved


def build_entrez_query(info: dict) -> str:
    """Build a single-taxonomy Entrez query for BLAST."""
    org_name = info["scientific_name"]
    return f'"{org_name}"[Organism]'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def trim_description(title: str, max_len: int = 120) -> str:
    """Keep only the first entry of a MULTISPECIES title and cap length."""
    first = title.split(" >")[0]
    return first if len(first) <= max_len else first[:max_len - 1] + "…"


def safe_filename(s: str) -> str:
    """Convert a sequence ID into a safe file-name stem."""
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


def merged_query_coverage(hsps, query_length: int) -> float:
    """
    Compute query coverage from merged non-overlapping HSP intervals.
    Handles multi-domain alignments where several HSPs each contribute
    coverage to different parts of the query.
    """
    intervals = sorted(
        (min(h.query_start, h.query_end), max(h.query_start, h.query_end))
        for h in hsps
    )
    cur_start, cur_end = intervals[0]
    merged_len = 0
    for start, end in intervals[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            merged_len += cur_end - cur_start + 1
            cur_start, cur_end = start, end
    merged_len += cur_end - cur_start + 1
    return merged_len / query_length * 100


def extract_accession(raw: str) -> str:
    """
    Normalise a BLAST accession field to a bare accession number.

    BLAST XML returns titles such as:
      'sp|P12345|PROT_BACSU Protein ... [Bacillus subtilis]'
    or plain accessions like 'XP_123456.1'.
    """
    parts = raw.split("|")
    if len(parts) >= 3:
        return parts[1]          # UniProt style: db|acc|entry
    if len(parts) == 2:
        return parts[-1]         # gb|acc  or similar
    return parts[0].split()[0]   # plain accession, drop description


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def run_blastp(seq: str, entrez_query: str, evalue: float, hitlist_size: int,
               db: str, max_retries: int, retry_delay: float):
    """Submit a BLASTP query to NCBI and return parsed BLAST records."""
    with tqdm(bar_format="  {desc} {elapsed}", desc="  Waiting for BLAST...", leave=False) as pbar:
        result_handle = retry_call(
            NCBIWWW.qblast,
            "blastp", db, seq,
            entrez_query=entrez_query,
            expect=evalue,
            hitlist_size=hitlist_size,
            alignments=hitlist_size,
            descriptions=hitlist_size,
            max_attempts=max_retries,
            base_delay=retry_delay,
            label="BLASTP",
        )
        pbar.set_description("  BLAST complete")
    records = list(NCBIXML.parse(result_handle))
    result_handle.close()
    return records


def filter_hits(blast_record, evalue_cutoff: float, coverage_cutoff: float,
                identity_cutoff: float, query_length: int) -> list[dict]:
    """
    Return one entry per unique accession that passes all filters.
    Coverage is computed from merged HSP intervals; E-value and identity
    are taken from the best (first) HSP.
    """
    passing: list[dict] = []
    seen: set[str] = set()

    for alignment in blast_record.alignments:
        best_hsp = alignment.hsps[0]

        if best_hsp.expect > evalue_cutoff:
            continue

        identity_pct = best_hsp.identities / best_hsp.align_length * 100
        if identity_pct < identity_cutoff:
            continue

        cov = merged_query_coverage(alignment.hsps, query_length)
        if cov < coverage_cutoff:
            continue

        acc = extract_accession(alignment.accession)
        if acc in seen:
            continue
        seen.add(acc)

        passing.append({
            "accession": acc,
            "title": alignment.title,
            "evalue": best_hsp.expect,
            "coverage": round(cov, 1),
            "identity": round(identity_pct, 1),
        })

    return passing


def fetch_sequences(accessions: list[str], max_retries: int, retry_delay: float,
                    batch_size: int = 200) -> dict:
    """
    Fetch FASTA sequences from NCBI Protein for a list of accessions.
    Returns {record.id: SeqRecord}.
    """
    sequences: dict = {}
    batches = range(0, len(accessions), batch_size)
    with tqdm(batches, desc="  Fetching sequences", unit="batch", leave=False) as pbar:
        for start in pbar:
            batch = accessions[start:start + batch_size]
            pbar.set_postfix(batch=f"{start // batch_size + 1}", fetched=len(sequences))
            try:
                handle = retry_call(
                    Entrez.efetch,
                    db="protein", id=",".join(batch), rettype="fasta", retmode="text",
                    max_attempts=max_retries, base_delay=retry_delay,
                    label=f"Entrez fetch batch {start // batch_size + 1}",
                )
                for rec in SeqIO.parse(handle, "fasta"):
                    sequences[rec.id] = rec
                handle.close()
            except Exception as exc:
                tqdm.write(
                    f"  [warning] Entrez fetch failed for batch "
                    f"{start // batch_size + 1}: {exc}",
                    file=sys.stderr,
                )
            if start + batch_size < len(accessions):
                time.sleep(0.4)
    return sequences


def match_sequence(acc: str, sequences: dict):
    """
    Look up a sequence by accession, tolerating partial matches
    (BLAST accessions and Entrez IDs occasionally differ slightly).
    """
    if acc in sequences:
        return sequences[acc]
    for key, rec in sequences.items():
        if key.startswith(acc) or acc.startswith(key) or acc in key:
            return rec
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    Entrez.email = args.email

    if args.seed is not None:
        random.seed(args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    checkpoint_path = outdir / "checkpoint.txt"
    completed = load_checkpoint(checkpoint_path) if args.resume else set()
    if completed:
        print(f"Resuming: {len(completed)} (query, taxonomy) pair(s) already done, skipping.")

    # Load query sequences
    query_records = list(SeqIO.parse(args.query, "fasta"))
    if not query_records:
        sys.exit(f"Error: no sequences found in {args.query}")

    # Resolve taxonomy names to TaxIDs and canonical names
    resolved = resolve_taxonomies(args.taxonomies, args.max_retries, args.retry_delay)

    tsv_path = outdir / "summary.tsv"
    tsv_columns = [
        "query_id", "query_length", "taxonomy", "accession", "description",
        "evalue", "identity_pct", "query_coverage_pct", "output_fasta",
    ]
    tsv_mode = "a" if (args.resume and tsv_path.exists()) else "w"

    total_pairs = len(query_records) * len(resolved)
    print(f"Queries       : {len(query_records)}")
    print(f"Taxonomies    : {len(resolved)}  ({', '.join(args.taxonomies)})")
    print(f"BLAST calls   : {total_pairs}  (one per query × taxonomy)")
    print(f"E-value       : {args.evalue}")
    print(f"Coverage      : {args.coverage}%")
    print(f"Identity      : {args.identity}%")
    print(f"Max seqs/tax  : {args.max_seqs}")
    print(f"BLAST pool    : {args.hitlist_size} hits per search")
    print(f"Output dir    : {outdir}")
    print(f"Summary TSV   : {tsv_path}")
    print()

    with open(tsv_path, tsv_mode, newline="") as tsv_fh:
        writer = csv.DictWriter(tsv_fh, fieldnames=tsv_columns, delimiter="\t")
        if tsv_mode == "w":
            writer.writeheader()

        query_bar = tqdm(query_records, desc="Queries", unit="query")
        for qrec in query_bar:
            qid = qrec.id
            qlen = len(qrec.seq)
            query_bar.set_description(f"Query: {qid}")

            # Skip query entirely if all its (query, taxonomy) pairs are done
            pending_taxes = [
                name for name in resolved
                if f"{qid}::{name}" not in completed
            ]
            if not pending_taxes:
                tqdm.write(f"\n=== Query: {qid} — all taxonomies done, skipping ===")
                continue

            tqdm.write(f"\n=== Query: {qid} ({qlen} aa) | {len(pending_taxes)} taxonomy search(es) ===")

            out_path = outdir / f"{safe_filename(qid)}.fasta"
            # Append to FASTA if it already has content from a previous (resumed) run
            fasta_has_content = out_path.exists() and out_path.stat().st_size > 0

            for tax_name in pending_taxes:
                tax_info = resolved[tax_name]
                checkpoint_key = f"{qid}::{tax_name}"
                entrez_q = build_entrez_query(tax_info)

                tqdm.write(f"  --- Taxonomy: {tax_name} ---")

                # 1. BLAST (scoped to this taxonomy only)
                try:
                    blast_records = run_blastp(
                        str(qrec.seq), entrez_q, args.evalue,
                        args.hitlist_size, args.db, args.max_retries, args.retry_delay,
                    )
                except Exception as exc:
                    tqdm.write(f"  [error] BLAST failed after retries: {exc}", file=sys.stderr)
                    continue

                if not blast_records or not blast_records[0].alignments:
                    tqdm.write("  No BLAST hits returned.")
                    save_checkpoint(checkpoint_path, checkpoint_key)
                    continue

                # 2. Filter
                passing = filter_hits(
                    blast_records[0], args.evalue, args.coverage, args.identity, qlen,
                )
                tqdm.write(f"  Hits passing filters : {len(passing)}")

                if not passing:
                    tqdm.write("  Nothing to write.")
                    save_checkpoint(checkpoint_path, checkpoint_key)
                    continue

                # 3. Random draw (up to max_seqs from this taxonomy's pool)
                max_seqs = args.taxonomy_max_seqs.get(tax_name, args.max_seqs)
                selected = random.sample(passing, min(max_seqs, len(passing)))
                tqdm.write(f"  Selected             : {len(selected)}")

                # 4. Fetch sequences
                accessions = [h["accession"] for h in selected]
                sequences = fetch_sequences(accessions, args.max_retries, args.retry_delay)

                # 5. Write FASTA and TSV rows
                fasta_mode = "a" if fasta_has_content else "w"
                written = 0
                with open(out_path, fasta_mode) as fh:
                    for hit in tqdm(selected, desc="  Writing", unit="seq", leave=False):
                        rec = match_sequence(hit["accession"], sequences)
                        if rec is None:
                            tqdm.write(
                                f"  [warning] Sequence not fetched: {hit['accession']}",
                                file=sys.stderr,
                            )
                            writer.writerow({
                                "query_id": qid,
                                "query_length": qlen,
                                "taxonomy": tax_name,
                                "accession": hit["accession"],
                                "description": trim_description(hit["title"]),
                                "evalue": f"{hit['evalue']:.2e}",
                                "identity_pct": hit["identity"],
                                "query_coverage_pct": hit["coverage"],
                                "output_fasta": "NOT_FETCHED",
                            })
                            continue
                        rec.description = (
                            f"evalue={hit['evalue']:.2e} "
                            f"id={hit['identity']}% "
                            f"qcov={hit['coverage']}% | {rec.description}"
                        )
                        SeqIO.write(rec, fh, "fasta")
                        writer.writerow({
                            "query_id": qid,
                            "query_length": qlen,
                            "taxonomy": tax_name,
                            "accession": hit["accession"],
                            "description": trim_description(hit["title"]),
                            "evalue": f"{hit['evalue']:.2e}",
                            "identity_pct": hit["identity"],
                            "query_coverage_pct": hit["coverage"],
                            "output_fasta": str(out_path),
                        })
                        written += 1

                tqdm.write(f"  Written {written} sequences -> {out_path}")
                fasta_has_content = True
                save_checkpoint(checkpoint_path, checkpoint_key)
                time.sleep(1)

    print(f"\nSummary written -> {tsv_path}")
    print("Done.")


if __name__ == "__main__":
    main()
