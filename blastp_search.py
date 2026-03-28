#!/usr/bin/env python3
"""
BLASTP search program with per-taxonomy search strategy.

Takes a list of taxonomy names and a FASTA file of query protein sequences,
runs one BLASTP search per (query, taxonomy) pair against NCBI, filters by
E-value, query coverage, and percent identity, and writes a random selection
of matching sequences to per-query FASTA output files.

Running one BLAST per taxonomy ensures that each taxon gets its own hit-list
pool, preventing heavily-annotated species from crowding out rarer ones.

Taxonomy searches for a given query are run in parallel (--workers) and
BLAST XML results are cached on disk to avoid redundant NCBI round-trips
when re-running with different filter parameters.

Requirements:
    pip install biopython tqdm
"""

import argparse
import csv
import hashlib
import io
import re
import random
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

VERSION = "1.2.1"

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
        help="File with one taxonomy name per line (can be combined with -t). "
             "Optionally add a tab-separated integer to override -n for that taxonomy.",
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
        "--workers", type=int, default=3,
        help="Number of taxonomy BLAST searches to run in parallel per query",
    )
    p.add_argument(
        "--no-cache", action="store_true", default=False,
        help="Disable BLAST result caching (always query NCBI)",
    )
    p.add_argument(
        "--cache-dir", metavar="DIR", default=None,
        help="Directory for BLAST XML cache (default: <outdir>/cache)",
    )
    p.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducible sequence selection. "
             "With --workers > 1 the draw order varies by thread scheduling; "
             "use --workers 1 for fully reproducible runs.",
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
    p.add_argument(
        "--no-multispecies", action="store_true", default=False,
        help="Exclude hits whose title contains 'MULTISPECIES'",
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
    return f'"{info["scientific_name"]}"[Organism]'


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
      'pdb|2I5L|X ...'
    or plain accessions like 'XP_123456.1'.
    """
    parts = raw.split("|")
    if len(parts) >= 3 and parts[0].lower() == "pdb":
        return f"{parts[1]}_{parts[2].split()[0]}"   # PDB: '2I5L_X'
    if len(parts) >= 3:
        return parts[1]          # UniProt style: db|acc|entry
    if len(parts) == 2:
        return parts[-1]         # gb|acc  or similar
    return parts[0].split()[0]   # plain accession, drop description


_PDB_RE = re.compile(r'^[0-9][A-Z0-9]{3}_[A-Z0-9]+$', re.IGNORECASE)


def is_pdb_accession(acc: str) -> bool:
    """Return True if acc looks like a PDB accession (e.g. '2I5L_X')."""
    return bool(_PDB_RE.match(acc))


def _fetch_single_pdb(acc: str, max_retries: int, retry_delay: float) -> dict:
    """
    Fetch a single PDB accession from NCBI by searching the protein database
    for its UID, then fetching the FASTA. Returns {record.id: SeqRecord}.
    """
    try:
        handle = retry_call(
            Entrez.esearch, db="protein", term=f"{acc}[accn]",
            max_attempts=max_retries, base_delay=retry_delay,
            label=f"esearch PDB {acc}",
        )
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            tqdm.write(f"  [warning] PDB accession not found in NCBI protein: {acc}",
                       file=sys.stderr)
            return {}
        uid = record["IdList"][0]
        handle = retry_call(
            Entrez.efetch, db="protein", id=uid, rettype="fasta", retmode="text",
            max_attempts=max_retries, base_delay=retry_delay,
            label=f"efetch PDB {acc}",
        )
        seqs = {rec.id: rec for rec in SeqIO.parse(handle, "fasta")}
        handle.close()
        return seqs
    except Exception as exc:
        tqdm.write(f"  [warning] Failed to fetch PDB accession {acc}: {exc}", file=sys.stderr)
        return {}


# ---------------------------------------------------------------------------
# BLAST with caching
# ---------------------------------------------------------------------------

def _blast_cache_key(seq: str, entrez_query: str, evalue: float,
                     hitlist_size: int, db: str) -> str:
    payload = f"{seq}|{entrez_query}|{evalue}|{hitlist_size}|{db}"
    return hashlib.md5(payload.encode()).hexdigest()


def run_blastp(seq: str, entrez_query: str, evalue: float, hitlist_size: int,
               db: str, max_retries: int, retry_delay: float,
               cache_dir: Path | None = None):
    """
    Submit a BLASTP query to NCBI and return parsed BLAST records.
    If cache_dir is set, check for a cached XML file first and save new
    results to the cache.
    """
    cache_file: Path | None = None
    if cache_dir is not None:
        key = _blast_cache_key(seq, entrez_query, evalue, hitlist_size, db)
        cache_file = cache_dir / f"{key}.xml"
        if cache_file.exists():
            tqdm.write(f"  [cache] {entrez_query[:60]}… — loading from cache")
            xml_str = cache_file.read_text()
            return list(NCBIXML.parse(io.StringIO(xml_str)))

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
    xml_str = result_handle.read()
    result_handle.close()

    if cache_file is not None:
        cache_file.write_text(xml_str)

    return list(NCBIXML.parse(io.StringIO(xml_str)))


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def filter_hits(blast_record, evalue_cutoff: float, coverage_cutoff: float,
                identity_cutoff: float, query_length: int,
                exclude_multispecies: bool = False) -> list[dict]:
    """
    Return one entry per unique accession that passes all filters.
    Coverage is computed from merged HSP intervals; E-value and identity
    are taken from the best (first) HSP.
    """
    passing: list[dict] = []
    seen: set[str] = set()

    for alignment in blast_record.alignments:
        if exclude_multispecies and "MULTISPECIES" in alignment.title:
            continue

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
    Regular accessions are fetched in batches; PDB accessions (e.g. '2I5L_X')
    are resolved individually via esearch + efetch.
    Returns {record.id: SeqRecord}.
    """
    regular = [a for a in accessions if not is_pdb_accession(a)]
    pdb     = [a for a in accessions if is_pdb_accession(a)]

    sequences: dict = {}

    # Batch-fetch regular accessions
    batches = range(0, len(regular), batch_size)
    for i, start in enumerate(batches):
        batch = regular[start:start + batch_size]
        try:
            handle = retry_call(
                Entrez.efetch,
                db="protein", id=",".join(batch), rettype="fasta", retmode="text",
                max_attempts=max_retries, base_delay=retry_delay,
                label=f"Entrez fetch batch {i + 1}",
            )
            for rec in SeqIO.parse(handle, "fasta"):
                sequences[rec.id] = rec
            handle.close()
        except Exception as exc:
            tqdm.write(
                f"  [warning] Entrez fetch failed for batch {i + 1}: {exc}",
                file=sys.stderr,
            )
        if start + batch_size < len(regular):
            time.sleep(0.4)

    # Individually resolve PDB accessions via esearch
    for acc in pdb:
        tqdm.write(f"  [PDB] fetching {acc} via esearch...")
        seqs = _fetch_single_pdb(acc, max_retries, retry_delay)
        sequences.update(seqs)
        time.sleep(0.34)

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
# Per-(query, taxonomy) worker
# ---------------------------------------------------------------------------

def process_taxonomy(
    qrec,
    tax_name: str,
    tax_info: dict,
    args,
    cache_dir: Path | None,
    rng: random.Random,
) -> tuple[str, list[tuple]]:
    """
    Run BLAST, filter, and fetch sequences for one (query, taxonomy) pair.
    Returns (tax_name, [(hit_dict, SeqRecord_or_None), ...]).
    Called from worker threads.
    """
    entrez_q = build_entrez_query(tax_info)
    qlen = len(qrec.seq)

    tqdm.write(f"  [{tax_name}] starting BLAST...")

    try:
        blast_records = run_blastp(
            str(qrec.seq), entrez_q, args.evalue,
            args.hitlist_size, args.db, args.max_retries, args.retry_delay,
            cache_dir=cache_dir,
        )
    except Exception as exc:
        tqdm.write(f"  [{tax_name}] BLAST failed after retries: {exc}", file=sys.stderr)
        return tax_name, []

    if not blast_records or not blast_records[0].alignments:
        tqdm.write(f"  [{tax_name}] no BLAST hits returned")
        return tax_name, []

    passing = filter_hits(
        blast_records[0], args.evalue, args.coverage, args.identity, qlen,
        exclude_multispecies=args.no_multispecies,
    )
    tqdm.write(f"  [{tax_name}] hits passing filters: {len(passing)}")

    if not passing:
        return tax_name, []

    max_seqs = args.taxonomy_max_seqs.get(tax_name, args.max_seqs)
    selected = rng.sample(passing, min(max_seqs, len(passing)))
    tqdm.write(f"  [{tax_name}] selected: {len(selected)}")

    sequences = fetch_sequences(
        [h["accession"] for h in selected], args.max_retries, args.retry_delay,
    )

    results = [(hit, match_sequence(hit["accession"], sequences)) for hit in selected]
    return tax_name, results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    Entrez.email = args.email

    base_seed = args.seed

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Set up cache directory
    cache_dir: Path | None = None
    if not args.no_cache:
        cache_dir = Path(args.cache_dir) if args.cache_dir else outdir / "cache"
        cache_dir.mkdir(parents=True, exist_ok=True)

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
    print(f"Workers       : {args.workers}  (parallel taxonomy searches per query)")
    print(f"Cache         : {'disabled' if cache_dir is None else str(cache_dir)}")
    print(f"E-value       : {args.evalue}")
    print(f"Coverage      : {args.coverage}%")
    print(f"Identity      : {args.identity}%")
    print(f"Max seqs/tax  : {args.max_seqs}")
    print(f"BLAST pool    : {args.hitlist_size} hits per search")
    print(f"Output dir    : {outdir}")
    print(f"Summary TSV   : {tsv_path}")
    print()

    write_lock = threading.Lock()

    with open(tsv_path, tsv_mode, newline="") as tsv_fh:
        writer = csv.DictWriter(tsv_fh, fieldnames=tsv_columns, delimiter="\t")
        if tsv_mode == "w":
            writer.writeheader()

        query_bar = tqdm(query_records, desc="Queries", unit="query")
        for qrec in query_bar:
            qid = qrec.id
            qlen = len(qrec.seq)
            query_bar.set_description(f"Query: {qid}")

            pending_taxes = [
                name for name in resolved
                if f"{qid}::{name}" not in completed
            ]
            if not pending_taxes:
                tqdm.write(f"\n=== Query: {qid} — all taxonomies done, skipping ===")
                continue

            tqdm.write(
                f"\n=== Query: {qid} ({qlen} aa) | "
                f"{len(pending_taxes)} taxonomy search(es) | "
                f"{min(args.workers, len(pending_taxes))} worker(s) ==="
            )

            out_path = outdir / f"{safe_filename(qid)}.fasta"

            # Build per-taxonomy RNGs so parallel draws are independent but
            # deterministic when --seed is given
            def make_rng(tax: str) -> random.Random:
                if base_seed is None:
                    return random.Random()
                return random.Random(base_seed ^ (hash(tax) & 0xFFFFFFFF))

            n_workers = min(args.workers, len(pending_taxes))
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(
                        process_taxonomy,
                        qrec, tax_name, resolved[tax_name],
                        args, cache_dir, make_rng(tax_name),
                    ): tax_name
                    for tax_name in pending_taxes
                }

                tax_bar = tqdm(
                    total=len(futures),
                    desc=f"  Taxonomies [{qid}]",
                    unit="tax",
                    leave=False,
                )
                for future in as_completed(futures):
                    tax_name = futures[future]
                    checkpoint_key = f"{qid}::{tax_name}"
                    try:
                        _, results = future.result()
                    except Exception as exc:
                        tqdm.write(f"  [{tax_name}] unexpected error: {exc}", file=sys.stderr)
                        tax_bar.update(1)
                        continue

                    # Write results under lock (FASTA + TSV + checkpoint)
                    with write_lock:
                        fasta_mode = "a" if (out_path.exists() and out_path.stat().st_size > 0) else "w"
                        written = 0
                        with open(out_path, fasta_mode) as fh:
                            for hit, rec in results:
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
                        tsv_fh.flush()
                        save_checkpoint(checkpoint_path, checkpoint_key)
                        tqdm.write(f"  [{tax_name}] written {written} sequences -> {out_path}")

                    tax_bar.update(1)
                    tax_bar.set_postfix(last=tax_name)
                tax_bar.close()

    print(f"\nSummary written -> {tsv_path}")
    print("Done.")


if __name__ == "__main__":
    main()
