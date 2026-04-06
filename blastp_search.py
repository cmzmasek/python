#!/usr/bin/env python3
# Copyright (C) 2026  Christian M. Zmasek
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
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
import collections
import csv
import hashlib
import io
import re
import random
import socket
import statistics
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed, TimeoutError as FutureTimeoutError
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

VERSION = "1.7.1"


@dataclass
class Hit:
    accession: str
    title: str
    evalue: float
    coverage: float
    identity: float

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
        "--api-key", default=None, metavar="KEY",
        help="NCBI API key (raises rate limit from 3 to 10 req/s; register at "
             "https://www.ncbi.nlm.nih.gov/account/)",
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
        help="Number of taxonomy BLAST searches to run in parallel per query "
             "(with --api-key you can safely raise this to 10)",
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
        "--blast-timeout", type=int, default=300,
        help="Socket timeout in seconds for NCBI connections (prevents silent hangs)",
    )
    p.add_argument(
        "--search-timeout", type=int, default=1200,
        help="Maximum total time in seconds allowed per BLAST search before retrying (default: 1200 = 20 min)",
    )
    p.add_argument(
        "--resume", action="store_true", default=False,
        help="Resume a previous run, skipping (query, taxonomy) pairs in checkpoint.txt",
    )
    p.add_argument(
        "--hit-counts", metavar="FILE", default=None,
        help="Write a continuously updated TSV matrix of hit counts "
             "(rows = taxonomies, columns = query sequences)",
    )
    p.add_argument(
        "--selected-hit-median-identity", metavar="FILE", default=None,
        help="Write a continuously updated TSV matrix of median percent identity "
             "of the selected hits (rows = taxonomies, columns = query sequences)",
    )
    p.add_argument(
        "--selected-hit-stats", metavar="FILE", default=None,
        help="Write a continuously updated TSV matrix of descriptive statistics "
             "(E-value range, identity range+mean, coverage range+mean) for the "
             "selected hits (rows = taxonomies, columns = query sequences)",
    )
    p.add_argument(
        "--selected-hit-majority-name", metavar="FILE", default=None,
        help="Write a continuously updated TSV matrix of the majority protein name "
             "among the selected hits (rows = taxonomies, columns = query sequences)",
    )
    p.add_argument(
        "--no-multispecies", action="store_true", default=False,
        help="Exclude hits whose title contains 'MULTISPECIES'",
    )
    p.add_argument(
        "--exclude-taxonomies", nargs="+", metavar="NAME", default=[],
        help="Taxonomy names to exclude from all BLAST searches, e.g. "
             "--exclude-taxonomies Coronaviridae Flaviviridae",
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
    Resolve taxonomy names to NCBI TaxIDs, canonical scientific names, and lineages.
    Returns {user_name: {"taxid": str|None, "scientific_name": str, "lineage": str}}.
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
            resolved[name] = {"taxid": None, "scientific_name": name, "lineage": name}
            continue

        if not ids:
            tqdm.write(f"  [warning] No TaxID found for '{name}', using as-is", file=sys.stderr)
            resolved[name] = {"taxid": None, "scientific_name": name, "lineage": name}
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
            lineage_raw = tax_records[0].get("Lineage", "")
            lineage = (lineage_raw.replace("; ", " > ") + " > " + sci_name
                       if lineage_raw else sci_name)
        except Exception as exc:
            tqdm.write(f"  [warning] Could not fetch canonical name for TaxID {taxid}: {exc}",
                       file=sys.stderr)
            sci_name = name
            lineage = name

        resolved[name] = {"taxid": taxid, "scientific_name": sci_name, "lineage": lineage}
        tqdm.write(f"  {name!r} -> TaxID {taxid} ({sci_name})")
        time.sleep(0.34)

    print()
    return resolved


def build_entrez_query(info: dict, exclude_taxonomies: Optional[list] = None) -> str:
    """Build a single-taxonomy Entrez query for BLAST, optionally excluding taxa."""
    query = f'"{info["scientific_name"]}"[Organism]'
    for name in (exclude_taxonomies or []):
        query += f' NOT "{name}"[Organism]'
    return query


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def trim_description(title: str, max_len: int = 120) -> str:
    """Keep only the first entry of a MULTISPECIES title and cap length."""
    first = title.split(" >")[0]
    return first if len(first) <= max_len else first[:max_len - 1] + "…"


def extract_protein_name(title: str) -> str:
    """
    Strip the accession prefix from a BLAST title, returning only the protein
    description.  Handles the common formats:
      sp|P12345|PROT_BACSU Description [Organism]  -> Description [Organism]
      gb|ABC123.1| Description [Organism]           -> Description [Organism]
      pdb|2I5L|X chain A                            -> chain A
      XP_123456.1 Description [Organism]            -> Description [Organism]
    """
    parts = title.split("|", 2)
    if len(parts) == 3:
        rest = parts[2]
        if rest.startswith(" "):
            # gb style: description follows the pipe directly
            return rest.strip()
        # sp/pdb style: first token is entry_name, skip it
        tokens = rest.split(" ", 1)
        return tokens[1].strip() if len(tokens) > 1 else rest.strip()
    if len(parts) == 2:
        return parts[1].strip()
    # plain accession: strip first whitespace-delimited token
    _, _, description = title.partition(" ")
    return description.strip()


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
               cache_dir: Optional[Path] = None, search_timeout: int = 1200):
    """
    Submit a BLASTP query to NCBI and return parsed BLAST records.
    If cache_dir is set, check for a cached XML file first and save new
    results to the cache.
    """
    cache_file: Optional[Path] = None
    if cache_dir is not None:
        key = _blast_cache_key(seq, entrez_query, evalue, hitlist_size, db)
        cache_file = cache_dir / f"{key}.xml"
        if cache_file.exists():
            tqdm.write(f"  [cache] {entrez_query[:60]}… — loading from cache")
            xml_str = cache_file.read_text()
            return list(NCBIXML.parse(io.StringIO(xml_str)))

    def _qblast():
        handle = NCBIWWW.qblast(
            "blastp", db, seq,
            entrez_query=entrez_query,
            expect=evalue,
            hitlist_size=hitlist_size,
            alignments=hitlist_size,
            descriptions=hitlist_size,
        )
        xml = handle.read()
        handle.close()
        # NCBI errors are returned as HTML/plain-text inside the handle instead
        # of raising — detect them here so retry_call can back off and retry.
        if "Error message from NCBI" in xml or "error code:" in xml.lower():
            raise RuntimeError(f"NCBI error response: {xml[:200].strip()}")
        return xml

    def _qblast_with_timeout():
        ex = ThreadPoolExecutor(max_workers=1)
        future = ex.submit(_qblast)
        try:
            result = future.result(timeout=search_timeout)
            ex.shutdown(wait=False)
            return result
        except FutureTimeoutError:
            ex.shutdown(wait=False)
            raise RuntimeError(
                f"BLAST search exceeded {search_timeout}s total time limit"
            )

    xml_str = retry_call(
        _qblast_with_timeout,
        max_attempts=max_retries,
        base_delay=retry_delay,
        label="BLASTP",
    )

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
    passing: dict[str, Hit] = {}

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
        if acc not in passing:
            passing[acc] = Hit(
                accession=acc,
                title=alignment.title,
                evalue=best_hsp.expect,
                coverage=round(cov, 1),
                identity=round(identity_pct, 1),
            )

    return list(passing.values())


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
def make_rng(tax: str, seed: Optional[int]) -> random.Random:
    """Return a per-taxonomy RNG seeded deterministically when seed is given."""
    if seed is None:
        return random.Random()
    return random.Random(seed ^ (hash(tax) & 0xFFFFFFFF))


# Per-(query, taxonomy) worker
# ---------------------------------------------------------------------------

def process_taxonomy(
    qrec,
    tax_name: str,
    tax_info: dict,
    args,
    cache_dir: Optional[Path],
    rng: random.Random,
) -> tuple[str, list[tuple], int]:
    """
    Run BLAST, filter, and fetch sequences for one (query, taxonomy) pair.
    Returns (tax_name, [(hit_dict, SeqRecord_or_None), ...], n_passing).
    n_passing is the count of hits passing all filters (before random draw).
    Called from worker threads.
    """
    entrez_q = build_entrez_query(tax_info, args.exclude_taxonomies)
    qlen = len(qrec.seq)

    jitter = rng.uniform(2.0, 5.0)
    tqdm.write(f"  [{tax_name}] waiting {jitter:.1f}s before BLAST...")
    time.sleep(jitter)
    tqdm.write(f"  [{tax_name}] starting BLAST...")

    try:
        blast_records = run_blastp(
            str(qrec.seq), entrez_q, args.evalue,
            args.hitlist_size, args.db, args.max_retries, args.retry_delay,
            cache_dir=cache_dir, search_timeout=args.search_timeout,
        )
    except Exception as exc:
        tqdm.write(f"  [{tax_name}] BLAST failed after retries: {exc}", file=sys.stderr)
        return tax_name, [], 0

    if not blast_records or not blast_records[0].alignments:
        tqdm.write(f"  [{tax_name}] no BLAST hits returned")
        return tax_name, [], 0

    passing = filter_hits(
        blast_records[0], args.evalue, args.coverage, args.identity, qlen,
        exclude_multispecies=args.no_multispecies,
    )
    tqdm.write(f"  [{tax_name}] hits passing filters: {len(passing)}")

    if not passing:
        return tax_name, [], 0

    max_seqs = args.taxonomy_max_seqs.get(tax_name, args.max_seqs)
    selected = rng.sample(passing, min(max_seqs, len(passing)))
    tqdm.write(f"  [{tax_name}] selected: {len(selected)}")

    sequences = fetch_sequences(
        [h.accession for h in selected], args.max_retries, args.retry_delay,
    )

    results = [(hit, match_sequence(hit.accession, sequences)) for hit in selected]
    return tax_name, results, len(passing)


# ---------------------------------------------------------------------------
# Hit-count matrix
# ---------------------------------------------------------------------------

def write_tsv_matrix(path: Path, data: dict, query_ids: list, tax_names: list,
                     lineages: dict) -> None:
    """Rewrite a TSV matrix (rows = taxa, columns = queries)."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["taxonomy", "lineage"] + query_ids)
        for tax in tax_names:
            row = [tax, lineages.get(tax, "")] + [data.get(tax, {}).get(qid, "") for qid in query_ids]
            writer.writerow(row)


def format_hit_stats(results: list) -> str:
    """
    Format descriptive statistics for the selected hits as a compact string.
    Example: 'n=5 E:1.2e-08–3.5e-06 id:80.0–95.3%(avg:87.1%) cov:88.0–100.0%(avg:94.2%)'
    """
    if not results:
        return "n=0"
    evalues   = [hit.evalue   for hit, _ in results]
    ids       = [hit.identity for hit, _ in results]
    coverages = [hit.coverage for hit, _ in results]
    e_min, e_max     = min(evalues),   max(evalues)
    id_min, id_max   = min(ids),       max(ids)
    cov_min, cov_max = min(coverages), max(coverages)
    id_median = round(statistics.median(ids),  1)
    cov_mean  = round(statistics.mean(coverages), 1)
    return (
        f"n={len(results)} "
        f"E:{e_min:.1e}\u2013{e_max:.1e} "
        f"id:{id_min}\u2013{id_max}%(med:{id_median}%) "
        f"cov:{cov_min}\u2013{cov_max}%(avg:{cov_mean}%)"
    )




def majority_name(results: list) -> str:
    """Return the most common trimmed protein name among the selected hits."""
    if not results:
        return ""
    names = [re.sub(r'\s*\[[^\]]*\]\s*$', '', extract_protein_name(trim_description(hit.title))).strip()
             for hit, _ in results]
    return collections.Counter(names).most_common(1)[0][0]




# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    socket.setdefaulttimeout(args.blast_timeout)

    base_seed = args.seed

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Set up cache directory
    cache_dir: Optional[Path] = None
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

    hit_counts_path      = Path(args.hit_counts)                    if args.hit_counts                    else None
    median_identity_path = Path(args.selected_hit_median_identity)  if args.selected_hit_median_identity  else None
    hit_stats_path       = Path(args.selected_hit_stats)            if args.selected_hit_stats            else None
    majority_name_path   = Path(args.selected_hit_majority_name)    if args.selected_hit_majority_name    else None

    query_ids = [rec.id for rec in query_records]
    lineages = {name: resolved[name]["lineage"] for name in args.taxonomies}
    hit_counts: dict = {tax: {} for tax in args.taxonomies}
    median_identity: dict = {tax: {} for tax in args.taxonomies}
    hit_stats: dict = {tax: {} for tax in args.taxonomies}
    majority_names: dict = {tax: {} for tax in args.taxonomies}

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
    print(f"NCBI API key  : {'yes' if args.api_key else 'no (max 3 req/s)'}")
    print(f"BLAST timeout : {args.blast_timeout}s (socket) / {args.search_timeout}s (per search)")
    print(f"Cache         : {'disabled' if cache_dir is None else str(cache_dir)}")
    print(f"E-value       : {args.evalue}")
    print(f"Coverage      : {args.coverage}%")
    print(f"Identity      : {args.identity}%")
    if args.exclude_taxonomies:
        print(f"Excluding     : {', '.join(args.exclude_taxonomies)}")
    print(f"Max seqs/tax  : {args.max_seqs}")
    print(f"BLAST pool    : {args.hitlist_size} hits per search")
    print(f"Output dir    : {outdir}")
    print(f"Summary TSV   : {tsv_path}")
    print(f"Hit counts    : {args.hit_counts or 'disabled'}")
    print(f"Median id     : {args.selected_hit_median_identity or 'disabled'}")
    print(f"Hit stats     : {args.selected_hit_stats or 'disabled'}")
    print(f"Majority name : {args.selected_hit_majority_name or 'disabled'}")
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

            n_workers = min(args.workers, len(pending_taxes))
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(
                        process_taxonomy,
                        qrec, tax_name, resolved[tax_name],
                        args, cache_dir, make_rng(tax_name, base_seed),
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
                        _, results, n_passing = future.result()
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
                                        f"  [warning] Sequence not fetched: {hit.accession}",
                                        file=sys.stderr,
                                    )
                                    writer.writerow({
                                        "query_id": qid,
                                        "query_length": qlen,
                                        "taxonomy": tax_name,
                                        "accession": hit.accession,
                                        "description": trim_description(hit.title),
                                        "evalue": f"{hit.evalue:.2e}",
                                        "identity_pct": hit.identity,
                                        "query_coverage_pct": hit.coverage,
                                        "output_fasta": "NOT_FETCHED",
                                    })
                                    continue
                                rec.description = (
                                    f"evalue={hit.evalue:.2e} "
                                    f"id={hit.identity}% "
                                    f"qcov={hit.coverage}% | {rec.description}"
                                )
                                SeqIO.write(rec, fh, "fasta")
                                writer.writerow({
                                    "query_id": qid,
                                    "query_length": qlen,
                                    "taxonomy": tax_name,
                                    "accession": hit.accession,
                                    "description": trim_description(hit.title),
                                    "evalue": f"{hit.evalue:.2e}",
                                    "identity_pct": hit.identity,
                                    "query_coverage_pct": hit.coverage,
                                    "output_fasta": str(out_path),
                                })
                                written += 1
                        tsv_fh.flush()
                        save_checkpoint(checkpoint_path, checkpoint_key)
                        hit_counts[tax_name][qid] = n_passing
                        if hit_counts_path:
                            write_tsv_matrix(hit_counts_path, hit_counts, query_ids, args.taxonomies, lineages)
                        identities = [hit.identity for hit, _ in results]
                        median_identity[tax_name][qid] = round(statistics.median(identities), 1) if identities else 0
                        if median_identity_path:
                            write_tsv_matrix(median_identity_path, median_identity, query_ids, args.taxonomies, lineages)
                        hit_stats[tax_name][qid] = format_hit_stats(results)
                        if hit_stats_path:
                            write_tsv_matrix(hit_stats_path, hit_stats, query_ids, args.taxonomies, lineages)
                        majority_names[tax_name][qid] = majority_name(results)
                        if majority_name_path:
                            write_tsv_matrix(majority_name_path, majority_names, query_ids, args.taxonomies, lineages)
                        tqdm.write(f"  [{tax_name}] written {written} sequences -> {out_path}")

                    tax_bar.update(1)
                    tax_bar.set_postfix(last=tax_name)
                tax_bar.close()

    print(f"\nSummary written -> {tsv_path}")
    if hit_counts_path:
        print(f"Hit counts written -> {hit_counts_path}")
    if median_identity_path:
        print(f"Median identity written -> {median_identity_path}")
    if hit_stats_path:
        print(f"Hit stats written -> {hit_stats_path}")
    if majority_name_path:
        print(f"Majority names written -> {majority_name_path}")
    print("Done.")


if __name__ == "__main__":
    main()
