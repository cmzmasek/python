#!/usr/bin/env python3
"""
fasta_relabel.py — Relabel FASTA description lines via UniProt / NCBI lookups.

Each sequence's accession is detected automatically (or forced with --db), then
the description line is rebuilt from a user-supplied format string such as:

    "{organism} | {id} | {name} | {nucl_id}"
    "{lineage}; {organism}; {id}"

A detailed report on labeling coverage, sequence-length statistics, and the
distribution of organisms / proteins is always generated after processing.
"""

import argparse
import dataclasses
import json
import math
import re
import statistics
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, Optional

VERSION = "1.0.0"

# ---------------------------------------------------------------------------
# Rate-limiting helpers
# ---------------------------------------------------------------------------

_last_request_time: float = 0.0
NCBI_INTERVAL = 0.34  # ≤3 req/s without an API key
UNIPROT_INTERVAL = 0.1


def _throttle(interval: float) -> None:
    global _last_request_time
    wait = interval - (time.monotonic() - _last_request_time)
    if wait > 0:
        time.sleep(wait)
    _last_request_time = time.monotonic()


def _get(url: str, timeout: int = 30) -> str:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "fasta_relabel/1.0 (bioinformatics tool)"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------

def parse_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (header_without_>, sequence) pairs from a FASTA file."""
    header: Optional[str] = None
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:]
                parts = []
            elif line:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def write_fasta(
        path: Path, records: list[tuple[str, str]], wrap: int = 60, overwrite: bool = True
) -> bool:
    """Write FASTA records to *path*; return False (and skip) if file exists and overwrite=False."""
    if path.exists() and not overwrite:
        print(
            f"  Skipping {path}: file already exists (use --overwrite to replace)",
            file=sys.stderr,
        )
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i: i + wrap] + "\n")
    return True


# ---------------------------------------------------------------------------
# Accession detection
# ---------------------------------------------------------------------------

# UniProt: 6- or 10-character accession per the UniProt spec
_UNIPROT_RE = re.compile(
    r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b"
)
# NCBI RefSeq with two-letter prefix+underscore: NP_, XP_, WP_, YP_, AP_, …
_REFSEQ_RE = re.compile(r"\b([A-Z]{2}_\d{6,9}(?:\.\d+)?)\b")
# Bare GenBank protein accession: 1–2 letters + 5–8 digits
_GENBANK_RE = re.compile(r"\b([A-Z]{1,2}\d{5,8}(?:\.\d+)?)\b")
# Legacy gi| pipe-format
_GI_RE = re.compile(r"gi\|(\d+)")


def detect_id(header: str, db_hint: str) -> tuple[str, str]:
    """Return (db, accession) inferred from a FASTA header line."""
    first_token = header.split()[0]

    # UniProt pipe format: sp|ACC|ENTRYNAME or tr|ACC|ENTRYNAME
    m = re.match(r"(?:sp|tr)\|([A-Z0-9]+)\|", first_token)
    if m and db_hint in ("auto", "uniprot"):
        return "uniprot", m.group(1)

    # Bare UniProt accession as the first token
    m = _UNIPROT_RE.match(first_token)
    if m and db_hint in ("auto", "uniprot"):
        return "uniprot", m.group(1)

    if db_hint in ("auto", "ncbi"):
        m = _GI_RE.search(header)
        if m:
            return "ncbi", m.group(1)
        m = _REFSEQ_RE.search(first_token)
        if m:
            return "ncbi", m.group(1)
        m = _GENBANK_RE.match(first_token)
        if m:
            return "ncbi", m.group(1)

    fallback = first_token.split("|")[-1]
    return (db_hint if db_hint != "auto" else "ncbi"), fallback


# ---------------------------------------------------------------------------
# Nucleotide ID helpers
# ---------------------------------------------------------------------------

# Molecule-type preference orders for --nucl-type
_NUCL_TYPE_PREF: dict[str, list[str]] = {
    "genomic": ["Genomic_DNA", "Genomic_RNA", "mRNA", "Transcribed_RNA", "Other_RNA"],
    "mrna": ["mRNA", "Transcribed_RNA", "Other_RNA", "Genomic_DNA", "Genomic_RNA"],
    # "any" = first EMBL entry regardless of type (handled separately)
}


def _nucl_id_from_uniprot(data: dict, nucl_type: str) -> tuple[str, list[str]]:
    """
    Return (chosen_nucl_id, all_nucl_ids) from UniProt EMBL cross-references.

    nucl_type: 'genomic' | 'mrna' | 'any'
    chosen_nucl_id is a single accession selected by nucl_type preference.
    all_nucl_ids is the complete list regardless of preference.
    """
    xrefs = data.get("uniProtKBCrossReferences", [])
    embl_refs = [x for x in xrefs if x.get("database") == "EMBL"]
    all_ids = [x["id"] for x in embl_refs if x.get("id")]

    if not embl_refs:
        return "", []

    if nucl_type == "any":
        return embl_refs[0]["id"], all_ids

    def mol_type(xref: dict) -> str:
        for prop in xref.get("properties", []):
            if prop.get("key") == "MoleculeType":
                return prop.get("value", "")
        return ""

    # Build a dict: molecule_type → first accession of that type
    by_type: dict[str, str] = {}
    for ref in embl_refs:
        mt = mol_type(ref)
        if mt not in by_type:
            by_type[mt] = ref["id"]

    pref_order = _NUCL_TYPE_PREF.get(nucl_type, [])
    for pref in pref_order:
        if pref in by_type:
            return by_type[pref], all_ids

    return embl_refs[0]["id"], all_ids  # ultimate fallback


def _nucl_id_from_ncbi_xml(gbseq: ET.Element) -> str:
    """
    Extract the nucleotide accession from a GenBank CDS 'coded_by' qualifier.

    coded_by examples:
      NM_000558.5:1..429
      complement(AL096865.18:c45836..45393)
      join(NM_000558.5:1..200,NM_000558.5:300..429)
    """
    for feat in gbseq.findall(".//GBFeature"):
        if feat.findtext("GBFeature_key") != "CDS":
            continue
        for qual in feat.findall(".//GBQualifier"):
            if qual.findtext("GBQualifier_name") == "coded_by":
                coded_by = qual.findtext("GBQualifier_value") or ""
                # Strip wrapper functions so the accession is the first token
                stripped = re.sub(r"(?:complement|join|order)\(", "", coded_by)
                m = re.search(r"([A-Za-z]{1,2}_?[0-9][A-Za-z0-9_.]*[0-9]):", stripped)
                if m:
                    return m.group(1)
    return ""


# ---------------------------------------------------------------------------
# UniProt query
# ---------------------------------------------------------------------------

def _query_uniprot(accession: str, nucl_type: str) -> dict:
    _throttle(UNIPROT_INTERVAL)
    data = json.loads(_get(f"https://rest.uniprot.org/uniprotkb/{accession}.json"))

    result: dict = {
        "id": accession,
        "protein_name": "",
        "gene": "",
        "organism": "",
        "taxid": "",
        "lineage": "",
        "entry_name": data.get("uniProtkbId", accession),
        "nucl_id": "",
        "nucl_ids": [],  # full list, stored for JSON report
    }

    # Protein name
    pd_ = data.get("proteinDescription", {})
    rec = pd_.get("recommendedName") or (pd_.get("submissionNames") or [{}])[0]
    result["protein_name"] = (rec.get("fullName") or {}).get("value", "")

    # Gene name
    genes = data.get("genes", [])
    if genes:
        result["gene"] = (genes[0].get("geneName") or {}).get("value", "")

    # Organism + lineage
    org = data.get("organism", {})
    result["organism"] = org.get("scientificName", "")
    result["taxid"] = str(org.get("taxonId", ""))
    result["lineage"] = "; ".join(org.get("lineage", []))

    # Nucleotide / genome ID from EMBL cross-references
    nucl_id, all_ids = _nucl_id_from_uniprot(data, nucl_type)
    result["nucl_id"] = nucl_id
    result["nucl_ids"] = all_ids

    return result


# ---------------------------------------------------------------------------
# NCBI query
# ---------------------------------------------------------------------------

def _query_ncbi(
        accession: str,
        email: Optional[str],
        api_key: Optional[str],
) -> dict:
    _throttle(NCBI_INTERVAL)
    params: dict = {"db": "protein", "id": accession, "rettype": "gb", "retmode": "xml"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            + urllib.parse.urlencode(params)
    )
    return _parse_ncbi_xml(_get(url), accession)


def _parse_ncbi_xml(xml_text: str, accession: str) -> dict:
    result: dict = {
        "id": accession,
        "protein_name": "",
        "gene": "",
        "organism": "",
        "taxid": "",
        "lineage": "",
        "entry_name": accession,
        "nucl_id": "",
        "nucl_ids": [],
    }
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return result

    gbseq = root.find(".//GBSeq")
    if gbseq is None:
        return result

    av = gbseq.findtext("GBSeq_accession-version")
    if av:
        result["entry_name"] = av

    result["organism"] = gbseq.findtext(".//GBSeq_organism") or ""
    result["lineage"] = gbseq.findtext(".//GBSeq_taxonomy") or ""

    for feat in gbseq.findall(".//GBFeature"):
        key = feat.findtext("GBFeature_key")
        quals = {
            q.findtext("GBQualifier_name"): q.findtext("GBQualifier_value")
            for q in feat.findall(".//GBQualifier")
        }
        if key == "source":
            xref = quals.get("db_xref", "") or ""
            if xref.startswith("taxon:"):
                result["taxid"] = xref.split(":", 1)[1]
        if key == "CDS":
            if not result["protein_name"] and quals.get("product"):
                result["protein_name"] = quals["product"]
            if not result["gene"] and quals.get("gene"):
                result["gene"] = quals["gene"]
        if key == "gene" and not result["gene"] and quals.get("gene"):
            result["gene"] = quals["gene"]

    if not result["protein_name"]:
        defn = gbseq.findtext("GBSeq_definition") or ""
        result["protein_name"] = re.sub(r"\s*\[.+?\]$", "", defn).strip()

    # Nucleotide accession from CDS coded_by qualifier
    nucl_id = _nucl_id_from_ncbi_xml(gbseq)
    result["nucl_id"] = nucl_id
    result["nucl_ids"] = [nucl_id] if nucl_id else []

    return result


# ---------------------------------------------------------------------------
# In-memory cache (avoid re-querying the same accession in a batch)
# ---------------------------------------------------------------------------

_cache: dict[str, dict] = {}


def fetch_info(
        accession: str,
        db: str,
        email: Optional[str],
        api_key: Optional[str],
        nucl_type: str,
        verbose: bool,
) -> tuple[Optional[dict], str]:
    """Return (info_dict_or_None, error_string)."""
    # Include nucl_type in the cache key so changing --nucl-type always
    # produces the correct result (UniProt re-picked; NCBI unaffected).
    key = f"{db}:{accession}:{nucl_type}"
    if key in _cache:
        return _cache[key], ""
    try:
        if verbose:
            print(f"  [{db}] querying {accession} …", file=sys.stderr)
        if db == "uniprot":
            info = _query_uniprot(accession, nucl_type)
        else:
            info = _query_ncbi(accession, email, api_key)
        _cache[key] = info
        return info, ""
    except urllib.error.HTTPError as exc:
        msg = f"HTTP {exc.code} {exc.reason}"
        print(f"  Warning: {msg} for {accession}", file=sys.stderr)
        return None, msg
    except Exception as exc:  # noqa: BLE001
        msg = str(exc)
        print(f"  Warning: could not fetch {accession}: {msg}", file=sys.stderr)
        return None, msg


# ---------------------------------------------------------------------------
# Header formatting
# ---------------------------------------------------------------------------

_ALIASES: dict[str, str] = {
    "id": "id",
    "accession": "id",
    "acc": "id",
    "name": "protein_name",
    "protein_name": "protein_name",
    "protein": "protein_name",
    "organism": "organism",
    "taxon": "organism",
    "species": "organism",
    "sci_name": "organism",
    "gene": "gene",
    "lineage": "lineage",
    "taxonomy": "lineage",
    "taxid": "taxid",
    "tax_id": "taxid",
    "entry": "entry_name",
    "entry_name": "entry_name",
    # Nucleotide / genome ID
    "nucl_id": "nucl_id",
    "genome_id": "nucl_id",
    "nucleotide_id": "nucl_id",
    "nucl_acc": "nucl_id",
    "coding_seq": "nucl_id",
}


def format_header(fmt: str, info: dict, use_underscores: bool = True) -> str:
    out = fmt
    for alias, field in _ALIASES.items():
        out = out.replace(f"{{{alias}}}", info.get(field, "") or "")
    out = out.strip()
    if use_underscores:
        out = re.sub(r"\s+", "_", out)
    return out


# ---------------------------------------------------------------------------
# Per-record result
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class RecordResult:
    source_file: str  # input filename (basename)
    original_header: str
    accession: str
    db_used: str  # 'uniprot' | 'ncbi'
    seq_len: int
    relabeled: bool  # False if lookup failed
    new_header: str
    organism: str = ""
    protein_name: str = ""
    gene: str = ""
    taxid: str = ""
    lineage: str = ""
    nucl_id: str = ""  # chosen nucleotide/genome accession
    nucl_ids: dataclasses.field(default_factory=list) = dataclasses.field(  # type: ignore[assignment]
        default_factory=list
    )  # full list of all linked nucleotide accessions
    error: str = ""  # non-empty when lookup failed


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def _n50(lengths: list[int]) -> int:
    """Sequence N50: smallest L s.t. sequences ≥ L account for ≥50% of total."""
    if not lengths:
        return 0
    total = sum(lengths)
    cumsum = 0
    for ln in sorted(lengths, reverse=True):
        cumsum += ln
        if cumsum >= total / 2:
            return ln
    return lengths[-1]


def _histogram(values: list[int], n_bins: int = 10) -> list[tuple[int, int, int]]:
    """Return [(bin_lo, bin_hi, count), …] with equal-width bins."""
    if not values:
        return []
    lo, hi = min(values), max(values)
    if lo == hi:
        return [(lo, hi, len(values))]
    bin_size = max(1, math.ceil((hi - lo + 1) / n_bins))
    bins = [0] * n_bins
    for v in values:
        idx = min((v - lo) // bin_size, n_bins - 1)
        bins[idx] += 1
    return [(lo + i * bin_size, lo + (i + 1) * bin_size - 1, cnt) for i, cnt in enumerate(bins)]


def _bar(fraction: float, width: int = 28) -> str:
    return "█" * round(fraction * width)


def _fmt_int(n: int) -> str:
    return f"{n:,}"


def _fmt_float(f: float, dp: int = 1) -> str:
    return f"{f:,.{dp}f}"


# ---------------------------------------------------------------------------
# Report building
# ---------------------------------------------------------------------------

RULE_HEAVY = "═" * 72
RULE_LIGHT = "─" * 72
COL = 72


def _center(text: str) -> str:
    return text.center(COL)


def _section(title: str) -> str:
    return f"\n{RULE_LIGHT}\n {title}\n{RULE_LIGHT}"


def build_report_text(
        results: list[RecordResult],
        fmt: str,
        db: str,
        nucl_type: str,
        top_n: int,
        input_files: list[Path],
        elapsed: float,
) -> str:
    lines: list[str] = []
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    # ── Header ──────────────────────────────────────────────────────────────
    lines += [
        RULE_HEAVY,
        _center(f"fasta_relabel {VERSION}  •  Run Report"),
        _center(now),
        RULE_HEAVY,
        f"  Format    : {fmt}",
        f"  Database  : {db}",
        f"  Nucl type : {nucl_type}",
        f"  Files     : {len(input_files)} processed",
        f"  Elapsed   : {elapsed:.1f} s",
    ]

    if not results:
        lines += ["", "  (no records found)", "", RULE_HEAVY]
        return "\n".join(lines)

    total = len(results)
    relabeled = sum(1 for r in results if r.relabeled)
    failed = total - relabeled
    unique_acc = len({r.accession for r in results})
    with_nucl = sum(1 for r in results if r.relabeled and r.nucl_id)
    lengths = [r.seq_len for r in results]

    # ── Labeling coverage ───────────────────────────────────────────────────
    lines += [_section("LABELING COVERAGE")]
    lines += [
        f"  {'Total sequences':<38} {_fmt_int(total):>6}",
        f"  {'Successfully relabeled':<38} {_fmt_int(relabeled):>6}  "
        f"({relabeled / total * 100:.1f} %)",
        f"  {'Kept original header (lookup failed)':<38} {_fmt_int(failed):>6}  "
        f"({failed / total * 100:.1f} %)",
        f"  {'Unique accessions queried':<38} {_fmt_int(unique_acc):>6}",
        f"  {'Cache hits':<38} {_fmt_int(total - unique_acc):>6}",
        f"  {'With nucleotide / genome ID':<38} {_fmt_int(with_nucl):>6}  "
        f"({with_nucl / total * 100:.1f} %)"
        if total else "",
        f"  {'Without nucleotide / genome ID':<38} {_fmt_int(relabeled - with_nucl):>6}  "
        f"({(relabeled - with_nucl) / total * 100:.1f} %)"
        if total else "",
    ]

    # ── Sequence length statistics ──────────────────────────────────────────
    lines += [_section("SEQUENCE LENGTH STATISTICS")]
    lmin = min(lengths)
    lmax = max(lengths)
    lmean = statistics.mean(lengths)
    lmedian = statistics.median(lengths)
    lstdev = statistics.stdev(lengths) if len(lengths) > 1 else 0.0
    ln50 = _n50(lengths)
    ltotal = sum(lengths)
    lines += [
        f"  {'Minimum':<22} {_fmt_int(lmin):>8} residues",
        f"  {'Maximum':<22} {_fmt_int(lmax):>8} residues",
        f"  {'Mean':<22} {_fmt_float(lmean):>8} residues",
        f"  {'Median':<22} {_fmt_float(lmedian):>8} residues",
        f"  {'Std deviation':<22} {_fmt_float(lstdev):>8} residues",
        f"  {'N50':<22} {_fmt_int(ln50):>8} residues",
        f"  {'Total residues':<22} {_fmt_int(ltotal):>8}",
    ]

    hist = _histogram(lengths, n_bins=10)
    if hist:
        hi_count = max(cnt for _, _, cnt in hist)
        w_lo = len(str(max(lo for lo, _, _ in hist)))
        w_hi = len(str(max(hi for _, hi, _ in hist)))
        lines.append("")
        lines.append("  Length histogram:")
        for lo, hi, cnt in hist:
            bar = _bar(cnt / hi_count, width=24) if hi_count else ""
            pct = cnt / total * 100
            lines.append(
                f"  [{lo:{w_lo},} – {hi:{w_hi},}]  {bar:<24}  "
                f"{_fmt_int(cnt):>5}  ({pct:5.1f} %)"
            )

    # ── Organism distribution ────────────────────────────────────────────────
    org_counts = Counter(r.organism or "(unknown)" for r in results if r.relabeled)
    _render_distribution(lines, "ORGANISM DISTRIBUTION", org_counts, total, top_n)

    # ── Protein distribution ─────────────────────────────────────────────────
    prot_counts = Counter(r.protein_name or "(unknown)" for r in results if r.relabeled)
    _render_distribution(lines, "PROTEIN DISTRIBUTION", prot_counts, total, top_n)

    # ── Nucleotide ID coverage per organism ──────────────────────────────────
    if with_nucl:
        lines += [_section("NUCLEOTIDE ID COVERAGE  (relabeled sequences only)")]
        # For each organism: how many have a nucl_id vs not
        org_nucl: dict[str, list[bool]] = {}
        for r in results:
            if r.relabeled:
                key = r.organism or "(unknown)"
                org_nucl.setdefault(key, []).append(bool(r.nucl_id))
        # Sort by organism name
        lines.append(
            f"  {'Organism':<38} {'w/ nucl_id':>10} {'total':>7} {'coverage':>9}"
        )
        lines.append("  " + "─" * 66)
        for org in sorted(org_nucl):
            flags = org_nucl[org]
            w = sum(flags)
            t = len(flags)
            cov = w / t * 100
            bar = _bar(w / t, width=10)
            truncated = org[:37] + "…" if len(org) > 38 else org
            lines.append(
                f"  {truncated:<38} {_fmt_int(w):>10} {_fmt_int(t):>7}  "
                f"{bar:<10}  {cov:5.1f} %"
            )

    # ── Per-file summary ─────────────────────────────────────────────────────
    lines += [_section("PER-FILE SUMMARY")]
    by_file: dict[str, list[RecordResult]] = {}
    for r in results:
        by_file.setdefault(r.source_file, []).append(r)

    col_w = max(len(f) for f in by_file) + 2
    col_w = max(col_w, 24)
    lines.append(
        f"  {'File':<{col_w}} {'Records':>8} {'Relabeled':>10} "
        f"{'Failed':>8} {'w/NuclID':>9} {'Min':>7} {'Max':>7} {'Mean':>7}"
    )
    lines.append("  " + "─" * (col_w + 62))
    for fname, recs in sorted(by_file.items()):
        n = len(recs)
        rel = sum(1 for r in recs if r.relabeled)
        fail = n - rel
        wn = sum(1 for r in recs if r.nucl_id)
        lens = [r.seq_len for r in recs]
        lines.append(
            f"  {fname:<{col_w}} {_fmt_int(n):>8} {_fmt_int(rel):>10} "
            f"{_fmt_int(fail):>8} {_fmt_int(wn):>9} "
            f"{_fmt_int(min(lens)):>7} {_fmt_int(max(lens)):>7} "
            f"{_fmt_float(statistics.mean(lens)):>7}"
        )

    # ── Failed accessions ─────────────────────────────────────────────────────
    failed_recs = [r for r in results if not r.relabeled]
    if failed_recs:
        lines += [_section(
            f"FAILED ACCESSIONS  ({len(failed_recs)} records kept with original header)"
        )]
        for r in failed_recs:
            err = f"  {r.error}" if r.error else ""
            lines.append(f"  {r.accession:<30}  {r.source_file}{err}")

    lines += ["", RULE_HEAVY, ""]
    return "\n".join(lines)


def _render_distribution(
        lines: list[str],
        title: str,
        counts: Counter,
        total: int,
        top_n: int,
) -> None:
    unique = len(counts)
    shown = min(top_n, unique)
    suffix = f"top {shown} of {_fmt_int(unique)}" if unique > top_n else f"{_fmt_int(unique)}"
    lines += [_section(f"{title}  ({suffix} unique)")]
    if not counts:
        lines.append("  (none)")
        return
    top = counts.most_common(top_n)
    max_count = top[0][1]
    label_w = min(40, max(len(label) for label, _ in top) + 1)
    for label, cnt in top:
        bar = _bar(cnt / max_count, width=22) if max_count else ""
        pct = cnt / total * 100
        truncated = label[:label_w - 1] + "…" if len(label) > label_w else label
        lines.append(
            f"  {truncated:<{label_w}}  {bar:<22}  {_fmt_int(cnt):>6}  ({pct:5.1f} %)"
        )


# ---------------------------------------------------------------------------
# TSV report (one row per record)
# ---------------------------------------------------------------------------

_TSV_FIELDS = [
    "source_file", "accession", "db_used", "relabeled",
    "seq_len", "organism", "taxid", "protein_name", "gene",
    "nucl_id", "nucl_ids", "lineage", "new_header", "error",
]


def build_report_tsv(results: list[RecordResult]) -> str:
    rows = ["\t".join(_TSV_FIELDS)]
    for r in results:
        row: list[str] = []
        for f in _TSV_FIELDS:
            val = getattr(r, f, "")
            if isinstance(val, list):
                val = ";".join(val)
            row.append(str(val))
        rows.append("\t".join(row))
    return "\n".join(rows) + "\n"


# ---------------------------------------------------------------------------
# JSON report
# ---------------------------------------------------------------------------

def build_report_json(
        results: list[RecordResult],
        fmt: str,
        db: str,
        nucl_type: str,
        input_files: list[Path],
        elapsed: float,
) -> str:
    total = len(results)
    relabeled = sum(1 for r in results if r.relabeled)
    lengths = [r.seq_len for r in results]
    with_nucl = sum(1 for r in results if r.relabeled and r.nucl_id)
    org_counts = Counter(r.organism or "" for r in results if r.relabeled)
    prot_counts = Counter(r.protein_name or "" for r in results if r.relabeled)

    doc = {
        "meta": {
            "version": VERSION,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "format_string": fmt,
            "database": db,
            "nucl_type": nucl_type,
            "files_processed": len(input_files),
            "elapsed_seconds": round(elapsed, 2),
        },
        "coverage": {
            "total_sequences": total,
            "relabeled": relabeled,
            "failed": total - relabeled,
            "unique_accessions": len({r.accession for r in results}),
            "with_nucl_id": with_nucl,
            "without_nucl_id": relabeled - with_nucl,
        },
        "length_stats": {
            "min": min(lengths) if lengths else 0,
            "max": max(lengths) if lengths else 0,
            "mean": round(statistics.mean(lengths), 2) if lengths else 0,
            "median": float(statistics.median(lengths)) if lengths else 0,
            "stdev": round(statistics.stdev(lengths), 2) if len(lengths) > 1 else 0,
            "n50": _n50(lengths),
            "total_residues": sum(lengths),
        },
        "organism_distribution": dict(org_counts.most_common()),
        "protein_distribution": dict(prot_counts.most_common()),
        "per_file": {
            fname: {
                "records": len(recs),
                "relabeled": sum(1 for r in recs if r.relabeled),
                "failed": sum(1 for r in recs if not r.relabeled),
                "with_nucl_id": sum(1 for r in recs if r.nucl_id),
                "length_min": min(r.seq_len for r in recs),
                "length_max": max(r.seq_len for r in recs),
                "length_mean": round(statistics.mean(r.seq_len for r in recs), 1),
            }
            for fname, recs in _group_by_file(results).items()
        },
        "records": [
            {**dataclasses.asdict(r), "nucl_ids": r.nucl_ids}
            for r in results
        ],
    }
    return json.dumps(doc, indent=2, ensure_ascii=False)


def _group_by_file(results: list[RecordResult]) -> dict[str, list[RecordResult]]:
    out: dict[str, list[RecordResult]] = {}
    for r in results:
        out.setdefault(r.source_file, []).append(r)
    return out


# ---------------------------------------------------------------------------
# Report dispatch helper + file-safety helper
# ---------------------------------------------------------------------------

# Maps --report-format value to the file extension used for per-file reports
_REPORT_EXT: dict[str, str] = {"text": ".txt", "tsv": ".tsv", "json": ".json"}


def _render_report(
        results: list[RecordResult],
        fmt: str,
        db: str,
        nucl_type: str,
        top_n: int,
        input_files: list[Path],
        elapsed: float,
        report_format: str,
) -> str:
    """Dispatch to the correct report builder and return the rendered string."""
    if report_format == "tsv":
        return build_report_tsv(results)
    if report_format == "json":
        return build_report_json(results, fmt, db, nucl_type, input_files, elapsed)
    return build_report_text(results, fmt, db, nucl_type, top_n, input_files, elapsed)


def _safe_write(path: Path, content: str, overwrite: bool) -> bool:
    """Write *content* to *path*; return True on success, False if skipped."""
    if path.exists() and not overwrite:
        print(
            f"  Skipping {path}: file already exists (use --overwrite to replace)",
            file=sys.stderr,
        )
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return True


# ---------------------------------------------------------------------------
# File processing
# ---------------------------------------------------------------------------

def process_file(
        input_path: Path,
        output_path: Path,
        fmt: str,
        db: str,
        email: Optional[str],
        api_key: Optional[str],
        nucl_type: str,
        use_underscores: bool,
        overwrite: bool,
        verbose: bool,
        dry_run: bool,
) -> list[RecordResult]:
    records = list(parse_fasta(input_path))
    if not records:
        print(f"Warning: no records found in {input_path}", file=sys.stderr)
        return []

    if verbose:
        print(f"\n{input_path}  ({len(records)} records)", file=sys.stderr)

    results: list[RecordResult] = []
    new_records: list[tuple[str, str]] = []

    for header, seq in records:
        detected_db, accession = detect_id(header, db)
        info, error = fetch_info(accession, detected_db, email, api_key, nucl_type, verbose)

        if info is None:
            res = RecordResult(
                source_file=input_path.name,
                original_header=header,
                accession=accession,
                db_used=detected_db,
                seq_len=len(seq),
                relabeled=False,
                new_header=header,
                error=error,
            )
            new_records.append((header, seq))
        else:
            info["id"] = accession
            new_header = format_header(fmt, info, use_underscores)
            if verbose:
                nucl_note = f"  nucl_id={info['nucl_id']!r}" if info.get("nucl_id") else ""
                print(f"  {accession!r} → {new_header!r}{nucl_note}", file=sys.stderr)
            res = RecordResult(
                source_file=input_path.name,
                original_header=header,
                accession=accession,
                db_used=detected_db,
                seq_len=len(seq),
                relabeled=True,
                new_header=new_header,
                organism=info.get("organism", ""),
                protein_name=info.get("protein_name", ""),
                gene=info.get("gene", ""),
                taxid=info.get("taxid", ""),
                lineage=info.get("lineage", ""),
                nucl_id=info.get("nucl_id", ""),
                nucl_ids=info.get("nucl_ids", []),
            )
            new_records.append((new_header, seq))

        results.append(res)

    if dry_run:
        for h, _ in new_records:
            print(f">{h}")
    else:
        if write_fasta(output_path, new_records, overwrite=overwrite):
            print(f"Wrote {output_path}", file=sys.stderr)

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fasta_relabel",
        description=f"fasta_relabel {VERSION} — Relabel FASTA headers with taxonomy/protein info from UniProt or NCBI.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
FORMAT FIELDS  (use any alias inside curly braces)
  {id}             Accession / sequence identifier
  {name}           Protein name  [aliases: protein_name, protein]
  {organism}       Scientific organism name  [aliases: taxon, species, sci_name]
  {gene}           Gene name
  {lineage}        Full taxonomic lineage (semicolon-separated)
  {taxid}          NCBI Taxonomy ID  [alias: tax_id]
  {entry_name}     UniProt entry name or GenBank accession+version  [alias: entry]
  {nucl_id}        Nucleotide / genome accession linked to this protein
                   [aliases: genome_id, nucleotide_id, nucl_acc, coding_seq]
                   Source: EMBL cross-ref (UniProt) or CDS coded_by (NCBI)

NUCLEOTIDE ID TYPE  (--nucl-type)
  genomic  Prefer Genomic_DNA accessions, fall back to mRNA/other  (default)
  mrna     Prefer mRNA / cDNA accessions, fall back to genomic
  any      Use the first cross-reference listed regardless of molecule type
  Note: for NCBI records the CDS coded_by qualifier is used directly and is
  unaffected by --nucl-type.

REPORT FORMATS
  text  Human-readable tables with ASCII bar charts (default → stderr)
  tsv   One row per sequence — import directly into Excel / R / pandas
  json  Fully structured JSON with per-record data and aggregate statistics

EXAMPLES
  fasta_relabel sequences.fasta
  fasta_relabel sequences.fasta --format "{organism} | {id} | {name}"
  fasta_relabel sequences.fasta --format "{organism} | {id} | {name} | {nucl_id}"
  fasta_relabel sequences.fasta --format "{id}|{nucl_id}|{name}" --nucl-type mrna
  fasta_relabel sequences.fasta --format "{id} {name}" --db ncbi --email me@example.com
  fasta_relabel seqs/ --format "{lineage}; {organism}" --output renamed/
  fasta_relabel a.fasta b.fasta --in-place --format "{id}|{taxon}|{name}|{genome_id}"
  fasta_relabel sequences.fasta --dry-run --format "{organism} {id} {name}"
  fasta_relabel sequences.fasta --report summary.txt
  fasta_relabel sequences.fasta --report data.tsv --report-format tsv
  fasta_relabel sequences.fasta --report data.json --report-format json
""",
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {VERSION}",
    )

    # ── Input / output ───────────────────────────────────────────────────────
    parser.add_argument(
        "input", nargs="+",
        help="FASTA file(s) or directory/directories to process",
    )
    parser.add_argument(
        "--format", "-f",
        default="{id} {name} [{organism}]",
        dest="fmt",
        metavar="FMT",
        help="Header format string (default: '{id} {name} [{organism}]')",
    )
    parser.add_argument(
        "--db",
        choices=["auto", "uniprot", "ncbi"],
        default="auto",
        help="Database to query — auto-detects from accession format (default: auto)",
    )
    parser.add_argument(
        "--output", "-o",
        metavar="DIR",
        help="Write output files into DIR (default: same dir as input + --suffix)",
    )
    parser.add_argument(
        "--suffix",
        default="_relabeled",
        help="Filename suffix when --output is not given (default: _relabeled)",
    )
    parser.add_argument(
        "--in-place", action="store_true",
        help="Relabel input files in place (mutually exclusive with --output)",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Allow overwriting existing output files and report files "
             "(default: skip and warn if a destination already exists)",
    )
    parser.add_argument(
        "--ext",
        default=".fasta,.fa,.fna,.faa,.ffn",
        metavar="EXTS",
        help="Comma-separated extensions for directory scans "
             "(default: .fasta,.fa,.fna,.faa,.ffn)",
    )

    # ── Database options ─────────────────────────────────────────────────────
    parser.add_argument(
        "--email",
        metavar="ADDR",
        help="E-mail for NCBI E-utilities (strongly recommended for batch jobs)",
    )
    parser.add_argument(
        "--ncbi-api-key",
        metavar="KEY",
        help="NCBI API key — raises the rate limit from 3 to 10 requests/second",
    )

    # ── Nucleotide ID options ─────────────────────────────────────────────────
    parser.add_argument(
        "--nucl-type",
        choices=["genomic", "mrna", "any"],
        default="genomic",
        metavar="TYPE",
        help="Which nucleotide cross-reference to use when multiple exist "
             "(genomic | mrna | any, default: genomic). "
             "Applies to UniProt EMBL cross-refs only; "
             "NCBI uses the CDS coded_by qualifier directly.",
    )

    # ── Report options ───────────────────────────────────────────────────────
    parser.add_argument(
        "--report", "-r",
        metavar="FILE",
        help="Write the run report to FILE instead of stderr (use '-' for stdout)",
    )
    parser.add_argument(
        "--report-format",
        choices=["text", "tsv", "json"],
        default="text",
        help="Report format: text (default), tsv, or json",
    )
    parser.add_argument(
        "--report-top",
        type=int,
        default=20,
        metavar="N",
        help="Show top N organisms/proteins in the text report (default: 20)",
    )
    parser.add_argument(
        "--no-report", action="store_true",
        help="Suppress both the aggregate run report and all per-file reports",
    )
    parser.add_argument(
        "--no-per-file-report", action="store_true",
        help="Do not write individual report files alongside each output FASTA "
             "(per-file reports are written by default)",
    )

    # ── Misc ─────────────────────────────────────────────────────────────────
    parser.add_argument(
        "--spaces", action="store_true",
        help="Keep spaces in output headers (default: replace all whitespace "
             "with underscores so each header is a single token)",
    )
    parser.add_argument(
        "--dry-run", "-n", action="store_true",
        help="Print renamed headers to stdout; do not write any files",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Show per-record progress on stderr",
    )

    return parser


def resolve_output(
        input_path: Path,
        output_dir: Optional[str],
        suffix: str,
        in_place: bool,
) -> Path:
    if in_place:
        return input_path
    if output_dir:
        return Path(output_dir) / input_path.name
    return input_path.with_stem(input_path.stem + suffix)


def main(argv: Optional[list[str]] = None) -> None:
    # Accept single-dash convenience aliases for --help and --version
    if argv is None:
        argv = sys.argv[1:]
    argv = [
        "--help" if a == "-help" else
        "--version" if a == "-version" else a
        for a in argv
    ]

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.in_place and args.output:
        parser.error("--in-place and --output are mutually exclusive")

    if args.ncbi_api_key:
        global NCBI_INTERVAL
        NCBI_INTERVAL = 0.1  # ≤10 req/s with an API key

    extensions = {e.strip() for e in args.ext.split(",")}

    # Collect input files
    input_files: list[Path] = []
    for raw in args.input:
        p = Path(raw)
        if p.is_dir():
            for ext in extensions:
                input_files.extend(sorted(p.rglob(f"*{ext}")))
        elif p.is_file():
            input_files.append(p)
        else:
            print(f"Error: {raw!r} is not a file or directory", file=sys.stderr)
            sys.exit(1)

    if not input_files:
        print("Error: no input files found", file=sys.stderr)
        sys.exit(1)

    # --in-place always implies overwrite (the user explicitly asked for it)
    overwrite = args.overwrite or args.in_place

    t_start = time.monotonic()
    all_results: list[RecordResult] = []
    skipped_fasta = 0

    for input_path in input_files:
        output_path = resolve_output(input_path, args.output, args.suffix, args.in_place)

        # ── Overwrite guard for output FASTA ────────────────────────────────
        if not args.dry_run and output_path.exists() and not overwrite:
            print(
                f"Skipping {input_path.name}: {output_path} already exists "
                f"(use --overwrite to replace)",
                file=sys.stderr,
            )
            skipped_fasta += 1
            continue

        file_results = process_file(
            input_path=input_path,
            output_path=output_path,
            fmt=args.fmt,
            db=args.db,
            email=args.email,
            api_key=args.ncbi_api_key,
            nucl_type=args.nucl_type,
            use_underscores=not args.spaces,
            overwrite=overwrite,
            verbose=args.verbose,
            dry_run=args.dry_run,
        )
        all_results.extend(file_results)

        # ── Per-file report ─────────────────────────────────────────────────
        if (
                file_results
                and not args.no_report
                and not args.no_per_file_report
                and not args.dry_run
        ):
            t_file = time.monotonic() - t_start
            per_report = _render_report(
                file_results, args.fmt, args.db, args.nucl_type,
                args.report_top, [input_path], t_file, args.report_format,
            )
            ext = _REPORT_EXT[args.report_format]
            per_report_path = output_path.with_name(output_path.stem + "_report" + ext)
            if _safe_write(per_report_path, per_report, overwrite):
                print(f"Per-file report → {per_report_path}", file=sys.stderr)

    elapsed = time.monotonic() - t_start

    if skipped_fasta:
        print(
            f"\nWarning: {skipped_fasta} file(s) skipped — output already exists.",
            file=sys.stderr,
        )

    # ── Aggregate report ─────────────────────────────────────────────────────
    if not args.no_report and all_results:
        report_text = _render_report(
            all_results, args.fmt, args.db, args.nucl_type,
            args.report_top, input_files, elapsed, args.report_format,
        )

        if args.report == "-":
            print(report_text)
        elif args.report:
            if _safe_write(Path(args.report), report_text, overwrite):
                print(f"Aggregate report → {args.report}", file=sys.stderr)
        else:
            print(report_text, file=sys.stderr)


if __name__ == "__main__":
    main()
