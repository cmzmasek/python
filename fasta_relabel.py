#!/usr/bin/env python3
"""
fasta_relabel.py — Relabel FASTA description lines via UniProt / NCBI lookups.

Each sequence's accession is detected automatically (or forced with --db), then
the description line is rebuilt from a user-supplied format string such as:

    "{organism} | {id} | {name} | {nucl_id}"
    "{lineage}; {organism}; {id}"
    "{organism}|{id}|{name}|{ec}|{go}"

A detailed report on labeling coverage, sequence-length statistics, and the
distribution of organisms / proteins is always generated after processing.
"""

import argparse
import dataclasses
import json
import math
import random
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


VERSION = "1.0.1"


# ---------------------------------------------------------------------------
# Rate-limiting helpers
# ---------------------------------------------------------------------------

NCBI_INTERVAL     = 0.34    # ≤3 req/s without an API key
UNIPROT_INTERVAL  = 0.10
ENSEMBL_INTERVAL  = 0.15    # conservative; Ensembl REST allows ~6-7 req/s unauthenticated
INTERPRO_INTERVAL = 0.20    # EBI servers: ~5 req/s to be safe
PDB_INTERVAL      = 0.10    # RCSB is generous; 10 req/s is well within limits

_last_ncbi_request:     float = 0.0
_last_uniprot_request:  float = 0.0
_last_ensembl_request:  float = 0.0
_last_interpro_request: float = 0.0
_last_pdb_request:      float = 0.0

# Default HTTP retry settings (overridden by --retries)
HTTP_RETRIES: int = 3
HTTP_BACKOFF: float = 1.0   # initial sleep seconds; doubles on each retry


def _throttle(interval: float, last_t: float) -> float:
    """Sleep if the required interval has not yet elapsed; return the new timestamp."""
    wait = interval - (time.monotonic() - last_t)
    if wait > 0:
        time.sleep(wait)
    return time.monotonic()


def _get(url: str, timeout: int = 30) -> str:
    """Fetch *url* with automatic retry on transient errors.

    Retries up to ``HTTP_RETRIES`` times with exponential back-off starting at
    ``HTTP_BACKOFF`` seconds.  HTTP 429 / 5xx and network errors are retried;
    HTTP 4xx (except 429) are raised immediately as non-retryable.
    """
    headers = {"User-Agent": f"fasta_relabel/{VERSION} (bioinformatics tool)"}
    last_exc: Exception = RuntimeError("no attempts made")
    for attempt in range(max(HTTP_RETRIES, 1)):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as exc:
            if exc.code in (429, 500, 502, 503, 504) and attempt < HTTP_RETRIES - 1:
                wait = HTTP_BACKOFF * (2 ** attempt)
                # Honour Retry-After header when the server supplies one
                retry_after = exc.headers.get("Retry-After", "")
                try:
                    wait = max(wait, float(retry_after))
                except (ValueError, TypeError):
                    pass
                time.sleep(wait)
                last_exc = exc
                continue
            raise
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            if attempt < HTTP_RETRIES - 1:
                time.sleep(HTTP_BACKOFF * (2 ** attempt))
                last_exc = exc
                continue
            raise
    raise last_exc


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


def _check_overwrite(path: Path, overwrite: bool) -> bool:
    """Return True if the caller should proceed; print a warning and return False otherwise.

    The guard only applies to regular files — devices such as ``/dev/null`` are
    always treated as writable regardless of the *overwrite* flag.
    """
    if path.is_file() and not overwrite:
        print(
            f"  Skipping {path}: file already exists (use --overwrite to replace)",
            file=sys.stderr,
        )
        return False
    return True


def _write_fasta_record(fh, header: str, seq: str, wrap: int = 60) -> None:
    """Write a single FASTA record to an already-open file handle."""
    fh.write(f">{header}\n")
    for i in range(0, len(seq), wrap):
        fh.write(seq[i : i + wrap] + "\n")


def write_fasta(
    path: Path, records: list[tuple[str, str]], wrap: int = 60, overwrite: bool = True
) -> bool:
    """Write FASTA records to *path*; return False (and skip) if file exists and overwrite=False."""
    if not _check_overwrite(path, overwrite):
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for header, seq in records:
            _write_fasta_record(fh, header, seq, wrap)
    return True


# ---------------------------------------------------------------------------
# Sequence quality helpers
# ---------------------------------------------------------------------------

# Ambiguity character sets for amino-acid and nucleotide sequences
_AMBIGUOUS_AA   = frozenset("BXZJbxzj")
_AMBIGUOUS_NUCL = frozenset("NRYSWKMBDHVnryswkmbdhv")
_CORE_NUCL      = frozenset("ACGTUacgtu")


def _detect_molecule_type(seq: str) -> str:
    """Return ``'nucl'`` if ≥85 % of non-gap characters are ACGTU, else ``'prot'``."""
    if not seq:
        return "prot"
    core = [c for c in seq if c not in "-. \t"]
    if not core:
        return "prot"
    return "nucl" if sum(1 for c in core if c in _CORE_NUCL) / len(core) >= 0.85 else "prot"


def _ambiguous_ratio(seq: str, amb_chars: frozenset) -> float:
    """Return the fraction of characters in *seq* that appear in *amb_chars*."""
    if not seq:
        return 0.0
    return sum(1 for c in seq if c in amb_chars) / len(seq)


def apply_sequence_filters(
    records: list[tuple[str, str]],
    source_file: str,
    min_len: Optional[int],
    max_len: Optional[int],
    max_ambiguous_ratio: Optional[float],
    dedupe_sequence: bool,
    sample: Optional[int],
    molecule_type: str,          # 'auto' | 'prot' | 'nucl'
) -> tuple[list[tuple[str, str]], list["RecordResult"]]:
    """Apply pre-lookup filters to *records*.

    Order: length → ambiguity → deduplication → random sample.
    Returns ``(kept, dropped)`` where *dropped* entries have ``filtered=True``.
    """
    # Determine ambiguity character set once for the whole batch
    if max_ambiguous_ratio is not None:
        mol = molecule_type
        if mol == "auto" and records:
            mol = _detect_molecule_type(records[0][1])
        amb_chars = _AMBIGUOUS_NUCL if mol == "nucl" else _AMBIGUOUS_AA
    else:
        amb_chars = _AMBIGUOUS_AA  # unused, but keeps the variable defined

    kept: list[tuple[str, str]] = []
    dropped: list[RecordResult] = []
    seen_seqs: set[str] = set()

    for header, seq in records:
        seq_upper = seq.upper()
        reason = ""

        if min_len is not None and len(seq) < min_len:
            reason = f"too short ({len(seq)} < {min_len})"
        elif max_len is not None and len(seq) > max_len:
            reason = f"too long ({len(seq)} > {max_len})"
        elif max_ambiguous_ratio is not None and _ambiguous_ratio(seq, amb_chars) > max_ambiguous_ratio:
            ratio = _ambiguous_ratio(seq, amb_chars)
            reason = f"ambiguous ratio {ratio:.3f} > {max_ambiguous_ratio:.3f}"
        elif dedupe_sequence and seq_upper in seen_seqs:
            reason = "duplicate sequence"

        if reason:
            dropped.append(_filtered_result(source_file, header, seq, reason))
        else:
            seen_seqs.add(seq_upper)
            kept.append((header, seq))

    # Random subsample (applied last so N refers to post-filter count)
    if sample is not None and sample < len(kept):
        chosen_idx = set(random.sample(range(len(kept)), sample))
        new_kept: list[tuple[str, str]] = []
        for i, (header, seq) in enumerate(kept):
            if i in chosen_idx:
                new_kept.append((header, seq))
            else:
                dropped.append(_filtered_result(source_file, header, seq, "not sampled"))
        kept = new_kept

    return kept, dropped


def _filtered_result(
    source_file: str, header: str, seq: str, reason: str
) -> "RecordResult":
    """Build a minimal ``RecordResult`` for a record dropped by a pre-lookup filter."""
    return RecordResult(
        source_file=source_file,
        original_header=header,
        accession=header.split()[0] if header else "",
        db_used="",
        seq_len=len(seq),
        relabeled=False,
        new_header=header,
        filtered=True,
        filter_reason=reason,
    )


def _compile_patterns(patterns: list[str]) -> list[re.Pattern]:
    """Compile a list of user-supplied strings to case-insensitive regex patterns.

    Invalid regex patterns are automatically escaped and treated as literal
    substrings so the user never gets a hard crash from a typo.
    """
    compiled: list[re.Pattern] = []
    for pat in patterns:
        try:
            compiled.append(re.compile(pat, re.IGNORECASE))
        except re.error:
            compiled.append(re.compile(re.escape(pat), re.IGNORECASE))
    return compiled


def _post_lookup_filter_reason(
    res: "RecordResult",
    org_include: list[re.Pattern],
    org_exclude: list[re.Pattern],
    taxid_include: set[str],
    taxid_exclude: set[str],
    lineage_include: list[re.Pattern],
    lineage_exclude: list[re.Pattern],
) -> str:
    """Return a non-empty reason string if *res* should be dropped, else ``''``.

    Records with no organism/taxid/lineage (lookup failed) are never filtered.
    """
    org     = res.organism or ""
    taxid   = res.taxid    or ""
    lineage = res.lineage  or ""

    if taxid and taxid_exclude and taxid in taxid_exclude:
        return f"taxid excluded: {taxid}"
    if taxid and taxid_include and taxid not in taxid_include:
        return f"taxid not included: {taxid}"
    if org and org_exclude and any(p.search(org) for p in org_exclude):
        return f"organism excluded: {org}"
    if org and org_include and not any(p.search(org) for p in org_include):
        return f"organism not included: {org}"
    if lineage and lineage_exclude and any(p.search(lineage) for p in lineage_exclude):
        return f"lineage excluded: {lineage}"
    if lineage and lineage_include and not any(p.search(lineage) for p in lineage_include):
        return f"lineage not included: {lineage}"
    return ""


def apply_post_lookup_filters(
    results: list["RecordResult"],
    new_records: list[tuple[str, str]],
    org_include: list[re.Pattern],
    org_exclude: list[re.Pattern],
    taxid_include: set[str],
    taxid_exclude: set[str],
    lineage_include: list[re.Pattern],
    lineage_exclude: list[re.Pattern],
) -> tuple[list["RecordResult"], list[tuple[str, str]], list["RecordResult"]]:
    """Apply all post-lookup filters in one pass.

    Returns ``(kept_results, kept_records, dropped_results)``.
    """
    any_filter = (
        org_include or org_exclude
        or taxid_include or taxid_exclude
        or lineage_include or lineage_exclude
    )
    if not any_filter:
        return results, new_records, []

    kept_results:   list[RecordResult]       = []
    kept_records:   list[tuple[str, str]]    = []
    dropped_results: list[RecordResult]      = []

    for res, rec in zip(results, new_records):
        reason = _post_lookup_filter_reason(
            res, org_include, org_exclude,
            taxid_include, taxid_exclude,
            lineage_include, lineage_exclude,
        )
        if reason:
            dropped_results.append(dataclasses.replace(res, filtered=True, filter_reason=reason))
        else:
            kept_results.append(res)
            kept_records.append(rec)

    return kept_results, kept_records, dropped_results


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
# UniProt pipe-format embedded anywhere in a string: sp|ACC|… or tr|ACC|…
_UNIPROT_PIPE_RE = re.compile(r"(?:sp|tr)\|([A-Z0-9]+)\|")
# Ensembl: gene (G), transcript (T), protein (P), exon (E) — optional species infix
_ENSEMBL_RE = re.compile(r"\b(ENS[A-Z]{0,4}[GTPE]\d{11}(?:\.\d+)?)\b")
# InterPro domain/family entry
_INTERPRO_RE = re.compile(r"\b(IPR\d{6})\b")
# PDB: digit + 3 alphanumeric chars, optional underscore + entity/chain number
_PDB_RE = re.compile(r"\b([0-9][A-Za-z0-9]{3}(?:_[0-9]+)?)\b")


def _match_token(token: str, db_hint: str) -> Optional[tuple[str, str]]:
    """Try to extract a recognised accession from a single field string.

    Returns ``(db, accession)`` on success, ``None`` if the token contains no
    known accession pattern.  Priority order (most-specific first):

    UniProt pipe → InterPro → Ensembl → UniProt bare → gi| → RefSeq →
    GenBank → PDB
    """
    token = token.strip()
    if not token:
        return None

    if db_hint in ("auto", "uniprot"):
        m = _UNIPROT_PIPE_RE.search(token)
        if m:
            return "uniprot", m.group(1)

    if db_hint in ("auto", "interpro"):
        m = _INTERPRO_RE.search(token)
        if m:
            return "interpro", m.group(1)

    if db_hint in ("auto", "ensembl"):
        m = _ENSEMBL_RE.search(token)
        if m:
            return "ensembl", m.group(1)

    if db_hint in ("auto", "uniprot"):
        m = _UNIPROT_RE.search(token)
        if m:
            return "uniprot", m.group(1)

    if db_hint in ("auto", "ncbi"):
        m = _GI_RE.search(token)
        if m:
            return "ncbi", m.group(1)
        m = _REFSEQ_RE.search(token)
        if m:
            return "ncbi", m.group(1)
        m = _GENBANK_RE.search(token)
        if m:
            return "ncbi", m.group(1)

    # PDB last in auto: digit-first pattern is distinctive but broad
    if db_hint in ("auto", "pdb"):
        m = _PDB_RE.match(token)
        if m:
            return "pdb", m.group(1).upper()

    return None


def detect_id(
    header: str,
    db_hint: str,
    id_delimiter: Optional[str] = None,
    id_field: Optional[int] = None,
) -> tuple[str, str]:
    """Return ``(db, accession)`` extracted from a FASTA header line.

    Parameters
    ----------
    header:
        The raw header text (without the leading ``>``).
    db_hint:
        One of ``"auto"``, ``"uniprot"``, or ``"ncbi"``.
    id_delimiter:
        Character(s) used to split *header* into fields before searching.
        ``None`` means split on any whitespace (the default).
    id_field:
        0-based index of the field that holds the accession.  When given, only
        that field is inspected; when omitted every field is tried in order.
    """
    # Resolve the fallback db once rather than duplicating the ternary below
    fallback_db = db_hint if db_hint != "auto" else "ncbi"

    # Split the header into candidate fields
    tokens: list[str] = (
        header.split(id_delimiter) if id_delimiter is not None else header.split()
    )

    # ── Pinned field: the user told us exactly where the ID lives ────────────
    if id_field is not None:
        try:
            candidate = tokens[id_field].strip()
        except IndexError:
            candidate = tokens[0].strip() if tokens else header
        hit = _match_token(candidate, db_hint)
        if hit:
            return hit
        # Use the field verbatim even if it doesn't look like a known accession
        return fallback_db, candidate

    # ── Scan every field in order ────────────────────────────────────────────
    for token in tokens:
        hit = _match_token(token, db_hint)
        if hit:
            return hit

    # ── Fallback: first field, pipe-stripped ─────────────────────────────────
    first = (tokens[0].strip() if tokens else header).split("|")[-1]
    return fallback_db, first


# ---------------------------------------------------------------------------
# Nucleotide ID helpers
# ---------------------------------------------------------------------------

# Molecule-type preference orders for --nucl-type
_NUCL_TYPE_PREF: dict[str, list[str]] = {
    "genomic": ["Genomic_DNA", "Genomic_RNA", "mRNA", "Transcribed_RNA", "Other_RNA"],
    "mrna":    ["mRNA", "Transcribed_RNA", "Other_RNA", "Genomic_DNA", "Genomic_RNA"],
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



# ---------------------------------------------------------------------------
# Shared info-dict factory
# ---------------------------------------------------------------------------

def _empty_info(accession: str, entry_name: str = "") -> dict:
    """Return a zeroed info dict with the canonical set of keys."""
    return {
        "id": accession,
        "protein_name": "",
        "gene": "",
        "organism": "",
        "taxid": "",
        "lineage": "",
        "entry_name": entry_name or accession,
        "nucl_id": "",
        "nucl_ids": [],
        "go_terms": [],    # list of "GO:XXXXXXX name" strings
        "ec_numbers": [],  # list of EC number strings, e.g. ["3.4.21.4"]
    }


# ---------------------------------------------------------------------------
# UniProt query
# ---------------------------------------------------------------------------

def _query_uniprot(accession: str, nucl_type: str) -> dict:
    global _last_uniprot_request
    _last_uniprot_request = _throttle(UNIPROT_INTERVAL, _last_uniprot_request)
    data = json.loads(_get(f"https://rest.uniprot.org/uniprotkb/{accession}.json"))

    result = _empty_info(accession, data.get("uniProtkbId", accession))

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

    # EC numbers from recommended / alternative / submission names
    ec_set: list[str] = []
    for name_block in (
        [pd_.get("recommendedName") or {}]
        + (pd_.get("alternativeNames") or [])
        + (pd_.get("submissionNames") or [])
    ):
        for ec in name_block.get("ecNumbers") or []:
            v = ec.get("value", "")
            if v and v not in ec_set:
                ec_set.append(v)
    result["ec_numbers"] = ec_set

    # GO terms from cross-references — stored as "X:GO:XXXXXXX term_name"
    # where X is F (molecular function), P (biological process), C (cellular component)
    go_terms: list[str] = []
    for xref in data.get("uniProtKBCrossReferences") or []:
        if xref.get("database") != "GO":
            continue
        go_id = xref.get("id", "")
        term_name = ""
        category = ""
        for prop in xref.get("properties") or []:
            if prop.get("key") == "GoTerm":
                # value is like "F:serine-type endopeptidase activity"
                raw = prop.get("value", "")
                if ":" in raw:
                    category, term_name = raw.split(":", 1)
                    term_name = term_name.strip()
                else:
                    term_name = raw
                break
        prefix = f"{category}:" if category else ""
        entry = f"{prefix}{go_id} {term_name}".strip() if term_name else f"{prefix}{go_id}".strip()
        if entry:
            go_terms.append(entry)
    result["go_terms"] = go_terms

    return result


# ---------------------------------------------------------------------------
# NCBI query
# ---------------------------------------------------------------------------

def _query_ncbi(
    accession: str,
    email: Optional[str],
    api_key: Optional[str],
) -> dict:
    global _last_ncbi_request
    _last_ncbi_request = _throttle(NCBI_INTERVAL, _last_ncbi_request)
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
    result = _empty_info(accession)
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
            # Extract nucleotide accession from coded_by — same pass, no second walk
            if not result["nucl_id"] and quals.get("coded_by"):
                coded_by = quals["coded_by"]
                stripped = re.sub(r"(?:complement|join|order)\(", "", coded_by)
                m = re.search(r"([A-Za-z]{1,2}_?[0-9][A-Za-z0-9_.]*[0-9]):", stripped)
                if m:
                    result["nucl_id"] = m.group(1)
            # EC numbers
            for ec_val in feat.findall(".//GBQualifier"):
                if ec_val.findtext("GBQualifier_name") == "EC_number":
                    ec = ec_val.findtext("GBQualifier_value") or ""
                    if ec and ec not in result["ec_numbers"]:
                        result["ec_numbers"].append(ec)
            # GO annotations (GO_function, GO_process, GO_component)
            _GO_QUAL_CAT = {"GO_function": "F", "GO_process": "P", "GO_component": "C"}
            for go_q in feat.findall(".//GBQualifier"):
                qname = go_q.findtext("GBQualifier_name") or ""
                if qname.startswith("GO_"):
                    qval = go_q.findtext("GBQualifier_value") or ""
                    # value format: "term name [GO:XXXXXXX]"
                    go_m = re.search(r"\[(GO:\d+)\]", qval)
                    go_id = go_m.group(1) if go_m else ""
                    term = re.sub(r"\s*\[GO:\d+\]", "", qval).strip()
                    category = _GO_QUAL_CAT.get(qname, "")
                    prefix = f"{category}:" if category else ""
                    entry = f"{prefix}{go_id} {term}".strip() if go_id else f"{prefix}{term}".strip()
                    if entry and entry not in result["go_terms"]:
                        result["go_terms"].append(entry)
        if key == "gene" and not result["gene"] and quals.get("gene"):
            result["gene"] = quals["gene"]

    if not result["protein_name"]:
        defn = gbseq.findtext("GBSeq_definition") or ""
        result["protein_name"] = re.sub(r"\s*\[.+?\]$", "", defn).strip()

    if result["nucl_id"]:
        result["nucl_ids"] = [result["nucl_id"]]

    return result


# ---------------------------------------------------------------------------
# Ensembl query
# ---------------------------------------------------------------------------

# Small per-session caches so walking protein → transcript → gene only costs
# one extra API call per unique parent ID, not one per sequence.
_ensembl_obj_cache: dict[str, dict] = {}
# species string (e.g. "homo_sapiens") → (taxid, lineage_str)
_ensembl_tax_cache: dict[str, tuple[str, str]] = {}


def _ensembl_lookup(eid: str) -> dict:
    """Fetch one Ensembl /lookup/id call; results are cached by ID."""
    if eid in _ensembl_obj_cache:
        return _ensembl_obj_cache[eid]
    global _last_ensembl_request
    _last_ensembl_request = _throttle(ENSEMBL_INTERVAL, _last_ensembl_request)
    url = (
        f"https://rest.ensembl.org/lookup/id/{urllib.parse.quote(eid)}"
        "?expand=0&content-type=application/json"
    )
    data = json.loads(_get(url))
    _ensembl_obj_cache[eid] = data
    return data


def _ensembl_taxonomy(species: str) -> tuple[str, str]:
    """Return ``(taxid, lineage_str)`` for an Ensembl species name like ``homo_sapiens``.

    Results are cached per species so a batch of same-organism sequences only
    pays for two extra API calls once.
    """
    if species in _ensembl_tax_cache:
        return _ensembl_tax_cache[species]
    sci_name = species.replace("_", " ").capitalize()
    taxid = lineage = ""
    try:
        global _last_ensembl_request
        _last_ensembl_request = _throttle(ENSEMBL_INTERVAL, _last_ensembl_request)
        nodes = json.loads(
            _get(
                f"https://rest.ensembl.org/taxonomy/name/{urllib.parse.quote(sci_name)}"
                "?content-type=application/json"
            )
        )
        if nodes:
            taxid = str(nodes[0].get("id", ""))
        if taxid:
            _last_ensembl_request = _throttle(ENSEMBL_INTERVAL, _last_ensembl_request)
            cls = json.loads(
                _get(
                    f"https://rest.ensembl.org/taxonomy/classification/{taxid}"
                    "?content-type=application/json"
                )
            )
            lineage = "; ".join(
                c["name"]
                for c in reversed(cls)
                if c.get("rank") not in ("no rank", None, "")
            )
    except Exception:  # noqa: BLE001
        pass
    _ensembl_tax_cache[species] = (taxid, lineage)
    return taxid, lineage


def _query_ensembl(accession: str) -> dict:
    """Fetch gene/protein metadata from the Ensembl REST API.

    Walks up the object hierarchy (Translation → Transcript → Gene) so that
    the gene description and display name are always populated regardless of
    which Ensembl ID type is supplied.
    """
    global _last_ensembl_request
    _last_ensembl_request = _throttle(ENSEMBL_INTERVAL, _last_ensembl_request)

    data = _ensembl_lookup(accession)
    obj_type = data.get("object_type", "")
    species = data.get("species", "")

    result = _empty_info(accession, data.get("id", accession))
    result["organism"] = species.replace("_", " ").capitalize()

    # Walk up the hierarchy to the Gene object which carries the description
    gene_data = data
    if obj_type == "Translation":
        transcript_id = data.get("Parent", "")
        if transcript_id:
            transcript_data = _ensembl_lookup(transcript_id)
            gene_id = transcript_data.get("Parent", "")
            gene_data = _ensembl_lookup(gene_id) if gene_id else transcript_data
            # Expose parent transcript as the nucleotide accession
            result["nucl_id"] = transcript_id
            result["nucl_ids"] = [transcript_id]
    elif obj_type == "Transcript":
        gene_id = data.get("Parent", "")
        if gene_id:
            gene_data = _ensembl_lookup(gene_id)
        result["nucl_id"] = accession
        result["nucl_ids"] = [accession]

    desc = gene_data.get("description", "") or ""
    result["protein_name"] = re.sub(r"\s*\[Source:[^\]]*\]", "", desc).strip()
    result["gene"] = gene_data.get("display_name", "") or data.get("display_name", "")

    if species:
        result["taxid"], result["lineage"] = _ensembl_taxonomy(species)

    # GO terms via /xrefs/id endpoint (best-effort; silently skipped on error)
    try:
        _last_ensembl_request = _throttle(ENSEMBL_INTERVAL, _last_ensembl_request)
        xrefs = json.loads(_get(
            f"https://rest.ensembl.org/xrefs/id/{urllib.parse.quote(accession)}"
            "?content-type=application/json"
        ))
        go_terms: list[str] = []
        for xref in xrefs:
            if xref.get("dbname") == "GO":
                go_id = xref.get("primary_id", "")
                go_name = xref.get("display_id", "") or xref.get("description", "")
                entry = f"{go_id} {go_name}".strip() if go_name and go_name != go_id else go_id
                if entry and entry not in go_terms:
                    go_terms.append(entry)
        result["go_terms"] = go_terms
    except Exception:  # noqa: BLE001
        pass

    return result


# ---------------------------------------------------------------------------
# InterPro query
# ---------------------------------------------------------------------------

def _query_interpro(accession: str) -> dict:
    """Fetch entry metadata from the EBI InterPro REST API.

    InterPro entries represent protein families / domains rather than
    individual sequences, so organism / lineage fields are left empty.
    The entry type (e.g. ``Domain``, ``Family``) is appended to the name.
    """
    global _last_interpro_request
    _last_interpro_request = _throttle(INTERPRO_INTERVAL, _last_interpro_request)
    url = (
        f"https://www.ebi.ac.uk/interpro/api/entry/interpro/{accession.upper()}/"
        "?format=json"
    )
    data = json.loads(_get(url))
    meta = data.get("metadata", {})

    result = _empty_info(accession.upper(), meta.get("accession", accession.upper()))
    name_block = meta.get("name", {})
    entry_name = name_block.get("name", "") or name_block.get("short", "")
    entry_type = meta.get("type", "")
    result["protein_name"] = f"{entry_name} [{entry_type}]" if entry_type else entry_name
    result["gene"] = name_block.get("short", "")
    # GO terms from InterPro metadata
    go_terms: list[str] = []
    for go in meta.get("go_terms") or []:
        go_id = go.get("identifier", "")
        go_name = go.get("name", "")
        category = (go.get("category") or {}).get("code", "")
        prefix = f"{category}:" if category else ""
        entry = f"{prefix}{go_id} {go_name}".strip() if go_name else f"{prefix}{go_id}".strip()
        if entry:
            go_terms.append(entry)
    result["go_terms"] = go_terms
    # Organism / lineage are not meaningful for InterPro entries
    return result


# ---------------------------------------------------------------------------
# PDB query
# ---------------------------------------------------------------------------

def _query_pdb(accession: str) -> dict:
    """Fetch structure metadata from the RCSB PDB REST API.

    Accession may be a bare 4-character PDB ID (``4HHB``) or include an
    entity number (``4HHB_1``).  When no entity is given, entity ``1`` is
    used.  Organism and gene data are taken from the first source organism of
    the polymer entity.
    """
    global _last_pdb_request
    _last_pdb_request = _throttle(PDB_INTERVAL, _last_pdb_request)

    parts = accession.upper().split("_", 1)
    pdb_id = parts[0]
    entity_id = parts[1] if len(parts) > 1 else "1"

    result = _empty_info(accession.upper(), f"{pdb_id}_{entity_id}")

    # Polymer-entity endpoint contains description, organism, and gene
    entity_url = (
        f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
    )
    try:
        entity = json.loads(_get(entity_url))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            # Fall back to entry-level title only
            entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            entry = json.loads(_get(entry_url))
            result["protein_name"] = (entry.get("struct") or {}).get("title", "")
            return result
        raise

    pe = entity.get("rcsb_polymer_entity") or {}
    result["protein_name"] = pe.get("pdbx_description", "") or pe.get("pdbx_fragment", "")

    sources = entity.get("rcsb_entity_source_organism") or []
    if sources:
        src = sources[0]
        result["organism"] = src.get("ncbi_scientific_name", "")
        taxid = src.get("ncbi_taxonomy_id")
        result["taxid"] = str(taxid) if taxid else ""
        gene_names = src.get("rcsb_gene_name") or []
        if gene_names:
            result["gene"] = gene_names[0].get("value", "")

    # EC numbers from combined enzyme class list
    ec_numbers: list[str] = []
    for ec_entry in pe.get("rcsb_enzyme_class_combined") or []:
        ec = ec_entry.get("ec", "")
        if ec and ec not in ec_numbers:
            ec_numbers.append(ec)
    result["ec_numbers"] = ec_numbers

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
        elif db == "ensembl":
            info = _query_ensembl(accession)
        elif db == "interpro":
            info = _query_interpro(accession)
        elif db == "pdb":
            info = _query_pdb(accession)
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
    # GO terms
    "go": "go_terms",
    "go_terms": "go_terms",
    "go_ids": "go_terms",
    # EC numbers
    "ec": "ec_numbers",
    "ec_numbers": "ec_numbers",
    "ec_number": "ec_numbers",
}


# Cache: format string → list of (placeholder, canonical_field) for only the aliases
# actually present in that string.  Built once per unique format string.
_fmt_substitutions: dict[str, list[tuple[str, str]]] = {}


_GO_CAT_PREFIX_RE = re.compile(r"^[A-Z]:")


def _strip_go_prefix(term: str) -> str:
    """Remove the leading ``"X:"`` category prefix from a stored GO term string."""
    return _GO_CAT_PREFIX_RE.sub("", term, count=1)


def format_header(
    fmt: str,
    info: dict,
    use_underscores: bool = True,
    max_go: Optional[int] = None,
    max_ec: Optional[int] = None,
    go_category: Optional[set] = None,
) -> str:
    if fmt not in _fmt_substitutions:
        _fmt_substitutions[fmt] = [
            (f"{{{alias}}}", field)
            for alias, field in _ALIASES.items()
            if f"{{{alias}}}" in fmt
        ]
    out = fmt
    for placeholder, field in _fmt_substitutions[fmt]:
        val = info.get(field, "") or ""
        if isinstance(val, list):
            if field == "go_terms":
                if go_category is not None:
                    # Keep only terms whose leading category code is in the requested set.
                    # Terms with no prefix (unknown category) are excluded when filtering.
                    val = [t for t in val if t[:2] in {f"{c}:" for c in go_category}]
                if max_go is not None:
                    val = val[:max_go]
                val = [_strip_go_prefix(t) for t in val]
            elif field == "ec_numbers" and max_ec is not None:
                val = val[:max_ec]
            val = ";".join(val)
        out = out.replace(placeholder, val)
    out = out.strip()
    if use_underscores:
        out = re.sub(r"\s+", "_", out)
    return out


# ---------------------------------------------------------------------------
# Per-record result
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class RecordResult:
    source_file: str        # input filename (basename)
    original_header: str
    accession: str
    db_used: str            # 'uniprot' | 'ncbi'
    seq_len: int
    relabeled: bool         # False if lookup failed
    new_header: str
    organism: str = ""
    protein_name: str = ""
    gene: str = ""
    taxid: str = ""
    lineage: str = ""
    nucl_id: str = ""                                           # chosen nucleotide/genome accession
    nucl_ids: list = dataclasses.field(default_factory=list)    # full list of linked nucleotide accessions
    go_terms: list = dataclasses.field(default_factory=list)    # GO term strings e.g. ["GO:0004252 serine-type endopeptidase activity"]
    ec_numbers: list = dataclasses.field(default_factory=list)  # EC numbers e.g. ["3.4.21.4"]
    error: str = ""                                             # non-empty when lookup failed
    filtered: bool = False                                      # True when dropped by a sequence filter
    filter_reason: str = ""                                     # human-readable reason for filtering


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

_POST_LOOKUP_PREFIXES = ("organism", "taxid", "lineage")


def _is_post_lookup_reason(reason: str) -> bool:
    return any(reason.startswith(p) for p in _POST_LOOKUP_PREFIXES)


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

    # ── Filters ──────────────────────────────────────────────────────────────
    filtered = [r for r in results if r.filtered]
    active   = [r for r in results if not r.filtered]
    if filtered:
        pre_filtered  = [r for r in filtered if not _is_post_lookup_reason(r.filter_reason)]
        post_filtered = [r for r in filtered if _is_post_lookup_reason(r.filter_reason)]
        lines += [_section(f"FILTERS  ({len(filtered)} removed, {len(active)} kept)")]
        if pre_filtered:
            lines.append("  Pre-lookup (sequence):")
            pre_counts: Counter = Counter(r.filter_reason for r in pre_filtered)
            for reason, cnt in pre_counts.most_common():
                lines.append(f"    {reason:<46} {_fmt_int(cnt):>6}")
        if post_filtered:
            lines.append("  Post-lookup (organism / taxid / lineage):")
            post_counts: Counter = Counter(r.filter_reason for r in post_filtered)
            for reason, cnt in post_counts.most_common():
                lines.append(f"    {reason:<46} {_fmt_int(cnt):>6}")

    total = len(active)
    relabeled = sum(1 for r in active if r.relabeled)
    failed = total - relabeled
    unique_acc = len({r.accession for r in active})
    with_nucl = sum(1 for r in active if r.relabeled and r.nucl_id)
    lengths = [r.seq_len for r in active]

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
    org_counts = Counter(r.organism or "(unknown)" for r in active if r.relabeled)
    _render_distribution(lines, "ORGANISM DISTRIBUTION", org_counts, total, top_n)

    # ── Protein distribution ─────────────────────────────────────────────────
    prot_counts = Counter(r.protein_name or "(unknown)" for r in active if r.relabeled)
    _render_distribution(lines, "PROTEIN DISTRIBUTION", prot_counts, total, top_n)

    # ── EC number distribution ────────────────────────────────────────────────
    ec_all = [ec for r in active if r.relabeled for ec in r.ec_numbers]
    if ec_all:
        ec_counts: Counter = Counter(ec_all)
        _render_distribution(lines, "EC NUMBER DISTRIBUTION", ec_counts, total, top_n)

    # ── GO term distribution ──────────────────────────────────────────────────
    go_all = [gt for r in active if r.relabeled for gt in r.go_terms]
    if go_all:
        go_counts: Counter = Counter(go_all)
        _render_distribution(lines, "GO TERM DISTRIBUTION", go_counts, total, top_n)

    # ── Nucleotide ID coverage per organism ──────────────────────────────────
    if with_nucl:
        lines += [_section("NUCLEOTIDE ID COVERAGE  (relabeled sequences only)")]
        # For each organism: how many have a nucl_id vs not
        org_nucl: dict[str, list[bool]] = {}
        for r in active:
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
    by_file = _group_by_file(active)

    col_w = max(len(f) for f in by_file) + 2 if by_file else 24
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
    failed_recs = [r for r in active if not r.relabeled]
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
    "nucl_id", "nucl_ids", "lineage", "go_terms", "ec_numbers",
    "new_header", "error", "filtered", "filter_reason",
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
    active = [r for r in results if not r.filtered]
    filtered = [r for r in results if r.filtered]
    total = len(active)
    relabeled = sum(1 for r in active if r.relabeled)
    lengths = [r.seq_len for r in active]
    with_nucl = sum(1 for r in active if r.relabeled and r.nucl_id)
    org_counts = Counter(r.organism or "" for r in active if r.relabeled)
    prot_counts = Counter(r.protein_name or "" for r in active if r.relabeled)
    ec_counts = Counter(ec for r in active if r.relabeled for ec in r.ec_numbers)
    go_counts = Counter(gt for r in active if r.relabeled for gt in r.go_terms)

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
        "filters": {
            "filtered_count": len(filtered),
            "kept_count": total,
            "pre_lookup": dict(
                Counter(
                    r.filter_reason for r in filtered
                    if not _is_post_lookup_reason(r.filter_reason)
                ).most_common()
            ),
            "post_lookup": dict(
                Counter(
                    r.filter_reason for r in filtered
                    if _is_post_lookup_reason(r.filter_reason)
                ).most_common()
            ),
        },
        "coverage": {
            "total_sequences": total,
            "relabeled": relabeled,
            "failed": total - relabeled,
            "unique_accessions": len({r.accession for r in active}),
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
        "ec_number_distribution": dict(ec_counts.most_common()),
        "go_term_distribution": dict(go_counts.most_common()),
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
            for fname, recs in _group_by_file(active).items()
        },
        "records": [dataclasses.asdict(r) for r in results],
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
    if not _check_overwrite(path, overwrite):
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return True


# ---------------------------------------------------------------------------
# Progress bar
# ---------------------------------------------------------------------------

class ProgressBar:
    """Minimal in-place terminal progress bar written to stderr.

    Pass ``total=0`` for *streaming mode*: a simple counter is shown instead
    of a bar and percentage (total is unknown).  Automatically disabled when
    stderr is not a TTY or *enabled* is False.  Uses ``\\r`` to overwrite the
    current line; call :meth:`clear` before any other stderr print.
    """

    BAR_WIDTH = 26

    def __init__(self, total: int, desc: str = "", enabled: bool = True) -> None:
        self.total = total          # 0 means "unknown / streaming"
        self.current = 0
        self._desc = (desc[:22] + "…") if len(desc) > 23 else desc.ljust(23)
        self._streaming = (total == 0)
        self._enabled = enabled and sys.stderr.isatty()
        self._line_len = 0

    def update(self, n: int = 1) -> None:
        """Advance the counter by *n* and redraw."""
        self.current += n
        if self._enabled:
            self._draw()

    def _draw(self) -> None:
        if self._streaming:
            line = f"  {self._desc} {self.current:,} records"
        else:
            frac = self.current / max(self.total, 1)
            filled = round(self.BAR_WIDTH * frac)
            bar = "█" * filled + "░" * (self.BAR_WIDTH - filled)
            line = (
                f"  {self._desc} [{bar}] "
                f"{self.current}/{self.total} ({frac * 100:.0f}%)"
            )
        self._line_len = len(line)
        print(f"\r{line}", end="", file=sys.stderr, flush=True)

    def clear(self) -> None:
        """Erase the bar line so a normal ``print`` can follow cleanly."""
        if self._enabled:
            print(f"\r{' ' * (self._line_len + 2)}\r", end="", file=sys.stderr, flush=True)
            self._enabled = False


# ---------------------------------------------------------------------------
# Processing configuration
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class ProcessConfig:
    """All per-run settings passed down to process_file and fetch_info."""
    fmt: str
    db: str
    email: Optional[str]
    api_key: Optional[str]
    nucl_type: str
    use_underscores: bool
    overwrite: bool
    id_delimiter: Optional[str]
    id_field: Optional[int]
    verbose: bool
    dry_run: bool
    show_progress: bool = True
    min_len: Optional[int] = None
    max_len: Optional[int] = None
    max_ambiguous_ratio: Optional[float] = None
    dedupe_sequence: bool = False
    sample: Optional[int] = None
    molecule_type: str = "auto"          # 'auto' | 'prot' | 'nucl'
    organism_include: list = dataclasses.field(default_factory=list)   # compiled re.Patterns
    organism_exclude: list = dataclasses.field(default_factory=list)   # compiled re.Patterns
    taxid_include: set = dataclasses.field(default_factory=set)        # str taxids
    taxid_exclude: set = dataclasses.field(default_factory=set)        # str taxids
    lineage_include: list = dataclasses.field(default_factory=list)    # compiled re.Patterns
    lineage_exclude: list = dataclasses.field(default_factory=list)    # compiled re.Patterns
    streaming: bool = False
    max_go: Optional[int] = None         # cap on GO terms included in header
    max_ec: Optional[int] = None         # cap on EC numbers included in header
    go_category: Optional[set] = None    # e.g. {"F", "P", "C"}; None = all categories


# ---------------------------------------------------------------------------
# File processing
# ---------------------------------------------------------------------------

def _run_lookups(
    records: list[tuple[str, str]],
    cfg: ProcessConfig,
    bar: "ProgressBar",
) -> list[tuple[str, str, Optional[dict], str]]:
    """Run DB lookups for all records, returning ``(db, accession, info, error)`` tuples."""
    out: list[tuple[str, str, Optional[dict], str]] = []
    for header, _seq in records:
        detected_db, accession = detect_id(header, cfg.db, cfg.id_delimiter, cfg.id_field)
        info, error = fetch_info(
            accession, detected_db, cfg.email, cfg.api_key, cfg.nucl_type, cfg.verbose
        )
        out.append((detected_db, accession, info, error))
        bar.update()
    return out


def process_file(
    input_path: Path,
    output_path: Path,
    cfg: ProcessConfig,
) -> list[RecordResult]:
    if cfg.streaming:
        return _process_file_streaming(input_path, output_path, cfg)

    records = list(parse_fasta(input_path))
    if not records:
        print(f"Warning: no records found in {input_path}", file=sys.stderr)
        return []

    # ── Pre-lookup filters (length / ambiguity / dedupe / sample) ────────
    any_pre_filter = (
        cfg.min_len is not None
        or cfg.max_len is not None
        or cfg.max_ambiguous_ratio is not None
        or cfg.dedupe_sequence
        or cfg.sample is not None
    )
    filtered_results: list[RecordResult] = []
    if any_pre_filter:
        records, filtered_results = apply_sequence_filters(
            records,
            source_file=input_path.name,
            min_len=cfg.min_len,
            max_len=cfg.max_len,
            max_ambiguous_ratio=cfg.max_ambiguous_ratio,
            dedupe_sequence=cfg.dedupe_sequence,
            sample=cfg.sample,
            molecule_type=cfg.molecule_type,
        )
        if cfg.verbose and filtered_results:
            print(
                f"  Pre-lookup filters: removed {len(filtered_results)} record(s) "
                f"({len(records)} remain)",
                file=sys.stderr,
            )

    if cfg.verbose:
        print(f"\n{input_path}  ({len(records)} records)", file=sys.stderr)

    # Progress bar: active only when verbose is off and stderr is a TTY
    bar = ProgressBar(
        total=len(records),
        desc=input_path.name,
        enabled=cfg.show_progress and not cfg.verbose and not cfg.dry_run,
    )

    # ── Phase 1: DB lookups ───────────────────────────────────────────────
    lookups = _run_lookups(records, cfg, bar)

    bar.clear()  # erase bar line before any final print

    # ── Phase 2: Build RecordResult list in original order ────────────────
    results: list[RecordResult] = []
    new_records: list[tuple[str, str]] = []

    for (header, seq), (detected_db, accession, info, error) in zip(records, lookups):
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
            new_header = format_header(cfg.fmt, info, cfg.use_underscores, cfg.max_go, cfg.max_ec, cfg.go_category)
            if cfg.verbose:
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
                go_terms=info.get("go_terms", []),
                ec_numbers=info.get("ec_numbers", []),
            )
            new_records.append((new_header, seq))
        results.append(res)

    # ── Post-lookup filters (organism / taxid / lineage) ─────────────────
    post_dropped: list[RecordResult] = []
    results, new_records, post_dropped = apply_post_lookup_filters(
        results, new_records,
        cfg.organism_include, cfg.organism_exclude,
        cfg.taxid_include,    cfg.taxid_exclude,
        cfg.lineage_include,  cfg.lineage_exclude,
    )
    if cfg.verbose and post_dropped:
        print(
            f"  Post-lookup filters: removed {len(post_dropped)} record(s) "
            f"({len(results)} remain)",
            file=sys.stderr,
        )

    if cfg.dry_run:
        for h, _ in new_records:
            print(f">{h}")
    else:
        if write_fasta(output_path, new_records, overwrite=cfg.overwrite):
            print(f"Wrote {output_path}", file=sys.stderr)

    # Filtered records are appended at the end so they appear in reports but
    # are excluded from the output FASTA.
    return results + filtered_results + post_dropped


def _process_file_streaming(
    input_path: Path,
    output_path: Path,
    cfg: ProcessConfig,
) -> list[RecordResult]:
    """Process a FASTA file in streaming mode.

    Sequences are consumed one at a time so only a small working set is held
    in memory.  Output is written incrementally to avoid buffering all records.

    Behavioural differences from batch mode:
    - The progress bar shows a plain counter (total is unknown).
    - ``--sample N`` uses Vitter's reservoir sampling (Algorithm R).
    - ``--dedupe-sequence`` stores a hash per unique sequence rather than the
      full string, trading a negligible collision risk for lower memory usage.
    """
    if cfg.verbose:
        print(f"\n{input_path}  (streaming)", file=sys.stderr)

    bar = ProgressBar(
        total=0,   # 0 → streaming / count-only mode
        desc=input_path.name,
        enabled=cfg.show_progress and not cfg.verbose and not cfg.dry_run,
    )

    # ── Molecule-type detection (peek at first record if needed) ──────────
    mol = cfg.molecule_type
    amb_chars: frozenset = _AMBIGUOUS_AA  # resolved on first record if mol=='auto'
    mol_resolved = (mol != "auto")

    # ── Per-record filter state ───────────────────────────────────────────
    seen_hashes: set[int] = set()          # for --dedupe-sequence (hash-based)

    # ── Reservoir for --sample (Vitter's Algorithm R) ─────────────────────
    # Each entry: (out_record, RecordResult) where out_record=(header, seq)
    reservoir: list[tuple[tuple[str, str], RecordResult]] = []
    stream_total = 0          # total records that reached the reservoir stage

    results:          list[RecordResult] = []
    filtered_results: list[RecordResult] = []

    # Open output file early so we can write incrementally (skipped in dry-run)
    out_fh = None
    if not cfg.dry_run:
        if not _check_overwrite(output_path, cfg.overwrite):
            return []
        output_path.parent.mkdir(parents=True, exist_ok=True)
        out_fh = open(output_path, "w")

    try:
        for header, seq in parse_fasta(input_path):

            # ── Resolve molecule type from first sequence ─────────────────
            if not mol_resolved:
                mol = _detect_molecule_type(seq)
                amb_chars = _AMBIGUOUS_NUCL if mol == "nucl" else _AMBIGUOUS_AA
                mol_resolved = True

            # ── Pre-lookup filters ────────────────────────────────────────
            seq_upper = seq.upper()
            reason = ""
            if cfg.min_len is not None and len(seq) < cfg.min_len:
                reason = f"too short ({len(seq)} < {cfg.min_len})"
            elif cfg.max_len is not None and len(seq) > cfg.max_len:
                reason = f"too long ({len(seq)} > {cfg.max_len})"
            elif cfg.max_ambiguous_ratio is not None and _ambiguous_ratio(seq, amb_chars) > cfg.max_ambiguous_ratio:
                ratio = _ambiguous_ratio(seq, amb_chars)
                reason = f"ambiguous ratio {ratio:.3f} > {cfg.max_ambiguous_ratio:.3f}"
            elif cfg.dedupe_sequence:
                h = hash(seq_upper)
                if h in seen_hashes:
                    reason = "duplicate sequence"
                else:
                    seen_hashes.add(h)

            if reason:
                filtered_results.append(_filtered_result(input_path.name, header, seq, reason))
                bar.update()
                continue

            # ── DB lookup ─────────────────────────────────────────────────
            detected_db, accession = detect_id(header, cfg.db, cfg.id_delimiter, cfg.id_field)
            info, error = fetch_info(
                accession, detected_db, cfg.email, cfg.api_key, cfg.nucl_type, cfg.verbose
            )

            if info is None:
                res = RecordResult(
                    source_file=input_path.name, original_header=header,
                    accession=accession, db_used=detected_db,
                    seq_len=len(seq), relabeled=False, new_header=header, error=error,
                )
                out_rec: tuple[str, str] = (header, seq)
            else:
                info["id"] = accession
                new_header = format_header(cfg.fmt, info, cfg.use_underscores, cfg.max_go, cfg.max_ec, cfg.go_category)
                if cfg.verbose:
                    nucl_note = f"  nucl_id={info['nucl_id']!r}" if info.get("nucl_id") else ""
                    print(f"  {accession!r} → {new_header!r}{nucl_note}", file=sys.stderr)
                res = RecordResult(
                    source_file=input_path.name, original_header=header,
                    accession=accession, db_used=detected_db,
                    seq_len=len(seq), relabeled=True, new_header=new_header,
                    organism=info.get("organism", ""), protein_name=info.get("protein_name", ""),
                    gene=info.get("gene", ""), taxid=info.get("taxid", ""),
                    lineage=info.get("lineage", ""), nucl_id=info.get("nucl_id", ""),
                    nucl_ids=info.get("nucl_ids", []),
                    go_terms=info.get("go_terms", []),
                    ec_numbers=info.get("ec_numbers", []),
                )
                out_rec = (new_header, seq)

            # ── Post-lookup filters ───────────────────────────────────────
            post_reason = _post_lookup_filter_reason(
                res,
                cfg.organism_include, cfg.organism_exclude,
                cfg.taxid_include,    cfg.taxid_exclude,
                cfg.lineage_include,  cfg.lineage_exclude,
            )
            if post_reason:
                filtered_results.append(
                    dataclasses.replace(res, filtered=True, filter_reason=post_reason)
                )
                bar.update()
                continue

            # ── Reservoir sampling (Vitter's Algorithm R) ─────────────────
            if cfg.sample is not None:
                stream_total += 1
                if len(reservoir) < cfg.sample:
                    reservoir.append((out_rec, res))
                else:
                    j = random.randrange(stream_total)
                    if j < cfg.sample:
                        # Evict the displaced entry
                        evicted_rec, evicted_res = reservoir[j]
                        filtered_results.append(
                            dataclasses.replace(evicted_res, filtered=True, filter_reason="not sampled")
                        )
                        reservoir[j] = (out_rec, res)
                    else:
                        filtered_results.append(
                            dataclasses.replace(res, filtered=True, filter_reason="not sampled")
                        )
            else:
                # Write immediately
                if cfg.dry_run:
                    print(f">{out_rec[0]}")
                elif out_fh is not None:
                    _write_fasta_record(out_fh, *out_rec)
                results.append(res)

            bar.update()

    finally:
        bar.clear()
        if out_fh is not None:
            out_fh.close()

    # ── Flush reservoir ───────────────────────────────────────────────────
    if cfg.sample is not None:
        if not cfg.dry_run and reservoir:
            # Re-open to write the reservoir (file was closed above)
            with open(output_path, "a") as fh:
                for out_rec, res in reservoir:
                    _write_fasta_record(fh, *out_rec)
        for out_rec, res in reservoir:
            if cfg.dry_run:
                print(f">{out_rec[0]}")
            results.append(res)

    if not cfg.dry_run:
        print(f"Wrote {output_path}", file=sys.stderr)

    return results + filtered_results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fasta_relabel",
        description=(
            f"fasta_relabel {VERSION} — Relabel FASTA headers with taxonomy/protein info "
            "from UniProt, NCBI, Ensembl, InterPro, or PDB."
        ),
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
  {ec}             EC numbers, semicolon-separated  [aliases: ec_numbers, ec_number]
                   Source: UniProt recommendedName, NCBI CDS EC_number, PDB enzyme class
  {go}             GO terms, semicolon-separated  [aliases: go_terms, go_ids]
                   Source: UniProt xrefs, NCBI GO_* qualifiers, InterPro metadata,
                   Ensembl xrefs.  Filter by category with --go-category F/P/C.

NUCLEOTIDE ID TYPE  (--nucl-type)
  genomic  Prefer Genomic_DNA accessions, fall back to mRNA/other  (default)
  mrna     Prefer mRNA / cDNA accessions, fall back to genomic
  any      Use the first cross-reference listed regardless of molecule type
  Note: for NCBI records the CDS coded_by qualifier is used directly and is
  unaffected by --nucl-type.  For Ensembl proteins, {nucl_id} is the parent
  transcript ID (ENST…).  Not applicable for InterPro or PDB.

DATABASES  (--db)
  auto      Auto-detect from accession pattern (default)
  uniprot   Swiss-Prot / TrEMBL  — e.g. P69905, sp|P69905|HBA_HUMAN
  ncbi      GenBank / RefSeq     — e.g. NP_000549.1, AAH01234
  ensembl   Ensembl gene/transcript/protein  — e.g. ENSG00000157764, ENSP00000000233
  interpro  InterPro domain/family entries   — e.g. IPR000001
  pdb       RCSB PDB structures              — e.g. 4HHB, 4HHB_1

REPORT FORMATS
  text  Human-readable tables with ASCII bar charts (default → stderr)
  tsv   One row per sequence — import directly into Excel / R / pandas
  json  Fully structured JSON with per-record data and aggregate statistics

ID EXTRACTION
  By default every whitespace-separated field of the header is scanned left-to-right
  for a recognised accession pattern.  Use --id-delimiter and --id-field to override:

  --id-delimiter "|"            scan each pipe-separated field
  --id-delimiter "|" --id-field 1   use exactly the second pipe field
  --id-delimiter ";"            scan semicolon-separated fields
  --id-field 2                  use the third whitespace field directly

  Example headers and matching options:
    >sp|P69905|HBA_HUMAN …       auto-detected without any flags
    >P69905 extra text …         auto-detected without any flags
    >gene_A|P69905|Homo …        --id-delimiter "|" (or auto-scan finds it)
    >gene_A|P69905|Homo …        --id-delimiter "|" --id-field 1  (explicit)
    >custom_prefix NP_001290.1 … auto-scan finds it, or --id-field 1

EXAMPLES
  fasta_relabel sequences.fasta
  fasta_relabel sequences.fasta --format "{organism} | {id} | {name}"
  fasta_relabel sequences.fasta --format "{organism} | {id} | {name} | {nucl_id}"
  fasta_relabel sequences.fasta --format "{id}|{nucl_id}|{name}" --nucl-type mrna
  fasta_relabel sequences.fasta --format "{organism}|{id}|{name}|{ec}|{go}"
  fasta_relabel sequences.fasta --format "{id}|{ec}" --db uniprot
  fasta_relabel sequences.fasta --format "{id} {name}" --db ncbi --email me@example.com
  fasta_relabel seqs/ --format "{lineage}; {organism}" --output renamed/
  fasta_relabel a.fasta b.fasta --in-place --format "{id}|{taxon}|{name}|{genome_id}"
  fasta_relabel sequences.fasta --dry-run --format "{organism} {id} {name}"
  fasta_relabel sequences.fasta --report summary.txt
  fasta_relabel sequences.fasta --report data.tsv --report-format tsv
  fasta_relabel sequences.fasta --report data.json --report-format json
  fasta_relabel custom.fasta --id-delimiter "|" --format "{organism}|{id}|{name}"
  fasta_relabel custom.fasta --id-delimiter "|" --id-field 2 --format "{id} {name}"
  fasta_relabel genes.fasta --db ensembl --format "{organism}|{id}|{name}|{nucl_id}"
  fasta_relabel domains.fasta --db interpro --format "{id} {name}"
  fasta_relabel structures.fasta --db pdb --format "{organism}|{id}|{name}"
  fasta_relabel sequences.fasta --retries 5
  fasta_relabel sequences.fasta --min-len 50 --max-len 2000
  fasta_relabel sequences.fasta --max-ambiguous-ratio 0.05
  fasta_relabel sequences.fasta --dedupe-sequence
  fasta_relabel sequences.fasta --min-len 30 --dedupe-sequence --max-ambiguous-ratio 0.1
  fasta_relabel sequences.fasta --sample 500
  fasta_relabel sequences.fasta --molecule-type nucl --max-ambiguous-ratio 0.05
  fasta_relabel sequences.fasta --organism-include "Homo sapiens" "Mus musculus"
  fasta_relabel sequences.fasta --organism-exclude synthetic artificial
  fasta_relabel sequences.fasta --organism-include "^Homo" --organism-exclude "virus"
  fasta_relabel sequences.fasta --taxid-include 9606 10090
  fasta_relabel sequences.fasta --taxid-exclude 32630
  fasta_relabel sequences.fasta --lineage-include Mammalia
  fasta_relabel sequences.fasta --lineage-exclude Viridae
  fasta_relabel big.fasta --streaming
  fasta_relabel big.fasta --streaming --sample 1000 --lineage-include Vertebrata
  fasta_relabel sequences.fasta --config run.toml
  fasta_relabel sequences.fasta --config run.toml --dry-run  # CLI overrides config

CONFIG FILE  (--config FILE)
  Any long option can appear as a key in a TOML file.  Use underscores
  instead of hyphens in key names.  CLI arguments always take precedence.

  Example run.toml ─────────────────────────────────────────────────────────
    format              = "{organism}|{id}|{name}|{nucl_id}"
    db                  = "auto"
    email               = "me@example.com"

    # Pre-lookup filters
    min_len             = 50
    max_len             = 5000
    max_ambiguous_ratio = 0.1
    dedupe_sequence     = true

    # Post-lookup filters
    organism_include    = ["Homo sapiens", "Mus musculus"]
    taxid_exclude       = ["32630"]
    lineage_include     = ["Vertebrata"]

    # Output / network
    suffix              = "_clean"
    streaming           = true
    retries             = 5
  ──────────────────────────────────────────────────────────────────────────
  Requires Python 3.11+ for full TOML support (tomllib).  A built-in
  minimal parser is used automatically on older versions and handles the
  subset shown above.
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
        "--id-delimiter",
        metavar="DELIM",
        default=None,
        help="Character(s) used to split each description line into fields when "
             "searching for the accession (default: any whitespace).  Common "
             "choices: '|' for pipe-delimited headers, ',' for CSV-style.  "
             "Combine with --id-field to pin the exact field position.",
    )
    parser.add_argument(
        "--id-field",
        metavar="N",
        type=int,
        default=None,
        help="0-based index of the field that contains the accession after "
             "splitting by --id-delimiter (default: search every field for a "
             "recognised accession pattern).  Example: '1' selects the second "
             "field of a pipe-delimited header such as 'prefix|ACC|suffix'.",
    )
    parser.add_argument(
        "--db",
        choices=["auto", "uniprot", "ncbi", "ensembl", "interpro", "pdb"],
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

    # ── Network options ───────────────────────────────────────────────────────
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        metavar="N",
        help="Maximum number of HTTP retry attempts for transient errors "
             "(429, 5xx, network failures).  Uses exponential back-off starting "
             "at 1 s (default: 3).",
    )
    # ── Sequence filters ─────────────────────────────────────────────────────
    parser.add_argument(
        "--min-len",
        type=int,
        default=None,
        metavar="N",
        help="Drop sequences shorter than N residues before any DB lookup",
    )
    parser.add_argument(
        "--max-len",
        type=int,
        default=None,
        metavar="N",
        help="Drop sequences longer than N residues before any DB lookup",
    )
    parser.add_argument(
        "--max-ambiguous-ratio",
        type=float,
        default=None,
        metavar="F",
        help="Drop sequences where the fraction of ambiguous / undetermined "
             "amino-acid characters (B, X, Z, J) exceeds F (0.0 – 1.0).  "
             "Example: --max-ambiguous-ratio 0.1 removes any sequence with "
             "more than 10%% ambiguous residues.",
    )
    parser.add_argument(
        "--dedupe-sequence",
        action="store_true",
        help="Keep only the first occurrence of each unique sequence; "
             "exact duplicates (case-insensitive) are dropped silently.",
    )
    parser.add_argument(
        "--sample",
        type=int,
        default=None,
        metavar="N",
        help="Randomly subsample N sequences per file after all other "
             "pre-lookup filters have been applied.  Uses uniform sampling "
             "without replacement.  If N ≥ the number of passing sequences "
             "all are kept.",
    )
    parser.add_argument(
        "--max-go",
        type=int,
        default=None,
        metavar="N",
        help="Include at most N GO terms in the {go} placeholder (terms are "
             "kept in the order returned by the database).  Default: all.",
    )
    parser.add_argument(
        "--max-ec",
        type=int,
        default=None,
        metavar="N",
        help="Include at most N EC numbers in the {ec} placeholder.  "
             "Default: all.",
    )
    parser.add_argument(
        "--go-category",
        nargs="+",
        metavar="CAT",
        default=None,
        help="Restrict {go} output to one or more GO categories: "
             "F (molecular function), P (biological process), "
             "C (cellular component).  Applied before --max-go.  "
             "GO terms with unknown category (e.g. from Ensembl) are "
             "excluded when this filter is active.  "
             "Example: --go-category F P",
    )
    parser.add_argument(
        "--molecule-type",
        choices=["auto", "prot", "nucl"],
        default="auto",
        metavar="TYPE",
        help="Molecule type used for --max-ambiguous-ratio: 'auto' (detect "
             "from first sequence), 'prot' (amino-acid ambiguity: B/X/Z/J), "
             "or 'nucl' (nucleotide IUPAC ambiguity: N/R/Y/S/W/K/M/B/D/H/V).  "
             "Default: auto.",
    )
    parser.add_argument(
        "--organism-include",
        nargs="+",
        default=[],
        metavar="PATTERN",
        help="After DB lookup, keep only sequences whose organism name matches "
             "at least one of the given patterns (case-insensitive regex).  "
             "Records with an unknown organism (lookup failed) are always kept.  "
             "Example: --organism-include 'Homo sapiens' 'Mus musculus'",
    )
    parser.add_argument(
        "--organism-exclude",
        nargs="+",
        default=[],
        metavar="PATTERN",
        help="After DB lookup, drop sequences whose organism name matches any "
             "of the given patterns (case-insensitive regex).  "
             "Example: --organism-exclude 'synthetic' 'artificial'",
    )
    parser.add_argument(
        "--taxid-include",
        nargs="+",
        default=[],
        metavar="TAXID",
        help="After DB lookup, keep only sequences whose NCBI taxon ID is in "
             "the given list.  Exact integer match; records with no taxid "
             "(lookup failed) are always kept.  "
             "Example: --taxid-include 9606 10090  (human + mouse)",
    )
    parser.add_argument(
        "--taxid-exclude",
        nargs="+",
        default=[],
        metavar="TAXID",
        help="After DB lookup, drop sequences whose NCBI taxon ID is in the "
             "given list.  Example: --taxid-exclude 9606",
    )
    parser.add_argument(
        "--lineage-include",
        nargs="+",
        default=[],
        metavar="PATTERN",
        help="After DB lookup, keep only sequences whose full taxonomic lineage "
             "matches at least one pattern (case-insensitive regex).  Matches "
             "against the semicolon-separated lineage string.  "
             "Example: --lineage-include Mammalia Aves  (keep mammals or birds)",
    )
    parser.add_argument(
        "--lineage-exclude",
        nargs="+",
        default=[],
        metavar="PATTERN",
        help="After DB lookup, drop sequences whose lineage matches any "
             "pattern.  Example: --lineage-exclude Viridae  (drop all viruses)",
    )
    parser.add_argument(
        "--streaming",
        action="store_true",
        help="Process sequences one at a time without loading the full file "
             "into memory.  Useful for very large FASTA files.  "
             "Uses Vitter's reservoir sampling for --sample and hash-based "
             "deduplication for --dedupe-sequence.  "
             "The progress bar switches to a counter display (no percentage).",
    )

    # ── Config file ──────────────────────────────────────────────────────────
    parser.add_argument(
        "--config",
        metavar="FILE",
        help="Load default settings from a TOML file.  Command-line arguments "
             "always take precedence over config file values.  "
             "See CONFIG FILE section below for the expected format.",
    )

    # ── Misc ─────────────────────────────────────────────────────────────────
    parser.add_argument(
        "--spaces", action="store_true",
        help="Keep spaces in output headers (default: replace all whitespace "
             "with underscores so each header is a single token)",
    )
    parser.add_argument(
        "--no-progress", action="store_true",
        help="Disable the per-file progress bar (auto-hidden when stderr is not "
             "a TTY or --verbose is active)",
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


# ---------------------------------------------------------------------------
# Config file loading  (TOML)
# ---------------------------------------------------------------------------

# Try the stdlib tomllib module (Python 3.11+); fall back to a minimal parser.
try:
    import tomllib as _tomllib
    def _toml_loads(text: str) -> dict:
        return _tomllib.loads(text)
except ImportError:
    _tomllib = None  # type: ignore[assignment]
    def _toml_loads(text: str) -> dict:  # type: ignore[misc]
        return _parse_toml_minimal(text)


def _strip_inline_comment(s: str) -> str:
    """Remove a trailing ``# comment`` from a value string, respecting strings."""
    in_str: Optional[str] = None
    for i, c in enumerate(s):
        if in_str:
            if c == in_str:
                in_str = None
        elif c in ('"', "'"):
            in_str = c
        elif c == "#":
            return s[:i].rstrip()
    return s


def _split_toml_array(inner: str) -> list[str]:
    """Split a comma-separated TOML array interior, respecting quoted strings."""
    items: list[str] = []
    depth = 0
    in_str: Optional[str] = None
    buf: list[str] = []
    for c in inner:
        if in_str:
            buf.append(c)
            if c == in_str:
                in_str = None
        elif c in ('"', "'"):
            in_str = c
            buf.append(c)
        elif c == "[":
            depth += 1
            buf.append(c)
        elif c == "]":
            depth -= 1
            buf.append(c)
        elif c == "," and depth == 0:
            item = "".join(buf).strip()
            if item:
                items.append(item)
            buf = []
        else:
            buf.append(c)
    item = "".join(buf).strip()
    if item:
        items.append(item)
    return items


def _parse_toml_value(s: str):
    """Parse a single TOML scalar or array value string."""
    s = s.strip()
    # Quoted string
    if len(s) >= 2 and s[0] in ('"', "'") and s[-1] == s[0]:
        return s[1:-1]
    # Array (possibly after multiline join)
    if s.startswith("[") and s.endswith("]"):
        inner = s[1:-1].strip()
        return [_parse_toml_value(item.strip()) for item in _split_toml_array(inner) if item.strip()]
    # Boolean
    if s == "true":
        return True
    if s == "false":
        return False
    # Integer
    try:
        return int(s)
    except ValueError:
        pass
    # Float
    try:
        return float(s)
    except ValueError:
        pass
    # Bare string (no surrounding quotes)
    return s


def _parse_toml_minimal(text: str) -> dict:
    """Minimal TOML parser supporting the subset used by fasta_relabel configs.

    Handles: bare strings, quoted strings, integers, floats, booleans,
    single-line arrays, multi-line arrays, inline comments, and section
    headers (section names are ignored — all keys are treated as top-level).
    """
    result: dict = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line or line.startswith("#"):
            continue
        # Section headers → ignored (keys from all sections are top-level)
        if line.startswith("["):
            continue
        if "=" not in line:
            continue
        key, _, rest = line.partition("=")
        key = key.strip()
        rest = _strip_inline_comment(rest.strip())
        # Multi-line array: opening bracket without closing bracket on same line
        if rest.startswith("[") and "]" not in rest:
            parts = [rest]
            while i < len(lines):
                cont = lines[i].strip()
                i += 1
                parts.append(cont)
                if "]" in cont:
                    break
            rest = " ".join(parts)
        result[key] = _parse_toml_value(rest)
    return result


def _flatten_toml(d: dict, _out: Optional[dict] = None) -> dict:
    """Recursively flatten a nested TOML dict so all keys are top-level.

    Nested section dicts produced by ``tomllib`` (e.g. ``{"filters":
    {"min_len": 50}}``) are merged into the root; section names are dropped.
    Later values win on key collision.
    """
    if _out is None:
        _out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            _flatten_toml(v, _out)
        else:
            _out[k] = v
    return _out


# Config keys that differ from their argparse dest name
_CONFIG_KEY_REMAP: dict[str, str] = {
    "format": "fmt",        # --format → dest "fmt"
    "ncbi_api_key": "ncbi_api_key",  # identity, but explicit for clarity
}

# List-valued options whose items must be strings (taxids may be written as
# bare integers in TOML; argparse always delivers them as strings from CLI)
_CONFIG_LIST_STR_FIELDS = frozenset({
    "taxid_include", "taxid_exclude",
    "organism_include", "organism_exclude",
    "lineage_include", "lineage_exclude",
    "input",
})


def _load_config(path: Path) -> dict:
    """Parse *path* as a TOML config file and return an argparse-defaults dict.

    Keys are mapped to their argparse ``dest`` names.  List fields that
    correspond to string options are normalised to ``list[str]``.  Unknown
    keys are passed through unchanged (argparse ignores extra namespace attrs).
    """
    text = path.read_text(encoding="utf-8")
    raw = _flatten_toml(_toml_loads(text))

    result: dict = {}
    for k, v in raw.items():
        dest = _CONFIG_KEY_REMAP.get(k, k)
        # Normalise list fields to list[str]
        if dest in _CONFIG_LIST_STR_FIELDS:
            if isinstance(v, list):
                v = [str(x) for x in v]
            else:
                v = [str(v)]
        result[dest] = v
    return result


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
        "--help"    if a == "-help"    else
        "--version" if a == "-version" else a
        for a in argv
    ]

    parser = build_parser()

    # ── Config file: applied as defaults BEFORE the full parse ───────────────
    # We scan argv manually so that --config works even when it appears anywhere
    # in the argument list, and so config defaults are set before parse_args().
    config_path: Optional[Path] = None
    for idx, tok in enumerate(argv):
        if tok == "--config" and idx + 1 < len(argv):
            config_path = Path(argv[idx + 1])
            break
        if tok.startswith("--config="):
            config_path = Path(tok[len("--config="):])
            break

    if config_path is not None:
        if not config_path.is_file():
            parser.error(f"--config: file not found: {config_path}")
        try:
            config_defaults = _load_config(config_path)
        except Exception as exc:
            parser.error(f"--config: could not parse {config_path}: {exc}")
        parser.set_defaults(**config_defaults)

    args = parser.parse_args(argv)

    if args.in_place and args.output:
        parser.error("--in-place and --output are mutually exclusive")

    if args.ncbi_api_key:
        global NCBI_INTERVAL
        NCBI_INTERVAL = 0.1  # ≤10 req/s with an API key

    global HTTP_RETRIES
    HTTP_RETRIES = args.retries

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

    if args.max_ambiguous_ratio is not None and not 0.0 <= args.max_ambiguous_ratio <= 1.0:
        parser.error("--max-ambiguous-ratio must be between 0.0 and 1.0")
    if args.sample is not None and args.sample < 1:
        parser.error("--sample must be a positive integer")
    if args.max_go is not None and args.max_go < 1:
        parser.error("--max-go must be a positive integer")
    if args.max_ec is not None and args.max_ec < 1:
        parser.error("--max-ec must be a positive integer")
    go_category: Optional[set] = None
    if args.go_category is not None:
        valid = {"F", "P", "C"}
        bad = {c.upper() for c in args.go_category} - valid
        if bad:
            parser.error(f"--go-category: invalid value(s) {bad}; use F, P, and/or C")
        go_category = {c.upper() for c in args.go_category}

    # Compile regex patterns early so bad patterns are caught before processing
    org_include_pats  = _compile_patterns(args.organism_include)
    org_exclude_pats  = _compile_patterns(args.organism_exclude)
    lin_include_pats  = _compile_patterns(args.lineage_include)
    lin_exclude_pats  = _compile_patterns(args.lineage_exclude)
    taxid_include_set = set(args.taxid_include)
    taxid_exclude_set = set(args.taxid_exclude)

    cfg = ProcessConfig(
        fmt=args.fmt,
        db=args.db,
        email=args.email,
        api_key=args.ncbi_api_key,
        nucl_type=args.nucl_type,
        use_underscores=not args.spaces,
        overwrite=overwrite,
        id_delimiter=args.id_delimiter,
        id_field=args.id_field,
        verbose=args.verbose,
        dry_run=args.dry_run,
        show_progress=not args.no_progress,
        min_len=args.min_len,
        max_len=args.max_len,
        max_ambiguous_ratio=args.max_ambiguous_ratio,
        dedupe_sequence=args.dedupe_sequence,
        sample=args.sample,
        molecule_type=args.molecule_type,
        organism_include=org_include_pats,
        organism_exclude=org_exclude_pats,
        taxid_include=taxid_include_set,
        taxid_exclude=taxid_exclude_set,
        lineage_include=lin_include_pats,
        lineage_exclude=lin_exclude_pats,
        streaming=args.streaming,
        max_go=args.max_go,
        max_ec=args.max_ec,
        go_category=go_category,
    )

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

        file_results = process_file(input_path, output_path, cfg)
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
