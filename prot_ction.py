#!/usr/bin/env python3
"""
Protein sequence conservation calculator.

Similarity metric: normalized inverted Hamming distance
    similarity = 1 - (mismatches / comparable_positions)

Modes (--mode):
  auto      Auto-detect: equal-length → MSA column comparison;
            unequal-length → sliding approach.  (default)
  aligned   Input is a pre-computed MSA.  All sequences must be the same
            length.  Gap characters (- .) are respected as alignment columns
            and skipped when computing similarity.  Per-position conservation
            is derived directly from the alignment columns.
  unaligned Treat all sequences as raw/ungapped, regardless of length.
            Any gap characters are stripped before processing.  The sliding
            approach is always used.

Taxonomy (--rank): fetches full taxonomic lineage from UniProt / NCBI for
    each accession, groups sequences by the chosen rank (e.g. genus, family),
    and reports the median pairwise similarity within and between groups.
    Sequences for which taxonomy cannot be retrieved are automatically dropped.

Input: FASTA file with UniProt (sp|P12345|GENE_HUMAN) or GenBank (NP_000001.1)
       identifiers in the definition line.
"""
from __future__ import annotations

VERSION = "1.0.0"

import re
import sys
import csv
import json
import time
import pathlib
import argparse
import statistics
import urllib.request
import urllib.error
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations
from xml.etree import ElementTree as ET


# ---------------------------------------------------------------------------
# Identifier parsing
# ---------------------------------------------------------------------------

def parse_identifier(header: str) -> tuple[str, str]:
    """
    Extract accession and database type from a FASTA header line.

    Supported formats:
      UniProt  >sp|P12345|GENE_HUMAN ...  or  >tr|A0A000|GENE_SPECIES ...
      GenBank  >NP_000001.1 ...  /  >XP_... / >gi|...|ref|NP_...|
      Plain    >MYPROT  (falls back to first whitespace-delimited token)

    Returns (accession, db_type) where db_type is 'uniprot', 'genbank', or 'unknown'.
    """
    header = header.lstrip(">").strip()

    # UniProt Swiss-Prot / TrEMBL
    m = re.match(r"(?:sp|tr)\|([A-Z0-9]+)\|\S+", header)
    if m:
        return m.group(1), "uniprot"

    # NCBI gi pipe format  gi|12345|ref|NP_000001.1|
    m = re.search(r"gi\|\d+\|\w+\|([A-Z_0-9.]+)\|?", header)
    if m:
        return m.group(1), "genbank"

    # Bare GenBank/RefSeq accession (e.g. NP_000001.1, XP_001234567.2, P12345)
    m = re.match(r"([A-Z]{1,3}[_]?\d+\.?\d*)", header)
    if m:
        return m.group(1), "genbank"

    # Fall back: first token
    return header.split()[0], "unknown"


# ---------------------------------------------------------------------------
# FASTA parser
# ---------------------------------------------------------------------------

def parse_fasta(path: str) -> list[tuple[str, str, str]]:
    """
    Parse a FASTA file.

    Returns list of (accession, db_type, sequence) tuples.
    Sequences are uppercased; ambiguous/gap characters are preserved.
    """
    records = []
    acc = db = None
    buf: list[str] = []

    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if acc is not None:
                    records.append((acc, db, "".join(buf).upper()))
                acc, db = parse_identifier(line)
                buf = []
            else:
                buf.append(line)

    if acc is not None:
        records.append((acc, db, "".join(buf).upper()))

    return records


# ---------------------------------------------------------------------------
# Taxonomy – data structures
# ---------------------------------------------------------------------------

NCBI_BASE    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
# NCBI rate limit: ≤3 req/s without key, ≤10 req/s with key
_NCBI_DELAY     = 0.34
_NCBI_DELAY_KEY = 0.11


@dataclass
class TaxNode:
    taxon_id: int
    name: str
    rank: str   # lowercase NCBI rank string, e.g. 'species', 'genus', 'no rank'


@dataclass
class TaxInfo:
    taxon_id: int
    scientific_name: str
    lineage: list[TaxNode] = field(default_factory=list)
    # lineage is ordered root → organism (organism itself is the last entry)

    def at_rank(self, rank: str) -> TaxNode | None:
        """Return the lineage node whose rank matches, or None."""
        rank = rank.strip().lower()
        for node in self.lineage:
            if node.rank == rank:
                return node
        return None

    def available_ranks(self) -> list[str]:
        return sorted({n.rank for n in self.lineage if n.rank != "no rank"})


# ---------------------------------------------------------------------------
# Taxonomy – persistent cache
# ---------------------------------------------------------------------------

_CACHE_VERSION = 1
_DEFAULT_CACHE  = "~/.protein_conservation_cache.json"


class TaxCache:
    """
    Persistent JSON cache for taxonomy lookups.

    Two independent lookup tables are stored:

      taxids   – protein accession (str) → NCBI taxon ID (int)
      taxonomy – NCBI taxon ID (int)     → TaxInfo (full lineage)

    Because a taxon ID is shared by many proteins, caching the taxonomy entry
    keyed on the taxon ID avoids redundant NCBI Taxonomy calls for proteins
    from the same organism.

    The cache is loaded from disk in ``__init__`` and written back only when
    ``save()`` is called and new entries have been added (dirty flag).  A safe
    write-then-rename strategy prevents corruption on interrupted runs.
    """

    def __init__(self, path: str) -> None:
        self._path  = pathlib.Path(path).expanduser()
        self._taxids:   dict[str, int]     = {}
        self._taxonomy: dict[int, TaxInfo] = {}
        self._dirty = False
        self._load()

    # ------------------------------------------------------------------
    # Serialisation helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _info_to_dict(info: TaxInfo) -> dict:
        return {
            "taxon_id":       info.taxon_id,
            "scientific_name": info.scientific_name,
            "lineage": [
                {"taxon_id": n.taxon_id, "name": n.name, "rank": n.rank}
                for n in info.lineage
            ],
        }

    @staticmethod
    def _dict_to_info(d: dict) -> TaxInfo:
        return TaxInfo(
            taxon_id=d["taxon_id"],
            scientific_name=d["scientific_name"],
            lineage=[
                TaxNode(taxon_id=n["taxon_id"], name=n["name"], rank=n["rank"])
                for n in d.get("lineage", [])
            ],
        )

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def _load(self) -> None:
        if not self._path.exists():
            return
        try:
            with self._path.open() as fh:
                data = json.load(fh)
            if data.get("version") != _CACHE_VERSION:
                return          # incompatible format – start fresh
            self._taxids = {k: int(v) for k, v in data.get("taxids", {}).items()}
            self._taxonomy = {
                int(k): self._dict_to_info(v)
                for k, v in data.get("taxonomy", {}).items()
            }
        except Exception:
            pass                # corrupt / unreadable – start fresh silently

    def save(self) -> None:
        """Flush new entries to disk; no-op when nothing changed."""
        if not self._dirty:
            return
        payload = {
            "version":  _CACHE_VERSION,
            "taxids":   self._taxids,
            "taxonomy": {
                str(k): self._info_to_dict(v)
                for k, v in self._taxonomy.items()
            },
        }
        try:
            self._path.parent.mkdir(parents=True, exist_ok=True)
            tmp = self._path.with_suffix(".tmp")
            with tmp.open("w") as fh:
                json.dump(payload, fh, indent=2)
            tmp.replace(self._path)
            self._dirty = False
        except Exception as exc:
            print(f"  Warning: could not write cache to {self._path}: {exc}",
                  file=sys.stderr)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_taxid(self, accession: str) -> int | None:
        return self._taxids.get(accession)

    def set_taxid(self, accession: str, taxid: int) -> None:
        if self._taxids.get(accession) != taxid:
            self._taxids[accession] = taxid
            self._dirty = True

    def get_taxonomy(self, taxid: int) -> TaxInfo | None:
        return self._taxonomy.get(taxid)

    def set_taxonomy(self, taxid: int, info: TaxInfo) -> None:
        if taxid not in self._taxonomy:
            self._taxonomy[taxid] = info
            self._dirty = True

    def stats(self) -> tuple[int, int]:
        """Return (n_accessions_cached, n_taxa_cached)."""
        return len(self._taxids), len(self._taxonomy)


# ---------------------------------------------------------------------------
# Taxonomy – API helpers
# ---------------------------------------------------------------------------

def _http_get(url: str, timeout: int = 20) -> str:
    req = urllib.request.Request(
        url, headers={"User-Agent": "protein_conservation/2.0 (github)"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def _ncbi_delay(api_key: str | None) -> None:
    time.sleep(_NCBI_DELAY_KEY if api_key else _NCBI_DELAY)


def _taxonomy_from_ncbi(taxid: int, api_key: str | None = None) -> TaxInfo | None:
    """
    Fetch full taxonomic lineage from the NCBI Taxonomy database.
    Returns a TaxInfo with the complete LineageEx plus the organism itself.
    """
    url = (f"{NCBI_BASE}/efetch.fcgi?db=taxonomy&id={taxid}&retmode=xml"
           + (f"&api_key={api_key}" if api_key else ""))
    _ncbi_delay(api_key)
    try:
        body = _http_get(url)
    except Exception:
        return None

    try:
        root = ET.fromstring(body)
    except ET.ParseError:
        return None

    taxon_el = root.find("Taxon")
    if taxon_el is None:
        return None

    sci_name = taxon_el.findtext("ScientificName", "").strip()
    rank     = taxon_el.findtext("Rank", "no rank").strip().lower()

    lineage: list[TaxNode] = []
    lineage_ex = taxon_el.find("LineageEx")
    if lineage_ex is not None:
        for node in lineage_ex.findall("Taxon"):
            lineage.append(TaxNode(
                taxon_id=int(node.findtext("TaxId", "0")),
                name=node.findtext("ScientificName", "").strip(),
                rank=node.findtext("Rank", "no rank").strip().lower(),
            ))
    lineage.append(TaxNode(taxon_id=taxid, name=sci_name, rank=rank))

    return TaxInfo(taxon_id=taxid, scientific_name=sci_name, lineage=lineage)


def _taxid_from_uniprot(accession: str) -> int | None:
    """
    Look up the NCBI taxon ID for a UniProt accession via the UniProt REST API.
    """
    url = f"{UNIPROT_BASE}/{accession}.json"
    try:
        body = _http_get(url)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise
    data = json.loads(body)
    return data.get("organism", {}).get("taxonId")


def _taxid_from_ncbi_protein(accession: str, api_key: str | None = None) -> int | None:
    """
    Look up the NCBI taxon ID for a GenBank/RefSeq protein accession.
    Parses the taxon:XXXXXX qualifier from the GenPept source feature.
    """
    url = (f"{NCBI_BASE}/efetch.fcgi?db=protein&id={accession}"
           f"&rettype=gp&retmode=xml"
           + (f"&api_key={api_key}" if api_key else ""))
    _ncbi_delay(api_key)
    try:
        body = _http_get(url)
    except Exception:
        return None

    try:
        root = ET.fromstring(body)
    except ET.ParseError:
        return None

    for qual in root.iter("GBQualifier"):
        if qual.findtext("GBQualifier_name") == "db_xref":
            val = qual.findtext("GBQualifier_value", "")
            if val.startswith("taxon:"):
                return int(val.split(":")[1])
    return None


def fetch_taxonomy(accession: str, db_type: str,
                   api_key: str | None = None,
                   cache: TaxCache | None = None) -> tuple[TaxInfo | None, bool]:
    """
    Fetch full taxonomy for one protein accession.

    Workflow
    --------
    UniProt  → UniProt REST API (taxon ID)  → NCBI Taxonomy (lineage)
    GenBank  → NCBI EUtils protein fetch    → NCBI Taxonomy (lineage)

    Both the accession→taxid mapping and the taxid→TaxInfo record are
    checked in *cache* before making any network request, and stored there
    after a successful fetch.

    Returns (TaxInfo | None, from_cache: bool).
    """
    try:
        # ── Step 1: resolve accession → taxon ID ──────────────────────────
        taxid: int | None = cache.get_taxid(accession) if cache else None
        taxid_from_cache  = taxid is not None

        if taxid is None:
            if db_type == "uniprot":
                taxid = _taxid_from_uniprot(accession)
            else:
                taxid = _taxid_from_ncbi_protein(accession, api_key)
            if taxid is None:
                return None, False
            if cache:
                cache.set_taxid(accession, taxid)

        # ── Step 2: resolve taxon ID → full lineage ────────────────────────
        info: TaxInfo | None = cache.get_taxonomy(taxid) if cache else None
        info_from_cache = info is not None

        if info is None:
            info = _taxonomy_from_ncbi(taxid, api_key)
            if info is None:
                return None, False
            if cache:
                cache.set_taxonomy(taxid, info)

        from_cache = taxid_from_cache and info_from_cache
        return info, from_cache

    except Exception:
        return None, False


def fetch_all_taxonomies(
    records: list[tuple[str, str, str]],
    api_key: str | None = None,
    cache: TaxCache | None = None,
) -> dict[str, TaxInfo]:
    """
    Fetch taxonomy for every record in *records*, printing progress.

    Cache hits are served instantly without any network call and are
    labelled ``[cached]`` in the progress output.  After all records are
    processed the cache is flushed to disk (if it is dirty).

    Returns a dict mapping accession → TaxInfo.
    Accessions for which taxonomy cannot be retrieved are omitted.
    """
    if cache:
        n_acc, n_tax = cache.stats()
        if n_acc or n_tax:
            print(f"  Cache: {n_acc} accession(s), {n_tax} taxon record(s) "
                  f"loaded from {cache._path}")

    tax_map: dict[str, TaxInfo] = {}
    n = len(records)
    n_hits = 0
    for idx, (acc, db, _seq) in enumerate(records, 1):
        print(f"  [{idx}/{n}] {acc} ({db})… ", end="", flush=True)
        info, from_cache = fetch_taxonomy(acc, db, api_key, cache)
        if info is None:
            print("not found – skipped")
        else:
            label = " [cached]" if from_cache else ""
            print(f"{info.scientific_name}{label}")
            tax_map[acc] = info
            if from_cache:
                n_hits += 1

    if cache:
        cache.save()
        if n_hits:
            print(f"  {n_hits}/{n} record(s) served from cache "
                  f"(no network call needed).")

    return tax_map


# ---------------------------------------------------------------------------
# Taxonomy – rank-based conservation
# ---------------------------------------------------------------------------

def group_by_rank(
    ids: list[str],
    tax_map: dict[str, TaxInfo],
    rank: str,
) -> dict[str, list[int]]:
    """
    Return {taxon_name: [sequence_indices]} grouped by the given rank.
    Sequences without taxonomy or without that rank are omitted.
    """
    groups: dict[str, list[int]] = {}
    for i, acc in enumerate(ids):
        info = tax_map.get(acc)
        if info is None:
            continue
        node = info.at_rank(rank)
        if node is None:
            continue
        groups.setdefault(node.name, []).append(i)
    return groups


def rank_conservation(
    ids: list[str],
    seqs: list[str],
    sim_mat: list[list[float]],
    tax_map: dict[str, TaxInfo],
    rank: str,
) -> dict:
    """
    Compute median pairwise similarity within and between taxa at *rank*.

    Returns a dict with keys:
      'rank'          – the requested rank string
      'groups'        – {taxon_name: {'members': [accessions],
                                      'intra_pairs': [(sim, id1, id2)],
                                      'intra_median': float | None}}
      'inter_pairs'   – [(sim, taxon_a, taxon_b, id1, id2)]
      'intra_median'  – overall median across all within-group pairs
      'inter_median'  – overall median across all between-group pairs
    """
    groups = group_by_rank(ids, tax_map, rank)

    group_data: dict[str, dict] = {}
    for taxon, indices in groups.items():
        intra = [
            (sim_mat[i][j], ids[i], ids[j])
            for i, j in combinations(indices, 2)
        ]
        group_data[taxon] = {
            "members": [ids[i] for i in indices],
            "intra_pairs": intra,
            "intra_median": statistics.median(s for s, *_ in intra) if intra else None,
        }

    taxon_names = list(groups.keys())
    inter_pairs: list[tuple] = []
    for ta, tb in combinations(taxon_names, 2):
        for i in groups[ta]:
            for j in groups[tb]:
                inter_pairs.append((sim_mat[i][j], ta, tb, ids[i], ids[j]))

    all_intra = [s for g in group_data.values() for s, *_ in g["intra_pairs"]]
    all_inter = [s for s, *_ in inter_pairs]

    return {
        "rank":         rank,
        "groups":       group_data,
        "inter_pairs":  inter_pairs,
        "intra_median": statistics.median(all_intra) if all_intra else None,
        "inter_median": statistics.median(all_inter) if all_inter else None,
    }


def permutation_test_rank(
    ids: list[str],
    sim_mat: list[list[float]],
    tax_map: dict[str, "TaxInfo"],
    rank: str,
    n_permutations: int = 1000,
    seed: "int | None" = None,
) -> dict:
    """
    Permutation test for rank-based conservation significance.

    Observed test statistic: T = intra_median − inter_median.
    Null distribution: shuffle group labels (preserving group sizes) and
    recompute T for each permutation.

    Returns a dict with keys:
      'observed_stat'  – observed T value (float | None)
      'p_value'        – fraction of permuted T ≥ observed T (float | None)
      'n_permutations' – number of permutations actually run
      'perm_mean'      – mean of permuted T values (float | None)
      'perm_stdev'     – stdev of permuted T values (float | None)
      'reason'         – non-empty string if test could not run
    """
    import random

    result = rank_conservation(ids, [], sim_mat, tax_map, rank)
    obs_intra = result["intra_median"]
    obs_inter = result["inter_median"]

    if obs_intra is None or obs_inter is None:
        return {
            "observed_stat": None, "p_value": None,
            "n_permutations": 0, "perm_mean": None, "perm_stdev": None,
            "reason": "cannot compute test statistic (need ≥2 groups with ≥2 members each)",
        }

    observed_stat = obs_intra - obs_inter

    # Build the assignment vector: list of group labels parallel to ids
    groups = result["groups"]
    label_vector: list[str] = [""] * len(ids)
    id_to_idx = {acc: i for i, acc in enumerate(ids)}
    for taxon, gdata in groups.items():
        for acc in gdata["members"]:
            if acc in id_to_idx:
                label_vector[id_to_idx[acc]] = taxon

    # Only keep indices that were assigned to a group
    valid_indices = [i for i, lbl in enumerate(label_vector) if lbl]
    if len(valid_indices) < 4:
        return {
            "observed_stat": observed_stat, "p_value": None,
            "n_permutations": 0, "perm_mean": None, "perm_stdev": None,
            "reason": "too few sequences for permutation test",
        }

    # Preserve group sizes
    group_sizes = [len(gdata["members"]) for gdata in groups.values()]
    group_labels = list(groups.keys())

    rng = random.Random(seed)
    perm_stats: list[float] = []

    for _ in range(n_permutations):
        # Shuffle valid_indices, then re-assign group labels by size
        shuffled = valid_indices[:]
        rng.shuffle(shuffled)
        perm_assignment: dict[int, str] = {}
        pos = 0
        for taxon, size in zip(group_labels, group_sizes):
            for idx in shuffled[pos:pos + size]:
                perm_assignment[idx] = taxon
            pos += size

        # Compute intra and inter medians for this permutation
        intra_sims: list[float] = []
        inter_sims: list[float] = []
        for ii, jj in combinations(valid_indices, 2):
            ta = perm_assignment.get(ii, "")
            tb = perm_assignment.get(jj, "")
            if not ta or not tb:
                continue
            if ta == tb:
                intra_sims.append(sim_mat[ii][jj])
            else:
                inter_sims.append(sim_mat[ii][jj])

        if not intra_sims or not inter_sims:
            continue
        perm_stats.append(
            statistics.median(intra_sims) - statistics.median(inter_sims)
        )

    if not perm_stats:
        return {
            "observed_stat": observed_stat, "p_value": None,
            "n_permutations": 0, "perm_mean": None, "perm_stdev": None,
            "reason": "no valid permutations produced",
        }

    p_value = sum(1 for t in perm_stats if t >= observed_stat) / len(perm_stats)
    perm_mean  = statistics.mean(perm_stats)
    perm_stdev = statistics.stdev(perm_stats) if len(perm_stats) > 1 else 0.0

    return {
        "observed_stat": observed_stat,
        "p_value":        p_value,
        "n_permutations": len(perm_stats),
        "perm_mean":      perm_mean,
        "perm_stdev":     perm_stdev,
        "reason":         "",
    }


def print_rank_conservation(result: dict, ids: list[str],
                             tax_map: dict[str, TaxInfo],
                             perm_result: "dict | None" = None) -> None:
    rank = result["rank"]
    groups = result["groups"]

    print(f"\nTaxonomic conservation at rank: {rank}")
    print("=" * 60)

    # Sequence → taxon mapping table
    print(f"\n  {'Accession':<20} {'Organism':<30} {rank.capitalize()}")
    print(f"  {'-'*20} {'-'*30} {'-'*20}")
    for acc in ids:
        info = tax_map.get(acc)
        if info is None:
            print(f"  {acc:<20} {'(no taxonomy)':<30} –")
            continue
        node = info.at_rank(rank)
        rank_name = node.name if node else f"(no {rank} rank)"
        print(f"  {acc:<20} {info.scientific_name:<30} {rank_name}")

    # Per-group within-rank median
    print(f"\n  Within-{rank} median similarity")
    print(f"  {'Taxon':<30} {'Members':>7} {'Pairs':>6} {'Median':>8}")
    print(f"  {'-'*30} {'-'*7} {'-'*6} {'-'*8}")
    for taxon, gdata in sorted(groups.items()):
        n_mem   = len(gdata["members"])
        n_pairs = len(gdata["intra_pairs"])
        med     = gdata["intra_median"]
        med_str = f"{med:.4f}" if med is not None else "  n/a"
        print(f"  {taxon:<30} {n_mem:>7} {n_pairs:>6} {med_str:>8}")

    # Between-group summary
    inter_pairs = result["inter_pairs"]
    if inter_pairs:
        print(f"\n  Between-{rank} median similarity")
        print(f"  {'Taxon A':<25} {'Taxon B':<25} {'Pairs':>6} {'Median':>8}")
        print(f"  {'-'*25} {'-'*25} {'-'*6} {'-'*8}")
        # Collect per-(ta,tb) medians
        inter_by_pair: dict[tuple, list[float]] = defaultdict(list)
        for sim, ta, tb, _id1, _id2 in inter_pairs:
            inter_by_pair[(ta, tb)].append(sim)
        for (ta, tb), sims in sorted(inter_by_pair.items()):
            med = statistics.median(sims)
            print(f"  {ta:<25} {tb:<25} {len(sims):>6} {med:>8.4f}")

    # Overall summary
    print(f"\n  Overall summary")
    intra_med = result["intra_median"]
    inter_med = result["inter_median"]
    print(f"  Intra-{rank} median:   "
          f"{intra_med:.4f}" if intra_med is not None else f"  Intra-{rank} median:   n/a")
    print(f"  Inter-{rank} median:   "
          f"{inter_med:.4f}" if inter_med is not None else f"  Inter-{rank} median:   n/a")
    if intra_med is not None and inter_med is not None:
        direction = "higher" if intra_med > inter_med else "lower"
        print(f"  Within-{rank} similarity is {direction} than between-{rank} "
              f"(Δ = {abs(intra_med - inter_med):.4f})")

    # Permutation test
    if perm_result is not None:
        print(f"\n  Permutation test (n={perm_result['n_permutations']})")
        reason = perm_result.get("reason", "")
        if reason:
            print(f"  Could not run test: {reason}")
        else:
            obs  = perm_result["observed_stat"]
            pval = perm_result["p_value"]
            pmn  = perm_result["perm_mean"]
            psd  = perm_result["perm_stdev"]
            print(f"  Observed T (intra − inter):  {obs:+.4f}")
            print(f"  Null mean ± stdev:           {pmn:+.4f} ± {psd:.4f}")
            sig = " *" if pval is not None and pval < 0.05 else ""
            pval_str = f"{pval:.4f}{sig}" if pval is not None else "n/a"
            print(f"  p-value (one-tailed):        {pval_str}")


# ---------------------------------------------------------------------------
# Taxonomy – heatmap
# ---------------------------------------------------------------------------

def _build_heatmap_matrix(
    result: dict,
) -> tuple[list[str], list[list[float]]]:
    """
    Build a symmetric N×N matrix of median similarities indexed by taxon name.

    Diagonal cell  [i][i] – intra-group median (NaN when the group has only
                            one member and therefore no pairs).
    Off-diagonal   [i][j] – inter-group median between group i and group j.
    """
    names = sorted(result["groups"].keys())
    N     = len(names)
    idx   = {n: i for i, n in enumerate(names)}
    nan   = float("nan")
    mat   = [[nan] * N for _ in range(N)]

    for name, gdata in result["groups"].items():
        i = idx[name]
        if gdata["intra_median"] is not None:
            mat[i][i] = gdata["intra_median"]

    inter_sims: dict[tuple[int, int], list[float]] = defaultdict(list)
    for sim, ta, tb, *_ in result["inter_pairs"]:
        i, j = idx[ta], idx[tb]
        inter_sims[(i, j)].append(sim)

    for (i, j), sims in inter_sims.items():
        med = statistics.median(sims)
        mat[i][j] = mat[j][i] = med

    return names, mat


def plot_rank_heatmap(result: dict, path: str) -> None:
    """
    Save a colour-coded N×N heatmap of median pairwise similarity
    between (and within) taxon groups at the selected rank.

    Requires matplotlib (``pip install matplotlib``).  The output format
    is inferred from the file extension (.png, .pdf, .svg, …).

    Layout
    ------
    * Rows and columns are taxon names, sorted alphabetically.
    * Each cell is annotated with the numeric median (3 decimal places).
    * Diagonal cells are framed with a black border to distinguish
      intra-group similarity from between-group similarity.
    * Cells with no pairs (single-member groups on the diagonal) are
      rendered in neutral grey and labelled "n/a".
    * A shared colour scale runs from 0 (red) through yellow to 1 (green),
      matching the direction "higher = more conserved".
    """
    try:
        import matplotlib
        matplotlib.use("Agg")               # file-only backend; safe for CLI
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import numpy as np
    except ImportError:
        print(
            "  Warning: matplotlib is not installed – heatmap skipped.\n"
            "  Install it with:  pip install matplotlib",
            file=sys.stderr,
        )
        return

    rank  = result["rank"]
    names, mat_lists = _build_heatmap_matrix(result)
    N     = len(names)

    if N < 2:
        print("  Note: heatmap requires at least 2 taxon groups – skipped.")
        return

    mat = np.array(mat_lists, dtype=float)   # NaN preserved

    # ── figure size scales with number of groups ──────────────────────────
    cell_px   = max(1.0, min(1.8, 10.0 / N))
    fig_w     = max(5.0, N * cell_px + 3.0)
    fig_h     = max(4.0, N * cell_px + 2.0)
    fig, ax   = plt.subplots(figsize=(fig_w, fig_h))

    cmap = plt.get_cmap("RdYlGn").copy()
    cmap.set_bad(color="#d3d3d3")            # NaN cells → light grey

    im = ax.imshow(mat, cmap=cmap, vmin=0.0, vmax=1.0, aspect="auto")

    # ── colour bar ────────────────────────────────────────────────────────
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Median pairwise similarity", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # ── cell text ─────────────────────────────────────────────────────────
    fs = max(6, min(10, int(cell_px * 5.5)))
    for i in range(N):
        for j in range(N):
            val = mat[i, j]
            if np.isnan(val):
                txt      = "n/a"
                txt_col  = "#888888"
                bold     = False
            else:
                txt      = f"{val:.3f}"
                # White text on dark ends of RdYlGn; black on the light middle
                txt_col  = "white" if val < 0.22 or val > 0.82 else "black"
                bold     = (i == j)
            ax.text(
                j, i, txt,
                ha="center", va="center",
                fontsize=fs,
                color=txt_col,
                fontweight="bold" if bold else "normal",
            )

    # ── thick border on diagonal (intra-group) cells ──────────────────────
    for k in range(N):
        rect = plt.Rectangle(
            (k - 0.5, k - 0.5), 1, 1,
            linewidth=2.2, edgecolor="black", facecolor="none",
        )
        ax.add_patch(rect)

    # ── tick labels ───────────────────────────────────────────────────────
    label_fs = max(7, min(10, int(cell_px * 5.5)))
    ax.set_xticks(range(N))
    ax.set_yticks(range(N))
    ax.set_xticklabels(names, rotation=40, ha="right", fontsize=label_fs,
                       style="italic")
    ax.set_yticklabels(names, fontsize=label_fs, style="italic")

    # ── title ─────────────────────────────────────────────────────────────
    ax.set_title(
        f"Median pairwise similarity  ·  grouped by {rank}\n"
        f"diagonal = within-{rank}   ·   off-diagonal = between-{rank}",
        fontsize=10, pad=14,
    )

    # ── legend for n/a cells ─────────────────────────────────────────────
    na_patch = mpatches.Patch(
        facecolor="#d3d3d3", edgecolor="#999999",
        label="n/a  (single-member group)",
    )
    ax.legend(
        handles=[na_patch],
        loc="upper left", bbox_to_anchor=(1.18, 1.02),
        fontsize=8, framealpha=0.85,
    )

    plt.tight_layout()
    try:
        plt.savefig(path, dpi=150, bbox_inches="tight")
    except ValueError as exc:
        print(f"  Error saving heatmap: {exc}", file=sys.stderr)
        plt.close(fig)
        return
    plt.close(fig)
    print(f"  Heatmap written to: {path}")


def plot_dendrogram(
    ids: list[str],
    sim_mat: list[list[float]],
    path: str,
    method: str = "average",
) -> None:
    """
    Save a hierarchical-clustering dendrogram of the sequences.

    Distance = 1 − similarity.  Linkage method is *average* (UPGMA) by
    default; any method accepted by ``scipy.cluster.hierarchy.linkage``
    is valid (single, complete, ward, …).

    Requires matplotlib and scipy (``pip install matplotlib scipy``).
    The output format is inferred from the file extension.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "  Warning: matplotlib is not installed – dendrogram skipped.\n"
            "  Install it with:  pip install matplotlib scipy",
            file=sys.stderr,
        )
        return
    try:
        from scipy.spatial.distance import squareform
        from scipy.cluster.hierarchy import linkage, dendrogram
    except ImportError:
        print(
            "  Warning: scipy is not installed – dendrogram skipped.\n"
            "  Install it with:  pip install scipy",
            file=sys.stderr,
        )
        return

    n = len(ids)
    if n < 2:
        print("  Note: dendrogram requires at least 2 sequences – skipped.")
        return

    # Build condensed distance matrix (upper triangle, row-major)
    import numpy as np
    dist_sq = np.array([[1.0 - sim_mat[i][j] for j in range(n)]
                        for i in range(n)], dtype=float)
    np.fill_diagonal(dist_sq, 0.0)
    # clip any floating-point negatives
    dist_sq = np.clip(dist_sq, 0.0, None)
    condensed = squareform(dist_sq, checks=False)

    Z = linkage(condensed, method=method)

    fig_h = max(4.0, 1.8 + n * 0.35)
    fig, ax = plt.subplots(figsize=(8, fig_h))

    dendrogram(
        Z,
        labels=ids,
        orientation="left",
        ax=ax,
        color_threshold=0.7 * max(Z[:, 2]) if len(Z) else 0,
        leaf_font_size=max(7, min(11, int(120 / n))),
    )

    ax.set_xlabel("Distance  (1 − similarity)", fontsize=9)
    ax.set_title(
        f"Hierarchical clustering  ·  {method} linkage  ·  "
        f"{n} sequences",
        fontsize=10, pad=10,
    )
    ax.tick_params(axis="y", labelsize=max(7, min(10, int(110 / n))))
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    try:
        plt.savefig(path, dpi=150, bbox_inches="tight")
    except ValueError as exc:
        print(f"  Error saving dendrogram: {exc}", file=sys.stderr)
        plt.close(fig)
        return
    plt.close(fig)
    print(f"  Dendrogram written to: {path}")


def export_phyloxml(
    ids: list[str],
    sim_mat: list[list[float]],
    path: str,
    method: str = "average",
) -> None:
    """
    Export the hierarchical clustering tree as a PhyloXML file.

    Uses the same distance matrix (1 − similarity) and linkage method as
    ``plot_dendrogram``.  Branch lengths represent the difference in
    clustering height between a node and its parent (i.e. the additional
    distance accumulated at each merge).

    Requires scipy (``pip install scipy``).  The output is a valid
    PhyloXML 1.10 document readable by tools such as Archaeopteryx,
    FigTree, Biopython, and iTOL.
    """
    try:
        from scipy.spatial.distance import squareform
        from scipy.cluster.hierarchy import linkage
        import numpy as np
    except ImportError:
        print(
            "  Warning: scipy is not installed – PhyloXML export skipped.\n"
            "  Install it with:  pip install scipy",
            file=sys.stderr,
        )
        return

    n = len(ids)
    if n < 2:
        print("  Note: PhyloXML export requires at least 2 sequences – skipped.")
        return

    dist_sq = np.array([[1.0 - sim_mat[i][j] for j in range(n)]
                        for i in range(n)], dtype=float)
    np.fill_diagonal(dist_sq, 0.0)
    dist_sq = np.clip(dist_sq, 0.0, None)
    Z = linkage(squareform(dist_sq, checks=False), method=method)

    # ── recursively build <clade> elements ──────────────────────────────────
    def _clade(node_id: int, parent_height: float) -> ET.Element:
        elem = ET.Element("clade")
        if node_id < n:
            bl = parent_height          # leaf height is 0
            ET.SubElement(elem, "name").text = ids[node_id]
        else:
            row    = Z[node_id - n]
            height = float(row[2])
            bl     = parent_height - height
            elem.append(_clade(int(row[0]), height))
            elem.append(_clade(int(row[1]), height))
        if bl > 1e-10:
            ET.SubElement(elem, "branch_length").text = f"{bl:.8f}"
        return elem

    root_id     = n + len(Z) - 1
    root_height = float(Z[-1][2])

    # ── assemble PhyloXML document ───────────────────────────────────────────
    phyloxml = ET.Element("phyloxml")
    phyloxml.set("xmlns",              "http://www.phyloxml.org")
    phyloxml.set("xmlns:xsi",          "http://www.w3.org/2001/XMLSchema-instance")
    phyloxml.set("xsi:schemaLocation", (
        "http://www.phyloxml.org "
        "http://www.phyloxml.org/1.10/phyloxml.xsd"
    ))

    phylogeny = ET.SubElement(phyloxml, "phylogeny")
    phylogeny.set("rooted", "true")
    ET.SubElement(phylogeny, "name").text = (
        f"Hierarchical clustering ({method} linkage)"
    )
    ET.SubElement(phylogeny, "description").text = (
        "Branch lengths are clustering-height differences (distance = 1 − similarity)."
    )
    phylogeny.append(_clade(root_id, root_height))

    # indent if Python ≥ 3.9, otherwise write flat
    try:
        ET.indent(phyloxml, space="  ")
    except AttributeError:
        pass

    out = pathlib.Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write(ET.tostring(phyloxml, encoding="unicode"))
        fh.write("\n")
    print(f"  PhyloXML written to: {path}")


# ---------------------------------------------------------------------------
# Scoring matrices
# ---------------------------------------------------------------------------

def _parse_matrix_table(raw: str) -> dict[str, dict[str, int]]:
    """
    Parse a whitespace-delimited matrix table (header row + data rows).
    The first token on each data row is the row amino-acid label.
    """
    lines  = [l for l in raw.strip().splitlines() if l.strip()]
    header = lines[0].split()
    mat: dict[str, dict[str, int]] = {}
    for row in lines[1:]:
        parts = row.split()
        aa    = parts[0]
        mat[aa] = {col: int(val) for col, val in zip(header, parts[1:])}
    return mat


# BLOSUM62 – Henikoff & Henikoff (1992), NCBI standard values
_BLOSUM62 = _parse_matrix_table("""\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4""")

# PAM250 – Dayhoff, Schwartz & Orcutt (1978)
_PAM250 = _parse_matrix_table("""\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  2 -2  0  0 -2  0  0  1 -1 -1 -2 -1 -1 -3  1  1  1 -6 -3  0
R -2  6  0 -1 -3  1 -1 -3  2 -2 -3  3  0 -4  0  0 -1  2 -4 -2
N  0  0  2  2 -4  1  1  0  2 -2 -3  1 -2 -3  0  1  0 -4 -2 -2
D  0 -1  2  4 -5  2  3  1  1 -2 -4  0 -3 -6 -1  0  0 -7 -4 -2
C -2 -3 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3  0 -2 -8  0 -2
Q  0  1  1  2 -5  4  2 -1  3 -2 -2  1 -1 -5  0 -1 -1 -5 -4 -2
E  0 -1  1  3 -5  2  4  0  1 -2 -3  0 -2 -5 -1  0  0 -7 -4 -2
G  1 -3  0  1 -3 -1  0  5 -2 -3 -4 -2 -3 -5  0  1  0 -7 -5 -1
H -1  2  2  1 -3  3  1 -2  6 -2 -2  0 -2 -2  0 -1 -1 -3  0 -2
I -1 -2 -2 -2 -2 -2 -2 -3 -2  5  2 -2  2  1 -2 -1  0 -5 -1  4
L -2 -3 -3 -4 -6 -2 -3 -4 -2  2  6 -3  4  2 -3 -3 -2 -2 -1  2
K -1  3  1  0 -5  1  0 -2  0 -2 -3  5  0 -5 -1  0  0 -3 -4 -2
M -1  0 -2 -3 -5 -1 -2 -3 -2  2  4  0  6  0 -2 -2 -1 -4 -2  2
F -3 -4 -3 -6 -4 -5 -5 -5 -2  1  2 -5  0  9 -5 -3 -3  0  7 -1
P  1  0  0 -1 -3  0 -1  0  0 -2 -3 -1 -2 -5  6  1  0 -6 -5 -1
S  1  0  1  0  0 -1  0  1 -1 -1 -3  0 -2 -3  1  2  1 -2 -3 -1
T  1 -1  0  0 -2 -1  0  0 -1  0 -2  0 -1 -3  0  1  3 -5 -3  0
W -6  2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4  0 -6 -2 -5 17  0 -6
Y -3 -4 -2 -4  0 -4 -4 -5  0 -1 -1 -4 -2  7 -5 -3 -3  0 10 -2
V  0 -2 -2 -2 -2 -2 -2 -1 -2  4  2 -2  2 -1 -1 -1  0 -6 -2  4""")

_BUILTIN_MATRICES: dict[str, dict[str, dict[str, int]]] = {
    "blosum62": _BLOSUM62,
    "pam250":   _PAM250,
}

# All names accepted by --scoring; others need Biopython installed.
SCORING_CHOICES = [
    "hamming",
    "blosum45", "blosum62", "blosum80", "blosum90",
    "pam30", "pam70", "pam120", "pam250",
]


def load_matrix(name: str) -> dict[str, dict[str, int]] | None:
    """
    Return a substitution matrix as a nested dict, or None for 'hamming'.

    Resolution order
    ----------------
    1. Biopython (Bio.Align.substitution_matrices) — covers all names.
    2. Built-in embedded tables — blosum62, pam250.
    3. ValueError with install hint for anything else.
    """
    if name == "hamming":
        return None

    try:
        from Bio.Align import substitution_matrices as _bsm
        bmat  = _bsm.load(name.upper())
        alpha = bmat.alphabet
        return {a: {b: int(bmat[a, b]) for b in alpha} for a in alpha}
    except Exception:
        pass

    if name in _BUILTIN_MATRICES:
        return _BUILTIN_MATRICES[name]

    raise ValueError(
        f"Scoring matrix '{name}' is not built in and Biopython is not installed.\n"
        f"  pip install biopython          (then all --scoring options work)\n"
        f"  Built-in matrices: {', '.join(_BUILTIN_MATRICES)}"
    )


def _matrix_min(matrix: dict[str, dict[str, int]]) -> int:
    """Global minimum value across the entire substitution matrix."""
    return min(v for row in matrix.values() for v in row.values())


def _norm_score(a: str, b: str,
                matrix: dict[str, dict[str, int]],
                mat_min: int) -> float | None:
    """
    Normalized substitution score for residue pair (a, b) → [0, 1].

    Normalization anchor: the best self-score for this position pair,
    max(mat[a][a], mat[b][b]), maps to 1.0; mat_min maps to 0.0.

    Returns None when either residue is absent from the matrix
    (unknown characters such as X, B, Z are silently skipped).
    """
    a_row = matrix.get(a)
    b_row = matrix.get(b)
    if a_row is None or b_row is None:
        return None
    s = a_row.get(b)
    if s is None:
        return None
    s_ref = max(a_row.get(a, mat_min), b_row.get(b, mat_min))
    denom = s_ref - mat_min
    if denom <= 0:
        return 0.0
    return max(0.0, (s - mat_min) / denom)


# ---------------------------------------------------------------------------
# Similarity metric
# ---------------------------------------------------------------------------

GAP = frozenset("-.")


def _hamming_sim_window(seq1: str, seq2: str) -> float:
    """
    Normalized inverted Hamming distance for two same-length strings.

    Gap characters in either position are skipped (not counted as mismatches
    or as comparable positions).

        similarity = 1 - mismatches / comparable_positions

    Returns 0.0 when there are no comparable positions.
    """
    comparable = 0
    mismatches = 0
    for a, b in zip(seq1, seq2):
        if a in GAP or b in GAP:
            continue
        comparable += 1
        if a != b:
            mismatches += 1
    if comparable == 0:
        return 0.0
    return 1.0 - mismatches / comparable


def _score_window(seq1: str, seq2: str,
                  matrix: "dict | None" = None,
                  mat_min: int = 0) -> float:
    """
    Score a same-length window using the chosen scoring scheme.

    Hamming (matrix is None): normalized inverted Hamming distance.
    Substitution matrix: mean per-position normalized score → [0, 1].
    Positions where either residue is a gap or absent from the matrix
    are skipped; returns 0.0 when no comparable positions remain.
    """
    if matrix is None:
        return _hamming_sim_window(seq1, seq2)

    scores: list[float] = []
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a in GAP or b in GAP:
            continue
        ns = _norm_score(a, b, matrix, mat_min)
        if ns is not None:
            scores.append(ns)
    return sum(scores) / len(scores) if scores else 0.0


def sliding_best_similarity(seq1: str, seq2: str,
                             matrix: "dict | None" = None,
                             mat_min: int = 0) -> tuple[float, int]:
    """
    Slide the shorter sequence along the longer and return the best similarity.

    The shorter sequence is positioned at every valid offset where it fits
    entirely inside the longer sequence (offsets 0 … len(long) - len(short)).
    At each offset the score is computed by ``_score_window``, which uses
    either the Hamming metric or the supplied substitution matrix.

    Returns (best_similarity, best_offset) where best_offset is the 0-based
    start position of the *shorter* sequence within the *longer* sequence.

    When both sequences have the same length there is only one offset (0).
    """
    if len(seq1) >= len(seq2):
        long_, short = seq1, seq2
    else:
        long_, short = seq2, seq1

    L, S = len(long_), len(short)
    if S == 0:
        return 0.0, 0

    best_sim = -1.0
    best_off = 0
    for off in range(L - S + 1):
        sim = _score_window(short, long_[off: off + S], matrix, mat_min)
        if sim > best_sim:
            best_sim = sim
            best_off = off

    return best_sim, best_off


def hamming_similarity(seq1: str, seq2: str) -> float:
    """Convenience wrapper (Hamming only); returns only the similarity score."""
    sim, _ = sliding_best_similarity(seq1, seq2)
    return sim


# ---------------------------------------------------------------------------
# Conservation calculations
# ---------------------------------------------------------------------------

def pairwise_matrix(seqs: list[str],
                    matrix: "dict | None" = None,
                    mat_min: int = 0) -> tuple[list[list[float]], list[list[int]]]:
    """
    Compute all pairwise sliding-best similarities.

    Returns
    -------
    sim_mat : n×n float matrix  – similarity scores (Hamming or substitution)
    off_mat : n×n int matrix    – best_offset[i][j] is the start position of
                                  the shorter of (seqs[i], seqs[j]) within the
                                  longer.  Diagonal entries are 0.
    """
    n = len(seqs)
    sim_mat = [[0.0] * n for _ in range(n)]
    off_mat = [[0]   * n for _ in range(n)]
    for i in range(n):
        sim_mat[i][i] = 1.0
    for i, j in combinations(range(n), 2):
        sim, off = sliding_best_similarity(seqs[i], seqs[j], matrix, mat_min)
        sim_mat[i][j] = sim_mat[j][i] = sim
        off_mat[i][j] = off_mat[j][i] = off
    return sim_mat, off_mat


def best_offsets_vs_reference(seqs: list[str], ref_idx: int,
                               matrix: "dict | None" = None,
                               mat_min: int = 0) -> list[int]:
    """
    For each sequence find its best offset relative to the reference sequence.

    The reference sequence is always the longest (or the user-chosen index).
    Sequences longer than the reference are NOT slid; they receive offset 0
    and are truncated to the reference length when the virtual alignment is
    built (unusual edge case – the reference should be the longest).

    Returns a list of offsets, one per sequence (reference offset is 0).
    """
    ref = seqs[ref_idx]
    offsets: list[int] = []
    for i, seq in enumerate(seqs):
        if i == ref_idx:
            offsets.append(0)
        else:
            _, off = sliding_best_similarity(ref, seq, matrix, mat_min)
            offsets.append(off)
    return offsets


def build_virtual_alignment(seqs: list[str], ref_idx: int,
                             offsets: list[int]) -> list[str]:
    """
    Place each sequence at its best offset within the reference and pad with
    gap characters so every string is the same length as the reference.

    Sequences shorter than (offset + their length) are right-padded; sequences
    longer than the reference are truncated (should not happen when ref is the
    longest sequence).
    """
    ref_len = len(seqs[ref_idx])
    aligned: list[str] = []
    for i, seq in enumerate(seqs):
        off = offsets[i]
        padded = "-" * off + seq[: ref_len - off]
        padded = padded + "-" * (ref_len - len(padded))  # right-pad if needed
        aligned.append(padded)
    return aligned


def henikoff_weights(seqs: list[str]) -> list[float]:
    """
    Position-based sequence weights (Henikoff & Henikoff, 1994).

    At each alignment column the contribution of sequence *i* is
    ``1 / (r × n_i)``, where *r* is the number of distinct residue types
    at that column (gaps excluded) and *n_i* is the count of sequences
    sharing the same residue.  Fully conserved and all-gap columns
    contribute nothing.

    Weights are normalized to sum to *N* (number of sequences) so that an
    average weight equals 1.0 — values above 1 flag under-represented
    sequences, values below 1 flag over-represented ones.
    """
    n    = len(seqs)
    _GAP = {"-", "."}
    raw  = [0.0] * n

    for pos in range(len(seqs[0])):
        col_idx = [i for i in range(n) if seqs[i][pos] not in _GAP]
        if len(col_idx) < 2:
            continue
        counts: dict[str, int] = {}
        for i in col_idx:
            aa = seqs[i][pos].upper()
            counts[aa] = counts.get(aa, 0) + 1
        r = len(counts)
        if r == 1:          # fully conserved — no discriminating power
            continue
        for i in col_idx:
            aa = seqs[i][pos].upper()
            raw[i] += 1.0 / (r * counts[aa])

    total = sum(raw)
    if total == 0.0:        # all positions conserved or all gaps → equal weights
        return [1.0] * n
    return [w * n / total for w in raw]


def per_position_conservation(seqs: list[str],
                               matrix: "dict | None" = None,
                               mat_min: int = 0,
                               weights: "list[float] | None" = None) -> list[float]:
    """
    Per-column conservation for an aligned (equal-length) set of sequences.

    Hamming (matrix is None): (weighted) fraction of matching pairs.
    Substitution matrix: (weighted) mean normalized pairwise score.
    Gap characters are skipped.  All-gap columns score 0.0; single-residue
    columns score 1.0.

    When *weights* is supplied (e.g. from ``henikoff_weights``), each pair
    (i, j) is weighted by ``weights[i] × weights[j]``; otherwise all pairs
    are weighted equally.
    """
    length = len(seqs[0])
    _GAP   = {"-", "."}
    w      = weights if weights is not None else [1.0] * len(seqs)
    scores: list[float] = []

    for pos in range(length):
        idx = [i for i in range(len(seqs)) if seqs[i][pos] not in _GAP]
        if len(idx) < 2:
            scores.append(1.0 if len(idx) == 1 else 0.0)
            continue

        if matrix is None:
            num = den = 0.0
            for i, j in combinations(idx, 2):
                a, b  = seqs[i][pos].upper(), seqs[j][pos].upper()
                ww    = w[i] * w[j]
                num  += ww * (1.0 if a == b else 0.0)
                den  += ww
            scores.append(num / den if den else 1.0)
        else:
            num = den = 0.0
            for i, j in combinations(idx, 2):
                a, b = seqs[i][pos].upper(), seqs[j][pos].upper()
                ns   = _norm_score(a, b, matrix, mat_min)
                if ns is not None:
                    ww   = w[i] * w[j]
                    num += ww * ns
                    den += ww
            scores.append(num / den if den else 0.0)

    return scores


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _col(text: str, width: int, right: bool = False) -> str:
    return f"{text:>{width}}" if right else f"{text:<{width}}"


def print_sequence_table(records: list[tuple[str, str, str]]) -> None:
    print("\nInput sequences")
    print("=" * 50)
    print(f"{'Accession':<20} {'Database':<12} {'Length':>8}")
    print("-" * 42)
    for acc, db, seq in records:
        print(f"{acc:<20} {db:<12} {len(seq):>8}")


def print_pairwise_matrix(ids: list[str], sim_mat: list[list[float]],
                          off_mat: list[list[int]] | None = None) -> None:
    n = len(ids)
    W = 9
    id_w = max(len(i) for i in ids)

    print("\nPairwise similarity matrix (normalized inverted Hamming, best sliding offset)")
    print("=" * 70)
    header = " " * (id_w + 2) + "  ".join(f"{i[:W]:>{W}}" for i in ids)
    print(header)
    for i, row_id in enumerate(ids):
        row = f"{row_id:<{id_w}}  " + "  ".join(f"{sim_mat[i][j]:>{W}.4f}"
                                                  for j in range(n))
        print(row)

    if off_mat is not None:
        print()
        print("Best offset matrix (start position of shorter seq within longer)")
        print("-" * 70)
        header2 = " " * (id_w + 2) + "  ".join(f"{i[:W]:>{W}}" for i in ids)
        print(header2)
        for i, row_id in enumerate(ids):
            row = f"{row_id:<{id_w}}  " + "  ".join(f"{off_mat[i][j]:>{W}}"
                                                      for j in range(n))
            print(row)


def print_pairwise_summary(ids: list[str], seqs: list[str],
                           sim_mat: list[list[float]],
                           off_mat: list[list[int]]) -> None:
    pairs = [
        (sim_mat[i][j], off_mat[i][j], ids[i], ids[j], len(seqs[i]), len(seqs[j]))
        for i, j in combinations(range(len(ids)), 2)
    ]
    sims = [s for s, *_ in pairs]

    print("\nPairwise similarity summary")
    print("=" * 50)
    print(f"  Pairs:    {len(sims)}")
    print(f"  Mean:     {sum(sims)/len(sims):.4f}")
    print(f"  Min:      {min(sims):.4f}")
    print(f"  Max:      {max(sims):.4f}")

    pairs_sorted = sorted(pairs, reverse=True)
    show = min(3, len(pairs_sorted))

    def _offset_note(off: int, la: int, lb: int) -> str:
        """Describe offset: which seq is shorter and where it was placed."""
        if la == lb:
            return ""
        shorter_id, longer_len = ("A", lb) if la <= lb else ("B", la)
        return f"  [shorter={shorter_id}, offset={off}, window={min(la,lb)}/{longer_len}]"

    print(f"\n  Most similar pair(s):")
    for sim, off, a, b, la, lb in pairs_sorted[:show]:
        print(f"    {a} <-> {b}:  {sim:.4f}{_offset_note(off, la, lb)}")

    if len(pairs_sorted) > show:
        print(f"\n  Least similar pair(s):")
        for sim, off, a, b, la, lb in pairs_sorted[-show:]:
            print(f"    {a} <-> {b}:  {sim:.4f}{_offset_note(off, la, lb)}")


def print_position_summary(scores: list[float], seqs: list[str], top_n: int) -> None:
    fully = sum(1 for s in scores if s >= 1.0)
    variable = sum(1 for s in scores if s < 1.0)
    mean = sum(scores) / len(scores)

    print("\nPer-position conservation summary")
    print("=" * 50)
    print(f"  Alignment length:           {len(scores)}")
    print(f"  Mean conservation:          {mean:.4f}")
    print(f"  Fully conserved positions:  {fully}")
    print(f"  Variable positions:         {variable}")

    indexed = sorted(enumerate(scores), key=lambda x: x[1], reverse=True)
    GAP = {"-", "."}

    def residue_str(pos: int) -> str:
        return "".join(s[pos] for s in seqs if s[pos] not in GAP)

    print(f"\n  Top {top_n} most conserved positions")
    print(f"  {'Pos':>6}  {'Score':>7}  Residues")
    print(f"  {'---':>6}  {'-----':>7}  --------")
    for pos, score in indexed[:top_n]:
        print(f"  {pos+1:>6}  {score:>7.4f}  {residue_str(pos)}")

    if len(indexed) > top_n:
        print(f"\n  Top {top_n} least conserved positions")
        print(f"  {'Pos':>6}  {'Score':>7}  Residues")
        print(f"  {'---':>6}  {'-----':>7}  --------")
        for pos, score in indexed[-top_n:]:
            print(f"  {pos+1:>6}  {score:>7.4f}  {residue_str(pos)}")


def _guard(path: "str | pathlib.Path", overwrite: bool) -> None:
    """Raise FileExistsError if *path* exists and *overwrite* is False."""
    if not overwrite and pathlib.Path(path).exists():
        raise FileExistsError(
            f"File already exists (use --overwrite to replace): {path}"
        )


def write_matrix_csv(path: str, ids: list[str], mat: list[list[float]]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + ids)
        for i, row_id in enumerate(ids):
            writer.writerow([row_id] + [f"{mat[i][j]:.6f}" for j in range(len(ids))])
    print(f"  Matrix written to: {path}")


def write_conservation_csv(path: str, scores: list[float]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["position", "conservation"])
        for i, s in enumerate(scores, 1):
            writer.writerow([i, f"{s:.6f}"])
    print(f"  Per-position conservation written to: {path}")


# ---------------------------------------------------------------------------
# Multi-file input helpers
# ---------------------------------------------------------------------------

_FASTA_EXTS = frozenset({".fasta", ".fa", ".faa", ".fna", ".fas"})


def expand_input_paths(
    paths: list[str], recursive: bool = False
) -> list[pathlib.Path]:
    """
    Expand a mixed list of file paths and directory paths to a deduplicated,
    sorted list of FASTA file Paths.

    Directories are scanned for files whose extension is in _FASTA_EXTS.
    With *recursive=True* the scan descends into sub-directories.
    Unrecognised paths are skipped with a warning.
    """
    found: list[pathlib.Path] = []
    for raw in paths:
        p = pathlib.Path(raw)
        if p.is_dir():
            pat = "**/*" if recursive else "*"
            hits = sorted(
                h for ext in _FASTA_EXTS
                for h in p.glob(f"{pat}{ext}")
                if h.is_file()
            )
            if not hits:
                print(
                    f"  Warning: directory '{raw}' contains no FASTA files "
                    f"({', '.join(_FASTA_EXTS)}).",
                    file=sys.stderr,
                )
            found.extend(hits)
        elif p.is_file():
            found.append(p)
        else:
            print(
                f"  Warning: '{raw}' is not a file or directory – skipped.",
                file=sys.stderr,
            )

    # Deduplicate by resolved absolute path, preserving encounter order
    seen: set[pathlib.Path] = set()
    result: list[pathlib.Path] = []
    for p in found:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            result.append(p)
    return result


def _resolve_out_path(
    flag_val: str | None,       # None → disabled; "" → auto-name; other → explicit
    stem: str,                  # input file stem, e.g. "proteins"
    suffix: str,                # e.g. "_matrix"
    ext: str,                   # with dot, e.g. ".csv"
    input_parent: pathlib.Path,
    out_dir: pathlib.Path | None,
    multi_file: bool,
) -> pathlib.Path | None:
    """
    Resolve an output path for one file in either single- or multi-file mode.

    Rules
    -----
    flag not given  → return None (output disabled)
    flag="" (const) → auto-name: {base}/{stem}{suffix}{ext}
    flag=explicit   → single-file: use as-is
                    → multi-file: if it's a directory use it as base;
                                  if it looks like a filename, warn + auto-name
    """
    if flag_val is None:
        return None

    auto_name = f"{stem}{suffix}{ext}"
    default_base = out_dir if out_dir is not None else input_parent

    if flag_val == "":
        return default_base / auto_name

    p = pathlib.Path(flag_val)

    if multi_file:
        # Is the given value a directory (existing or trailing sep)?
        is_dir_path = p.is_dir() or flag_val.endswith(("/", "\\"))
        if is_dir_path:
            return p / auto_name
        if p.suffix:
            # Looks like a specific filename → would be clobbered; warn
            print(
                f"  Note: '{flag_val}' is a filename, not a directory. "
                f"In multi-file mode outputs are auto-named; "
                f"writing to '{default_base / auto_name}' instead.",
                file=sys.stderr,
            )
            return default_base / auto_name
        # No extension → treat as directory path
        return p / auto_name

    # Single-file mode: honour the explicit path
    return p


def _resolve_heatmap_path(
    flag_val: str | None,
    stem: str,
    input_parent: pathlib.Path,
    out_dir: pathlib.Path | None,
    multi_file: bool,
) -> pathlib.Path | None:
    """
    Like _resolve_out_path but the output extension comes from the flag value
    itself.  Defaults to .png when auto-naming.
    """
    if flag_val is None:
        return None

    default_base = out_dir if out_dir is not None else input_parent
    _IMG_EXTS = frozenset({".png", ".pdf", ".svg", ".eps", ".jpg", ".jpeg"})

    if flag_val == "":
        return default_base / f"{stem}_heatmap.png"

    p = pathlib.Path(flag_val)

    if multi_file:
        if p.suffix.lower() in _IMG_EXTS:
            # Treat the suffix as the desired format, auto-name the stem
            return default_base / f"{stem}_heatmap{p.suffix}"
        if p.is_dir() or not p.suffix:
            return p / f"{stem}_heatmap.png"
        print(
            f"  Note: '{flag_val}' treated as output directory for heatmaps.",
            file=sys.stderr,
        )
        return p / f"{stem}_heatmap.png"

    return p


def _resolve_dendrogram_path(
    flag_val: str | None,
    stem: str,
    input_parent: pathlib.Path,
    out_dir: pathlib.Path | None,
    multi_file: bool,
) -> pathlib.Path | None:
    """Like _resolve_heatmap_path but for dendrograms."""
    if flag_val is None:
        return None

    default_base = out_dir if out_dir is not None else input_parent
    _IMG_EXTS = frozenset({".png", ".pdf", ".svg", ".eps", ".jpg", ".jpeg"})

    if flag_val == "":
        return default_base / f"{stem}_dendrogram.png"

    p = pathlib.Path(flag_val)

    if multi_file:
        if p.suffix.lower() in _IMG_EXTS:
            return default_base / f"{stem}_dendrogram{p.suffix}"
        if p.is_dir() or not p.suffix:
            return p / f"{stem}_dendrogram.png"
        return p / f"{stem}_dendrogram.png"

    return p


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--version", action="version", version=f"%(prog)s {VERSION}",
    )
    # ── input ──────────────────────────────────────────────────────────────
    p.add_argument(
        "fasta", nargs="+",
        metavar="PATH",
        help=(
            "One or more FASTA files and/or directories to analyse.  "
            "Directories are scanned for files with extensions "
            + ", ".join(sorted(_FASTA_EXTS)) + ".  "
            "Multiple paths are processed in sequence and a cross-file "
            "summary is printed at the end."
        ),
    )
    p.add_argument(
        "--recursive", action="store_true",
        help="Scan directories recursively for FASTA files.",
    )
    # ── analysis ───────────────────────────────────────────────────────────
    p.add_argument(
        "--matrix", action="store_true",
        help="Print the full pairwise similarity and offset matrices to stdout",
    )
    p.add_argument(
        "--top", type=int, default=10, metavar="N",
        help="Number of most/least conserved positions to display (default: 10)",
    )
    p.add_argument(
        "--mode", choices=["auto", "aligned", "unaligned"], default="auto",
        help=(
            "Sequence input mode.  "
            "'auto' (default): equal-length → MSA column comparison, "
            "unequal-length → sliding approach.  "
            "'aligned': treat input as a pre-computed MSA; all sequences must "
            "share the same length; gap characters (- .) mark alignment columns "
            "and are skipped in similarity calculations.  "
            "'unaligned': strip any gap characters first, then always use the "
            "sliding approach regardless of length."
        ),
    )
    p.add_argument(
        "--scoring",
        choices=SCORING_CHOICES,
        default="hamming",
        metavar="MATRIX",
        help=(
            "Scoring scheme for pairwise similarity and per-position "
            "conservation.  "
            "'hamming' (default): normalized inverted Hamming distance "
            "(binary match/mismatch, gap-aware).  "
            "Substitution matrix names (blosum45, blosum62, blosum80, "
            "blosum90, pam30, pam70, pam120, pam250): mean per-position "
            "normalized substitution score, mapped to [0, 1] so that a "
            "perfect self-alignment scores 1.0.  "
            "'blosum62' and 'pam250' are built-in; all others additionally "
            "require Biopython (pip install biopython).  "
            f"Choices: {', '.join(SCORING_CHOICES)}."
        ),
    )
    # ── output files ───────────────────────────────────────────────────────
    p.add_argument(
        "--out-dir", metavar="DIR",
        help=(
            "Directory for all file outputs.  In multi-file mode output names "
            "are derived automatically from each input file's stem.  "
            "Created if it does not exist."
        ),
    )
    p.add_argument(
        "--out-matrix", metavar="FILE", nargs="?", const="",
        help=(
            "Write pairwise similarity matrix to a CSV file.  "
            "FILE is required in single-file mode.  "
            "Omit FILE (or use with --out-dir) to auto-name as "
            "{stem}_matrix.csv."
        ),
    )
    p.add_argument(
        "--out-conservation", metavar="FILE", nargs="?", const="",
        help=(
            "Write per-position conservation scores to a CSV file.  "
            "Omit FILE to auto-name as {stem}_conservation.csv."
        ),
    )
    # ── taxonomy options ───────────────────────────────────────────────────
    p.add_argument(
        "--rank", metavar="RANK",
        help=(
            "Taxonomic rank at which to group sequences and compute median "
            "pairwise conservation (e.g. species, genus, family, order, class, "
            "phylum, superkingdom).  Requires network access to UniProt / NCBI."
        ),
    )
    p.add_argument(
        "--ncbi-api-key", metavar="KEY",
        help="NCBI API key (raises rate limit from 3 to 10 requests/second).",
    )
    p.add_argument(
        "--cache-file", metavar="FILE", default=_DEFAULT_CACHE,
        help=(
            f"JSON file used to cache taxonomy lookups between runs "
            f"(default: {_DEFAULT_CACHE}).  Pass an empty string to disable caching."
        ),
    )
    p.add_argument(
        "--no-cache", action="store_true",
        help="Disable the taxonomy cache entirely (always fetch from network).",
    )
    p.add_argument(
        "--heatmap", metavar="FILE", nargs="?", const="",
        help=(
            "Save a similarity heatmap.  Format is inferred from the extension "
            "(.png, .pdf, .svg, …); defaults to .png when auto-naming.  "
            "Requires matplotlib.  Omit FILE to auto-name as {stem}_heatmap.png."
        ),
    )
    p.add_argument(
        "--dendrogram", metavar="FILE", nargs="?", const="",
        help=(
            "Save a hierarchical-clustering dendrogram of the sequences.  "
            "Distance = 1 − similarity; linkage method is UPGMA (average).  "
            "Format is inferred from the extension (.png, .pdf, .svg, …); "
            "defaults to .png when auto-naming.  "
            "Requires matplotlib and scipy (pip install matplotlib scipy).  "
            "Omit FILE to auto-name as {stem}_dendrogram.png."
        ),
    )
    p.add_argument(
        "--dendrogram-method",
        metavar="METHOD",
        default="average",
        help=(
            "Linkage method for the dendrogram (default: average / UPGMA).  "
            "Any method accepted by scipy.cluster.hierarchy.linkage is valid: "
            "single, complete, average, weighted, ward, centroid, median.  "
            "Also controls the tree topology in --out-phyloxml."
        ),
    )
    p.add_argument(
        "--out-phyloxml", metavar="FILE", nargs="?", const="",
        help=(
            "Export the hierarchical clustering tree as a PhyloXML file.  "
            "The same linkage method as --dendrogram-method is used.  "
            "Requires scipy (pip install scipy).  "
            "Omit FILE to auto-name as {stem}_dendrogram.xml."
        ),
    )
    p.add_argument(
        "--out-rank-csv", metavar="FILE", nargs="?", const="",
        help=(
            "Write per-group rank-conservation table to a CSV file.  "
            "Omit FILE to auto-name as {stem}_rank_{rank}.csv."
        ),
    )
    p.add_argument(
        "--weights", action="store_true",
        help=(
            "Apply Henikoff & Henikoff (1994) position-based sequence weights "
            "when computing per-position conservation.  Sequences that are "
            "over-represented in the input are down-weighted; rare sequences "
            "are up-weighted.  Has no effect unless the virtual alignment can "
            "be constructed (requires ≥2 sequences)."
        ),
    )
    p.add_argument(
        "--permutations", type=int, default=0, metavar="N",
        help=(
            "Run a permutation test for rank conservation significance.  "
            "N random shuffles of group labels are performed and the "
            "one-tailed p-value is the fraction of shuffles with "
            "T = intra_median − inter_median ≥ observed T.  "
            "Requires --rank.  Typical values: 1000 or 10000.  "
            "Default: 0 (test disabled)."
        ),
    )
    p.add_argument(
        "--overwrite", action="store_true",
        help=(
            "Allow output files to be overwritten if they already exist.  "
            "By default the tool skips any output file that is already present "
            "and prints a warning, leaving the existing file untouched."
        ),
    )
    p.add_argument(
        "--report", metavar="FILE",
        help=(
            "Write a self-contained HTML report combining all results "
            "(sequence tables, coloured similarity matrices, per-position "
            "conservation plots, and taxonomy heatmaps).  "
            "The format is inferred from the extension: "
            ".html (default) produces HTML only; "
            ".pdf additionally converts to PDF via WeasyPrint "
            "(pip install weasyprint) and saves the HTML alongside.  "
            "Respects --out-dir for relative paths."
        ),
    )
    return p


def _write_rank_csv(path: pathlib.Path, result: dict) -> None:
    rank   = result["rank"]
    groups = result["groups"]
    rows: list[dict] = []
    for taxon, gdata in sorted(groups.items()):
        med = gdata["intra_median"]
        rows.append({
            rank: taxon,
            "members": len(gdata["members"]),
            "intra_pairs": len(gdata["intra_pairs"]),
            "intra_median_similarity": f"{med:.6f}" if med is not None else "",
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Rank conservation written to: {path}")


# ---------------------------------------------------------------------------
# Per-file analysis  (called once per FASTA in both single and multi-file mode)
# ---------------------------------------------------------------------------

def _analyse_fasta(
    args:             argparse.Namespace,
    fasta_path:       pathlib.Path,
    out_matrix:       pathlib.Path | None,
    out_conservation: pathlib.Path | None,
    out_rank_csv:     pathlib.Path | None,
    heatmap_path:     pathlib.Path | None,
    dendrogram_path:  pathlib.Path | None,
    phyloxml_path:    pathlib.Path | None,
    cache:            "TaxCache | None",
) -> dict:
    """
    Run the full conservation analysis for one FASTA file.

    Returns a summary dict consumed by the multi-file table.
    On error the dict contains 'error' with a message and all numeric
    fields set to None.
    """
    summary: dict = {
        "file": fasta_path,
        "n_seqs": 0, "n_kept": 0, "mode": "–",
        "sim_median": None, "intra_median": None, "inter_median": None,
        "outputs": [], "error": None,
    }

    try:
        # ── load ──────────────────────────────────────────────────────────
        print(f"Reading: {fasta_path}")
        records = parse_fasta(str(fasta_path))
        if not records:
            summary["error"] = "no sequences found"
            return summary
        summary["n_seqs"] = len(records)

        # ── taxonomy (optional) ───────────────────────────────────────────
        tax_map: dict[str, TaxInfo] = {}
        if args.rank:
            print(f"\nFetching taxonomy for {len(records)} sequence(s)…")
            tax_map = fetch_all_taxonomies(
                records, api_key=args.ncbi_api_key, cache=cache
            )
            before  = len(records)
            records = [r for r in records if r[0] in tax_map]
            dropped = before - len(records)
            if dropped:
                print(f"  Dropped {dropped} sequence(s) with no taxonomy.")
            if len(records) < 2:
                summary["error"] = "fewer than 2 sequences have retrievable taxonomy"
                return summary

        summary["n_kept"] = len(records)
        print_sequence_table(records)

        if len(records) < 2:
            summary["error"] = "fewer than 2 sequences"
            return summary

        ids  = [r[0] for r in records]
        seqs = [r[2] for r in records]

        # ── mode / gap handling ───────────────────────────────────────────
        mode = args.mode
        gap_re = re.compile(r"[-.]")

        if mode == "unaligned":
            stripped   = [gap_re.sub("", s) for s in seqs]
            n_stripped = sum(1 for o, s in zip(seqs, stripped) if len(o) != len(s))
            if n_stripped:
                print(f"\nMode: unaligned – gap characters stripped from "
                      f"{n_stripped} sequence(s).")
            seqs = stripped

        lengths = [len(s) for s in seqs]

        if mode == "aligned":
            if len(set(lengths)) != 1:
                summary["error"] = (
                    f"--mode aligned requires equal lengths "
                    f"(found {sorted(set(lengths))})"
                )
                return summary
            has_gaps = any(gap_re.search(s) for s in seqs)
            if not has_gaps:
                print(
                    "\nMode: aligned – no gap characters (- .) found; "
                    "treating sequences as column-aligned."
                )
            else:
                print(
                    f"\nMode: aligned – MSA columns used directly "
                    f"({lengths[0]} positions, gap characters skipped)."
                )
            is_aligned = True

        elif mode == "unaligned":
            is_aligned = False
            print(
                f"\nMode: unaligned – sliding-window approach "
                f"(lengths {min(lengths)}–{max(lengths)} aa)."
            )

        else:  # auto
            is_aligned = len(set(lengths)) == 1
            if is_aligned:
                print(
                    "\nMode: auto – equal-length sequences; "
                    "treating as MSA (standard Hamming similarity)."
                )
            else:
                print(
                    f"\nMode: auto – unequal lengths "
                    f"({min(lengths)}–{max(lengths)} aa); "
                    "using sliding-window approach."
                )

        summary["mode"] = (
            "aligned"   if is_aligned  else "unaligned"
        )

        # ── load scoring matrix ───────────────────────────────────────────
        scoring = getattr(args, "scoring", "hamming")
        try:
            score_matrix = load_matrix(scoring)
        except ValueError as exc:
            summary["error"] = str(exc)
            return summary
        mat_min = _matrix_min(score_matrix) if score_matrix else 0

        # ── pairwise similarity ───────────────────────────────────────────
        label = scoring.upper() if scoring != "hamming" else "Hamming"
        print(f"\nScoring: {label}")
        print("Computing pairwise similarities…")
        sim_mat, off_mat = pairwise_matrix(seqs, score_matrix, mat_min)
        print_pairwise_summary(ids, seqs, sim_mat, off_mat)

        all_pairs = [
            sim_mat[i][j]
            for i, j in combinations(range(len(ids)), 2)
        ]
        summary["sim_median"] = statistics.median(all_pairs) if all_pairs else None

        if args.matrix:
            print_pairwise_matrix(
                ids, sim_mat, off_mat if not is_aligned else None
            )

        if out_matrix is not None:
            out_matrix.parent.mkdir(parents=True, exist_ok=True)
            try:
                _guard(out_matrix, args.overwrite)
                write_matrix_csv(str(out_matrix), ids, sim_mat)
                summary["outputs"].append(str(out_matrix))
            except FileExistsError as exc:
                print(f"  Skipped: {exc}", file=sys.stderr)

        dend_method = getattr(args, "dendrogram_method", "average")

        if dendrogram_path is not None:
            dendrogram_path.parent.mkdir(parents=True, exist_ok=True)
            try:
                _guard(dendrogram_path, args.overwrite)
                plot_dendrogram(ids, sim_mat, str(dendrogram_path),
                                method=dend_method)
                summary["outputs"].append(str(dendrogram_path))
            except FileExistsError as exc:
                print(f"  Skipped: {exc}", file=sys.stderr)

        if phyloxml_path is not None:
            phyloxml_path.parent.mkdir(parents=True, exist_ok=True)
            try:
                _guard(phyloxml_path, args.overwrite)
                export_phyloxml(ids, sim_mat, str(phyloxml_path),
                                method=dend_method)
                summary["outputs"].append(str(phyloxml_path))
            except FileExistsError as exc:
                print(f"  Skipped: {exc}", file=sys.stderr)

        # ── per-position conservation ─────────────────────────────────────
        if is_aligned:
            print("\nComputing per-position conservation…")
            aligned_seqs = seqs
        else:
            ref_idx = lengths.index(max(lengths))
            print(
                f"\nBuilding virtual alignment using '{ids[ref_idx]}' as "
                f"reference (length {lengths[ref_idx]} aa)…"
            )
            offsets = best_offsets_vs_reference(seqs, ref_idx,
                                                   score_matrix, mat_min)
            for i, (seq_id, off) in enumerate(zip(ids, offsets)):
                if i != ref_idx:
                    print(f"  {seq_id}: placed at offset {off}"
                          f"  (covers positions {off+1}–{off+len(seqs[i])})")
            aligned_seqs = build_virtual_alignment(seqs, ref_idx, offsets)

        seq_weights: "list[float] | None" = None
        if getattr(args, "weights", False):
            seq_weights = henikoff_weights(aligned_seqs)
            print(f"\nHenikoff sequence weights: "
                  f"{', '.join(f'{w:.3f}' for w in seq_weights)}")

        scores = per_position_conservation(aligned_seqs, score_matrix, mat_min,
                                           weights=seq_weights)
        print_position_summary(scores, aligned_seqs, args.top)

        if out_conservation is not None:
            out_conservation.parent.mkdir(parents=True, exist_ok=True)
            try:
                _guard(out_conservation, args.overwrite)
                write_conservation_csv(str(out_conservation), scores)
                summary["outputs"].append(str(out_conservation))
            except FileExistsError as exc:
                print(f"  Skipped: {exc}", file=sys.stderr)

        # ── rank-based taxonomy conservation ──────────────────────────────
        rank_result: dict | None = None
        perm_result: dict | None = None
        if args.rank:
            rank_result = rank_conservation(
                ids, seqs, sim_mat, tax_map, args.rank
            )
            n_perm = getattr(args, "permutations", 0)
            if n_perm > 0:
                print(f"\nRunning permutation test ({n_perm} permutations)…")
                perm_result = permutation_test_rank(
                    ids, sim_mat, tax_map, args.rank,
                    n_permutations=n_perm,
                )
            print_rank_conservation(rank_result, ids, tax_map, perm_result)

            summary["intra_median"] = rank_result["intra_median"]
            summary["inter_median"] = rank_result["inter_median"]

            if heatmap_path is not None:
                try:
                    _guard(heatmap_path, args.overwrite)
                    plot_rank_heatmap(rank_result, str(heatmap_path))
                    summary["outputs"].append(str(heatmap_path))
                except FileExistsError as exc:
                    print(f"  Skipped: {exc}", file=sys.stderr)

            missing = [
                acc for acc in ids
                if acc in tax_map and tax_map[acc].at_rank(args.rank) is None
            ]
            if missing:
                avail = sorted({
                    r for acc in ids if acc in tax_map
                    for r in tax_map[acc].available_ranks()
                })
                print(
                    f"\n  Note: {len(missing)} sequence(s) have no "
                    f"'{args.rank}' rank in their lineage and were excluded."
                )
                print(f"  Available ranks in this dataset: {', '.join(avail)}")

            if out_rank_csv is not None and rank_result["groups"]:
                try:
                    _guard(out_rank_csv, args.overwrite)
                    _write_rank_csv(out_rank_csv, rank_result)
                    summary["outputs"].append(str(out_rank_csv))
                except FileExistsError as exc:
                    print(f"  Skipped: {exc}", file=sys.stderr)

        elif heatmap_path is not None:
            print(
                "  Note: --heatmap requires --rank to be set; "
                "heatmap skipped for this file.",
                file=sys.stderr,
            )

        # ── collect data for HTML/PDF report ──────────────────────────────
        summary["_report"] = {
            "records":             records,
            "ids":                 ids,
            "sim_mat":             sim_mat,
            "off_mat":             off_mat,
            "conservation_scores": scores,
            "is_aligned":          is_aligned,
            "scoring":             scoring,
            "rank_result":         rank_result,
            "tax_map":             tax_map,
            "perm_result":         perm_result,
            "seq_weights":         seq_weights,
            "dendrogram_method":   getattr(args, "dendrogram_method", "average"),
        }

    except Exception as exc:          # noqa: BLE001 – catch-all for per-file resilience
        summary["error"] = str(exc)

    return summary


# ---------------------------------------------------------------------------
# Multi-file summary table
# ---------------------------------------------------------------------------

def _print_multifile_summary(results: list[dict], rank: str | None) -> None:
    """Print a consolidated one-row-per-file summary after all files are done."""
    print("\n" + "━" * 72)
    print("Multi-file summary")
    print("━" * 72)

    # Column widths
    fw = max(len(str(r["file"].name)) for r in results)
    fw = max(fw, 20)

    hdr_parts = [
        f"{'File':<{fw}}",
        f"{'Seqs':>5}",
        f"{'Kept':>5}",
        f"{'Mode':<10}",
        f"{'Sim med':>8}",
    ]
    if rank:
        hdr_parts += [f"{'Intra':>8}", f"{'Inter':>8}"]
    hdr_parts.append(f"  {'Status'}")
    print("  " + "  ".join(hdr_parts))
    print("  " + "  ".join(
        "-" * len(h.strip()) if h.strip() else "-" * len(h) for h in hdr_parts
    ))

    for r in results:
        name = r["file"].name
        if r["error"]:
            status = f"ERROR: {r['error']}"
            row = (
                f"  {name:<{fw}}  {r['n_seqs']:>5}  {r['n_kept']:>5}  "
                f"{'–':<10}  {'–':>8}"
            )
            if rank:
                row += f"  {'–':>8}  {'–':>8}"
            row += f"  {status}"
        else:
            sm  = f"{r['sim_median']:.4f}"  if r["sim_median"]  is not None else "–"
            inm = f"{r['intra_median']:.4f}" if r["intra_median"] is not None else "–"
            inr = f"{r['inter_median']:.4f}" if r["inter_median"] is not None else "–"
            row = (
                f"  {name:<{fw}}  {r['n_seqs']:>5}  {r['n_kept']:>5}  "
                f"{r['mode']:<10}  {sm:>8}"
            )
            if rank:
                row += f"  {inm:>8}  {inr:>8}"
            n_out = len(r["outputs"])
            row += f"  OK{f' ({n_out} file(s) written)' if n_out else ''}"
        print(row)

    n_ok  = sum(1 for r in results if not r["error"])
    n_err = len(results) - n_ok
    print(f"\n  {n_ok} file(s) completed successfully"
          + (f", {n_err} with errors" if n_err else "") + ".")


# ---------------------------------------------------------------------------
# HTML / PDF report
# ---------------------------------------------------------------------------

def _h(s: object) -> str:
    """HTML-escape any value."""
    return (str(s)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;"))


def _sim_color(v: float) -> str:
    """HSL background colour for a similarity value in [0, 1]."""
    try:
        if v != v:          # NaN guard
            return "#d8d8d8"
        v = max(0.0, min(1.0, float(v)))
        return f"hsl({int(v * 120)},60%,88%)"
    except Exception:
        return "#d8d8d8"


def _conservation_plot_b64(scores: list[float], title: str) -> str | None:
    """Render a per-position conservation line chart; return base64 PNG or None."""
    try:
        import base64, io
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    fig, ax = plt.subplots(figsize=(11, 2.8))
    x = list(range(1, len(scores) + 1))
    ax.fill_between(x, scores, alpha=0.25, color="#2c6fad")
    ax.plot(x, scores, color="#1a4e8c", linewidth=0.9, rasterized=True)
    if scores:
        med = statistics.median(scores)
        ax.axhline(med, color="#c0392b", linestyle="--", linewidth=1,
                   alpha=0.8, label=f"median = {med:.3f}")
        ax.legend(fontsize=8, loc="upper right")
    ax.set_xlim(1, max(len(scores), 1))
    ax.set_ylim(-0.02, 1.08)
    ax.set_xlabel("Position", fontsize=9)
    ax.set_ylabel("Conservation", fontsize=9)
    ax.set_title(title, fontsize=10, pad=6)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, alpha=0.25, linewidth=0.5)
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def _heatmap_b64(rank_result: dict) -> str | None:
    """Render the rank heatmap to an embedded base64 PNG string."""
    import os, tempfile, base64
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        tmp = f.name
    try:
        plot_rank_heatmap(rank_result, tmp)
        with open(tmp, "rb") as f:
            data = f.read()
        return base64.b64encode(data).decode() if data else None
    except Exception:
        return None
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def _dendrogram_b64(ids: list[str], sim_mat: list[list[float]],
                    method: str = "average") -> str | None:
    """Render a dendrogram to an embedded base64 PNG string."""
    import os, tempfile, base64
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        tmp = f.name
    try:
        plot_dendrogram(ids, sim_mat, tmp, method=method)
        with open(tmp, "rb") as f:
            data = f.read()
        return base64.b64encode(data).decode() if data else None
    except Exception:
        return None
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


_REPORT_CSS = """\
:root{--blue:#2c6fad;--dark:#1a1a2e;--border:#dde3ed;--bg2:#f8f9fc}
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif;
     font-size:14px;line-height:1.6;color:#333;max-width:1280px;
     margin:0 auto;padding:2rem;background:#fff}
h1{font-size:1.75rem;color:var(--dark);border-bottom:3px solid var(--blue);
   padding-bottom:.5rem;margin-bottom:.5rem}
h2{font-size:1.2rem;color:var(--dark);border-bottom:1px solid var(--border);
   padding-bottom:.3rem;margin:2rem 0 .7rem}
h3{font-size:1rem;color:#0f3460;margin:1.2rem 0 .4rem}
.meta{color:#666;font-size:.88rem;margin:.3rem 0 1rem}
.toc{background:var(--bg2);border:1px solid var(--border);border-radius:8px;
     padding:1.2rem 1.6rem;margin-bottom:2rem}
.toc h2{margin-top:0;border:none;padding:0}
.toc ol{padding-left:1.4rem;margin-top:.5rem}
.toc li{margin:.2rem 0}
.toc a{color:var(--blue);text-decoration:none}
.toc a:hover{text-decoration:underline}
.file-section{border:1px solid var(--border);border-radius:8px;
              padding:1.4rem 1.8rem;margin:1.5rem 0}
.error-box{background:#fef0f0;border:1px solid #f5c6cb;border-radius:6px;
           color:#c0392b;padding:.8rem 1rem;margin:.5rem 0}
table{border-collapse:collapse;margin:.8rem 0;font-size:.88rem;width:100%}
th{background:var(--blue);color:#fff;padding:.45rem .8rem;text-align:left;
   font-weight:600;white-space:nowrap}
td{padding:.35rem .8rem;border-bottom:1px solid var(--border)}
tr:last-child td{border-bottom:none}
tbody tr:hover td{background:#f0f4ff}
.matrix-wrap{overflow-x:auto;margin:.8rem 0}
.matrix th,.matrix td{text-align:center;padding:.35rem .55rem;
                       font-size:.82rem;white-space:nowrap}
.matrix th:first-child,.matrix td:first-child{text-align:left;background:var(--bg2);
                                               font-weight:600}
.bold{font-weight:700}
img.plot{max-width:100%;height:auto;margin:.6rem 0;display:block;
         border:1px solid var(--border);border-radius:4px}
.badge{display:inline-block;padding:.15rem .5rem;border-radius:4px;
       font-size:.78rem;font-weight:600;margin-left:.4rem;vertical-align:middle}
.badge-ok{background:#d4edda;color:#155724}
.badge-err{background:#f8d7da;color:#721c24}
@media print{
  .toc{page-break-after:always}
  .file-section{page-break-inside:avoid}
  a{color:inherit;text-decoration:none}
}"""


def _html_summary_table(results: list[dict], rank: str | None) -> str:
    rank_cols = (
        f"<th>Intra-{_h(rank)}</th><th>Inter-{_h(rank)}</th>" if rank else ""
    )
    rows = []
    for i, r in enumerate(results, 1):
        name = _h(r["file"].name)
        href = f"#file-{i}"
        if r["error"]:
            rows.append(
                f'<tr><td><a href="{href}">{name}</a></td>'
                f'<td>{r["n_seqs"]}</td><td>{r["n_kept"]}</td>'
                f'<td>–</td><td>–</td>'
                + ("<td>–</td><td>–</td>" if rank else "")
                + f'<td><span class="badge badge-err">ERROR</span> {_h(r["error"])}</td></tr>'
            )
        else:
            sm  = f'{r["sim_median"]:.4f}'   if r["sim_median"]  is not None else "–"
            inm = f'{r["intra_median"]:.4f}' if r["intra_median"] is not None else "–"
            inr = f'{r["inter_median"]:.4f}' if r["inter_median"] is not None else "–"
            rank_cells = f"<td>{inm}</td><td>{inr}</td>" if rank else ""
            rows.append(
                f'<tr><td><a href="{href}">{name}</a></td>'
                f'<td>{r["n_seqs"]}</td><td>{r["n_kept"]}</td>'
                f'<td>{_h(r.get("mode", "–"))}</td><td>{sm}</td>'
                f'{rank_cells}'
                f'<td><span class="badge badge-ok">OK</span></td></tr>'
            )
    return (
        '<table>'
        '<thead><tr><th>File</th><th>Seqs</th><th>Kept</th>'
        f'<th>Mode</th><th>Sim median</th>{rank_cols}<th>Status</th>'
        '</tr></thead>'
        f'<tbody>{"".join(rows)}</tbody>'
        '</table>'
    )


def _html_seq_table(records: list) -> str:
    rows = "".join(
        f'<tr><td class="bold">{_h(acc)}</td><td>{_h(db)}</td>'
        f'<td>{len(seq):,}</td></tr>'
        for acc, db, seq in records
    )
    return (
        '<h3>Sequences</h3>'
        '<table><thead><tr><th>Accession</th><th>Database</th>'
        '<th>Length (aa)</th></tr></thead>'
        f'<tbody>{rows}</tbody></table>'
    )


def _html_sim_matrix(ids: list[str], sim_mat: list[list[float]],
                      off_mat: list[list[int]], is_aligned: bool) -> str:
    th_cells = "".join(f"<th>{_h(aid)}</th>" for aid in ids)
    rows = []
    for i, row_id in enumerate(ids):
        cells = [f'<td class="bold">{_h(row_id)}</td>']
        for j in range(len(ids)):
            v      = sim_mat[i][j]
            bg     = _sim_color(v)
            val_s  = "1.000" if i == j else f"{v:.3f}"
            off_s  = (f'<br><small style="color:#888">off={off_mat[i][j]}</small>'
                      if not is_aligned and i != j else "")
            cells.append(
                f'<td style="background:{bg}">'
                f'<span class="bold">{val_s}</span>{off_s}</td>'
            )
        rows.append(f'<tr>{"".join(cells)}</tr>')
    note = ("" if is_aligned else
            '<p style="font-size:.82rem;color:#666;margin-top:.3rem">'
            'off = best sliding offset of the shorter sequence</p>')
    return (
        '<h3>Pairwise Similarity Matrix</h3>'
        '<div class="matrix-wrap"><table class="matrix">'
        f'<thead><tr><th></th>{th_cells}</tr></thead>'
        f'<tbody>{"".join(rows)}</tbody>'
        f'</table></div>{note}'
    )


def _html_rank_section(rank_result: dict, tax_map: dict,
                        ids: list[str], rank: str,
                        perm_result: "dict | None" = None) -> str:
    groups    = rank_result["groups"]
    intra_med = rank_result["intra_median"]
    inter_med = rank_result["inter_median"]

    # Sequence → taxon mapping
    seq_rows = []
    for acc in ids:
        info = tax_map.get(acc)
        if info is None:
            seq_rows.append(
                f'<tr><td class="bold">{_h(acc)}</td>'
                f'<td><em>(no taxonomy)</em></td><td>–</td></tr>'
            )
        else:
            node  = info.at_rank(rank)
            rname = _h(node.name) if node else f"<em>(no {_h(rank)} rank)</em>"
            seq_rows.append(
                f'<tr><td class="bold">{_h(acc)}</td>'
                f'<td><em>{_h(info.scientific_name)}</em></td>'
                f'<td>{rname}</td></tr>'
            )

    # Within-group medians
    intra_rows = []
    for taxon, gdata in sorted(groups.items()):
        med = gdata["intra_median"]
        bg  = _sim_color(med) if med is not None else "#e0e0e0"
        intra_rows.append(
            f'<tr><td><em>{_h(taxon)}</em></td>'
            f'<td>{len(gdata["members"])}</td>'
            f'<td>{len(gdata["intra_pairs"])}</td>'
            f'<td style="background:{bg}"><span class="bold">'
            f'{"n/a" if med is None else f"{med:.4f}"}</span></td></tr>'
        )

    # Between-group medians
    inter_by_pair: dict = defaultdict(list)
    for sim, ta, tb, *_ in rank_result["inter_pairs"]:
        inter_by_pair[(ta, tb)].append(sim)
    inter_rows = []
    for (ta, tb), sims in sorted(inter_by_pair.items()):
        med = statistics.median(sims)
        inter_rows.append(
            f'<tr><td><em>{_h(ta)}</em></td><td><em>{_h(tb)}</em></td>'
            f'<td>{len(sims)}</td>'
            f'<td style="background:{_sim_color(med)}">'
            f'<span class="bold">{med:.4f}</span></td></tr>'
        )

    intra_s = f"{intra_med:.4f}" if intra_med is not None else "n/a"
    inter_s = f"{inter_med:.4f}" if inter_med is not None else "n/a"
    diff_note = ""
    if intra_med is not None and inter_med is not None:
        direction = "higher" if intra_med > inter_med else "lower"
        diff_note = (
            f'<p style="margin-top:.5rem;font-size:.88rem">'
            f'Within-{_h(rank)} median is <strong>{direction}</strong> than '
            f'between-{_h(rank)} (Δ = {abs(intra_med - inter_med):.4f})</p>'
        )

    inter_section = ""
    if inter_rows:
        inter_section = (
            f'<h3>Between-{_h(rank)} Median Similarity</h3>'
            f'<table><thead><tr><th>{_h(rank.capitalize())} A</th>'
            f'<th>{_h(rank.capitalize())} B</th>'
            f'<th>Pairs</th><th>Median</th></tr></thead>'
            f'<tbody>{"".join(inter_rows)}</tbody></table>'
        )

    heatmap_img = ""
    hb64 = _heatmap_b64(rank_result)
    if hb64:
        heatmap_img = (
            f'<h3>Similarity Heatmap (by {_h(rank)})</h3>'
            f'<img class="plot" src="data:image/png;base64,{hb64}" '
            f'alt="Rank heatmap">'
        )

    # Permutation test block
    perm_section = ""
    if perm_result is not None:
        reason = perm_result.get("reason", "")
        n_perm = perm_result["n_permutations"]
        if reason:
            perm_section = (
                f'<p style="margin-top:.6rem;font-size:.88rem">'
                f'<strong>Permutation test (n={n_perm}):</strong> {_h(reason)}</p>'
            )
        else:
            obs  = perm_result["observed_stat"]
            pval = perm_result["p_value"]
            pmn  = perm_result["perm_mean"]
            psd  = perm_result["perm_stdev"]
            sig_badge = (' <span style="color:#c00;font-weight:bold">*</span>'
                         if pval is not None and pval < 0.05 else "")
            pval_str  = f"{pval:.4f}" if pval is not None else "n/a"
            obs_str   = f"{obs:+.4f}" if obs is not None else "n/a"
            perm_section = (
                f'<p style="margin-top:.6rem;font-size:.88rem">'
                f'<strong>Permutation test (n={n_perm}):</strong> &nbsp; '
                f'Observed T = <strong>{_h(obs_str)}</strong> &nbsp; '
                f'Null mean ± sd = {_h(f"{pmn:+.4f}")} ± {_h(f"{psd:.4f}")} &nbsp; '
                f'p = <strong>{_h(pval_str)}</strong>{sig_badge}</p>'
            )

    return (
        f'<h3>Taxonomic Conservation — grouped by {_h(rank)}</h3>'
        f'<table><thead><tr><th>Accession</th><th>Organism</th>'
        f'<th>{_h(rank.capitalize())}</th></tr></thead>'
        f'<tbody>{"".join(seq_rows)}</tbody></table>'

        f'<h3>Within-{_h(rank)} Median Similarity</h3>'
        f'<table><thead><tr><th>{_h(rank.capitalize())}</th>'
        f'<th>Members</th><th>Pairs</th><th>Median</th></tr></thead>'
        f'<tbody>{"".join(intra_rows)}</tbody></table>'

        + inter_section
        + (f'<p style="margin-top:.6rem">'
           f'<strong>Overall intra-{_h(rank)} median:</strong> {intra_s} &nbsp; '
           f'<strong>inter-{_h(rank)} median:</strong> {inter_s}</p>')
        + diff_note
        + perm_section
        + heatmap_img
    )


def _html_file_section(summary: dict, file_id: str,
                        idx: int, rank: str | None) -> str:
    name  = _h(summary["file"].name)
    rpt   = summary.get("_report")
    err   = summary.get("error")
    badge = ('<span class="badge badge-err">ERROR</span>'
             if err else '<span class="badge badge-ok">OK</span>')
    header = f'<h2 id="{file_id}">[{idx}] {name} {badge}</h2>'

    if err:
        return (f'<div class="file-section">{header}'
                f'<div class="error-box">{_h(err)}</div></div>')
    if rpt is None:
        return (f'<div class="file-section">{header}'
                f'<div class="error-box">No analysis data available.</div></div>')

    records     = rpt["records"]
    ids         = rpt["ids"]
    sim_mat     = rpt["sim_mat"]
    off_mat     = rpt["off_mat"]
    scores      = rpt["conservation_scores"]
    is_aligned  = rpt["is_aligned"]
    scoring     = rpt["scoring"]
    rank_result = rpt.get("rank_result")
    tax_map     = rpt.get("tax_map") or {}

    scoring_label = scoring.upper() if scoring != "hamming" else "Hamming"
    mode_label    = "aligned (MSA)" if is_aligned else "unaligned (sliding)"
    meta = (f'<p class="meta">'
            f'<strong>{len(records)}</strong> sequences &nbsp;·&nbsp; '
            f'Mode: <strong>{_h(mode_label)}</strong> &nbsp;·&nbsp; '
            f'Scoring: <strong>{_h(scoring_label)}</strong></p>')

    # Conservation plot (or text fallback)
    b64 = _conservation_plot_b64(
        scores, f"Per-Position Conservation – {summary['file'].name}"
    )
    if b64:
        cons_block = (
            '<h3>Per-Position Conservation</h3>'
            f'<img class="plot" src="data:image/png;base64,{b64}" '
            f'alt="Conservation plot">'
        )
    elif scores:
        med = statistics.median(scores)
        cons_block = (
            '<h3>Per-Position Conservation</h3>'
            f'<p>median={med:.4f} &nbsp; min={min(scores):.4f} &nbsp; '
            f'max={max(scores):.4f} '
            f'<em>(install matplotlib for the plot)</em></p>'
        )
    else:
        cons_block = ""

    rank_block = ""
    if rank and rank_result:
        perm_result = rpt.get("perm_result")
        rank_block = _html_rank_section(rank_result, tax_map, ids, rank, perm_result)

    dend_block = ""
    dend_method = rpt.get("dendrogram_method", "average")
    db64 = _dendrogram_b64(ids, sim_mat, method=dend_method)
    if db64:
        dend_block = (
            f'<h3>Hierarchical Clustering Dendrogram '
            f'<span style="font-weight:normal;font-size:.85em">'
            f'({_h(dend_method)} linkage)</span></h3>'
            f'<img class="plot" src="data:image/png;base64,{db64}" '
            f'alt="Dendrogram">'
        )

    return (
        f'<div class="file-section">'
        + header + meta
        + _html_seq_table(records)
        + _html_sim_matrix(ids, sim_mat, off_mat, is_aligned)
        + cons_block
        + dend_block
        + rank_block
        + '</div>'
    )


def _render_html_report(results: list[dict], args,
                         html_path: pathlib.Path) -> None:
    """Write a self-contained HTML report combining all analysed files."""
    import datetime
    run_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    scoring  = getattr(args, "scoring", "hamming")
    rank     = getattr(args, "rank", None)
    mode     = getattr(args, "mode", "auto")
    n_files  = len(results)

    scoring_label = scoring.upper() if scoring != "hamming" else "Hamming"
    rank_note = (f' &nbsp;·&nbsp; Rank: <strong>{_h(rank)}</strong>'
                 if rank else "")

    toc_items = "\n".join(
        f'<li><a href="#file-{i}">{_h(r["file"].name)}'
        + (' &nbsp;<span style="color:#c0392b">✗</span>' if r["error"] else "")
        + "</a></li>"
        for i, r in enumerate(results, 1)
    )
    toc = (
        '<div class="toc"><h2>Contents</h2><ol>'
        '<li><a href="#summary">Summary table</a></li>'
        f'{toc_items}'
        '</ol></div>'
    )

    file_sections = "\n".join(
        _html_file_section(r, f"file-{i}", i, rank)
        for i, r in enumerate(results, 1)
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Protein Conservation Report</title>
<style>
{_REPORT_CSS}
</style>
</head>
<body>
<h1>Protein Conservation Report</h1>
<p class="meta">
  Generated: {run_date} &nbsp;·&nbsp;
  {n_files} file(s) &nbsp;·&nbsp;
  Scoring: <strong>{_h(scoring_label)}</strong> &nbsp;·&nbsp;
  Mode: <strong>{_h(mode)}</strong>{rank_note}
</p>
{toc}
<div id="summary"><h2>Summary</h2>
{_html_summary_table(results, rank)}
</div>
{file_sections}
</body>
</html>"""

    _guard(html_path, getattr(args, "overwrite", False))
    html_path.write_text(html, encoding="utf-8")
    print(f"  HTML report written to: {html_path}")


def _render_report(results: list[dict], args, report_arg: str) -> None:
    """
    Generate a self-contained HTML report and, if a .pdf path is given,
    also convert it via WeasyPrint (pip install weasyprint).

    Extension rules
    ---------------
    .html (or no extension) → HTML file only
    .pdf                    → HTML saved as <stem>.html, then converted to PDF
    """
    p = pathlib.Path(report_arg)
    if not p.suffix:
        p = p.with_suffix(".html")

    want_pdf  = p.suffix.lower() == ".pdf"
    html_path = p.with_suffix(".html") if want_pdf else p

    # Respect --out-dir for relative paths
    if getattr(args, "out_dir", None) and not p.is_absolute():
        base      = pathlib.Path(args.out_dir)
        html_path = base / html_path.name
        if want_pdf:
            p = base / p.name

    html_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        _render_html_report(results, args, html_path)
    except FileExistsError as exc:
        print(f"  Skipped: {exc}", file=sys.stderr)
        return

    if want_pdf:
        try:
            _guard(p, getattr(args, "overwrite", False))
            from weasyprint import HTML as _WPHTML
            _WPHTML(filename=str(html_path)).write_pdf(str(p))
            print(f"  PDF report written to: {p}")
        except ImportError:
            print(
                "  Note: WeasyPrint not installed – PDF output skipped.\n"
                "        pip install weasyprint    (then re-run for PDF output)\n"
                f"        Or open {html_path} in a browser and print-to-PDF.",
                file=sys.stderr,
            )
        except Exception as exc:
            print(f"  PDF generation failed: {exc}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    args = build_parser().parse_args()

    # ── expand input paths ─────────────────────────────────────────────────
    input_files = expand_input_paths(args.fasta, recursive=args.recursive)
    if not input_files:
        sys.exit("Error: no FASTA files found in the specified path(s).")

    multi = len(input_files) > 1
    if multi:
        print(f"Found {len(input_files)} FASTA file(s) to process.")

    # ── shared taxonomy cache ──────────────────────────────────────────────
    cache: TaxCache | None = None
    if args.rank and not args.no_cache and args.cache_file:
        cache = TaxCache(args.cache_file)

    # ── out-dir ────────────────────────────────────────────────────────────
    out_dir: pathlib.Path | None = None
    if args.out_dir:
        out_dir = pathlib.Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

    # ── process each file ──────────────────────────────────────────────────
    results: list[dict] = []
    for idx, fasta_path in enumerate(input_files, 1):
        if multi:
            bar = "━" * 70
            print(f"\n{bar}")
            print(f"[{idx}/{len(input_files)}]  {fasta_path}")
            print(bar)

        stem   = fasta_path.stem
        parent = fasta_path.parent

        # Resolve output paths for this file
        out_matrix = _resolve_out_path(
            args.out_matrix, stem, "_matrix", ".csv", parent, out_dir, multi
        )
        out_cons = _resolve_out_path(
            args.out_conservation, stem, "_conservation", ".csv",
            parent, out_dir, multi
        )
        rank_suffix = f"_rank_{args.rank}" if args.rank else "_rank"
        out_rank = _resolve_out_path(
            args.out_rank_csv, stem, rank_suffix, ".csv", parent, out_dir, multi
        )
        heatmap = _resolve_heatmap_path(
            args.heatmap, stem, parent, out_dir, multi
        )
        dendrogram_p = _resolve_dendrogram_path(
            args.dendrogram, stem, parent, out_dir, multi
        )
        phyloxml_p = _resolve_out_path(
            args.out_phyloxml, stem, "_dendrogram", ".xml", parent, out_dir, multi
        )

        summary = _analyse_fasta(
            args, fasta_path,
            out_matrix, out_cons, out_rank, heatmap, dendrogram_p, phyloxml_p,
            cache,
        )
        results.append(summary)

        if summary["error"] and not multi:
            sys.exit(f"Error: {summary['error']}")

    # ── cross-file summary ─────────────────────────────────────────────────
    if multi:
        _print_multifile_summary(results, args.rank)

    # ── HTML / PDF report ──────────────────────────────────────────────────
    if args.report:
        _render_report(results, args, args.report)

    # Flush cache once after all files
    if cache is not None:
        cache.save()


if __name__ == "__main__":
    main()
