#!/usr/bin/env python3
"""
Protein Pairwise Distance Calculator and Visualizer
Supports:
  • PAM     (40, 80, 120, 160, 200, 250)
  • BLOSUM  (30, 45, 50, 62, 80, 90)
  • Classical (Hamming / p-distance, Jukes-Cantor, Kimura)
  • Cosine  (one-hot, physicochemical, BLOSUM62-embedding, k-mer frequency)
  • Optional TSV label replacement file (2-column: original_id → display_label)
  • Heatmap output as PNG, PDF, SVG, or all formats simultaneously
"""

import argparse
import csv
import datetime
import json
import os
import sys
import warnings
from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Callable

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import (
    dendrogram, linkage, optimal_leaf_ordering, leaves_list, fcluster,
)
from scipy.spatial.distance import squareform, cosine as scipy_cosine
from scipy.stats import pearsonr, spearmanr

warnings.filterwarnings("ignore")

# ══════════════════════════════════════════════════════════════════════════════
#  OUTPUT FORMAT MANAGER
# ══════════════════════════════════════════════════════════════════════════════

SUPPORTED_FORMATS = ["png", "pdf", "svg", "all"]


class FigureExporter:
    """
    Handles saving matplotlib figures to one or more output formats.

    Supported formats : png, pdf, svg, all (= png + pdf + svg)

    PDF behaviour
    ─────────────
    • Single method  → one PDF file per figure.
    • Multi-page PDF → enabled by passing a shared PdfPages object via
                       begin_pdf() / end_pdf(); all subsequent save() calls
                       append pages to the same file.
    • Metadata       → author, creation date, and figure title are embedded
                       in PDF metadata automatically.

    Usage (single figure)
    ─────────────────────
        exporter = FigureExporter(fmt="pdf", dpi=150)
        exporter.save(fig, "output/heatmap", title="My Heatmap")

    Usage (multi-page PDF)
    ──────────────────────
        exporter = FigureExporter(fmt="pdf")
        exporter.begin_multipage_pdf("output/report.pdf")
        exporter.save(fig1, title="Figure 1")
        exporter.save(fig2, title="Figure 2")
        exporter.end_multipage_pdf()
    """

    def __init__(
            self,
            fmt: str = "png",
            dpi: int = 150,
            transparent: bool = False,
            facecolor: str = "white",
    ):
        fmt = fmt.lower()
        if fmt not in SUPPORTED_FORMATS:
            raise ValueError(
                f"Unsupported format '{fmt}'. "
                f"Choose from: {SUPPORTED_FORMATS}"
            )
        self.fmt = fmt
        self.dpi = dpi
        self.transparent = transparent
        self.facecolor = facecolor

        # Multi-page PDF state
        self._pdf_pages: Optional[PdfPages] = None
        self._pdf_path: Optional[str] = None
        self._page_count: int = 0

    # ── multi-page PDF context ────────────────────────────────────────────────

    def begin_multipage_pdf(self, path: str) -> None:
        """Open a PdfPages object for multi-page PDF output."""
        if self._pdf_pages is not None:
            raise RuntimeError("A multi-page PDF is already open. Call end_multipage_pdf() first.")
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        metadata = {
            "Author": "protein_distances pipeline",
            "Creator": "matplotlib / PdfPages",
            "Subject": "Protein pairwise distance analysis",
            "CreationDate": datetime.datetime.now(),
        }
        self._pdf_pages = PdfPages(path, metadata=metadata)
        self._pdf_path = path
        self._page_count = 0
        print(f"  ✓  Multi-page PDF opened → {path}")

    def end_multipage_pdf(self) -> None:
        """Close and finalise the multi-page PDF."""
        if self._pdf_pages is None:
            return
        self._pdf_pages.close()
        print(f"  ✓  Multi-page PDF closed ({self._page_count} page(s)) → {self._pdf_path}")
        self._pdf_pages = None
        self._pdf_path = None
        self._page_count = 0

    # context-manager support
    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.end_multipage_pdf()

    # ── main save method ──────────────────────────────────────────────────────

    def save(
            self,
            fig: plt.Figure,
            base_path: Optional[str] = None,
            title: Optional[str] = None,
    ) -> List[str]:
        """
        Save *fig* to disk.

        Parameters
        ──────────
        fig        : matplotlib Figure to save.
        base_path  : File path WITHOUT extension. Required unless a multi-page
                     PDF is active (in which case the figure is appended there).
        title      : Embedded in PDF page metadata when using PdfPages.

        Returns
        ───────
        List of file paths that were written.
        """
        saved: List[str] = []

        # If multi-page PDF is open, always append there regardless of fmt
        if self._pdf_pages is not None:
            self._pdf_pages.savefig(
                fig,
                dpi=self.dpi,
                bbox_inches="tight",
                facecolor=self.facecolor,
                transparent=self.transparent,
                metadata={"Title": title or ""} if title else None,
            )
            self._page_count += 1
            saved.append(f"{self._pdf_path} (page {self._page_count})")
            return saved

        # ── standalone save ───────────────────────────────────────────────────
        if base_path is None:
            raise ValueError("base_path is required when not using multi-page PDF mode.")

        Path(base_path).parent.mkdir(parents=True, exist_ok=True)

        formats_to_write = (
            ["png", "pdf", "svg"] if self.fmt == "all"
            else [self.fmt]
        )

        for ext in formats_to_write:
            out = f"{base_path}.{ext}"
            kwargs: Dict = dict(
                dpi=self.dpi,
                bbox_inches="tight",
                facecolor=self.facecolor,
                transparent=self.transparent,
            )
            if ext == "pdf":
                # Embed metadata in standalone PDF
                kwargs["metadata"] = self._pdf_metadata(title or Path(base_path).name)
            fig.savefig(out, format=ext, **kwargs)
            saved.append(out)

        return saved

    def save_and_report(
            self,
            fig: plt.Figure,
            base_path: Optional[str],
            label: str,
            title: Optional[str] = None,
    ) -> List[str]:
        """save() + print confirmation lines."""
        paths = self.save(fig, base_path, title=title)
        for p in paths:
            print(f"  ✓  {label} → {p}")
        return paths

    # ── helpers ───────────────────────────────────────────────────────────────

    @staticmethod
    def _pdf_metadata(title: str) -> Dict[str, str]:
        return {
            "Title": title,
            "Author": "protein_distances pipeline",
            "Creator": "matplotlib",
            "Subject": "Protein pairwise distance heatmap",
            "CreationDate": datetime.datetime.now().strftime("%Y%m%d%H%M%S"),
        }

    @property
    def extension(self) -> str:
        """Primary extension (first format if 'all')."""
        return "png" if self.fmt == "all" else self.fmt

    def extensions(self) -> List[str]:
        return ["png", "pdf", "svg"] if self.fmt == "all" else [self.fmt]

    def describe(self) -> str:
        exts = self.extensions()
        return f"Format: {self.fmt.upper()}  →  {' + '.join(e.upper() for e in exts)}"


# ══════════════════════════════════════════════════════════════════════════════
#  TSV LABEL MANAGER
# ══════════════════════════════════════════════════════════════════════════════

class LabelManager:
    """
    Loads an optional two-column TSV file mapping original sequence IDs
    to user-friendly display labels.

    TSV format
    ──────────
    # comment
    original_id<TAB>display_label
    SeqID_001<TAB>Homo sapiens hemoglobin alpha

    Rules
    ─────
    • Lines starting with '#' are skipped.
    • Header row auto-detected if first token is id/seq_id/original/name/label.
    • IDs not present in the TSV keep their original name.
    • Duplicate original IDs → last entry wins (warning printed).
    • BOM characters stripped transparently.
    """

    _HEADER_TOKENS = {"id", "seq_id", "original", "sequence_id", "name", "label", "original_id"}

    def __init__(self, tsv_path: Optional[str] = None, delimiter: str = "\t"):
        self._map: Dict[str, str] = {}
        self._reverse: Dict[str, str] = {}
        self._path = tsv_path
        self._loaded = False

        if tsv_path:
            self._load(tsv_path, delimiter)

    # ── public ────────────────────────────────────────────────────────────────

    def replace(self, original_id: str) -> str:
        return self._map.get(original_id.strip(), original_id.strip())

    def replace_all(self, ids: List[str]) -> List[str]:
        return [self.replace(i) for i in ids]

    def original(self, display_label: str) -> str:
        return self._reverse.get(display_label, display_label)

    @property
    def loaded(self) -> bool:
        return self._loaded

    @property
    def n_mappings(self) -> int:
        return len(self._map)

    def summary(self) -> str:
        if not self._loaded:
            return "  LabelManager: no TSV loaded – using original IDs."
        lines = [f"  LabelManager: {self.n_mappings} label(s) from '{self._path}'"]
        for orig, lbl in self._map.items():
            lines.append(f"    {orig!r:40s} → {lbl!r}")
        return "\n".join(lines)

    def validate_against_alignment(self, alignment: MultipleSeqAlignment) -> None:
        if not self._loaded:
            return
        aln_ids = {rec.id for rec in alignment}
        missing = set(self._map) - aln_ids
        if missing:
            print("  ⚠  TSV IDs not found in alignment (ignored):")
            for m in sorted(missing):
                print(f"       '{m}'")
        matched = len(self._map) - len(missing)
        print(f"  ✓  TSV: {matched}/{len(self._map)} IDs matched.")

    # ── private ───────────────────────────────────────────────────────────────

    def _load(self, path: str, delimiter: str) -> None:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Label TSV not found: '{path}'")

        n_loaded = n_skipped = n_dup = 0
        with open(path, "r", encoding="utf-8-sig") as fh:
            for lineno, raw in enumerate(fh, 1):
                line = raw.rstrip("\n\r")
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                parts = [p.strip() for p in line.split(delimiter)]
                if lineno == 1 and parts[0].lower() in self._HEADER_TOKENS:
                    continue
                if len(parts) < 2:
                    print(f"  ⚠  TSV line {lineno}: <2 columns – skipped.")
                    n_skipped += 1
                    continue
                orig, disp = parts[0], parts[1]
                if orig in self._map:
                    print(f"  ⚠  TSV line {lineno}: duplicate '{orig}' – overwriting.")
                    n_dup += 1
                self._map[orig] = disp
                self._reverse[disp] = orig
                n_loaded += 1

        if n_loaded:
            self._loaded = True
            print(f"  ✓  TSV loaded: {n_loaded} mappings "
                  f"({n_skipped} skipped, {n_dup} duplicates)")
        else:
            print(f"  ⚠  TSV '{path}' contained no valid mappings.")


# ══════════════════════════════════════════════════════════════════════════════
#  AMINO ACID PROPERTIES
# ══════════════════════════════════════════════════════════════════════════════

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}

_RAW_PHYSICOCHEMICAL: Dict[str, List[float]] = {
    "A": [1.8, 88.6, 0.00, 0.0, 0.0, 0.0, 0.0],
    "C": [2.5, 108.5, 1.48, 8.3, 0.0, 0.0, 0.0],
    "D": [-3.5, 111.1, 49.70, 3.9, 0.0, -1.0, 0.0],
    "E": [-3.5, 138.4, 49.90, 4.1, 0.0, -1.0, 0.0],
    "F": [2.8, 189.9, 0.35, 0.0, 0.0, 0.0, 1.0],
    "G": [-0.4, 60.1, 0.00, 0.0, 0.0, 0.0, 0.0],
    "H": [-3.2, 153.2, 51.60, 6.0, 0.0, 0.1, 1.0],
    "I": [4.5, 166.7, 0.15, 0.0, 0.0, 0.0, 0.0],
    "K": [-3.9, 168.6, 49.50, 10.5, 0.0, 1.0, 0.0],
    "L": [3.8, 166.7, 0.45, 0.0, 0.0, 0.0, 0.0],
    "M": [1.9, 162.9, 1.43, 0.0, 0.0, 0.0, 0.0],
    "N": [-3.5, 114.1, 3.38, 0.0, 0.0, 0.0, 0.0],
    "P": [-1.6, 112.7, 1.58, 0.0, 1.0, 0.0, 0.0],
    "Q": [-3.5, 143.8, 3.53, 0.0, 0.0, 0.0, 0.0],
    "R": [-4.5, 173.4, 52.00, 12.5, 0.0, 1.0, 0.0],
    "S": [-0.8, 89.0, 1.67, 0.0, 0.0, 0.0, 0.0],
    "T": [-0.7, 116.1, 1.66, 0.0, 0.0, 0.0, 0.0],
    "V": [4.2, 140.0, 0.13, 0.0, 0.0, 0.0, 0.0],
    "W": [-0.9, 227.8, 2.10, 0.0, 0.0, 0.0, 1.0],
    "Y": [-1.3, 193.6, 1.61, 10.1, 0.0, 0.0, 1.0],
}


def _normalise_columns(raw: Dict[str, List[float]]) -> Dict[str, np.ndarray]:
    mat = np.array([raw[aa] for aa in AA_ORDER], dtype=float)
    lo, hi = mat.min(0), mat.max(0)
    d = hi - lo;
    d[d == 0] = 1.0
    return {aa: (mat[i] - lo) / d for i, aa in enumerate(AA_ORDER)}


PHYSICOCHEMICAL: Dict[str, np.ndarray] = _normalise_columns(_RAW_PHYSICOCHEMICAL)


def _one_hot(aa: str) -> np.ndarray:
    v = np.zeros(20, dtype=float)
    i = AA_INDEX.get(aa.upper(), -1)
    if i >= 0: v[i] = 1.0
    return v


def _build_blosum62_embedding() -> Dict[str, np.ndarray]:
    try:
        mat = substitution_matrices.load("BLOSUM62")
        return {aa: np.array([float(mat[aa][o]) for o in AA_ORDER], dtype=float)
                for aa in AA_ORDER}
    except Exception:
        return {aa: _one_hot(aa) for aa in AA_ORDER}


BLOSUM62_EMBEDDING: Dict[str, np.ndarray] = _build_blosum62_embedding()


# ══════════════════════════════════════════════════════════════════════════════
#  SEQUENCE ENCODER
# ══════════════════════════════════════════════════════════════════════════════

class SequenceEncoder:
    SCHEMES = ["one_hot", "physicochemical", "blosum_embed", "kmer_freq"]

    def __init__(self, scheme: str = "physicochemical", k: int = 2):
        if scheme not in self.SCHEMES:
            raise ValueError(f"Unknown scheme '{scheme}'.")
        self.scheme = scheme
        self.k = k
        self._kmer_vocab: Optional[List[str]] = None

    def encode(self, seq: str) -> np.ndarray:
        seq = seq.upper()
        if self.scheme == "one_hot":
            return self._pos_encode(seq, _one_hot, 20)
        elif self.scheme == "physicochemical":
            return self._pos_encode(seq, lambda aa: PHYSICOCHEMICAL.get(aa, np.zeros(7)), 7)
        elif self.scheme == "blosum_embed":
            return self._pos_encode(seq, lambda aa: BLOSUM62_EMBEDDING.get(aa, np.zeros(20)), 20)
        elif self.scheme == "kmer_freq":
            return self._kmer_encode(seq)

    def build_vocab(self, sequences: List[str]) -> None:
        kmers: set = set()
        for s in sequences:
            clean = s.upper().replace("-", "")
            for i in range(len(clean) - self.k + 1):
                kmers.add(clean[i:i + self.k])
        self._kmer_vocab = sorted(kmers)

    @staticmethod
    def _pos_encode(seq: str, fn: Callable, dim: int) -> np.ndarray:
        return np.concatenate([np.zeros(dim) if aa == "-" else fn(aa) for aa in seq])

    def _kmer_encode(self, seq: str) -> np.ndarray:
        clean = seq.replace("-", "")
        if self._kmer_vocab is None:
            self.build_vocab([seq])
        counts = Counter(clean[i:i + self.k] for i in range(len(clean) - self.k + 1))
        total = max(sum(counts.values()), 1)
        return np.array([counts.get(km, 0) / total for km in self._kmer_vocab], dtype=float)


# ══════════════════════════════════════════════════════════════════════════════
#  COSINE DISTANCE ENGINE
# ══════════════════════════════════════════════════════════════════════════════

class CosineDistanceEngine:
    def __init__(self, scheme: str = "physicochemical", k: int = 2, normalise: bool = True):
        self.scheme = scheme
        self.k = k
        self.normalise = normalise
        self.encoder = SequenceEncoder(scheme=scheme, k=k)

    def compute_matrix(self, alignment: MultipleSeqAlignment) -> np.ndarray:
        seqs = [str(rec.seq) for rec in alignment]
        if self.scheme == "kmer_freq":
            self.encoder.build_vocab(seqs)
        vectors = np.array([self.encoder.encode(s) for s in seqs], dtype=float)
        if self.normalise and self.scheme != "kmer_freq":
            norms = np.linalg.norm(vectors, axis=1, keepdims=True)
            norms[norms == 0] = 1.0
            vectors /= norms
        n = len(vectors)
        D = np.zeros((n, n), dtype=float)
        for i, j in combinations(range(n), 2):
            u, v = vectors[i], vectors[j]
            nu, nv = np.linalg.norm(u), np.linalg.norm(v)
            d = 1.0 if (nu == 0 or nv == 0) else float(scipy_cosine(u, v))
            D[i, j] = D[j, i] = max(0.0, min(1.0, d))
        return D


# ══════════════════════════════════════════════════════════════════════════════
#  METHOD REGISTRY
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class MatrixInfo:
    name: str
    family: str
    biopython_key: str
    description: str
    typical_identity: str


MATRIX_REGISTRY: Dict[str, MatrixInfo] = {
    "PAM40": MatrixInfo("PAM40", "PAM", "PAM40", "PAM40  – very close (~85–90% id)", "~85–90%"),
    "PAM80": MatrixInfo("PAM80", "PAM", "PAM80", "PAM80  – close (~75–85% id)", "~75–85%"),
    "PAM120": MatrixInfo("PAM120", "PAM", "PAM120", "PAM120 – moderate (~65–75% id)", "~65–75%"),
    "PAM160": MatrixInfo("PAM160", "PAM", "PAM160", "PAM160 – distant (~50–65% id)", "~50–65%"),
    "PAM200": MatrixInfo("PAM200", "PAM", "PAM200", "PAM200 – distant (~40–55% id)", "~40–55%"),
    "PAM250": MatrixInfo("PAM250", "PAM", "PAM250", "PAM250 – very distant (<40% id)", "<40%"),
    "BLOSUM30": MatrixInfo("BLOSUM30", "BLOSUM", "BLOSUM30", "BLOSUM30 – very distant (<30% id)", "<30%"),
    "BLOSUM45": MatrixInfo("BLOSUM45", "BLOSUM", "BLOSUM45", "BLOSUM45 – distant (~30–40% id)", "~30–40%"),
    "BLOSUM50": MatrixInfo("BLOSUM50", "BLOSUM", "BLOSUM50", "BLOSUM50 – moderate (~40–50% id)", "~40–50%"),
    "BLOSUM62": MatrixInfo("BLOSUM62", "BLOSUM", "BLOSUM62", "BLOSUM62 – general (~50–62% id)", "~50–62%"),
    "BLOSUM80": MatrixInfo("BLOSUM80", "BLOSUM", "BLOSUM80", "BLOSUM80 – close (~80% id)", "~80%"),
    "BLOSUM90": MatrixInfo("BLOSUM90", "BLOSUM", "BLOSUM90", "BLOSUM90 – very close (>90% id)", ">90%"),
    "hamming": MatrixInfo("hamming", "classical", "", "p-distance", "any"),
    "jukes_cantor": MatrixInfo("jukes_cantor", "classical", "", "Jukes-Cantor 20-state", "any"),
    "kimura": MatrixInfo("kimura", "classical", "", "Kimura empirical distance", "any"),
    "cosine_onehot": MatrixInfo("cosine_onehot", "cosine", "", "Cosine – one-hot positions", "any"),
    "cosine_physico": MatrixInfo("cosine_physico", "cosine", "", "Cosine – physicochemical (7D)", "any"),
    "cosine_blosum": MatrixInfo("cosine_blosum", "cosine", "", "Cosine – BLOSUM62 embedding", "any"),
    "cosine_kmer": MatrixInfo("cosine_kmer", "cosine", "", "Cosine – k-mer frequencies", "any"),
}

PAM_METHODS = [k for k, v in MATRIX_REGISTRY.items() if v.family == "PAM"]
BLOSUM_METHODS = [k for k, v in MATRIX_REGISTRY.items() if v.family == "BLOSUM"]
CLASSICAL_METHODS = [k for k, v in MATRIX_REGISTRY.items() if v.family == "classical"]
COSINE_METHODS = [k for k, v in MATRIX_REGISTRY.items() if v.family == "cosine"]
ALL_METHODS = list(MATRIX_REGISTRY.keys())

COSINE_SCHEME_MAP = {
    "cosine_onehot": "one_hot",
    "cosine_physico": "physicochemical",
    "cosine_blosum": "blosum_embed",
    "cosine_kmer": "kmer_freq",
}

FAMILY_COLORS = {
    "PAM": "#E07B39",
    "BLOSUM": "#4A90D9",
    "classical": "#5BAD72",
    "cosine": "#9B59B6",
}

COLORMAPS = {
    "PAM": "YlOrRd",
    "BLOSUM": "YlGnBu",
    "classical": "BuGn",
    "cosine": "PuRd",
    "default": "viridis",
}


def family_of(m: str) -> str:
    return MATRIX_REGISTRY.get(m, MatrixInfo("?", "classical", "", "", "")).family


# ══════════════════════════════════════════════════════════════════════════════
#  SUBSTITUTION MATRIX CACHE
# ══════════════════════════════════════════════════════════════════════════════

class SubstitutionMatrixCache:
    _cache: Dict[str, object] = {}

    @classmethod
    def get(cls, name: str):
        if name not in cls._cache:
            cls._cache[name] = substitution_matrices.load(name)
        return cls._cache[name]

    @classmethod
    def score(cls, name: str, aa1: str, aa2: str) -> float:
        try:
            return float(cls.get(name)[aa1.upper()][aa2.upper()])
        except (KeyError, IndexError):
            return 0.0


# ══════════════════════════════════════════════════════════════════════════════
#  CLASSICAL DISTANCE FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def _comparable_pairs(s1: str, s2: str):
    return [(a, b) for a, b in zip(s1, s2) if not (a == "-" and b == "-")]


def hamming_distance(s1: str, s2: str) -> float:
    p = _comparable_pairs(str(s1), str(s2))
    return 0.0 if not p else sum(1 for a, b in p if a != b) / len(p)


def jukes_cantor_distance(s1: str, s2: str, n: int = 20) -> float:
    p = hamming_distance(s1, s2)
    v = 1.0 - (n / (n - 1)) * p
    return np.inf if v <= 0 else -((n - 1) / n) * np.log(v)


def kimura_distance(s1: str, s2: str) -> float:
    p = hamming_distance(s1, s2)
    v = 1.0 - p - 0.2 * p ** 2
    return np.inf if v <= 0 else -np.log(v)


def substitution_matrix_distance(
        s1: str, s2: str, matrix_name: str, metric: str = "normalized_score"
) -> float:
    pairs = _comparable_pairs(str(s1), str(s2))
    if not pairs: return 1.0
    if metric == "normalized_score":
        obs = sum(SubstitutionMatrixCache.score(matrix_name, a, b) for a, b in pairs)
        ss1 = sum(SubstitutionMatrixCache.score(matrix_name, a, a) for a, b in pairs)
        ss2 = sum(SubstitutionMatrixCache.score(matrix_name, b, b) for a, b in pairs)
        den = (ss1 + ss2) / 2.0
        return max(0.0, 1.0 - obs / den) if den else 1.0
    elif metric == "log_odds":
        sc = [SubstitutionMatrixCache.score(matrix_name, a, b) for a, b in pairs]
        slf = [SubstitutionMatrixCache.score(matrix_name, a, a) for a, b in pairs]
        mx = max(slf) if slf else 1.0
        return max(0.0, 1.0 - np.mean(sc) / mx) if mx else 1.0
    elif metric == "percent_similarity":
        pos = sum(1 for a, b in pairs if SubstitutionMatrixCache.score(matrix_name, a, b) > 0)
        return 1.0 - pos / len(pairs)
    else:
        raise ValueError(f"Unknown sub_metric '{metric}'")


# ══════════════════════════════════════════════════════════════════════════════
#  DISTANCE MATRIX BUILDER
# ══════════════════════════════════════════════════════════════════════════════

def calculate_distance_matrix(
        alignment: MultipleSeqAlignment,
        method: str,
        sub_metric: str = "normalized_score",
        cosine_k: int = 2,
) -> np.ndarray:
    n = len(alignment)

    if method in COSINE_METHODS:
        scheme = COSINE_SCHEME_MAP[method]
        print(f"\n  Cosine | encoding='{scheme}'"
              + (f" k={cosine_k}" if scheme == "kmer_freq" else ""))
        D = CosineDistanceEngine(scheme=scheme, k=cosine_k).compute_matrix(alignment)
        up = D[np.triu_indices(n, k=1)]
        print(f"  ✓  min={up.min():.4f}  max={up.max():.4f}  mean={up.mean():.4f}")
        return D

    total = n * (n - 1) // 2
    print(f"\n  Computing {total} pairs  [method={method}] …")
    D = np.zeros((n, n), dtype=float)
    for cnt, (i, j) in enumerate(combinations(range(n), 2), 1):
        if method == "hamming":
            d = hamming_distance(str(alignment[i].seq), str(alignment[j].seq))
        elif method == "jukes_cantor":
            d = jukes_cantor_distance(str(alignment[i].seq), str(alignment[j].seq))
        elif method == "kimura":
            d = kimura_distance(str(alignment[i].seq), str(alignment[j].seq))
        else:
            key = MATRIX_REGISTRY[method].biopython_key
            d = substitution_matrix_distance(
                str(alignment[i].seq), str(alignment[j].seq), key, sub_metric)
        D[i, j] = D[j, i] = np.nan if (np.isinf(d) or np.isnan(d)) else d
        if cnt % max(1, total // 10) == 0:
            print(f"    {cnt}/{total}  ({100 * cnt / total:.0f}%)")

    finite = D[np.isfinite(D) & (D > 0)]
    if len(finite):
        fmax = finite.max()
        nan_n = int(np.isnan(D).sum()) // 2
        if nan_n:
            print(f"  ⚠  {nan_n} saturated → clamped to {fmax:.4f}")
            D = np.where(np.isnan(D), fmax, D)
        up = D[np.triu_indices(n, k=1)]
        print(f"  ✓  min={up.min():.4f}  max={up.max():.4f}  mean={up.mean():.4f}")
    return D


# ══════════════════════════════════════════════════════════════════════════════
#  MSA PARSING
# ══════════════════════════════════════════════════════════════════════════════

FORMATS = ["fasta", "clustal", "phylip", "stockholm", "nexus", "phylip-relaxed"]


def parse_alignment(filepath: str, fmt: Optional[str] = None) -> MultipleSeqAlignment:
    order = ([fmt] + [f for f in FORMATS if f != fmt]) if fmt else FORMATS
    for f in order:
        try:
            aln = AlignIO.read(filepath, f)
            print(f"  ✓  Parsed as '{f}'")
            return aln
        except Exception:
            continue
    raise ValueError(f"Cannot parse '{filepath}'. Tried: {order}")


def print_alignment_summary(
        aln: MultipleSeqAlignment,
        label_mgr: LabelManager,
) -> Tuple[int, int]:
    n, L = len(aln), aln.get_alignment_length()
    print(f"\n{'─' * 70}")
    print(f"  Sequences  : {n}   |   Alignment length : {L} columns")
    print(f"  {'Original ID':<35}  {'Display Label':<28}  Gaps%")
    print(f"  {'─' * 35}  {'─' * 28}  {'─' * 5}")
    for rec in aln:
        gp = str(rec.seq).count("-") / L * 100
        label = label_mgr.replace(rec.id)
        marker = " ◀" if label != rec.id else ""
        print(f"  {rec.id[:34]:<35}  {label[:27]:<28}  {gp:5.1f}%{marker}")
    print(f"{'─' * 70}\n")
    if n < 2: raise ValueError("Need ≥ 2 sequences.")
    return n, L


# ══════════════════════════════════════════════════════════════════════════════
#  VISUAL HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def shorten(s: str, n: int = 28) -> str:
    return s if len(s) <= n else s[:n - 1] + "…"


def pick_text_color(rgba) -> str:
    r, g, b = rgba[:3]
    return "white" if 0.299 * r + 0.587 * g + 0.114 * b < 0.55 else "black"


def _draw_cluster_boxes(ax, order, cluster_map):
    n = len(order)
    reord = [cluster_map[order[i]] for i in range(n)]
    unique = sorted(set(reord))
    palette = plt.cm.get_cmap("Set2")(np.linspace(0, 1, max(len(unique), 1)))
    for ci, cl_id in enumerate(unique):
        idxs = [i for i, c in enumerate(reord) if c == cl_id]
        if not idxs: continue
        lo, hi = min(idxs), max(idxs)
        ax.add_patch(plt.Rectangle(
            (lo - 0.5, lo - 0.5), hi - lo + 1, hi - lo + 1,
            fill=False, edgecolor=palette[ci], linewidth=2.5, linestyle="--"))


# ══════════════════════════════════════════════════════════════════════════════
#  PLOT 1 – HEATMAP + DENDROGRAM
# ══════════════════════════════════════════════════════════════════════════════

def plot_heatmap_with_dendrogram(
        dist_matrix: np.ndarray,
        labels: List[str],
        original_ids: List[str],
        method: str = "BLOSUM62",
        linkage_method: str = "average",
        colormap: Optional[str] = None,
        exporter: Optional[FigureExporter] = None,
        base_path: Optional[str] = None,
        title: Optional[str] = None,
        figsize: Optional[Tuple] = None,
        n_clusters: int = 0,
) -> Tuple[plt.Figure, np.ndarray, List[int]]:
    n = len(labels)
    short = [shorten(lb) for lb in labels]
    fam = family_of(method)
    cmap = colormap or COLORMAPS.get(fam, "viridis")

    # Clustering
    condensed = np.clip(squareform(dist_matrix), 0, None)
    Z = linkage(condensed, method=linkage_method)
    try:
        Z = optimal_leaf_ordering(Z, condensed)
    except Exception:
        pass
    order = leaves_list(Z)
    reord = dist_matrix[np.ix_(order, order)]
    reord_labels = [short[i] for i in order]
    reord_origs = [original_ids[i] for i in order]

    thr = 0.3 * Z[:, 2].max()
    cl_ids = (fcluster(Z, t=n_clusters, criterion="maxclust")
              if n_clusters > 0
              else fcluster(Z, t=thr, criterion="distance"))
    cluster_map = {order[i]: cl_ids[i] for i in range(n)}

    # Figure
    auto = max(10, n * 0.6)
    fig = plt.figure(figsize=figsize or (auto + 5, auto + 2.5), facecolor="white")
    gs = GridSpec(3, 4, figure=fig,
                  width_ratios=[0.9, 0.12, 5.5, 0.35],
                  height_ratios=[0.9, 0.10, 5.5],
                  hspace=0.02, wspace=0.02)
    ax_td = fig.add_subplot(gs[0, 2])
    ax_ld = fig.add_subplot(gs[2, 0])
    ax_hm = fig.add_subplot(gs[2, 2])
    ax_cb = fig.add_subplot(gs[2, 3])
    ax_in = fig.add_subplot(gs[0, 0]);
    ax_in.axis("off")

    dend_thr = 0.65 * Z[:, 2].max()

    dendrogram(Z, ax=ax_td, color_threshold=dend_thr,
               above_threshold_color="#888", orientation="top", no_labels=True)
    ax_td.set_xlim(0, n * 10);
    ax_td.axis("off");
    ax_td.set_facecolor("#fafafa")

    dendrogram(Z, ax=ax_ld, color_threshold=dend_thr,
               above_threshold_color="#888", orientation="left", no_labels=True)
    ax_ld.set_ylim(0, n * 10);
    ax_ld.invert_yaxis()
    ax_ld.axis("off");
    ax_ld.set_facecolor("#fafafa")

    cmap_obj = plt.cm.get_cmap(cmap)
    vmax = reord.max() or 1.0
    im = ax_hm.imshow(reord, aspect="auto", cmap=cmap_obj,
                      vmin=0, vmax=vmax, interpolation="nearest")

    fs = max(5, min(11, 130 // n))
    ax_hm.set_xticks(range(n));
    ax_hm.set_yticks(range(n))
    ax_hm.set_xticklabels(reord_labels, rotation=45, ha="right",
                          fontsize=fs, fontfamily="monospace")
    ax_hm.set_yticklabels(reord_labels, fontsize=fs, fontfamily="monospace")
    ax_hm.set_xticks(np.arange(-0.5, n), minor=True)
    ax_hm.set_yticks(np.arange(-0.5, n), minor=True)
    ax_hm.grid(which="minor", color="white", linewidth=0.3, alpha=0.5)
    ax_hm.tick_params(which="minor", length=0);
    ax_hm.tick_params(length=2, pad=2)

    if n <= 25:
        cfs = max(4, min(8, 85 // n))
        for ii in range(n):
            for jj in range(n):
                v = reord[ii, jj]
                rgba = cmap_obj(v / vmax)
                ax_hm.text(jj, ii, f"{v:.3f}", ha="center", va="center",
                           fontsize=cfs, color=pick_text_color(rgba),
                           fontweight="bold" if ii == jj else "normal")

    for k in range(n):
        ax_hm.add_patch(plt.Rectangle((k - 0.5, k - 0.5), 1, 1,
                                      fill=False, edgecolor="#555", linewidth=1.2))
    _draw_cluster_boxes(ax_hm, order, cluster_map)

    cbar = plt.colorbar(im, cax=ax_cb)
    cbar.set_label("Distance", fontsize=10, labelpad=6)
    cbar.ax.tick_params(labelsize=8)

    info = MATRIX_REGISTRY.get(method)
    lbl_changed = any(l != o for l, o in zip(labels, original_ids))
    lbl_note = "custom labels ✓" if lbl_changed else "original IDs"
    # Include output format in info box
    fmt_note = exporter.describe() if exporter else "PNG"
    txt = (f"Method  : {method}\n"
           f"Family  : {fam}\n"
           f"Linkage : {linkage_method}\n"
           f"n seqs  : {n}\n"
           f"Labels  : {lbl_note}\n"
           f"Output  : {fmt_note}\n"
           + (f"Id range: {info.typical_identity}" if info else ""))
    ax_in.text(0.05, 0.95, txt, transform=ax_in.transAxes,
               fontsize=7.5, va="top", fontfamily="monospace",
               bbox=dict(boxstyle="round", facecolor="#f0f4ff",
                         edgecolor=FAMILY_COLORS.get(fam, "#aaa"), alpha=0.9))

    ttl = title or (f"Protein Pairwise Distance  ·  {method}  ·  "
                    f"Linkage: {linkage_method}  ·  n={n}")
    fig.suptitle(ttl, fontsize=13, fontweight="bold", y=0.998)
    plt.tight_layout(rect=[0, 0, 1, 0.995])

    if exporter and base_path:
        exporter.save_and_report(fig, base_path, label="Heatmap", title=ttl)

    return fig, reord, list(order)


# ══════════════════════════════════════════════════════════════════════════════
#  PLOT 2 – COSINE ENCODING COMPARISON
# ══════════════════════════════════════════════════════════════════════════════

def plot_cosine_encoding_comparison(
        alignment: MultipleSeqAlignment,
        labels: List[str],
        exporter: Optional[FigureExporter] = None,
        base_path: Optional[str] = None,
) -> plt.Figure:
    schemes = [
        ("one_hot", "One-Hot (20D)"),
        ("physicochemical", "Physicochemical (7D)"),
        ("blosum_embed", "BLOSUM62 Embedding (20D)"),
        ("kmer_freq", "k-mer Frequency (k=2)"),
    ]
    n = len(labels)
    short = [shorten(lb, 18) for lb in labels]
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), facecolor="white")
    axes = axes.flatten()
    for ax, (scheme, title) in zip(axes, schemes):
        D = CosineDistanceEngine(scheme=scheme, k=2).compute_matrix(alignment)
        im = ax.imshow(D, cmap="PuRd", vmin=0, vmax=1,
                       aspect="auto", interpolation="nearest")
        fs = max(5, min(9, 120 // n))
        ax.set_xticks(range(n));
        ax.set_yticks(range(n))
        ax.set_xticklabels(short, rotation=45, ha="right", fontsize=fs, fontfamily="monospace")
        ax.set_yticklabels(short, fontsize=fs, fontfamily="monospace")
        ax.set_title(title, fontsize=11, fontweight="bold", color="#9B59B6")
        plt.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
        if n <= 15:
            for ii in range(n):
                for jj in range(n):
                    rgba = plt.cm.PuRd(D[ii, jj])
                    ax.text(jj, ii, f"{D[ii, jj]:.2f}", ha="center", va="center",
                            fontsize=max(4, 55 // n), color=pick_text_color(rgba))
    ttl = "Cosine Distance: Encoding Scheme Comparison"
    fig.suptitle(ttl, fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    if exporter and base_path:
        exporter.save_and_report(fig, base_path, label="Cosine schemes", title=ttl)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  PLOT 3 – MDS MAP
# ══════════════════════════════════════════════════════════════════════════════

def plot_mds_sequence_map(
        all_matrices: Dict[str, np.ndarray],
        labels: List[str],
        methods: List[str],
        exporter: Optional[FigureExporter] = None,
        base_path: Optional[str] = None,
) -> plt.Figure:
    try:
        from sklearn.manifold import MDS
        def mds2d(D):
            return MDS(n_components=2, dissimilarity="precomputed",
                       random_state=42, n_init=4).fit_transform(D)
    except ImportError:
        def mds2d(D):
            n = D.shape[0];
            H = np.eye(n) - np.ones((n, n)) / n
            B = -0.5 * H @ (D ** 2) @ H
            vals, vecs = np.linalg.eigh(B)
            idx = np.argsort(-vals)
            return vecs[:, idx][:, :2] * np.sqrt(np.maximum(vals[idx][:2], 0))

    n_m = len(methods);
    short = [shorten(lb, 14) for lb in labels]
    fig, axes = plt.subplots(1, n_m, figsize=(5.5 * n_m, 5), facecolor="white", squeeze=False)
    for ax, method in zip(axes[0], methods):
        fam = family_of(method);
        D = all_matrices[method]
        coords = mds2d(np.clip(D, 0, None))
        Z = linkage(squareform(np.clip(D, 0, None)), method="average")
        cl = fcluster(Z, t=0.3 * Z[:, 2].max(), criterion="distance")
        col = plt.cm.tab10(np.linspace(0, 1, max(cl)))
        for i, (x, y) in enumerate(coords):
            ax.scatter(x, y, color=col[cl[i] - 1], s=120, zorder=3,
                       edgecolors="white", linewidths=0.8)
            ax.annotate(short[i], (x, y), textcoords="offset points",
                        xytext=(5, 5), fontsize=7, fontfamily="monospace", color="#333")
        ax.set_title(method, fontsize=11, fontweight="bold",
                     color=FAMILY_COLORS.get(fam, "#333"))
        ax.set_xlabel("MDS-1", fontsize=9);
        ax.set_ylabel("MDS-2", fontsize=9)
        ax.grid(alpha=0.25);
        ax.set_facecolor("#fafafa")
    ttl = "2-D MDS Sequence Map"
    fig.suptitle(ttl, fontsize=13, fontweight="bold")
    plt.tight_layout()
    if exporter and base_path:
        exporter.save_and_report(fig, base_path, label="MDS map", title=ttl)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  PLOT 4 – MULTI-METHOD COMPARISON
# ══════════════════════════════════════════════════════════════════════════════

def plot_multi_method_comparison(
        matrices: Dict[str, np.ndarray],
        labels: List[str],
        exporter: Optional[FigureExporter] = None,
        base_path: Optional[str] = None,
) -> plt.Figure:
    methods = list(matrices.keys());
    n_m = len(methods)
    n = len(labels);
    short = [shorten(lb, 16) for lb in labels]
    cell = max(3.0, 0.35 * n)
    fig, axes = plt.subplots(1, n_m, figsize=(cell * n_m + 2, cell + 2.5), facecolor="#f8f8f8")
    if n_m == 1: axes = [axes]
    for ax, method in zip(axes, methods):
        mat = matrices[method];
        fam = family_of(method)
        cmap = COLORMAPS.get(fam, "viridis")
        im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=mat.max(),
                       aspect="auto", interpolation="nearest")
        ax.set_title(method, fontsize=9, fontweight="bold",
                     color=FAMILY_COLORS.get(fam, "#333"), pad=4)
        fs = max(4, min(8, 100 // n))
        ax.set_xticks(range(n))
        ax.set_xticklabels(short, rotation=45, ha="right", fontsize=fs, fontfamily="monospace")
        ax.set_yticks(range(n))
        ax.set_yticklabels(short, fontsize=fs, fontfamily="monospace")
        if n <= 15:
            cm = plt.cm.get_cmap(cmap);
            vm = mat.max() or 1
            for ii in range(n):
                for jj in range(n):
                    rgba = cm(mat[ii, jj] / vm)
                    ax.text(jj, ii, f"{mat[ii, jj]:.2f}", ha="center", va="center",
                            fontsize=max(3, 55 // n), color=pick_text_color(rgba))
        plt.colorbar(im, ax=ax, shrink=0.8, pad=0.03)
    handles = [mpatches.Patch(color=c, label=fam)
               for fam, c in FAMILY_COLORS.items()
               if any(family_of(m) == fam for m in methods)]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               fontsize=8, framealpha=0.7, bbox_to_anchor=(0.5, -0.01))
    ttl = "Multi-Method Distance Comparison"
    fig.suptitle(ttl, fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    if exporter and base_path:
        exporter.save_and_report(fig, base_path, label="Method comparison", title=ttl)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  PLOT 5 – DISTANCE DISTRIBUTIONS
# ══════════════════════════════════════════════════════════════════════════════

def plot_distance_distribution(
        matrices: Dict[str, np.ndarray],
        labels: List[str],
        exporter: Optional[FigureExporter] = None,
        base_path: Optional[str] = None,
) -> plt.Figure:
    n = len(labels);
    methods = list(matrices.keys())
    short = [shorten(lb, 16) for lb in labels]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), facecolor="white")
    ax = axes[0]
    for method, mat in matrices.items():
        fam = family_of(method);
        col = FAMILY_COLORS.get(fam, "steelblue")
        up = mat[np.triu_indices(n, k=1)]
        ax.hist(up, bins=min(30, max(5, len(up) // 3)), alpha=0.45,
                color=col, edgecolor="white", label=f"{method} (μ={up.mean():.3f})")
        ax.axvline(up.mean(), color=col, linestyle="--", linewidth=1.5)
    ax.set_xlabel("Distance", fontsize=11);
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("Pairwise Distance Distribution", fontsize=12, fontweight="bold")
    ax.legend(fontsize=7, framealpha=0.85);
    ax.grid(alpha=0.3)
    ax2 = axes[1];
    x = np.arange(n);
    bw = 0.8 / max(len(methods), 1)
    for mi, (method, mat) in enumerate(matrices.items()):
        fam = family_of(method);
        col = FAMILY_COLORS.get(fam, "steelblue")
        means = np.array([np.mean(np.concatenate([mat[i, :i], mat[i, i + 1:]])) for i in range(n)])
        off = (mi - len(methods) / 2 + 0.5) * bw
        ax2.bar(x + off, means, width=bw * 0.88, color=col, alpha=0.72, label=method)
    fs = max(6, min(10, 110 // n))
    ax2.set_xticks(x)
    ax2.set_xticklabels(short, rotation=45, ha="right", fontsize=fs, fontfamily="monospace")
    ax2.set_ylabel("Mean Distance to Others", fontsize=11)
    ax2.set_title("Per-Sequence Mean Distance", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=7, framealpha=0.85);
    ax2.grid(alpha=0.3, axis="y")
    plt.tight_layout()
    ttl = "Distance Distributions"
    if exporter and base_path:
        exporter.save_and_report(fig, base_path, label="Distributions", title=ttl)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  PLOT 6 – METHOD CORRELATION
# ══════════════════════════════════════════════════════════════════════════════

def plot_method_correlation(
        matrices: Dict[str, np.ndarray],
        exporter: Optional[FigureExporter] = None,
        base_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    methods = list(matrices.keys());
    n_m = len(methods)
    if n_m < 2: print("  ⚠  Need ≥ 2 methods for correlation."); return None
    n_seq = list(matrices.values())[0].shape[0]
    idx = np.triu_indices(n_seq, k=1)
    vecs = {m: matrices[m][idx] for m in methods}
    fig, axes = plt.subplots(n_m, n_m, figsize=(3.5 * n_m, 3.5 * n_m), facecolor="white")
    for i, mi in enumerate(methods):
        for j, mj in enumerate(methods):
            ax = axes[i][j];
            di, dj = vecs[mi], vecs[mj]
            fi = family_of(mi)
            if i == j:
                ax.hist(di, bins=20, color=FAMILY_COLORS.get(fi, "steelblue"),
                        edgecolor="white", alpha=0.8)
                ax.set_title(mi, fontsize=8, fontweight="bold",
                             color=FAMILY_COLORS.get(fi, "#333"))
            elif i < j:
                ax.scatter(dj, di, s=7, alpha=0.45, color=FAMILY_COLORS.get(fi, "steelblue"))
                r, _ = pearsonr(di, dj);
                rs, _ = spearmanr(di, dj)
                ax.set_title(f"r={r:.3f}  ρ={rs:.3f}", fontsize=7)
            else:
                res = di - dj
                ax.axhline(0, color="red", linewidth=1, linestyle="--")
                ax.scatter(dj, res, s=5, alpha=0.35, color=FAMILY_COLORS.get(fi, "gray"))
                ax.set_title(f"{mi}−{mj}", fontsize=7)
            if j == 0: ax.set_ylabel(mi, fontsize=7)
            if i == n_m - 1: ax.set_xlabel(mj, fontsize=7)
            ax.tick_params(labelsize=5)
    ttl = "Method Correlation Matrix"
    fig.suptitle(ttl, fontsize=10, fontweight="bold")
    plt.tight_layout()
    if exporter and base_path:
        exporter.save_and_report(fig, base_path, label="Correlation matrix", title=ttl)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  MULTI-PAGE PDF REPORT
# ══════════════════════════════════════════════════════════════════════════════

def build_multipage_pdf_report(
        all_figures: List[Tuple[plt.Figure, str]],  # [(fig, title), …]
        pdf_path: str,
        pipeline_meta: Dict,
) -> None:
    """
    Combine every figure produced by the pipeline into a single, structured
    multi-page PDF report with a cover page.

    Parameters
    ──────────
    all_figures   : list of (figure, page_title) pairs in order.
    pdf_path      : destination path for the combined PDF.
    pipeline_meta : dict with keys: input_file, methods, linkage, n_seqs,
                    label_tsv, output_prefix, timestamp.
    """
    Path(pdf_path).parent.mkdir(parents=True, exist_ok=True)

    metadata = {
        "Title": "Protein Distance Analysis Report",
        "Author": "protein_distances pipeline",
        "Subject": f"Analysis of {pipeline_meta.get('input_file', 'unknown')}",
        "CreationDate": datetime.datetime.now(),
    }

    with PdfPages(pdf_path, metadata=metadata) as pdf:

        # ── Cover page ────────────────────────────────────────────────────────
        cover = plt.figure(figsize=(11, 8.5), facecolor="white")
        ax = cover.add_subplot(111);
        ax.axis("off")

        # Title block
        ax.text(0.5, 0.88, "Protein Pairwise Distance Analysis",
                ha="center", va="center", fontsize=22, fontweight="bold",
                transform=ax.transAxes, color="#1a1a2e")
        ax.text(0.5, 0.80, "Automated Report",
                ha="center", va="center", fontsize=14,
                transform=ax.transAxes, color="#4a4a6a")

        # Horizontal rule
        ax.axhline(0.76, xmin=0.05, xmax=0.95, color="#4A90D9", linewidth=2)

        # Metadata table
        ts = pipeline_meta.get("timestamp", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        rows = [
            ("Input file", pipeline_meta.get("input_file", "–")),
            ("Label TSV", pipeline_meta.get("label_tsv", "(none)")),
            ("Sequences", str(pipeline_meta.get("n_seqs", "–"))),
            ("Methods", ", ".join(pipeline_meta.get("methods", []))),
            ("Linkage", pipeline_meta.get("linkage", "–")),
            ("Sub-metric", pipeline_meta.get("sub_metric", "–")),
            ("Output prefix", pipeline_meta.get("output_prefix", "–")),
            ("Figures", str(len(all_figures))),
            ("Generated", ts),
        ]
        y = 0.70
        for key, val in rows:
            ax.text(0.10, y, f"{key}:", ha="left", va="center",
                    fontsize=10, fontweight="bold", transform=ax.transAxes, color="#333")
            ax.text(0.38, y, val, ha="left", va="center",
                    fontsize=10, transform=ax.transAxes, color="#555",
                    fontfamily="monospace")
            y -= 0.062

        # Figure index
        ax.axhline(y - 0.01, xmin=0.05, xmax=0.95, color="#ccc", linewidth=1)
        y -= 0.05
        ax.text(0.10, y, "Contents:", ha="left", va="center",
                fontsize=10, fontweight="bold", transform=ax.transAxes, color="#333")
        y -= 0.05
        for page_n, (_, pg_title) in enumerate(all_figures, start=2):
            ax.text(0.12, y, f"  p.{page_n:02d}   {pg_title}",
                    ha="left", va="center", fontsize=9,
                    transform=ax.transAxes, color="#444", fontfamily="monospace")
            y -= 0.045
            if y < 0.04: break  # avoid running off page

        # Footer
        ax.text(0.5, 0.02, f"Generated by protein_distances pipeline  ·  {ts}",
                ha="center", va="bottom", fontsize=8,
                transform=ax.transAxes, color="#aaa")

        pdf.savefig(cover, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(cover)

        # ── Content pages ─────────────────────────────────────────────────────
        for fig, pg_title in all_figures:
            pdf.savefig(fig, dpi=150, bbox_inches="tight", facecolor="white")

    n_pages = len(all_figures) + 1  # +1 for cover
    print(f"  ✓  Multi-page PDF report ({n_pages} pages) → {pdf_path}")


# ══════════════════════════════════════════════════════════════════════════════
#  MATRIX I/O
# ══════════════════════════════════════════════════════════════════════════════

def save_distance_matrix(
        D: np.ndarray,
        display_labels: List[str],
        original_ids: List[str],
        base: str,
) -> None:
    n = len(display_labels)
    with open(base + ".phy", "w") as f:
        f.write(f"   {n}\n")
        for i, orig in enumerate(original_ids):
            f.write(f"{orig[:10].ljust(10)} "
                    + " ".join(f"{D[i, j]:.6f}" for j in range(n)) + "\n")
    print(f"  ✓  PHYLIP          → {base}.phy")

    with open(base + ".csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + display_labels)
        for i, lb in enumerate(display_labels):
            w.writerow([lb] + [f"{D[i, j]:.6f}" for j in range(n)])
    print(f"  ✓  CSV             → {base}.csv")

    with open(base + "_id_map.tsv", "w") as f:
        f.write("original_id\tdisplay_label\n")
        for o, d in zip(original_ids, display_labels):
            f.write(f"{o}\t{d}\n")
    print(f"  ✓  ID map TSV      → {base}_id_map.tsv")

    up = D[np.triu_indices(n, k=1)]
    with open(base + "_summary.json", "w") as f:
        json.dump({
            "n_sequences": n,
            "min_distance": float(up.min()),
            "max_distance": float(up.max()),
            "mean_distance": float(up.mean()),
            "std_distance": float(up.std()),
            "sequences": [{"original_id": o, "display_label": d}
                          for o, d in zip(original_ids, display_labels)],
        }, f, indent=2)
    print(f"  ✓  JSON summary    → {base}_summary.json")


# ══════════════════════════════════════════════════════════════════════════════
#  DEMO GENERATOR
# ══════════════════════════════════════════════════════════════════════════════

def generate_demo_alignment(
        n_seqs: int = 12,
        aln_len: int = 200,
        n_clusters: int = 3,
        seed: int = 42,
) -> Tuple[str, str]:
    np.random.seed(seed)
    AA = list("ACDEFGHIKLMNPQRSTVWY")
    per = [n_seqs // n_clusters] * n_clusters
    per[-1] += n_seqs - sum(per)
    cl_names = [chr(65 + c) for c in range(n_clusters)]
    organisms = [
        "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Gallus gallus",
        "Danio rerio", "Xenopus laevis", "Drosophila melanogaster",
        "Caenorhabditis elegans", "Saccharomyces cerevisiae", "Arabidopsis thaliana",
        "Escherichia coli", "Bacillus subtilis", "Streptomyces coelicolor",
        "Methanocaldococcus jannaschii", "Sulfolobus solfataricus",
    ]
    records = [];
    tsv_rows = [];
    org_idx = 0
    for c in range(n_clusters):
        ancestor = "".join(np.random.choice(AA, aln_len))
        for s in range(per[c]):
            mut = list(ancestor)
            for pos in np.random.choice(aln_len, int(aln_len * np.random.uniform(0.05, 0.15)), replace=False):
                mut[pos] = np.random.choice(AA)
            for pos in np.random.choice(aln_len, int(aln_len * np.random.uniform(0.02, 0.05)), replace=False):
                mut[pos] = "-"
            seq_id = f"Cluster{cl_names[c]}_{s + 1:02d}"
            org = organisms[org_idx % len(organisms)]
            records.append(SeqRecord(Seq("".join(mut)), id=seq_id, description=""))
            tsv_rows.append((seq_id, f"{org} | Cluster-{cl_names[c]}"))
            org_idx += 1

    fasta_path = "demo_protein_alignment.fasta"
    tsv_path = "demo_labels.tsv"
    with open(fasta_path, "w") as f:
        for rec in records: f.write(f">{rec.id}\n{rec.seq}\n")
    with open(tsv_path, "w") as f:
        f.write("original_id\tdisplay_label\n")
        for o, d in tsv_rows: f.write(f"{o}\t{d}\n")
    print(f"  ✓  Demo MSA    ({n_seqs} seqs, {aln_len} cols) → {fasta_path}")
    print(f"  ✓  Demo labels ({len(tsv_rows)} mappings)       → {tsv_path}")
    return fasta_path, tsv_path


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def run_pipeline(
        input_file: str,
        methods: List[str],
        label_tsv: Optional[str] = None,
        linkage_method: str = "average",
        sub_metric: str = "normalized_score",
        cosine_k: int = 2,
        fmt: Optional[str] = None,
        output_prefix: str = "protein_dist",
        output_format: str = "png",  # png | pdf | svg | all
        dpi: int = 150,
        pdf_report: bool = True,  # build combined PDF report?
        show_plots: bool = True,
        save_matrices: bool = True,
        compare_methods: bool = True,
        run_mds: bool = True,
) -> Tuple[Dict[str, np.ndarray], LabelManager]:
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print(f"\n{'═' * 65}")
    print(f"  Protein Distance Analysis Pipeline")
    print(f"{'═' * 65}")
    print(f"  Input         : {input_file}")
    print(f"  Label TSV     : {label_tsv or '(none)'}")
    print(f"  Methods       : {', '.join(methods)}")
    print(f"  Linkage       : {linkage_method}")
    print(f"  Sub-metric    : {sub_metric}")
    print(f"  k (cosine)    : {cosine_k}")
    print(f"  Output format : {output_format.upper()}")
    print(f"  DPI           : {dpi}")
    print(f"  PDF report    : {'yes' if pdf_report else 'no'}")

    # Exporter
    exporter = FigureExporter(fmt=output_format, dpi=dpi)
    print(f"  {exporter.describe()}")

    # ── 1. Labels ─────────────────────────────────────────────────────────────
    print(f"\n[1/7]  Loading labels …")
    label_mgr = LabelManager(label_tsv)
    print(label_mgr.summary())

    # ── 2. Parse ──────────────────────────────────────────────────────────────
    print(f"\n[2/7]  Parsing alignment …")
    alignment = parse_alignment(input_file, fmt=fmt)
    label_mgr.validate_against_alignment(alignment)
    n_seqs, _ = print_alignment_summary(alignment, label_mgr)
    original_ids = [rec.id for rec in alignment]
    display_labels = label_mgr.replace_all(original_ids)

    # ── 3. Compute ─────────────────────────────────────────────────────────────
    print(f"\n[3/7]  Computing distance matrices …")
    all_matrices: Dict[str, np.ndarray] = {}
    for method in methods:
        info = MATRIX_REGISTRY.get(method)
        print(f"\n  ── {method}  ({info.description if info else ''}) ──")
        all_matrices[method] = calculate_distance_matrix(
            alignment, method, sub_metric=sub_metric, cosine_k=cosine_k,
        )

    # ── 4. Save matrices ───────────────────────────────────────────────────────
    if save_matrices:
        print(f"\n[4/7]  Saving matrices …")
        for method, D in all_matrices.items():
            save_distance_matrix(D, display_labels, original_ids,
                                 f"{output_prefix}_{method}")
    else:
        print(f"\n[4/7]  Skipping matrix save.")

    # ── 5-7. Figures ───────────────────────────────────────────────────────────
    # Collect all (figure, title) pairs for the PDF report
    all_figures: List[Tuple[plt.Figure, str]] = []

    print(f"\n[5/7]  Heatmaps …")
    for method, D in all_matrices.items():
        ttl = f"Heatmap – {method}"
        base = f"{output_prefix}_{method}_heatmap"
        fig, _, _ = plot_heatmap_with_dendrogram(
            D, display_labels, original_ids,
            method=method, linkage_method=linkage_method,
            exporter=exporter, base_path=base, title=ttl,
        )
        all_figures.append((fig, ttl))

    print(f"\n[6/7]  Distribution & comparison plots …")
    fig_dist = plot_distance_distribution(
        all_matrices, display_labels,
        exporter=exporter, base_path=f"{output_prefix}_distributions",
    )
    all_figures.append((fig_dist, "Distance Distributions"))

    if compare_methods and len(methods) > 1:
        fig_cmp = plot_multi_method_comparison(
            all_matrices, display_labels,
            exporter=exporter, base_path=f"{output_prefix}_method_comparison",
        )
        all_figures.append((fig_cmp, "Multi-Method Comparison"))

        fig_cor = plot_method_correlation(
            all_matrices,
            exporter=exporter, base_path=f"{output_prefix}_method_correlation",
        )
        if fig_cor: all_figures.append((fig_cor, "Method Correlation"))

    print(f"\n[7/7]  Extra visualisations …")
    if any(m in COSINE_METHODS for m in methods):
        fig_cos = plot_cosine_encoding_comparison(
            alignment, display_labels,
            exporter=exporter, base_path=f"{output_prefix}_cosine_schemes",
        )
        all_figures.append((fig_cos, "Cosine Encoding Comparison"))

    if run_mds and len(methods) >= 1:
        fig_mds = plot_mds_sequence_map(
            all_matrices, display_labels,
            methods=methods[:min(4, len(methods))],
            exporter=exporter, base_path=f"{output_prefix}_mds_map",
        )
        all_figures.append((fig_mds, "2-D MDS Sequence Map"))

    # ── Combined PDF report ────────────────────────────────────────────────────
    if pdf_report:
        meta = dict(
            input_file=input_file,
            label_tsv=label_tsv or "(none)",
            methods=methods,
            linkage=linkage_method,
            sub_metric=sub_metric,
            n_seqs=n_seqs,
            output_prefix=output_prefix,
            timestamp=timestamp,
        )
        build_multipage_pdf_report(
            all_figures,
            pdf_path=f"{output_prefix}_REPORT.pdf",
            pipeline_meta=meta,
        )

    if show_plots:
        plt.show()

    print(f"\n{'═' * 65}")
    print(f"     Done!   Timestamp: {timestamp}")
    print(f"{'═' * 65}\n")
    return all_matrices, label_mgr


# ══════════════════════════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════════════════════════

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="protein_distances",
        description="Protein pairwise distances – PAM, BLOSUM, Classical, Cosine | PDF/PNG/SVG output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Output formats  (--output-format)
──────────────────────────────────
  png   – raster PNG  (default, 150 DPI)
  pdf   – vector PDF  (publication-quality, embedded metadata)
  svg   – vector SVG  (editable in Inkscape / Illustrator)
  all   – PNG + PDF + SVG simultaneously

PDF notes
─────────
  • Each figure is also saved as an individual PDF (when --output-format pdf/all).
  • A combined multi-page PDF report (cover page + all figures) is always
    written as  <prefix>_REPORT.pdf  unless --no-pdf-report is passed.
  • PDF metadata (title, author, date) is embedded automatically.

TSV label file  (--labels)
───────────────────────────
  # comment
  original_id<TAB>display_label
  SeqID_001<TAB>Homo sapiens hemoglobin alpha

Method families
───────────────
  PAM       : {', '.join(PAM_METHODS)}
  BLOSUM    : {', '.join(BLOSUM_METHODS)}
  Classical : {', '.join(CLASSICAL_METHODS)}
  Cosine    : {', '.join(COSINE_METHODS)}

Examples
────────
  # Demo (auto-generates MSA + TSV + PDF report)
  ppwd_vis.py --demo --output-format pdf

  # PDF output with custom labels
  ppwd_vis.py -i aln.fasta --labels labels.tsv \\
      -m BLOSUM62 kimura --output-format pdf

  # All formats simultaneously
  ppwd_vis.py -i aln.fasta -m PAM120 cosine_physico \\
      --output-format all --dpi 300

  # Single SVG per figure, no combined report
  ppwd_vis.py -i aln.fasta -m BLOSUM62 \\
      --output-format svg --no-pdf-report
        """
    )
    p.add_argument("-i", "--input", type=str, help="MSA file path")
    p.add_argument("--labels", type=str, default=None, metavar="TSV",
                   help="Optional 2-column TSV: original_id → display_label")
    p.add_argument("-m", "--methods", type=str, nargs="+",
                   default=["BLOSUM62"], metavar="METHOD")
    p.add_argument("-l", "--linkage", type=str, default="average",
                   choices=["single", "complete", "average", "ward", "weighted"])
    p.add_argument("--sub-metric", type=str, default="normalized_score",
                   choices=["normalized_score", "log_odds", "percent_similarity"])
    p.add_argument("--cosine-k", type=int, default=2)
    p.add_argument("-f", "--format", type=str, default=None,
                   help="Alignment format: fasta|clustal|phylip|stockholm|nexus")
    p.add_argument("-o", "--output", type=str, default="protein_distances",
                   help="Output file prefix (default: protein_distances)")
    p.add_argument("--output-format", type=str, default="png",
                   choices=SUPPORTED_FORMATS,
                   help="Figure output format: png | pdf | svg | all  (default: png)")
    p.add_argument("--dpi", type=int, default=150,
                   help="Resolution for PNG output (default: 150)")
    p.add_argument("--no-pdf-report", action="store_true",
                   help="Skip the combined multi-page PDF report")
    p.add_argument("--demo", action="store_true")
    p.add_argument("--all-pam", action="store_true")
    p.add_argument("--all-blosum", action="store_true")
    p.add_argument("--all-cosine", action="store_true")
    p.add_argument("--no-show", action="store_true")
    p.add_argument("--no-save-matrix", action="store_true")
    p.add_argument("--no-compare", action="store_true")
    p.add_argument("--no-mds", action="store_true")
    p.add_argument("--list-methods", action="store_true")
    p.add_argument("--validate-tsv", type=str, default=None, metavar="TSV")
    return p


def list_methods() -> None:
    print(f"\n{'─' * 72}")
    print(f"  {'Method':<20}  {'Family':<10}  Description")
    print(f"{'─' * 72}")
    for name, info in MATRIX_REGISTRY.items():
        print(f"  {name:<20}  {info.family:<10}  {info.description}")
    print(f"{'─' * 72}\n")


# ══════════════════════════════════════════════════════════════════════════════
#  ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    args = build_parser().parse_args()

    if args.list_methods:
        list_methods();
        sys.exit(0)

    if args.validate_tsv:
        if not args.input:
            print("  ✗  --validate-tsv requires -i/--input");
            sys.exit(1)
        aln = parse_alignment(args.input, args.format)
        mgr = LabelManager(args.validate_tsv)
        mgr.validate_against_alignment(aln)
        print(mgr.summary());
        sys.exit(0)

    methods = list(args.methods)
    if args.all_pam:
        methods = PAM_METHODS + [m for m in methods if m not in PAM_METHODS]
    if args.all_blosum:
        methods = BLOSUM_METHODS + [m for m in methods if m not in BLOSUM_METHODS]
    if args.all_cosine:
        methods = COSINE_METHODS + [m for m in methods if m not in COSINE_METHODS]

    for m in methods:
        if m not in MATRIX_REGISTRY:
            print(f"  ✗  Unknown method '{m}'. Use --list-methods.");
            sys.exit(1)

    label_tsv = args.labels
    if args.demo or not args.input:
        print("\n  🧬  DEMO mode …")
        input_file, label_tsv = generate_demo_alignment(n_seqs=12, n_clusters=3)
        if set(methods) == {"BLOSUM62"}:
            methods = ["PAM120", "BLOSUM62", "cosine_physico", "cosine_kmer", "kimura"]
    else:
        if not os.path.exists(args.input):
            print(f"  ✗  File not found: '{args.input}'");
            sys.exit(1)
        input_file = args.input

    run_pipeline(
        input_file=input_file,
        methods=methods,
        label_tsv=label_tsv,
        linkage_method=args.linkage,
        sub_metric=args.sub_metric,
        cosine_k=args.cosine_k,
        fmt=args.format,
        output_prefix=args.output,
        output_format=args.output_format,
        dpi=args.dpi,
        pdf_report=not args.no_pdf_report,
        show_plots=not args.no_show,
        save_matrices=not args.no_save_matrix,
        compare_methods=not args.no_compare,
        run_mds=not args.no_mds,
    )
