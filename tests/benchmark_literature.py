#!/usr/bin/env python3
"""
Literature-based benchmark harness for pdifFinder.

Compares pdifFinder output against manually curated pdif site coordinates
from published literature to measure recall, precision, and positional accuracy.

Skeleton — actual comparison/metric/mock logic lives in Tasks 8-12.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import statistics
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from unittest.mock import patch

from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

POSITION_TOLERANCE: int = 5


@dataclass
class ReferenceEntry:
    """One curated pdif site from the literature reference CSV."""
    plasmid_accession: str
    pdif_start: int
    pdif_end: int
    pdif_sequence: str
    strand: str
    orientation: str
    source_pmid: str
    source_description: str
    curator_notes: str
    curation_date: str


@dataclass
class PdifOutput:
    """Parsed pdifFinder output for a single plasmid run (stub)."""
    plasmid_accession: str
    pdif_sites: List[dict] = field(default_factory=list)
    amr_genes: List[dict] = field(default_factory=list)
    pdif_modules: List[dict] = field(default_factory=list)
    return_code: int = 0
    stdout: str = ""
    stderr: str = ""


@dataclass
class ComparisonResult:
    """Result of comparing pdifFinder output against literature references for one plasmid."""
    plasmid_accession: str = ""
    matched_pairs: List[dict] = field(default_factory=list)
    missed_references: List[ReferenceEntry] = field(default_factory=list)
    unvalidated_detections: List[dict] = field(default_factory=list)
    failure_categories: Dict[str, int] = field(default_factory=dict)


@dataclass
class FailureCategory:
    """Categorisation of why a reference site was not matched (stub)."""
    category: str
    detail: str = ""


@dataclass
class Metrics:
    """Aggregate benchmark metrics across all plasmids."""
    detection_rate: float = 0.0
    precision: float = 0.0
    position_median_error: float = -1.0
    position_mean_error: float = -1.0
    position_max_error: float = -1.0
    fpr_per_kb: float = -1.0
    per_plasmid_breakdown: List[dict] = field(default_factory=list)
    failure_category_summary: Dict[str, int] = field(default_factory=dict)
    total_reference_sites: int = 0
    total_matched_sites: int = 0
    total_missed_sites: int = 0
    total_unvalidated: int = 0

    @property
    def recall(self) -> float:
        """Alias for detection_rate (backward compat)."""
        return self.detection_rate

    @property
    def f1(self) -> float:
        """F1 score computed from detection_rate (recall) and precision."""
        if self.detection_rate + self.precision == 0:
            return 0.0
        return 2 * (self.precision * self.detection_rate) / (self.precision + self.detection_rate)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class BenchmarkResult:
    """Top-level result of a full benchmark run."""
    mode: str
    timestamp: str
    plasmid_filter: Optional[str]
    comparisons: List[ComparisonResult] = field(default_factory=list)
    failures: List[str] = field(default_factory=list)
    metrics: Metrics = field(default_factory=Metrics)
    output_dir: str = ""

    def to_dict(self) -> dict:
        return {
            "mode": self.mode,
            "timestamp": self.timestamp,
            "plasmid_filter": self.plasmid_filter,
            "comparisons": [asdict(c) for c in self.comparisons],
            "failures": self.failures,
            "metrics": self.metrics.to_dict(),
            "output_dir": self.output_dir,
        }


class LiteratureBenchmark:
    """Benchmark harness that compares pdifFinder output against published literature.

    Supports two modes:
      - **pdif_only**: Tests pdif site detection in isolation by mocking AMR gene
        detection (``findResistanceGene``). Does not require blastn. Useful for
        quick validation of the core sliding-window algorithm.
      - **full_pipeline**: End-to-end test with real blastn-based AMR gene
        detection. Requires blastn in PATH. Prefers GenBank input when available.

    Usage::

        python -m tests.benchmark_literature --mode pdif_only
        python -m tests.benchmark_literature --mode full_pipeline

    Results include detection rate (recall), precision, positional error, and
    per-plasmid breakdown. See ``compute_metrics()`` for threshold targets.
    """

    def __init__(
        self,
        mode: str,
        output_dir: str,
        plasmid_filter: Optional[str] = None,
    ) -> None:
        self.mode = mode
        self.output_dir = output_dir
        self.plasmid_filter = plasmid_filter

        self._base_dir = Path(__file__).resolve().parent
        self._benchmark_data = self._base_dir / "benchmark_data"
        self._reference_csv = self._benchmark_data / "literature_reference.csv"
        self._plasmids_dir = self._benchmark_data / "plasmids"

        self.log = logging.getLogger("LiteratureBenchmark")
        self._seeds_cache: Optional[List[Tuple[str, str]]] = None

    @staticmethod
    def _non_comment_lines(filepath: Path):
        """Yield lines from *filepath* that are not blank and do not start with #."""
        with open(filepath, "r") as fh:
            for raw in fh:
                stripped = raw.strip()
                if stripped and not stripped.startswith("#"):
                    yield stripped

    def load_reference_data(self) -> List[ReferenceEntry]:
        """Read literature_reference.csv and return list of ReferenceEntry."""
        entries: List[ReferenceEntry] = []
        if not self._reference_csv.exists():
            self.log.warning("Reference CSV not found: %s", self._reference_csv)
            return entries

        reader = csv.DictReader(self._non_comment_lines(self._reference_csv))
        for row in reader:
            try:
                entry = ReferenceEntry(
                    plasmid_accession=row.get("plasmid_accession", "").strip(),
                    pdif_start=int(row.get("pdif_start", 0)),
                    pdif_end=int(row.get("pdif_end", 0)),
                    pdif_sequence=row.get("pdif_sequence", "").strip(),
                    strand=row.get("strand", "+").strip(),
                    orientation=row.get("orientation", "").strip(),
                    source_pmid=row.get("source_pmid", "").strip(),
                    source_description=row.get("source_description", "").strip(),
                    curator_notes=row.get("curator_notes", "").strip(),
                    curation_date=row.get("curation_date", "").strip(),
                )
                entries.append(entry)
            except (ValueError, KeyError) as exc:
                self.log.warning("Skipping malformed row: %s | %s", exc, row)

        self.log.info("Loaded %d reference entries from %s", len(entries), self._reference_csv.name)
        return entries

    def load_plasmid_sequences(self) -> Dict[str, str]:
        """Read FASTA files from tests/benchmark_data/plasmids/ → {accession: sequence}."""
        try:
            from Bio import SeqIO
        except ImportError:
            self.log.error("BioPython is required (pip install biopython)")
            return {}

        sequences: Dict[str, str] = {}
        if not self._plasmids_dir.exists():
            self.log.warning("Plasmids directory not found: %s", self._plasmids_dir)
            return sequences

        for fasta_path in sorted(self._plasmids_dir.glob("*.fasta")):
            if fasta_path.stem == "all_plasmids":
                continue
            try:
                record = SeqIO.read(str(fasta_path), "fasta")
                accession = record.id.split()[0]
                sequences[accession] = str(record.seq)
            except Exception as exc:
                self.log.warning("Failed to read FASTA %s: %s", fasta_path.name, exc)

        self.log.info("Loaded %d plasmid sequences from %s", len(sequences), self._plasmids_dir)
        return sequences

    def run_pdif_finder(self, fasta_path: str, output_dir: str) -> PdifOutput:
        """Run pdifFinder on a single FASTA file, dispatching by mode."""
        self.log.info("run_pdif_finder(%s, %s) mode=%s", fasta_path, output_dir, self.mode)
        acc = Path(fasta_path).stem.replace("_", ".")

        if self.mode == "pdif_only":
            return self._run_pdif_only(fasta_path, output_dir)
        elif self.mode == "full_pipeline":
            return self._run_full_pipeline(fasta_path, output_dir)
        else:
            self.log.warning("Unknown mode '%s' — returning empty PdifOutput", self.mode)
            return PdifOutput(plasmid_accession=acc)

    def _run_pdif_only(self, fasta_path: str, output_dir: str) -> PdifOutput:
        """Run pdifFinder's pdif detection with mocked AMR gene detection.

        Uses ``unittest.mock.patch`` on ``PdifFinder.pdifFinder.findResistanceGene``
        to bypass the AMR guard at ``singleThread()`` line 797 so that pdif site
        detection runs without requiring real resistance gene BLAST hits.

        Returns a ``PdifOutput`` with parsed results from ``parse_pdif_output()``.
        """
        import PdifFinder.pdifFinder as pf
        from tests.helpers import parse_pdif_output

        scripts_dir = os.path.dirname(os.path.abspath(pf.__file__))
        acc = Path(fasta_path).stem.replace("_", ".")

        os.makedirs(output_dir, exist_ok=True)
        tmp_dir = os.path.join(output_dir, "tmp")
        os.makedirs(tmp_dir, exist_ok=True)
        for sub in ("pairSearch", "seedSearch"):
            os.makedirs(os.path.join(tmp_dir, sub), exist_ok=True)

        # make_sort() will be called on AMRgene1.txt at line 798; since our mock
        # does not write that file we create an empty placeholder to avoid a
        # FileNotFoundError inside make_sort.
        amr_path = os.path.join(output_dir, "AMRgene1.txt")
        if not os.path.exists(amr_path):
            Path(amr_path).touch()
        # changepdifname reads AMRgene1.txt and writes AMRgene.txt only when
        # there are lines to process. With an empty AMRgene1.txt the output
        # file is never created, which crashes Getpdifmoduleseq. Touch it here
        # so the empty-file guard in Getpdifmoduleseq handles it gracefully.
        amr_gene_path = os.path.join(output_dir, "AMRgene.txt")
        if not os.path.exists(amr_gene_path):
            Path(amr_gene_path).touch()

        mock_find_rg = patch(
            "PdifFinder.pdifFinder.findResistanceGene",
            return_value=["1-1000"],
        )

        try:
            mock_find_rg.start()
            pf.singleThread(
                inFile=fasta_path,
                outdir=output_dir,
                blastnPath="/usr/bin/blastn",
                scriptsDir=scripts_dir,
                circle="true",
            )
            parsed = parse_pdif_output(output_dir)
        finally:
            mock_find_rg.stop()

        return PdifOutput(
            plasmid_accession=acc,
            pdif_sites=parsed.get("pdif_sites", []),
            amr_genes=parsed.get("amr_genes", []),
            pdif_modules=parsed.get("pdif_modules", []),
        )

    @staticmethod
    def _circular_distance(a: int, b: int, length: int) -> int:
        """Shortest distance between two positions on a circular molecule."""
        d = abs(a - b)
        return min(d, length - d)

    def _compare_plasmid(
        self, references: List[ReferenceEntry], detected: PdifOutput,
        plasmid_length: int = 0, is_circular: bool = False,
    ) -> ComparisonResult:
        matched_pairs: List[dict] = []
        missed_references: List[ReferenceEntry] = []
        used_detected: set = set()

        for ref in references:
            ref_start = ref.pdif_start
            ref_end = ref.pdif_end
            ref_seq = ref.pdif_sequence.upper()
            found = False

            for det_idx, det_site in enumerate(detected.pdif_sites):
                if det_idx in used_detected:
                    continue
                det_start = int(det_site.get("start", 0))
                det_end = int(det_site.get("end", 0))
                det_seq = det_site.get("pdif_site", "").upper()

                if not det_seq or det_start <= 0 or det_end <= 0:
                    continue

                if is_circular and plasmid_length > 0:
                    start_dist = self._circular_distance(ref_start, det_start, plasmid_length)
                    end_dist = self._circular_distance(ref_end, det_end, plasmid_length)
                    pos_match = (start_dist <= POSITION_TOLERANCE and end_dist <= POSITION_TOLERANCE)
                else:
                    pos_match = (abs(det_start - ref_start) <= POSITION_TOLERANCE
                                 and abs(det_end - ref_end) <= POSITION_TOLERANCE)

                if not pos_match:
                    continue

                rev_comp = str(Seq(det_seq).reverse_complement())
                if det_seq == ref_seq or rev_comp == ref_seq:
                    matched_pairs.append(dict(
                        ref_start=ref_start, ref_end=ref_end,
                        det_index=det_idx, det_start=det_start, det_end=det_end,
                    ))
                    used_detected.add(det_idx)
                    found = True
                    break

            if not found:
                missed_references.append(ref)

        unvalidated = [
            dict(det_index=i)
            for i in range(len(detected.pdif_sites))
            if i not in used_detected
        ]

        return ComparisonResult(
            plasmid_accession=references[0].plasmid_accession if references else "unknown",
            matched_pairs=matched_pairs,
            missed_references=missed_references,
            unvalidated_detections=unvalidated,
            failure_categories={},
        )

    def compare_results(
        self, reference: ReferenceEntry, detected: PdifOutput,
        plasmid_length: int = 0, is_circular: bool = False,
    ) -> ComparisonResult:
        """Compare one reference pdif site against detected pdif sites.

        Matches detected sites within ``POSITION_TOLERANCE`` bp of the
        reported reference coordinates, supporting circular coordinate
        wrapping and strand-agnostic (reverse-complement) sequence matching.
        """
        ref_start = reference.pdif_start
        ref_end = reference.pdif_end
        ref_seq = reference.pdif_sequence.upper()

        matched_pairs: List[dict] = []
        used_detected: set = set()

        for det_idx, det_site in enumerate(detected.pdif_sites):
            det_start = int(det_site.get("start", 0))
            det_end = int(det_site.get("end", 0))
            det_seq = det_site.get("pdif_site", "").upper()

            if not det_seq or det_start <= 0 or det_end <= 0:
                continue

            if is_circular and plasmid_length > 0:
                start_dist = self._circular_distance(
                    ref_start, det_start, plasmid_length)
                end_dist = self._circular_distance(
                    ref_end, det_end, plasmid_length)
                pos_match = (
                    start_dist <= POSITION_TOLERANCE
                    and end_dist <= POSITION_TOLERANCE
                )
            else:
                pos_match = (
                    abs(det_start - ref_start) <= POSITION_TOLERANCE
                    and abs(det_end - ref_end) <= POSITION_TOLERANCE
                )

            if not pos_match:
                continue

            rev_comp = str(Seq(det_seq).reverse_complement())
            if det_seq == ref_seq or rev_comp == ref_seq:
                matched_pairs.append(dict(
                    ref_index=0,
                    det_index=det_idx,
                    ref_start=ref_start,
                    ref_end=ref_end,
                    det_start=det_start,
                    det_end=det_end,
                ))
                used_detected.add(det_idx)

        unvalidated = [
            dict(det_index=i)
            for i in range(len(detected.pdif_sites))
            if i not in used_detected
        ]

        if matched_pairs:
            return ComparisonResult(
                plasmid_accession=reference.plasmid_accession,
                matched_pairs=matched_pairs,
                missed_references=[],
                unvalidated_detections=unvalidated,
                failure_categories={},
            )

        return ComparisonResult(
            plasmid_accession=reference.plasmid_accession,
            matched_pairs=[],
            missed_references=[reference],
            unvalidated_detections=unvalidated,
            failure_categories={},
        )

    # ------------------------------------------------------------------
    # seed database helpers (Task 11)
    # ------------------------------------------------------------------

    def _load_seeds(self) -> List[Tuple[str, str]]:
        """Parse redundant.seed.fa → list of (xerC, xerD) 11-mer pairs.

        Each seed entry is 28 bp: xerC[0:11]  spacer[11:17]  xerD[17:28].
        Cached on first call so every failure-categorisation reuses the
        same parsed list.
        """
        if self._seeds_cache is not None:
            return self._seeds_cache

        seeds: List[Tuple[str, str]] = []
        seed_path = (
            self._base_dir.parent / "PdifFinder" / "data" / "redundant.seed.fa"
        )
        if not seed_path.exists():
            self.log.warning("Seed database not found: %s", seed_path)
            self._seeds_cache = seeds
            return seeds

        current_seq = ""
        for raw in seed_path.read_text().splitlines():
            line = raw.strip()
            if line.startswith(">"):
                if len(current_seq) == 28:
                    seeds.append((current_seq[0:11], current_seq[17:28]))
                current_seq = ""
            elif line:
                current_seq += line

        if len(current_seq) == 28:
            seeds.append((current_seq[0:11], current_seq[17:28]))

        self._seeds_cache = seeds
        self.log.debug("Loaded %d seeds from redundant.seed.fa", len(seeds))
        return seeds

    @staticmethod
    def _hamming_distance(seq1: str, seq2: str) -> int:
        return sum(1 for a, b in zip(seq1, seq2) if a != b)

    # ------------------------------------------------------------------
    # failure categorisation (Task 11)
    # ------------------------------------------------------------------

    def categorize_failure(
        self,
        reference_entry: ReferenceEntry,
        detected_list: List[dict],
        amr_found: bool = True,
    ) -> str:
        """Categorise why a reference pdif site was missed.

        Priority order (first matching wins):
          1. NOT_SOUGHT    – AMR not found so pdifFinder never searched
          2. NO_SEED_MATCH – neither xerC nor xerD matches any seed
          3. NO_PAIR       – halves match seeds but not as a single-entry pair
          4. WRONG_POSITION – seeds & pair OK but pdifFinder placed site
                               >5 bp from the literature coordinate (default)
        """
        # ── 1. NOT_SOUGHT ──────────────────────────────────────────
        if not amr_found:
            return "NOT_SOUGHT"

        pdif_seq = reference_entry.pdif_sequence
        if len(pdif_seq) < 28:
            return "WRONG_POSITION"

        xerC_ref = pdif_seq[0:11].upper()
        xerD_ref = pdif_seq[17:28].upper()

        seeds = self._load_seeds()
        if not seeds:
            return "WRONG_POSITION"

        # ── 2. NO_SEED_MATCH ───────────────────────────────────────
        # Check whether ANY seed's xerC matches *and* ANY seed's xerD
        # matches (they can be from different entries at this stage).
        # pdifFinder uses two threshold passes: (3,2) then (2,3).
        xerC_matched = False
        xerD_matched = False
        threshold_pairs = [(3, 2), (2, 3)]

        for seed_xc, seed_xd in seeds:
            for th_c, th_d in threshold_pairs:
                if self._hamming_distance(seed_xc, xerC_ref) <= th_c:
                    xerC_matched = True
                if self._hamming_distance(seed_xd, xerD_ref) <= th_d:
                    xerD_matched = True
            if xerC_matched and xerD_matched:
                break

        if not xerC_matched or not xerD_matched:
            return "NO_SEED_MATCH"

        # ── 3. NO_PAIR ─────────────────────────────────────────────
        # Simplified pairing check: both halves must match the *same*
        # seed entry.  If xerC and xerD each match a seed individually
        # but no single seed matches both, the pairing step would fail.
        same_seed_matched = False
        for seed_xc, seed_xd in seeds:
            for th_c, th_d in threshold_pairs:
                if (self._hamming_distance(seed_xc, xerC_ref) <= th_c
                        and self._hamming_distance(seed_xd, xerD_ref) <= th_d):
                    same_seed_matched = True
                    break
            if same_seed_matched:
                break

        if not same_seed_matched:
            return "NO_PAIR"

        # ── 4. WRONG_POSITION (default) ────────────────────────────
        return "WRONG_POSITION"

    def compute_metrics(
        self,
        all_comparisons: List[ComparisonResult],
        total_plasmid_kb: float = 0,
        has_negative_controls: bool = False,
    ) -> Metrics:
        """Compute benchmark metrics from comparison results.

        Produces four primary metrics:
          - **Detection Rate (Recall)**: fraction of reference sites matched
          - **Precision**: fraction of detected sites that match a reference
          - **Position Error**: median/mean/max bp deviation on matched sites
          - **FPR per kb**: false positives per kilobase (only with negative controls)

        Also computes per-plasmid breakdown and failure-category summary.

        Target thresholds (from decisions.md):
          Detection Rate >= 0.80, Precision >= 0.70,
          Position Median Error <= 5 bp, FPR <= 0.01/kb
        """

        # --- 1. Detection Rate (Recall) ---
        total_reference_sites = sum(
            len(c.missed_references) + len(c.matched_pairs) for c in all_comparisons
        )
        total_matched_sites = sum(len(c.matched_pairs) for c in all_comparisons)
        total_missed_sites = sum(len(c.missed_references) for c in all_comparisons)
        total_unvalidated = sum(len(c.unvalidated_detections) for c in all_comparisons)

        detection_rate = (
            total_matched_sites / total_reference_sites
            if total_reference_sites > 0
            else 0.0
        )

        # --- 2. Precision (excluding unvalidated detections) ---
        total_detected = total_matched_sites + total_unvalidated
        precision = (
            total_matched_sites / total_detected if total_detected > 0 else 0.0
        )

        # --- 3. Position Accuracy (only on matched pairs) ---
        position_errors: List[float] = []
        for c in all_comparisons:
            for pair in c.matched_pairs:
                ref_start = int(pair.get("ref_start", 0))
                ref_end = int(pair.get("ref_end", 0))
                det_start = int(pair.get("det_start", 0))
                det_end = int(pair.get("det_end", 0))
                error = max(abs(det_start - ref_start), abs(det_end - ref_end))
                position_errors.append(float(error))

        if position_errors:
            position_median_error = statistics.median(position_errors)
            position_mean_error = statistics.mean(position_errors)
            position_max_error = max(position_errors)
        else:
            position_median_error = -1.0
            position_mean_error = -1.0
            position_max_error = -1.0

        # --- 4. False Positive Rate (only for negative controls) ---
        if has_negative_controls and total_plasmid_kb > 0:
            false_positives = sum(
                len(c.unvalidated_detections) for c in all_comparisons
                if len(c.missed_references) == 0 and len(c.matched_pairs) == 0
            )
            fpr_per_kb = false_positives / total_plasmid_kb
        else:
            fpr_per_kb = -1.0

        # --- 5. Per-plasmid breakdown ---
        per_plasmid: Dict[str, dict] = {}
        for c in all_comparisons:
            acc = c.plasmid_accession
            if acc not in per_plasmid:
                per_plasmid[acc] = {
                    "plasmid_accession": acc,
                    "reference_count": 0,
                    "matched_count": 0,
                    "missed_count": 0,
                    "unvalidated_count": 0,
                    "detection_rate": 0.0,
                    "failure_categories": {},
                }
            entry = per_plasmid[acc]
            entry["reference_count"] += len(c.missed_references) + len(c.matched_pairs)
            entry["matched_count"] += len(c.matched_pairs)
            entry["missed_count"] += len(c.missed_references)
            entry["unvalidated_count"] += len(c.unvalidated_detections)

        # Compute per-plasmid detection rates
        for entry in per_plasmid.values():
            ref_count = entry["reference_count"]
            entry["detection_rate"] = (
                entry["matched_count"] / ref_count if ref_count > 0 else 0.0
            )

        # Merge failure categories per plasmid
        for c in all_comparisons:
            acc = c.plasmid_accession
            for cat, count in c.failure_categories.items():
                per_plasmid[acc]["failure_categories"][cat] = (
                    per_plasmid[acc]["failure_categories"].get(cat, 0) + count
                )

        per_plasmid_breakdown = list(per_plasmid.values())

        # --- 6. Failure category summary ---
        failure_summary: Dict[str, int] = {}
        for c in all_comparisons:
            for cat, count in c.failure_categories.items():
                failure_summary[cat] = failure_summary.get(cat, 0) + count

        return Metrics(
            detection_rate=detection_rate,
            precision=precision,
            position_median_error=position_median_error,
            position_mean_error=position_mean_error,
            position_max_error=position_max_error,
            fpr_per_kb=fpr_per_kb,
            per_plasmid_breakdown=per_plasmid_breakdown,
            failure_category_summary=failure_summary,
            total_reference_sites=total_reference_sites,
            total_matched_sites=total_matched_sites,
            total_missed_sites=total_missed_sites,
            total_unvalidated=total_unvalidated,
        )

    def _run_full_pipeline(self, fasta_path: str, output_dir: str) -> PdifOutput:
        """Run pdifFinder end-to-end with real blastn-based AMR gene detection.

        Unlike ``_run_pdif_only``, this mode does NOT mock any pdifFinder internals.
        It calls ``singleThread()`` with actual blastn, requiring ``blastn`` to be
        installed and available in the system PATH. Checks blastn availability at
        runtime and returns an empty ``PdifOutput`` with ``stderr="blastn_unavailable"``
        if it is missing.

        Prefers GenBank input when a ``.gb`` or ``.gbk`` file exists for the plasmid
        in ``tests/benchmark_data/plasmids/``, falling back to FASTA otherwise.
        """
        import PdifFinder.pdifFinder as pf

        acc = Path(fasta_path).stem.replace("_", ".")

        blastnPath = pf.check_dependencies()
        if not blastnPath:
            self.log.warning(
                "blastn not available — skipping full_pipeline for %s", acc,
            )
            return PdifOutput(
                plasmid_accession=acc,
                stderr="blastn_unavailable",
            )

        scriptsDir = os.path.dirname(pf.__file__)

        plasmid_outdir = os.path.join(output_dir, acc.replace(".", "_"))
        pf.makeOutdir(plasmid_outdir)

        circle = "true"
        circularSeq = "true"

        gb_path = self._plasmids_dir / f"{acc}.gb"
        gbk_path = self._plasmids_dir / f"{acc}.gbk"

        if gb_path.exists():
            in_file, is_circular = pf.getSeqFromGenbankFile(str(gb_path), plasmid_outdir)
            circularSeq = "true" if is_circular else "false"
        elif gbk_path.exists():
            in_file, is_circular = pf.getSeqFromGenbankFile(str(gbk_path), plasmid_outdir)
            circularSeq = "true" if is_circular else "false"
        else:
            in_file = pf.getSeqFromFastaFile(fasta_path, plasmid_outdir)

        try:
            pf.singleThread(
                in_file, plasmid_outdir, blastnPath, scriptsDir, circle, circularSeq,
            )
        except Exception as exc:
            self.log.error(
                "pdifFinder singleThread crashed for %s: %s", acc, exc,
            )
            return PdifOutput(
                plasmid_accession=acc,
                return_code=1,
                stderr=str(exc),
            )

        from tests.helpers import parse_pdif_output

        parsed = parse_pdif_output(plasmid_outdir)

        return PdifOutput(
            plasmid_accession=acc,
            pdif_sites=parsed.get("pdif_sites", []),
            amr_genes=parsed.get("amr_genes", []),
            pdif_modules=parsed.get("pdif_modules", []),
        )

    def _filter_entries(self, entries: List[ReferenceEntry]) -> List[ReferenceEntry]:
        if self.plasmid_filter:
            filtered = [e for e in entries if e.plasmid_accession == self.plasmid_filter]
            self.log.info(
                "Plasmid filter '%s': %d / %d entries kept",
                self.plasmid_filter, len(filtered), len(entries),
            )
            return filtered
        return entries

    def run(self) -> BenchmarkResult:
        self.log.info("%s", "=" * 60)
        self.log.info("Literature Benchmark — mode=%s", self.mode)
        self.log.info("%s", "=" * 60)

        references = self.load_reference_data()
        references = self._filter_entries(references)
        if not references:
            self.log.warning("No reference entries to benchmark (empty or all filtered)")

        plasmids = self.load_plasmid_sequences()
        self.log.info("Plasmids available: %s", list(plasmids.keys()))

        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = Path(self.output_dir) / stamp
        run_dir.mkdir(parents=True, exist_ok=True)

        all_comparisons: List[ComparisonResult] = []
        all_failures: List[str] = []

        target_accessions = {r.plasmid_accession for r in references}
        for acc in sorted(target_accessions):
            if acc not in plasmids:
                self.log.warning("No plasmid sequence for %s — skipping", acc)
                continue

            tmp_fasta = run_dir / f"{acc.replace('.', '_')}.fasta"
            tmp_fasta.write_text(f">{acc}\n{plasmids[acc]}\n")
            # pdifFinder internally hardcodes 'inputFile.fasta' for angularPlasmid;
            # write a copy so angularPlasmid can find the sequence.
            input_fasta = run_dir / "inputFile.fasta"
            if not input_fasta.exists():
                input_fasta.symlink_to(tmp_fasta.name)

            pdif_out = self.run_pdif_finder(str(tmp_fasta), str(run_dir))

            # Collect ALL reference entries for this plasmid
            plasmid_refs = [r for r in references if r.plasmid_accession == acc]
            if not plasmid_refs:
                continue

            # Compare ALL reference sites against ALL detected sites at once
            # (not one-at-a-time, which inflated unvalidated counts)
            comp = self._compare_plasmid(plasmid_refs, pdif_out)
            all_comparisons.append(comp)

            # Categorize each missed reference
            for ref in comp.missed_references:
                cat = self.categorize_failure(ref, pdif_out.pdif_sites)
                all_failures.append(cat)

        metrics = self.compute_metrics(all_comparisons)

        result = BenchmarkResult(
            mode=self.mode,
            timestamp=datetime.now().isoformat(),
            plasmid_filter=self.plasmid_filter,
            comparisons=all_comparisons,
            failures=all_failures,
            metrics=metrics,
            output_dir=str(run_dir),
        )

        result_path = run_dir / "benchmark_result.json"
        result_path.write_text(json.dumps(result.to_dict(), indent=2, default=str))
        self.log.info("Benchmark result written to %s", result_path)

        self.log.info("Metrics: precision=%.3f recall=%.3f f1=%.3f",
                       metrics.precision, metrics.recall, metrics.f1)
        self.log.info("Result output_dir: %s", run_dir)

        return result


def _setup_logging(output_dir: str, verbose: bool = False) -> logging.Logger:
    """Configure dual logging: INFO→console, DEBUG→file."""
    logger = logging.getLogger("LiteratureBenchmark")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, "benchmark.log")
    fh = logging.FileHandler(log_file, mode="a")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)-7s] %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(logging.Formatter("%(levelname)-7s | %(message)s"))
    logger.addHandler(ch)

    logger.debug("Logging initialised — file: %s", log_file)
    return logger


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Literature-based benchmark for pdifFinder",
    )
    parser.add_argument(
        "--mode",
        choices=["pdif_only", "full_pipeline"],
        required=True,
        help="Benchmark mode: pdif_only (just pdif site detection) or full_pipeline",
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parent / "benchmark_data" / "results"),
        help="Output directory for benchmark results (default: tests/benchmark_data/results)",
    )
    parser.add_argument(
        "--plasmid",
        default=None,
        help="Filter benchmark to a single plasmid accession (e.g. KY984047.1)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable DEBUG-level console output",
    )
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    _setup_logging(args.output_dir, verbose=args.verbose)

    bench = LiteratureBenchmark(
        mode=args.mode,
        output_dir=args.output_dir,
        plasmid_filter=args.plasmid,
    )
    bench.run()


if __name__ == "__main__":
    main()
