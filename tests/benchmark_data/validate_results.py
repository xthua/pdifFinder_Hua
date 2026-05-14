#!/usr/bin/env python3
"""
validate_results.py — Validate benchmark results JSON against schema.

Validates two known output formats:
  1. BenchmarkResult  — direct output of BenchmarkResult.to_dict() (e.g. pdif_only_results.json)
  2. ComparisonSummary — side-by-side comparison output (e.g. full_pipeline_results.json)

Exits with code 0 if valid, 1 with error messages if invalid.

Usage:
    python tests/benchmark_data/validate_results.py <results.json>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List


# ---------------------------------------------------------------------------
# Schema (loaded from results_schema.json)
# ---------------------------------------------------------------------------

SCHEMA_PATH = Path(__file__).resolve().parent / "results" / "results_schema.json"


def _load_schema() -> dict:
    """Load the JSON Schema from the sibling schema file."""
    if not SCHEMA_PATH.exists():
        print(f"ERROR: Schema file not found: {SCHEMA_PATH}", file=sys.stderr)
        sys.exit(1)
    with SCHEMA_PATH.open(encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Type-checking helpers
# ---------------------------------------------------------------------------

Errors = List[str]


def _check_type(value: Any, expected: type | tuple, path: str, errors: Errors) -> None:
    """Append error to *errors* if *value* is not of *expected* type."""
    if not isinstance(value, expected):
        types_str = getattr(expected, "__name__", str(expected))
        errors.append(
            f"{path}: expected {types_str}, got {type(value).__name__}"
        )


def _check_optional_type(value: Any, expected: type | tuple, path: str, errors: Errors) -> None:
    """Type-check a value that may be None."""
    if value is not None:
        _check_type(value, expected, path, errors)


def _check_required_keys(obj: dict, required: List[str], path: str, errors: Errors) -> None:
    """Check that all *required* keys exist in *obj*."""
    for key in required:
        if key not in obj:
            errors.append(f"{path}: missing required key '{key}'")


def _check_extra_keys(obj: dict, allowed: set, path: str, errors: Errors) -> None:
    """Report any keys in *obj* not in *allowed*."""
    extra = set(obj.keys()) - allowed
    if extra:
        errors.append(f"{path}: unexpected extra key(s): {sorted(extra)}")


# ---------------------------------------------------------------------------
# Validator for BenchmarkResult format  (matches results_schema.json)
# ---------------------------------------------------------------------------

VALID_MODES = {"pdif_only", "full_pipeline"}

COMPARISON_REQUIRED = [
    "plasmid_accession", "matched_pairs", "missed_references",
    "unvalidated_detections", "failure_categories",
]
MATCHED_PAIR_REQUIRED = [
    "ref_index", "det_index", "ref_start", "ref_end",
    "det_start", "det_end",
]
UNVALIDATED_REQUIRED = ["det_index"]
REFERENCE_ENTRY_REQUIRED = [
    "plasmid_accession", "pdif_start", "pdif_end", "pdif_sequence",
    "strand", "orientation", "source_pmid", "source_description",
    "curator_notes", "curation_date",
]
METRICS_REQUIRED = [
    "detection_rate", "precision", "position_median_error",
    "position_mean_error", "position_max_error", "fpr_per_kb",
    "per_plasmid_breakdown", "failure_category_summary",
    "total_reference_sites", "total_matched_sites",
    "total_missed_sites", "total_unvalidated",
]
PER_PLASMID_REQUIRED = [
    "plasmid_accession", "reference_count", "matched_count",
    "missed_count", "unvalidated_count", "detection_rate",
    "failure_categories",
]
TOP_LEVEL_REQUIRED = [
    "mode", "timestamp", "plasmid_filter", "comparisons",
    "failures", "metrics", "output_dir",
]


def _validate_per_plasmid_entry(entry: Any, path: str, errors: Errors) -> None:
    """Validate one per-plasmid breakdown entry."""
    _check_type(entry, dict, path, errors)
    if not isinstance(entry, dict):
        return
    _check_required_keys(entry, PER_PLASMID_REQUIRED, path, errors)
    _check_extra_keys(entry, set(PER_PLASMID_REQUIRED), path, errors)

    for key, typ in [
        ("plasmid_accession", str),
        ("reference_count", int),
        ("matched_count", int),
        ("missed_count", int),
        ("unvalidated_count", int),
        ("detection_rate", (int, float)),
        ("failure_categories", dict),
    ]:
        if key in entry:
            _check_type(entry[key], typ, f"{path}.{key}", errors)
    # detection_rate in [0, 1]
    if isinstance(entry.get("detection_rate"), (int, float)):
        dr = entry["detection_rate"]
        if dr < 0.0 or dr > 1.0:
            errors.append(f"{path}.detection_rate: value {dr} out of range [0.0, 1.0]")


def _validate_matched_pair(pair: Any, path: str, errors: Errors) -> None:
    _check_type(pair, dict, path, errors)
    if not isinstance(pair, dict):
        return
    _check_required_keys(pair, MATCHED_PAIR_REQUIRED, path, errors)
    _check_extra_keys(pair, set(MATCHED_PAIR_REQUIRED), path, errors)
    for key in ("ref_index", "det_index", "ref_start", "ref_end", "det_start", "det_end"):
        if key in pair:
            _check_type(pair[key], int, f"{path}.{key}", errors)
            if isinstance(pair[key], int) and pair[key] < 0:
                errors.append(f"{path}.{key}: negative value {pair[key]}")


def _validate_unvalidated(det: Any, path: str, errors: Errors) -> None:
    _check_type(det, dict, path, errors)
    if not isinstance(det, dict):
        return
    _check_required_keys(det, UNVALIDATED_REQUIRED, path, errors)
    _check_extra_keys(det, set(UNVALIDATED_REQUIRED), path, errors)
    if "det_index" in det:
        _check_type(det["det_index"], int, f"{path}.det_index", errors)
        if isinstance(det["det_index"], int) and det["det_index"] < 0:
            errors.append(f"{path}.det_index: negative value {det['det_index']}")


def _validate_reference_entry(ref: Any, path: str, errors: Errors) -> None:
    _check_type(ref, dict, path, errors)
    if not isinstance(ref, dict):
        return
    _check_required_keys(ref, REFERENCE_ENTRY_REQUIRED, path, errors)
    _check_extra_keys(ref, set(REFERENCE_ENTRY_REQUIRED), path, errors)
    for key, typ in [
        ("plasmid_accession", str),
        ("pdif_start", int),
        ("pdif_end", int),
        ("pdif_sequence", str),
        ("strand", str),
        ("orientation", str),
        ("source_pmid", str),
        ("source_description", str),
        ("curator_notes", str),
        ("curation_date", str),
    ]:
        if key in ref:
            _check_type(ref[key], typ, f"{path}.{key}", errors)


def _validate_comparison(comp: Any, path: str, errors: Errors) -> None:
    """Validate one ComparisonResult."""
    _check_type(comp, dict, path, errors)
    if not isinstance(comp, dict):
        return
    _check_required_keys(comp, COMPARISON_REQUIRED, path, errors)
    _check_extra_keys(comp, set(COMPARISON_REQUIRED), path, errors)

    # Validate nested arrays
    for idx, pair in enumerate(comp.get("matched_pairs", [])):
        _validate_matched_pair(pair, f"{path}.matched_pairs[{idx}]", errors)
    for idx, ref in enumerate(comp.get("missed_references", [])):
        _validate_reference_entry(ref, f"{path}.missed_references[{idx}]", errors)
    for idx, det in enumerate(comp.get("unvalidated_detections", [])):
        _validate_unvalidated(det, f"{path}.unvalidated_detections[{idx}]", errors)

    # failure_categories must be dict[str, int]
    fc = comp.get("failure_categories", {})
    _check_type(fc, dict, f"{path}.failure_categories", errors)
    if isinstance(fc, dict):
        for k, v in fc.items():
            _check_type(k, str, f"{path}.failure_categories", errors)
            _check_type(v, int, f"{path}.failure_categories.{k}", errors)
            if isinstance(v, int) and v < 0:
                errors.append(f"{path}.failure_categories.{k}: negative count {v}")


def _validate_metrics(metrics: Any, path: str, errors: Errors) -> None:
    """Validate the Metrics sub-object."""
    _check_type(metrics, dict, path, errors)
    if not isinstance(metrics, dict):
        return
    _check_required_keys(metrics, METRICS_REQUIRED, path, errors)
    _check_extra_keys(metrics, set(METRICS_REQUIRED), path, errors)

    for key, typ in [
        ("detection_rate", (int, float)),
        ("precision", (int, float)),
        ("position_median_error", (int, float)),
        ("position_mean_error", (int, float)),
        ("position_max_error", (int, float)),
        ("fpr_per_kb", (int, float)),
        ("per_plasmid_breakdown", list),
        ("failure_category_summary", dict),
        ("total_reference_sites", int),
        ("total_matched_sites", int),
        ("total_missed_sites", int),
        ("total_unvalidated", int),
    ]:
        if key in metrics:
            _check_type(metrics[key], typ, f"{path}.{key}", errors)

    # Rate / precision in [0, 1]
    for rate_key in ("detection_rate", "precision"):
        if isinstance(metrics.get(rate_key), (int, float)):
            v = metrics[rate_key]
            if v < 0.0 or v > 1.0:
                errors.append(f"{path}.{rate_key}: value {v} out of range [0.0, 1.0]")

    # Integer counts >= 0
    for count_key in ("total_reference_sites", "total_matched_sites",
                      "total_missed_sites", "total_unvalidated"):
        if isinstance(metrics.get(count_key), int) and metrics[count_key] < 0:
            errors.append(f"{path}.{count_key}: negative value {metrics[count_key]}")

    # Validate per_plasmid_breakdown
    ppb = metrics.get("per_plasmid_breakdown", [])
    _check_type(ppb, list, f"{path}.per_plasmid_breakdown", errors)
    if isinstance(ppb, list):
        for idx, entry in enumerate(ppb):
            _validate_per_plasmid_entry(entry, f"{path}.per_plasmid_breakdown[{idx}]", errors)

    # Validate failure_category_summary
    fcs = metrics.get("failure_category_summary", {})
    _check_type(fcs, dict, f"{path}.failure_category_summary", errors)
    if isinstance(fcs, dict):
        for k, v in fcs.items():
            _check_type(v, int, f"{path}.failure_category_summary.{k}", errors)
            if isinstance(v, int) and v < 0:
                errors.append(f"{path}.failure_category_summary.{k}: negative count {v}")


def validate_benchmark_result(data: dict, path: str = "root") -> Errors:
    """Validate a BenchmarkResult (pdif_only style). Return list of error strings."""
    errors: Errors = []

    _check_type(data, dict, path, errors)
    if not isinstance(data, dict):
        return errors

    _check_required_keys(data, TOP_LEVEL_REQUIRED, path, errors)
    _check_extra_keys(data, set(TOP_LEVEL_REQUIRED), path, errors)

    # mode
    if "mode" in data:
        _check_type(data["mode"], str, f"{path}.mode", errors)
        if isinstance(data["mode"], str) and data["mode"] not in VALID_MODES:
            errors.append(f"{path}.mode: unexpected value '{data['mode']}' (expected {VALID_MODES})")

    # timestamp
    if "timestamp" in data:
        _check_type(data["timestamp"], str, f"{path}.timestamp", errors)

    # plasmid_filter
    if "plasmid_filter" in data:
        if data["plasmid_filter"] is not None:
            _check_type(data["plasmid_filter"], str, f"{path}.plasmid_filter", errors)

    # comparisons
    if "comparisons" in data:
        _check_type(data["comparisons"], list, f"{path}.comparisons", errors)
        if isinstance(data["comparisons"], list):
            for idx, comp in enumerate(data["comparisons"]):
                _validate_comparison(comp, f"{path}.comparisons[{idx}]", errors)

    # failures
    if "failures" in data:
        _check_type(data["failures"], list, f"{path}.failures", errors)
        if isinstance(data["failures"], list):
            for idx, f in enumerate(data["failures"]):
                _check_type(f, str, f"{path}.failures[{idx}]", errors)

    # metrics
    if "metrics" in data:
        _validate_metrics(data["metrics"], f"{path}.metrics", errors)

    # output_dir
    if "output_dir" in data:
        _check_type(data["output_dir"], str, f"{path}.output_dir", errors)

    return errors


# ---------------------------------------------------------------------------
# Validator for ComparisonSummary format  (full_pipeline_results.json style)
# ---------------------------------------------------------------------------

SUMMARY_REQUIRED = ["precision", "recall", "f1", "total_reference_sites",
                    "total_matched_sites", "total_missed_sites",
                    "total_unvalidated", "plasmids_processed", "note"]

PER_PLASMID_POSITIVE_REQUIRED = [
    "accession", "length", "amr_genes", "amr_gene_list",
    "pdif_only_detection_rate", "full_pipeline_detection_rate",
    "matched_sites_both_modes", "reference_sites",
]
PER_PLASMID_NEGATIVE_REQUIRED = [
    "accession", "length", "pdif_only_detection_rate", "reference_sites",
]

COMPARISON_TOP_LEVEL_REQUIRED = [
    "title", "timestamp", "blastn_available", "amr_database",
    "plasmids_total", "amr_positive_plasmids", "amr_negative_plasmids",
    "summary", "per_plasmid", "key_insights",
]


def _validate_summary_block(block: Any, path: str, errors: Errors) -> None:
    """Validate one mode summary block (pdif_only / full_pipeline)."""
    _check_type(block, dict, path, errors)
    if not isinstance(block, dict):
        return
    _check_required_keys(block, SUMMARY_REQUIRED, path, errors)

    for key, typ in [
        ("precision", (int, float)),
        ("recall", (int, float)),
        ("f1", (int, float)),
        ("total_reference_sites", int),
        ("total_matched_sites", int),
        ("total_missed_sites", int),
        ("total_unvalidated", int),
        ("plasmids_processed", int),
        ("note", str),
    ]:
        if key in block:
            _check_type(block[key], typ, f"{path}.{key}", errors)


def _validate_per_plasmid_positive(entry: Any, path: str, errors: Errors) -> None:
    """Validate one per-plasmid amr_positive entry."""
    _check_type(entry, dict, path, errors)
    if not isinstance(entry, dict):
        return
    _check_required_keys(entry, PER_PLASMID_POSITIVE_REQUIRED, path, errors)

    for key, typ in [
        ("accession", str),
        ("length", int),
        ("amr_genes", int),
        ("amr_gene_list", list),
        ("pdif_only_detection_rate", (int, float)),
        ("full_pipeline_detection_rate", (int, float)),
        ("matched_sites_both_modes", int),
        ("reference_sites", int),
    ]:
        if key in entry:
            _check_type(entry[key], typ, f"{path}.{key}", errors)


def _validate_per_plasmid_negative(entry: Any, path: str, errors: Errors) -> None:
    """Validate one per-plasmid amr_negative entry."""
    _check_type(entry, dict, path, errors)
    if not isinstance(entry, dict):
        return
    _check_required_keys(entry, PER_PLASMID_NEGATIVE_REQUIRED, path, errors)

    for key, typ in [
        ("accession", str),
        ("length", int),
        ("pdif_only_detection_rate", (int, float)),
        ("reference_sites", int),
    ]:
        if key in entry:
            _check_type(entry[key], typ, f"{path}.{key}", errors)


def validate_comparison_summary(data: dict, path: str = "root") -> Errors:
    """Validate a ComparisonSummary (full_pipeline style). Return list of error strings."""
    errors: Errors = []

    _check_type(data, dict, path, errors)
    if not isinstance(data, dict):
        return errors

    _check_required_keys(data, COMPARISON_TOP_LEVEL_REQUIRED, path, errors)

    for key, typ in [
        ("title", str),
        ("timestamp", str),
        ("blastn_available", bool),
        ("amr_database", str),
        ("plasmids_total", int),
        ("amr_positive_plasmids", int),
        ("amr_negative_plasmids", int),
        ("key_insights", list),
    ]:
        if key in data:
            _check_type(data[key], typ, f"{path}.{key}", errors)

    # summary
    summary = data.get("summary", {})
    _check_type(summary, dict, f"{path}.summary", errors)
    if isinstance(summary, dict):
        for mode_key in ("pdif_only", "full_pipeline"):
            if mode_key in summary:
                _validate_summary_block(
                    summary[mode_key], f"{path}.summary.{mode_key}", errors
                )
            else:
                errors.append(f"{path}.summary: missing '{mode_key}'")

    # per_plasmid
    per_plasmid = data.get("per_plasmid", {})
    _check_type(per_plasmid, dict, f"{path}.per_plasmid", errors)
    if isinstance(per_plasmid, dict):
        # amr_positive
        pos = per_plasmid.get("amr_positive", [])
        _check_type(pos, list, f"{path}.per_plasmid.amr_positive", errors)
        if isinstance(pos, list):
            for idx, entry in enumerate(pos):
                _validate_per_plasmid_positive(
                    entry, f"{path}.per_plasmid.amr_positive[{idx}]", errors
                )
        # amr_negative
        neg = per_plasmid.get("amr_negative", [])
        _check_type(neg, list, f"{path}.per_plasmid.amr_negative", errors)
        if isinstance(neg, list):
            for idx, entry in enumerate(neg):
                _validate_per_plasmid_negative(
                    entry, f"{path}.per_plasmid.amr_negative[{idx}]", errors
                )

    return errors


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def validate_file(json_path: str) -> bool:
    """Validate the JSON results file. Returns True if valid, False otherwise."""
    path = Path(json_path)

    if not path.exists():
        print(f"ERROR: File not found: {path}")
        return False

    try:
        with path.open(encoding="utf-8") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"ERROR: Invalid JSON in {path.name}: {e}")
        return False

    _check_type(data, dict, "root", [])
    if not isinstance(data, dict):
        print(f"ERROR: {path.name}: top-level value is not a JSON object")
        return False

    # Detect format: BenchmarkResult has 'metrics' key, ComparisonSummary has 'summary' key
    if "metrics" in data:
        print(f"  Detected format: BenchmarkResult (metrics-based)")
        errors = validate_benchmark_result(data)
    elif "summary" in data and "title" in data:
        print(f"  Detected format: ComparisonSummary (side-by-side)")
        errors = validate_comparison_summary(data)
    else:
        errors = [f"Unrecognised format: top-level keys = {sorted(data.keys())}"]
        print(f"ERROR: Unrecognised result format")
        for e in errors:
            print(f"  {e}")
        return False

    if errors:
        print(f"\nFound {len(errors)} validation error(s):")
        for err in errors:
            print(f"  {err}")
        return False
    else:
        print(f"All validations passed (0 errors).")
        return True


def main() -> int:
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <results.json> [<results.json> ...]")
        return 1

    all_passed = True
    for json_path in sys.argv[1:]:
        fname = Path(json_path).name
        print(f"\n{'=' * 60}")
        print(f"Validating: {fname}")
        print(f"{'=' * 60}")
        if not validate_file(json_path):
            all_passed = False

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
