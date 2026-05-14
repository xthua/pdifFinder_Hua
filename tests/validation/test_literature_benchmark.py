#!/usr/bin/env python3
"""
Literature benchmark validation tests.

Verifies benchmark results against expected thresholds. All tests use the
``@pytest.mark.benchmark`` marker so they are excluded from the default
test suite (run with ``-m benchmark``).

Tests 4–7 depend on cached ``pdif_only_results.json`` and are skipped if
the file does not exist.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from tests.benchmark_data.validate_reference import validate as validate_reference_csv
from tests.benchmark_literature import LiteratureBenchmark

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BENCHMARK_DATA = PROJECT_ROOT / "tests" / "benchmark_data"
RESULTS_DIR = BENCHMARK_DATA / "results"
PDF_ONLY_RESULTS = RESULTS_DIR / "pdif_only_results.json"
REFERENCE_CSV = BENCHMARK_DATA / "literature_reference.csv"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def pdif_only_results_path() -> Path:
    """Path to cached pdif-only benchmark results."""
    return PDF_ONLY_RESULTS


@pytest.fixture
def pdif_only_results() -> dict | None:
    """Load cached pdif-only benchmark results as a dict, or ``None``."""
    if PDF_ONLY_RESULTS.exists():
        with PDF_ONLY_RESULTS.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    return None


# ---------------------------------------------------------------------------
# Test 1 — Reference CSV validation
# ---------------------------------------------------------------------------


@pytest.mark.benchmark
def test_literature_reference_csv_valid() -> None:
    """Validate ``literature_reference.csv`` against the schema."""
    assert REFERENCE_CSV.exists(), f"Reference CSV not found: {REFERENCE_CSV}"
    is_valid = validate_reference_csv(str(REFERENCE_CSV))
    assert is_valid, "Literature reference CSV failed validation"


# ---------------------------------------------------------------------------
# Test 2 — Benchmark harness initialisation
# ---------------------------------------------------------------------------


@pytest.mark.benchmark
def test_benchmark_harness_runs(tmp_path: Path) -> None:
    """Verify ``LiteratureBenchmark`` initialises without error."""
    bench = LiteratureBenchmark(mode="pdif_only", output_dir=str(tmp_path))
    assert bench is not None
    assert bench.mode == "pdif_only"
    assert bench._reference_csv == REFERENCE_CSV


# ---------------------------------------------------------------------------
# Test 3 — Results file existence
# ---------------------------------------------------------------------------


@pytest.mark.benchmark
def test_pdif_only_results_exist() -> None:
    """Verify ``pdif_only_results.json`` exists on disk.

    Fails with a clear message directing the user to run the benchmark first.
    """
    assert PDF_ONLY_RESULTS.exists(), (
        f"pdif_only_results.json not found at {PDF_ONLY_RESULTS}. "
        "Run the literature benchmark first to generate results."
    )


# ---------------------------------------------------------------------------
# Test 4 — Detection rate threshold
# ---------------------------------------------------------------------------


@pytest.mark.benchmark
@pytest.mark.skipif(
    not PDF_ONLY_RESULTS.exists(),
    reason="pdif_only_results.json not available — run literature benchmark first",
)
def test_detection_rate_above_threshold(pdif_only_results: dict) -> None:
    """Verify ``detection_rate`` >= 0.80 in benchmark results."""
    metrics = pdif_only_results.get("metrics", {})
    # Support both the current Metrics dataclass key and the legacy key.
    dr = metrics.get("detection_rate", metrics.get("recall"))
    assert dr is not None, (
        "Neither 'detection_rate' nor 'recall' found in metrics"
    )
    assert dr >= 0.80, f"Detection rate {dr:.3f} < 0.80 threshold"


# ---------------------------------------------------------------------------
# Test 5 — Precision threshold
# ---------------------------------------------------------------------------


@pytest.mark.benchmark
@pytest.mark.skipif(
    not PDF_ONLY_RESULTS.exists(),
    reason="pdif_only_results.json not available — run literature benchmark first",
)
def test_precision_above_threshold(pdif_only_results: dict) -> None:
    """Verify ``precision`` >= 0.70 in benchmark results."""
    metrics = pdif_only_results.get("metrics", {})
    precision = metrics.get("precision")
    assert precision is not None, "'precision' not found in metrics"
    assert precision >= 0.70, f"Precision {precision:.3f} < 0.70 threshold"


# ---------------------------------------------------------------------------
# Test 6 — Position median error threshold
# ---------------------------------------------------------------------------


@pytest.mark.benchmark
@pytest.mark.skipif(
    not PDF_ONLY_RESULTS.exists(),
    reason="pdif_only_results.json not available — run literature benchmark first",
)
def test_position_median_error(pdif_only_results: dict) -> None:
    """Verify ``position_median_error`` <= 5 bp in benchmark results."""
    metrics = pdif_only_results.get("metrics", {})
    err = metrics.get("position_median_error")
    assert err is not None, "'position_median_error' not found in metrics"
    assert err <= 5, f"Position median error {err} > 5 bp threshold"


# ---------------------------------------------------------------------------
# Test 7 — Results JSON structure
# ---------------------------------------------------------------------------


@pytest.mark.benchmark
@pytest.mark.skipif(
    not PDF_ONLY_RESULTS.exists(),
    reason="pdif_only_results.json not available — run literature benchmark first",
)
def test_results_json_structure(pdif_only_results: dict) -> None:
    """Verify all required top-level and metrics fields exist."""
    # ---- Top-level fields ----
    for field in ("mode", "timestamp", "metrics", "comparisons", "failures", "output_dir"):
        assert field in pdif_only_results, (
            f"Missing top-level field '{field}'"
        )

    # ---- Metrics fields ----
    metrics = pdif_only_results["metrics"]
    required_metric_fields = (
        "detection_rate",
        "precision",
        "position_median_error",
        "position_mean_error",
        "position_max_error",
        "total_reference_sites",
        "total_matched_sites",
        "total_missed_sites",
        "total_unvalidated",
    )
    for field in required_metric_fields:
        assert field in metrics, f"Missing metrics field '{field}'"
