
## 2026-05-10: Plan has wrong accession for pS30-1

### Problem
The plan specifies GenBank accession KU987654 for Blackwell & Hall (2017) pS30-1, but KU987654 is a plant ITS sequence (Prangos denticulata, 587bp), not a bacterial plasmid.

### Root cause
The plan was written with an incorrect accession number. The correct accession is KY617771.

### Resolution
Identified the correct accession (KY617771) by:
1. Web search for "Blackwell Hall 2017 pS30-1 GenBank accession"
2. Paper text: "The sequence of pS30-1 has been submitted to GenBank under accession number KY617771"
3. NCBI eutils validation - KY617771 returns pS30-1 18,234bp plasmid

### Impact
None - downstream tasks should reference KY617771.1, not KU987654. KY617771.1 is already present in plasmid_accessions.json.

---

## pdif_only Benchmark Run Issues — 2026-05-10

### 1. angularPlasmid hardcoded filename
**Status**: WORKED AROUND (symlink)
pdifFinder `angularPlasmid()` at line 85 hardcodes `fastaFile = outdir + '/inputFile.fasta'`. 
The benchmark writes FASTA as `{accession}_1.fasta`. This causes FileNotFoundError.
**Fix applied**: Create symlink `inputFile.fasta -> {acc}_1.fasta` in run directory (benchmark_literature.py:706-708).

### 2. changepdifname doesn't create AMRgene.txt for empty input
**Status**: WORKED AROUND (pre-touch)
`changepdifname()` reads `AMRgene1.txt` and only writes `AMRgene.txt` inside the loop body.
When `AMRgene1.txt` is empty (mocked pdif_only mode), the output file is never created.
`Getpdifmoduleseq()` then crashes with FileNotFoundError for `AMRgene.txt`.
**Fix applied**: Pre-touch `AMRgene.txt` alongside `AMRgene1.txt` (benchmark_literature.py:276-279).

### 3. failure_category_summary is empty
**Status**: OPEN
All per-plasmid entries have `failure_categories: {}` despite 12 missed references.
The `categorize_failure()` method in the benchmark code may not be producing categories,
or the categories aren't being propagated into the results. Need investigation.

### 4. fpr_per_kb returns -1.0
**Status**: OPEN
The false-positive-rate-per-kb metric returns -1.0, which appears to be a sentinel value
indicating the metric wasn't computed. Likely because plasmid lengths aren't passed to
the metric computation in pdif_only mode, or the FPR formula requires data not available
in this mode.

### 5. High false positive rate (low precision)
**Status**: OPEN
Precision is 1.76% with 1733 unvalidated detections across 12 plasmids.
This suggests pdifFinder's seed-based detection produces many false positives.
Some may be genuine pdif-like sequences not in the literature reference,
but the scale (1733 vs 31 true positives) indicates poor specificity.
## F2 — Code Quality Review (2026-05-10)

### Files Reviewed
| # | File | Lines | LSP | Verdict |
|---|------|-------|-----|---------|
| 1 | `tests/benchmark_literature.py` | 855 | CLEAN | 2 minor issues |
| 2 | `tests/helpers.py` | 372 | CLEAN | 2 minor issues |
| 3 | `tests/validation/test_literature_benchmark.py` | 189 | 1 err (Pyright pytest import) | 1 minor + 1 config issue |
| 4 | `tests/benchmark_data/validate_reference.py` | 172 | CLEAN | CLEAN |
| 5 | `tests/benchmark_data/validate_results.py` | 528 | CLEAN | 1 minor issue |

### Issues Found (all Low severity)

1. **benchmark_literature.py:71-75** — `FailureCategory` dataclass defined but NEVER instantiated. `categorize_failure()` returns `str`, making this dataclass dead code.
2. **benchmark_literature.py:2,50,72** — Stale docstrings: "Skeleton", "stub" (code is fully implemented, not skeleton)
3. **helpers.py:11** — `from pathlib import Path` — unused import
4. **helpers.py:334,349** — Bare `except Exception` in `cleanup_temp_files()` and `cleanup_temp_dirs()`
5. **validate_results.py:30-36** — `_load_schema()` defined but never called in validation flow

### Configuration
- **pytest.ini**: `benchmark` marker IS registered (line 21) but pytest 9.0.2 emits warnings anyway — likely `[tool:pytest]` vs `[pytest]` header mismatch. Non-blocking.

### parse_pdif_output Fix Assessment
- **CLEAN**: `tests/helpers.py:165-244` correctly parses 9-column `pdif_site.txt` format
- All 9 fields extracted: sequence_id, pdif_name, start, end, xerC, spacer, xerD, orientation, pdif_db_name
- Backward-compatible: `pdif_site` computed as `xerC + spacer + xerD`
- AMR gene parsing (7-col) also correct
- pdif module parsing handles variable-width columns

### Test Results
- **5/7 passed**, 2 failed (THRESHOLD assertions, NOT code bugs):
  - `test_detection_rate_above_threshold`: 0.721 < 0.80 target
  - `test_precision_above_threshold`: 0.018 < 0.70 target
- These reflect pdifFinder's actual performance against the benchmark data, not test harness defects.

### No Issues Found
- No todo/fixme/hack stubs in any file
- No empty catch blocks
- No over-abstraction
- `validate_reference.py`: completely clean
