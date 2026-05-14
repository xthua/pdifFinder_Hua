# Test Suite
## OVERVIEW
Pytest suite for PdifFinder — 15+ files across unit, integration, validation, benchmark. Coverage-enabled, custom markers, literature-based ground truth (43 curated pdif sites, 3 studies).

## STRUCTURE
```
tests/
├── conftest.py              # Shared fixtures (test_data_dir, sample_fasta/gb paths)
├── helpers.py               # Seq gen, pdifFinder runner, output parser (371 lines)
├── test_framework.py        # Smoke: pytest + import check
├── test_bug_fixes.py        # Source-level regressions (inspect.getsource)
├── test_performance.py      # Speedup >=1.2x threshold
├── test_thread_safety.py    # Lock validation via source inspection
├── test_whole_plasmid.py    # Whole-plasmid scanning + findFeatureEndSeq
├── validate_output.py       # Output existence + column format (605 lines)
├── validate_scanning.py     # Algorithm coverage: basic/edge/large
├── benchmark.py             # Performance benchmark CLI
├── benchmark_literature.py  # Literature benchmark, 2 modes (923 lines)
├── integration/
│   └── test_pdif_finder_pipeline.py
├── validation/
│   ├── test_against_real_pdif.py     # DR + FPR (xfail for design limits)
│   └── test_literature_benchmark.py  # DR>=0.80, precision>=0.70
├── benchmark_data/          # Reference CSV, cached results, scripts
└── test_data/               # Synthetic sequences: positive, negative, edge_case
```
## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| Run pdifFinder | `helpers.py:run_pdif_finder()` | Subprocess wrapper |
| Parse output | `helpers.py:parse_pdif_output()` | 9-col pdif_site.txt + 7-col AMRgene.txt |
| Create sequences | `helpers.py:create_sequence_with_pdif_sites()` | Injects known pdif sites |
| Shared fixtures | `conftest.py` | test_data_dir, sample_fasta/gb |
| Benchmark harness | `benchmark_literature.py:LiteratureBenchmark` | pdif_only / full_pipeline |
| Ground truth | `benchmark_data/literature_reference.csv` | 43 sites, 3 studies |
| Known limitations | `test_bug_fixes.py`, `validation/test_against_real_pdif.py` | xfail markers |
| Standalone test | `helpers.py` __main__ block | Quick utility check |
## CONVENTIONS
- `test_` prefix on functions; `Test` prefix on classes.
- Registered markers (pytest.ini): `slow`, `integration`, `unit`, `data`, `db`, `benchmark`. Register new BEFORE use.
- `performance` marker used in code but NOT registered — avoid until registered.
- Coverage: `--cov=PdifFinder` (prod only). No `--cov-fail-under`.
- Benchmark: `pytest -m benchmark` (run separately).
- Mocking: `unittest.mock.patch` for `findResistanceGene()`. Source inspection for bug fix verification.
- Two test styles: pytest-discoverable (tests/) + root-level ad-hoc scripts (not discovered).
## ANTI-PATTERNS
- Do NOT add unregistered pytest markers — add to pytest.ini first.
- Do NOT use `performance` marker until registered.
- Do NOT delete xfail markers — they document known limitations.
- Do NOT move root-level test scripts into tests/ — they are ad-hoc debug tools, not CI tests.
