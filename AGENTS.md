# PdifFinder Knowledge Base

**Generated:** 2026-05-11
**Commit:** d2131af
**Branch:** master

## OVERVIEW

PdifFinder v1.1.1 — Python bioinformatics CLI for annotating antimicrobial resistance (AMR) genes, 28bp pdif recombination sites, and pdif-ARG modules in bacterial plasmids. Uses sliding-window algorithm + BLAST-based AMR annotation. GPLv3.

## STRUCTURE

```
pdifFinder_Hua/
├── PdifFinder/             # Installable package (camelCase — deviation from PEP 8)
│   ├── pdifFinder.py       # ⬅ Core: CLI, BLAST, sliding window, threading (1126 lines)
│   ├── angularPlasmid.py   # Circular plasmid HTML/SVG visualization
│   ├── pdifmodulecharts.py # pdif-ARG module SVG chart
│   ├── echarts.py          # ECharts-based feature rendering
│   ├── Getpdifmoduleseq.py # DEAD CODE — hardcoded Windows paths, not used
│   ├── AMRDB/              # BLAST nucleotide DB for resistance genes
│   └── data/               # Reference data (redundant.seed.fa, pdifdatabase.fasta)
├── tests/                  # Pytest suite (15+ files, coverage-enabled)
│   ├── helpers.py          # Core test utilities (371 lines)
│   ├── conftest.py         # Shared fixtures
│   ├── benchmark_literature.py # Literature benchmark (923 lines)
│   ├── integration/        # End-to-end pipeline tests
│   ├── validation/         # Real-pdif validation + literature benchmark tests
│   └── benchmark_data/     # Cached results, reference CSV, plasmid fixtures
├── docs/                   # Algorithm design and analysis docs
├── scripts/                # extract_plasmid_accessions.py
├── setup.py                # Legacy setuptools (no pyproject.toml)
├── requirements.txt        # Pinned: biopython==1.78, matplotlib, numpy, pandas
├── pytest.ini              # Coverage + custom markers
├── debug_*.py              # Ad-hoc debug scripts (project root)
├── test_*.py               # Ad-hoc test scripts (project root — duplicates tests/ patterns)
└── profile_scanning.py     # cProfile perf profiling script
```

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| CLI entry / arg parsing | `PdifFinder/pdifFinder.py:main()` L1065 | Console script: `pdifFinder` |
| Sliding window pdif detection | `PdifFinder/pdifFinder.py:findMatchFragmentThread()` | 28bp window, XerC[0:11], XerD[17:28] |
| Threading / locks | `PdifFinder/pdifFinder.py` L19-20 | `fragment_lock`, `pair_lock` |
| BLAST AMR gene detection | `PdifFinder/pdifFinder.py:checkResistanceGenePos()` | Subprocess `blastn`, requires PATH |
| Circular plasmid visualization | `PdifFinder/angularPlasmid.py` | Also runnable standalone |
| Module SVG charts | `PdifFinder/pdifmodulecharts.py` | Gene color mapping in `data/genecolor.txt` |
| Test utilities | `tests/helpers.py` | `run_pdif_finder()`, `parse_pdif_output()`, cleanup |
| Literature benchmark | `tests/benchmark_literature.py` | `LiteratureBenchmark` class, pdif_only & full_pipeline modes |
| Reference data | `tests/benchmark_data/literature_reference.csv` | 43 curated pdif sites, 3 published studies |
| Seed database | `PdifFinder/data/redundant.seed.fa` | Seed pairs used by algorithm |

## CODE MAP

| Symbol | Type | Location | Refs | Role |
|--------|------|----------|------|------|
| `main()` | function | `PdifFinder/pdifFinder.py:1065` | Entry point | CLI orchestrator |
| `arg_parse()` | function | `PdifFinder/pdifFinder.py:24` | — | Argument parsing |
| `findPdif()` | function | `PdifFinder/pdifFinder.py` | — | Top-level pdif detection pipeline |
| `findMatchFragmentThread()` | function | `PdifFinder/pdifFinder.py` | Threaded | Sliding window seed matching |
| `findPossiblePairThread()` | function | `PdifFinder/pdifFinder.py` | Threaded | Pair detection (XerC+XerD) |
| `angularPlasmid()` | function | `PdifFinder/angularPlasmid.py` | Imported | Circular plasmid HTML gen |
| `pdifmoduleechart()` | function | `PdifFinder/pdifmodulecharts.py` | Imported | Module SVG gen |
| `LiteratureBenchmark` | class | `tests/benchmark_literature.py` | Test | Full benchmark harness |
| `run_pdif_finder()` | function | `tests/helpers.py:126` | Test utility | Subprocess pdifFinder runner |
| `parse_pdif_output()` | function | `tests/helpers.py:165` | Test utility | Output parser (6 file types) |

## CONVENTIONS

**Naming**: INCONSISTENT. Snake_case (`arg_parse`) and camelCase (`findPdif`, `inFile`, `seqList`) mixed in same file. No enforced standard.

**Git commits**: Conventional commits — `type(scope): description`. Types: `fix`, `perf`, `docs`, `build`. Example: `fix(pdifFinder): repair XerC/XerD scanning algorithm`

**Thread safety**: All threaded file writes MUST use `with fragment_lock:` or `with pair_lock:`. See `PdifFinder/pdifFinder.py` L19-20.

**Pytest markers**: All custom markers MUST be registered in `pytest.ini`. Currently: `slow`, `integration`, `unit`, `data`, `db`, `benchmark`. Note: `performance` marker used in code but NOT registered — will warn.

**Dependencies**: Pinned to EXACT versions in both `setup.py` and `requirements.txt`. Never loosen pins without testing.

**External dependency**: `blastn` (BLAST+ 2.10.1+) required in PATH. NOT a pip dependency.

**No type hints**: Zero type annotations. No mypy. Do NOT add type hints without team consensus — this is a deliberate style choice.

**String formatting**: `%` formatting and `+` concatenation. NO f-strings (Python 3.5/3.6 compat).

**Package layout**: Flat (not src-layout). Package dir is `PdifFinder/` (camelCase — historical, do NOT rename).

## ANTI-PATTERNS (THIS PROJECT)

- **DO NOT** add f-strings — breaks Python 3.5/3.6 compat declared in setup.py
- **DO NOT** remove `fragment_lock` or `pair_lock` — thread safety depends on them
- **DO NOT** touch `findMatchFragmentThreadOriginal()` — preserved for benchmark comparison
- **DO NOT** rename `PdifFinder/` to lowercase — breaks console_scripts entry point
- **DO NOT** add `pyproject.toml` without also updating setup.py — dual-config maintenance required
- **DO NOT** use `shell=True` in subprocess outside BLAST calls — security concern
- **AVOID** bare `except:` — found at `pdifFinder.py:849`, use specific exceptions
- **TEST FILES**: `test_*.py` must use `test_` function prefix for pytest discovery. Root-level ad-hoc scripts are NOT pytest-discoverable (they use `if __name__` guards)

## UNIQUE STYLES

- **CamelCase package**: `PdifFinder` (not `pdif_finder`). Historical, matches PyPI name.
- **Mixed naming**: Both camelCase and snake_case coexist. Match surrounding code style when editing a file.
- **Dual README**: `README.md` (GitHub) + `README.rst` (PyPI). Update both.
- **Commented-out ISfinder code**: Left intentionally in `angularPlasmid.py` and `echarts.py`. Do NOT remove — may be un-commented in future release.
- **Source-code-level tests**: `test_bug_fixes.py` and `test_thread_safety.py` use `inspect.getsource()` to verify source code, not runtime behavior.

## COMMANDS

```bash
# Install
pip install .

# Run (requires BLAST+ in PATH)
pdifFinder -i input.fasta -o output_dir

# Tests
pytest                              # All non-benchmark tests
pytest -m benchmark                 # Literature benchmark tests
pytest tests/validate_scanning.py   # Single file

# Benchmarks
python -m tests.benchmark_literature --mode pdif_only
python -m tests.benchmark --compare-before-after

# Profile
python profile_scanning.py

# Lint (add when configured — currently none)
# python -m ruff check PdifFinder/
```

## NOTES

- **BLAST+ must be in PATH** — `blastn` invoked via subprocess. Not pip-installable.
- **Circular genome scanning**: Uses modulo indexing. `--circular-seq=auto|true|false` flag.
- **Mismatch thresholds**: Hardcoded — XerC: 3, XerD: 2, Pair CD: 4, Pair DC: 6. Not configurable.
- **Performance**: ~45% speedup from early-exit strategies. Deferred: Boyer-Moore, seed clustering, bitmask.
- **Open bugs**: bare `except:`, `shell=True`, empty `failure_category_summary`, `fpr_per_kb` -1.0 sentinel.
- **iOS/debug scripts at root**: `debug_*.py`, `test_*.py`, `profile_scanning.py` are ad-hoc dev tools. Not part of pytest suite.
