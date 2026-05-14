# PdifFinder Package

## OVERVIEW
Core algorithm, CLI, threading, BLAST orchestration, and visualization. 6 Python modules, 1565 lines total.

## STRUCTURE
```
PdifFinder/
├── pdifFinder.py           # Core: CLI, sliding window, BLAST, threading (1126 lines)
├── angularPlasmid.py       # Circular plasmid HTML/SVG viz (199 lines, standalone-capable)
├── pdifmodulecharts.py     # pdif-ARG module SVG charts (120 lines)
├── echarts.py              # ECharts feature rendering (120 lines, partially disabled)
├── Getpdifmoduleseq.py     # DEAD CODE — hardcoded Windows paths (12 lines)
├── AMRDB/                  # sequences.fasta + BLAST .n* index files
└── data/                   # redundant.seed.fa, pdifdatabase.fasta, genecolor.txt
```

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Seed parsing (4-entry groups) | `pdifFinder.py:207-212` | XerC[0:11], XerD[17:28], 4-entry layout |
| Circular vs linear limit | `pdifFinder.py:447` | `limit = seq_len if circular else seq_len - 27` |
| BLAST identity/coverage defaults | `pdifFinder.py:613-614` | `default_ident=90`, `default_cov=60` |
| Original thread (benchmark only) | `pdifFinder.py:395` | `findMatchFragmentThreadOriginal()` |
| Angular plasmid standalone | `angularPlasmid.py:199-200` | `if __name__ == "__main__"` |
| ISfinder commented-out refs | `angularPlasmid.py:16`, `echarts.py:9` | Do NOT remove |
| Dead code | `Getpdifmoduleseq.py` | Hardcoded `E:/Python/...`. Ignore entirely. |

## CONVENTIONS

- Mixed snake_case/camelCase: `arg_parse()` vs `findPdif()`, `inFile`. Match surrounding style.
- Threaded file writes inside `with fragment_lock:` / `with pair_lock:` only.
- Hardcoded thresholds: XerC=3, XerD=2, pair CD=4, DC=6. Do NOT change without domain justification.
- Circular: modulo indexing (`i % seq_len`). Linear: stop at `len(seq) - 27`.
- BLAST defaults: identity=90, coverage=60. Change only with domain knowledge.
- No f-strings, no type hints (Python 3.5 compat). Use `%` formatting.
- Flat layout (camelCase dir). Renaming breaks `setup.py` console_scripts.
- `findMatchFragmentThreadOriginal()` preserved for benchmark comparison only. DO NOT call in production.

## ANTI-PATTERNS

- DO NOT add f-strings — breaks Python 3.5/3.6 compat
- DO NOT remove `fragment_lock` / `pair_lock` — thread safety
- DO NOT change 4-entry seed parsing — algorithm depends on it
- DO NOT modify `findMatchFragmentThreadOriginal()` — benchmark only
- DO NOT add code to `Getpdifmoduleseq.py` — dead code, ignore
- DO NOT add significant code to `pdifFinder.py` — 1126 lines, needs decomposition
- DO NOT remove commented-out ISfinder refs — may be resurrected
- DO NOT use `shell=True` outside BLAST subprocess calls
