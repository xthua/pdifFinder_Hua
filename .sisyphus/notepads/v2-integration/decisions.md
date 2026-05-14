# V2 Scanner Integration Decisions

## 2026-05-14

### Parse pdif_site1.txt as fallback
`parse_old_pdif_output()` first tries `pdif_site.txt` (9-column format after `changepdifname/finalFilter`), then falls back to `pdif_site1.txt` (7-column raw format). This avoids depending on `pdif_pair.txt` being created.

### Direct singleThread call (not subprocess)
Old scanner is invoked via `pf.singleThread()` with mocked `findResistanceGene` instead of running `pdifFinder --inFile ...` as subprocess. This is needed because the old scanner's `main()` creates nested directories and uses `split()`/`process()` flow that doesn't work with single FASTA files.

### 2bp dedup window (as specified)
Deduplication uses 2bp merge window as specified in task. Noted that this is insufficient for the v2 scanner's current sensitivity level — most candidates are 3-4bp apart. Score-based filtering or wider windows would be needed for practical use.

### Module-level PWM cache
`_get_pwms()` caches PWMs at module level to avoid rebuilding on every scan. Uses `load_pwms()` from existing `xer_pwm.json`.
