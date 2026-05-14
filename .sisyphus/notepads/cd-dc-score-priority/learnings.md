# cd-dc-score-priority Learnings

## Changes Made
1. **findMatchFragmentThread** (line 488): Added `total_mm = mismatches_c + mismatches_d` and appended `|{score}` to output format
2. **findPdif parsing** (lines 256-264): Changed `posOrientDict` from `dict[int, set[str]]` to `dict[int, dict[str, int]]`, stores score via `min(score, existing)`
3. **findPdif orientation selection** (lines 278-282): Replaced CD-hardcoded priority with `min(orientations, key=orientations.get)` — picks orientation with lowest mismatch score

## Test Impact
- Output format changed from `pos|orient` to `pos|orient|score` — broke `test_whole_plasmid.py` which parsed with `int(line.strip())`
- Fixed all 3 occurrences in `test_whole_plasmid.py` to use `int(line.strip().split('|')[0])`
- Also fixed `positions_orig` and `positions_opt` variants

## Pre-existing Issues (not caused by this change)
- `test_detection_rate_above_threshold` (0.721 < 0.80) — cached benchmark data issue
- `test_precision_above_threshold` (0.018 < 0.70) — cached benchmark data issue
- `shutil.rmtree` OSError on container overlayfs — blocks `-g` flag CLI testing

## Verification
- All 29 non-benchmark tests pass
- All 6 test_whole_plasmid tests pass
- LSP diagnostics: only pre-existing pyright warnings + false positives on `mismatches_c`/`mismatches_d` (guaranteed bound when `passed=True`)
- Direct unit tests confirm min-score orientation selection works correctly
