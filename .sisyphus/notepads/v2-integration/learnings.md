# V2 Scanner Integration Learnings

## 2026-05-14

### Old pdifFinder Runner
- Requires `findResistanceGene` mocked via `unittest.mock.patch` to return `["1-1000"]` for pdif detection to run
- Requires placeholder `AMRgene1.txt` and `AMRgene.txt` files to avoid crashes in `make_sort()` and `changepdifname()`
- Must use `circle='true'` and call `pf.singleThread()` directly (not subprocess)
- Output is in `pdif_site1.txt` (7+ column tab-separated), not `pdif_site.txt` (which requires `finalFilter` to run)
- `pdif_pair.txt` may not be created if no valid pairs found; parse `pdif_site1.txt` as fallback

### V2 Scanner
- `scan_sequence()` returns list of dicts with canonicalized orientation, total_score, pdif_type
- PWM thresholds (0.6 XerC, 0.5 XerD) produce **many false positives** — on a 9kb plasmid, ~763-1906 sites
- DC proportion consistently ~50% (random pairing of independently-detected arms)
- Score range: ~9.2 to ~17.8. Real sites cluster at high scores (14+)
- 2bp dedup window is insufficient for the combinatorial explosion; larger windows or score filtering needed

### Comparison
- Old scanner: 8-17 sites per plasmid (seed-pair based, more conservative)
- New scanner: 763-1906 sites per plasmid (independent arm detection, combinatorial pairing)
- Old DC proportion: 20-50% (varies by plasmid)
- New DC proportion: ~48-51% (consistently near 50%)
- 3 new Hall plasmids (NC_010402, CP008708, FM210331) downloaded successfully but old scanner found no pdif pairs

### Next Steps
- V2 scanner needs PWM threshold tuning (try 0.8/0.7 for XerC/XerD) OR score-based filtering
- Wider dedup window (10-20bp) to collapse overlapping candidates into single sites
- Consider greedy best-pair selection instead of all-pairs combinatorial assembly
