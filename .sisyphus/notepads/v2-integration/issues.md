# V2 Scanner Integration Issues

## 2026-05-14

### Resolved
1. **Old scanner returns 0 sites** — Fixed by mocking `findResistanceGene` and providing placeholder files
2. **Module import fails without PYTHONPATH** — Script must be run with `PYTHONPATH=.` from project root
3. **Hall plasmid download** — NC_010402, CP008708, FM210331 downloaded successfully via NCBI efetch

### Unresolved
1. **V2 scanner over-detection** — PWM thresholds (0.6/0.5) cause combinatorial explosion. On 9kb plasmids: 763-1906 sites. Comparison numbers are inflated but DC proportion measurement (~50%) is still valid.
2. **Old scanner fails on Hall plasmids** — NC_010402 (p2ABAYE), CP008708 (p2AB5075), FM210331 (pABV01) found seed matches (8-9) but no valid pairs. May need different mock or genuine BLAST hits.
