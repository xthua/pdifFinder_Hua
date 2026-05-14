# pdif_v2_scanner.py - Implementation Notes

## Decisions

### PWM Structure
- Each PWM is `{nt: [11 probabilities]}` where nt in {A,C,G,T}
- Pseudocount 0.01 added per cell before column normalization: `p(nt,i) = (count + 0.01) / (N + 0.04)`
- Degenerate XerD PWM adds 3 variant 11-mers with weight 0.5 each to the frequency table before normalization
- Saved as JSON at `PdifFinder/data/xer_pwm.json`

### Arm Scanner (scan_xer_arm)
- XerC thresholds: PWM >= 0.6 * max_score OR hamming <= 2
- XerD thresholds: PWM >= 0.5 * max_score OR hamming <= 3
- Max possible PWM score computed as sum of column maxes

### Pair Assembly (assemble_pdif_sites)
- Separation between arm start positions: 16-22 bp (requirements specify this range)
- Both CD (XerC-then-XerD) and DC (XerD-then-XerC) orientations assembled
- O(n) using position-indexed dictionaries

### Dual-Strand Scanning
- Reverse complement: scan revcomp strand with same PWMs, then map coordinates back
- Forward mapping: fwd_start = seq_len - (rev_start + 28)
- Orientation flips: CD on revcomp → DC on forward, DC on revcomp → CD on forward

### Canonical Normalization
- Canonical form is ALWAYS 5'-XerC-spacer-XerD-3' (CD orientation)
- DC sites are reverse-complemented to produce canonical representation
- Both raw and canonical orientation tracked

### Scoring
- total_score = xerC_pwm + xerD_pwm + spacer_bonus(0.1 for 6bp) + paired_bonus(0.2) + context_bonus(0.1 within 100bp)
- Paired_bonus always 0.2 since both arms must be independently detected to form a pair

### Degenerate Rescue
- Only activates when feature_boundaries is non-empty
- Rescans positions where XerC score >= 0.6*max but no pair formed
- Requires position within 200bp of a feature boundary
- Uses degenerate XerD PWM with threshold 0.4*max OR hamming <= 3

## Output Format (per scan_sequence)
```
{'contig': str, 'start': int, 'end': int, 'strand': '+'|'-',
 'raw_orientation': 'CD'|'DC', 'canonical_orientation': 'CD',
 'xerC_seq': 11bp, 'spacer_seq': 6bp, 'xerD_seq': 11bp,
 'total_score': float, 'pdif_type': 'canonical'|'degenerate',
 'rescue_flag': 'standard'|'degenerate'}
```
