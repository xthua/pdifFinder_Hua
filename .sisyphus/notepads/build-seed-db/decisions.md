# build-seed-db Decisions

## Orientation Detection: Structural vs. Scoring
**Decision**: Use structural TTAT-based orientation detection instead of consensus scoring.
**Rationale**: XerD arms are characterized by "TTAT" at the first 4 positions (positions 1-4 of the 11bp arm). This is a strong, unambiguous signal. Scoring-based approaches (comparing XerC/XerD consensus) failed because some XerD arms resemble XerC consensus patterns and vice versa, leading to incorrect flips.

## RC Dedup Strategy
**Decision**: After orientation and standard dedup, remove reverse-complement duplicates.
**Rationale**: Some sites have both the forward and RC sequences in the source data (e.g., cameranesi CSV gives forward strand, literature_reference gives feature strand). After orientation, both can appear CD with different CRs. Keeping both creates redundancy without adding biological value. Keep the one with higher validation status.

## CR Grouping Threshold: hd ≤ 1 (not ≤ 2)
**Decision**: Use Hamming distance ≤ 1 for merging CR groups.
**Rationale**: With 6bp CR strings, hd ≤ 2 creates excessive transitive chaining (16 different CRs collapsing into one group). With hd ≤ 1, groups are more biologically coherent while still allowing single-base variants to merge.

## Consensus: Majority Rule with Pseudocount
**Decision**: Use majority-rule consensus with pseudocount 0.1 for PWM generation.
**Rationale**: Standard bioinformatics practice — pseudocounts prevent zero probabilities and allow detection of rare variants while majority rule produces clean consensus sequences for single-group seeds.

## GenBank Parsing
**Decision**: Use Bio.SeqIO to parse GenBank features, extracting sequences from features with /note containing "XerC" or "XerD".
**Rationale**: Properly handles complement locations via BioPython's feature.extract(), which automatically reverse-complements when location.strand == -1.
