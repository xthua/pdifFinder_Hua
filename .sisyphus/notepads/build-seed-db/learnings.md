# build-seed-db Learnings

## Pipeline: build_seed_db.py

### Data Collection
- 6 source types: NAR2025 (2), Cameranesi2018 (17), Blackwell2017 (8), GenBank annotations (17), SeedDB (41), LiteratureReference (47)
- 132 raw → 53 unique → 43 after RC-dedup → 46 final (43 + 3 novel consensus)
- GenBank and LiteratureReference significantly overlap with Cameranesi/Blackwell CSV data
- SeedDB entries mostly redundant with literature-curated data (41 → 4 unique)

### Orientation Key Learning
- **TTAT at positions 17-20 is the definitive CD-orientation signal**: XerD arms almost always start with TTAT (positions 1-4 of the 11bp XerD arm)
- Scoring-based orientation (comparing XerC/XerD consensus patterns) failed because XerC and XerD arms can score similarly on both sequences and their RCs
- The structural rule "if right side has TTAT and left doesn't → CD; if left has TTAT and right doesn't → flip" handles 50/53 cases correctly
- 3 edge cases (AGTTCGTATA..., GTATAAGGTG..., AGTTGTAATA...) lack TTAT entirely — left as-is

### RC Dedup
- 10 sequence pairs were reverse complements of each other
- Both orientations can appear valid CD (TTAT on right) when the spacer region is palindromic-ish
- Added RC dedup step after orientation to remove these duplicates

### CR Grouping
- Changed hd threshold from ≤2 to ≤1 after discovering transitive chaining across 6bp strings
- With hd≤2: 16 CRs merged into one group (chain: GGTGTA→CGTGTA→...→CAACCA)
- With hd≤1: 14 distinct groups, more biologically meaningful
- Group representatives: GGTGTA (10 members), CAGCCA (7), GAGATT (5), AATCTC (5), etc.

### CSV Parsing Issue
- literature_reference.csv has comment lines starting with `#` — csv.DictReader treats them as data
- Fixed by pre-filtering lines with `not ln.startswith("#")` before CSV parsing

### Output Files
- xer_seeds_v1.fa: 46 CD seeds → 46 DC seeds (verified RC relationship)
- xer_seeds_v1.pwm.json: 14 groups with per-position nucleotide frequencies
- build_report.json: Complete provenance tracking and summary stats
