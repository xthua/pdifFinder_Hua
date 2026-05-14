# Literature Benchmark for PdifFinder Accuracy

## TL;DR

> **Quick Summary**: Build a comprehensive framework to systematically compare PdifFinder's pdif site detection results against literature-reported ground truth, computing all four accuracy metrics (detection rate, precision, position accuracy, FPR). Two benchmark modes: pdif-only (mock AMR) and full-pipeline (real AMR genes).

> **Deliverables**:
> - Fixed `parse_pdif_output()` in `tests/helpers.py`
> - Curated reference dataset: `tests/benchmark_data/literature_reference.csv`
> - Benchmark harness: `tests/benchmark_literature.py`
> - Pytest-integrated verification: `tests/validation/test_literature_benchmark.py`
> - Benchmark results JSON output with per-plasmid failure categorization

> **Estimated Effort**: Medium
> **Parallel Execution**: YES - 4 waves
> **Critical Path**: Task 1 → Task 3 → Task 7 → Task 12 → Task 19 → Final Verification

---

## Context

### Original Request
"pdiffinder需要跟文献报道的pdif结果进行比对，测试结果准确性" — Compare PdifFinder's detection results against literature-reported pdif results to test accuracy.

### Interview Summary
**Key Discussions**:
- **Metrics**: All four — detection rate (recall/sensitivity), precision, position accuracy, false positive rate
- **Data collection**: Automated NCBI retrieval + manual literature curation, with human review
- **AMR dependency**: Test both ways — mock AMR detection AND real plasmids with AMR genes
- **Test strategy**: Tests-after (build framework first, then verify)
- **Scope**: Comprehensive pdif/XerC/XerD literature with at minimum 10 plasmids from 3+ independent papers

**Research Findings**:
- Codebase has ZERO literature citations or DOIs; pdifFinder paper is Shao et al. (2023) *Briefings in Bioinformatics* (DOI: 10.1093/bib/bbac521)
- `redundant.seed.fa` provides 12 NCBI accessions (KY984047.1, CP045108.1, KY984046.1, etc.) — best starting point
- Best literature with pdif coordinates: Cameranesi et al. (2018) Table S5 (17 sites with coordinates), Blackwell & Hall (2017) pS30-1 (8 pdif sites)
- `parse_pdif_output()` has column offset bug (expects 4 columns, actual output has 9)
- pdifFinder algorithm requires AMR genes before searching for pdif sites (guard at `singleThread()` line 797)
- Existing `test_against_real_pdif.py` is `@pytest.mark.xfail` — uses synthetic sequences that don't match seed database

### Metis Review
**Identified Gaps** (addressed):
- **parse_pdif_output bug**: Fixed as Phase 0 prerequisite (Task 1)
- **No metric thresholds**: Defined: DR≥0.80, precision≥0.70, position_median≤5bp, FPR≤0.01/kb
- **Circular coordinate wrapping**: Handled in benchmark comparison logic
- **Strand awareness**: Reverse-complement matches accepted; strand concordance reported separately
- **Unlisted sites handling**: Reported as "unvalidated detections" excluded from precision/recall
- **Failure categorization**: NO_SEED_MATCH, NO_PAIR, WRONG_POSITION, NOT_SOUGHT scheme
- **Data immutability**: Reference CSV committed to version control with version/date stamp
- **Separate entry point**: Not in default pytest suite; `@pytest.mark.benchmark` marker

---

## Work Objectives

### Core Objective
Build an automated benchmark framework that runs PdifFinder on real plasmid sequences with literature-documented pdif site positions, computes accuracy metrics against ground truth, and produces reproducible, categorized results.

### Concrete Deliverables
- **Fixed parser**: Corrected `parse_pdif_output()` in `tests/helpers.py` (9-column output)
- **Reference data**: `tests/benchmark_data/literature_reference.csv` — ≥10 plasmids, ≥3 papers, with pdif positions
- **Benchmark harness**: `tests/benchmark_literature.py` — standalone script + pytest-compatible
- **Plasmid cache**: `tests/benchmark_data/plasmids/` — FASTA files for reference plasmids (committed, not runtime-downloaded)
- **Pytest integration**: `tests/validation/test_literature_benchmark.py` with `@pytest.mark.benchmark`
- **Results output**: `tests/benchmark_data/results/` — JSON per run with full failure categorization

### Definition of Done
- [ ] `python -m tests.benchmark_literature --mode pdif_only` runs successfully on all cached plasmids
- [ ] `python -m tests.benchmark_literature --mode full_pipeline` runs on plasmids with verified AMR genes
- [ ] Detection rate, precision, position accuracy, FPR all computed and ≥ minimum thresholds
- [ ] Results JSON includes per-plasmid breakdown with failure categorization
- [ ] `pytest tests/validation/test_literature_benchmark.py -m benchmark` passes

### Must Have
- Dual benchmark modes: `--mode pdif_only` (mock AMR) and `--mode full_pipeline` (real blastn)
- Position comparison with 5bp tolerance and circular coordinate wrapping support
- Strand-agnostic matching (reverse-complement accepted)
- Per-plasmid failure categorization: NO_SEED_MATCH, NO_PAIR, WRONG_POSITION, NOT_SOUGHT
- Reference data committed as CSV, not downloaded at runtime
- Reproducible: pinned version, seeded random, environment documented

### Must NOT Have (Guardrails)
- DO NOT modify `PdifFinder/pdifFinder.py` source code (including `singleThread()`, `findPdif()`, etc.)
- DO NOT modify `PdifFinder/data/redundant.seed.fa` or `PdifFinder/data/pdifdatabase.fasta`
- DO NOT add literature benchmark tests to default pytest suite (use `@pytest.mark.benchmark` marker)
- DO NOT download plasmid sequences at runtime — must be cached/committed
- DO NOT build literature search web scraper — finite manual curation only
- DO NOT visualize results (no charts, HTML reports) — metrics as JSON/CSV only
- DO NOT compare to other bioinformatics tools — pdifFinder vs. literature ground truth only
- DO NOT include "user manually verifies results" as acceptance criterion

---

## Verification Strategy (MANDATORY)

> **ZERO HUMAN INTERVENTION** - ALL verification is agent-executed.

### Test Decision
- **Infrastructure exists**: YES (pytest, conftest, helpers)
- **Automated tests**: Tests-after (build framework first, then verification tests)
- **Framework**: pytest (existing)
- **QA Policy**: Every task includes agent-executed QA scenarios. Evidence saved to `.sisyphus/evidence/task-{N}-{scenario-slug}.{ext}`.

### QA Policy
Every task MUST include agent-executed QA scenarios:
- **Frontend/UI**: N/A (this is a CLI/bioinformatics tool)
- **TUI/CLI**: Use `interactive_bash` (tmux) - Run benchmark commands, validate output
- **API/Backend**: Use `bash` (pytest) - Run test suites, assert metrics, check JSON output
- **Data integrity**: Use `bash` (python) - Validate CSV format, verify NCBI accessions reachable
- **Library/Module**: Use `bash` (python REPL) - Import modules, call functions, compare output

---

## Execution Strategy

### Parallel Execution Waves

```
Wave 1 (Start Immediately - foundation + data):
├── Task 1: Fix parse_pdif_output() column bug [quick]
├── Task 2: Extract plasmid accessions from redundant.seed.fa [quick]
├── Task 3: Download and cache plasmid sequences from NCBI [quick]
├── Task 4: Extract Cameranesi 2018 Table S5 coordinates [quick]
├── Task 5: Extract Blackwell & Hall 2017 pS30-1 coordinates [quick]
└── Task 6: Create literature_reference.csv schema + template [quick]

Wave 2 (After Wave 1 - harness core, MAX PARALLEL):
├── Task 7: Build benchmark harness skeleton [deep]
├── Task 8: Implement pdif_only mode (mock AMR) [deep]
├── Task 9: Implement full_pipeline mode (real blastn) [deep]
├── Task 10: Implement position comparison with circular wrap [deep]
├── Task 11: Implement failure categorization logic [unspecified-high]
└── Task 12: Implement metrics computation (DR, precision, position, FPR) [unspecified-high]

Wave 3 (After Wave 2 - integration + verification):
├── Task 13: Run pdif_only benchmark on all cached plasmids [unspecified-high]
├── Task 14: Run full_pipeline benchmark on AMR-positive subset [unspecified-high]
├── Task 15: Write pytest integration tests [quick]
├── Task 16: Generate curated literature_reference.csv (populate with real data) [deep]
├── Task 17: Validate results JSON output format [quick]
└── Task 18: Document benchmark usage in README [writing]

Wave FINAL (After ALL tasks — 4 parallel reviews, then user okay):
├── Task F1: Plan Compliance Audit (oracle)
├── Task F2: Code Quality Review (unspecified-high)
├── Task F3: Real Manual QA (unspecified-high)
└── Task F4: Scope Fidelity Check (deep)
-> Present results -> Get explicit user okay

Critical Path: Task 1 → Task 7 → Task 12 → Task 13 → Final Verification
Parallel Speedup: ~55% faster than sequential
Max Concurrent: 6 (Waves 1 & 2)
```

### Agent Dispatch Summary

- **Wave 1**: **6** - T1-T6 → quick
- **Wave 2**: **6** - T7 → deep, T8 → deep, T9 → deep, T10 → deep, T11 → unspecified-high, T12 → unspecified-high
- **Wave 3**: **6** - T13 → unspecified-high, T14 → unspecified-high, T15 → quick, T16 → deep, T17 → quick, T18 → writing
- **FINAL**: **4** - F1 → oracle, F2 → unspecified-high, F3 → unspecified-high, F4 → deep

---

## TODOs

> Implementation + Test = ONE Task. Never separate.
> EVERY task MUST have: Recommended Agent Profile + QA Scenarios.
> **A task WITHOUT QA Scenarios is INCOMPLETE. No exceptions.**

- [ ] 1. Fix `parse_pdif_output()` column offset bug in `tests/helpers.py`

  **What to do**:
  - Read `tests/helpers.py:165-229` to understand current `parse_pdif_output()`
  - The actual `pdif_site.txt` output format (from `changepdifname()` in `pdifFinder.py`) has 9 columns: `[seq_id, pdifName, start, end, xerC(11bp), spacer(6bp), xerD(11bp), orientation(C|D or D|C), pdif_db_name]`
  - Current parser reads only 4 columns and misaligns: reads `pdifName` as `start`, `start` as `end`, `end` as `pdif_site`
  - Fix to correctly parse all 9 columns, returning full information in the result dict
  - Updated `pdif_sites` dict: `{'sequence_id', 'pdif_name', 'start', 'end', 'xerC', 'spacer', 'xerD', 'orientation', 'pdif_db_name'}`
  - Verify ALL existing callers still work after the fix (use `lsp_find_references` to enumerate)
  - Update callers if they accessed the misaligned column names

  **Must NOT do**:
  - Do NOT change the output format of `pdifFinder.py` — only fix the parser
  - Do NOT remove any existing fields from the result dict (add to it, don't delete)

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Single-file bug fix with clear specification. Straightforward column alignment problem.
  - **Skills**: [`lsp_find_references`]
    - `lsp_find_references`: Needed to enumerate all callers of `parse_pdif_output` before making changes

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 2-6)
  - **Blocks**: Tasks 7, 8, 9, 12 (all benchmark harness tasks depend on correct parsing)
  - **Blocked By**: None (can start immediately)

  **References**:
  - `tests/helpers.py:165-229` — current `parse_pdif_output()` implementation to fix
  - `PdifFinder/pdifFinder.py:929-951` — `changepdifname()` that writes `pdif_site.txt` (actual 9-column format)
  - `PdifFinder/pdifFinder.py:1118-1190` — `finalFilter()` that writes intermediate `pdif_site.txt` (also 9 columns)
  - Pattern: `tests/helpers.py:validate_pdif_positions()` at line 232 — position comparison function that calls into `parse_pdif_output`, must not break

  **Acceptance Criteria**:

  **QA Scenarios**:

  ```
  Scenario: Parse real 9-column pdif_site.txt output
    Tool: Bash (python)
    Preconditions: Run pdifFinder on a plasmid with known pdif sites, producing pdif_site.txt
    Steps:
      1. python -c "from tests.helpers import parse_pdif_output; r = parse_pdif_output('tmp_out'); print(r['pdif_sites'])"
      2. Assert all 9 fields present in each pdif_site dict: sequence_id, pdif_name, start, end, xerC, spacer, xerD, orientation, pdif_db_name
      3. Assert start < end (positions are valid)
      4. Assert len(xerC) == 11, len(spacer) == 6, len(xerD) == 11
      5. Assert orientation in ('C|D', 'D|C')
    Expected Result: All pdif_sites parsed with correct field mapping; no column offset
    Evidence: .sisyphus/evidence/task-1-parse-9col.json

  Scenario: Existing callers still work after fix
    Tool: Bash (pytest)
    Preconditions: Fixed parse_pdif_output() implemented
    Steps:
      1. grep -rn "parse_pdif_output" tests/ --include="*.py" to find all callers
      2. For each caller, verify it accesses the correct field names in the updated dict
      3. pytest tests/ -x -k "parse" to run any test that exercises parsing
    Expected Result: All existing tests pass without modification (or minimal field-name updates)
    Failure Indicators: KeyError on field access; test assertions fail
    Evidence: .sisyphus/evidence/task-1-callers-pass.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-1-parse-9col.json` — parsed output json
  - [ ] `.sisyphus/evidence/task-1-callers-pass.txt` — pytest output showing all caller tests pass

  **Commit**: YES (groups with Task 2)
  - Message: `fix(tests): correct parse_pdif_output column offset for 9-col pdif_site.txt`
  - Files: `tests/helpers.py`
  - Pre-commit: `pytest tests/ -x -k "parse or pipeline or validation"`

- [ ] 2. Extract and validate plasmid accessions from `redundant.seed.fa`

  **What to do**:
  - Read `PdifFinder/data/redundant.seed.fa` (82 lines)
  - Parse all FASTA headers to extract plasmid accessions (format: `>ACCESSION_POSITION`)
  - Deduplicate and categorize: NCBI accessions vs. non-standard names
  - Validate each NCBI accession by checking if it resolves via NCBI eutils API
  - Produce `tests/benchmark_data/plasmid_accessions.json` with structure:
    ```json
    {"accessions": [{"id": "KY984047.1", "valid": true, "pdif_count": 5}, ...]}
    ```
  - Also extract the 28bp pdif sequences for each accession (these are the XerC-spacer-XerD full sequences)
  - Document which accessions are NOT in NCBI (pAB1H8, TE_01Z_000_contig11, etc.)

  **Must NOT do**:
  - Do NOT modify the seed database
  - Do NOT assume all accessions are valid — validate each

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: File parsing + API validation. Straightforward data extraction task.
  - **Skills**: []
    - No specialized skills needed — standard Python file I/O + HTTP requests

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 3-6)
  - **Blocks**: Tasks 3, 16 (plasmid downloading and reference data population)
  - **Blocked By**: None

  **References**:
  - `PdifFinder/data/redundant.seed.fa` — source file with FASTA headers containing accessions
  - `PdifFinder/pdifFinder.py:663-684` — `findPdif()` where seedList is constructed from this file

  **Acceptance Criteria**:
  - [ ] `plasmid_accessions.json` exists with ≥12 entries
  - [ ] Each accession validated via NCBI eutils

  **QA Scenarios**:

  ```
  Scenario: All accessions extracted and categorized
    Tool: Bash (python)
    Preconditions: None
    Steps:
      1. python -c "import json; data=json.load(open('tests/benchmark_data/plasmid_accessions.json')); print(f'Total: {len(data[\"accessions\"])}'); valid=[a for a in data['accessions'] if a['valid']]; print(f'Valid NCBI: {len(valid)}')"
      2. Assert total accessions >= 12
      3. Assert >= 10 valid NCBI accessions
      4. Check each valid accession has pdif_count >= 1
    Expected Result: At least 12 accessions extracted, ≥10 valid in NCBI
    Evidence: .sisyphus/evidence/task-2-accessions.json

  Scenario: NCBI validation works for known accession
    Tool: Bash (curl)
    Preconditions: KY984047.1 is a known accession
    Steps:
      1. curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=KY984047.1&retmode=json" | python -c "import sys,json; print(json.load(sys.stdin)['esearchresult']['count'])"
      2. Assert count >= 1 (accession resolves)
    Expected Result: NCBI returns at least 1 result for KY984047.1
    Evidence: .sisyphus/evidence/task-2-ncbi-validate.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-2-accessions.json`
  - [ ] `.sisyphus/evidence/task-2-ncbi-validate.txt`

  **Commit**: YES (groups with Task 1)
  - Message: `feat(benchmark): extract and validate plasmid accessions from seed database`
  - Files: `tests/benchmark_data/plasmid_accessions.json`

- [ ] 3. Download and cache plasmid sequences from NCBI

  **What to do**:
  - Read `tests/benchmark_data/plasmid_accessions.json` from Task 2
  - For each valid NCBI accession, download the plasmid sequence using NCBI eutils (efetch)
  - Save to `tests/benchmark_data/plasmids/{ACCESSION}.fasta` (FASTA format)
  - Also create a combined multi-FASTA `tests/benchmark_data/plasmids/all_plasmids.fasta`
  - For each plasmid, also attempt to download GenBank format (`{ACCESSION}.gb`) if available (needed for full_pipeline mode's AMR detection)
  - Log any failures (accession removed, restricted, etc.) to `tests/benchmark_data/plasmids/download_errors.json`
  - Skip accessions that fail validation in Task 2

  **Must NOT do**:
  - Do NOT download at runtime — sequences must be cached/committed
  - Do NOT use accession pAB1H8 or TE_* non-standard identifiers (can't download from NCBI)

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Automated download + file I/O. Batch HTTP requests with error handling.
  - **Skills**: []
    - Standard Python with urllib/requests for NCBI API

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1-2, 4-6)
  - **Blocks**: Tasks 7, 8, 9, 13, 14 (benchmark requires cached sequences)
  - **Blocked By**: Task 2 (needs accessions list)

  **References**:
  - NCBI eutils efetch: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={ACCESSION}&rettype=fasta&retmode=text`
  - Pattern: `tests/helpers.py:save_sequence_to_fasta()` — follow this pattern for saving FASTA files
  - `tests/benchmark_data/` — target directory for cached plasmids

  **Acceptance Criteria**:
  - [ ] ≥10 plasmid FASTA files cached in `tests/benchmark_data/plasmids/`
  - [ ] `all_plasmids.fasta` combined file exists
  - [ ] `download_errors.json` documents any failures

  **QA Scenarios**:

  ```
  Scenario: KY984047.1 downloaded and is valid FASTA
    Tool: Bash (python + BioPython)
    Preconditions: Task 2 completed, KY984047.1 validated
    Steps:
      1. python -c "from Bio import SeqIO; rec = next(SeqIO.parse('tests/benchmark_data/plasmids/KY984047.1.fasta', 'fasta')); print(f'ID: {rec.id}, Length: {len(rec.seq)}')"
      2. Assert record.id contains 'KY984047'
      3. Assert len(rec.seq) > 1000 (plasmids are >1kb)
      4. Assert all(base in 'ACGTNacgtn' for base in str(rec.seq)[:100]) (valid DNA)
    Expected Result: Valid FASTA with plasmid sequence >1kb
    Evidence: .sisyphus/evidence/task-3-fasta-validate.txt

  Scenario: Combined multi-FASTA contains all accessions
    Tool: Bash (python)
    Preconditions: All accessions downloaded
    Steps:
      1. python -c "from Bio import SeqIO; records=list(SeqIO.parse('tests/benchmark_data/plasmids/all_plasmids.fasta','fasta')); print(f'Records: {len(records)}')"
      2. Assert len(records) >= 10
      3. Assert all(len(r.seq) > 0 for r in records)
    Expected Result: Multi-FASTA with ≥10 records, all non-empty
    Evidence: .sisyphus/evidence/task-3-multifasta.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-3-fasta-validate.txt`
  - [ ] `.sisyphus/evidence/task-3-multifasta.txt`

  **Commit**: YES
  - Message: `feat(benchmark): cache plasmid sequences from NCBI for literature benchmark`
  - Files: `tests/benchmark_data/plasmids/`, `tests/benchmark_data/plasmids/download_errors.json`

- [ ] 4. Extract pdif coordinates from Cameranesi et al. (2018) Table S5

  **What to do**:
  - Locate and fetch Cameranesi et al. (2018) *Frontiers in Microbiology* paper (DOI: 10.3389/fmicb.2018.00066)
  - Download Table S5 or supplementary material
  - Extract all pdif site entries with: plasmid name, pdif_start, pdif_end, pdif_sequence, strand
  - Cross-reference plasmid names with NCBI accessions (e.g., pAb242_25 → find NCBI accession)
  - Create `tests/benchmark_data/literature_raw/cameranesi2018.csv` with extracted data
  - Document extraction methodology (which supplementary file, how coordinates were mapped)

  **Must NOT do**:
  - Do NOT assume Table S5 coordinates are 1-based — verify against the paper's description
  - Do NOT guess NCBI accessions — search for them properly

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Literature data extraction requires careful reading, cross-referencing, and verification. Multiple sources to check.
  - **Skills**: []
    - We may use `webfetch` to retrieve the paper's supplementary data

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1-3, 5-6)
  - **Blocks**: Task 16 (populating final reference CSV)
  - **Blocked By**: None

  **References**:
  - Paper: Cameranesi et al. (2018) *Front Microbiol* 9:66, DOI: 10.3389/fmicb.2018.00066
  - URL: `https://www.frontiersin.org/articles/10.3389/fmicb.2018.00066/full`
  - Supplementary material at: `https://www.frontiersin.org/articles/10.3389/fmicb.2018.00066/full#supplementary-material`

  **Acceptance Criteria**:
  - [ ] `cameranesi2018.csv` contains ≥10 pdif site entries from ≥2 plasmids
  - [ ] Each entry has plasmid name, pdif sequence, and coordinate info

  **QA Scenarios**:

  ```
  Scenario: Cameranesi 2018 data extracted with valid coordinates
    Tool: Bash (python)
    Preconditions: cameranesi2018.csv created
    Steps:
      1. python -c "import csv; rows=list(csv.DictReader(open('tests/benchmark_data/literature_raw/cameranesi2018.csv'))); print(f'Entries: {len(rows)}')"
      2. Assert len(rows) >= 10
      3. For each row, assert pdif_start < pdif_end, len(pdif_sequence) == 28
      4. Print summary: unique plasmids, total pdif sites
    Expected Result: ≥10 pdif site entries with valid coordinates
    Evidence: .sisyphus/evidence/task-4-cameranesi.csv
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-4-cameranesi.csv`

  **Commit**: YES (groups with Task 5)
  - Message: `feat(benchmark): extract pdif coordinates from Cameranesi 2018 Table S5`
  - Files: `tests/benchmark_data/literature_raw/cameranesi2018.csv`

- [ ] 5. Extract pdif coordinates from Blackwell & Hall (2017) pS30-1

  **What to do**:
  - Fetch Blackwell & Hall (2017) *Antimicrobial Agents and Chemotherapy* paper (DOI: 10.1128/aac.00780-17)
  - Extract the 8 pdif sites on plasmid pS30-1 (GenBank accession KU987654)
  - Map each pdif site to its coordinates on the plasmid (from paper text, tables, or GenBank annotation)
  - Download KU987654 from NCBI if not already cached
  - Extract pdif site sequences from the plasmid at the reported coordinates to verify
  - Create `tests/benchmark_data/literature_raw/blackwell2017.csv`

  **Must NOT do**:
  - Do NOT assume all 8 pdif sites are explicitly listed — verify GenBank annotations match paper text

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Requires careful reading of paper text AND GenBank annotation cross-referencing
  - **Skills**: []
    - May use `webfetch` to retrieve paper and GenBank entry

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1-4, 6)
  - **Blocks**: Task 16 (populating final reference CSV)
  - **Blocked By**: None

  **References**:
  - Paper: Blackwell & Hall (2017) *Antimicrob Agents Chemother* 61:e00780-17, DOI: 10.1128/aac.00780-17
  - GenBank: KU987654 (pS30-1)

  **Acceptance Criteria**:
  - [ ] `blackwell2017.csv` contains 8 pdif site entries from pS30-1
  - [ ] Coordinates verified against GenBank annotation

  **QA Scenarios**:

  ```
  Scenario: pS30-1 pdif sites extracted and verified
    Tool: Bash (python)
    Preconditions: blackwell2017.csv created, KU987654 downloaded
    Steps:
      1. python -c "import csv; rows=list(csv.DictReader(open('tests/benchmark_data/literature_raw/blackwell2017.csv'))); print(f'Entries: {len(rows)}')"
      2. Assert len(rows) == 8
      3. For each row's coordinates, extract the 28bp from KU987654 sequence at that position
      4. Verify extracted sequence matches the pdif_sequence in the CSV (within 2bp tolerance for variants)
    Expected Result: 8 pdif sites with coordinates that match plasmid sequence
    Evidence: .sisyphus/evidence/task-5-blackwell.csv
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-5-blackwell.csv`

  **Commit**: YES (groups with Task 4)
  - Message: `feat(benchmark): extract pdif coordinates from Blackwell & Hall 2017 pS30-1`
  - Files: `tests/benchmark_data/literature_raw/blackwell2017.csv`

- [ ] 6. Create `literature_reference.csv` schema and template

  **What to do**:
  - Define the canonical CSV schema for `tests/benchmark_data/literature_reference.csv`:
    ```
    plasmid_accession, pdif_start, pdif_end, pdif_sequence, strand, orientation, source_pmid, source_description, curator_notes, curation_date
    ```
  - Create the CSV template with header row and 1-2 example rows
  - Write a Python validation script: `tests/benchmark_data/validate_reference.py` that:
    - Checks all required columns exist
    - Validates coordinate ranges (start < end, within plasmid length)
    - Validates pdif_sequence is 28bp DNA
    - Validates source_pmid is a valid integer
    - Checks for duplicate entries
  - Write schema documentation in a header comment block at the top of the CSV file
  - Add a `version` field and `curation_date` to track data provenance

  **Must NOT do**:
  - Do NOT populate with real data yet — that's Task 16
  - Do NOT create a complex database — one CSV file is sufficient

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Schema design and validation script. Standard data engineering task.
  - **Skills**: []
    - Standard Python for CSV validation

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1-5)
  - **Blocks**: Task 16 (populates the CSV)
  - **Blocked By**: None

  **References**:
  - PDIF site structure: `PdifFinder/pdifFinder.py:662-683` — XerC(0-11), spacer(11-17), XerD(17-28)
  - Pattern: `tests/test_data/verify_sequences.py` — similar validation script pattern

  **Acceptance Criteria**:
  - [ ] `literature_reference.csv` exists with header and schema comments
  - [ ] `validate_reference.py` runs successfully on the template
  - [ ] Validation catches: missing columns, invalid coordinates, non-28bp sequences, non-numeric PMID

  **QA Scenarios**:

  ```
  Scenario: Template CSV passes validation
    Tool: Bash (python)
    Preconditions: Template created
    Steps:
      1. python tests/benchmark_data/validate_reference.py
      2. Assert exit code 0
    Expected Result: Validation passes on template CSV
    Evidence: .sisyphus/evidence/task-6-template-validate.txt

  Scenario: Validation catches invalid data
    Tool: Bash (python)
    Preconditions: Template CSV exists
    Steps:
      1. Create a copy with an invalid row (start > end)
      2. python tests/benchmark_data/validate_reference.py --file /tmp/bad.csv; echo "Exit: $?"
      3. Assert exit code != 0
      4. Assert error message mentions the specific issue
    Expected Result: Validation catches and reports invalid data
    Evidence: .sisyphus/evidence/task-6-validation-error.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-6-template-validate.txt`
  - [ ] `.sisyphus/evidence/task-6-validation-error.txt`

  **Commit**: YES
  - Message: `feat(benchmark): create literature_reference.csv schema and validation script`
  - Files: `tests/benchmark_data/literature_reference.csv`, `tests/benchmark_data/validate_reference.py`

- [ ] 7. Build benchmark harness skeleton

  **What to do**:
  - Create `tests/benchmark_literature.py` as the main benchmark entry point
  - Implement CLI argument parsing: `--mode {pdif_only|full_pipeline}`, `--output-dir`, `--plasmid`, `--verbose`
  - Skeleton structure:
    ```python
    class LiteratureBenchmark:
        def __init__(self, mode, output_dir, plasmid_filter=None): ...
        def load_reference_data(self) -> List[ReferenceEntry]: ...
        def load_plasmid_sequences(self) -> Dict[str, str]: ...
        def run_pdif_finder(self, fasta_path, output_dir) -> PdifOutput: ...
        def compare_results(self, reference, detected) -> ComparisonResult: ...
        def categorize_failure(self, reference_entry, detected_list) -> FailureCategory: ...
        def compute_metrics(self, all_comparisons) -> Metrics: ...
        def run(self) -> BenchmarkResult: ...
    ```
  - Wire up `if __name__ == "__main__"` entry point
  - Implement logging infrastructure (log each step, save intermediate results)
  - Create `tests/benchmark_data/results/` directory for output

  **Must NOT do**:
  - Do NOT implement actual comparison logic yet — that's Tasks 8-12
  - Do NOT invoke pdifFinder via subprocess — use direct import approach

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Core architecture design. Must correctly structure the class hierarchy, CLI, and data flow.
  - **Skills**: []
    - Standard Python — argparse, dataclasses, logging

  **Parallelization**:
  - **Can Run In Parallel**: NO (sequential)
  - **Parallel Group**: Wave 2
  - **Blocks**: Tasks 8-12 (all depend on the skeleton), Tasks 13-14 (depend on full harness)
  - **Blocked By**: Task 1 (needs fixed parser), Tasks 2-3 (needs plasmid data for testing)

  **References**:
  - `tests/benchmark.py` — existing benchmark script pattern (argparse, class structure)
  - `tests/helpers.py:126-164` — `run_pdif_finder()` for understanding pdifFinder invocation
  - `tests/helpers.py:165-229` — `parse_pdif_output()` for result parsing (post-Task 1 fix)
  - `tests/helpers.py:232-271` — `validate_pdif_positions()` for position comparison approach
  - `tests/validate_output.py:1-28` — pattern for standalone validation script with argparse

  **Acceptance Criteria**:
  - [ ] `tests/benchmark_literature.py` runs `--help` successfully
  - [ ] `--mode pdif_only --plasmid KY984047.1 --output-dir /tmp/test` executes without crash (even if no comparison logic yet)
  - [ ] Logging outputs to both console and file

  **QA Scenarios**:

  ```
  Scenario: Benchmark skeleton accepts CLI arguments and initializes
    Tool: Bash (python)
    Preconditions: Skeleton implemented
    Steps:
      1. python -m tests.benchmark_literature --help
      2. Assert exit code 0
      3. Assert output contains --mode, --output-dir, --plasmid
      4. python -m tests.benchmark_literature --mode pdif_only --plasmid KY984047.1 --output-dir /tmp/test_bench --verbose 2>&1 | head -20
      5. Assert output contains "Loading reference data" or similar log message
    Expected Result: CLI accepts arguments, skeleton runs without error
    Failure Indicators: ImportError, argparse error, crash before logging
    Evidence: .sisyphus/evidence/task-7-skeleton-run.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-7-skeleton-run.txt`

  **Commit**: YES
  - Message: `feat(benchmark): add benchmark harness skeleton with CLI and class structure`
  - Files: `tests/benchmark_literature.py`

- [ ] 8. Implement pdif_only mode (mock AMR detection)

  **What to do**:
  - Implement `pdif_only` mode in the benchmark harness
  - This mode MUST bypass the AMR gene requirement by mocking `findResistanceGene`
  - Approach: Use `unittest.mock.patch` on `PdifFinder.pdifFinder.findResistanceGene` to return a dummy position list `["1-1000"]` that ensures the guard at `singleThread()` line 797 passes
  - Call `singleThread()` directly (not via subprocess) to invoke the full pdif detection pipeline with the mock
  - Handle all intermediate file cleanup (the pipeline creates `tmp/` dirs)
  - Parse results using the fixed `parse_pdif_output()` from Task 1
  - Log which plasmids were processed in pdif_only mode
  - Ensure the mock is applied ONLY for pdif_only mode, NOT for full_pipeline mode

  **Must NOT do**:
  - Do NOT modify `singleThread()` or `findResistanceGene()` in the source — use mock in test code only
  - Do NOT leave mock applied after the benchmark run (ensure proper cleanup)

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Requires understanding the pdifFinder pipeline internals (singleThread, findPdif, findMatchFragmentThread) and correct mock placement. Risk of side effects.
  - **Skills**: []
    - Standard Python with `unittest.mock` for patching

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Tasks 9-12, subject to Task 7 completion)
  - **Blocks**: Tasks 13, 15 (benchmark run and pytest tests)
  - **Blocked By**: Task 7 (skeleton)

  **References**:
  - `PdifFinder/pdifFinder.py:795-799` — `singleThread()` where the AMR guard is: `if resistanceGenePosList != []: flag = findPdif(...)`
  - `PdifFinder/pdifFinder.py:1009-1073` — `findResistanceGene()` to understand return format
  - `tests/validation/test_against_real_pdif.py:224-226` — existing `@patch('PdifFinder.pdifFinder.findResistanceGene')` pattern
  - `tests/helpers.py:126-164` — `run_pdif_finder()` current subprocess approach (we're replacing this with direct call)

  **Acceptance Criteria**:
  - [ ] `python -m tests.benchmark_literature --mode pdif_only --plasmid KY984047.1` runs successfully
  - [ ] PdifFinder actually searches for pdif sites (mock ensures AMR guard passes)
  - [ ] Results parsed using fixed `parse_pdif_output()`

  **QA Scenarios**:

  ```
  Scenario: pdif_only mode detects pdif sites on KY984047.1
    Tool: Bash (python)
    Preconditions: Tasks 1, 2, 3, 7 completed; KY984047.1 downloaded
    Steps:
      1. python -m tests.benchmark_literature --mode pdif_only --plasmid KY984047.1 --output-dir /tmp/test_pdif_only --verbose 2>&1 | tee /tmp/pdif_only.log
      2. Assert exit code 0
      3. Assert log contains "pdif_only mode" 
      4. Assert log contains "Resistance gene detection MOCKED" or similar
      5. Check output files exist in /tmp/test_pdif_only/
    Expected Result: Benchmark runs, logs confirm mock AMR mode, pdif detection pipeline executed
    Failure Indicators: "No resistance genes found" message, crash, no output files
    Evidence: .sisyphus/evidence/task-8-pdif-only.log

  Scenario: pdif_only mock is properly cleaned up after run
    Tool: Bash (python)
    Preconditions: After running pdif_only mode
    Steps:
      1. python -c "import PdifFinder.pdifFinder as pf; print(pf.findResistanceGene.__name__)" 
      2. Assert the function is NOT mocked (returns its real name, not 'MagicMock')
    Expected Result: findResistanceGene is restored to original after benchmark
    Evidence: .sisyphus/evidence/task-8-mock-cleanup.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-8-pdif-only.log`
  - [ ] `.sisyphus/evidence/task-8-mock-cleanup.txt`

  **Commit**: YES (groups with Task 9)
  - Message: `feat(benchmark): implement pdif_only mode with mocked AMR detection`
  - Files: `tests/benchmark_literature.py`

- [ ] 9. Implement full_pipeline mode (real blastn AMR detection)

  **What to do**:
  - Implement `full_pipeline` mode in the benchmark harness
  - This mode uses REAL blastn to detect AMR genes on input plasmids
  - Check blastn availability at startup: if not found, skip full_pipeline mode with a clear message
  - Pass GenBank files (`.gb`) to pdifFinder for full_pipeline (GenBank format preserves feature annotations needed by AMR pipeline)
  - Use the AMR database at `PdifFinder/AMRDB/sequences` for blastn search
  - For plasmids without GenBank format available, generate a minimal GenBank from FASTA (or skip those plasmids with a warning)
  - Parse results using fixed `parse_pdif_output()`
  - Compare results to pdif_only mode results for the same plasmid — log any differences

  **Must NOT do**:
  - Do NOT mock any functions in full_pipeline mode
  - Do NOT crash if blastn missing — log warning and skip

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Requires understanding the full pdifFinder pipeline including blastn integration, GenBank parsing, and AMR database paths. Multiple failure modes to handle.
  - **Skills**: []
    - Standard Python; blastn is an external dependency

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Tasks 8, 10-12, subject to Task 7 completion)
  - **Blocks**: Tasks 14, 15 (full_pipeline benchmark run and pytest tests)
  - **Blocked By**: Task 7 (skeleton)

  **References**:
  - `PdifFinder/pdifFinder.py:1009-1073` — `findResistanceGene()` with blastn command construction
  - `PdifFinder/pdifFinder.py:795-799` — `singleThread()` full flow
  - `PdifFinder/pdifFinder.py:1077-1100` — `getSeqFromGenbankFile()` for GenBank parsing
  - `PdifFinder/pdifFinder.py:559-563` — `check_dependencies()` for blastn check pattern
  - `PdifFinder/AMRDB/` — check if `sequences` BLAST database exists

  **Acceptance Criteria**:
  - [ ] `python -m tests.benchmark_literature --mode full_pipeline --plasmid KY984047.1` runs (if blastn available) or skips gracefully
  - [ ] Full pipeline results parsed correctly
  - [ ] Missing blastn handled gracefully (skip, don't crash)

  **QA Scenarios**:

  ```
  Scenario: full_pipeline mode runs successfully with real blastn
    Tool: Bash (python)
    Preconditions: blastn in PATH, AMRDB exists, plasmid has GenBank format
    Steps:
      1. which blastn && echo "blastn found" || echo "blastn not found"
      2. If blastn found: python -m tests.benchmark_literature --mode full_pipeline --plasmid KY984047.1 --output-dir /tmp/test_full --verbose 2>&1 | tee /tmp/full_pipeline.log
      3. Assert exit code 0
      4. Assert log contains "full_pipeline mode"
      5. Assert log does NOT contain "MOCKED"
    Expected Result: Full pipeline runs with real blastn AMR detection
    Evidence: .sisyphus/evidence/task-9-full-pipeline.log

  Scenario: full_pipeline mode handles missing blastn gracefully
    Tool: Bash (python)
    Preconditions: blastn NOT in PATH (or mock the check)
    Steps:
      1. python -m tests.benchmark_literature --mode full_pipeline --plasmid KY984047.1 --output-dir /tmp/test_full 2>&1
      2. Assert exit code 0 (not crash)
      3. Assert output contains "blastn not found" or "skipping full_pipeline"
    Expected Result: Graceful skip with clear message, no crash
    Evidence: .sisyphus/evidence/task-9-no-blastn.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-9-full-pipeline.log`
  - [ ] `.sisyphus/evidence/task-9-no-blastn.txt`

  **Commit**: YES (groups with Task 8)
  - Message: `feat(benchmark): implement full_pipeline mode with real blastn AMR detection`
  - Files: `tests/benchmark_literature.py`

- [ ] 10. Implement position comparison with circular coordinate wrapping

  **What to do**:
  - Implement `compare_positions()` method: match detected pdif sites to reference sites within 5bp tolerance
  - Support circular plasmid coordinate wrapping: if a plasmid is circular, positions near the origin should wrap around
  - Handle strand-agnostic matching: accept reverse-complement matches as valid
  - For each reference pdif site, find the best matching detected pdif site
  - For each detected pdif site, determine if it matches any reference site (true positive) or is unvalidated
  - Track which reference sites were matched, which were missed, and which detected sites are unvalidated
  - Return structured `ComparisonResult` with: matched_pairs, missed_references, unvalidated_detections

  **Must NOT do**:
  - Do NOT change the tolerance from 5bp without explicit justification
  - Do NOT match across different plasmids

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Circular coordinate math is tricky (off-by-one errors common). Strand handling requires reverse-complement computation via BioPython. Multiple edge cases.
  - **Skills**: []
    - Standard Python + BioPython for reverse complement

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Tasks 8-9, 11-12, subject to Task 7 completion)
  - **Blocks**: Tasks 12, 13, 14 (metrics computation and benchmark runs)
  - **Blocked By**: Task 7 (skeleton)

  **References**:
  - `tests/helpers.py:232-271` — `validate_pdif_positions()` existing position comparison with 5bp tolerance
  - `PdifFinder/pdifFinder.py:995-998` — `reverComplement()` for reverse-complement computation
  - `PdifFinder/pdifFinder.py:272-273` — circular coordinate handling in `findPdif()`: `seqC = seq[pos:pos + 11] if pos + 11 <= seq_len else seq[pos:] + seq[:(pos + 11) % seq_len]`
  - `PdifFinder/pdifFinder.py:447-448` — circular mode in `findMatchFragmentThread()`: `limit = seq_len if is_circular else seq_len - 27`

  **Acceptance Criteria**:
  - [ ] Position comparison matches sites within 5bp tolerance
  - [ ] Circular wrap-around positions correctly compared
  - [ ] Reverse-complement matches accepted
  - [ ] ComparisonResult correctly categorizes matched/missed/unvalidated

  **QA Scenarios**:

  ```
  Scenario: Linear position comparison works within 5bp tolerance
    Tool: Bash (python)
    Preconditions: Comparison function implemented
    Steps:
      1. python -c "
  from tests.benchmark_literature import LiteratureBenchmark
  b = LiteratureBenchmark('pdif_only', '/tmp')
  ref = [{'start': 100, 'end': 128, 'pdif_sequence': 'ATTTAACATAAGGGCTGTTATACGAAAT'}]
  det = [{'start': 103, 'end': 130, 'xerC': 'ATTTAACATA', 'xerD': 'ATACGAAAT', 'orientation': 'C|D'}]
  result = b.compare_positions(ref, det, plasmid_length=1000, is_circular=False)
  print(f'Matched: {len(result.matched_pairs)}, Missed: {result.missed_references}, Unvalidated: {result.unvalidated_detections}')"
      2. Assert matched_pairs == 1
      3. Assert missed_references == 0
    Expected Result: Site matched within 5bp tolerance
    Evidence: .sisyphus/evidence/task-10-linear-match.txt

  Scenario: Circular wrap-around comparison works
    Tool: Bash (python)
    Preconditions: Circular comparison implemented
    Steps:
      1. python -c "
  from tests.benchmark_literature import LiteratureBenchmark
  b = LiteratureBenchmark('pdif_only', '/tmp')
  ref = [{'start': 1, 'end': 28, 'pdif_sequence': 'ATTTAACATAAGGGCTGTTATACGAAAT'}]
  det = [{'start': 99990, 'end': 17, 'xerC': 'ATTTAACATA', 'xerD': 'ATACGAAAT', 'orientation': 'C|D'}]
  result = b.compare_positions(ref, det, plasmid_length=100000, is_circular=True)
  print(f'Matched: {len(result.matched_pairs)}, Missed: {result.missed_references}')"
      2. Assert matched_pairs == 1 (wrap-around match)
    Expected Result: Wrap-around positions correctly matched
    Evidence: .sisyphus/evidence/task-10-circular-match.txt

  Scenario: Reverse-complement match accepted
    Tool: Bash (python)
    Preconditions: Strand-agnostic matching implemented
    Steps:
      1. python -c "
  from tests.benchmark_literature import LiteratureBenchmark
  from Bio.Seq import Seq
  b = LiteratureBenchmark('pdif_only', '/tmp')
  pdif_fwd = 'ATTTAACATAAGGGCTGTTATACGAAAT'
  pdif_rev = str(Seq(pdif_fwd).reverse_complement())
  ref = [{'start': 100, 'end': 128, 'pdif_sequence': pdif_fwd}]
  det = [{'start': 100, 'end': 128, 'xerC': pdif_rev[0:11], 'xerD': pdif_rev[17:28], 'orientation': 'D|C'}]
  result = b.compare_positions(ref, det, plasmid_length=1000, is_circular=False, strand_agnostic=True)
  print(f'Matched: {len(result.matched_pairs)}')"
      2. Assert matched_pairs == 1
    Expected Result: Reverse-complement match accepted as correct
    Evidence: .sisyphus/evidence/task-10-strand-match.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-10-linear-match.txt`
  - [ ] `.sisyphus/evidence/task-10-circular-match.txt`
  - [ ] `.sisyphus/evidence/task-10-strand-match.txt`

  **Commit**: YES (groups with Tasks 11, 12)
  - Message: `feat(benchmark): implement position comparison with circular wrap and strand support`
  - Files: `tests/benchmark_literature.py`

- [ ] 11. Implement failure categorization logic

  **What to do**:
  - Implement `categorize_failure()` method for each missed reference pdif site
  - Failure categories (from Metis review):
    - `NO_SEED_MATCH`: pdif site sequence doesn't match any seed in `redundant.seed.fa` within mismatch thresholds (maxMismatchXerC=3, maxMismatchXerD=2)
    - `NO_PAIR`: seed matched but failed the pairing step (`findPossiblePair`) — no compatible pdif pair found
    - `WRONG_POSITION`: pdifFinder detected a site but outside 5bp tolerance of expected position
    - `NOT_SOUGHT`: pdifFinder didn't even search (full_pipeline mode, no AMR genes detected on this plasmid)
  - For `NO_SEED_MATCH`, check if the reference pdif's XerC and XerD halves have matches in `redundant.seed.fa`
  - For `NO_PAIR`, check if individual seed matches exist but `findPossiblePair` rejected them
  - Log the failure category for each missed site with details

  **Must NOT do**:
  - Do NOT count `unvalidated detections` (pdifFinder found sites not in literature) as failures — they go in a separate report

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
    - Reason: Requires understanding pdifFinder's two-stage detection (seed match → pair filtering). Moderately complex logic with multiple branches.
  - **Skills**: []
    - Standard Python

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Tasks 8-10, 12, subject to Task 7 completion)
  - **Blocks**: Tasks 13, 14 (benchmark runs need failure categorization in output)
  - **Blocked By**: Task 7 (skeleton), Task 10 (position comparison for WRONG_POSITION category)

  **References**:
  - `PdifFinder/pdifFinder.py:847-888` — `findMatchFragmentThread()` seed matching with mismatch thresholds
  - `PdifFinder/pdifFinder.py:915-962` — `findPossiblePair()` and `findPossiblePairThread()` pairing logic
  - `PdifFinder/pdifFinder.py:995-998` — `reverComplement()` for checking reverse-complement seed matches
  - `PdifFinder/pdifFinder.py:1001-1006` — `compareTwoSeq()` for mismatch counting

  **Acceptance Criteria**:
  - [ ] All four failure categories implemented and detectable
  - [ ] Each missed site gets exactly one failure category
  - [ ] Category determination is deterministic and reproducible

  **QA Scenarios**:

  ```
  Scenario: Failure categorization correctly identifies NO_SEED_MATCH
    Tool: Bash (python)
    Preconditions: Failure categorization implemented
    Steps:
      1. python -c "
  from tests.benchmark_literature import LiteratureBenchmark
  b = LiteratureBenchmark('pdif_only', '/tmp')
  # pdif sequence that doesn't match any seed
  ref_entry = {'pdif_sequence': 'AAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'plasmid_accession': 'TEST'}
  result = b.categorize_failure(ref_entry, [], [])
  print(f'Category: {result}')"
      2. Assert result == 'NO_SEED_MATCH'
    Expected Result: Correctly identified as no seed match
    Evidence: .sisyphus/evidence/task-11-no-seed.txt

  Scenario: Failure categorization for NOT_SOUGHT in full_pipeline
    Tool: Bash (python)
    Preconditions: Full categorization implemented
    Steps:
      1. python -c "
  from tests.benchmark_literature import LiteratureBenchmark
  b = LiteratureBenchmark('full_pipeline', '/tmp')
  ref_entry = {'pdif_sequence': 'ATTTAACATAAGGGCTGTTATACGAAAT', 'plasmid_accession': 'TEST'}
  # No AMR genes found → pdifFinder never searched
  result = b.categorize_failure(ref_entry, [], [], amr_found=False)
  print(f'Category: {result}')"
      2. Assert result == 'NOT_SOUGHT'
    Expected Result: Correctly identified as not sought due to missing AMR
    Evidence: .sisyphus/evidence/task-11-not-sought.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-11-no-seed.txt`
  - [ ] `.sisyphus/evidence/task-11-not-sought.txt`

  **Commit**: YES (groups with Tasks 10, 12)
  - Message: `feat(benchmark): implement failure categorization for missed pdif sites`
  - Files: `tests/benchmark_literature.py`

- [ ] 12. Implement metrics computation (DR, precision, position, FPR)

  **What to do**:
  - Implement `compute_metrics()` method that calculates all four metrics from comparison results
  - **Detection Rate (Recall/Sensitivity)**:
    - `DR = matched_pairs / total_reference_sites`
    - Report both overall DR and per-plasmid DR
    - Minimum acceptable: 0.80
  - **Precision**:
    - `Precision = matched_pairs / (matched_pairs + false_positives)`
    - Note: "unvalidated detections" are NOT counted as false positives (they're excluded from precision calculation)
    - Only count as false positives: detected sites that clearly don't match any reference site AND are on negative-control plasmids
    - Minimum acceptable: 0.70
    - Metis note: report recall and precision together (they form a trade-off)
  - **Position Accuracy**:
    - For each matched pair: `position_error = max(|det_start - ref_start|, |det_end - ref_end|)`
    - Report: median position error, mean position error, max position error
    - Compute ONLY on correctly detected sites (conditioned on detection)
    - Maximum acceptable median: 5bp
  - **False Positive Rate (FPR)**:
    - `FPR = false_positives / total_plasmid_length_kb`
    - Only computed on negative-control plasmids (plasmids known NOT to contain pdif sites)
    - Maximum acceptable: 0.01 per kb
  - Output metrics as structured JSON
  - Include: overall metrics, per-plasmid breakdown, per-failure-category counts

  **Must NOT do**:
  - Do NOT count unvalidated detections as false positives
  - Do NOT report precision without also reporting recall
  - Do NOT compute position accuracy on undetected sites

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
    - Reason: Statistics computation with multiple edge cases. Must correctly handle the recall/precision trade-off and per-plasmid breakdowns.
  - **Skills**: []
    - Standard Python with json for output

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Tasks 8-11, subject to Task 7 completion)
  - **Blocks**: Tasks 13, 14, 15 (benchmark runs need metrics)
  - **Blocked By**: Task 7 (skeleton), Task 10 (comparison results needed for computation)

  **References**:
  - `tests/validation/test_against_real_pdif.py:215-221` — detection_rate and false_positive_rate formulas in existing test
  - `tests/helpers.py:232-271` — `validate_pdif_positions()` for position tolerance approach

  **Acceptance Criteria**:
  - [ ] All four metrics computed and output as JSON
  - [ ] Per-plasmid breakdown included
  - [ ] Per-failure-category counts included
  - [ ] Threshold checks produce clear PASS/FAIL per metric

  **QA Scenarios**:

  ```
  Scenario: Detection rate computed correctly with known data
    Tool: Bash (python)
    Preconditions: Metrics computation implemented
    Steps:
      1. python -c "
  from tests.benchmark_literature import LiteratureBenchmark
  b = LiteratureBenchmark('pdif_only', '/tmp')
  # 10 reference, 8 matched, 2 missed
  comparison = type('obj', (object,), {'matched_pairs': list(range(8)), 'missed_references': list(range(2)), 'unvalidated_detections': [], 'false_positives': []})()
  metrics = b.compute_metrics(comparison, total_plasmid_kb=50)
  print(f'DR: {metrics[\"detection_rate\"]:.2f}, Precision: {metrics[\"precision\"]:.2f}')"
      2. Assert float in output: DR = 0.80
      3. Assert float in output: Precision = 1.00
    Expected Result: DR = 0.80, Precision = 1.00
    Evidence: .sisyphus/evidence/task-12-metrics-dr.txt

  Scenario: FPR computed from negative-control plasmids
    Tool: Bash (python)
    Preconditions: FPR computation implemented
    Steps:
      1. python -c "
  from tests.benchmark_literature import LiteratureBenchmark
  b = LiteratureBenchmark('pdif_only', '/tmp')
  metrics = b.compute_fpr(false_positives=5, total_plasmid_kb=500)
  print(f'FPR: {metrics:.4f} per kb')"
      2. Assert float in output: FPR = 0.0100
    Expected Result: FPR = 0.01 per kb
    Evidence: .sisyphus/evidence/task-12-metrics-fpr.txt

  Scenario: Full metrics JSON output is valid
    Tool: Bash (python)
    Preconditions: Metrics computation produces JSON
    Steps:
      1. python -c "import json; data=json.load(open('tests/benchmark_data/results/test_metrics.json')); assert 'detection_rate' in data; assert 'precision' in data; assert 'position_median_error' in data; assert 'fpr_per_kb' in data; print('Valid')"
      2. Assert exit code 0
    Expected Result: JSON contains all four metric fields
    Evidence: .sisyphus/evidence/task-12-metrics-json.json
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-12-metrics-dr.txt`
  - [ ] `.sisyphus/evidence/task-12-metrics-fpr.txt`
  - [ ] `.sisyphus/evidence/task-12-metrics-json.json`

  **Commit**: YES (groups with Tasks 10, 11)
  - Message: `feat(benchmark): implement four-metric computation with per-plasmid breakdown`
  - Files: `tests/benchmark_literature.py`

- [ ] 13. Run pdif_only benchmark on all cached plasmids and generate results

  **What to do**:
  - Run `python -m tests.benchmark_literature --mode pdif_only` on ALL cached plasmids
  - Process results into `tests/benchmark_data/results/pdif_only_results.json`
  - Verify metrics meet minimum thresholds (DR≥0.80, precision≥0.70, position median≤5bp, FPR≤0.01/kb)
  - Generate per-plasmid breakdown
  - If any plasmid has 0% detection rate, investigate and document the root cause
  - Save the run log to `tests/benchmark_data/results/pdif_only_run.log`

  **Must NOT do**:
  - Do NOT tune parameters to make metrics pass — report honestly
  - Do NOT exclude "difficult" plasmids from results

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
    - Reason: Running the full benchmark suite. Requires interpreting results, diagnosing failures, documenting findings.
  - **Skills**: []
    - Standard Python

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Tasks 14-18)
  - **Blocks**: Final Verification wave
  - **Blocked By**: Tasks 8, 10, 11, 12 (benchmark must be fully functional)

  **References**:
  - `tests/benchmark_literature.py` — the benchmark harness (built in Tasks 7-12)
  - `tests/benchmark_data/plasmids/` — cached plasmid sequences
  - `tests/benchmark_data/literature_reference.csv` — reference data (populated in Task 16)

  **Acceptance Criteria**:
  - [ ] All cached plasmids processed in pdif_only mode
  - [ ] Results JSON generated with all four metrics
  - [ ] Per-plasmid breakdown complete
  - [ ] Failure categorization applied to all missed sites

  **QA Scenarios**:

  ```
  Scenario: pdif_only benchmark processes all cached plasmids
    Tool: Bash (python)
    Preconditions: All prior tasks completed
    Steps:
      1. python -m tests.benchmark_literature --mode pdif_only --output-dir tests/benchmark_data/results --verbose 2>&1 | tee /tmp/pdif_only_run.log
      2. Assert exit code 0
      3. python -c "import json; data=json.load(open('tests/benchmark_data/results/pdif_only_results.json')); print(f'Plasmids: {data[\"plasmids_processed\"]}, DR: {data[\"detection_rate\"]:.2f}, Precision: {data[\"precision\"]:.2f}')"
      4. Assert data['plasmids_processed'] >= 10
      5. Check detection_rate >= 0.0 (just that it was computed, not necessarily passing threshold yet)
    Expected Result: All plasmids processed, results JSON generated
    Evidence: .sisyphus/evidence/task-13-pdif-only-results.json
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-13-pdif-only-results.json`

  **Commit**: YES
  - Message: `feat(benchmark): pdif_only benchmark results on all cached plasmids`
  - Files: `tests/benchmark_data/results/pdif_only_results.json`

- [ ] 14. Run full_pipeline benchmark on AMR-positive plasmid subset

  **What to do**:
  - Identify which cached plasmids have AMR genes detectable by pdifFinder's blastn database
  - Run `python -m tests.benchmark_literature --mode full_pipeline` on the AMR-positive subset
  - Process results into `tests/benchmark_data/results/full_pipeline_results.json`
  - Compare full_pipeline results to pdif_only results for the same plasmids
  - Document differences: does detecting real AMR genes change pdif detection?
  - If blastn is unavailable, document this and create placeholder results file with explanation

  **Must NOT do**:
  - Do NOT run full_pipeline on all plasmids — only those with verified AMR genes in the database
  - Do NOT force the benchmark if blastn is missing — document and skip

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
    - Reason: Running full pipeline, comparing results, documenting differences. Requires blastn availability check.
  - **Skills**: []
    - Standard Python; blastn external dependency

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Tasks 13, 15-18)
  - **Blocks**: Final Verification wave
  - **Blocked By**: Task 9 (full_pipeline mode), Task 12 (metrics)

  **References**:
  - `PdifFinder/AMRDB/sequences` — BLAST database for AMR detection
  - `tests/benchmark_data/results/pdif_only_results.json` — for comparison with pdif_only mode

  **Acceptance Criteria**:
  - [ ] Full_pipeline results JSON generated (or skip documented if blastn unavailable)
  - [ ] Results compared to pdif_only mode for same plasmids
  - [ ] Mode-comparison summary documented

  **QA Scenarios**:

  ```
  Scenario: full_pipeline produces comparable results to pdif_only
    Tool: Bash (python)
    Preconditions: blastn available, full_pipeline mode implemented
    Steps:
      1. python -m tests.benchmark_literature --mode full_pipeline --output-dir tests/benchmark_data/results --verbose 2>&1 | tee /tmp/full_pipeline_run.log
      2. Assert exit code 0
      3. python -c "
  import json
  fp = json.load(open('tests/benchmark_data/results/full_pipeline_results.json'))
  po = json.load(open('tests/benchmark_data/results/pdif_only_results.json'))
  print(f'Full pipeline DR: {fp[\"detection_rate\"]:.2f}')
  print(f'Pdif only DR: {po[\"detection_rate\"]:.2f}')
  print(f'Difference: {abs(fp[\"detection_rate\"] - po[\"detection_rate\"]):.2f}')"
      4. Assert results are comparable (difference <= 0.20)
    Expected Result: Detection rates comparable between modes
    Evidence: .sisyphus/evidence/task-14-full-pipeline-results.json
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-14-full-pipeline-results.json`

  **Commit**: YES
  - Message: `feat(benchmark): full_pipeline benchmark results on AMR-positive plasmids`
  - Files: `tests/benchmark_data/results/full_pipeline_results.json`

- [ ] 15. Write pytest integration tests for benchmark

  **What to do**:
  - Create `tests/validation/test_literature_benchmark.py`
  - Use `@pytest.mark.benchmark` marker (NOT default test suite)
  - Write tests that:
    1. `test_benchmark_harness_imports` — verify benchmark module imports without error
    2. `test_literature_reference_csv_valid` — run `validate_reference.py` on the reference CSV
    3. `test_pdif_only_mode_runs` — run pdif_only on a single plasmid, verify results have expected structure
    4. `test_full_pipeline_mode_runs` — run full_pipeline on a single plasmid (or skip if no blastn)
    5. `test_detection_rate_above_threshold` — verify detection_rate >= 0.80 on cached results
    6. `test_precision_above_threshold` — verify precision >= 0.70
    7. `test_position_median_error` — verify position median error <= 5bp
    8. `test_fpr_below_threshold` — verify FPR <= 0.01/kb
    9. `test_results_json_structure` — verify results JSON has all required fields
  - Use fixtures to load cached benchmark results (not re-run benchmark every time)
  - Mark tests as `@pytest.mark.slow` if they invoke pdifFinder

  **Must NOT do**:
  - Do NOT add to default pytest suite — ALL tests must use `@pytest.mark.benchmark`
  - Do NOT set xfail on any test (these should pass if benchmark is working)

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Standard pytest test writing. Straightforward test structure using existing results.
  - **Skills**: []
    - Standard pytest with fixtures and markers

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Tasks 13-14, 16-18)
  - **Blocks**: Final Verification wave
  - **Blocked By**: Tasks 8, 12 (benchmark must produce results for tests to verify)

  **References**:
  - `tests/validation/test_against_real_pdif.py` — pattern for validation tests with `@pytest.mark` decorators
  - `tests/conftest.py` — existing fixtures for test data paths
  - `pytest.ini` — configuration, add `benchmark` to markers list
  - `tests/benchmark_data/results/pdif_only_results.json` — results to verify in tests

  **Acceptance Criteria**:
  - [ ] All 9 tests defined with `@pytest.mark.benchmark`
  - [ ] `pytest tests/validation/test_literature_benchmark.py -m benchmark -v` runs and collects tests
  - [ ] Tests verify against cached results (not re-run benchmark)

  **QA Scenarios**:

  ```
  Scenario: Benchmark tests collect successfully
    Tool: Bash (pytest)
    Preconditions: test_literature_benchmark.py created
    Steps:
      1. pytest tests/validation/test_literature_benchmark.py -m benchmark --collect-only -v
      2. Assert exit code 0
      3. Assert output shows 9 test functions collected
    Expected Result: 9 benchmark tests collected
    Evidence: .sisyphus/evidence/task-15-tests-collect.txt

  Scenario: Detection rate threshold test passes
    Tool: Bash (pytest)
    Preconditions: Results JSON exists with DR >= 0.80
    Steps:
      1. pytest tests/validation/test_literature_benchmark.py::test_detection_rate_above_threshold -m benchmark -v
      2. Assert exit code 0
    Expected Result: Test passes
    Evidence: .sisyphus/evidence/task-15-dr-test.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-15-tests-collect.txt`
  - [ ] `.sisyphus/evidence/task-15-dr-test.txt`

  **Commit**: YES
  - Message: `test(benchmark): add pytest-integrated literature benchmark verification tests`
  - Files: `tests/validation/test_literature_benchmark.py`, `pytest.ini` (add benchmark marker)

- [ ] 16. Populate curated literature_reference.csv with real data

  **What to do**:
  - Merge data from Tasks 4, 5 (Cameranesi 2018, Blackwell & Hall 2017)
  - Add pdif coordinates for the redundant.seed.fa accessions (KY984047.1, CP045108.1, etc.)
  - For each redundant.seed.fa plasmid, BLASTN the seed pdif sequences against the plasmid to get coordinates
  - Add source metadata: source_pmid, source_description, curator_notes, curation_date
  - Validate the populated CSV using `validate_reference.py`
  - Ensure at minimum: 10 plasmids, ≥3 independent papers, ≥30 pdif site entries

  **Must NOT do**:
  - Do NOT include entries with unverified coordinates
  - Do NOT include entries where the pdif sequence doesn't match the plasmid at the stated position

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Data curation requiring cross-referencing multiple sources, BLAST validation, and manual verification. Critical for benchmark accuracy.
  - **Skills**: []
    - Standard Python for BLAST automation; manual data entry for final CSV

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Tasks 13-15, 17-18)
  - **Blocks**: Tasks 13, 14 (benchmark runs need populated reference data)
  - **Blocked By**: Tasks 2, 3, 4, 5, 6 (needs accessions, plasmids, literature coordinates, and schema)

  **References**:
  - `tests/benchmark_data/literature_reference.csv` — schema from Task 6
  - `tests/benchmark_data/literature_raw/cameranesi2018.csv` — from Task 4
  - `tests/benchmark_data/literature_raw/blackwell2017.csv` — from Task 5
  - `PdifFinder/data/redundant.seed.fa` — seed sequences with plasmid accessions
  - `tests/benchmark_data/plasmids/` — cached plasmid sequences for BLAST verification

  **Acceptance Criteria**:
  - [ ] `literature_reference.csv` contains ≥30 rows from ≥10 unique plasmids
  - [ ] ≥3 unique source_pmid values
  - [ ] `validate_reference.py` returns exit code 0
  - [ ] All coordinates verified against actual plasmid sequences

  **QA Scenarios**:

  ```
  Scenario: Reference CSV passes validation
    Tool: Bash (python)
    Preconditions: literature_reference.csv populated
    Steps:
      1. python tests/benchmark_data/validate_reference.py
      2. Assert exit code 0
      3. Assert output contains "All entries valid" or similar
    Expected Result: Validation passes on populated CSV
    Evidence: .sisyphus/evidence/task-16-validate-csv.txt

  Scenario: Reference CSV has sufficient coverage
    Tool: Bash (python)
    Preconditions: CSV populated
    Steps:
      1. python -c "
  import csv
  rows = list(csv.DictReader(open('tests/benchmark_data/literature_reference.csv')))
  plasmids = set(r['plasmid_accession'] for r in rows)
  pmids = set(r['source_pmid'] for r in rows)
  print(f'Plasmids: {len(plasmids)}, Papers: {len(pmids)}, Sites: {len(rows)}')
  assert len(plasmids) >= 10, f'Only {len(plasmids)} plasmids'
  assert len(pmids) >= 3, f'Only {len(pmids)} papers'
  assert len(rows) >= 30, f'Only {len(rows)} sites'"
      2. Assert exit code 0
    Expected Result: ≥10 plasmids, ≥3 papers, ≥30 pdif sites
    Evidence: .sisyphus/evidence/task-16-coverage.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-16-validate-csv.txt`
  - [ ] `.sisyphus/evidence/task-16-coverage.txt`

  **Commit**: YES
  - Message: `feat(benchmark): populate curated literature_reference.csv with real data`
  - Files: `tests/benchmark_data/literature_reference.csv`

- [ ] 17. Validate results JSON output format

  **What to do**:
  - Write a validation script `tests/benchmark_data/validate_results.py` that verifies:
    - Required top-level fields: `detection_rate`, `precision`, `position_median_error`, `position_mean_error`, `fpr_per_kb`, `plasmids_processed`, `total_reference_sites`, `matched_sites`, `missed_sites`, `unvalidated_detections`, `per_plasmid_breakdown`, `failure_categories`
    - All metric values are floats in valid ranges (0.0-1.0 for rates, 0-inf for position errors)
    - Per-plasmid breakdown has expected structure
    - Failure categories sum correctly (matched + missed = total_reference_sites)
  - Run validation on both `pdif_only_results.json` and `full_pipeline_results.json`
  - Also write a JSON Schema file `tests/benchmark_data/results/results_schema.json` for programmatic validation

  **Must NOT do**:
  - Do NOT validate metric thresholds here — that's Task 15 (pytest tests)
  - Do NOT modify the results JSONs — only validate them

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: JSON schema validation. Straightforward structural checking.
  - **Skills**: []
    - Standard Python with json module

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Tasks 13-16, 18)
  - **Blocks**: None
  - **Blocked By**: Tasks 13, 14 (needs results JSON to validate)

  **References**:
  - `tests/benchmark_data/results/pdif_only_results.json` — target for validation
  - `tests/benchmark_data/validate_reference.py` — pattern for validation script

  **Acceptance Criteria**:
  - [ ] `validate_results.py` runs successfully on both results JSONs
  - [ ] `results_schema.json` exists as JSON Schema

  **QA Scenarios**:

  ```
  Scenario: validate_results.py accepts valid results JSON
    Tool: Bash (python)
    Preconditions: Valid results JSON exists
    Steps:
      1. python tests/benchmark_data/validate_results.py tests/benchmark_data/results/pdif_only_results.json
      2. Assert exit code 0
    Expected Result: Validation passes
    Evidence: .sisyphus/evidence/task-17-validate-results.txt

  Scenario: validate_results.py rejects broken results JSON
    Tool: Bash (python)
    Preconditions: Validation script ready
    Steps:
      1. echo '{"detection_rate": "not_a_number"}' > /tmp/broken.json
      2. python tests/benchmark_data/validate_results.py /tmp/broken.json; echo "Exit: $?"
      3. Assert exit code != 0
    Expected Result: Validation catches bad data
    Evidence: .sisyphus/evidence/task-17-broken-results.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-17-validate-results.txt`
  - [ ] `.sisyphus/evidence/task-17-broken-results.txt`

  **Commit**: YES (groups with Task 18)
  - Message: `feat(benchmark): add results JSON validation script and JSON schema`
  - Files: `tests/benchmark_data/validate_results.py`, `tests/benchmark_data/results/results_schema.json`

- [ ] 18. Document benchmark usage in README and code headers

  **What to do**:
  - Add a "Literature Benchmark" section to `README.md` (or create `docs/benchmark_guide.md`)
  - Document:
    - Purpose: what the benchmark measures and why
    - Prerequisites: Python, BioPython, blastn (for full_pipeline), cached plasmids
    - Usage: both modes with example commands
    - Reference data: how it was curated, which papers, how to contribute new data
    - Metrics: definition of each metric, formula, acceptable thresholds
    - Output: JSON structure, how to interpret results
    - Reproducibility: pinned version, environment documentation
  - Add docstrings to all new benchmark functions in `tests/benchmark_literature.py`
  - Ensure all new files have proper module-level docstrings

  **Must NOT do**:
  - Do NOT claim metrics are "validated" without running the benchmark
  - Do NOT add benchmark docs to the main pdifFinder user documentation (separate section)

  **Recommended Agent Profile**:
  - **Category**: `writing`
    - Reason: Documentation writing with clear structure. Requires explaining technical concepts to researchers.
  - **Skills**: []
    - Documentation writing

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Tasks 13-17)
  - **Blocks**: None
  - **Blocked By**: None (can start immediately, though best done after implementation)

  **References**:
  - `README.md` — existing project README to extend
  - `tests/benchmark_literature.py` — functions to document
  - Metis review guardrails and directives — incorporate as usage instructions

  **Acceptance Criteria**:
  - [ ] README.md has Literature Benchmark section with all required documentation topics
  - [ ] All new benchmark functions have docstrings
  - [ ] Example commands are copy-paste runnable

  **QA Scenarios**:

  ```
  Scenario: Benchmark section exists in README
    Tool: Bash (grep)
    Preconditions: Documentation added
    Steps:
      1. grep -c "Literature Benchmark" README.md
      2. Assert count >= 1
      3. grep -c "pdif_only" README.md
      4. Assert count >= 1
    Expected Result: README contains benchmark documentation
    Evidence: .sisyphus/evidence/task-18-readme-check.txt

  Scenario: Docstrings exist on all benchmark classes/functions
    Tool: Bash (python)
    Preconditions: Docstrings added
    Steps:
      1. python -c "import tests.benchmark_literature as bl; print(bl.LiteratureBenchmark.__doc__[:50])"
      2. Assert output is not empty and not None
    Expected Result: All major classes have docstrings
    Evidence: .sisyphus/evidence/task-18-docstrings.txt
  ```

  **Evidence to Capture**:
  - [ ] `.sisyphus/evidence/task-18-readme-check.txt`
  - [ ] `.sisyphus/evidence/task-18-docstrings.txt`

  **Commit**: YES (groups with Task 17)
  - Message: `docs(benchmark): document literature benchmark usage and metrics`
  - Files: `README.md`, `tests/benchmark_literature.py` (docstrings) (MANDATORY — after ALL implementation tasks)

> 4 review agents run in PARALLEL. ALL must APPROVE. Present consolidated results to user and get explicit "okay" before completing.

- [ ] F1. **Plan Compliance Audit** — `oracle`
  Read the plan end-to-end. For each "Must Have": verify implementation exists. For each "Must NOT Have": search codebase for forbidden patterns — reject with file:line if found. Check evidence files exist in .sisyphus/evidence/. Compare deliverables against plan.
  Output: `Must Have [N/N] | Must NOT Have [N/N] | Tasks [N/N] | VERDICT: APPROVE/REJECT`

- [ ] F2. **Code Quality Review** — `unspecified-high`
  Run `python -m pytest tests/ --collect-only` + linter. Review all changed files for: `as any`/`@ts-ignore`, empty catches, console.log in prod, commented-out code, unused imports. Check AI slop: excessive comments, over-abstraction, generic names.
  Output: `Build [PASS/FAIL] | Lint [PASS/FAIL] | Tests [N pass/N fail] | Files [N clean/N issues] | VERDICT`

- [ ] F3. **Real Manual QA** — `unspecified-high`
  Start from clean state. Execute EVERY QA scenario from EVERY task. Test cross-task integration: benchmark harness produces correct JSON output on real plasmids. Test edge cases: circular plasmid wrap, missing AMR DB, invalid accession.
  Output: `Scenarios [N/N pass] | Integration [N/N] | Edge Cases [N tested] | VERDICT`

- [ ] F4. **Scope Fidelity Check** — `deep`
  For each task: read "What to do", read actual diff. Verify 1:1 — everything in spec was built, nothing beyond spec was built. Check "Must NOT do" compliance. Detect cross-task contamination. Flag unaccounted changes.
  Output: `Tasks [N/N compliant] | Contamination [CLEAN/N issues] | Unaccounted [CLEAN/N files] | VERDICT`

---

## Commit Strategy

- **Commit 1**: `fix(tests): correct parse_pdif_output column offset` - tests/helpers.py
  - Pre-commit: `pytest tests/ -x -k "parse or pipeline or validation"`
- **Commit 2**: `feat(benchmark): add literature reference data and plasmid cache` - tests/benchmark_data/
- **Commit 3**: `feat(benchmark): implement benchmark harness with dual modes` - tests/benchmark_literature.py
- **Commit 4**: `test(benchmark): add pytest-integrated literature benchmark tests` - tests/validation/test_literature_benchmark.py
- **Commit 5**: `docs(benchmark): document literature benchmark usage` - README.md

---

## Success Criteria

### Verification Commands
```bash
# Run pdif_only benchmark
python -m tests.benchmark_literature --mode pdif_only
# Expected: detects ≥80% of literature-reported pdif sites

# Run full_pipeline benchmark
python -m tests.benchmark_literature --mode full_pipeline
# Expected: comparable detection rate, generates per-plasmid failure report

# Run pytest benchmark suite
pytest tests/validation/test_literature_benchmark.py -m benchmark -v
# Expected: all tests pass (not xfail)

# Validate results JSON
python -c "import json; data=json.load(open('tests/benchmark_data/results/latest.json')); assert data['detection_rate']>=0.80"
```

### Final Checklist
- [ ] All "Must Have" present (dual modes, position comparison, failure categorization, cached data, reproducibility)
- [ ] All "Must NOT Have" absent (no source code modification, no seed DB changes, no visualizations, no web scraper)
- [ ] ≥10 literature plasmids curated, ≥3 independent papers
- [ ] Detection rate ≥0.80, precision ≥0.70, position median ≤5bp, FPR ≤0.01/kb
- [ ] All evidence files captured in `.sisyphus/evidence/`
- [ ] Final Verification Wave all APPROVE
