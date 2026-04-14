# Fix XerC and XerD Annotation Algorithm in PdifFinder

## TL;DR

> **Quick Summary**: Fix critical bug causing pdif site detection to only check first position of sequences, improve algorithm to scan entire plasmids, add test framework, and fix all identified code issues.
> 
> **Deliverables**: 
> - Fixed `pdifFinder.py` with corrected scanning algorithm
> - Comprehensive test suite with pytest framework
> - Updated documentation if needed
> - Validation scripts for pdif site detection
> 
> **Estimated Effort**: Medium (algorithm fixes + test framework)
> **Parallel Execution**: YES - 3 waves (setup, core fixes, testing/validation)
> **Critical Path**: Test framework → Bug fixes → Whole plasmid scanning → Integration tests

---

## Context

### Original Request
User asked: "XerC和XerD的标注是否存在问题？" (Are there issues with XerC and XerD annotations?)

### Interview Summary
**Key Discussions**:
- User suspected algorithm had defects - confirmed true by analysis
- Critical bug found: line 393 `for i in range(1):` only checks first position, should scan entire sequence
- Additional issues: variable name conflicts, spelling errors, logical issues
- Biological validation shows basic approach (28bp, XerC[0:11], XerD[17:28], mismatch thresholds) is scientifically sound
- User wants "comprehensive fix and improvement" not just minimal bug fix
- User wants to set up test framework (pytest/unittest) since none exists
- User decided: Scan entire plasmids (not just feature-end fragments)

**Research Findings**:
- **librarian**: Confirmed biological correctness - pdif sites are 28bp with XerC (0:11), 6bp spacer, XerD (17:28). Mismatch thresholds (XerC:3, XerD:2) are reasonable.
- **explore**: Found critical bug and additional code issues. Also identified algorithm limitation: only scans 28bp fragments around feature ends (`TAA`, `TGT`, `CAT`), not whole plasmids.

### Metis Review
**Identified Gaps** (addressed):
- **Scanning scope**: Current algorithm scans only feature-end fragments → User chose "scan entire plasmids"
- **Performance vs accuracy**: Will implement efficient scanning algorithm (not just brute-force fix)
- **Test data**: Will create synthetic test sequences with known pdif sites
- **Algorithm improvements**: Beyond bug fixes, will optimize scanning logic

---

## Work Objectives

### Core Objective
Fix the XerC/XerD annotation algorithm to correctly detect pdif sites throughout entire plasmid sequences, not just at the first position or feature ends.

### Concrete Deliverables
1. Fixed `pdifFinder.py` with corrected scanning algorithm
2. Test framework (pytest) with comprehensive test suite
3. Synthetic test sequences (positive/negative test cases)
4. Validation scripts to verify algorithm correctness
5. Updated documentation if needed

### Definition of Done
- [ ] `pytest tests/` runs all tests successfully (100% pass)
- [ ] Algorithm detects pdif sites at known positions in synthetic test sequences
- [ ] Algorithm scans entire plasmid sequences (not just first position or feature ends)
- [ ] All identified bugs (lines 393, 417, 474, etc.) are fixed
- [ ] Code quality improvements applied (spelling, variable conflicts, redundancy)

### Must Have
1. Fix critical bug: line 393 `for i in range(1):` → `for i in range(len(seq) - 27):`
2. Fix variable name conflict at line 417
3. Fix spelling: `maxMistach` → `maxMismatch`
4. Fix mathematical error at line 474
5. Remove redundant comparisons (lines 517-518)
6. Implement whole plasmid scanning (not just feature-end fragments)
7. Create test framework with comprehensive test suite
8. Thread-safe file operations

### Must NOT Have (Guardrails)
1. No changes to biological foundation (28bp pattern, XerC/XerD positions, mismatch thresholds) without explicit approval
2. No new features beyond algorithm fixes (e.g., new visualization, additional databases)
3. No over-engineering or excessive abstraction
4. No refactoring of unrelated code
5. No extensive documentation bloat beyond necessary comments

---

## Verification Strategy (MANDATORY)

> **ZERO HUMAN INTERVENTION** - ALL verification is agent-executable. No exceptions.

### Test Decision
- **Infrastructure exists**: NO (will be created)
- **Automated tests**: YES (TDD approach)
- **Framework**: pytest with synthetic test sequences
- **TDD approach**: Each fix will be test-driven: RED (failing test) → GREEN (minimal fix) → REFACTOR

### QA Policy
Every task MUST include agent-executed QA scenarios. Evidence saved to `.sisyphus/evidence/task-{N}-{scenario-slug}.{ext}`.

- **Python code**: Use pytest to run unit tests, capture output
- **CLI execution**: Use Bash to run pdifFinder commands, validate output files
- **Algorithm validation**: Use custom validation scripts to compare expected vs actual results

---

## Execution Strategy

### Parallel Execution Waves

> Maximize throughput by grouping independent tasks into parallel waves.
> Each wave completes before the next begins.
> Target: 5-8 tasks per wave.

```
Wave 1 (Foundation & Test Setup - Start Immediately):
├── Task 1: Set up pytest test framework
├── Task 2: Create synthetic test sequences
├── Task 3: Create validation utilities
├── Task 4: Analyze current algorithm limitations
├── Task 5: Design improved scanning algorithm
└── Task 6: Set up evidence directory structure

Wave 2 (Core Bug Fixes - After Wave 1):
├── Task 7: Fix critical bug (line 393 loop range)
├── Task 8: Fix variable name conflicts
├── Task 9: Fix spelling errors and code quality issues
├── Task 10: Fix mathematical error (line 474)
├── Task 11: Remove redundant comparisons
└── Task 12: Implement thread-safe file operations

Wave 3 (Algorithm Improvement - After Wave 2):
├── Task 13: Implement whole plasmid scanning
├── Task 14: Optimize scanning algorithm for performance
├── Task 15: Add comprehensive test coverage
├── Task 16: Create integration tests
├── Task 17: Validate against real pdif sequences
└── Task 18: Update documentation if needed

Wave FINAL (Verification - After ALL tasks):
├── Task F1: Plan compliance audit (oracle)
├── Task F2: Code quality review (unspecified-high)
├── Task F3: Real manual QA (unspecified-high)
└── Task F4: Scope fidelity check (deep)
-> Present results -> Get explicit user okay
```

### Dependency Matrix

- **1-6**: - - 7-18, 1
- **7**: 1, 4, 5 - 13, 14, 15, 2
- **8-12**: 1 - 13, 14, 15, 2
- **13**: 7, 8-12 - 16, 17, 3
- **14**: 7, 8-12 - 16, 17, 3
- **15**: 7, 8-12 - 16, 17, 3
- **16**: 13, 14, 15 - 18, 4
- **17**: 13, 14, 15 - 18, 4
- **18**: 16, 17 - F1-F4, FINAL

### Agent Dispatch Summary

- **1**: **6** - T1 → `quick`, T2 → `quick`, T3 → `quick`, T4 → `deep`, T5 → `deep`, T6 → `quick`
- **2**: **6** - T7 → `deep`, T8 → `quick`, T9 → `quick`, T10 → `deep`, T11 → `quick`, T12 → `deep`
- **3**: **6** - T13 → `deep`, T14 → `deep`, T15 → `unspecified-high`, T16 → `deep`, T17 → `unspecified-high`, T18 → `writing`
- **FINAL**: **4** - F1 → `oracle`, F2 → `unspecified-high`, F3 → `unspecified-high`, F4 → `deep`

---

## TODOs

> Implementation + Test = ONE Task. Never separate.
> EVERY task MUST have: Recommended Agent Profile + Parallelization info + QA Scenarios.
> **A task WITHOUT QA Scenarios is INCOMPLETE. No exceptions.**

- [x] 1. Set up pytest test framework

  **What to do**:
  - Install pytest and related dependencies (pytest-cov, pytest-mock)
  - Create `tests/` directory structure
  - Create `tests/conftest.py` with common fixtures
  - Create basic test configuration (`pytest.ini` or `setup.cfg`)
  - Verify pytest works by running `python -m pytest tests/ --collect-only`

  **Must NOT do**:
  - Do not modify existing production code
  - Do not add excessive dependencies beyond pytest ecosystem

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Simple setup task requiring standard Python tooling knowledge
  - **Skills**: None needed for basic pytest setup
  - **Skills Evaluated but Omitted**:
    - `git-master`: Not needed for test framework setup

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 2-6)
  - **Blocks**: Tasks 7-18 (all core fixes depend on test framework)
  - **Blocked By**: None (can start immediately)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `setup.py:59-64` - Dependencies list (check if pytest should be added)
  - Current project structure for `tests/` directory placement

  **External References** (libraries and frameworks):
  - Official docs: `https://docs.pytest.org/en/stable/getting-started.html` - Pytest setup and configuration

  **WHY Each Reference Matters**:
  - `setup.py` shows existing dependencies pattern to follow
  - Pytest docs provide correct configuration syntax

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] pytest installed and importable
  - [ ] `python -m pytest --version` → shows pytest version
  - [ ] `tests/` directory exists with `conftest.py`
  - [ ] `python -m pytest tests/ --collect-only` → shows 0 tests collected (expected at this stage)

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify pytest installation
    Tool: Bash
    Preconditions: Python environment active
    Steps:
      1. Run `python -c "import pytest; print(pytest.__version__)"`
      2. Check output contains version number (e.g., "7.4" or higher)
    Expected Result: pytest version printed without ImportError
    Failure Indicators: ImportError or "ModuleNotFoundError: No module named 'pytest'"
    Evidence: .sisyphus/evidence/task-1-pytest-version.txt

  Scenario: Verify test directory structure
    Tool: Bash
    Preconditions: In project root directory
    Steps:
      1. Run `ls -la tests/`
      2. Check `tests/conftest.py` exists
      3. Check `tests/__init__.py` exists (optional but recommended)
    Expected Result: tests/ directory contains at least conftest.py
    Failure Indicators: tests/ directory missing or empty
    Evidence: .sisyphus/evidence/task-1-test-dir.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-1-pytest-version.txt`: Output of pytest version check
  - [ ] `task-1-test-dir.txt`: Directory listing of tests/

  **Commit**: YES
  - Message: `build(tests): set up pytest test framework`
  - Files: `tests/`, `pytest.ini`/`setup.cfg`, `requirements-test.txt` (if created)
  - Pre-commit: None (first commit)

- [x] 2. Create synthetic test sequences

  **What to do**:
  - Create `tests/test_data/` directory
  - Generate synthetic DNA sequences (1000-5000bp) with known pdif sites at specific positions
  - Create positive test: sequence with 3 pdif sites at positions 100-127, 500-527, 1000-1027
  - Create negative test: sequence with no pdif sites
  - Create edge case test: sequence with partial/near-miss pdif sites
  - Save as FASTA files with descriptive names

  **Must NOT do**:
  - Do not use real patient/genomic data (synthetic only)
  - Do not create overly complex sequences (>10,000bp unless needed)

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Creating test data requires basic bioinformatics knowledge
  - **Skills**: None needed for synthetic sequence generation
  - **Skills Evaluated but Omitted**:
    - `playwright`: Not relevant for sequence generation

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 3-6)
  - **Blocks**: Tasks 15-17 (test coverage and validation depend on test data)
  - **Blocked By**: None (can start immediately)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `data/pdifdatabase.fasta:2-30` - Example pdif site sequences (28bp)
  - `data/redundant.seed.fa:1-30` - Example seed sequences

  **API/Type References** (contracts to implement against):
  - FASTA format: sequence header lines start with ">", sequence lines

  **Test References** (testing patterns to follow):
  - Need to match pdifFinder input format (FASTA)

  **External References** (libraries and frameworks):
  - Biopython for sequence generation: `from Bio.Seq import Seq`

  **WHY Each Reference Matters**:
  - `pdifdatabase.fasta` shows correct 28bp pdif site sequences to embed
  - FASTA format is required input for pdifFinder
  - Biopython provides reliable sequence generation

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] `tests/test_data/` directory exists
  - [ ] At least 3 test FASTA files created (positive, negative, edge)
  - [ ] Positive test contains exactly 3 pdif sites at specified positions
  - [ ] Files pass basic validation (can be read by `Bio.SeqIO`)

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Validate test sequences format
    Tool: Bash (Python script)
    Preconditions: In project root directory
    Steps:
      1. Run `python -c "from Bio import SeqIO; records = list(SeqIO.parse('tests/test_data/positive_test.fasta', 'fasta')); print(f'Read {len(records)} sequences')"`
      2. Check output shows "Read 1 sequences" (or more if multi-FASTA)
      3. Run `python -c "from Bio import SeqIO; record = next(SeqIO.parse('tests/test_data/positive_test.fasta', 'fasta')); print(f'Sequence length: {len(record.seq)}')"`
    Expected Result: Sequences readable, lengths > 1000bp
    Failure Indicators: Bio.SeqIO parse errors, sequences too short
    Evidence: .sisyphus/evidence/task-2-sequence-validation.txt

  Scenario: Verify pdif sites in positive test
    Tool: Bash (Python script)
    Preconditions: Positive test sequence created
    Steps:
      1. Extract positions 100-127, 500-527, 1000-1027 from sequence
      2. Compare extracted 28bp sequences to known pdif patterns
      3. Check each is 28bp and contains XerC[0:11] and XerD[17:28] pattern
    Expected Result: All 3 positions contain valid pdif sites (28bp with XerC/XerD)
    Failure Indicators: Wrong positions, incorrect sequences, missing pdif sites
    Evidence: .sisyphus/evidence/task-2-pdif-validation.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-2-sequence-validation.txt`: Sequence reading validation output
  - [ ] `task-2-pdif-validation.txt`: Pdif site verification output

  **Commit**: YES (group with Task 1)
  - Message: `test(data): add synthetic test sequences`
  - Files: `tests/test_data/*.fasta`
  - Pre-commit: `python -c "from Bio import SeqIO; SeqIO.parse('tests/test_data/positive_test.fasta', 'fasta')"`

- [x] 3. Create validation utilities

  **What to do**:
  - Create `tests/validate_scanning.py`: Script to verify algorithm scans entire sequence
  - Create `tests/validate_output.py`: Script to compare pdifFinder output to expected results
  - Create `tests/helpers.py`: Common test utilities (sequence generation, position checking)
  - Add command-line interface for validation scripts

  **Must NOT do**:
  - Do not modify pdifFinder.py production code
  - Do not create overly complex validation (keep simple and focused)

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Creating utility scripts requires basic Python skills
  - **Skills**: None needed for utility scripts
  - **Skills Evaluated but Omitted**:
    - `playwright`: Not relevant for validation scripts

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1-2, 4-6)
  - **Blocks**: Tasks 15-17 (test execution depends on validation utilities)
  - **Blocked By**: None (can start immediately)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - Existing Python files in project for coding style
  - `pdifFinder.py` for import patterns and structure

  **Test References** (testing patterns to follow):
  - Need to output PASS/FAIL with clear messages
  - Should handle edge cases gracefully

  **WHY Each Reference Matters**:
  - Existing code shows project's Python style conventions
  - Clear PASS/FAIL output enables automated validation

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] `tests/validate_scanning.py` exists and can be imported
  - [ ] `tests/validate_output.py` exists and can be imported
  - [ ] `tests/helpers.py` exists with utility functions
  - [ ] Scripts run without syntax errors

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Test validation script imports
    Tool: Bash
    Preconditions: In project root directory
    Steps:
      1. Run `python -c "import tests.validate_scanning; print('Import successful')"`
      2. Run `python -c "import tests.validate_output; print('Import successful')"`
      3. Run `python -c "import tests.helpers; print('Import successful')"`
    Expected Result: All imports succeed without errors
    Failure Indicators: ImportError, SyntaxError, ModuleNotFoundError
    Evidence: .sisyphus/evidence/task-3-import-test.txt

  Scenario: Test validation script execution
    Tool: Bash
    Preconditions: Validation scripts created
    Steps:
      1. Run `python tests/validate_scanning.py --help`
      2. Check output shows usage/help message
      3. Run `python tests/validate_output.py --help`
      4. Check output shows usage/help message
    Expected Result: Both scripts show help messages (or run without error)
    Failure Indicators: Script crashes, no help/usage shown
    Evidence: .sisyphus/evidence/task-3-script-execution.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-3-import-test.txt`: Import test output
  - [ ] `task-3-script-execution.txt`: Script execution output

  **Commit**: YES (group with Tasks 1-2)
  - Message: `test(utils): add validation utilities`
  - Files: `tests/validate_scanning.py`, `tests/validate_output.py`, `tests/helpers.py`
  - Pre-commit: `python -c "import tests.validate_scanning; import tests.validate_output; import tests.helpers"`

- [x] 4. Analyze current algorithm limitations

  **What to do**:
  - Read and understand `findMatchFragmentThread` function in detail
  - Document current limitations: only scans feature-end fragments, not whole plasmids
  - Analyze the bug: `for i in range(1):` only checks first position
  - Document variable name conflicts and other code quality issues
  - Create analysis report in `docs/algorithm_analysis.md`

  **Must NOT do**:
  - Do not modify production code yet
  - Do not change algorithm behavior

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Requires deep understanding of algorithm and bug analysis
  - **Skills**: None needed for analysis
  - **Skills Evaluated but Omitted**:
    - `git-master`: Not needed for analysis phase

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1-3, 5-6)
  - **Blocks**: Tasks 5, 7-14 (design and fixes depend on analysis)
  - **Blocked By**: None (can start immediately)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `pdifFinder.py:383-424` - `findMatchFragmentThread` function
  - `pdifFinder.py:208-248` - Thread dispatch and seed list processing
  - `pdifFinder.py:209-216` - Seed sequence extraction

  **Test References** (testing patterns to follow):
  - Current behavior needs to be documented before changes

  **WHY Each Reference Matters**:
  - `findMatchFragmentThread` is the core algorithm with the bug
  - Thread dispatch shows how algorithm is called
  - Seed extraction shows XerC/XerD pattern

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Analysis report created: `docs/algorithm_analysis.md`
  - [ ] Report documents: bug location, algorithm limitations, code quality issues
  - [ ] Report includes specific line numbers and descriptions

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify analysis report exists
    Tool: Bash
    Preconditions: In project root directory
    Steps:
      1. Check `docs/algorithm_analysis.md` file exists
      2. Read first 10 lines to confirm content
    Expected Result: Analysis report file exists with content
    Failure Indicators: File missing or empty
    Evidence: .sisyphus/evidence/task-4-analysis-exists.txt

  Scenario: Verify analysis covers key points
    Tool: Bash (grep)
    Preconditions: Analysis report created
    Steps:
      1. Run `grep -i "range(1)" docs/algorithm_analysis.md`
      2. Run `grep -i "variable.*conflict" docs/algorithm_analysis.md`
      3. Run `grep -i "feature.*end.*fragment" docs/algorithm_analysis.md`
    Expected Result: All grep commands find matches (analysis covers these issues)
    Failure Indicators: One or more grep commands return no matches
    Evidence: .sisyphus/evidence/task-4-analysis-content.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-4-analysis-exists.txt`: File existence check
  - [ ] `task-4-analysis-content.txt`: Grep results showing analysis coverage

  **Commit**: YES (group with Tasks 1-3)
  - Message: `docs: analyze algorithm limitations and bugs`
  - Files: `docs/algorithm_analysis.md`
  - Pre-commit: None

- [x] 5. Design improved scanning algorithm

  **What to do**:
  - Design algorithm to scan entire plasmid sequences (not just feature ends)
  - Consider performance optimizations (Boyer-Moore, sliding window)
  - Design solution for whole plasmid scanning while maintaining compatibility
  - Create design document: `docs/algorithm_design.md`
  - Include pseudocode for improved algorithm

  **Must NOT do**:
  - Do not implement design yet (design only)
  - Do not change biological parameters (28bp, XerC/XerD positions, mismatch thresholds)

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Algorithm design requires deep problem-solving
  - **Skills**: None needed for design
  - **Skills Evaluated but Omitted**:
    - `playwright`: Not relevant for algorithm design

  **Parallelization**:
  - **Can Run In Parallel**: YES (but depends on Task 4 analysis)
  - **Parallel Group**: Wave 1 (with Tasks 1-4, 6)
  - **Blocks**: Tasks 7, 13-14 (implementation depends on design)
  - **Blocked By**: Task 4 (needs analysis to inform design)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `pdifFinder.py:383-424` - Current algorithm to improve
  - `pdifFinder.py:527-570` - BM (Boyer-Moore) function for pattern matching

  **External References** (libraries and frameworks):
  - Boyer-Moore algorithm literature for efficient string searching
  - Biopython sequence searching patterns

  **WHY Each Reference Matters**:
  - Current algorithm shows what needs to be modified
  - BM function shows existing efficient searching capability in codebase
  - Boyer-Moore literature informs performance optimization

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Design document created: `docs/algorithm_design.md`
  - [ ] Design includes: whole plasmid scanning approach, performance considerations
  - [ ] Design includes pseudocode for improved algorithm
  - [ ] Design maintains backward compatibility where possible

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify design document exists
    Tool: Bash
    Preconditions: In project root directory
    Steps:
      1. Check `docs/algorithm_design.md` file exists
      2. Check file size > 1KB (meaningful content)
    Expected Result: Design document exists with content
    Failure Indicators: File missing or very small (< 500 bytes)
    Evidence: .sisyphus/evidence/task-5-design-exists.txt

  Scenario: Verify design covers whole plasmid scanning
    Tool: Bash (grep)
    Preconditions: Design document created
    Steps:
      1. Run `grep -i "whole.*plasmid\|entire.*sequence" docs/algorithm_design.md`
      2. Run `grep -i "sliding.*window\|scan.*all" docs/algorithm_design.md`
    Expected Result: Grep finds matches (design covers whole plasmid scanning)
    Failure Indicators: No matches found
    Evidence: .sisyphus/evidence/task-5-design-content.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-5-design-exists.txt`: File existence and size check
  - [ ] `task-5-design-content.txt`: Grep results showing design coverage

  **Commit**: YES (group with Tasks 1-4)
  - Message: `docs: design improved scanning algorithm`
  - Files: `docs/algorithm_design.md`
  - Pre-commit: None

- [x] 6. Set up evidence directory structure

  **What to do**:
  - Create `.sisyphus/evidence/` directory if not exists
  - Create subdirectories: `task-evidence/`, `validation/`, `screenshots/`
  - Create README explaining evidence structure
  - Ensure proper permissions for file creation

  **Must NOT do**:
  - Do not modify existing evidence files
  - Do not create excessive directory structure

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Simple directory setup task
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - All skills: Not needed for directory creation

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1-5)
  - **Blocks**: All tasks with QA scenarios (need evidence directory)
  - **Blocked By**: None (can start immediately)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - Existing `.sisyphus/` directory structure

  **WHY Each Reference Matters**:
  - Follows existing project organization patterns

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] `.sisyphus/evidence/` directory exists
  - [ ] Subdirectories created: `task-evidence/`, `validation/`, `screenshots/`
  - [ ] Directory writable by current user

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify evidence directory structure
    Tool: Bash
    Preconditions: In project root directory
    Steps:
      1. Run `ls -la .sisyphus/evidence/`
      2. Check `task-evidence/`, `validation/`, `screenshots/` directories exist
      3. Run `touch .sisyphus/evidence/test.txt && rm .sisyphus/evidence/test.txt`
    Expected Result: Directories exist and are writable
    Failure Indicators: Directories missing, permission errors
    Evidence: .sisyphus/evidence/task-6-directory-check.txt

  Scenario: Test evidence file creation
    Tool: Bash
    Preconditions: Evidence directories created
    Steps:
      1. Run `echo "test" > .sisyphus/evidence/task-evidence/test-file.txt`
      2. Check file created: `cat .sisyphus/evidence/task-evidence/test-file.txt`
      3. Cleanup: `rm .sisyphus/evidence/task-evidence/test-file.txt`
    Expected Result: File created, contains "test", can be removed
    Failure Indicators: Permission errors, file not created
    Evidence: .sisyphus/evidence/task-6-write-test.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-6-directory-check.txt`: Directory listing and permissions check
  - [ ] `task-6-write-test.txt`: File creation test output

  **Commit**: NO (evidence directory is temporary, not committed to repo)
  - Message: N/A
  - Files: N/A
  - Pre-commit: N/A

- [x] 7. Fix critical bug (line 393 loop range)

  **What to do**:
  - Locate `findMatchFragmentThread` function in `pdifFinder.py`
  - Change line 393 from `for i in range(1):` to `for i in range(len(seq) - 27):`
  - Ensure loop iterates over all possible starting positions for 28bp pdif sites
  - Update related variables: `initPos = i`, `endPos = i + 11`
  - Verify boundary conditions: `(endPos + 17) <= len(seq)` should check `i + 28 <= len(seq)`
  - Test the fix with synthetic test sequences

  **Must NOT do**:
  - Do not change biological parameters (28bp, XerC/XerD positions, mismatch thresholds)
  - Do not modify unrelated parts of the function

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Requires careful understanding of algorithm and boundary conditions
  - **Skills**: None needed for code modification
  - **Skills Evaluated but Omitted**:
    - `git-master`: Not needed for this specific fix

  **Parallelization**:
  - **Can Run In Parallel**: NO (depends on Tasks 4-5 analysis and design)
  - **Parallel Group**: Wave 2 (with Tasks 8-12)
  - **Blocks**: Tasks 13-14 (algorithm improvements depend on this fix)
  - **Blocked By**: Tasks 4-5 (need analysis and design completed)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `pdifFinder.py:383-424` - `findMatchFragmentThread` function
  - `pdifFinder.py:393` - Current buggy line: `for i in range(1):`
  - `pdifFinder.py:396` - `initPos = i` and `endPos = i + 11`

  **Test References** (testing patterns to follow):
  - Use synthetic test sequences from Task 2
  - Use validation utilities from Task 3

  **WHY Each Reference Matters**:
  - `findMatchFragmentThread` is the function to modify
  - Line 393 is the exact location of the bug
  - Lines 396 show how `initPos` and `endPos` are used

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Line 393 changed to `for i in range(len(seq) - 27):`
  - [ ] `python -m pytest tests/test_bug_fixes.py::test_loop_range` passes
  - [ ] Algorithm detects pdif sites at all positions in synthetic test sequence

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify line 393 fix
    Tool: Bash (grep)
    Preconditions: pdifFinder.py modified
    Steps:
      1. Run `grep -n "for i in range" pdifFinder/pdifFinder.py | grep "393"`
      2. Check output contains "for i in range(len(seq) - 27):"
    Expected Result: Line 393 shows correct loop range
    Failure Indicators: Line 393 unchanged or shows different range
    Evidence: .sisyphus/evidence/task-7-line-393-fix.txt

  Scenario: Test algorithm detects pdif sites at multiple positions
    Tool: Bash (Python script)
    Preconditions: Synthetic test sequence with pdif sites at positions 100, 500, 1000
    Steps:
      1. Run `python tests/validate_scanning.py --test-multiple-positions`
      2. Check output contains "All 3 pdif sites detected"
    Expected Result: Algorithm detects all 3 pdif sites
    Failure Indicators: Algorithm misses some pdif sites
    Evidence: .sisyphus/evidence/task-7-multi-position-test.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-7-line-393-fix.txt`: Grep output showing fixed line
  - [ ] `task-7-multi-position-test.txt`: Test output showing pdif site detection

  **Commit**: YES
  - Message: `fix(pdifFinder): correct loop range in findMatchFragmentThread`
  - Files: `pdifFinder/pdifFinder.py`
  - Pre-commit: `python -m pytest tests/test_bug_fixes.py::test_loop_range`

- [x] 8. Fix variable name conflicts

  **What to do**:
  - Locate `findMatchFragmentThread` function in `pdifFinder.py`
  - Change line 417 variable name from `k` to `m` (or another unused name)
  - Line 417: `for k in range(11):` → `for m in range(11):`
  - Update variable usage inside loop: `if tempSeqD[m] != seqd[m]:`
  - Ensure no other variable conflicts exist in the function

  **Must NOT do**:
  - Do not rename outer loop variable `k` at line 387
  - Do not modify algorithm logic beyond variable renaming

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Simple variable renaming task
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - All skills: Not needed for variable renaming

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 7, 9-12 in Wave 2)
  - **Parallel Group**: Wave 2 (with Tasks 7, 9-12)
  - **Blocks**: Tasks 13-14 (algorithm improvements)
  - **Blocked By**: Task 4 (analysis should identify this issue)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `pdifFinder.py:383-424` - `findMatchFragmentThread` function
  - `pdifFinder.py:387` - Outer loop: `for k in range(start, end, 2):`
  - `pdifFinder.py:417` - Inner loop conflict: `for k in range(11):`

  **WHY Each Reference Matters**:
  - Shows the variable conflict between outer and inner loops
  - Line 417 is the exact location to fix

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Line 417 changed to use different variable name (e.g., `m`)
  - [ ] Inner loop uses new variable name consistently
  - [ ] `python -m pytest tests/test_bug_fixes.py::test_variable_conflict` passes

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify variable conflict fix
    Tool: Bash (grep)
    Preconditions: pdifFinder.py modified
    Steps:
      1. Run `grep -n "for .* in range(11):" pdifFinder/pdifFinder.py`
      2. Check output does NOT contain "for k in range(11):"
      3. Check output contains "for m in range(11):" (or other non-conflicting name)
    Expected Result: No variable conflict (inner loop uses different variable than outer loop)
    Failure Indicators: Both loops use same variable name `k`
    Evidence: .sisyphus/evidence/task-8-variable-fix.txt

  Scenario: Test algorithm still works after variable rename
    Tool: Bash (Python script)
    Preconditions: Synthetic test sequence
    Steps:
      1. Run `python tests/validate_output.py --regression-test`
      2. Check output contains "Variable rename: PASS"
    Expected Result: Algorithm functions correctly after variable rename
    Failure Indicators: Algorithm broken due to variable rename
    Evidence: .sisyphus/evidence/task-8-regression-test.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-8-variable-fix.txt`: Grep output showing no conflict
  - [ ] `task-8-regression-test.txt`: Regression test output

  **Commit**: YES (group with Task 7)
  - Message: `fix(pdifFinder): resolve variable name conflict in findMatchFragmentThread`
  - Files: `pdifFinder/pdifFinder.py`
  - Pre-commit: `python -m pytest tests/test_bug_fixes.py::test_variable_conflict`

- [x] 9. Fix spelling errors and code quality issues

  **What to do**:
  - Locate `maxMistachXerC` and `maxMistachXerD` variables in `pdifFinder.py`
  - Change to `maxMismatchXerC` and `maxMismatchXerD` (correct spelling)
  - Update all references to these variables (lines 391-392, 404, 412, 420)
  - Check for other spelling errors in the function
  - Improve code comments if needed

  **Must NOT do**:
  - Do not change variable values (3 and 2 thresholds)
  - Do not modify algorithm logic

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Simple spelling correction
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - All skills: Not needed for spelling fixes

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 7-8, 10-12 in Wave 2)
  - **Parallel Group**: Wave 2 (with Tasks 7-8, 10-12)
  - **Blocks**: Tasks 13-14 (algorithm improvements)
  - **Blocked By**: Task 4 (analysis should identify this issue)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `pdifFinder.py:391` - `maxMistachXerC = 3`
  - `pdifFinder.py:392` - `maxMistachXerD = 2`
  - `pdifFinder.py:404` - `if tempMistachNumberCLeft > maxMistachXerC:`
  - `pdifFinder.py:412` - `if tempMistachNumberC > maxMistachXerC:`
  - `pdifFinder.py:420` - `if tempMistachNumberD > maxMistachXerD:`

  **WHY Each Reference Matters**:
  - Shows all locations where spelling needs correction

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] All instances of `maxMistach` changed to `maxMismatch`
  - [ ] Variables retain original values (3 and 2)
  - [ ] `python -m pytest tests/test_bug_fixes.py::test_spelling_fixes` passes

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify spelling corrections
    Tool: Bash (grep)
    Preconditions: pdifFinder.py modified
    Steps:
      1. Run `grep -n "maxMistach" pdifFinder/pdifFinder.py`
      2. Check output is empty (no misspelled instances remain)
      3. Run `grep -n "maxMismatch" pdifFinder/pdifFinder.py`
      4. Check output shows 5 instances (391, 392, 404, 412, 420)
    Expected Result: All misspellings corrected, correct spelling used
    Failure Indicators: `maxMistach` still found, or `maxMismatch` not found
    Evidence: .sisyphus/evidence/task-9-spelling-fix.txt

  Scenario: Test algorithm still works after spelling fix
    Tool: Bash (Python script)
    Preconditions: Synthetic test sequence
    Steps:
      1. Run `python tests/validate_output.py --regression-test`
      2. Check output contains "Spelling fix: PASS"
    Expected Result: Algorithm functions correctly after spelling corrections
    Failure Indicators: Algorithm broken due to variable rename
    Evidence: .sisyphus/evidence/task-9-regression-test.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-9-spelling-fix.txt`: Grep output showing spelling corrections
  - [ ] `task-9-regression-test.txt`: Regression test output

  **Commit**: YES (group with Tasks 7-8)
  - Message: `style(pdifFinder): fix spelling of maxMismatch variables`
  - Files: `pdifFinder/pdifFinder.py`
  - Pre-commit: `python -m pytest tests/test_bug_fixes.py::test_spelling_fixes`

- [x] 10. Fix mathematical error (line 474)

  **What to do**:
  - Locate line 474 in `pdifFinder.py`: `restNumber = batchLength - 50 * batchNumber`
  - Change `50` to `batch` (variable name): `restNumber = batchLength - batch * batchNumber`
  - Verify the calculation: `restNumber` should be the remainder after dividing `batchLength` into batches of size `batch`
  - Test with sample values to ensure correct calculation

  **Must NOT do**:
  - Do not change algorithm logic beyond this correction
  - Do not modify unrelated parts of the function

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Simple mathematical correction
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - All skills: Not needed for mathematical correction

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 7-9, 11-12 in Wave 2)
  - **Parallel Group**: Wave 2 (with Tasks 7-9, 11-12)
  - **Blocks**: None (code quality fix, doesn't block other tasks)
  - **Blocked By**: Task 4 (analysis should identify this issue)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `pdifFinder.py:474` - Line with mathematical error
  - `pdifFinder.py:470-476` - Context around the error

  **WHY Each Reference Matters**:
  - Line 474 is the exact location of the error
  - Context shows how `restNumber` is used

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Line 474 changed to `restNumber = batchLength - batch * batchNumber`
  - [ ] `python -m pytest tests/test_bug_fixes.py::test_math_error` passes
  - [ ] Calculation produces correct remainder for sample inputs

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify mathematical error fix
    Tool: Bash (grep)
    Preconditions: pdifFinder.py modified
    Steps:
      1. Run `grep -n "restNumber = batchLength - " pdifFinder/pdifFinder.py`
      2. Check output contains "restNumber = batchLength - batch * batchNumber"
      3. Verify `50` is not in the line
    Expected Result: Line shows correct formula with `batch` not `50`
    Failure Indicators: Line unchanged or still contains `50`
    Evidence: .sisyphus/evidence/task-10-math-fix.txt

  Scenario: Test calculation correctness
    Tool: Bash (Python script)
    Preconditions: Mathematical fix applied
    Steps:
      1. Run `python tests/validate_math.py --test-batch-calculation`
      2. Check output contains "Batch calculation: PASS"
    Expected Result: Calculation produces correct remainders
    Failure Indicators: Incorrect remainder calculation
    Evidence: .sisyphus/evidence/task-10-calculation-test.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-10-math-fix.txt`: Grep output showing fixed line
  - [ ] `task-10-calculation-test.txt`: Calculation test output

  **Commit**: YES (group with Tasks 7-9)
  - Message: `fix(pdifFinder): correct mathematical error in batch calculation`
  - Files: `pdifFinder/pdifFinder.py`
  - Pre-commit: `python -m pytest tests/test_bug_fixes.py::test_math_error`

- [x] 11. Remove redundant comparisons

  **What to do**:
  - Locate lines 517-518 in `pdifFinder.py`
  - Identify redundant comparisons: `misMatch3 = misMatch2` and `misMatch4 = misMatch1`
  - Remove these redundant assignments if they serve no purpose
  - Verify removal doesn't break any downstream logic
  - Check for other redundant code in the same function

  **Must NOT do**:
  - Do not remove non-redundant comparisons
  - Do not modify algorithm logic

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Simple code cleanup
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - All skills: Not needed for code cleanup

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 7-10, 12 in Wave 2)
  - **Parallel Group**: Wave 2 (with Tasks 7-10, 12)
  - **Blocks**: None (code quality fix)
  - **Blocked By**: Task 4 (analysis should identify this issue)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `pdifFinder.py:517-518` - Redundant comparison lines
  - `pdifFinder.py:515-520` - Context around redundant code

  **WHY Each Reference Matters**:
  - Lines 517-518 are the exact locations of redundant code
  - Context shows if these variables are used later

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Redundant lines 517-518 removed or commented out
  - [ ] `python -m pytest tests/test_bug_fixes.py::test_redundant_code` passes
  - [ ] Algorithm functions correctly without redundant code

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify redundant code removal
    Tool: Bash (grep)
    Preconditions: pdifFinder.py modified
    Steps:
      1. Run `grep -n "misMatch3 = misMatch2" pdifFinder/pdifFinder.py`
      2. Run `grep -n "misMatch4 = misMatch1" pdifFinder/pdifFinder.py`
    Expected Result: Both grep commands return no matches (lines removed)
    Failure Indicators: Lines still present
    Evidence: .sisyphus/evidence/task-11-redundant-removal.txt

  Scenario: Test algorithm after redundant code removal
    Tool: Bash (Python script)
    Preconditions: Redundant code removed
    Steps:
      1. Run `python tests/validate_output.py --regression-test`
      2. Check output contains "Redundant code removal: PASS"
    Expected Result: Algorithm functions correctly
    Failure Indicators: Algorithm broken after removal
    Evidence: .sisyphus/evidence/task-11-regression-test.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-11-redundant-removal.txt`: Grep output showing lines removed
  - [ ] `task-11-regression-test.txt`: Regression test output

  **Commit**: YES (group with Tasks 7-10)
  - Message: `refactor(pdifFinder): remove redundant comparisons`
  - Files: `pdifFinder/pdifFinder.py`
  - Pre-commit: `python -m pytest tests/test_bug_fixes.py::test_redundant_code`

- [x] 12. Implement thread-safe file operations

  **What to do**:
  - Analyze file operations in `findMatchFragmentThread` and related functions
  - Identify multiple threads writing to same files without locking
  - Implement thread-safe file writing (e.g., using `threading.Lock`, queue, or separate files)
  - Ensure no data corruption or race conditions
  - Test with multi-threaded execution

  **Must NOT do**:
  - Do not remove multi-threading (keep performance benefit)
  - Do not significantly reduce algorithm performance

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Requires understanding of threading and file I/O
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - `git-master`: Not needed for threading fixes

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 7-11 in Wave 2)
  - **Parallel Group**: Wave 2 (with Tasks 7-11)
  - **Blocks**: Tasks 13-18 (algorithm improvements should use thread-safe operations)
  - **Blocked By**: Task 4 (analysis should identify this issue)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `pdifFinder.py:423-424` - File write operation in thread: `with open(outfile, 'a') as w:`
  - `pdifFinder.py:245-248` - Thread start/join pattern
  - `pdifFinder.py:272-280` - Other file operations that may need thread safety

  **External References** (libraries and frameworks):
  - Python `threading.Lock` documentation for thread safety
  - Best practices for thread-safe file I/O in Python

  **WHY Each Reference Matters**:
  - Line 423 shows the file write that needs protection
  - Thread start/join shows concurrency pattern
  - Other file operations may also need protection

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Thread-safe file operations implemented (e.g., using `Lock`)
  - [ ] `python -m pytest tests/test_thread_safety.py` passes
  - [ ] No data corruption in multi-threaded tests

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify thread-safe file operations
    Tool: Bash (Python script)
    Preconditions: Thread-safe implementation applied
    Steps:
      1. Run `python tests/test_thread_safety.py --test-file-writes`
      2. Check output contains "Thread-safe file operations: PASS"
    Expected Result: No data corruption in concurrent file writes
    Failure Indicators: Missing data, corrupted writes, race conditions
    Evidence: .sisyphus/evidence/task-12-thread-safety-test.txt

  Scenario: Test algorithm with multiple threads
    Tool: Bash (Python script)
    Preconditions: Thread-safe implementation
    Steps:
      1. Run `python tests/validate_output.py --test-concurrent`
      2. Check output contains "Concurrent execution: PASS"
    Expected Result: Algorithm works correctly with multiple threads
    Failure Indicators: Crashes, incorrect results with concurrency
    Evidence: .sisyphus/evidence/task-12-concurrent-test.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-12-thread-safety-test.txt`: Thread safety test output
  - [ ] `task-12-concurrent-test.txt`: Concurrent execution test output

  **Commit**: YES (group with Tasks 7-11)
  - Message: `fix(pdifFinder): implement thread-safe file operations`
  - Files: `pdifFinder/pdifFinder.py`
  - Pre-commit: `python -m pytest tests/test_thread_safety.py`

- [x] 13. Implement whole plasmid scanning

  **What to do**:
  - Modify algorithm to scan entire plasmid sequences instead of just feature-end fragments
  - Update `findFeatureEndSeq` or bypass its fragment extraction
  - Ensure algorithm processes full sequences from input files
  - Maintain compatibility with existing output formats
  - Test with full-length synthetic plasmids

  **Must NOT do**:
  - Do not break existing functionality for feature-end scanning
  - Do not significantly reduce performance without optimization

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Major algorithm change requiring careful design implementation
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - `playwright`: Not relevant for algorithm changes

  **Parallelization**:
  - **Can Run In Parallel**: NO (depends on Tasks 7-12 fixes)
  - **Parallel Group**: Wave 3 (with Tasks 14-18)
  - **Blocks**: Tasks 16-17 (integration tests and validation)
  - **Blocked By**: Tasks 7-12 (core fixes must be complete)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `pdifFinder.py:220` - `findFeatureEndSeq` function call
  - `pdifFinder.py:527-570` - BM function for efficient searching
  - `docs/algorithm_design.md` - Design document from Task 5

  **WHY Each Reference Matters**:
  - `findFeatureEndSeq` is the current fragment extraction to modify
  - BM function shows efficient searching pattern
  - Design document guides implementation

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Algorithm processes full plasmid sequences
  - [ ] `python -m pytest tests/test_whole_plasmid.py` passes
  - [ ] Detects pdif sites throughout entire sequences

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify whole plasmid scanning
    Tool: Bash (Python script)
    Preconditions: Whole plasmid scanning implemented
    Steps:
      1. Run `python tests/validate_scanning.py --test-whole-plasmid`
      2. Check output contains "Whole plasmid scanning: PASS"
    Expected Result: Algorithm scans entire sequences
    Failure Indicators: Only scans fragments or misses sites
    Evidence: .sisyphus/evidence/task-13-whole-plasmid-test.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-13-whole-plasmid-test.txt`: Whole plasmid scanning test output

  **Commit**: YES
  - Message: `feat(pdifFinder): implement whole plasmid scanning`
  - Files: `pdifFinder/pdifFinder.py`
  - Pre-commit: `python -m pytest tests/test_whole_plasmid.py`

- [x] 14. Optimize scanning algorithm for performance

  **What to do**:
  - Implement performance optimizations (Boyer-Moore, sliding window optimizations)
  - Profile algorithm to identify bottlenecks
  - Reduce unnecessary computations and memory usage
  - Ensure optimizations don't affect accuracy
  - Benchmark before/after performance

  **Must NOT do**:
  - Do not introduce bugs while optimizing
  - Do not change biological matching logic

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Performance optimization requires algorithmic expertise
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - All skills: Not needed for optimization

  **Parallelization**:
  - **Can Run In Parallel**: NO (depends on Task 13)
  - **Parallel Group**: Wave 3 (with Tasks 13, 15-18)
  - **Blocks**: Tasks 16-17 (integration and validation)
  - **Blocked By**: Task 13 (whole plasmid scanning must work first)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `pdifFinder.py:527-570` - Existing BM function for pattern matching
  - Performance profiling tools (cProfile, timeit)

  **External References** (libraries and frameworks):
  - Boyer-Moore algorithm literature
  - Python optimization best practices

  **WHY Each Reference Matters**:
  - BM function shows existing optimization pattern
  - Profiling tools identify bottlenecks

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Performance improved by at least 20% on benchmark sequences
  - [ ] `python -m pytest tests/test_performance.py` passes
  - [ ] Accuracy maintained (no regressions)

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify performance improvement
    Tool: Bash (Python script)
    Preconditions: Optimizations implemented
    Steps:
      1. Run `python tests/benchmark.py --compare-before-after`
      2. Check output contains "Performance improvement: YES"
    Expected Result: Measurable performance improvement
    Failure Indicators: Performance degraded or unchanged
    Evidence: .sisyphus/evidence/task-14-performance-test.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-14-performance-test.txt`: Performance benchmark results

  **Commit**: YES
  - Message: `perf(pdifFinder): optimize scanning algorithm performance`
  - Files: `pdifFinder/pdifFinder.py`
  - Pre-commit: `python -m pytest tests/test_performance.py`

- [ ] 15. Add comprehensive test coverage

  **What to do**:
  - Create comprehensive test suite covering all bug fixes and improvements
  - Add unit tests for each fixed bug (loop range, variable conflict, spelling, math error, redundant code)
  - Add integration tests for whole plasmid scanning and performance optimizations
  - Add edge case tests (empty sequences, no pdif sites, multiple pdif sites)
  - Ensure test coverage > 80% for modified code

  **Must NOT do**:
  - Do not write tests for unrelated code
  - Do not create flaky or non-deterministic tests

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
    - Reason: Comprehensive test suite requires thorough understanding of all changes
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - All skills: Not needed for test writing

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 13-14, 16-18 in Wave 3)
  - **Parallel Group**: Wave 3 (with Tasks 13-14, 16-18)
  - **Blocks**: Tasks 16-17 (integration tests part of this task)
  - **Blocked By**: Tasks 7-14 (need all fixes to test)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - Existing test patterns in project (none, but follow pytest conventions)
  - Test structure from Tasks 1-3

  **WHY Each Reference Matters**:
  - Need to follow consistent test patterns

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Test coverage > 80% for modified `pdifFinder.py` functions
  - [ ] `python -m pytest tests/ --cov=pdifFinder --cov-report=term-missing` shows adequate coverage
  - [ ] All tests pass

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Verify test coverage
    Tool: Bash
    Preconditions: Comprehensive test suite created
    Steps:
      1. Run `python -m pytest tests/ --cov=pdifFinder --cov-report=term-missing`
      2. Check coverage report shows >80% for modified functions
    Expected Result: Adequate test coverage achieved
    Failure Indicators: Coverage below 80%
    Evidence: .sisyphus/evidence/task-15-coverage.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-15-coverage.txt`: Test coverage report

  **Commit**: YES
  - Message: `test(pdifFinder): add comprehensive test coverage`
  - Files: `tests/`
  - Pre-commit: `python -m pytest tests/`

- [ ] 16. Create integration tests

  **What to do**:
  - Create integration tests that test the full pdifFinder pipeline
  - Test end-to-end: FASTA input → pdifFinder execution → output validation
  - Test with various input types (single sequence, multi-FASTA, GenBank)
  - Verify output files (pdif_site.txt, pdif_pair.txt) are correctly generated
  - Test cross-platform compatibility

  **Must NOT do**:
  - Do not test external dependencies (BLAST) unless necessary
  - Do not create slow integration tests (> 1 minute each)

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Integration testing requires understanding of full pipeline
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - All skills: Not needed for integration tests

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 13-15, 17-18 in Wave 3)
  - **Parallel Group**: Wave 3 (with Tasks 13-15, 17-18)
  - **Blocks**: Task 17 (validation uses integration tests)
  - **Blocked By**: Tasks 13-15 (need core functionality and tests)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - pdifFinder command-line interface usage
  - Expected output formats from README

  **WHY Each Reference Matters**:
  - Need to match expected CLI behavior and output formats

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Integration tests pass for all input types
  - [ ] `python -m pytest tests/integration/` passes
  - [ ] End-to-end pipeline works correctly

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Run integration tests
    Tool: Bash
    Preconditions: Integration tests created
    Steps:
      1. Run `python -m pytest tests/integration/ -v`
      2. Check all tests pass
    Expected Result: All integration tests pass
    Failure Indicators: Any integration test fails
    Evidence: .sisyphus/evidence/task-16-integration-tests.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-16-integration-tests.txt`: Integration test output

  **Commit**: YES
  - Message: `test(pdifFinder): add integration tests`
  - Files: `tests/integration/`
  - Pre-commit: `python -m pytest tests/integration/`

- [ ] 17. Validate against real pdif sequences

  **What to do**:
  - Test algorithm against real pdif sequences from `data/pdifdatabase.fasta`
  - Verify algorithm detects known pdif sites correctly
  - Compare results with expected annotations (if available)
  - Test edge cases and boundary conditions
  - Ensure no false positives/negatives on known data

  **Must NOT do**:
  - Do not modify reference database
  - Do not assume perfect accuracy (allow for mismatch thresholds)

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
    - Reason: Validation requires careful comparison with ground truth
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - All skills: Not needed for validation

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 13-16, 18 in Wave 3)
  - **Parallel Group**: Wave 3 (with Tasks 13-16, 18)
  - **Blocks**: Task 18 (documentation may include validation results)
  - **Blocked By**: Tasks 13-16 (need working algorithm and tests)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `data/pdifdatabase.fasta` - Reference pdif sequences
  - `data/redundant.seed.fa` - Seed sequences

  **WHY Each Reference Matters**:
  - Real pdif sequences provide ground truth for validation

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Algorithm detects >90% of known pdif sites in reference database
  - [ ] False positive rate <5% on negative control sequences
  - [ ] `python -m pytest tests/validation/` passes

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Validate against reference database
    Tool: Bash (Python script)
    Preconditions: Validation tests created
    Steps:
      1. Run `python tests/validate_against_reference.py`
      2. Check output contains "Validation PASS: >90% detection, <5% false positives"
    Expected Result: Algorithm meets accuracy thresholds
    Failure Indicators: Accuracy below thresholds
    Evidence: .sisyphus/evidence/task-17-validation.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-17-validation.txt`: Validation results

  **Commit**: YES
  - Message: `test(pdifFinder): validate against real pdif sequences`
  - Files: `tests/validation/`
  - Pre-commit: `python -m pytest tests/validation/`

- [ ] 18. Update documentation if needed

  **What to do**:
  - Review and update documentation based on changes
  - Update README if algorithm behavior changes significantly
  - Add documentation for new features (whole plasmid scanning)
  - Ensure usage examples are still accurate
  - Update inline code comments if needed

  **Must NOT do**:
  - Do not create extensive documentation for minor changes
  - Do not remove existing useful documentation

  **Recommended Agent Profile**:
  - **Category**: `writing`
    - Reason: Documentation writing requires clear communication
  - **Skills**: None needed
  - **Skills Evaluated but Omitted**:
    - All skills: Not needed for documentation

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 13-17 in Wave 3)
  - **Parallel Group**: Wave 3 (with Tasks 13-17)
  - **Blocks**: None (final task)
  - **Blocked By**: Tasks 13-17 (need to know what changed)

  **References** (CRITICAL - Be Exhaustive):

  **Pattern References** (existing code to follow):
  - `README.md` - Existing documentation
  - `docs/` directory if exists

  **WHY Each Reference Matters**:
  - Need to maintain consistent documentation style

  **Acceptance Criteria**:

  **If TDD (tests enabled):**
  - [ ] Documentation updated if needed
  - [ ] README accurately reflects algorithm capabilities
  - [ ] No outdated information

  **QA Scenarios (MANDATORY - task is INCOMPLETE without these):**

  ```
  Scenario: Check documentation updates
    Tool: Bash (grep)
    Preconditions: Documentation reviewed/updated
    Steps:
      1. Run `grep -i "whole plasmid\|entire sequence" README.md`
      2. Check if documentation mentions scanning improvements
    Expected Result: Documentation reflects changes (if changes made)
    Failure Indicators: Documentation contradicts implementation
    Evidence: .sisyphus/evidence/task-18-docs-check.txt
  ```

  **Evidence to Capture**:
  - [ ] `task-18-docs-check.txt`: Documentation check output

  **Commit**: YES
  - Message: `docs: update documentation for algorithm improvements`
  - Files: `README.md`, `docs/`
  - Pre-commit: None

---

## Final Verification Wave (MANDATORY — after ALL implementation tasks)

> 4 review agents run in PARALLEL. ALL must APPROVE. Present consolidated results to user and get explicit "okay" before completing.

- [ ] F1. **Plan Compliance Audit** — `oracle`
  Read the plan end-to-end. For each "Must Have": verify implementation exists (read file, curl endpoint, run command). For each "Must NOT Have": search codebase for forbidden patterns — reject with file:line if found. Check evidence files exist in .sisyphus/evidence/. Compare deliverables against plan.
  Output: `Must Have [N/N] | Must NOT Have [N/N] | Tasks [N/N] | VERDICT: APPROVE/REJECT`

- [ ] F2. **Code Quality Review** — `unspecified-high`
  Run `python -m pytest tests/` + check for `as any`/`@ts-ignore`, empty catches, console.log in prod, commented-out code, unused imports. Check AI slop: excessive comments, over-abstraction, generic names (data/result/item/temp).
  Output: `Tests [N pass/N fail] | Files [N clean/N issues] | VERDICT`

- [ ] F3. **Real Manual QA** — `unspecified-high`
  Execute EVERY QA scenario from EVERY task — follow exact steps, capture evidence. Test cross-task integration. Test edge cases: empty sequences, sequences with no pdif sites, sequences with multiple pdif sites.
  Output: `Scenarios [N/N pass] | Integration [N/N] | Edge Cases [N tested] | VERDICT`

- [ ] F4. **Scope Fidelity Check** — `deep`
  For each task: read "What to do", read actual diff (git log/diff). Verify 1:1 — everything in spec was built (no missing), nothing beyond spec was built (no creep). Check "Must NOT do" compliance.
  Output: `Tasks [N/N compliant] | Contamination [CLEAN/N issues] | Unaccounted [CLEAN/N files] | VERDICT`

---

## Commit Strategy

- **1**: `fix(pdifFinder): repair XerC/XerD scanning algorithm` - pdifFinder.py, tests/
- **2**: `test(pdifFinder): add comprehensive test suite` - tests/, test_data/

---

## Success Criteria

### Verification Commands
```bash
# Run all tests
python -m pytest tests/ -v
# Expected: All tests pass (0 failures)

# Run pdifFinder on synthetic test sequence
pdifFinder -i tests/test_data/positive_test.fasta -o test_output
# Expected: pdif_site.txt contains exactly 3 pdif sites at positions 100-127, 500-527, 1000-1027

# Validate algorithm scans entire sequence
python tests/validate_scanning.py
# Expected: "All positions scanned: PASS"
```

### Final Checklist
- [ ] All "Must Have" present
- [ ] All "Must NOT Have" absent
- [ ] All tests pass
- [ ] Algorithm scans entire plasmids
- [ ] Critical bugs fixed
- [ ] Test framework established
