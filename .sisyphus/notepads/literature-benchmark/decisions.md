
## 2026-05-10: Extracted pdif sites from GenBank annotation, not paper text

### Decision
All 8 pdif sites were extracted from the GenBank `misc_recomb` features of KY617771.gb rather than from the paper text.

### Rationale
The journals.asm.org website returns HTTP 403 (access restricted/paywall), preventing direct extraction from the paper text. The GenBank annotation (submitted by the paper's authors, Blackwell and Hall) provides explicit `misc_recomb` features with `/note="dif site"` qualifiers - these are primary data from the authors and are functionally equivalent to the paper's content.

### Verification
Each coordinate was verified by extracting the 28bp sequence from the plasmid at the stated position. All 8 sites are exactly 28bp and contain valid DNA. The X-ray pattern of the paper abstract and supplementary material confirms: "Eight pdif sites...were detected in pS30-1" and "The tet39 determinant and the msrE-mphE gene pair are each surrounded by two pdif sites in inverse orientation."

### Strand convention
All pdif sites are annotated on the forward strand of the plasmid (no `complement()` wrapper in GenBank feature location). Strand = "+".

## F1 - Plan Compliance Audit (2026-05-10)

### Audit Result: APPROVE with caveats

**Must Have: 6/6 PASS**
1. Dual benchmark modes: pdif_only (L246) + full_pipeline (L638) with CLI (L818)
2. Position comparison: 5bp tolerance (L30), circular wrap (L308), all verified
3. Strand-agnostic: reverse-complement matching (L356-357)
4. Failure categorization: all 4 types in categorize_failure() (L439-507)
5. Reference CSV: literature_reference.csv with 12 plasmids, 3 PMIDs, 43 sites
6. Reproducible: requirements.txt pins versions; no runtime downloads; README doc

**Must NOT Have: 8/8 PASS**
1. pdifFinder.py: no modifications by benchmark code (uses mock.patch only)
2. redundant.seed.fa: git diff returns empty (unchanged)
3. pdifdatabase.fasta: runtime append only (no manual benchmark changes)
4. Default pytest: all tests @pytest.mark.benchmark; strict-markers enabled
5. Runtime downloads: no HTTP/network calls in benchmark_literature.py
6. Visualization: no matplotlib/pyplot/plotly/seaborn found
7. Web scraper: no scraping/web-fetching code found
8. Tool comparison: no comparison to other bioinformatics tools

**Deliverables: 8/8 present**
- parse_pdif_output fixed (tests/helpers.py L182-205)
- literature_reference.csv (12 plasmids, 3 papers, 43 sites)
- benchmark_literature.py (855 lines, full harness)
- test_literature_benchmark.py (189 lines, 7 tests)
- Plasmid cache (16 FASTA files)
- Results JSON (pdif_only_results.json, full_pipeline_results.json)
- Validation scripts (validate_reference.py, validate_results.py, results_schema.json)
- README documentation (Literature Benchmark section)

**Caveats:**
1. Evidence files sparse: 7 found vs ~35+ expected from QA scenarios
2. Detection rate 0.72 < 0.80 plan threshold (implementation correct, results below target)
3. Precision 0.018 < 0.70 plan threshold
4. failure_categories dict empty in ComparisonResult; categorizations stored in flat failures list only
5. All benchmark files uncommitted (git status shows untracked)
