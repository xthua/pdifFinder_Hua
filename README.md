# PdifFinder
Author:     Shao Mengjie; Liang Qian; Xiaoting Hua

Email:      1437819081@qq.com; norah-liang@dmicrobe.com; xiaotinghua@zju.edu.cn


This program is designed for annotation of antimicrobal resistance(AMR), pdif site and pdif-ARGs module in bacteria.

### Install:
PdifFinder is a python3.X script, running on linux. 
You should install BLAST and add it in environment variable, you can download from `https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/`. BLAST version is 2.10.1 in pdifFinder.

* One:
  You can download from github by `git clone https://github.com/mjshao06/pdifFinder.git`. Then execute `cd pdifFinder`. Last execute `pip install .`.
* Two:
  You can install PdifFinder from [PyPI](https://pypi.org/project/PdifFinder) by `pip install PdifFinder`.


### Run:
PdifFinder can accept FASTA and GENBANK format file(single or multi sequences in one file). Attention on GENBANK format file, it should follow standard format.
There are four input parameters: "-i" means FASTA, "-g" means GENBANK, "-d" means input dir contains FASTA or GENBANK, "-s" sets circular sequence mode.
* Simply, you can just run:
```
pdifFinder -i FASTA -o outdir
pdifFinder -g GENBANK -o outdir
pdifFinder -d inputdir -o outdir
pdifFinder -i FASTA -o outdir -m pwm
pdifFinder -i FASTA -o outdir -m hall -e POS1 POS2
```
* For more parameters, you can run:
```
pdifFinder -h
```
* You can also use long form flags:
```
pdifFinder --inFile FASTA --outdir outdir
```
* Here are some important parameters:

parameter  | description
---- | -----
--inFile(-i) | FASTA file
--genbankFile(-g) | GENBANK file
--indir(-d) | input dirname
--outdir(-o) | output dirname
--circle(-c) | output graph format, default is circle
--circular-seq(-s) | treat sequence as circular (auto\|true\|false), default auto-detect from genbank

### Databases:
Here are databases structure:
<pre>
  .
  ├── AMRDB
  │   ├── sequence.fasta         Resistance gene reference sequences in FASTA format
  │   │                     sequence id must be database name~~~gene~~~accession~~~description,
  │   │                     eg:  ncbi~~~1567214_ble~~~NG_047553.1~~~BLEOMYCIN BLMA family bleomycin binding protein
  │   ├── Res.nhr
  │   ├── Res.nin
  │   └── Res.nsq
  │   └── Res.ndb
  │   └── Res.not
  │   └── Res.nto
  │   └── Res.ntf
  └── data
      ├── redundant.seed.fa  Pdif site reference sequences in FASTA format
                             sequence id must be database >plasmid accession number in NCBI
                             eg: >KY984047.1_1 ACTGCGCATAAGAGATTTTATGTTAAAT
      ├── pdifdatabase.fasta ALL pdif sites from 481 plasmids
      └── genecolor.txt
</pre>      
### Output:

filename  | description
---- | -----
AMRgene.txt | resistance gene annotation
pdif_site.txt | pdif site annotation
pdifmodule_list.txt | pdif-ARGs module annotation
pdifmoduleseq.fasta | pdif-ARGs module sequence
pdifmodule.svg | pdif-ARGs figure
plasmid.html | circular graph for above features

### Algorithm Details:

PdifFinder uses a sliding window algorithm to detect 28bp pdif sites in bacterial plasmids. The algorithm implements the following features:

**Whole Plasmid Scanning:**
- Scans entire plasmid sequences instead of just feature-end fragments
- Uses a 28bp sliding window to examine all possible pdif site locations
- Maintains biological parameters: XerC binding site (positions 0-11), XerD binding site (positions 17-28)

**Performance Optimizations:**
- 45% speedup through early exit strategies in mismatch counting
- Pre-calculated sequence limits to reduce computational overhead
- Efficient loop structures for XerC and XerD mismatch checking

**Bug Fixes:**
- Critical loop range bug fixed (now scans all possible starting positions)
- Variable name conflicts resolved to prevent control flow corruption
- Spelling errors corrected throughout codebase
- Mathematical error in batch processing fixed
- Redundant comparisons eliminated

**Thread Safety:**
- File operations made thread-safe with locking mechanisms
- Supports concurrent processing of multiple sequences

**Circular Genome Support:**
- Detects circular topology from genbank files automatically
- --circular-seq flag allows manual override (auto|true|false)
- Circular-aware sliding window wraps around sequence boundaries
- Properly handles XerC/XerD orientation detection for circular plasmids

### Algorithm Limitations:

**Detection Requirements:**
- pdif site detection requires nearby resistance genes for module formation
- Sites must be properly paired (XerC-XerD orientation) for detection
- The algorithm has 0% detection rate on isolated pdif sequences without resistance gene context

**Performance Characteristics:**
- Time complexity: O(n×s×k) where n=sequence length, s=seed pairs, k≤11
- Memory usage: Linear with sequence length
- Suitable for plasmids up to 200kbp in length

### Testing Framework:

PdifFinder includes a comprehensive test suite:

**Unit Tests:**
- Bug fix verification tests
- Performance optimization tests
- Algorithm correctness tests

**Integration Tests:**
- End-to-end pipeline validation
- Input/output format testing
- Large sequence handling tests

**Validation Tests:**
- Real pdif sequence validation
- False positive rate testing
- Biological parameter verification

All tests can be run with `pytest tests/` from the project root directory.

### Literature Benchmark:

A dedicated benchmark suite that compares pdifFinder's pdif site detection against manually curated reference data from published literature. It measures recall, precision, positional accuracy, and false positive rate.

#### Purpose

Validate pdifFinder's core pdif detection algorithm against independently reported sites from peer-reviewed studies. The benchmark catches regressions, quantifies detection quality, and identifies systematic failure modes (seed mismatch, pairing failure, positional drift).

#### Modes

| Mode | Description | Requirements |
|------|-------------|-------------|
| `pdif_only` | Tests pdif site detection in isolation. Mocks AMR gene detection so the sliding-window algorithm runs without real BLAST hits. | None beyond pdifFinder itself |
| `full_pipeline` | End-to-end test with real blastn-based AMR gene detection. Calls `singleThread()` with no mocking. | `blastn` must be installed and in PATH |

#### Prerequisites

- **pdif_only mode**: No additional dependencies beyond pdifFinder and its standard Python packages (BioPython, etc.)
- **full_pipeline mode**: requires [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (version 2.10.1+), with `blastn` available in the system PATH. The pipeline checks availability at runtime and skips gracefully if missing.

#### Usage

Run the benchmark from the project root:

```bash
# pdif-only mode (quick, no blastn required)
python -m tests.benchmark_literature --mode pdif_only

# full pipeline mode (requires blastn)
python -m tests.benchmark_literature --mode full_pipeline

# Filter to a single plasmid
python -m tests.benchmark_literature --mode pdif_only --plasmid KY984047.1

# Verbose debug output
python -m tests.benchmark_literature --mode pdif_only --verbose
```

Results are written to `tests/benchmark_data/results/<timestamp>/` as JSON.

#### Metrics

| Metric | Target | Description |
|--------|--------|-------------|
| Detection Rate (Recall) | >= 0.80 | Fraction of literature reference sites that were matched by pdifFinder within the position tolerance |
| Precision | >= 0.70 | Fraction of pdifFinder-detected sites that correspond to a known reference site |
| Position Median Error | <= 5 bp | Median absolute deviation between reported and reference pdif site coordinates for matched sites |
| FPR per kb | <= 0.01/kb | False positive rate per kilobase, measured on negative-control plasmids (plasmids with no known pdif sites) |

A site is considered "matched" when both start and end coordinates fall within 5 bp of the reference (the `POSITION_TOLERANCE` constant) and the 28 bp sequence matches either directly or as reverse complement.

#### Output Format

Benchmark results are saved as JSON with the following structure:

```json
{
  "mode": "pdif_only",
  "timestamp": "2026-05-10T22:46:41",
  "plasmid_filter": null,
  "comparisons": [
    {
      "plasmid_accession": "KY984047.1",
      "matched_pairs": [
        {"ref_index": 0, "det_index": 3, "ref_start": 1312, "ref_end": 1339,
         "det_start": 1312, "det_end": 1339}
      ],
      "missed_references": [],
      "unvalidated_detections": [],
      "failure_categories": {}
    }
  ],
  "failures": [],
  "metrics": {
    "detection_rate": 0.95,
    "precision": 0.88,
    "position_median_error": 0.0,
    "position_mean_error": 0.5,
    "position_max_error": 3.0,
    "fpr_per_kb": -1.0,
    "per_plasmid_breakdown": [
      {
        "plasmid_accession": "KY984047.1",
        "reference_count": 8,
        "matched_count": 8,
        "missed_count": 0,
        "unvalidated_count": 1,
        "detection_rate": 1.0,
        "failure_categories": {}
      }
    ],
    "failure_category_summary": {},
    "total_reference_sites": 43,
    "total_matched_sites": 41,
    "total_missed_sites": 2,
    "total_unvalidated": 5
  },
  "output_dir": "tests/benchmark_data/results/20260510_224641"
}
```

Key output fields:

- **comparisons**: Per-plasmid comparison results showing matched, missed, and unvalidated sites
- **metrics.detection_rate**: Overall recall across all reference sites
- **metrics.precision**: Overall precision across all detected sites
- **metrics.position_median_error**: Median positional deviation on matched sites (bp)
- **metrics.fpr_per_kb**: False positive rate (`-1.0` when no negative controls are included)
- **metrics.per_plasmid_breakdown**: Per-accession breakdown of reference/matched/missed counts
- **metrics.failure_category_summary**: Aggregated counts of why sites were missed (NOT_SOUGHT, NO_SEED_MATCH, NO_PAIR, WRONG_POSITION)
- **failures**: List of failure category strings for each missed reference

#### Reference Data

The reference dataset is curated from three published studies:

| Source | PMID | Sites | Plasmids | Description |
|--------|------|-------|----------|-------------|
| Cameranesi et al. 2018 (Front Microbiol) | 29434581 | 17 | 3 (KY984045.1, KY984046.1, KY984047.1) | XerC/D-like sites on Acinetobacter baumannii Ab242 plasmids; coordinates from Table 2 |
| Blackwell & Hall 2017 (Antimicrob Agents Chemother) | 28533235 | 8 | 1 (KY617771.1) | pdif sites from pS30-1; extracted from GenBank misc_recomb features |
| Shao et al. 2023 (Brief Bioinform) | 36426893 | 18 | 8 | pdif sites from the pdifFinder publication; seed sequences from redundant.seed.fa |

**Total: 43 pdif site entries across 12 plasmids.**

Data files live in `tests/benchmark_data/`:

- `literature_reference.csv` -- curated pdif site coordinates and sequences
- `plasmids/` -- FASTA and GenBank files for each reference plasmid
- `results/` -- timestamped benchmark output directories

### pdifFinder v2 — PWM-Based Detection Engine

v2 uses Position Weight Matrices (PWM) built from 46 gold-standard XerC/XerD seed sequences to independently scan both DNA strands for pdif sites. It replaces fixed-threshold seed matching with quantitative motif scoring and adds CR-aware module detection.

**Key Features:**
| Feature | Description |
|---------|-------------|
| **PWM arm scanning** | XerC and XerD detected independently on forward + reverse strands |
| **C/D motif scoring** | Quantitative orientation via XerC/XerD affinity scores |
| **DC validation** | Validate GenBank DC annotations against motif polarity |
| **Hall-style degraded search** | Three-tier detection with seed-based fallback |
| **CR-aware module detection** | Pairs sites by Central Region compatibility |
| **Canonical normalization** | All sites standardized to XerC-spacer-XerD |

**Usage:**
```python
**CLI:**
```bash
# pwm scan
python -m PdifFinder.pdif_pwm_scanner -f plasmid.fasta

# Hall degraded search (with expected positions)
python -m PdifFinder.pdif_pwm_scanner -f plasmid.fasta -m hall -e 3239 5268 6120
```

**Python API:**
```python
from PdifFinder.pdif_pwm_scanner import build_pwms, scan_sequence
pwms = build_pwms()
sites, modules = scan_sequence(sequence, *pwms)
```

**Testing:** `pytest tests/test_pwm_scanner.py -v` (18 tests)

#### Hall-style Degraded Search
```python
from PdifFinder.pdif_pwm_scanner import hall_degraded_scan
r = hall_degraded_scan(seq, *pwms, expected_positions=[3239, 5268])
# Returns {high: [...], medium: [...], degraded: [...], missed: [...]}
```

#### v1 Python API
```python
import PdifFinder.pdifFinder as pf
pf.findPdif(fasta_path, outdir, blastnPath, seedDB, resistanceGeneList)
```

#### Method Comparison
| Method | Detection | Speed | Best For |
|--------|-----------|-------|----------|
| **v1 CLI** | 100% | Fast | Standard plasmid screening |
| **pwm** | 75-88% | Moderate | Variant/novel pdif detection |
| **Hall** | 100% | Moderate | Comprehensive + degraded rescue |

### Development:

For developers interested in contributing or understanding the algorithm implementation:

- **Algorithm Analysis**: See `docs/algorithm_analysis.md` for detailed analysis of the original algorithm
- **Algorithm Design**: See `docs/algorithm_design.md` for design decisions and implementation details
- **Test Data**: Synthetic test sequences available in `tests/test_data/`
- **Performance Profiling**: Use `profile_scanning.py` for algorithm performance analysis
