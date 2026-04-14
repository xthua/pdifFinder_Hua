# Algorithm Design - Improved pdif Site Scanning

**Date:** Tue Apr 14 2026  
**Task:** 5. Design improved scanning algorithm  
**File:** `pdifFinder/pdifFinder.py`

## Executive Summary

This document outlines the design for improving the pdifFinder algorithm to:
1. Scan **entire plasmid sequences** instead of just feature-end fragments
2. Fix the critical `range(1)` bug that limits scanning to first position
3. Maintain biological correctness (28bp pattern, XerC[0:11], XerD[17:28], mismatch thresholds)
4. Optimize performance for whole plasmid scanning
5. Preserve backward compatibility where possible

## 1. Current Limitations (From Analysis)

### 1.1 Algorithm Flow
```
findPdif() → findFeatureEndSeq() → fragments → findMatchFragmentThread()
```

### 1.2 Key Issues
1. **Fragment limitation**: Only scans 28bp fragments around `TAA`, `TGT`, `CAT` patterns
2. **Critical bug**: `for i in range(1):` only checks first position of each fragment
3. **Performance**: Current naive scanning is O(n×m) where n=sequence length, m=seed patterns
4. **Thread safety**: Multiple threads write to same files without synchronization

## 2. Design Goals

### 2.1 Primary Goals
1. **Whole plasmid scanning**: Detect pdif sites anywhere in plasmid sequences
2. **Correctness**: Fix `range(1)` bug to scan all positions
3. **Performance**: Efficient scanning of long sequences (up to 1Mbp)
4. **Accuracy**: Maintain biological parameters (28bp, mismatch thresholds)

### 2.2 Secondary Goals
1. **Backward compatibility**: Maintain existing output formats
2. **Thread safety**: Proper synchronization for concurrent execution
3. **Code quality**: Fix variable conflicts, spelling errors, redundant code
4. **Testability**: Enable comprehensive testing with synthetic sequences

## 3. Proposed Architecture

### 3.1 High-Level Overview
```
Input: FASTA/GENBANK file
↓
Parse sequences (Bio.SeqIO)
↓
For each sequence:
  - Extract full sequence (not fragments)
  - Scan for pdif sites using optimized algorithm
  - Apply XerC/XerD mismatch thresholds
  - Output results in existing format
```

### 3.2 Key Design Decisions

#### 3.2.1 Whole Plasmid Scanning
- **Replace `findFeatureEndSeq()`** with direct sequence extraction
- **Scan entire sequence** using sliding window of 28bp
- **Maintain XerC/XerD pattern matching** with existing seed database

#### 3.2.2 Performance Optimization
Two approaches considered:

**Option A: Enhanced Sliding Window**
- Simple sliding window across sequence
- Early exit when mismatch thresholds exceeded
- Optimized for 28bp window size

**Option B: Boyer-Moore Adaptation**
- Use existing BM implementation for pattern matching
- Match 11bp XerC/XerD regions separately
- Combine with mismatch tolerance

**Recommended: Hybrid Approach**
- Use sliding window for simplicity and maintainability
- Apply early exit optimization
- Consider BM for exact pattern matching in future optimization

#### 3.2.3 Thread Safety
- **File-level locking** using `threading.Lock`
- **Separate output files** per thread with post-merge
- **Queue-based writing** for concurrent operations

## 4. Detailed Algorithm Design

### 4.1 Sequence Processing Pipeline

```python
def process_sequence(sequence, seed_list, outdir):
    """
    Process a single sequence for pdif sites
    """
    results = []
    
    # Whole plasmid scanning (replace fragment extraction)
    for start_pos in range(len(sequence) - 27):  # 28bp window
        window = sequence[start_pos:start_pos + 28]
        
        # Check against all seed pairs
        for seed_c, seed_d in seed_pairs:
            if is_pdif_site(window, seed_c, seed_d):
                results.append({
                    'start': start_pos + 1,  # 1-based
                    'end': start_pos + 28,
                    'xerc_seq': window[0:11],
                    'spacer': window[11:17],
                    'xerd_seq': window[17:28]
                })
                break  # Found pdif site, move to next position
    
    return results
```

### 4.2 pdif Site Detection Algorithm

```python
def is_pdif_site(window, seed_c, seed_d, max_mismatch_c=3, max_mismatch_d=2):
    """
    Determine if a 28bp window contains a pdif site
    matching the given XerC/XerD seed patterns
    """
    xerc_window = window[0:11]    # Positions 0-10
    xerd_window = window[17:28]   # Positions 17-27
    
    # Check XerC region
    mismatch_c = count_mismatches(xerc_window, seed_c)
    if mismatch_c > max_mismatch_c:
        return False
    
    # Check XerD region  
    mismatch_d = count_mismatches(xerd_window, seed_d)
    if mismatch_d > max_mismatch_d:
        return False
    
    return True


def count_mismatches(seq1, seq2):
    """
    Count mismatches between two sequences of equal length
    """
    return sum(1 for a, b in zip(seq1, seq2) if a != b)
```

### 4.3 Seed Database Processing

```python
def load_seed_database(pdif_db):
    """
    Load seed patterns from pdif database
    Returns list of (xerc_seq, xerd_seq) pairs
    """
    seed_pairs = []
    for rec in SeqIO.parse(pdif_db, 'fasta'):
        seq = str(rec.seq)
        xerc = seq[0:11]   # First 11bp: XerC binding site
        xerd = seq[17:28]  # Last 11bp: XerD binding site
        seed_pairs.append((xerc, xerd))
    
    return seed_pairs
```

### 4.4 Performance Optimizations

#### 4.4.1 Early Exit Optimization
```python
def is_pdif_site_optimized(window, seed_c, seed_d, max_mismatch_c=3, max_mismatch_d=2):
    """
    Optimized version with early exit
    """
    xerc_window = window[0:11]
    
    # Early exit for XerC (more stringent check)
    mismatch_c = 0
    for i in range(11):
        if xerc_window[i] != seed_c[i]:
            mismatch_c += 1
            if mismatch_c > max_mismatch_c:
                return False  # Early exit
    
    xerd_window = window[17:28]
    
    # Check XerD
    mismatch_d = 0
    for i in range(11):
        if xerd_window[i] != seed_d[i]:
            mismatch_d += 1
            if mismatch_d > max_mismatch_d:
                return False  # Early exit
    
    return True
```

#### 4.4.2 Sliding Window Optimization
```python
def scan_sequence_optimized(sequence, seed_pairs):
    """
    Optimized scanning with window reuse
    """
    results = []
    seq_len = len(sequence)
    
    if seq_len < 28:
        return results  # Too short for pdif site
    
    # Pre-calculate windows for each position
    for start in range(seq_len - 27):
        window = sequence[start:start + 28]
        
        # Try seed pairs in order
        for seed_c, seed_d in seed_pairs:
            if is_pdif_site_optimized(window, seed_c, seed_d):
                results.append(start)
                break  # Found pdif site
        
        # Optional: Skip ahead if no match possible
        # (Advanced optimization for future)
    
    return results
```

## 5. Integration with Existing Code

### 5.1 Minimal Changes Approach

#### 5.1.1 Modify `findPdif()` function
```python
# Current (line 220):
seq, pos0List, seqList = findFeatureEndSeq(inFile, outdir)

# Proposed:
seq = extract_full_sequence(inFile)  # New function
seqList = [seq]  # Single element list for compatibility
pos0List = [0]   # Starting position offset
```

#### 5.1.2 New Helper Function
```python
def extract_full_sequence(inFile):
    """
    Extract full sequence from input file
    Maintains compatibility with existing code
    """
    for rec in SeqIO.parse(inFile, 'fasta'):
        return str(rec.seq).upper()
    return ""
```

### 5.2 Thread Safety Implementation

```python
import threading

# Global lock for file writing
file_lock = threading.Lock()

def write_result_thread_safe(outfile, result):
    """
    Thread-safe file writing
    """
    with file_lock:
        with open(outfile, 'a') as f:
            f.write(result + '\n')
```

### 5.3 Backward Compatibility

#### 5.3.1 Output Format
Maintain existing output files:
- `pdif_site.txt`: Same format (start end XerC spacer XerD)
- `pdif_pair.txt`: Same pairwise analysis
- All visualization files unchanged

#### 5.3.2 API Compatibility
- Command-line interface unchanged
- Input file formats unchanged (FASTA/GENBANK)
- Output directory structure unchanged

## 6. Performance Analysis

### 6.1 Complexity Analysis

| Algorithm | Time Complexity | Space Complexity | Notes |
|-----------|-----------------|------------------|-------|
| Current (fragments) | O(f × s × p) | O(1) | f=fragments, s=seed pairs, p=positions |
| Proposed (whole) | O(n × s × 11) | O(1) | n=sequence length, s=seed pairs |
| Optimized (early exit) | O(n × s × k) | O(1) | k≤11, early exit reduces comparisons |

Where:
- n: Sequence length (up to 1,000,000 bp)
- s: Seed pairs (~100-200 from database)
- f: Fragments (typically 10-100)

### 6.2 Expected Performance

| Sequence Length | Current Algorithm | Proposed Algorithm | Improvement |
|-----------------|------------------|-------------------|-------------|
| 10,000 bp | ~100 fragments × scanning | 10,000 positions × scanning | ~100x slower* |
| 100,000 bp | ~1,000 fragments × scanning | 100,000 positions × scanning | ~100x slower* |
| 1,000,000 bp | ~10,000 fragments × scanning | 1,000,000 positions × scanning | ~100x slower* |

*Note: Slower in raw operations but detects ALL pdif sites, not just fragment-localized ones.

### 6.3 Optimization Opportunities

1. **Seed clustering**: Group similar seed patterns to reduce comparisons
2. **Bitmask encoding**: Encode DNA as 2-bit values for faster comparison
3. **SIMD operations**: Use vectorized operations for mismatch counting
4. **Parallel scanning**: Divide sequence across threads (already implemented)

## 7. Implementation Plan

### 7.1 Phase 1: Core Bug Fixes (Wave 2)
1. Fix `range(1)` bug in `findMatchFragmentThread`
2. Resolve variable name conflicts
3. Fix spelling errors (`maxMistach` → `maxMismatch`)
4. Correct mathematical error (line 474)
5. Remove redundant comparisons
6. Implement thread-safe file operations

### 7.2 Phase 2: Whole Plasmid Scanning (Wave 3)
1. Replace `findFeatureEndSeq` with full sequence extraction
2. Update `findMatchFragmentThread` for whole sequence scanning
3. Implement early exit optimization
4. Add comprehensive test coverage
5. Performance benchmarking

### 7.3 Phase 3: Advanced Optimizations (Optional)
1. Implement Boyer-Moore based scanning
2. Add seed pattern clustering
3. Bitmask encoding for faster comparisons
4. Profile and optimize hotspots

## 8. Testing Strategy

### 8.1 Unit Tests
- Test `is_pdif_site` with known pdif sites
- Test edge cases (short sequences, no pdif sites)
- Test mismatch threshold logic

### 8.2 Integration Tests
- End-to-end testing with synthetic sequences
- Compare results with expected pdif site positions
- Verify output file formats

### 8.3 Performance Tests
- Benchmark scanning speed on various sequence lengths
- Compare before/after optimization
- Measure memory usage

### 8.4 Biological Validation
- Test with known pdif-containing plasmids
- Verify sensitivity and specificity
- Compare with existing biological knowledge

## 9. Risk Mitigation

### 9.1 Technical Risks
1. **Performance degradation**: Addressed by optimizations and benchmarking
2. **False positives/negatives**: Addressed by comprehensive testing
3. **Backward compatibility**: Addressed by maintaining output formats

### 9.2 Implementation Risks
1. **Complexity**: Incremental implementation with testing at each step
2. **Thread safety**: Use proven locking patterns
3. **Code quality**: Follow existing patterns and add documentation

### 9.3 Biological Risks
1. **Parameter changes**: Maintain existing mismatch thresholds (3 for XerC, 2 for XerD)
2. **Pattern recognition**: Use same 28bp pattern definition
3. **Validation**: Test with known biological examples

## 10. Conclusion

The proposed design addresses the fundamental limitations of the current pdifFinder algorithm while maintaining biological correctness and backward compatibility. The key improvements are:

1. **Whole plasmid scanning** instead of fragment-limited scanning
2. **Correct scanning range** (fixing the `range(1)` bug)
3. **Performance optimizations** for efficient large sequence processing
4. **Thread-safe implementation** for reliable concurrent execution

This design provides a clear roadmap for implementation in Tasks 7-18, with measurable success criteria and comprehensive testing strategy.

---

## Appendix A: Pseudocode Summary

### A.1 Main Algorithm
```
For each input sequence:
  1. Extract full sequence
  2. Load seed patterns from database
  3. For position = 0 to length-28:
      a. Extract 28bp window
      b. For each seed pair (XerC, XerD):
          i. Count mismatches in XerC region (0-10)
          ii. If mismatches > 3, continue to next seed
          iii. Count mismatches in XerD region (17-27)
          iv. If mismatches > 2, continue to next seed
          v. Record pdif site position
  4. Output results
```

### A.2 Optimized Version
```
For each input sequence:
  1. Extract full sequence
  2. Load seed patterns
  3. For position = 0 to length-28:
      a. window = sequence[position:position+28]
      b. xerc = window[0:11], xerd = window[17:28]
      c. For each seed pair (seed_c, seed_d):
          i. mismatches_c = 0
          ii. For i = 0 to 10:
              - If xerc[i] != seed_c[i]: mismatches_c++
              - If mismatches_c > 3: break (early exit)
          iii. If mismatches_c <= 3:
              - mismatches_d = 0
              - For i = 0 to 10:
                  - If xerd[i] != seed_d[i]: mismatches_d++
                  - If mismatches_d > 2: break (early exit)
              - If mismatches_d <= 2: record position
  4. Output results
```

## Appendix B: Compatibility Matrix

| Feature | Current Algorithm | Proposed Algorithm | Change Required |
|---------|------------------|-------------------|----------------|
| Input formats | FASTA, GENBANK | FASTA, GENBANK | None |
| Output files | pdif_site.txt, etc. | Same files | None (format preserved) |
| Scanning scope | Feature-end fragments | Whole plasmids | Major change |
| Performance | Fast (limited scope) | Slower (full scope) | Optimization needed |
| Accuracy | Misses non-fragment sites | Detects all sites | Improvement |
| Threading | Unsafe file writes | Thread-safe | Moderate change |
| Parameters | 28bp, XerC/XerD thresholds | Same parameters | None |

## Appendix C: Reference Implementation Snippets

### C.1 Sequence Extraction Replacement
```python
# Replace in findPdif() function:
# Current:
seq, pos0List, seqList = findFeatureEndSeq(inFile, outdir)

# New:
def extract_full_sequence(inFile):
    for rec in SeqIO.parse(inFile, 'fasta'):
        return str(rec.seq).upper()
    return ""

seq = extract_full_sequence(inFile)
seqList = [seq]  # Single sequence for compatibility
pos0List = [0]   # Starting at position 0
```

### C.2 Fixed Scanning Loop
```python
# In findMatchFragmentThread(), line 393:
# Current bug:
for i in range(1):  # Only checks position 0

# Fixed version:
for i in range(len(seq) - 27):  # Scan all positions for 28bp window
    # Ensure boundary check:
    if i + 28 <= len(seq):
        # ... existing logic with i as starting position
```

### C.3 Thread-Safe File Writing
```python
import threading

write_lock = threading.Lock()

def write_pdif_result(outfile, position):
    with write_lock:
        with open(outfile, 'a') as f:
            f.write(f"{position}\n")
```