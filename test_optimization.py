#!/usr/bin/env python3
"""
Test that optimized scanning algorithm produces identical results to original.
"""
import sys
sys.path.insert(0, '.')
import os
import tempfile
import shutil
from Bio import SeqIO
import PdifFinder.pdifFinder as pf

def load_seeds(pdif_db='PdifFinder/data/redundant.seed.fa'):
    """Load seed list as in findPdif."""
    seedList = []
    for rec in SeqIO.parse(pdif_db, 'fasta'):
        seq1 = str(rec.seq)
        XerC = seq1[0:11]
        XerD = seq1[17:28]
        seedList.append(XerC)
        seedList.append(XerD)
        seedList.append(XerD)
        seedList.append(XerC)
    return seedList

def create_test_sequence(length=200, pdif_site=None, position=100):
    """Create a random DNA sequence with optional pdif site."""
    import random
    random.seed(42)
    sequence = ''.join(random.choice('ACGT') for _ in range(length))
    if pdif_site:
        if position + len(pdif_site) > len(sequence):
            raise ValueError('Position out of range')
        seq_list = list(sequence)
        seq_list[position:position+len(pdif_site)] = list(pdif_site)
        sequence = ''.join(seq_list)
    return sequence

def run_function(func, seedList, seq, start1, end1):
    """Run a scanning function and return sorted positions."""
    outdir = tempfile.mkdtemp(prefix='pdif_test_')
    os.makedirs(outdir + '/tmp/seedSearch/0', exist_ok=True)
    num4 = 0
    pos0 = 0
    name = 1
    func(num4, pos0, outdir, name, seedList, seq, start1, end1)
    outfile = outdir + '/tmp/seedSearch/0/1'
    positions = []
    if os.path.exists(outfile):
        with open(outfile, 'r') as f:
            for line in f:
                positions.append(int(line.strip()))
    shutil.rmtree(outdir)
    return sorted(positions)

def test_single_seed_pair():
    """Test with a single seed pair."""
    seedList = load_seeds()
    # Use first seed pair
    seq = create_test_sequence(500)
    start1 = 1
    end1 = 1
    pos_original = run_function(pf.findMatchFragmentThreadOriginal, seedList, seq, start1, end1)
    pos_optimized = run_function(pf.findMatchFragmentThread, seedList, seq, start1, end1)
    assert pos_original == pos_optimized, f"Mismatch: original {pos_original}, optimized {pos_optimized}"
    print(f"✓ Single seed pair test passed (found {len(pos_original)} positions)")

def test_all_seeds_small_sequence():
    """Test with all seed pairs on a small sequence."""
    seedList = load_seeds()
    seq = create_test_sequence(300)
    start1 = 1
    end1 = len(seedList) // 2
    pos_original = run_function(pf.findMatchFragmentThreadOriginal, seedList, seq, start1, end1)
    pos_optimized = run_function(pf.findMatchFragmentThread, seedList, seq, start1, end1)
    assert pos_original == pos_optimized, f"Mismatch: original {len(pos_original)} positions, optimized {len(pos_optimized)} positions"
    print(f"✓ All seeds test passed (found {len(pos_original)} positions)")

def test_with_pdif_site():
    """Test with a known pdif site inserted."""
    seedList = load_seeds()
    # Take the first seed sequence as pdif site
    first_rec = next(SeqIO.parse('PdifFinder/data/redundant.seed.fa', 'fasta'))
    pdif_site = str(first_rec.seq)  # 28bp
    seq = create_test_sequence(1000, pdif_site, position=200)
    start1 = 1
    end1 = len(seedList) // 2
    pos_original = run_function(pf.findMatchFragmentThreadOriginal, seedList, seq, start1, end1)
    pos_optimized = run_function(pf.findMatchFragmentThread, seedList, seq, start1, end1)
    assert pos_original == pos_optimized, f"Mismatch with pdif site"
    # Expect at least position 200 found
    assert 200 in pos_original, f"Expected pdif site at position 200 not found in original"
    assert 200 in pos_optimized, f"Expected pdif site at position 200 not found in optimized"
    print(f"✓ Pdif site detection test passed (found {len(pos_original)} positions including position 200)")

def test_edge_cases():
    """Test edge cases: short sequence, no matches."""
    seedList = load_seeds()
    # Sequence shorter than 28bp
    seq = create_test_sequence(20)
    start1 = 1
    end1 = len(seedList) // 2
    pos_original = run_function(pf.findMatchFragmentThreadOriginal, seedList, seq, start1, end1)
    pos_optimized = run_function(pf.findMatchFragmentThread, seedList, seq, start1, end1)
    assert pos_original == pos_optimized == [], f"Should be empty for short sequence"
    print("✓ Short sequence test passed")
    # Sequence with no pdif sites (random)
    seq = create_test_sequence(500)
    pos_original = run_function(pf.findMatchFragmentThreadOriginal, seedList, seq, start1, end1)
    pos_optimized = run_function(pf.findMatchFragmentThread, seedList, seq, start1, end1)
    assert pos_original == pos_optimized, f"Mismatch for random sequence"
    print("✓ Random sequence test passed")

if __name__ == '__main__':
    print("Testing optimized scanning algorithm...")
    test_single_seed_pair()
    test_all_seeds_small_sequence()
    test_with_pdif_site()
    test_edge_cases()
    print("\nAll tests passed! Optimized algorithm matches original.")