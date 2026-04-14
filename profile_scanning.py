#!/usr/bin/env python3
"""
Profile the scanning algorithm to identify bottlenecks.
"""
import sys
sys.path.insert(0, '.')
import os
import tempfile
import time
import cProfile
import pstats
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

def create_test_sequence(length=2000, pdif_site=None, position=100):
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

def profile_single_seed_pair(seq_length=1000, seed_index=0):
    """Profile scanning with a single seed pair."""
    seedList = load_seeds()
    # Select first seed pair (indices 0,1)
    seq = create_test_sequence(seq_length)
    outdir = tempfile.mkdtemp(prefix='pdif_profile_')
    os.makedirs(outdir + '/tmp/seedSearch/0', exist_ok=True)
    
    num4 = 0
    pos0 = 0
    name = 1
    start1 = seed_index + 1  # 1-based indexing in original
    end1 = start1
    
    print(f'Profiling sequence length {seq_length}, seed pair {seed_index}')
    profiler = cProfile.Profile()
    profiler.enable()
    pf.findMatchFragmentThread(num4, pos0, outdir, name, seedList, seq, start1, end1)
    profiler.disable()
    
    stats = pstats.Stats(profiler).sort_stats('cumulative')
    stats.print_stats(20)
    
    # Cleanup
    import shutil
    shutil.rmtree(outdir)

def profile_all_seeds(seq_length=1000):
    """Profile scanning with all seed pairs."""
    seedList = load_seeds()
    seq = create_test_sequence(seq_length)
    outdir = tempfile.mkdtemp(prefix='pdif_profile_')
    os.makedirs(outdir + '/tmp/seedSearch/0', exist_ok=True)
    
    num4 = 0
    pos0 = 0
    name = 1
    start1 = 1
    end1 = len(seedList) // 2  # seedListPairLength
    
    print(f'Profiling sequence length {seq_length}, all {end1} seed pairs')
    start_time = time.time()
    pf.findMatchFragmentThread(num4, pos0, outdir, name, seedList, seq, start1, end1)
    elapsed = time.time() - start_time
    print(f'Total time: {elapsed:.2f}s')
    
    # Cleanup
    import shutil
    shutil.rmtree(outdir)
    return elapsed

def time_varying_lengths():
    """Measure time for different sequence lengths."""
    seedList = load_seeds()
    lengths = [100, 500, 1000, 2000, 5000, 10000]
    times = []
    for L in lengths:
        seq = create_test_sequence(L)
        outdir = tempfile.mkdtemp(prefix='pdif_profile_')
        os.makedirs(outdir + '/tmp/seedSearch/0', exist_ok=True)
        num4 = 0
        pos0 = 0
        name = 1
        start1 = 1
        end1 = len(seedList) // 2
        start = time.time()
        pf.findMatchFragmentThread(num4, pos0, outdir, name, seedList, seq, start1, end1)
        elapsed = time.time() - start
        times.append(elapsed)
        print(f'Length {L}: {elapsed:.3f}s')
        import shutil
        shutil.rmtree(outdir)
    return lengths, times

if __name__ == '__main__':
    print("=== Profiling pdifFinder scanning algorithm ===")
    # Profile single seed pair with cProfile
    profile_single_seed_pair(seq_length=1000)
    # Time all seeds with varying lengths
    lengths, times = time_varying_lengths()
    print("\nLength vs Time:")
    for L, t in zip(lengths, times):
        print(f'{L}\t{t:.3f}')