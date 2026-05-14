#!/usr/bin/env python3
"""
Test whole plasmid scanning functionality.
Focus on verifying that findFeatureEndSeq returns whole plasmid sequence
and that pdifFinder can detect pdif sites throughout entire plasmids.
"""

import os
import sys
import tempfile
from pathlib import Path
import pytest

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import PdifFinder.pdifFinder as pf
from Bio import SeqIO
from tests.helpers import (
    generate_test_sequence,
    save_sequence_to_fasta,
    cleanup_temp_files,
    cleanup_temp_dirs,
)


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


def test_find_feature_end_seq_returns_whole_plasmid():
    """Test that findFeatureEndSeq returns the entire plasmid sequence."""
    test_seq = generate_test_sequence(1000, seed=42)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">test_plasmid\n{test_seq}\n")
        fasta_file = f.name
    
    outdir = None
    try:
        outdir = tempfile.mkdtemp(prefix='pdif_test_whole_plasmid_')
        seq, pos0_list, seq_list = pf.findFeatureEndSeq(fasta_file, outdir)
        
        assert seq == test_seq
        assert len(seq_list) == 1
        assert seq_list[0] == test_seq
        assert pos0_list == [0]
    finally:
        cleanup_temp_files([fasta_file])
        if outdir:
            cleanup_temp_dirs([outdir])


def test_whole_plasmid_scanning_detects_pdif_sites():
    """Test that pdifFinder can detect pdif sites anywhere in the plasmid."""
    seed_list = load_seeds()
    if len(seed_list) < 2:
        pytest.skip("No seed sequences found in database")
    
    first_rec = next(SeqIO.parse('PdifFinder/data/redundant.seed.fa', 'fasta'))
    pdif_site = str(first_rec.seq)
    background = generate_test_sequence(1000, seed=123)
    position = 200
    test_seq = background[:position] + pdif_site + background[position + len(pdif_site):]
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">test_with_pdif\n{test_seq}\n")
        fasta_file = f.name
    
    outdir = None
    try:
        outdir = tempfile.mkdtemp(prefix='pdif_test_whole_plasmid_')
        tmp_dir = Path(outdir) / 'tmp' / 'seedSearch' / '0'
        tmp_dir.mkdir(parents=True, exist_ok=True)
        outfile = tmp_dir / '1'
        
        pf.findMatchFragmentThread(
            num4=0,
            pos0=0,
            outdir=outdir,
            name='1',
            seedList=seed_list,
            seq=test_seq,
            start1=1,
            end1=len(seed_list) // 2
        )
        
        if outfile.exists():
            with open(outfile, 'r') as f:
                positions = [int(line.strip().split('|')[0]) for line in f if line.strip()]
            found = any(abs(pos - position) <= 5 for pos in positions)
            assert found, f"Expected pdif site at position {position} not found. Found positions: {positions}"
        else:
            pytest.fail(f"No output file created at {outfile}")
    finally:
        cleanup_temp_files([fasta_file])


        if outdir:
            cleanup_temp_dirs([outdir])


def test_original_vs_optimized_consistency():
    """Test that original and optimized algorithms produce same results."""
    seed_list = load_seeds()
    if len(seed_list) < 2:
        pytest.skip("No seed sequences found in database")
    
    # Use a moderate length sequence
    test_seq = generate_test_sequence(500, seed=999)
    # Insert a known pdif site from first seed
    first_rec = next(SeqIO.parse('PdifFinder/data/redundant.seed.fa', 'fasta'))
    pdif_site = str(first_rec.seq)
    position = 100
    test_seq = test_seq[:position] + pdif_site + test_seq[position + len(pdif_site):]
    
    outdir_orig = tempfile.mkdtemp(prefix='pdif_test_orig_')
    outdir_opt = tempfile.mkdtemp(prefix='pdif_test_opt_')
    try:
        # Setup directories
        for outdir in [outdir_orig, outdir_opt]:
            tmp_dir = Path(outdir) / 'tmp' / 'seedSearch' / '0'
            tmp_dir.mkdir(parents=True, exist_ok=True)
        
        # Run original algorithm
        pf.findMatchFragmentThreadOriginal(
            num4=0,
            pos0=0,
            outdir=outdir_orig,
            name='1',
            seedList=seed_list,
            seq=test_seq,
            start1=1,
            end1=len(seed_list) // 2
        )
        # Run optimized algorithm
        pf.findMatchFragmentThread(
            num4=0,
            pos0=0,
            outdir=outdir_opt,
            name='1',
            seedList=seed_list,
            seq=test_seq,
            start1=1,
            end1=len(seed_list) // 2
        )
        
        # Compare outputs
        outfile_orig = Path(outdir_orig) / 'tmp' / 'seedSearch' / '0' / '1'
        outfile_opt = Path(outdir_opt) / 'tmp' / 'seedSearch' / '0' / '1'
        positions_orig = []
        positions_opt = []
        if outfile_orig.exists():
            with open(outfile_orig, 'r') as f:
                positions_orig = [int(line.strip().split('|')[0]) for line in f if line.strip()]
        if outfile_opt.exists():
            with open(outfile_opt, 'r') as f:
                positions_opt = [int(line.strip().split('|')[0]) for line in f if line.strip()]
        
        # They should match (allow for ordering differences)
        assert set(positions_orig) == set(positions_opt), \
            f"Original and optimized results differ: original {positions_orig}, optimized {positions_opt}"
        
    finally:
        cleanup_temp_dirs([outdir_orig, outdir_opt])


def test_edge_case_empty_sequence():
    """Test scanning with empty sequence (should handle gracefully)."""
    seed_list = load_seeds()
    if len(seed_list) < 2:
        pytest.skip("No seed sequences found in database")
    
    test_seq = ""
    outdir = tempfile.mkdtemp(prefix='pdif_test_empty_')
    try:
        tmp_dir = Path(outdir) / 'tmp' / 'seedSearch' / '0'
        tmp_dir.mkdir(parents=True, exist_ok=True)
        # Should not crash
        pf.findMatchFragmentThread(
            num4=0,
            pos0=0,
            outdir=outdir,
            name='1',
            seedList=seed_list,
            seq=test_seq,
            start1=1,
            end1=len(seed_list) // 2
        )
        outfile = tmp_dir / '1'
        # Expect empty output (or no file)
        if outfile.exists():
            with open(outfile, 'r') as f:
                content = f.read()
                assert content == "", f"Expected empty output, got {content}"
    finally:
        cleanup_temp_dirs([outdir])


def test_edge_case_short_sequence():
    """Test scanning with sequence shorter than 28bp."""
    seed_list = load_seeds()
    if len(seed_list) < 2:
        pytest.skip("No seed sequences found in database")
    
    test_seq = generate_test_sequence(20, seed=111)  # 20bp < 28
    outdir = tempfile.mkdtemp(prefix='pdif_test_short_')
    try:
        tmp_dir = Path(outdir) / 'tmp' / 'seedSearch' / '0'
        tmp_dir.mkdir(parents=True, exist_ok=True)
        pf.findMatchFragmentThread(
            num4=0,
            pos0=0,
            outdir=outdir,
            name='1',
            seedList=seed_list,
            seq=test_seq,
            start1=1,
            end1=len(seed_list) // 2
        )
        outfile = tmp_dir / '1'
        # Should produce no matches (since sequence too short for 28bp window)
        if outfile.exists():
            with open(outfile, 'r') as f:
                positions = [int(line.strip().split('|')[0]) for line in f if line.strip()]
                assert len(positions) == 0, f"Expected no matches for short sequence, got {positions}"
    finally:
        cleanup_temp_dirs([outdir])


def test_edge_case_no_pdif_sites():
    """Test scanning with sequence containing no pdif sites."""
    seed_list = load_seeds()
    if len(seed_list) < 2:
        pytest.skip("No seed sequences found in database")
    
    # Random sequence unlikely to contain pdif sites
    test_seq = generate_test_sequence(1000, seed=222)
    outdir = tempfile.mkdtemp(prefix='pdif_test_nopdif_')
    try:
        tmp_dir = Path(outdir) / 'tmp' / 'seedSearch' / '0'
        tmp_dir.mkdir(parents=True, exist_ok=True)
        pf.findMatchFragmentThread(
            num4=0,
            pos0=0,
            outdir=outdir,
            name='1',
            seedList=seed_list,
            seq=test_seq,
            start1=1,
            end1=len(seed_list) // 2
        )
        outfile = tmp_dir / '1'
        # Might produce some false positives due to mismatch thresholds, but we can accept some
        # For now just ensure no crash
        if outfile.exists():
            with open(outfile, 'r') as f:
                positions = [int(line.strip().split('|')[0]) for line in f if line.strip()]
                # Log number of false positives
                if len(positions) > 0:
                    print(f"Note: {len(positions)} potential false positives detected in random sequence")
    finally:
        cleanup_temp_dirs([outdir])


if __name__ == '__main__':
    test_find_feature_end_seq_returns_whole_plasmid()
    print("✓ findFeatureEndSeq test passed")
    test_whole_plasmid_scanning_detects_pdif_sites()
    print("✓ whole plasmid scanning test passed")
    print("All whole plasmid tests passed!")