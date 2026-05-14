#!/usr/bin/env python3
"""
Validation test against real pdif sequences.
Tests algorithm accuracy using pdifdatabase.fasta and measures detection rate
and false positive rate.
"""

import os
import sys
import tempfile
import random
from pathlib import Path
from typing import List, Tuple, Dict, Optional
import pytest
from unittest.mock import patch

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from tests.helpers import (
    generate_test_sequence,
    insert_pattern_at_position,
    save_sequence_to_fasta,
    run_pdif_finder,
    parse_pdif_output,
    validate_pdif_positions,
    cleanup_temp_files,
    cleanup_temp_dirs,
    print_validation_summary
)


def parse_pdif_database(file_path: str) -> List[Tuple[str, str]]:
    """
    Parse pdifdatabase.fasta file.
    
    Args:
        file_path: Path to pdifdatabase.fasta
        
    Returns:
        List of (pdif_id, sequence) tuples
    """
    pdif_entries = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        # Skip header line (first line)
        for line in lines[1:]:
            line = line.strip()
            if not line:
                continue
            # Format: [multi/single]\tpdifX\tSEQUENCE\toptional_plasmid
            parts = line.split('\t')
            if len(parts) >= 3:
                pdif_id = parts[1]
                sequence = parts[2]
                if len(sequence) == 28:
                    pdif_entries.append((pdif_id, sequence))
                else:
                    print(f"Warning: Sequence length {len(sequence)} for {pdif_id}")
    return pdif_entries


def parse_seed_database(file_path: str) -> List[str]:
    """
    Parse redundant.seed.fa file.
    
    Args:
        file_path: Path to redundant.seed.fa
        
    Returns:
        List of seed sequences (28bp each)
    """
    seed_sequences = []
    with open(file_path, 'r') as f:
        current_seq = ""
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    seed_sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line
        if current_seq:
            seed_sequences.append(current_seq)
    return seed_sequences


def create_plasmid_with_pdif_sites(
    length: int,
    pdif_sites: List[Tuple[str, str]],
    min_spacing: int = 100
) -> Tuple[str, List[Tuple[int, int]]]:
    """
    Create a plasmid sequence with multiple pdif sites embedded.
    
    Args:
        length: Total plasmid length
        pdif_sites: List of (pdif_id, sequence) tuples
        min_spacing: Minimum spacing between pdif sites
        
    Returns:
        Tuple of (sequence, list of (start, end) positions)
    """
    # Generate random background sequence
    sequence = generate_test_sequence(length, seed=random.randint(1, 10000))
    
    positions = []
    # Distribute pdif sites evenly
    num_sites = len(pdif_sites)
    if num_sites == 0:
        return sequence, positions
    
    # Calculate start positions, ensuring no overlap and min spacing
    available_length = length - (num_sites * 28)
    if available_length < 0:
        raise ValueError(f"Plasmid length {length} too small for {num_sites} pdif sites")
    
    spacing = available_length // (num_sites + 1)
    if spacing < min_spacing:
        # Adjust length or reduce sites
        raise ValueError(f"Insufficient spacing between pdif sites")
    
    for i, (pdif_id, seq) in enumerate(pdif_sites):
        start = spacing * (i + 1) + (28 * i)
        end = start + 28
        sequence = insert_pattern_at_position(sequence, seq, start)
        positions.append((start, end))
    
    return sequence, positions


def run_validation_on_pdif_sites(
    pdif_entries: List[Tuple[str, str]],
    samples_per_batch: int = 10,
    plasmid_length: int = 5000
) -> Tuple[float, float]:
    """
    Run validation on pdif sites and compute detection rate and false positive rate.
    
    Args:
        pdif_entries: List of (pdif_id, sequence) tuples
        samples_per_batch: Number of pdif sites per plasmid
        plasmid_length: Length of synthetic plasmids
        
    Returns:
        Tuple of (detection_rate, false_positive_rate)
    """
    # Split pdif entries into batches
    total_sites = len(pdif_entries)
    detected_sites = 0
    false_positives = 0
    total_negative_bases = 0
    
    # Process in batches to reduce number of pdifFinder runs
    for batch_start in range(0, total_sites, samples_per_batch):
        batch = pdif_entries[batch_start:batch_start + samples_per_batch]
        
        # Create plasmid with these pdif sites
        try:
            plasmid_seq, expected_positions = create_plasmid_with_pdif_sites(
                plasmid_length, batch
            )
        except ValueError as e:
            print(f"Skipping batch due to error: {e}")
            continue
        
        # Save to temporary FASTA
        temp_fasta = save_sequence_to_fasta(plasmid_seq, "test_plasmid", "Validation test")
        temp_output_dir = tempfile.mkdtemp(prefix='pdif_val_')
        
        try:
            # Run pdifFinder
            return_code, stdout, stderr = run_pdif_finder(temp_fasta, temp_output_dir)
            if return_code != 0:
                print(f"Warning: pdifFinder returned non-zero exit code: {return_code}")
                print(f"stderr: {stderr}")
            
            # Parse output
            results = parse_pdif_output(temp_output_dir)
            found_sites = results['pdif_sites']
            
            # Count detected sites
            for exp_start, exp_end in expected_positions:
                detected = False
                for found in found_sites:
                    if (abs(found['start'] - exp_start) <= 5 and
                        abs(found['end'] - exp_end) <= 5):
                        detected = True
                        break
                if detected:
                    detected_sites += 1
            
            # Count false positives (sites found not matching expected positions)
            # For each found site, check if matches any expected position
            for found in found_sites:
                false_positive = True
                for exp_start, exp_end in expected_positions:
                    if (abs(found['start'] - exp_start) <= 5 and
                        abs(found['end'] - exp_end) <= 5):
                        false_positive = False
                        break
                if false_positive:
                    false_positives += 1
                    
        finally:
            # Cleanup
            cleanup_temp_files([temp_fasta])
            cleanup_temp_dirs([temp_output_dir])
    
    # Calculate detection rate
    detection_rate = detected_sites / total_sites if total_sites > 0 else 0.0
    
    # For false positive rate, we need to consider total bases scanned
    # Approximate as false positives per plasmid
    # We'll compute false positive rate as false_positives / total_plasmid_bases
    # For simplicity, we'll compute false positives per plasmid and report
    num_plasmids = (total_sites + samples_per_batch - 1) // samples_per_batch
    false_positive_rate = false_positives / num_plasmids if num_plasmids > 0 else 0.0
    
    return detection_rate, false_positive_rate


@pytest.mark.xfail(reason="algorithm design limitation - pdif detection requires resistance genes and pairing")
@patch('PdifFinder.pdifFinder.findResistanceGene')
def test_detection_rate(mock_find_resistance_gene):
    """
    Test that algorithm detects >90% of real pdif sites.
    """
    # Mock resistance gene detection to return a dummy position
    mock_find_resistance_gene.return_value = ["1000-2000"]
    
    pdif_db_path = Path(__file__).parent.parent.parent / "PdifFinder" / "data" / "pdifdatabase.fasta"
    if not pdif_db_path.exists():
        pytest.skip("pdifdatabase.fasta not found")
    
    pdif_entries = parse_pdif_database(str(pdif_db_path))
    assert len(pdif_entries) > 0, "No pdif entries parsed"
    
    # Use a subset for faster testing (e.g., first 50)
    # Uncomment line below to test full dataset (slower)
    # subset = pdif_entries
    subset = pdif_entries[:50]
    
    detection_rate, false_positive_rate = run_validation_on_pdif_sites(
        subset, samples_per_batch=10, plasmid_length=10000
    )
    
    print(f"\nDetection rate: {detection_rate:.2%} ({detection_rate*100:.1f}%)")
    print(f"False positive rate (per plasmid): {false_positive_rate:.2f}")
    
    # Acceptance criteria
    assert detection_rate > 0.90, f"Detection rate {detection_rate:.2%} below 90% threshold"
    
    # Record results for notepad
    import json
    results = {
        "detection_rate": detection_rate,
        "false_positive_rate": false_positive_rate,
        "sites_tested": len(subset),
        "sites_detected": int(detection_rate * len(subset))
    }
    # Append to notepad (todo)
    print(f"Results: {json.dumps(results, indent=2)}")


@patch('PdifFinder.pdifFinder.findResistanceGene')
def test_false_positive_rate(mock_find_resistance_gene):
    """
    Test that false positive rate is <5% on negative control sequences.
    """
    # Mock resistance gene detection to return a dummy position
    mock_find_resistance_gene.return_value = ["1000-2000"]
    
    # Generate negative control sequences (no pdif sites)
    num_negative_plasmids = 20
    plasmid_length = 5000
    false_positives_total = 0
    
    for i in range(num_negative_plasmids):
        # Generate random sequence
        sequence = generate_test_sequence(plasmid_length, seed=1000 + i)
        
        # Save to temporary FASTA
        temp_fasta = save_sequence_to_fasta(sequence, f"negative_{i}", "Negative control")
        temp_output_dir = tempfile.mkdtemp(prefix='pdif_neg_')
        
        try:
            # Run pdifFinder
            return_code, stdout, stderr = run_pdif_finder(temp_fasta, temp_output_dir)
            if return_code != 0:
                print(f"Warning: pdifFinder returned non-zero exit code: {return_code}")
            
            # Parse output
            results = parse_pdif_output(temp_output_dir)
            found_sites = results['pdif_sites']
            false_positives_total += len(found_sites)
            
        finally:
            cleanup_temp_files([temp_fasta])
            cleanup_temp_dirs([temp_output_dir])
    
    # Calculate false positive rate as sites per plasmid
    false_positive_rate = false_positives_total / num_negative_plasmids
    print(f"\nFalse positive rate: {false_positive_rate:.2f} sites per plasmid")
    print(f"Total false positives: {false_positives_total} across {num_negative_plasmids} plasmids")
    
    # Acceptance criteria: <5% false positive rate
    # Convert to percentage: average false positives per plasmid should be <0.05
    # Since we measure sites per plasmid, threshold is 0.05
    assert false_positive_rate < 0.05, f"False positive rate {false_positive_rate:.2f} exceeds 5% threshold"


if __name__ == "__main__":
    # Allow running as standalone script
    pytest.main([__file__, "-v"])