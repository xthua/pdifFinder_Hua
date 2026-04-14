#!/usr/bin/env python3
"""
Common test utilities for pdifFinder validation.
Provides helper functions for sequence generation, position checking, and validation.
"""

import os
import sys
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple, Dict, Optional
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def generate_test_sequence(length: int, seed: int = 42) -> str:
    """
    Generate a random DNA sequence of given length.
    
    Args:
        length: Length of sequence to generate
        seed: Random seed for reproducibility
        
    Returns:
        Random DNA sequence string
    """
    random.seed(seed)
    return ''.join(random.choice('ACGT') for _ in range(length))


def insert_pattern_at_position(sequence: str, pattern: str, position: int) -> str:
    """
    Insert a pattern at the specified position in the sequence.
    
    Args:
        sequence: Original DNA sequence
        pattern: Pattern to insert
        position: 0-based start position
        
    Returns:
        Modified sequence with pattern inserted
    """
    if position + len(pattern) > len(sequence):
        raise ValueError(f"Position {position} + pattern length {len(pattern)} exceeds sequence length {len(sequence)}")
    
    seq_list = list(sequence)
    for i in range(len(pattern)):
        seq_list[position + i] = pattern[i]
    
    return ''.join(seq_list)


def create_sequence_with_pdif_sites(
    length: int = 2000,
    pdif_sites: Optional[List[Tuple[str, int]]] = None,
    seed: int = 42
) -> str:
    """
    Create a test sequence with pdif sites at specified positions.
    
    Args:
        length: Total sequence length
        pdif_sites: List of (pdif_site_sequence, position) tuples
        seed: Random seed for background sequence
        
    Returns:
        DNA sequence with pdif sites inserted
    """
    # Default pdif sites if none provided
    if pdif_sites is None:
        pdif_sites = [
            ("ATTTAACATAAGGGCTGTTATACGAAAT", 100),  # pdif1 at position 100
            ("AGTACATATAACAAAGATTATGTTAAAT", 500),  # pdif2 at position 500
            ("AATTAAAATACCTTCTGTTATGTGCAAC", 1000), # pdif3 at position 1000
        ]
    
    # Generate random background sequence
    sequence = generate_test_sequence(length, seed)
    
    # Insert pdif sites
    for pdif_site, position in pdif_sites:
        sequence = insert_pattern_at_position(sequence, pdif_site, position)
    
    return sequence


def save_sequence_to_fasta(
    sequence: str,
    sequence_id: str,
    description: str = "",
    output_path: Optional[str] = None
) -> str:
    """
    Save a sequence to a FASTA file.
    
    Args:
        sequence: DNA sequence string
        sequence_id: Sequence identifier
        description: Sequence description
        output_path: Output file path (creates temp file if None)
        
    Returns:
        Path to the created FASTA file
    """
    if output_path is None:
        # Create temporary file
        fd, temp_path = tempfile.mkstemp(suffix='.fasta', prefix='test_')
        os.close(fd)
        output_path = temp_path
    
    record = SeqRecord(
        Seq(sequence),
        id=sequence_id,
        description=description
    )
    
    with open(output_path, "w") as f:
        SeqIO.write(record, f, "fasta")
    
    return output_path


def run_pdif_finder(
    input_file: str,
    output_dir: Optional[str] = None,
    circle_format: str = "true"
) -> Tuple[int, str, str]:
    """
    Run pdifFinder on an input file.
    
    Args:
        input_file: Path to input FASTA/GenBank file
        output_dir: Output directory (creates temp dir if None)
        circle_format: Circle format parameter
        
    Returns:
        Tuple of (return_code, stdout, stderr)
    """
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='pdif_test_')
    
    # Build command
    cmd = [
        sys.executable, "-m", "PdifFinder.pdifFinder",
        "--inFile", input_file,
        "--outdir", output_dir,
        "--circle", circle_format
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "Command timed out after 5 minutes"


def parse_pdif_output(output_dir: str) -> Dict[str, List[Dict]]:
    """
    Parse pdifFinder output files.
    
    Args:
        output_dir: Directory containing pdifFinder output
        
    Returns:
        Dictionary with parsed results from each output file
    """
    results = {
        'pdif_sites': [],
        'amr_genes': [],
        'pdif_modules': []
    }
    
    # Parse pdif_site.txt
    pdif_site_file = os.path.join(output_dir, 'pdif_site.txt')
    if os.path.exists(pdif_site_file):
        with open(pdif_site_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        results['pdif_sites'].append({
                            'sequence_id': parts[0],
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'pdif_site': parts[3]
                        })
    
    # Parse AMRgene.txt
    amr_gene_file = os.path.join(output_dir, 'AMRgene.txt')
    if os.path.exists(amr_gene_file):
        with open(amr_gene_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        results['amr_genes'].append({
                            'sequence_id': parts[0],
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'gene_name': parts[3],
                            'accession': parts[4],
                            'description': parts[5]
                        })
    
    # Parse pdifmodule_list.txt
    pdif_module_file = os.path.join(output_dir, 'pdifmodule_list.txt')
    if os.path.exists(pdif_module_file):
        with open(pdif_module_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        results['pdif_modules'].append({
                            'sequence_id': parts[0],
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'module_type': parts[3],
                            'genes': parts[4]
                        })
    
    return results


def validate_pdif_positions(
    expected_positions: List[Tuple[int, int]],
    actual_results: List[Dict],
    tolerance: int = 5
) -> Tuple[bool, List[str]]:
    """
    Validate that pdif sites were found at expected positions.
    
    Args:
        expected_positions: List of (start, end) tuples for expected pdif sites
        actual_results: List of actual pdif site results from parse_pdif_output
        tolerance: Allowable position difference in bp
        
    Returns:
        Tuple of (success, list of error messages)
    """
    errors = []
    
    # Check each expected position
    for exp_start, exp_end in expected_positions:
        found = False
        for result in actual_results:
            act_start = result['start']
            act_end = result['end']
            
            # Check if position is within tolerance
            if (abs(act_start - exp_start) <= tolerance and 
                abs(act_end - exp_end) <= tolerance):
                found = True
                break
        
        if not found:
            errors.append(f"Expected pdif site at position {exp_start}-{exp_end} not found")
    
    # Check for unexpected findings
    if len(actual_results) > len(expected_positions):
        extra = len(actual_results) - len(expected_positions)
        errors.append(f"Found {extra} unexpected pdif site(s)")
    
    return len(errors) == 0, errors


def print_validation_summary(
    test_name: str,
    success: bool,
    errors: List[str],
    details: Optional[Dict] = None
) -> None:
    """
    Print a formatted validation summary.
    
    Args:
        test_name: Name of the test
        success: Whether validation passed
        errors: List of error messages
        details: Additional details to print
    """
    print("\n" + "="*60)
    print(f"VALIDATION: {test_name}")
    print("="*60)
    
    if success:
        print("✅ PASS: All checks passed")
    else:
        print("❌ FAIL: Validation errors found")
        for error in errors:
            print(f"  - {error}")
    
    if details:
        print("\nDetails:")
        for key, value in details.items():
            print(f"  {key}: {value}")
    
    print("="*60)


def cleanup_temp_files(file_paths: List[str]) -> None:
    """
    Clean up temporary files.
    
    Args:
        file_paths: List of file paths to remove
    """
    for path in file_paths:
        try:
            if os.path.exists(path):
                os.remove(path)
        except Exception as e:
            print(f"Warning: Could not remove {path}: {e}")


def cleanup_temp_dirs(dir_paths: List[str]) -> None:
    """
    Clean up temporary directories.
    
    Args:
        dir_paths: List of directory paths to remove
    """
    for path in dir_paths:
        try:
            if os.path.exists(path):
                import shutil
                shutil.rmtree(path)
        except Exception as e:
            print(f"Warning: Could not remove directory {path}: {e}")


if __name__ == "__main__":
    # Test the helper functions
    print("Testing helper functions...")
    
    # Test sequence generation
    seq = generate_test_sequence(100)
    print(f"Generated sequence length: {len(seq)}")
    
    # Test pattern insertion
    test_seq = "A" * 100
    pattern = "ATCG"
    modified = insert_pattern_at_position(test_seq, pattern, 10)
    print(f"Pattern inserted at position 10: {modified[10:14]}")
    
    # Test sequence with pdif sites
    pdif_seq = create_sequence_with_pdif_sites(length=1500)
    print(f"Sequence with pdif sites length: {len(pdif_seq)}")
    
    print("\nHelper functions test completed successfully!")