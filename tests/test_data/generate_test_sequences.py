#!/usr/bin/env python3
"""
Generate synthetic test sequences for pdifFinder testing.
Creates FASTA files with known pdif sites at specific positions.
"""

import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def generate_random_dna_sequence(length):
    """Generate a random DNA sequence of given length."""
    return ''.join(random.choice('ACGT') for _ in range(length))

def insert_pdif_site(sequence, position, pdif_site):
    """
    Insert a pdif site at the specified position in the sequence.
    pdif_site should be a 28bp string.
    """
    if position + 28 > len(sequence):
        raise ValueError(f"Position {position} + 28 exceeds sequence length {len(sequence)}")
    
    # Convert to list for modification
    seq_list = list(sequence)
    for i in range(28):
        seq_list[position + i] = pdif_site[i]
    
    return ''.join(seq_list)

def create_positive_test_sequence():
    """Create a sequence with 3 pdif sites at positions 100-127, 500-527, 1000-1027."""
    # Create a 1500bp random sequence
    seq_length = 1500
    sequence = generate_random_dna_sequence(seq_length)
    
    # Use actual pdif sites from the database
    pdif_sites = [
        "ATTTAACATAAGGGCTGTTATACGAAAT",  # pdif1 from database
        "AGTACATATAACAAAGATTATGTTAAAT",  # pdif2 from database  
        "AATTAAAATACCTTCTGTTATGTGCAAC",  # pdif3 from database
    ]
    
    # Insert pdif sites at specified positions
    positions = [100, 500, 1000]
    
    for pos, pdif in zip(positions, pdif_sites):
        sequence = insert_pdif_site(sequence, pos, pdif)
    
    return sequence

def create_negative_test_sequence():
    """Create a sequence with no pdif sites."""
    # Create a 2000bp random sequence
    seq_length = 2000
    sequence = generate_random_dna_sequence(seq_length)
    
    # Ensure no pdif sites by avoiding patterns that might match
    # We'll just use random sequence - probability of matching a specific 28bp pattern is extremely low
    return sequence

def create_edge_case_test_sequence():
    """Create a sequence with partial/near-miss pdif sites."""
    seq_length = 1800
    sequence = generate_random_dna_sequence(seq_length)
    
    # Create partial pdif sites (first 14bp only)
    partial_pdif = "ATTTAACATAAGGG"  # First half of pdif1
    
    # Insert partial sites at various positions
    partial_positions = [200, 600, 900]
    for pos in partial_positions:
        for i in range(14):
            if pos + i < len(sequence):
                sequence = sequence[:pos + i] + partial_pdif[i] + sequence[pos + i + 1:]
    
    # Create near-miss sites (1-2 mutations)
    near_miss_pdif = "ATTTAACATAAGGGCTGTTATACGAAAG"  # Last base changed from T to G
    
    # Insert near-miss sites
    near_miss_positions = [400, 800]
    for pos in near_miss_positions:
        if pos + 28 <= len(sequence):
            for i in range(28):
                sequence = sequence[:pos + i] + near_miss_pdif[i] + sequence[pos + i + 1:]
    
    return sequence

def save_fasta_file(filename, sequence_id, sequence, description=""):
    """Save a sequence as a FASTA file."""
    record = SeqRecord(
        Seq(sequence),
        id=sequence_id,
        description=description
    )
    
    with open(filename, "w") as f:
        SeqIO.write(record, f, "fasta")
    
    print(f"Saved {filename}")

def main():
    """Generate all test sequences."""
    # Set random seed for reproducibility
    random.seed(42)
    
    # Create output directory
    output_dir = "/work/tests/test_data"
    os.makedirs(output_dir, exist_ok=True)
    
    print("Generating synthetic test sequences...")
    
    # 1. Positive test sequence
    print("\n1. Generating positive test sequence...")
    positive_seq = create_positive_test_sequence()
    save_fasta_file(
        os.path.join(output_dir, "positive_test.fasta"),
        "test_positive",
        positive_seq,
        "Synthetic sequence with 3 pdif sites at positions 100-127, 500-527, 1000-1027"
    )
    
    # 2. Negative test sequence  
    print("\n2. Generating negative test sequence...")
    negative_seq = create_negative_test_sequence()
    save_fasta_file(
        os.path.join(output_dir, "negative_test.fasta"),
        "test_negative",
        negative_seq,
        "Synthetic sequence with no pdif sites"
    )
    
    # 3. Edge case test sequence
    print("\n3. Generating edge case test sequence...")
    edge_case_seq = create_edge_case_test_sequence()
    save_fasta_file(
        os.path.join(output_dir, "edge_case_test.fasta"),
        "test_edge_case",
        edge_case_seq,
        "Synthetic sequence with partial and near-miss pdif sites"
    )
    
    print("\nAll test sequences generated successfully!")
    
    # Print verification info
    print("\nVerification information:")
    print(f"Positive test sequence length: {len(positive_seq)} bp")
    print(f"Negative test sequence length: {len(negative_seq)} bp")
    print(f"Edge case test sequence length: {len(edge_case_seq)} bp")
    
    # Show pdif sites in positive test
    print("\nPdif sites in positive test (positions 100-127, 500-527, 1000-1027):")
    positions = [100, 500, 1000]
    for pos in positions:
        site = positive_seq[pos:pos+28]
        print(f"  Position {pos}-{pos+27}: {site}")

if __name__ == "__main__":
    main()