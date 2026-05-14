#!/usr/bin/env python3
import sys
import os
import tempfile
sys.path.insert(0, '.')
from tests.helpers import (
    generate_test_sequence,
    insert_pattern_at_position,
    save_sequence_to_fasta,
    run_pdif_finder,
    parse_pdif_output,
    cleanup_temp_files,
    cleanup_temp_dirs
)

# Use a seed sequence that should definitely match
seed_seq = "ACTGCGCATAAGAGATTTTATGTTAAAT"  # from redundant.seed.fa
plasmid_length = 5000
spacing = 1000
start = spacing
sequence = generate_test_sequence(plasmid_length, seed=123)
sequence = insert_pattern_at_position(sequence, seed_seq, start)

print(f"Created plasmid length {len(sequence)}")
print(f"Inserted seed at {start}-{start+28}")

temp_fasta = save_sequence_to_fasta(sequence, "test_plasmid", "Seed test")
temp_output_dir = tempfile.mkdtemp(prefix='pdif_debug_')
print(f"Output dir: {temp_output_dir}")

try:
    return_code, stdout, stderr = run_pdif_finder(temp_fasta, temp_output_dir)
    print(f"Return code: {return_code}")
    if stdout:
        print(f"Stdout (first 500 chars): {stdout[:500]}")
    if stderr:
        print(f"Stderr (first 500 chars): {stderr[:500]}")
    
    # Check intermediate files
    seed_search_dir = os.path.join(temp_output_dir, 'tmp', 'seedSearch')
    if os.path.exists(seed_search_dir):
        print(f"Seed search dir exists")
        for root, dirs, files in os.walk(seed_search_dir):
            for file in files:
                path = os.path.join(root, file)
                with open(path, 'r') as f:
                    content = f.read().strip()
                    if content:
                        print(f"File {path}: {content}")
    else:
        print(f"Seed search dir missing")
    
    # Check possible pdif site file
    possible_site_file = os.path.join(temp_output_dir, 'tmp', 'possible_pdif_site.txt')
    if os.path.exists(possible_site_file):
        with open(possible_site_file, 'r') as f:
            lines = f.readlines()
            print(f"Possible pdif sites ({len(lines)} lines):")
            for line in lines:
                print(f"  {line.strip()}")
    else:
        print(f"Possible pdif site file missing")
    
    # Check final output
    results = parse_pdif_output(temp_output_dir)
    print(f"Found {len(results['pdif_sites'])} pdif sites")
    for site in results['pdif_sites']:
        print(f"  Position: {site['start']}-{site['end']}, Sequence: {site['pdif_site']}")
    
finally:
    cleanup_temp_files([temp_fasta])
    cleanup_temp_dirs([temp_output_dir])