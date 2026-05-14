#!/usr/bin/env python3
import sys
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

# Use a known pdif site from database
pdif_site = "ATTTAACATAAGGGCTGTTATACGAAAT"  # pdif1
plasmid_length = 2000
spacing = 500
start = spacing
sequence = generate_test_sequence(plasmid_length, seed=42)
sequence = insert_pattern_at_position(sequence, pdif_site, start)

print(f"Created plasmid length {len(sequence)}")
print(f"Inserted pdif site at {start}-{start+28}")

temp_fasta = save_sequence_to_fasta(sequence, "test_plasmid", "Single pdif test")
temp_output_dir = tempfile.mkdtemp(prefix='pdif_test_')
print(f"Output dir: {temp_output_dir}")

try:
    return_code, stdout, stderr = run_pdif_finder(temp_fasta, temp_output_dir)
    print(f"Return code: {return_code}")
    if stderr:
        print(f"Stderr: {stderr[:200]}")
    
    results = parse_pdif_output(temp_output_dir)
    print(f"Found {len(results['pdif_sites'])} pdif sites")
    for site in results['pdif_sites']:
        print(f"  Position: {site['start']}-{site['end']}, Sequence: {site['pdif_site']}")
        if abs(site['start'] - start) <= 5:
            print("  -> MATCH expected position!")
finally:
    cleanup_temp_files([temp_fasta])
    cleanup_temp_dirs([temp_output_dir])