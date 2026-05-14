#!/usr/bin/env python3
import sys
import tempfile
import os
from unittest.mock import patch
sys.path.insert(0, '.')
from Bio import SeqIO
from tests.helpers import generate_test_sequence, run_pdif_finder, parse_pdif_output, cleanup_temp_files, cleanup_temp_dirs

# Mock findResistanceGene
with patch('PdifFinder.pdifFinder.findResistanceGene') as mock_find:
    mock_find.return_value = ["1000-2000"]
    
    # Load a known pdif site from database
    pdif_db = 'PdifFinder/data/redundant.seed.fa'
    first_rec = next(SeqIO.parse(pdif_db, 'fasta'))
    pdif_site = str(first_rec.seq)
    print(f"Seed sequence: {pdif_site}")
    
    # Create sequence with pdif site at known position
    background = generate_test_sequence(1000, seed=123)
    position = 200
    test_seq = background[:position] + pdif_site + background[position + len(pdif_site):]
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">test_with_pdif\n{test_seq}\n")
        fasta_file = f.name
    
    output_dir = tempfile.mkdtemp(prefix='pdif_test_seed_')
    try:
        # Run pdifFinder
        return_code, stdout, stderr = run_pdif_finder(
            input_file=fasta_file,
            output_dir=output_dir
        )
        
        print(f"Return code: {return_code}")
        print(f"Stdout:\n{stdout}")
        print(f"Stderr:\n{stderr}")
        
        # Parse output files to ensure they exist and have correct format
        results = parse_pdif_output(output_dir)
        
        print(f"pdif_sites found: {len(results['pdif_sites'])}")
        for site in results['pdif_sites']:
            print(f"  Position: {site['start']}-{site['end']}, Sequence: {site['pdif_site']}")
        
        # Check intermediate files
        seed_search_dir = os.path.join(output_dir, 'tmp', 'seedSearch')
        if os.path.exists(seed_search_dir):
            print(f"Seed search dir exists")
            for root, dirs, files in os.walk(seed_search_dir):
                for file in files:
                    path = os.path.join(root, file)
                    with open(path, 'r') as ff:
                        content = ff.read().strip()
                        if content:
                            print(f"File {path}: {content}")
        else:
            print(f"Seed search dir missing")
        
    finally:
        cleanup_temp_files([fasta_file])
        cleanup_temp_dirs([output_dir])