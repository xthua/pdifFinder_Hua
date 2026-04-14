#!/usr/bin/env python3
import sys
sys.path.insert(0, '.')
import os
import tempfile
import shutil
from Bio import SeqIO

# Import pdifFinder
import PdifFinder.pdifFinder as pf

# Create output directory
outdir = tempfile.mkdtemp(prefix='pdif_test_')
print(f'Output dir: {outdir}')
pf.makeOutdir(outdir)

# Positive test file
inFile = 'tests/test_data/positive_test.fasta'
pdifDB = 'PdifFinder/data/redundant.seed.fa'
blastnPath = ''  # not used
resistanceGenePosList = []

# Call findPdif
print('Calling findPdif...')
flag = pf.findPdif(inFile, outdir, blastnPath, pdifDB, resistanceGenePosList)
print(f'Return flag: {flag}')

# Check for output files
possible_pdif_site = outdir + '/tmp/possible_pdif_site.txt'
if os.path.exists(possible_pdif_site):
    with open(possible_pdif_site, 'r') as f:
        lines = f.readlines()
        print(f'Possible pdif sites ({len(lines)} lines):')
        for line in lines:
            print(line.strip())
else:
    print('No possible_pdif_site.txt')

# Check seedSearch directory
seed_search_dir = outdir + '/tmp/seedSearch/0'
if os.path.exists(seed_search_dir):
    for fname in os.listdir(seed_search_dir):
        with open(os.path.join(seed_search_dir, fname), 'r') as f:
            positions = [line.strip() for line in f if line.strip()]
            print(f'File {fname}: {len(positions)} positions')
            if positions:
                print('  ' + ', '.join(positions[:10]))

# Cleanup
shutil.rmtree(outdir)