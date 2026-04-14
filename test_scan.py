#!/usr/bin/env python3
import sys
sys.path.insert(0, '.')
import os
import tempfile
from Bio import SeqIO

# Build seedList as in findPdif
pdifDB = 'PdifFinder/data/redundant.seed.fa'
seedList = []
for rec in SeqIO.parse(pdifDB, 'fasta'):
    seq1 = str(rec.seq)
    XerC = seq1[0:11]
    XerD = seq1[17:28]
    seedList.append(XerC)
    seedList.append(XerD)
    seedList.append(XerD)
    seedList.append(XerC)
print(f'seedList length: {len(seedList)}')

# Create a simple exact match to the first seed pair
first_seq = next(SeqIO.parse(pdifDB, 'fasta'))
full_seq = str(first_seq.seq)  # 28bp
print(f'First seed sequence: {full_seq}')
xerc = full_seq[0:11]
xerd = full_seq[17:28]
print(f'XerC: {xerc}, XerD: {xerd}')

# Create a longer sequence with this pdif site at position 100
seq_length = 200
sequence = 'N' * seq_length
# insert at position 100
pos = 100
seq_list = list(sequence)
seq_list[pos:pos+28] = list(full_seq)
sequence = ''.join(seq_list)
print(f'Test sequence length: {len(sequence)}')

# Create temporary output directory
outdir = tempfile.mkdtemp(prefix='pdif_test_')
os.makedirs(outdir + '/tmp/seedSearch/0', exist_ok=True)

# Import pdifFinder module
import PdifFinder.pdifFinder as pf

# Call findMatchFragmentThread with appropriate arguments
num4 = 0
pos0 = 0
name = 1
seq = sequence  # whole plasmid
start1 = 1
end1 = 1  # only first seed pair
seedListPairLength = int(len(seedList) / 2)
print(f'seedListPairLength: {seedListPairLength}')
# We'll call with start1=1, end1=1 (first pair)
pf.findMatchFragmentThread(num4, pos0, outdir, name, seedList, seq, start1, end1)

# Check output file
outfile = outdir + '/tmp/seedSearch/0/1'
if os.path.exists(outfile):
    with open(outfile, 'r') as f:
        lines = f.readlines()
        print(f'Found {len(lines)} positions:')
        for line in lines:
            print(line.strip())
else:
    print('No output file created')

# Cleanup
import shutil
shutil.rmtree(outdir)