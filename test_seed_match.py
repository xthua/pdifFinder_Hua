#!/usr/bin/env python3
import sys
sys.path.insert(0, '.')
from Bio import SeqIO

# Load seed database
seed_file = 'PdifFinder/data/redundant.seed.fa'
seed_pairs = []
for rec in SeqIO.parse(seed_file, 'fasta'):
    seq = str(rec.seq)
    if len(seq) != 28:
        continue
    xerc = seq[0:11]
    xerd = seq[17:28]
    seed_pairs.append((xerc, xerd))
    # also add reverse orientation? they add both orientations
    seed_pairs.append((xerd, xerc))
print(f'Loaded {len(seed_pairs)} seed pairs')

# Synthetic pdif sites from pdifdatabase.fasta
synthetic = [
    ("ATTTAACATAAGGGCTGTTATACGAAAT", 100),
    ("AGTACATATAACAAAGATTATGTTAAAT", 500),
    ("AATTAAAATACCTTCTGTTATGTGCAAC", 1000),
]

def mismatch(seq1, seq2):
    return sum(1 for a,b in zip(seq1, seq2) if a != b)

for pdif, pos in synthetic:
    xerc = pdif[0:11]
    xerd = pdif[17:28]
    print(f'\nPDIF site at {pos}:')
    print(f'  XerC: {xerc}')
    print(f'  XerD: {xerd}')
    found = False
    for i, (sc, sd) in enumerate(seed_pairs):
        if mismatch(xerc, sc) <= 3 and mismatch(xerd, sd) <= 2:
            print(f'    matches seed pair {i}: XerC={sc} XerD={sd}')
            print(f'    mismatches: C={mismatch(xerc, sc)} D={mismatch(xerd, sd)}')
            found = True
            break
    if not found:
        print('    NO MATCHING SEED PAIR')

# Also check if any seed pair matches with reversed orientation (XerD/XerC)
print('\nChecking reversed orientation (XerD/XerC):')
for pdif, pos in synthetic:
    xerc = pdif[0:11]
    xerd = pdif[17:28]
    found = False
    for i, (sc, sd) in enumerate(seed_pairs):
        if mismatch(xerd, sc) <= 3 and mismatch(xerc, sd) <= 2:
            print(f'    PDIF {pos} matches reversed seed pair {i}')
            found = True
            break