#!/usr/bin/env python3
import sys
sys.path.insert(0, '.')
from Bio import SeqIO
import os

# Load seedList
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
seedListPairLength = len(seedList) // 2
print(f'seedListPairLength: {seedListPairLength}')

# Load positive test sequence
pos_fasta = 'tests/test_data/positive_test.fasta'
seq_record = next(SeqIO.parse(pos_fasta, 'fasta'))
seq = str(seq_record.seq).upper()
print(f'Sequence length: {len(seq)}')

# Expected positions
expected_positions = [100, 500, 1000]

# Define mismatch counting as in findMatchFragmentThread
def check_window(start, seq, seed_c, seed_d):
    # start is 0-based
    if start + 28 > len(seq):
        return False
    maxMismatchXerC = 3
    maxMismatchXerD = 2
    initPos = start
    endPos = start + 11
    tempSeqC = seq[initPos:endPos]
    tempMismatchNumberC = 0
    tempMismatchNumberD = 0
    tempMismatchNumberCLeft = 0
    # first 5 positions
    for l in range(5):
        if tempSeqC[l] != seed_c[l]:
            tempMismatchNumberCLeft += 1
    if tempMismatchNumberCLeft > maxMismatchXerC:
        return False
    else:
        for j in range(11):
            if tempSeqC[j] != seed_c[j]:
                tempMismatchNumberC += 1
    if tempMismatchNumberC > maxMismatchXerC:
        return False
    else:
        if (endPos + 17) <= len(seq):
            tempSeqD = seq[initPos + 17:endPos + 17]
            for m in range(11):
                if tempSeqD[m] != seed_d[m]:
                    tempMismatchNumberD += 1
            if tempMismatchNumberD > maxMismatchXerD:
                return False
            else:
                return True
    return False

# For each expected position, check each seed pair
for pos in expected_positions:
    print(f'\nChecking position {pos}:')
    window = seq[pos:pos+28]
    xerc = window[0:11]
    xerd = window[17:28]
    print(f'  XerC: {xerc}')
    print(f'  XerD: {xerd}')
    found = False
    for k in range(0, len(seedList), 2):
        seed_c = seedList[k]
        seed_d = seedList[k+1]
        if check_window(pos, seq, seed_c, seed_d):
            print(f'    matches seed pair {k//2}: C={seed_c} D={seed_d}')
            found = True
            break
    if not found:
        print('    NO MATCH')
        # Let's compute mismatches manually
        for k in range(0, len(seedList), 2):
            seed_c = seedList[k]
            seed_d = seedList[k+1]
            mismatch_c = sum(1 for a,b in zip(xerc, seed_c) if a != b)
            mismatch_d = sum(1 for a,b in zip(xerd, seed_d) if a != b)
            if mismatch_c <= 3 and mismatch_d <= 2:
                print(f'      but seed pair {k//2} passes simple mismatch: C={mismatch_c} D={mismatch_d}')
                # debug early exit
                # compute tempMismatchNumberCLeft
                tempMismatchNumberCLeft = sum(1 for a,b in zip(xerc[:5], seed_c[:5]) if a != b)
                print(f'        mismatches first5: {tempMismatchNumberCLeft}')
                # compute tempMismatchNumberC
                tempMismatchNumberC = mismatch_c
                print(f'        total XerC mismatches: {tempMismatchNumberC}')
                # check condition
                if tempMismatchNumberCLeft > 3:
                    print('        FAIL: first5 mismatches > 3')
                if tempMismatchNumberC > 3:
                    print('        FAIL: total XerC mismatches > 3')
                # check XerD region
                if (pos + 28) <= len(seq):
                    tempSeqD = seq[pos + 17: pos + 28]
                    mismatch_d = sum(1 for a,b in zip(tempSeqD, seed_d) if a != b)
                    print(f'        XerD mismatches: {mismatch_d}')
                else:
                    print('        XerD out of bounds')
                # check if (endPos + 17) <= len(seq)
                endPos = pos + 11
                if (endPos + 17) <= len(seq):
                    print('        XerD region within bounds')
                else:
                    print('        XerD region out of bounds')

# Also scan whole sequence for any matches
print('\n--- Scanning whole sequence for matches ---')
matches = []
for start in range(len(seq) - 27):
    for k in range(0, len(seedList), 2):
        if check_window(start, seq, seedList[k], seedList[k+1]):
            matches.append(start)
            break
print(f'Found {len(matches)} matches at positions: {matches}')
# compare with expected
for pos in expected_positions:
    if pos in matches:
        print(f'Position {pos} detected')
    else:
        print(f'Position {pos} MISSED')