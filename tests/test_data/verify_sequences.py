#!/usr/bin/env python3
"""
Verify that generated test sequences are readable by Bio.SeqIO.
"""

from Bio import SeqIO
import os

def verify_fasta_file(filename):
    """Verify a FASTA file can be read by Bio.SeqIO."""
    try:
        records = list(SeqIO.parse(filename, "fasta"))
        if len(records) == 0:
            return False, "No records found"
        
        record = records[0]
        return True, f"Success: ID='{record.id}', Length={len(record.seq)} bp"
    except Exception as e:
        return False, f"Error: {str(e)}"

def check_pdif_sites(sequence, expected_positions):
    """Check if pdif sites are at expected positions."""
    results = []
    for pos in expected_positions:
        if pos + 28 <= len(sequence):
            site = sequence[pos:pos+28]
            results.append((pos, site))
    return results

def main():
    """Verify all test sequences."""
    test_dir = "/work/tests/test_data"
    
    test_files = [
        ("positive_test.fasta", [100, 500, 1000]),
        ("negative_test.fasta", []),
        ("edge_case_test.fasta", []),
    ]
    
    print("Verifying test sequences with Bio.SeqIO...")
    print("=" * 60)
    
    all_passed = True
    
    for filename, expected_positions in test_files:
        filepath = os.path.join(test_dir, filename)
        print(f"\nVerifying: {filename}")
        print("-" * 40)
        
        # Check if file exists
        if not os.path.exists(filepath):
            print(f"  ❌ File not found: {filepath}")
            all_passed = False
            continue
        
        # Verify Bio.SeqIO can read it
        success, message = verify_fasta_file(filepath)
        if success:
            print(f"  ✅ {message}")
            
            # For positive test, check pdif sites
            if expected_positions:
                records = list(SeqIO.parse(filepath, "fasta"))
                sequence = str(records[0].seq)
                pdif_sites = check_pdif_sites(sequence, expected_positions)
                
                print(f"  Checking pdif sites at positions {expected_positions}:")
                for pos, site in pdif_sites:
                    print(f"    Position {pos}: {site}")
                    
                    # Verify it's a valid pdif site (28bp, only A,C,G,T)
                    if len(site) == 28 and all(base in 'ACGT' for base in site):
                        print(f"      ✅ Valid 28bp DNA sequence")
                    else:
                        print(f"      ❌ Invalid sequence")
                        all_passed = False
        else:
            print(f"  ❌ {message}")
            all_passed = False
    
    print("\n" + "=" * 60)
    if all_passed:
        print("✅ All test sequences verified successfully!")
    else:
        print("❌ Some test sequences failed verification.")
    
    return all_passed

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)