#!/usr/bin/env python3
"""
Validate that pdifFinder algorithm scans entire sequence correctly.
This script tests that the scanning algorithm doesn't miss regions and
properly identifies pdif sites at various positions.
"""

import argparse
import sys
import os
import tempfile
from pathlib import Path

# Add parent directory to path to import helpers
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from helpers import (
    create_sequence_with_pdif_sites,
    save_sequence_to_fasta,
    run_pdif_finder,
    parse_pdif_output,
    validate_pdif_positions,
    print_validation_summary,
    cleanup_temp_files,
    cleanup_temp_dirs
)


def test_scanning_algorithm_basic() -> bool:
    """
    Test basic scanning algorithm with pdif sites at standard positions.
    """
    print("Testing basic scanning algorithm...")
    
    # Create test sequence with pdif sites at positions 100, 500, 1000
    sequence = create_sequence_with_pdif_sites(
        length=1500,
        pdif_sites=[
            ("ATTTAACATAAGGGCTGTTATACGAAAT", 100),   # pdif1
            ("AGTACATATAACAAAGATTATGTTAAAT", 500),   # pdif2
            ("AATTAAAATACCTTCTGTTATGTGCAAC", 1000),  # pdif3
        ]
    )
    
    # Save to temporary FASTA file
    fasta_file = save_sequence_to_fasta(
        sequence,
        "test_scan_basic",
        "Test sequence for basic scanning validation"
    )
    
    # Create temporary output directory
    output_dir = tempfile.mkdtemp(prefix='pdif_scan_test_')
    
    try:
        # Run pdifFinder
        return_code, stdout, stderr = run_pdif_finder(fasta_file, output_dir)
        
        if return_code != 0:
            print(f"pdifFinder failed with return code {return_code}")
            print(f"stderr: {stderr}")
            return False
        
        # Parse results
        results = parse_pdif_output(output_dir)
        pdif_sites = results['pdif_sites']
        
        # Expected positions (1-based for pdifFinder output)
        expected_positions = [
            (101, 128),   # position 100-127 (0-based) -> 101-128 (1-based)
            (501, 528),   # position 500-527 -> 501-528
            (1001, 1028)  # position 1000-1027 -> 1001-1028
        ]
        
        # Validate
        success, errors = validate_pdif_positions(expected_positions, pdif_sites)
        
        # Print summary
        details = {
            "Sequence length": f"{len(sequence)} bp",
            "Expected pdif sites": len(expected_positions),
            "Found pdif sites": len(pdif_sites),
            "pdifFinder return code": return_code
        }
        
        print_validation_summary("Basic Scanning Algorithm", success, errors, details)
        
        return success
        
    finally:
        # Cleanup
        cleanup_temp_files([fasta_file])
        cleanup_temp_dirs([output_dir])


def test_scanning_algorithm_edge_positions() -> bool:
    """
    Test scanning algorithm with pdif sites at edge positions.
    """
    print("\nTesting scanning algorithm with edge positions...")
    
    # Create test sequence with pdif sites at very beginning and end
    sequence = create_sequence_with_pdif_sites(
        length=2000,
        pdif_sites=[
            ("ATTTAACATAAGGGCTGTTATACGAAAT", 0),     # At very beginning
            ("AGTACATATAACAAAGATTATGTTAAAT", 1972),  # At very end (2000-28)
        ]
    )
    
    # Save to temporary FASTA file
    fasta_file = save_sequence_to_fasta(
        sequence,
        "test_scan_edge",
        "Test sequence for edge position scanning validation"
    )
    
    # Create temporary output directory
    output_dir = tempfile.mkdtemp(prefix='pdif_scan_edge_')
    
    try:
        # Run pdifFinder
        return_code, stdout, stderr = run_pdif_finder(fasta_file, output_dir)
        
        if return_code != 0:
            print(f"pdifFinder failed with return code {return_code}")
            print(f"stderr: {stderr}")
            return False
        
        # Parse results
        results = parse_pdif_output(output_dir)
        pdif_sites = results['pdif_sites']
        
        # Expected positions (1-based for pdifFinder output)
        expected_positions = [
            (1, 28),       # position 0-27 -> 1-28
            (1973, 2000),  # position 1972-1999 -> 1973-2000
        ]
        
        # Validate
        success, errors = validate_pdif_positions(expected_positions, pdif_sites)
        
        # Print summary
        details = {
            "Sequence length": f"{len(sequence)} bp",
            "Expected pdif sites": len(expected_positions),
            "Found pdif sites": len(pdif_sites),
            "pdifFinder return code": return_code
        }
        
        print_validation_summary("Edge Position Scanning", success, errors, details)
        
        return success
        
    finally:
        # Cleanup
        cleanup_temp_files([fasta_file])
        cleanup_temp_dirs([output_dir])


def test_scanning_algorithm_close_sites() -> bool:
    """
    Test scanning algorithm with pdif sites close together.
    """
    print("\nTesting scanning algorithm with close pdif sites...")
    
    # Create test sequence with pdif sites very close together
    sequence = create_sequence_with_pdif_sites(
        length=1000,
        pdif_sites=[
            ("ATTTAACATAAGGGCTGTTATACGAAAT", 100),   # pdif1
            ("AGTACATATAACAAAGATTATGTTAAAT", 130),   # pdif2 (30bp after pdif1)
            ("AATTAAAATACCTTCTGTTATGTGCAAC", 160),   # pdif3 (30bp after pdif2)
        ]
    )
    
    # Save to temporary FASTA file
    fasta_file = save_sequence_to_fasta(
        sequence,
        "test_scan_close",
        "Test sequence for close pdif site scanning validation"
    )
    
    # Create temporary output directory
    output_dir = tempfile.mkdtemp(prefix='pdif_scan_close_')
    
    try:
        # Run pdifFinder
        return_code, stdout, stderr = run_pdif_finder(fasta_file, output_dir)
        
        if return_code != 0:
            print(f"pdifFinder failed with return code {return_code}")
            print(f"stderr: {stderr}")
            return False
        
        # Parse results
        results = parse_pdif_output(output_dir)
        pdif_sites = results['pdif_sites']
        
        # Expected positions (1-based for pdifFinder output)
        expected_positions = [
            (101, 128),   # position 100-127 -> 101-128
            (131, 158),   # position 130-157 -> 131-158
            (161, 188),   # position 160-187 -> 161-188
        ]
        
        # Validate
        success, errors = validate_pdif_positions(expected_positions, pdif_sites)
        
        # Print summary
        details = {
            "Sequence length": f"{len(sequence)} bp",
            "Expected pdif sites": len(expected_positions),
            "Found pdif sites": len(pdif_sites),
            "pdifFinder return code": return_code,
            "Site spacing": "30bp between sites"
        }
        
        print_validation_summary("Close Site Scanning", success, errors, details)
        
        return success
        
    finally:
        # Cleanup
        cleanup_temp_files([fasta_file])
        cleanup_temp_dirs([output_dir])


def test_scanning_algorithm_large_sequence() -> bool:
    """
    Test scanning algorithm with a large sequence.
    """
    print("\nTesting scanning algorithm with large sequence...")
    
    # Create a large test sequence (10k bp) with pdif sites
    sequence = create_sequence_with_pdif_sites(
        length=10000,
        pdif_sites=[
            ("ATTTAACATAAGGGCTGTTATACGAAAT", 1000),   # pdif1 at 1k
            ("AGTACATATAACAAAGATTATGTTAAAT", 5000),   # pdif2 at 5k
            ("AATTAAAATACCTTCTGTTATGTGCAAC", 9000),   # pdif3 at 9k
        ]
    )
    
    # Save to temporary FASTA file
    fasta_file = save_sequence_to_fasta(
        sequence,
        "test_scan_large",
        "Test sequence for large sequence scanning validation (10k bp)"
    )
    
    # Create temporary output directory
    output_dir = tempfile.mkdtemp(prefix='pdif_scan_large_')
    
    try:
        # Run pdifFinder
        return_code, stdout, stderr = run_pdif_finder(fasta_file, output_dir)
        
        if return_code != 0:
            print(f"pdifFinder failed with return code {return_code}")
            print(f"stderr: {stderr}")
            return False
        
        # Parse results
        results = parse_pdif_output(output_dir)
        pdif_sites = results['pdif_sites']
        
        # Expected positions (1-based for pdifFinder output)
        expected_positions = [
            (1001, 1028),   # position 1000-1027 -> 1001-1028
            (5001, 5028),   # position 5000-5027 -> 5001-5028
            (9001, 9028),   # position 9000-9027 -> 9001-9028
        ]
        
        # Validate
        success, errors = validate_pdif_positions(expected_positions, pdif_sites)
        
        # Print summary
        details = {
            "Sequence length": f"{len(sequence)} bp",
            "Expected pdif sites": len(expected_positions),
            "Found pdif sites": len(pdif_sites),
            "pdifFinder return code": return_code
        }
        
        print_validation_summary("Large Sequence Scanning", success, errors, details)
        
        return success
        
    finally:
        # Cleanup
        cleanup_temp_files([fasta_file])
        cleanup_temp_dirs([output_dir])


def run_all_scanning_tests() -> bool:
    """
    Run all scanning algorithm tests.
    
    Returns:
        True if all tests pass, False otherwise
    """
    print("="*60)
    print("VALIDATING SCANNING ALGORITHM")
    print("="*60)
    
    test_results = []
    
    # Run all tests
    test_results.append(("Basic Scanning", test_scanning_algorithm_basic()))
    test_results.append(("Edge Positions", test_scanning_algorithm_edge_positions()))
    test_results.append(("Close Sites", test_scanning_algorithm_close_sites()))
    test_results.append(("Large Sequence", test_scanning_algorithm_large_sequence()))
    
    # Print summary
    print("\n" + "="*60)
    print("SCANNING ALGORITHM VALIDATION SUMMARY")
    print("="*60)
    
    all_passed = True
    for test_name, passed in test_results:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{test_name}: {status}")
        if not passed:
            all_passed = False
    
    print("\n" + "="*60)
    if all_passed:
        print("✅ ALL SCANNING TESTS PASSED")
    else:
        print("❌ SOME SCANNING TESTS FAILED")
    print("="*60)
    
    return all_passed


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Validate pdifFinder scanning algorithm",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                    # Run all scanning tests
  %(prog)s --test basic       # Run only basic scanning test
  %(prog)s --test edge        # Run only edge position test
  %(prog)s --test close       # Run only close sites test
  %(prog)s --test large       # Run only large sequence test
        """
    )
    
    parser.add_argument(
        "--test",
        choices=["basic", "edge", "close", "large", "all"],
        default="all",
        help="Specific test to run (default: all)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose output"
    )
    
    args = parser.parse_args()
    
    # Set verbose mode
    if not args.verbose:
        # Redirect stdout to reduce output
        import io
        sys.stdout = io.StringIO()
    
    try:
        if args.test == "all":
            success = run_all_scanning_tests()
        elif args.test == "basic":
            success = test_scanning_algorithm_basic()
        elif args.test == "edge":
            success = test_scanning_algorithm_edge_positions()
        elif args.test == "close":
            success = test_scanning_algorithm_close_sites()
        elif args.test == "large":
            success = test_scanning_algorithm_large_sequence()
        else:
            print(f"Unknown test: {args.test}")
            return 1
        
        return 0 if success else 1
        
    except Exception as e:
        print(f"Error during validation: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())