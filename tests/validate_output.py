#!/usr/bin/env python3
"""
Validate pdifFinder output format and content.
This script compares pdifFinder output to expected results and checks output file formats.
"""

import argparse
import sys
import os
import tempfile
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional

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


def validate_output_files_exist(output_dir: str) -> Tuple[bool, List[str]]:
    """
    Check that all expected output files exist.
    
    Args:
        output_dir: Directory containing pdifFinder output
        
    Returns:
        Tuple of (success, list of missing files)
    """
    expected_files = [
        'pdif_site.txt',
        'AMRgene.txt',
        'pdifmodule_list.txt',
        'pdifmoduleseq.fasta',
        'pdifmodule.svg',
        'plasmid.html'
    ]
    
    missing_files = []
    for filename in expected_files:
        filepath = os.path.join(output_dir, filename)
        if not os.path.exists(filepath):
            missing_files.append(filename)
    
    return len(missing_files) == 0, missing_files


def validate_pdif_site_format(output_dir: str) -> Tuple[bool, List[str]]:
    """
    Validate pdif_site.txt file format.
    
    Args:
        output_dir: Directory containing pdifFinder output
        
    Returns:
        Tuple of (success, list of format errors)
    """
    errors = []
    pdif_site_file = os.path.join(output_dir, 'pdif_site.txt')
    
    if not os.path.exists(pdif_site_file):
        errors.append("pdif_site.txt file not found")
        return False, errors
    
    try:
        with open(pdif_site_file, 'r') as f:
            lines = f.readlines()
        
        # Check header
        if not lines or not lines[0].startswith('#'):
            errors.append("Missing header comment in pdif_site.txt")
        
        # Check data lines
        data_lines = [line for line in lines if line.strip() and not line.startswith('#')]
        for i, line in enumerate(data_lines, 1):
            parts = line.strip().split('\t')
            if len(parts) < 4:
                errors.append(f"Line {i}: Expected at least 4 tab-separated columns, got {len(parts)}")
                continue
            
            # Check position columns are integers
            try:
                start_pos = int(parts[1])
                end_pos = int(parts[2])
                
                if start_pos <= 0:
                    errors.append(f"Line {i}: Start position must be positive, got {start_pos}")
                if end_pos <= 0:
                    errors.append(f"Line {i}: End position must be positive, got {end_pos}")
                if start_pos >= end_pos:
                    errors.append(f"Line {i}: Start position ({start_pos}) must be less than end position ({end_pos})")
                
                # Check pdif site sequence length (should be 28bp)
                pdif_site = parts[3]
                if len(pdif_site) != 28:
                    errors.append(f"Line {i}: pdif site should be 28bp, got {len(pdif_site)}bp: {pdif_site}")
                
            except ValueError:
                errors.append(f"Line {i}: Position columns must be integers: {parts[1]}, {parts[2]}")
        
        return len(errors) == 0, errors
        
    except Exception as e:
        errors.append(f"Error reading pdif_site.txt: {e}")
        return False, errors


def validate_amr_gene_format(output_dir: str) -> Tuple[bool, List[str]]:
    """
    Validate AMRgene.txt file format.
    
    Args:
        output_dir: Directory containing pdifFinder output
        
    Returns:
        Tuple of (success, list of format errors)
    """
    errors = []
    amr_gene_file = os.path.join(output_dir, 'AMRgene.txt')
    
    if not os.path.exists(amr_gene_file):
        errors.append("AMRgene.txt file not found")
        return False, errors
    
    try:
        with open(amr_gene_file, 'r') as f:
            lines = f.readlines()
        
        # Check header
        if not lines or not lines[0].startswith('#'):
            errors.append("Missing header comment in AMRgene.txt")
        
        # Check data lines
        data_lines = [line for line in lines if line.strip() and not line.startswith('#')]
        for i, line in enumerate(data_lines, 1):
            parts = line.strip().split('\t')
            if len(parts) < 6:
                errors.append(f"Line {i}: Expected at least 6 tab-separated columns, got {len(parts)}")
                continue
            
            # Check position columns are integers
            try:
                start_pos = int(parts[1])
                end_pos = int(parts[2])
                
                if start_pos <= 0:
                    errors.append(f"Line {i}: Start position must be positive, got {start_pos}")
                if end_pos <= 0:
                    errors.append(f"Line {i}: End position must be positive, got {end_pos}")
                if start_pos >= end_pos:
                    errors.append(f"Line {i}: Start position ({start_pos}) must be less than end position ({end_pos})")
                
            except ValueError:
                errors.append(f"Line {i}: Position columns must be integers: {parts[1]}, {parts[2]}")
        
        return len(errors) == 0, errors
        
    except Exception as e:
        errors.append(f"Error reading AMRgene.txt: {e}")
        return False, errors


def validate_pdif_module_format(output_dir: str) -> Tuple[bool, List[str]]:
    """
    Validate pdifmodule_list.txt file format.
    
    Args:
        output_dir: Directory containing pdifFinder output
        
    Returns:
        Tuple of (success, list of format errors)
    """
    errors = []
    pdif_module_file = os.path.join(output_dir, 'pdifmodule_list.txt')
    
    if not os.path.exists(pdif_module_file):
        errors.append("pdifmodule_list.txt file not found")
        return False, errors
    
    try:
        with open(pdif_module_file, 'r') as f:
            lines = f.readlines()
        
        # Check header
        if not lines or not lines[0].startswith('#'):
            errors.append("Missing header comment in pdifmodule_list.txt")
        
        # Check data lines
        data_lines = [line for line in lines if line.strip() and not line.startswith('#')]
        for i, line in enumerate(data_lines, 1):
            parts = line.strip().split('\t')
            if len(parts) < 5:
                errors.append(f"Line {i}: Expected at least 5 tab-separated columns, got {len(parts)}")
                continue
            
            # Check position columns are integers
            try:
                start_pos = int(parts[1])
                end_pos = int(parts[2])
                
                if start_pos <= 0:
                    errors.append(f"Line {i}: Start position must be positive, got {start_pos}")
                if end_pos <= 0:
                    errors.append(f"Line {i}: End position must be positive, got {end_pos}")
                if start_pos >= end_pos:
                    errors.append(f"Line {i}: Start position ({start_pos}) must be less than end position ({end_pos})")
                
            except ValueError:
                errors.append(f"Line {i}: Position columns must be integers: {parts[1]}, {parts[2]}")
        
        return len(errors) == 0, errors
        
    except Exception as e:
        errors.append(f"Error reading pdifmodule_list.txt: {e}")
        return False, errors


def validate_sequence_files(output_dir: str) -> Tuple[bool, List[str]]:
    """
    Validate sequence output files (FASTA format).
    
    Args:
        output_dir: Directory containing pdifFinder output
        
    Returns:
        Tuple of (success, list of format errors)
    """
    errors = []
    
    # Check pdifmoduleseq.fasta
    fasta_file = os.path.join(output_dir, 'pdifmoduleseq.fasta')
    if os.path.exists(fasta_file):
        try:
            from Bio import SeqIO
            records = list(SeqIO.parse(fasta_file, "fasta"))
            if not records:
                errors.append("pdifmoduleseq.fasta is empty or invalid FASTA format")
        except Exception as e:
            errors.append(f"Error parsing pdifmoduleseq.fasta: {e}")
    else:
        errors.append("pdifmoduleseq.fasta file not found")
    
    return len(errors) == 0, errors


def validate_graphics_files(output_dir: str) -> Tuple[bool, List[str]]:
    """
    Validate graphics output files exist and are non-empty.
    
    Args:
        output_dir: Directory containing pdifFinder output
        
    Returns:
        Tuple of (success, list of errors)
    """
    errors = []
    
    # Check SVG file
    svg_file = os.path.join(output_dir, 'pdifmodule.svg')
    if os.path.exists(svg_file):
        if os.path.getsize(svg_file) == 0:
            errors.append("pdifmodule.svg file is empty")
    else:
        errors.append("pdifmodule.svg file not found")
    
    # Check HTML file
    html_file = os.path.join(output_dir, 'plasmid.html')
    if os.path.exists(html_file):
        if os.path.getsize(html_file) == 0:
            errors.append("plasmid.html file is empty")
        # Check it's actually HTML
        with open(html_file, 'r') as f:
            content = f.read(100)  # Read first 100 chars
            if '<html' not in content.lower() and '<!doctype' not in content.lower():
                errors.append("plasmid.html doesn't appear to be valid HTML")
    else:
        errors.append("plasmid.html file not found")
    
    return len(errors) == 0, errors


def test_output_format_validation() -> bool:
    """
    Test output format validation with a simple sequence.
    """
    print("Testing output format validation...")
    
    # Create a simple test sequence
    sequence = create_sequence_with_pdif_sites(
        length=1000,
        pdif_sites=[
            ("ATTTAACATAAGGGCTGTTATACGAAAT", 100),
        ]
    )
    
    # Save to temporary FASTA file
    fasta_file = save_sequence_to_fasta(
        sequence,
        "test_output_format",
        "Test sequence for output format validation"
    )
    
    # Create temporary output directory
    output_dir = tempfile.mkdtemp(prefix='pdif_output_test_')
    
    try:
        # Run pdifFinder
        return_code, stdout, stderr = run_pdif_finder(fasta_file, output_dir)
        
        if return_code != 0:
            print(f"pdifFinder failed with return code {return_code}")
            print(f"stderr: {stderr}")
            return False
        
        # Run all validations
        all_errors = []
        
        # 1. Check all files exist
        files_exist, missing_files = validate_output_files_exist(output_dir)
        if not files_exist:
            all_errors.append(f"Missing files: {', '.join(missing_files)}")
        
        # 2. Validate pdif_site.txt format
        pdif_format_ok, pdif_errors = validate_pdif_site_format(output_dir)
        if not pdif_format_ok:
            all_errors.extend([f"pdif_site.txt: {e}" for e in pdif_errors])
        
        # 3. Validate AMRgene.txt format
        amr_format_ok, amr_errors = validate_amr_gene_format(output_dir)
        if not amr_format_ok:
            all_errors.extend([f"AMRgene.txt: {e}" for e in amr_errors])
        
        # 4. Validate pdifmodule_list.txt format
        module_format_ok, module_errors = validate_pdif_module_format(output_dir)
        if not module_format_ok:
            all_errors.extend([f"pdifmodule_list.txt: {e}" for e in module_errors])
        
        # 5. Validate sequence files
        seq_files_ok, seq_errors = validate_sequence_files(output_dir)
        if not seq_files_ok:
            all_errors.extend([f"Sequence files: {e}" for e in seq_errors])
        
        # 6. Validate graphics files
        graphics_ok, graphics_errors = validate_graphics_files(output_dir)
        if not graphics_ok:
            all_errors.extend([f"Graphics files: {e}" for e in graphics_errors])
        
        # Print summary
        success = len(all_errors) == 0
        
        details = {
            "Sequence length": f"{len(sequence)} bp",
            "Output directory": output_dir,
            "pdifFinder return code": return_code,
            "Files checked": "pdif_site.txt, AMRgene.txt, pdifmodule_list.txt, pdifmoduleseq.fasta, pdifmodule.svg, plasmid.html"
        }
        
        print_validation_summary("Output Format Validation", success, all_errors, details)
        
        return success
        
    finally:
        # Cleanup
        cleanup_temp_files([fasta_file])
        cleanup_temp_dirs([output_dir])


def test_output_content_validation() -> bool:
    """
    Test output content validation with known pdif sites.
    """
    print("\nTesting output content validation...")
    
    # Create test sequence with known pdif sites
    sequence = create_sequence_with_pdif_sites(
        length=1500,
        pdif_sites=[
            ("ATTTAACATAAGGGCTGTTATACGAAAT", 100),
            ("AGTACATATAACAAAGATTATGTTAAAT", 500),
        ]
    )
    
    # Save to temporary FASTA file
    fasta_file = save_sequence_to_fasta(
        sequence,
        "test_output_content",
        "Test sequence for output content validation"
    )
    
    # Create temporary output directory
    output_dir = tempfile.mkdtemp(prefix='pdif_content_test_')
    
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
            (501, 528),   # position 500-527 -> 501-528
        ]
        
        # Validate content
        success, errors = validate_pdif_positions(expected_positions, pdif_sites)
        
        # Also check the actual pdif site sequences match
        if success and pdif_sites:
            expected_sites = ["ATTTAACATAAGGGCTGTTATACGAAAT", "AGTACATATAACAAAGATTATGTTAAAT"]
            for i, (result, expected) in enumerate(zip(pdif_sites, expected_sites)):
                if result['pdif_site'] != expected:
                    errors.append(f"Pdif site {i+1} sequence mismatch: expected {expected}, got {result['pdif_site']}")
                    success = False
        
        # Print summary
        details = {
            "Sequence length": f"{len(sequence)} bp",
            "Expected pdif sites": len(expected_positions),
            "Found pdif sites": len(pdif_sites),
            "pdifFinder return code": return_code
        }
        
        print_validation_summary("Output Content Validation", success, errors, details)
        
        return success
        
    finally:
        # Cleanup
        cleanup_temp_files([fasta_file])
        cleanup_temp_dirs([output_dir])


def test_empty_output_validation() -> bool:
    """
    Test output validation with sequence that should produce no results.
    """
    print("\nTesting empty output validation...")
    
    # Create sequence with no pdif sites (just random)
    from helpers import generate_test_sequence
    sequence = generate_test_sequence(length=1000, seed=123)
    
    # Save to temporary FASTA file
    fasta_file = save_sequence_to_fasta(
        sequence,
        "test_empty_output",
        "Test sequence that should produce no pdif sites"
    )
    
    # Create temporary output directory
    output_dir = tempfile.mkdtemp(prefix='pdif_empty_test_')
    
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
        
        # Should have no pdif sites
        success = len(pdif_sites) == 0
        errors = []
        
        if not success:
            errors.append(f"Expected no pdif sites, but found {len(pdif_sites)}")
        
        # Print summary
        details = {
            "Sequence length": f"{len(sequence)} bp",
            "Expected pdif sites": 0,
            "Found pdif sites": len(pdif_sites),
            "pdifFinder return code": return_code
        }
        
        print_validation_summary("Empty Output Validation", success, errors, details)
        
        return success
        
    finally:
        # Cleanup
        cleanup_temp_files([fasta_file])
        cleanup_temp_dirs([output_dir])


def run_all_output_tests() -> bool:
    """
    Run all output validation tests.
    
    Returns:
        True if all tests pass, False otherwise
    """
    print("="*60)
    print("VALIDATING OUTPUT FORMAT AND CONTENT")
    print("="*60)
    
    test_results = []
    
    # Run all tests
    test_results.append(("Output Format", test_output_format_validation()))
    test_results.append(("Output Content", test_output_content_validation()))
    test_results.append(("Empty Output", test_empty_output_validation()))
    
    # Print summary
    print("\n" + "="*60)
    print("OUTPUT VALIDATION SUMMARY")
    print("="*60)
    
    all_passed = True
    for test_name, passed in test_results:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{test_name}: {status}")
        if not passed:
            all_passed = False
    
    print("\n" + "="*60)
    if all_passed:
        print("✅ ALL OUTPUT TESTS PASSED")
    else:
        print("❌ SOME OUTPUT TESTS FAILED")
    print("="*60)
    
    return all_passed


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Validate pdifFinder output format and content",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                    # Run all output tests
  %(prog)s --test format      # Run only output format test
  %(prog)s --test content     # Run only output content test
  %(prog)s --test empty       # Run only empty output test
        """
    )
    
    parser.add_argument(
        "--test",
        choices=["format", "content", "empty", "all"],
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
            success = run_all_output_tests()
        elif args.test == "format":
            success = test_output_format_validation()
        elif args.test == "content":
            success = test_output_content_validation()
        elif args.test == "empty":
            success = test_empty_output_validation()
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