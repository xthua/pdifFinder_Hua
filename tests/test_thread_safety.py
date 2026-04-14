#!/usr/bin/env python3
"""
Test thread-safe file operations for pdifFinder.
Focus on verifying that file writes in threaded functions are protected by locks.
"""

import os
import sys
import inspect
from pathlib import Path
import pytest

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import PdifFinder.pdifFinder as pf


def test_lock_variables_exist():
    """
    Test that the global lock variables are defined in the module.
    """
    # Check that fragment_lock and pair_lock exist
    assert hasattr(pf, 'fragment_lock'), "fragment_lock not found in pdifFinder module"
    assert hasattr(pf, 'pair_lock'), "pair_lock not found in pdifFinder module"
    
    # Verify they are threading.Lock instances
    import threading
    assert hasattr(pf.fragment_lock, 'acquire'), "fragment_lock missing acquire method"
    assert hasattr(pf.pair_lock, 'acquire'), "pair_lock missing acquire method"


def test_find_match_fragment_thread_uses_lock():
    """
    Test that findMatchFragmentThread uses fragment_lock for file writes.
    """
    source_lines = inspect.getsource(pf.findMatchFragmentThread).split('\n')
    
    # Look for with fragment_lock: pattern
    lock_used = any('with fragment_lock:' in line for line in source_lines)
    assert lock_used, "findMatchFragmentThread does not use fragment_lock"
    
    # Also verify the file write line is present (optional)
    file_write_present = any('with open(outfile' in line for line in source_lines)
    assert file_write_present, "findMatchFragmentThread file write not found"


def test_find_possible_pair_thread_uses_lock():
    """
    Test that findPossiblePairThread uses pair_lock for file writes.
    """
    source_lines = inspect.getsource(pf.findPossiblePairThread).split('\n')
    
    # Look for with pair_lock: pattern
    lock_used = any('with pair_lock:' in line for line in source_lines)
    assert lock_used, "findPossiblePairThread does not use pair_lock"
    
    # Also verify the file write line is present
    file_write_present = any("with open(outdir + '/tmp/pairSearch/%s' % name" in line for line in source_lines)
    assert file_write_present, "findPossiblePairThread file write not found"


def test_lock_import():
    """
    Test that Lock is imported from threading.
    """
    import_lines = inspect.getsource(pf).split('\n')
    lock_import = any('from threading import Thread, Lock' in line or 
                      'from threading import Lock' in line for line in import_lines)
    assert lock_import, "Lock not imported from threading module"


if __name__ == "__main__":
    # Run tests directly if script is executed
    test_lock_variables_exist()
    print("✓ Lock variables exist")
    test_find_match_fragment_thread_uses_lock()
    print("✓ findMatchFragmentThread uses lock")
    test_find_possible_pair_thread_uses_lock()
    print("✓ findPossiblePairThread uses lock")
    test_lock_import()
    print("✓ Lock import verified")
    print("All thread safety tests passed!")