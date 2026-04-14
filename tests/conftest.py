"""
Pytest configuration file for PdifFinder tests.
This file is automatically loaded by pytest and provides shared fixtures and configuration.
"""

import os
import sys
import pytest

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


@pytest.fixture
def test_data_dir():
    """Fixture to provide path to test data directory."""
    return os.path.join(os.path.dirname(__file__), 'data')


@pytest.fixture
def sample_fasta_file():
    """Fixture to provide path to a sample FASTA file for testing."""
    return os.path.join(os.path.dirname(__file__), 'data', 'sample.fasta')


@pytest.fixture
def sample_genbank_file():
    """Fixture to provide path to a sample GenBank file for testing."""
    return os.path.join(os.path.dirname(__file__), 'data', 'sample.gb')