"""Test file to verify pytest framework is working."""

def test_pytest_working():
    """Simple test to verify pytest is working."""
    assert True


def test_import_pdif_finder():
    """Test that PdifFinder can be imported."""
    try:
        import PdifFinder
        assert True
    except ImportError:
        assert False, "Failed to import PdifFinder"