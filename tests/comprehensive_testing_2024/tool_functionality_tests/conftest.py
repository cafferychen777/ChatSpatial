"""
Pytest configuration for ChatSpatial tool functionality tests.
"""

import pytest
import warnings
import numpy as np
import scanpy as sc
from pathlib import Path
import tempfile
import shutil

# Configure scanpy settings for testing
sc.settings.verbosity = 0  # Suppress scanpy output during tests
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Configure numpy and warnings
np.random.seed(42)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

# Define test markers
def pytest_configure(config):
    """Configure pytest markers"""
    config.addinivalue_line("markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')")
    config.addinivalue_line("markers", "integration: marks tests as integration tests")
    config.addinivalue_line("markers", "requires_gpu: marks tests that require GPU")
    config.addinivalue_line("markers", "requires_r: marks tests that require R")
    config.addinivalue_line("markers", "memory_intensive: marks tests that use significant memory")


@pytest.fixture(scope="session")
def test_data_dir():
    """Create temporary directory for test data"""
    temp_dir = tempfile.mkdtemp(prefix="chatspatial_test_")
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)


@pytest.fixture(scope="session") 
def datasets_dir():
    """Path to the test datasets directory"""
    datasets_path = Path(__file__).parent.parent / "datasets"
    if not datasets_path.exists():
        pytest.skip(f"Test datasets directory not found: {datasets_path}")
    return datasets_path


@pytest.fixture(autouse=True)
def reset_scanpy_settings():
    """Reset scanpy settings before each test"""
    sc.settings.verbosity = 0
    yield
    # Cleanup after test
    sc.settings.verbosity = 0


@pytest.fixture(autouse=True)  
def suppress_warnings():
    """Suppress warnings during testing"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add markers based on test names"""
    for item in items:
        # Mark slow tests
        if "slow" in item.name or "performance" in item.name or "benchmark" in item.name:
            item.add_marker(pytest.mark.slow)
        
        # Mark integration tests
        if "integration" in item.name or "end_to_end" in item.name:
            item.add_marker(pytest.mark.integration)
        
        # Mark tests requiring specific dependencies
        if "gpu" in item.name:
            item.add_marker(pytest.mark.requires_gpu)
        
        if "r_" in item.name or "_r_" in item.name:
            item.add_marker(pytest.mark.requires_r)
            
        # Mark memory intensive tests
        if "large" in item.name or "memory" in item.name:
            item.add_marker(pytest.mark.memory_intensive)