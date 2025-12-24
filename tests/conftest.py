"""
pytest configuration and global fixtures

This file defines shared fixtures and configuration for all tests.
"""
import pytest
import numpy as np
from pathlib import Path


# ========== Pytest Configuration ==========

def pytest_configure(config):
    """Configure pytest markers"""
    config.addinivalue_line("markers", "unit: Unit tests (fast)")
    config.addinivalue_line("markers", "integration: Integration tests (slower)")
    config.addinivalue_line("markers", "slow: Slow tests (>5 seconds)")
    config.addinivalue_line("markers", "requires_r: Tests requiring R environment")
    config.addinivalue_line("markers", "requires_gpu: Tests requiring GPU")
    config.addinivalue_line("markers", "mcp: MCP protocol related tests")


# ========== Path Fixtures ==========

@pytest.fixture(scope="session")
def test_data_dir():
    """Return test data directory path"""
    return Path(__file__).parent / "fixtures" / "sample_data"


@pytest.fixture(scope="session")
def sample_visium_path(test_data_dir):
    """Return sample Visium data path"""
    path = test_data_dir / "visium_sample.h5ad"
    if not path.exists():
        pytest.skip(f"Test data not found: {path}. Run: python tests/fixtures/mock_adata.py")
    return str(path)


# ========== AnnData Fixtures ==========

@pytest.fixture
def mock_adata():
    """Create basic mock AnnData object (100 cells × 200 genes)"""
    from tests.fixtures.mock_adata import create_mock_adata
    return create_mock_adata(n_obs=100, n_vars=200)


@pytest.fixture
def mock_adata_with_spatial():
    """Create mock AnnData with spatial coordinates"""
    from tests.fixtures.mock_adata import create_mock_adata
    return create_mock_adata(n_obs=100, n_vars=200, add_spatial=True)


@pytest.fixture
def mock_adata_with_clusters():
    """Create mock AnnData with cluster labels"""
    from tests.fixtures.mock_adata import create_mock_adata
    return create_mock_adata(n_obs=100, n_vars=200, add_clusters=True)


@pytest.fixture
def mock_reference_adata():
    """Create reference dataset for annotation"""
    from tests.fixtures.mock_adata import create_mock_reference_adata
    return create_mock_reference_adata(n_obs=500, n_vars=200)


@pytest.fixture
def mock_velocity_adata():
    """Create AnnData with velocity layers"""
    from tests.fixtures.mock_adata import create_mock_velocity_adata
    return create_mock_velocity_adata(n_obs=100, n_vars=200)


# ========== Data Store Fixtures ==========

@pytest.fixture
def data_store(mock_adata):
    """Create test data_store dictionary"""
    return {"test_dataset": {"adata": mock_adata}}


@pytest.fixture
def data_store_with_spatial(mock_adata_with_spatial):
    """Create data_store with spatial data"""
    return {"test_dataset": {"adata": mock_adata_with_spatial}}


# ========== MCP Context Mock ==========

class MockContext:
    """Mock MCP Context for testing

    Records all info/warning/error messages for verification.
    """

    def __init__(self):
        self.messages = []
        self.info_messages = []
        self.warning_messages = []
        self.error_messages = []

    async def info(self, message: str):
        """Record info message"""
        self.messages.append(("info", message))
        self.info_messages.append(message)

    async def warning(self, message: str):
        """Record warning message"""
        self.messages.append(("warning", message))
        self.warning_messages.append(message)

    async def error(self, message: str):
        """Record error message"""
        self.messages.append(("error", message))
        self.error_messages.append(message)

    def get_messages(self, level=None):
        """Get messages list"""
        if level is None:
            return self.messages
        return [msg for lvl, msg in self.messages if lvl == level]

    def clear(self):
        """Clear all messages"""
        self.messages = []
        self.info_messages = []
        self.warning_messages = []
        self.error_messages = []


@pytest.fixture
def mock_context():
    """Create mock MCP context"""
    return MockContext()


# ========== Session-level Setup ==========

@pytest.fixture(scope="session", autouse=True)
def setup_test_environment():
    """Session-level test environment setup (runs automatically)"""
    # Set random seed for reproducibility
    np.random.seed(42)

    # Generate test datasets if they don't exist
    test_data_dir = Path(__file__).parent / "fixtures" / "sample_data"
    if not test_data_dir.exists() or not any(test_data_dir.iterdir()):
        print("\n📦 Generating test datasets...")
        from tests.fixtures.mock_adata import generate_sample_datasets
        generate_sample_datasets(str(test_data_dir))

    yield

    # Teardown (cleanup)
    pass


# ========== Utility Fixtures ==========

@pytest.fixture
def temp_output_dir(tmp_path):
    """Create temporary output directory"""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return str(output_dir)
