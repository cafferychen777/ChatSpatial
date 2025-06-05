"""
Test spatial domain identification functionality.
"""

import pytest
import numpy as np
import pandas as pd
import scanpy as sc
from unittest.mock import AsyncMock, MagicMock

# Import the modules we want to test
from chatspatial.models.data import SpatialDomainParameters
from chatspatial.models.analysis import SpatialDomainResult
from chatspatial.tools.spatial_domains import identify_spatial_domains


@pytest.fixture
def mock_adata():
    """Create a mock AnnData object for testing"""
    # Create synthetic spatial transcriptomics data
    n_obs = 100
    n_vars = 50
    
    # Create expression matrix
    X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars)).astype(float)
    
    # Create AnnData object
    adata = sc.AnnData(X)
    
    # Add spatial coordinates
    spatial_coords = np.random.uniform(0, 10, size=(n_obs, 2))
    adata.obsm['spatial'] = spatial_coords
    
    # Add gene names
    adata.var_names = [f"Gene_{i}" for i in range(n_vars)]
    
    # Add cell barcodes
    adata.obs_names = [f"Cell_{i}" for i in range(n_obs)]
    
    # Add highly variable genes
    adata.var['highly_variable'] = np.random.choice([True, False], size=n_vars, p=[0.3, 0.7])
    
    return adata


@pytest.fixture
def mock_data_store(mock_adata):
    """Create a mock data store"""
    return {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test",
            "n_cells": mock_adata.n_obs,
            "n_genes": mock_adata.n_vars
        }
    }


@pytest.mark.asyncio
async def test_identify_spatial_domains_clustering():
    """Test spatial domain identification using clustering methods"""
    # Create mock data
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))
    mock_adata.var['highly_variable'] = np.random.choice([True, False], size=30, p=[0.5, 0.5])
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    # Test leiden clustering
    params = SpatialDomainParameters(
        method="leiden",
        n_domains=5,
        resolution=0.5,
        include_image=False  # Skip visualization for testing
    )
    
    # Mock context
    mock_context = AsyncMock()
    
    # Run the function
    result = await identify_spatial_domains("test_data", mock_data_store, params, mock_context)
    
    # Verify the result
    assert isinstance(result, SpatialDomainResult)
    assert result.data_id == "test_data"
    assert result.method == "leiden"
    assert result.n_domains > 0
    assert result.domain_key.startswith("spatial_domains_leiden")
    assert len(result.domain_counts) == result.n_domains
    
    # Verify that domain labels were added to the data
    assert result.domain_key in mock_data_store["test_data"]["adata"].obs.columns


@pytest.mark.asyncio
async def test_identify_spatial_domains_invalid_method():
    """Test error handling for invalid methods"""
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    # Test invalid method
    params = SpatialDomainParameters(
        method="invalid_method",  # This should cause an error
        n_domains=5
    )
    
    mock_context = AsyncMock()
    
    # This should raise an error
    with pytest.raises(ValueError, match="Unsupported method"):
        await identify_spatial_domains("test_data", mock_data_store, params, mock_context)


@pytest.mark.asyncio
async def test_identify_spatial_domains_no_spatial_coords():
    """Test error handling when no spatial coordinates are available"""
    # Create data without spatial coordinates
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    # Note: no spatial coordinates added
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    params = SpatialDomainParameters(method="leiden", n_domains=5)
    mock_context = AsyncMock()
    
    # This should raise an error
    with pytest.raises(ValueError, match="No spatial coordinates found"):
        await identify_spatial_domains("test_data", mock_data_store, params, mock_context)


@pytest.mark.asyncio
async def test_identify_spatial_domains_dataset_not_found():
    """Test error handling when dataset is not found"""
    mock_data_store = {}  # Empty data store
    
    params = SpatialDomainParameters(method="leiden", n_domains=5)
    mock_context = AsyncMock()
    
    # This should raise an error
    with pytest.raises(ValueError, match="Dataset .* not found"):
        await identify_spatial_domains("nonexistent_data", mock_data_store, params, mock_context)


def test_spatial_domain_parameters_validation():
    """Test parameter validation"""
    # Test valid parameters
    params = SpatialDomainParameters(
        method="stagate",
        n_domains=7,
        stagate_alpha=0.5,
        resolution=0.8
    )
    assert params.method == "stagate"
    assert params.n_domains == 7
    assert params.stagate_alpha == 0.5
    
    # Test invalid n_domains (should be > 0)
    with pytest.raises(ValueError):
        SpatialDomainParameters(n_domains=0)
    
    # Test invalid n_domains (should be <= 50)
    with pytest.raises(ValueError):
        SpatialDomainParameters(n_domains=100)
    
    # Test invalid alpha (should be between 0 and 1)
    with pytest.raises(ValueError):
        SpatialDomainParameters(stagate_alpha=1.5)


def test_spatial_domain_result_creation():
    """Test SpatialDomainResult creation"""
    result = SpatialDomainResult(
        data_id="test_data",
        method="leiden",
        n_domains=5,
        domain_key="spatial_domains_leiden",
        domain_counts={"0": 20, "1": 15, "2": 10, "3": 8, "4": 7},
        statistics={"method": "leiden", "resolution": 0.5}
    )
    
    assert result.data_id == "test_data"
    assert result.method == "leiden"
    assert result.n_domains == 5
    assert len(result.domain_counts) == 5
    assert sum(result.domain_counts.values()) == 60


@pytest.mark.asyncio
async def test_identify_spatial_domains_with_refinement():
    """Test spatial domain identification with refinement"""
    # Create mock data with spatial structure
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    
    # Create spatial coordinates with some structure
    x_coords = np.concatenate([
        np.random.normal(2, 0.5, 25),  # First cluster
        np.random.normal(8, 0.5, 25)   # Second cluster
    ])
    y_coords = np.concatenate([
        np.random.normal(2, 0.5, 25),  # First cluster
        np.random.normal(8, 0.5, 25)   # Second cluster
    ])
    mock_adata.obsm['spatial'] = np.column_stack([x_coords, y_coords])
    mock_adata.var['highly_variable'] = np.random.choice([True, False], size=30, p=[0.5, 0.5])
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    # Test with refinement enabled
    params = SpatialDomainParameters(
        method="leiden",
        n_domains=3,
        refine_domains=True,
        include_image=False
    )
    
    mock_context = AsyncMock()
    
    # Run the function
    result = await identify_spatial_domains("test_data", mock_data_store, params, mock_context)
    
    # Verify that both original and refined domains were created
    assert result.domain_key in mock_data_store["test_data"]["adata"].obs.columns
    assert result.refined_domain_key is not None
    assert result.refined_domain_key in mock_data_store["test_data"]["adata"].obs.columns


if __name__ == "__main__":
    pytest.main([__file__])
