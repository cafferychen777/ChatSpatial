"""
Test spatial variable genes identification functionality.
"""

import pytest
import numpy as np
import pandas as pd
import scanpy as sc
from unittest.mock import AsyncMock, MagicMock, patch

# Import the modules we want to test
from chatspatial.models.data import SpatialVariableGenesParameters
from chatspatial.models.analysis import SpatialVariableGenesResult
from chatspatial.tools.spatial_variable_genes import identify_spatial_variable_genes


@pytest.fixture
def mock_adata():
    """Create a mock AnnData object for testing"""
    # Create synthetic spatial transcriptomics data
    n_obs = 100
    n_vars = 50
    
    # Create expression matrix with some spatial structure
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars)).astype(float)
    
    # Add some spatial structure to first few genes
    spatial_coords = np.random.uniform(0, 10, size=(n_obs, 2))
    for i in range(5):  # First 5 genes have spatial patterns
        # Create a spatial gradient
        spatial_effect = np.exp(-((spatial_coords[:, 0] - 5)**2 + (spatial_coords[:, 1] - 5)**2) / 4)
        X[:, i] = X[:, i] * (1 + 2 * spatial_effect)
    
    # Create AnnData object
    adata = sc.AnnData(X)
    adata.obsm['spatial'] = spatial_coords
    
    # Add gene names
    adata.var_names = [f"Gene_{i:03d}" for i in range(n_vars)]
    
    # Add cell barcodes
    adata.obs_names = [f"Spot_{i:04d}" for i in range(n_obs)]
    
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
async def test_identify_spatial_variable_genes_mock():
    """Test spatial variable genes identification with mocked SpatialDE"""
    # Create mock data
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    # Mock SpatialDE results
    mock_spatialde_results = pd.DataFrame({
        'g': [f'Gene_{i:03d}' for i in range(10)],
        'pval': np.random.uniform(0.001, 0.1, 10),
        'qval': np.random.uniform(0.001, 0.1, 10),
        'l': np.random.uniform(0.5, 2.0, 10),
        'FSV': np.random.uniform(0.1, 0.8, 10)
    })
    mock_spatialde_results = mock_spatialde_results.sort_values('qval')
    
    # Mock the SpatialDE functions
    with patch('chatspatial.tools.spatial_variable_genes.SPATIALDE_AVAILABLE', True), \
         patch('chatspatial.tools.spatial_variable_genes.spd') as mock_spd, \
         patch('chatspatial.tools.spatial_variable_genes.NaiveDE') as mock_naive:
        
        # Setup mocks
        mock_spd.run.return_value = mock_spatialde_results
        mock_naive.stabilize.return_value = pd.DataFrame(np.random.randn(30, 50))
        mock_naive.regress_out.return_value = pd.DataFrame(np.random.randn(30, 50))
        
        # Test parameters
        params = SpatialVariableGenesParameters(
            method="spatialde",
            n_top_genes=5,
            significance_threshold=0.05,
            include_image=False  # Skip visualization for testing
        )
        
        # Mock context
        mock_context = AsyncMock()
        
        # Run the function
        result = await identify_spatial_variable_genes("test_data", mock_data_store, params, mock_context)
        
        # Verify the result
        assert isinstance(result, SpatialVariableGenesResult)
        assert result.data_id == "test_data"
        assert result.method == "spatialde"
        assert result.n_tested_genes == len(mock_spatialde_results)
        assert len(result.top_genes) <= params.n_top_genes
        assert result.results_key.startswith("spatial_variable_genes_spatialde")


@pytest.mark.asyncio
async def test_identify_spatial_variable_genes_no_spatialde():
    """Test error handling when SpatialDE is not available"""
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    # Mock SpatialDE as not available
    with patch('chatspatial.tools.spatial_variable_genes.SPATIALDE_AVAILABLE', False):
        params = SpatialVariableGenesParameters(method="spatialde")
        mock_context = AsyncMock()
        
        # This should raise an ImportError
        with pytest.raises(RuntimeError, match="SpatialDE is not installed"):
            await identify_spatial_variable_genes("test_data", mock_data_store, params, mock_context)


@pytest.mark.asyncio
async def test_identify_spatial_variable_genes_no_spatial_coords():
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
    
    params = SpatialVariableGenesParameters(method="spatialde")
    mock_context = AsyncMock()
    
    # This should raise an error
    with pytest.raises(RuntimeError, match="No spatial coordinates found"):
        await identify_spatial_variable_genes("test_data", mock_data_store, params, mock_context)


@pytest.mark.asyncio
async def test_identify_spatial_variable_genes_dataset_not_found():
    """Test error handling when dataset is not found"""
    mock_data_store = {}  # Empty data store
    
    params = SpatialVariableGenesParameters(method="spatialde")
    mock_context = AsyncMock()
    
    # This should raise an error
    with pytest.raises(ValueError, match="Dataset .* not found"):
        await identify_spatial_variable_genes("nonexistent_data", mock_data_store, params, mock_context)


@pytest.mark.asyncio
async def test_identify_spatial_variable_genes_with_aeh():
    """Test spatial variable genes identification with AEH"""
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    # Mock SpatialDE results with significant genes
    mock_spatialde_results = pd.DataFrame({
        'g': [f'Gene_{i:03d}' for i in range(10)],
        'pval': [0.001, 0.002, 0.003, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5],
        'qval': [0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6],
        'l': np.random.uniform(0.5, 2.0, 10),
        'FSV': np.random.uniform(0.1, 0.8, 10)
    })
    
    # Mock AEH results
    mock_aeh_results = pd.DataFrame({
        'g': [f'Gene_{i:03d}' for i in range(5)],  # First 5 genes are significant
        'pattern': [0, 0, 1, 1, 2],
        'membership': [0.9, 0.8, 0.9, 0.7, 0.8]
    })
    
    mock_patterns = np.random.randn(50, 3)  # 3 patterns
    
    # Mock the SpatialDE functions
    with patch('chatspatial.tools.spatial_variable_genes.SPATIALDE_AVAILABLE', True), \
         patch('chatspatial.tools.spatial_variable_genes.spd') as mock_spd, \
         patch('chatspatial.tools.spatial_variable_genes.NaiveDE') as mock_naive:
        
        # Setup mocks
        mock_spd.run.return_value = mock_spatialde_results
        mock_spd.aeh.spatial_patterns.return_value = (mock_aeh_results, mock_patterns)
        mock_naive.stabilize.return_value = pd.DataFrame(np.random.randn(30, 50))
        mock_naive.regress_out.return_value = pd.DataFrame(np.random.randn(30, 50))
        
        # Test parameters with AEH enabled
        params = SpatialVariableGenesParameters(
            method="spatialde",
            n_top_genes=5,
            significance_threshold=0.05,
            perform_aeh=True,
            aeh_n_patterns=3,
            include_image=False
        )
        
        mock_context = AsyncMock()
        
        # Run the function
        result = await identify_spatial_variable_genes("test_data", mock_data_store, params, mock_context)
        
        # Verify AEH was performed
        assert result.aeh_performed is True
        assert result.n_patterns == 3
        assert result.aeh_patterns_key is not None
        assert result.aeh_membership_key is not None


def test_spatial_variable_genes_parameters_validation():
    """Test parameter validation"""
    # Test valid parameters
    params = SpatialVariableGenesParameters(
        method="spatialde",
        n_top_genes=100,
        significance_threshold=0.05,
        perform_aeh=True,
        aeh_n_patterns=5
    )
    assert params.method == "spatialde"
    assert params.n_top_genes == 100
    assert params.significance_threshold == 0.05
    
    # Test invalid n_top_genes (should be > 0)
    with pytest.raises(ValueError):
        SpatialVariableGenesParameters(n_top_genes=0)
    
    # Test invalid significance_threshold (should be between 0 and 1)
    with pytest.raises(ValueError):
        SpatialVariableGenesParameters(significance_threshold=1.5)
    
    # Test invalid aeh_n_patterns (should be > 0)
    with pytest.raises(ValueError):
        SpatialVariableGenesParameters(aeh_n_patterns=0)


def test_spatial_variable_genes_result_creation():
    """Test SpatialVariableGenesResult creation"""
    result = SpatialVariableGenesResult(
        data_id="test_data",
        method="spatialde",
        n_significant_genes=25,
        n_tested_genes=100,
        significance_threshold=0.05,
        top_genes=["Gene_001", "Gene_002", "Gene_003"],
        results_key="spatial_variable_genes_spatialde",
        statistics={"method": "spatialde", "min_qval": 0.001}
    )
    
    assert result.data_id == "test_data"
    assert result.method == "spatialde"
    assert result.n_significant_genes == 25
    assert result.n_tested_genes == 100
    assert len(result.top_genes) == 3
    assert result.top_genes[0] == "Gene_001"


@pytest.mark.asyncio
async def test_unsupported_method():
    """Test error handling for unsupported methods"""
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    # Test unsupported method
    params = SpatialVariableGenesParameters(method="spark")  # SPARK not implemented yet
    mock_context = AsyncMock()
    
    # This should raise a NotImplementedError
    with pytest.raises(RuntimeError, match="not yet implemented"):
        await identify_spatial_variable_genes("test_data", mock_data_store, params, mock_context)


if __name__ == "__main__":
    pytest.main([__file__])
