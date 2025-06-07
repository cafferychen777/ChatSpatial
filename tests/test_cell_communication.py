"""
Tests for LIANA+ cell communication analysis functionality
"""

import pytest
import numpy as np
import pandas as pd
import scanpy as sc
from unittest.mock import AsyncMock, patch, MagicMock

from chatspatial.models.data import CellCommunicationParameters
from chatspatial.tools.cell_communication import analyze_cell_communication
from chatspatial.models.analysis import CellCommunicationResult


@pytest.fixture
def mock_adata():
    """Create a mock AnnData object for testing"""
    # Create synthetic spatial transcriptomics data
    n_obs = 100
    n_vars = 50
    
    # Create expression matrix
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars)).astype(float)
    
    # Create AnnData object
    adata = sc.AnnData(X)
    
    # Add spatial coordinates
    spatial_coords = np.random.uniform(0, 10, size=(n_obs, 2))
    adata.obsm['spatial'] = spatial_coords
    
    # Add gene names (include some common ligands and receptors)
    gene_names = [
        'TGFB1', 'TGFBR1', 'TGFBR2', 'PDGFA', 'PDGFRA',
        'VEGFA', 'VEGFR1', 'VEGFR2', 'TNF', 'TNFRSF1A'
    ] + [f"Gene_{i:03d}" for i in range(10, n_vars)]
    adata.var_names = gene_names[:n_vars]
    
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
async def test_analyze_cell_communication_liana():
    """Test cell communication analysis with mocked LIANA+"""
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))

    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }

    # Mock LIANA+ functions
    with patch('chatspatial.tools.cell_communication.LIANA_AVAILABLE', True), \
         patch('chatspatial.tools.cell_communication.li') as mock_li:

        # Mock LIANA+ spatial neighbors
        mock_li.ut.spatial_neighbors.return_value = None

        # Mock LIANA+ bivariate analysis
        mock_lrdata = MagicMock()
        mock_lrdata.shape = (50, 10)  # 50 spots, 10 LR pairs
        mock_lrdata.var = pd.DataFrame({
            'morans': [0.5, 0.3, 0.7, 0.2, 0.8, 0.1, 0.6, 0.4, 0.9, 0.05],
            'cosine': [0.4, 0.2, 0.6, 0.1, 0.7, 0.05, 0.5, 0.3, 0.8, 0.02]
        }, index=[f'LR_pair_{i}' for i in range(10)])

        mock_li.mt.bivariate.return_value = mock_lrdata

        params = CellCommunicationParameters(
            method="liana",
            species="human",
            liana_local_metric="cosine",
            liana_global_metric="morans",
            liana_n_perms=50,
            include_image=False
        )

        mock_context = AsyncMock()

        result = await analyze_cell_communication("test_data", mock_data_store, params, mock_context)

        assert isinstance(result, CellCommunicationResult)
        assert result.method == "liana"
        assert result.species == "human"
        assert result.analysis_type == "spatial"
        assert result.n_lr_pairs == 10





@pytest.mark.asyncio
async def test_analyze_cell_communication_no_liana():
    """Test error handling when LIANA+ is not available"""
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))

    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }

    # Mock LIANA+ as not available
    with patch('chatspatial.tools.cell_communication.LIANA_AVAILABLE', False):
        params = CellCommunicationParameters(method="liana")
        mock_context = AsyncMock()

        # This should raise an ImportError
        with pytest.raises(ImportError, match="LIANA+ is not installed"):
            await analyze_cell_communication("test_data", mock_data_store, params, mock_context)


@pytest.mark.asyncio
async def test_analyze_cell_communication_invalid_method():
    """Test error handling for invalid method"""
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))

    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }

    # Test with invalid method - this should fail at parameter validation level
    with pytest.raises(ValueError):
        params = CellCommunicationParameters(method="invalid_method")


def test_cell_communication_parameters():
    """Test CellCommunicationParameters validation"""
    # Test valid parameters
    params = CellCommunicationParameters(
        method="liana",
        species="human",
        liana_local_metric="cosine",
        liana_global_metric="morans",
        liana_n_perms=100,
        min_cells=5
    )

    assert params.method == "liana"
    assert params.liana_local_metric == "cosine"
    assert params.liana_global_metric == "morans"

    # Test invalid global metric
    with pytest.raises(ValueError):
        CellCommunicationParameters(liana_global_metric="geary")  # Not supported


def test_cell_communication_result():
    """Test CellCommunicationResult creation"""
    result = CellCommunicationResult(
        data_id="test",
        method="liana",
        species="human",
        database="liana",
        n_lr_pairs=10,
        n_significant_pairs=5,
        analysis_type="spatial",
        statistics={"method": "liana", "local_metric": "cosine", "global_metric": "morans"}
    )

    assert result.method == "liana"
    assert result.n_lr_pairs == 10
    assert result.n_significant_pairs == 5
    assert result.analysis_type == "spatial"


@pytest.mark.asyncio
async def test_analyze_cell_communication_missing_data():
    """Test error handling for missing data"""
    mock_data_store = {}

    params = CellCommunicationParameters(method="liana")
    mock_context = AsyncMock()

    with pytest.raises(ValueError, match="Dataset test_data not found"):
        await analyze_cell_communication("test_data", mock_data_store, params, mock_context)


if __name__ == "__main__":
    pytest.main([__file__])
