"""
Test cell communication analysis functionality.
"""

import pytest
import numpy as np
import pandas as pd
import scanpy as sc
from unittest.mock import AsyncMock, MagicMock, patch

# Import the modules we want to test
from chatspatial.models.data import CellCommunicationParameters
from chatspatial.models.analysis import CellCommunicationResult
from chatspatial.tools.cell_communication import analyze_cell_communication


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
async def test_analyze_cell_communication_mock_commot():
    """Test cell communication analysis with mocked COMMOT"""
    # Create mock data
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))
    mock_adata.var_names = [f'Gene_{i}' for i in range(30)]
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    # Mock COMMOT functions
    mock_lr_db = pd.DataFrame({
        'ligand': ['TGFB1', 'PDGFA', 'VEGFA'],
        'receptor': ['TGFBR1_TGFBR2', 'PDGFRA', 'VEGFR1_VEGFR2'],
        'pathway': ['TGF_pathway', 'PDGF_pathway', 'VEGF_pathway']
    })
    
    with patch('chatspatial.tools.cell_communication.COMMOT_AVAILABLE', True), \
         patch('chatspatial.tools.cell_communication.ct') as mock_ct:
        
        # Setup mocks
        mock_ct.pp.ligand_receptor_database.return_value = mock_lr_db
        mock_ct.tl.spatial_communication.return_value = None
        
        # Mock the communication results
        mock_sender_data = np.random.uniform(0, 1, size=(50, 3))
        mock_receiver_data = np.random.uniform(0, 1, size=(50, 3))
        
        def mock_spatial_communication(adata, **kwargs):
            adata.obsm['commot-cellchat-sum-sender'] = mock_sender_data
            adata.obsm['commot-cellchat-sum-receiver'] = mock_receiver_data
        
        mock_ct.tl.spatial_communication.side_effect = mock_spatial_communication
        
        # Test parameters
        params = CellCommunicationParameters(
            method="commot",
            species="human",
            database="cellchat",
            include_image=False  # Skip visualization for testing
        )
        
        # Mock context
        mock_context = AsyncMock()
        
        # Run the function
        result = await analyze_cell_communication("test_data", mock_data_store, params, mock_context)
        
        # Verify the result
        assert isinstance(result, CellCommunicationResult)
        assert result.data_id == "test_data"
        assert result.method == "commot"
        assert result.species == "human"
        assert result.database == "cellchat"
        assert result.n_lr_pairs == len(mock_lr_db)
        assert result.commot_sender_key is not None
        assert result.commot_receiver_key is not None


@pytest.mark.asyncio
async def test_analyze_cell_communication_mock_spatialdm():
    """Test cell communication analysis with mocked SpatialDM"""
    # Create mock data
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))
    mock_adata.var_names = [f'Gene_{i}' for i in range(30)]
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    # Mock SpatialDM results
    mock_global_res = pd.DataFrame({
        'selected': [True, True, False, False, True],
        'z_pval': [0.01, 0.02, 0.15, 0.20, 0.03],
        'perm_pval': [0.02, 0.03, 0.18, 0.25, 0.04]
    }, index=['LR_pair_1', 'LR_pair_2', 'LR_pair_3', 'LR_pair_4', 'LR_pair_5'])
    
    with patch('chatspatial.tools.cell_communication.SPATIALDM_AVAILABLE', True), \
         patch('chatspatial.tools.cell_communication.sdm') as mock_sdm:
        
        # Setup mocks
        def mock_extract_lr(adata, species, min_cell):
            adata.uns['num_pairs'] = 5
        
        def mock_spatialdm_global(adata, **kwargs):
            adata.uns['global_res'] = mock_global_res
        
        mock_sdm.weight_matrix.return_value = None
        mock_sdm.extract_lr.side_effect = mock_extract_lr
        mock_sdm.spatialdm_global.side_effect = mock_spatialdm_global
        mock_sdm.sig_pairs.return_value = None
        mock_sdm.spatialdm_local.return_value = None
        mock_sdm.sig_spots.return_value = None
        
        # Test parameters
        params = CellCommunicationParameters(
            method="spatialdm",
            species="human",
            perform_global_analysis=True,
            perform_local_analysis=False,
            include_image=False
        )
        
        # Mock context
        mock_context = AsyncMock()
        
        # Run the function
        result = await analyze_cell_communication("test_data", mock_data_store, params, mock_context)
        
        # Verify the result
        assert isinstance(result, CellCommunicationResult)
        assert result.data_id == "test_data"
        assert result.method == "spatialdm"
        assert result.species == "human"
        assert result.n_lr_pairs == 5
        assert result.n_significant_pairs == 3  # 3 selected pairs
        assert result.global_results_key == "global_res"
        assert len(result.top_lr_pairs) > 0


@pytest.mark.asyncio
async def test_analyze_cell_communication_no_commot():
    """Test error handling when COMMOT is not available"""
    mock_adata = sc.AnnData(np.random.negative_binomial(5, 0.3, size=(50, 30)).astype(float))
    mock_adata.obsm['spatial'] = np.random.uniform(0, 10, size=(50, 2))
    
    mock_data_store = {
        "test_data": {
            "adata": mock_adata,
            "name": "Test Dataset",
            "type": "test"
        }
    }
    
    # Mock COMMOT as not available
    with patch('chatspatial.tools.cell_communication.COMMOT_AVAILABLE', False):
        params = CellCommunicationParameters(method="commot")
        mock_context = AsyncMock()
        
        # This should raise an ImportError
        with pytest.raises(RuntimeError, match="COMMOT is not installed"):
            await analyze_cell_communication("test_data", mock_data_store, params, mock_context)


@pytest.mark.asyncio
async def test_analyze_cell_communication_no_spatial_coords():
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
    
    params = CellCommunicationParameters(method="commot")
    mock_context = AsyncMock()
    
    # This should raise an error
    with pytest.raises(RuntimeError, match="No spatial coordinates found"):
        await analyze_cell_communication("test_data", mock_data_store, params, mock_context)


@pytest.mark.asyncio
async def test_analyze_cell_communication_dataset_not_found():
    """Test error handling when dataset is not found"""
    mock_data_store = {}  # Empty data store
    
    params = CellCommunicationParameters(method="commot")
    mock_context = AsyncMock()
    
    # This should raise an error
    with pytest.raises(ValueError, match="Dataset .* not found"):
        await analyze_cell_communication("nonexistent_data", mock_data_store, params, mock_context)


def test_cell_communication_parameters_validation():
    """Test parameter validation"""
    # Test valid parameters
    params = CellCommunicationParameters(
        method="commot",
        species="human",
        database="cellchat",
        commot_dis_thr=200.0,
        min_cells=3
    )
    assert params.method == "commot"
    assert params.species == "human"
    assert params.database == "cellchat"
    assert params.commot_dis_thr == 200.0
    
    # Test invalid min_cells (should be >= 0)
    with pytest.raises(ValueError):
        CellCommunicationParameters(min_cells=-1)
    
    # Test invalid distance threshold (should be > 0)
    with pytest.raises(ValueError):
        CellCommunicationParameters(commot_dis_thr=0.0)
    
    # Test invalid cutoff (should be between 0 and 1)
    with pytest.raises(ValueError):
        CellCommunicationParameters(spatialdm_cutoff=1.5)


def test_cell_communication_result_creation():
    """Test CellCommunicationResult creation"""
    result = CellCommunicationResult(
        data_id="test_data",
        method="commot",
        species="human",
        database="cellchat",
        n_lr_pairs=10,
        n_significant_pairs=5,
        top_lr_pairs=["TGFB1_TGFBR1", "PDGFA_PDGFRA", "VEGFA_VEGFR1"],
        statistics={"method": "commot", "distance_threshold": 200.0}
    )
    
    assert result.data_id == "test_data"
    assert result.method == "commot"
    assert result.species == "human"
    assert result.database == "cellchat"
    assert result.n_lr_pairs == 10
    assert result.n_significant_pairs == 5
    assert len(result.top_lr_pairs) == 3
    assert result.top_lr_pairs[0] == "TGFB1_TGFBR1"


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
    params = CellCommunicationParameters(method="cellphonedb")  # Not implemented yet
    mock_context = AsyncMock()
    
    # This should raise a NotImplementedError
    with pytest.raises(RuntimeError, match="not yet implemented"):
        await analyze_cell_communication("test_data", mock_data_store, params, mock_context)


if __name__ == "__main__":
    pytest.main([__file__])
