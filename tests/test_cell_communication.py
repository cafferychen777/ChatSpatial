"""
Comprehensive tests for cell_communication.py

This module tests cell-cell communication analysis functionality for spatial transcriptomics data.
"""

import pytest
import numpy as np
import pandas as pd
import anndata as ad
from unittest.mock import Mock, patch, AsyncMock, MagicMock
from scipy import sparse
import warnings

from chatspatial.tools.cell_communication import (
    analyze_cell_communication,
    LIANA_AVAILABLE,
    _analyze_communication_liana,
    _detect_species_from_genes,
    _get_liana_resource_name
)
from chatspatial.models.data import CellCommunicationParameters
from chatspatial.models.analysis import CellCommunicationResult


class TestCellCommunicationAnalysis:
    """Test cell communication analysis functions"""
    
    @pytest.fixture
    def mock_spatial_adata(self):
        """Create mock spatial AnnData object with cell types"""
        n_obs = 200
        n_vars = 1000
        
        # Create expression matrix
        X = sparse.random(n_obs, n_vars, density=0.3, random_state=42).tocsr()
        
        obs = pd.DataFrame({
            'sample': ['S1'] * n_obs,
            'cell_type': np.random.choice(['T_cell', 'B_cell', 'Macrophage', 'Fibroblast'], n_obs)
        }, index=[f"cell_{i}" for i in range(n_obs)])
        
        var = pd.DataFrame({
            'gene_name': [f"Gene_{i}" for i in range(n_vars)],
            'gene_symbol': [f"Gene_{i}" for i in range(n_vars)]
        }, index=[f"Gene_{i}" for i in range(n_vars)])
        
        adata = ad.AnnData(X=X, obs=obs, var=var)
        # Add spatial coordinates
        adata.obsm['spatial'] = np.random.rand(n_obs, 2) * 100
        
        return adata
    
    @pytest.fixture
    def mock_data_store(self, mock_spatial_adata):
        """Create mock data store"""
        return {
            'spatial_data': {'adata': mock_spatial_adata, 'name': 'spatial_data'}
        }
    
    @pytest.fixture
    def mock_context(self):
        """Create mock context"""
        context = AsyncMock()
        context.info = AsyncMock()
        context.warning = AsyncMock()
        context.error = AsyncMock()
        return context
    
    @pytest.mark.asyncio
    async def test_analyze_communication_dataset_not_found(self, mock_data_store, mock_context):
        """Test error when dataset not found"""
        params = CellCommunicationParameters(method="liana")
        
        with pytest.raises(ValueError, match="Dataset nonexistent not found"):
            await analyze_cell_communication("nonexistent", mock_data_store, params, mock_context)
    
    @pytest.mark.asyncio
    async def test_analyze_communication_no_spatial_coords(self, mock_data_store, mock_context):
        """Test error when spatial coordinates missing"""
        # Remove spatial coordinates
        del mock_data_store['spatial_data']['adata'].obsm['spatial']
        
        params = CellCommunicationParameters(method="liana")
        
        with pytest.raises(RuntimeError, match="No spatial coordinates found"):
            await analyze_cell_communication("spatial_data", mock_data_store, params, mock_context)
    
    @pytest.mark.asyncio
    async def test_analyze_communication_invalid_method(self, mock_data_store, mock_context):
        """Test error with invalid method"""
        # CellCommunicationParameters validates the method, so we need to bypass
        params = CellCommunicationParameters(method="liana")
        params.method = "invalid_method"
        
        with pytest.raises(RuntimeError, match="Unsupported method"):
            await analyze_cell_communication("spatial_data", mock_data_store, params, mock_context)
    
    @pytest.mark.asyncio
    async def test_analyze_communication_liana_not_installed(self, mock_data_store, mock_context):
        """Test error when LIANA not installed"""
        params = CellCommunicationParameters(method="liana")
        
        with patch('chatspatial.tools.cell_communication.LIANA_AVAILABLE', False):
            with pytest.raises(RuntimeError, match="LIANA"):
                await analyze_cell_communication("spatial_data", mock_data_store, params, mock_context)
    
    @pytest.mark.asyncio
    async def test_analyze_communication_liana_success(self, mock_data_store, mock_context):
        """Test successful LIANA analysis"""
        params = CellCommunicationParameters(
            method="liana",
            cell_type_key="cell_type",
            min_cells=5
        )
        
        # Mock LIANA analysis result
        mock_result_data = {
            "n_lr_pairs": 100,
            "n_significant_pairs": 25,
            "global_results_key": "liana_res",
            "top_lr_pairs": [
                "CXCL12_CXCR4",
                "CCL5_CCR5"
            ],
            "local_analysis_performed": False,
            "analysis_type": "basic"
        }
        
        with patch('chatspatial.tools.cell_communication.LIANA_AVAILABLE', True):
            with patch('chatspatial.tools.cell_communication._analyze_communication_liana', return_value=mock_result_data) as mock_analyze:
                result = await analyze_cell_communication("spatial_data", mock_data_store, params, mock_context)
                
                assert isinstance(result, CellCommunicationResult)
                assert result.method == "liana"
                assert result.n_lr_pairs == 100
                assert result.n_significant_pairs == 25
                assert len(result.top_lr_pairs) == 2
                assert result.database == "liana"
                
                # Check that analyze was called
                mock_analyze.assert_called_once()
    
    @pytest.mark.asyncio
    async def test_data_normalization(self, mock_data_store, mock_context):
        """Test that data is normalized if not already"""
        params = CellCommunicationParameters(method="liana")
        
        # Ensure log1p is not in uns
        if 'log1p' in mock_data_store['spatial_data']['adata'].uns:
            del mock_data_store['spatial_data']['adata'].uns['log1p']
        
        with patch('chatspatial.tools.cell_communication.LIANA_AVAILABLE', True):
            with patch('chatspatial.tools.cell_communication._analyze_communication_liana') as mock_analyze:
                mock_analyze.return_value = {
                    "n_lr_pairs": 0,
                    "n_significant_pairs": 0,
                    "analysis_type": "basic"
                }
                
                with patch('scanpy.pp.normalize_total') as mock_normalize:
                    with patch('scanpy.pp.log1p') as mock_log1p:
                        await analyze_cell_communication("spatial_data", mock_data_store, params, mock_context)
                        
                        # Check normalization was called
                        mock_normalize.assert_called_once()
                        mock_log1p.assert_called_once()


class TestHelperFunctions:
    """Test helper functions"""
    
    def test_detect_species_from_genes(self):
        """Test species detection from gene names"""
        # Create mock data with human genes
        adata = ad.AnnData(X=np.random.rand(10, 5))
        adata.var_names = ['GAPDH', 'ACTB', 'TP53', 'EGFR', 'BRAF']
        
        species = _detect_species_from_genes(adata)
        assert species == 'human'
        
        # Create mock data with mouse genes
        adata.var_names = ['Gapdh', 'Actb', 'Tp53', 'Egfr', 'Braf']
        
        species = _detect_species_from_genes(adata)
        assert species == 'mouse'
    
    def test_get_liana_resource_name(self):
        """Test LIANA resource name selection"""
        # Test human resources
        assert _get_liana_resource_name('human', 'consensus') == 'consensus'
        assert _get_liana_resource_name('human', 'cellphonedb') == 'cellphonedb'
        assert _get_liana_resource_name('human', 'invalid') == 'invalid'  # Returns as-is for human
        
        # Test mouse resources
        assert _get_liana_resource_name('mouse', 'consensus') == 'mouseconsensus'
        assert _get_liana_resource_name('mouse', 'cellphonedb') == 'cellphonedb'  # Available for mouse


class TestLianaAnalysis:
    """Test LIANA-specific analysis"""
    
    @pytest.mark.asyncio
    async def test_analyze_communication_liana_basic(self):
        """Test basic LIANA analysis workflow"""
        # Create mock data
        n_obs = 100
        n_vars = 500
        
        X = np.random.rand(n_obs, n_vars)
        
        obs = pd.DataFrame({
            'cell_type': np.random.choice(['T_cell', 'B_cell', 'Macrophage'], n_obs)
        })
        
        var = pd.DataFrame({
            'gene_name': [f"Gene_{i}" for i in range(n_vars)]
        }, index=[f"Gene_{i}" for i in range(n_vars)])
        
        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.obsm['spatial'] = np.random.rand(n_obs, 2)
        # Add spatial connectivity to bypass computation
        adata.obsp['spatial_connectivities'] = sparse.random(n_obs, n_obs, density=0.1)
        
        params = CellCommunicationParameters(
            method="liana",
            cell_type_key="cell_type",
            analysis_type="basic",
            perform_spatial_analysis=False  # Force cluster analysis
        )
        
        mock_context = AsyncMock()
        
        # Mock LIANA
        mock_li = MagicMock()
        
        # Mock method functions
        mock_mt = MagicMock()
        mock_rank_aggregate = MagicMock()
        mock_mt.rank_aggregate = mock_rank_aggregate
        mock_li.mt = mock_mt
        # Add ut module with spatial_neighbors method
        mock_ut = MagicMock()
        mock_spatial_neighbors = MagicMock()
        mock_ut.spatial_neighbors = mock_spatial_neighbors
        mock_li.ut = mock_ut
        
        # Setup resource
        mock_get_resource = MagicMock()
        mock_li.get_resource = mock_get_resource
        
        # Mock results
        adata.uns['liana_res'] = pd.DataFrame({
            'ligand_complex': ['CXCL12', 'CCL5'],
            'receptor_complex': ['CXCR4', 'CCR5'],
            'source': ['Fibroblast', 'T_cell'],
            'target': ['T_cell', 'Macrophage'],
            'magnitude_rank': [1, 2],
            'specificity_rank': [1, 2]
        })
        
        # Patch the liana import inside _analyze_communication_liana
        with patch.dict('sys.modules', {'liana': mock_li}):
            with patch('chatspatial.tools.cell_communication._detect_species_from_genes', return_value='human'):
                result = await _analyze_communication_liana(adata, params, mock_context)
            
            assert result['n_lr_pairs'] == 2
            assert result['analysis_type'] == 'cluster'
            assert result['liana_results_key'] == 'liana_res'


class TestIntegrationScenarios:
    """Test complete communication analysis workflows"""
    
    @pytest.mark.asyncio
    async def test_complete_liana_workflow(self):
        """Test complete LIANA workflow"""
        # Create realistic test data
        n_obs = 200
        n_vars = 1000
        
        # Create expression data with some structure
        X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars))
        
        # Create cell types with spatial organization
        cell_types = []
        for i in range(n_obs):
            if i < 50:
                cell_types.append('T_cell')
            elif i < 100:
                cell_types.append('B_cell')
            elif i < 150:
                cell_types.append('Macrophage')
            else:
                cell_types.append('Fibroblast')
        
        obs = pd.DataFrame({
            'cell_type': cell_types
        })
        
        var = pd.DataFrame({
            'gene_name': [f"Gene_{i}" for i in range(n_vars)]
        }, index=[f"Gene_{i}" for i in range(n_vars)])
        
        adata = ad.AnnData(X=X, obs=obs, var=var)
        
        # Add spatial coordinates with some clustering
        spatial = np.zeros((n_obs, 2))
        for i in range(n_obs):
            if i < 50:  # T cells clustered in one area
                spatial[i] = np.random.normal([20, 20], 5, 2)
            elif i < 100:  # B cells in another
                spatial[i] = np.random.normal([80, 20], 5, 2)
            elif i < 150:  # Macrophages
                spatial[i] = np.random.normal([20, 80], 5, 2)
            else:  # Fibroblasts spread out
                spatial[i] = np.random.uniform(0, 100, 2)
        
        adata.obsm['spatial'] = spatial
        
        data_store = {
            'test_data': {'adata': adata, 'name': 'test_data'}
        }
        
        context = AsyncMock()
        
        params = CellCommunicationParameters(
            method="liana",
            cell_type_key="cell_type",
            min_cells=10,
            threshold=0.05,
            analysis_type="basic"
        )
        
        # Mock LIANA
        with patch('chatspatial.tools.cell_communication.LIANA_AVAILABLE', True):
            with patch('chatspatial.tools.cell_communication._analyze_communication_liana') as mock_analyze:
                mock_analyze.return_value = {
                    "n_lr_pairs": 50,
                    "n_significant_pairs": 15,
                    "global_results_key": "liana_res",
                    "top_lr_pairs": [
                        "CXCL12_CXCR4"
                    ],
                    "local_analysis_performed": False,
                    "analysis_type": "basic",
                    "patterns_identified": True
                }
                
                result = await analyze_cell_communication("test_data", data_store, params, context)
                
                assert isinstance(result, CellCommunicationResult)
                assert result.n_lr_pairs == 50
                assert result.n_significant_pairs == 15
                assert result.patterns_identified == True
                assert len(result.top_lr_pairs) == 1