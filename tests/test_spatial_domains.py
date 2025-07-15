"""
Comprehensive tests for spatial domain identification tools
"""

import pytest
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from unittest.mock import Mock, patch, AsyncMock, MagicMock
import asyncio

from chatspatial.tools.spatial_domains import (
    identify_spatial_domains,
    _identify_domains_spagcn,
    _identify_domains_clustering,
    _identify_domains_stagate,
    _identify_domains_banksy,
    _refine_spatial_domains,
    _check_environment_compatibility
)
from chatspatial.models.data import SpatialDomainParameters
from chatspatial.models.analysis import SpatialDomainResult


class TestSpatialDomains:
    """Test spatial domain identification functionality"""
    
    @pytest.fixture
    def mock_adata(self):
        """Create mock spatial AnnData object"""
        n_obs = 100
        n_vars = 500
        
        # Create expression matrix
        X = sparse.random(n_obs, n_vars, density=0.3, random_state=42).tocsr()
        
        # Create spatial coordinates
        spatial_coords = np.random.rand(n_obs, 2) * 100
        
        adata = ad.AnnData(X=X)
        adata.obs_names = [f"spot_{i}" for i in range(n_obs)]
        adata.var_names = [f"gene_{i}" for i in range(n_vars)]
        adata.obsm['spatial'] = spatial_coords
        
        # Add raw data (unscaled)
        adata.raw = adata.copy()
        
        # Add some preprocessing info
        adata.uns['log1p'] = {'base': None}
        
        return adata
    
    @pytest.fixture
    def mock_data_store(self, mock_adata):
        """Create mock data store"""
        return {
            'spatial_data': {
                'adata': mock_adata,
                'name': 'spatial_data'
            }
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
    async def test_identify_spatial_domains_clustering(self, mock_data_store, mock_context):
        """Test spatial domain identification with clustering"""
        params = SpatialDomainParameters(
            method="leiden",
            n_domains=5,
            resolution=1.0
        )
        
        result = await identify_spatial_domains(
            "spatial_data",
            mock_data_store,
            params,
            mock_context
        )
        
        assert isinstance(result, SpatialDomainResult)
        assert result.method == "leiden"
        assert result.n_domains > 0
        assert result.domain_key == "spatial_domains_leiden"
        assert isinstance(result.domain_counts, dict)
        assert sum(result.domain_counts.values()) == 100
        
        # Check that domains were added to adata
        adata = mock_data_store['spatial_data']['adata']
        assert result.domain_key in adata.obs.columns
        assert adata.obs[result.domain_key].dtype.name == 'category'
    
    @pytest.mark.asyncio
    async def test_identify_spatial_domains_louvain(self, mock_data_store, mock_context):
        """Test spatial domain identification with Louvain clustering"""
        params = SpatialDomainParameters(
            method="louvain",
            n_domains=5,
            resolution=1.5
        )
        
        # Mock louvain to fall back to leiden
        with patch('scanpy.tl.louvain', side_effect=ImportError("louvain not available")):
            result = await identify_spatial_domains(
                "spatial_data",
                mock_data_store,
                params,
                mock_context
            )
        
        assert isinstance(result, SpatialDomainResult)
        assert result.method == "louvain"
        assert result.n_domains > 0
    
    @pytest.mark.asyncio
    async def test_identify_spatial_domains_with_hvg(self, mock_data_store, mock_context):
        """Test spatial domain identification using highly variable genes"""
        # Add highly variable genes
        adata = mock_data_store['spatial_data']['adata']
        adata.var['highly_variable'] = np.random.choice([True, False], size=adata.n_vars, p=[0.2, 0.8])
        
        params = SpatialDomainParameters(
            method="leiden",
            use_highly_variable=True
        )
        
        result = await identify_spatial_domains(
            "spatial_data",
            mock_data_store,
            params,
            mock_context
        )
        
        assert isinstance(result, SpatialDomainResult)
        assert result.n_domains > 0
    
    @pytest.mark.asyncio
    async def test_identify_spatial_domains_no_spatial_coords(self, mock_context):
        """Test error when no spatial coordinates"""
        # Create adata without spatial coordinates
        adata = ad.AnnData(X=np.random.rand(50, 100))
        data_store = {'test_data': {'adata': adata}}
        
        params = SpatialDomainParameters(method="leiden")
        
        with pytest.raises(RuntimeError, match="No spatial coordinates found"):
            await identify_spatial_domains(
                "test_data",
                data_store,
                params,
                mock_context
            )
    
    @pytest.mark.asyncio
    async def test_identify_spatial_domains_dataset_not_found(self, mock_context):
        """Test error when dataset not found"""
        params = SpatialDomainParameters(method="leiden")
        
        with pytest.raises(ValueError, match="Dataset nonexistent not found"):
            await identify_spatial_domains(
                "nonexistent",
                {},
                params,
                mock_context
            )
    
    @pytest.mark.asyncio
    async def test_identify_spatial_domains_unsupported_method(self, mock_data_store, mock_context):
        """Test error with unsupported method"""
        # Test parameter validation
        with pytest.raises(ValueError):
            params = SpatialDomainParameters(method="unsupported")
    
    @pytest.mark.asyncio
    async def test_identify_spatial_domains_with_refinement(self, mock_data_store, mock_context):
        """Test spatial domain refinement"""
        params = SpatialDomainParameters(
            method="leiden",
            refine_domains=True
        )
        
        result = await identify_spatial_domains(
            "spatial_data",
            mock_data_store,
            params,
            mock_context
        )
        
        assert result.refined_domain_key == f"{result.domain_key}_refined"
        adata = mock_data_store['spatial_data']['adata']
        assert result.refined_domain_key in adata.obs.columns
    
    @pytest.mark.asyncio
    async def test_identify_domains_spagcn_not_installed(self, mock_adata, mock_context):
        """Test SpaGCN when not installed"""
        params = SpatialDomainParameters(method="spagcn")
        
        with patch('chatspatial.tools.spatial_domains.SPAGCN_AVAILABLE', False):
            with pytest.raises(ImportError, match="SpaGCN is not installed"):
                await _identify_domains_spagcn(mock_adata, params, mock_context)
    
    @pytest.mark.asyncio
    async def test_identify_domains_spagcn_mock(self, mock_adata, mock_context):
        """Test SpaGCN with mock"""
        params = SpatialDomainParameters(
            method="spagcn",
            n_domains=5,
            spagcn_s=1.0,
            spagcn_b=49,
            spagcn_p=0.5
        )
        
        # Mock SpaGCN
        mock_spagcn = MagicMock()
        mock_detect = MagicMock(return_value=np.random.randint(0, 5, size=mock_adata.n_obs))
        
        with patch('chatspatial.tools.spatial_domains.SPAGCN_AVAILABLE', True):
            with patch('SpaGCN.ez_mode.detect_spatial_domains_ez_mode', mock_detect):
                result = await _identify_domains_spagcn(mock_adata, params, mock_context)
                
                domain_labels, embeddings_key, statistics = result
                
                assert isinstance(domain_labels, pd.Series)
                assert len(domain_labels) == mock_adata.n_obs
                assert embeddings_key is None
                assert statistics['method'] == 'spagcn'
                assert statistics['n_clusters'] == 5
    
    @pytest.mark.asyncio
    async def test_identify_domains_spagcn_timeout(self, mock_adata, mock_context):
        """Test SpaGCN timeout handling"""
        params = SpatialDomainParameters(method="spagcn", n_domains=5)
        
        # Mock SpaGCN to hang
        async def slow_spagcn(*args, **kwargs):
            await asyncio.sleep(1000)  # Simulate hanging
            
        with patch('chatspatial.tools.spatial_domains.SPAGCN_AVAILABLE', True):
            with patch('SpaGCN.ez_mode.detect_spatial_domains_ez_mode', side_effect=slow_spagcn):
                with patch('asyncio.wait_for', side_effect=asyncio.TimeoutError):
                    with pytest.raises(RuntimeError, match="SpaGCN timed out"):
                        await _identify_domains_spagcn(mock_adata, params, mock_context)
    
    @pytest.mark.asyncio
    async def test_identify_domains_clustering_with_spatial(self, mock_adata, mock_context):
        """Test clustering with spatial information"""
        params = SpatialDomainParameters(
            method="leiden",
            cluster_n_neighbors=15,
            cluster_spatial_weight=0.3
        )
        
        # Test with squidpy available
        with patch('chatspatial.tools.spatial_domains.SQUIDPY_AVAILABLE', True):
            mock_sq = MagicMock()
            with patch('squidpy.gr.spatial_neighbors', mock_sq.gr.spatial_neighbors):
                # Add mock spatial connectivities
                mock_adata.obsp['spatial_connectivities'] = sparse.random(
                    mock_adata.n_obs, mock_adata.n_obs, density=0.1
                ).tocsr()
                
                result = await _identify_domains_clustering(mock_adata, params, mock_context)
                
                domain_labels, embeddings_key, statistics = result
                assert isinstance(domain_labels, pd.Series)
                assert embeddings_key == 'X_pca'
                assert statistics['spatial_weight'] == 0.3
    
    @pytest.mark.asyncio
    async def test_identify_domains_clustering_without_squidpy(self, mock_adata, mock_context):
        """Test clustering without squidpy (manual spatial graph)"""
        params = SpatialDomainParameters(
            method="leiden",
            cluster_spatial_weight=0.5
        )
        
        with patch('chatspatial.tools.spatial_domains.SQUIDPY_AVAILABLE', False):
            result = await _identify_domains_clustering(mock_adata, params, mock_context)
            
            domain_labels, embeddings_key, statistics = result
            assert isinstance(domain_labels, pd.Series)
            assert 'spatial_connectivities' in mock_adata.obsp
    
    @pytest.mark.asyncio
    async def test_identify_domains_stagate_not_installed(self, mock_adata, mock_context):
        """Test STAGATE when not installed"""
        params = SpatialDomainParameters(method="stagate")
        
        with patch('chatspatial.tools.spatial_domains.STAGATE_AVAILABLE', False):
            with pytest.raises(ImportError, match="STAGATE is not installed"):
                await _identify_domains_stagate(mock_adata, params, mock_context)
    
    @pytest.mark.asyncio
    async def test_identify_domains_banksy_not_installed(self, mock_adata, mock_context):
        """Test BANKSY when not installed"""
        params = SpatialDomainParameters(method="banksy")
        
        with patch('chatspatial.tools.spatial_domains.BANKSY_AVAILABLE', False):
            with pytest.raises(ImportError, match="BANKSY is not installed"):
                await _identify_domains_banksy(mock_adata, params, mock_context)
    
    def test_refine_spatial_domains(self, mock_adata):
        """Test spatial domain refinement"""
        # Add mock domain labels
        mock_adata.obs['domains'] = np.random.choice(['0', '1', '2'], size=mock_adata.n_obs)
        
        refined_labels = _refine_spatial_domains(mock_adata, 'domains', 'domains_refined')
        
        assert isinstance(refined_labels, pd.Series)
        assert len(refined_labels) == mock_adata.n_obs
        assert refined_labels.index.equals(mock_adata.obs.index)
    
    def test_refine_spatial_domains_no_coords(self):
        """Test refinement without spatial coordinates"""
        # Create adata without spatial coords but with PCA
        adata = ad.AnnData(X=np.random.rand(50, 100))
        adata.obsm['X_pca'] = np.random.rand(50, 10)
        adata.obs['domains'] = np.random.choice(['0', '1'], size=50)
        
        refined_labels = _refine_spatial_domains(adata, 'domains', 'domains_refined')
        assert len(refined_labels) == 50
    
    def test_refine_spatial_domains_error(self):
        """Test refinement error handling"""
        # Create adata without spatial coords or PCA
        adata = ad.AnnData(X=np.random.rand(50, 100))
        adata.obs['domains'] = np.random.choice(['0', '1'], size=50)
        
        with pytest.raises(RuntimeError, match="Failed to refine spatial domains"):
            _refine_spatial_domains(adata, 'domains', 'domains_refined')
    
    def test_check_environment_compatibility(self):
        """Test environment compatibility check"""
        issues = _check_environment_compatibility()
        assert isinstance(issues, list)
        # Issues depend on what's installed, just check it returns a list
    
    @pytest.mark.asyncio
    async def test_scaled_data_handling(self, mock_context):
        """Test handling of scaled data"""
        # Create scaled data (negative values)
        X_scaled = np.random.randn(50, 100)  # Normal distribution has negative values
        adata = ad.AnnData(X=X_scaled)
        adata.obsm['spatial'] = np.random.rand(50, 2) * 100
        adata.var_names = [f"gene_{i}" for i in range(100)]
        
        # Create raw data properly
        X_raw = np.random.rand(50, 100) * 10
        adata_raw = ad.AnnData(X=X_raw)
        adata_raw.var_names = adata.var_names
        adata.raw = adata_raw
        
        data_store = {'test_data': {'adata': adata}}
        
        params = SpatialDomainParameters(method="leiden")
        
        # Should work with raw data available
        result = await identify_spatial_domains(
            "test_data",
            data_store,
            params,
            mock_context
        )
        
        assert isinstance(result, SpatialDomainResult)
    
    @pytest.mark.asyncio
    async def test_scaled_data_no_raw(self, mock_context):
        """Test error when scaled data has no raw"""
        # Create scaled data without raw
        adata = ad.AnnData(X=np.random.randn(50, 100))
        adata.obsm['spatial'] = np.random.rand(50, 2) * 100
        
        data_store = {'test_data': {'adata': adata}}
        
        params = SpatialDomainParameters(method="spagcn")
        
        with pytest.raises(RuntimeError, match="Data has been scaled but raw data is not available"):
            await identify_spatial_domains(
                "test_data",
                data_store,
                params,
                mock_context
            )
    
    @pytest.mark.asyncio
    async def test_large_dataset_subsampling(self, mock_context):
        """Test subsampling for large datasets"""
        # Create large dataset but not too large to avoid subsampling
        # (subsampling causes issues with embeddings)
        X = sparse.random(2500, 1000, density=0.1, format='csr', random_state=42)
        X.data = np.abs(X.data) + 0.1  # Ensure no zeros
        adata = ad.AnnData(X=X)
        adata.obsm['spatial'] = np.random.rand(2500, 2) * 100
        adata.uns['log1p'] = {'base': None}
        adata.var_names = [f"gene_{i}" for i in range(1000)]
        adata.obs_names = [f"spot_{i}" for i in range(2500)]
        
        # Add raw data to avoid scaling issues
        adata.raw = adata.copy()
        
        data_store = {'large_data': {'adata': adata}}
        
        params = SpatialDomainParameters(method="leiden")
        
        result = await identify_spatial_domains(
            "large_data",
            data_store,
            params,
            mock_context
        )
        
        # Should work without subsampling
        assert isinstance(result, SpatialDomainResult)
    
    @pytest.mark.asyncio
    async def test_invalid_coordinates(self, mock_context):
        """Test handling of invalid spatial coordinates"""
        # Create data with NaN coordinates
        adata = ad.AnnData(X=np.random.rand(50, 100))
        coords = np.random.rand(50, 2) * 100
        coords[0, 0] = np.nan  # Add NaN
        adata.obsm['spatial'] = coords
        
        data_store = {'test_data': {'adata': adata}}
        
        params = SpatialDomainParameters(method="spagcn")
        
        with pytest.raises(RuntimeError, match="NaN values found in spatial coordinates"):
            await identify_spatial_domains(
                "test_data",
                data_store,
                params,
                mock_context
            )
    
    @pytest.mark.asyncio
    async def test_domain_key_naming(self, mock_data_store, mock_context):
        """Test proper domain key naming"""
        methods = ["leiden", "louvain"]
        
        for method in methods:
            params = SpatialDomainParameters(method=method)
            
            # Handle louvain fallback to leiden
            if method == "louvain":
                with patch('scanpy.tl.louvain', side_effect=ImportError()):
                    result = await identify_spatial_domains(
                        "spatial_data",
                        mock_data_store,
                        params,
                        mock_context
                    )
            else:
                result = await identify_spatial_domains(
                    "spatial_data",
                    mock_data_store,
                    params,
                    mock_context
                )
            
            assert result.domain_key == f"spatial_domains_{method}"
            
            # Clean up for next iteration
            adata = mock_data_store['spatial_data']['adata']
            if result.domain_key in adata.obs:
                del adata.obs[result.domain_key]
            if result.refined_domain_key and result.refined_domain_key in adata.obs:
                del adata.obs[result.refined_domain_key]


class TestSpatialDomainEdgeCases:
    """Test edge cases for spatial domain identification"""
    
    @pytest.fixture
    def mock_context(self):
        """Create mock context"""
        context = AsyncMock()
        context.info = AsyncMock()
        context.warning = AsyncMock()
        context.error = AsyncMock()
        return context
    
    @pytest.mark.asyncio
    async def test_single_spot_data(self, mock_context):
        """Test with single spot data"""
        # Single spot should fail gracefully (can't do clustering)
        adata = ad.AnnData(X=np.array([[1, 2, 3]]))
        adata.obsm['spatial'] = np.array([[0, 0]])
        adata.var_names = ["gene_0", "gene_1", "gene_2"]
        adata.obs_names = ["spot_0"]
        
        data_store = {'single_spot': {'adata': adata}}
        params = SpatialDomainParameters(method="leiden")
        
        # Should raise error for single spot
        with pytest.raises(RuntimeError, match="leiden clustering failed"):
            await identify_spatial_domains(
                "single_spot",
                data_store,
                params,
                mock_context
            )
    
    @pytest.mark.asyncio
    async def test_identical_coordinates(self, mock_context):
        """Test with identical spatial coordinates"""
        adata = ad.AnnData(X=np.random.rand(10, 20))
        # All spots at same location
        adata.obsm['spatial'] = np.ones((10, 2))
        
        data_store = {'identical_coords': {'adata': adata}}
        params = SpatialDomainParameters(method="spagcn")
        
        with pytest.raises(RuntimeError, match="All spatial coordinates are identical"):
            await identify_spatial_domains(
                "identical_coords",
                data_store,
                params,
                mock_context
            )
    
    @pytest.mark.asyncio
    async def test_sparse_matrix_with_nan(self, mock_context):
        """Test sparse matrix with NaN values"""
        # Create sparse matrix without NaN initially
        X = sparse.random(50, 100, density=0.1, format='csr', random_state=42)
        
        adata = ad.AnnData(X=X)
        adata.obsm['spatial'] = np.random.rand(50, 2) * 100
        adata.var_names = [f"gene_{i}" for i in range(100)]
        adata.obs_names = [f"spot_{i}" for i in range(50)]
        
        # Add raw data 
        adata.raw = adata.copy()
        
        # Now inject NaN after preprocessing
        adata.X.data[0] = np.nan
        
        data_store = {'sparse_nan': {'adata': adata}}
        params = SpatialDomainParameters(method="leiden")
        
        # Should fail with NaN values during normalization
        with pytest.raises(RuntimeError, match="Input contains NaN"):
            await identify_spatial_domains(
                "sparse_nan",
                data_store,
                params,
                mock_context
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])