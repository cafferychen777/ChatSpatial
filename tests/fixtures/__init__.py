"""
Test fixtures for ChatSpatial MCP Server
Provides standardized test data and utilities for all test layers
"""

import pytest
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from typing import Dict, Any
from anndata import AnnData


@pytest.fixture(scope="session")
def standard_adata():
    """Standard test dataset with full preprocessing pipeline"""
    np.random.seed(42)
    n_cells, n_genes = 500, 1000
    
    # Generate synthetic expression data
    X = np.random.poisson(3, (n_cells, n_genes)).astype(np.float32)
    
    # Create AnnData object
    adata = AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell_{i:03d}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"Gene_{i:03d}" for i in range(n_genes)])
    )
    
    # Add spatial coordinates (4 distinct clusters)
    centers = np.array([[0, 0], [15, 0], [0, 15], [15, 15]])
    coords = []
    clusters = []
    
    for i in range(n_cells):
        cluster_idx = i % 4
        center = centers[cluster_idx]
        coord = center + np.random.normal(0, 3, 2)
        coords.append(coord)
        clusters.append(f"Cluster_{cluster_idx}")
    
    adata.obsm['spatial'] = np.array(coords)
    adata.obs['leiden'] = pd.Categorical(clusters)
    adata.obs['batch'] = pd.Categorical(['Batch1', 'Batch2'] * (n_cells // 2))
    adata.obs['cell_type'] = pd.Categorical(['TypeA', 'TypeB', 'TypeC'] * (n_cells // 3 + 1))[:n_cells]
    
    # Full preprocessing pipeline
    sc.pp.highly_variable_genes(adata, n_top_genes=100)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    
    # Add neighborhood enrichment results
    n_clusters = len(adata.obs['leiden'].cat.categories)
    enrichment_matrix = np.random.normal(0, 2, (n_clusters, n_clusters))
    enrichment_matrix = (enrichment_matrix + enrichment_matrix.T) / 2
    np.fill_diagonal(enrichment_matrix, np.abs(np.diagonal(enrichment_matrix)) + 2)
    adata.uns['leiden_nhood_enrichment'] = {'zscore': enrichment_matrix}
    
    return adata


@pytest.fixture(scope="session")
def minimal_adata():
    """Minimal test dataset with only basic data"""
    np.random.seed(123)
    n_cells, n_genes = 100, 200
    
    X = np.random.poisson(2, (n_cells, n_genes)).astype(np.float32)
    
    adata = AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell_{i:03d}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"Gene_{i:03d}" for i in range(n_genes)])
    )
    
    # Only spatial coordinates, no preprocessing
    x_coords = np.random.normal(0, 5, n_cells)
    y_coords = np.random.normal(0, 5, n_cells)
    adata.obsm['spatial'] = np.column_stack([x_coords, y_coords])
    
    return adata


@pytest.fixture(scope="session")
def multi_batch_adata():
    """Multi-batch dataset for integration testing"""
    np.random.seed(456)
    n_cells, n_genes = 300, 500
    
    X = np.random.poisson(4, (n_cells, n_genes)).astype(np.float32)
    
    adata = AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell_{i:03d}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"Gene_{i:03d}" for i in range(n_genes)])
    )
    
    # Spatial coordinates with batch effects
    batch_labels = []
    coords = []
    
    for i in range(n_cells):
        batch = f"Batch_{i % 3 + 1}"
        batch_labels.append(batch)
        
        # Add batch-specific spatial offset
        offset = np.array([0, 0]) if batch == "Batch_1" else \
                 np.array([20, 0]) if batch == "Batch_2" else \
                 np.array([10, 20])
        
        coord = offset + np.random.normal(0, 4, 2)
        coords.append(coord)
    
    adata.obsm['spatial'] = np.array(coords)
    adata.obs['batch'] = pd.Categorical(batch_labels)
    adata.obs['cell_type'] = pd.Categorical(['TypeA', 'TypeB', 'TypeC', 'TypeD'] * (n_cells // 4 + 1))[:n_cells]
    
    # Basic preprocessing
    sc.pp.highly_variable_genes(adata, n_top_genes=50)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    return adata


@pytest.fixture(scope="session")
def rna_velocity_adata():
    """Dataset with RNA velocity layers"""
    np.random.seed(789)
    n_cells, n_genes = 200, 300
    
    # Spliced and unspliced counts
    spliced = np.random.poisson(5, (n_cells, n_genes)).astype(np.float32)
    unspliced = np.random.poisson(2, (n_cells, n_genes)).astype(np.float32)
    
    adata = AnnData(
        X=spliced,
        obs=pd.DataFrame(index=[f"Cell_{i:03d}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"Gene_{i:03d}" for i in range(n_genes)])
    )
    
    # Add velocity layers
    adata.layers['spliced'] = spliced
    adata.layers['unspliced'] = unspliced
    
    # Spatial coordinates
    x_coords = np.random.normal(0, 8, n_cells)
    y_coords = np.random.normal(0, 8, n_cells)
    adata.obsm['spatial'] = np.column_stack([x_coords, y_coords])
    
    # Basic processing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=30)
    
    return adata


@pytest.fixture(scope="session")
def corrupted_adata():
    """Dataset with various data quality issues for robustness testing"""
    np.random.seed(999)
    n_cells, n_genes = 150, 250
    
    X = np.random.poisson(3, (n_cells, n_genes)).astype(np.float32)
    
    # Introduce NaN and inf values
    X[10:15, 50:55] = np.nan
    X[20:25, 100:105] = np.inf
    X[30:35, 150:155] = -np.inf
    
    adata = AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell_{i:03d}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"Gene_{i:03d}" for i in range(n_genes)])
    )
    
    # Spatial coordinates with some missing values
    coords = np.random.normal(0, 6, (n_cells, 2))
    coords[50:55] = np.nan  # Missing spatial coordinates for some cells
    adata.obsm['spatial'] = coords
    
    # Inconsistent metadata
    adata.obs['leiden'] = pd.Categorical(['A', 'B', 'C'] * (n_cells // 3) + ['A'] * (n_cells % 3))
    adata.obs['batch'] = pd.Categorical(['Batch1'] * 100 + [''] * 50)  # Empty batch labels
    
    return adata


@pytest.fixture
def mock_data_store(standard_adata, minimal_adata, multi_batch_adata):
    """Mock data store for testing"""
    return {
        "standard": {"adata": standard_adata},
        "minimal": {"adata": minimal_adata},
        "multi_batch": {"adata": multi_batch_adata}
    }


@pytest.fixture
def test_image_path():
    """Path for saving test images"""
    test_dir = Path(__file__).parent
    test_image_dir = test_dir / "output_images"
    test_image_dir.mkdir(exist_ok=True)
    return test_image_dir


class TestUtils:
    """Utility functions for testing"""
    
    @staticmethod
    def assert_valid_image(image_obj):
        """Assert that an object is a valid Image"""
        from mcp.types import Image
        assert isinstance(image_obj, Image)
        assert image_obj.mimeType in ["image/png", "image/jpeg"]
        assert len(image_obj.data) > 0
    
    @staticmethod
    def assert_valid_error_response(error_obj):
        """Assert that an object is a valid error response"""
        assert hasattr(error_obj, 'isError')
        assert error_obj.isError is True
        assert hasattr(error_obj, 'content')
    
    @staticmethod
    def count_cells_in_cluster(adata, cluster_key, cluster_value):
        """Count cells in a specific cluster"""
        return sum(adata.obs[cluster_key] == cluster_value)
    
    @staticmethod
    def get_feature_values(adata, feature):
        """Get feature values from adata"""
        if feature in adata.var_names:
            return adata[:, feature].X.toarray().flatten()
        elif feature in adata.obs.columns:
            return adata.obs[feature].values
        else:
            return None


@pytest.fixture
def test_utils():
    """Test utilities fixture"""
    return TestUtils()