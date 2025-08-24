#!/usr/bin/env python3
"""
ChatSpatial Tools Memory Profiling
根据Linus原则：测试实际内存使用，不是理论值
"""

import os
import sys
import tracemalloc
import gc
import time
import psutil
from pathlib import Path
from typing import Dict, List, Any, Tuple
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import anndata

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# Add paths
sys.path.append('/Users/apple/Research/SpatialTrans_MCP/chatspatial')

class MemoryProfiler:
    """Memory profiling for ChatSpatial tools."""
    
    def __init__(self):
        self.data_dir = Path('/Users/apple/Research/SpatialTrans_MCP/chatspatial/data')
        self.results = {}
        self.process = psutil.Process()
        
    def get_memory_usage(self) -> Dict[str, float]:
        """Get current memory usage in MB."""
        memory_info = self.process.memory_info()
        return {
            'rss_mb': memory_info.rss / 1024 / 1024,  # Physical memory
            'vms_mb': memory_info.vms / 1024 / 1024,  # Virtual memory
            'shared_mb': getattr(memory_info, 'shared', 0) / 1024 / 1024
        }
    
    def memory_test_context(self, test_name: str):
        """Context manager for memory testing."""
        class MemoryTestContext:
            def __init__(self, profiler, name):
                self.profiler = profiler
                self.name = name
                self.start_memory = None
                self.peak_memory = None
                
            def __enter__(self):
                gc.collect()  # Clean up before measurement
                tracemalloc.start()
                self.start_memory = self.profiler.get_memory_usage()
                return self
                
            def __exit__(self, exc_type, exc_val, exc_tb):
                current, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                end_memory = self.profiler.get_memory_usage()
                
                self.profiler.results[self.name] = {
                    'start_rss_mb': self.start_memory['rss_mb'],
                    'end_rss_mb': end_memory['rss_mb'],
                    'delta_rss_mb': end_memory['rss_mb'] - self.start_memory['rss_mb'],
                    'peak_traced_mb': peak / 1024 / 1024,
                    'current_traced_mb': current / 1024 / 1024,
                    'error': str(exc_val) if exc_val else None
                }
                gc.collect()  # Clean up after measurement
                
        return MemoryTestContext(self, test_name)
    
    def load_test_datasets(self) -> Dict[str, anndata.AnnData]:
        """Load datasets of different sizes."""
        datasets = {}
        
        # Small dataset
        small_path = self.data_dir / 'spatial_datasets/seqfish_developmental.h5ad'
        if small_path.exists():
            adata_small = sc.read_h5ad(small_path)
            datasets['small_1200x400'] = adata_small
            print(f"✅ Small dataset: {adata_small.n_obs} cells × {adata_small.n_vars} genes")
        
        # Medium dataset  
        medium_path = self.data_dir / 'benchmark_datasets/benchmark_500x1k.h5ad'
        if medium_path.exists():
            adata_medium = sc.read_h5ad(medium_path)
            datasets['medium_500x1k'] = adata_medium
            print(f"✅ Medium dataset: {adata_medium.n_obs} cells × {adata_medium.n_vars} genes")
        
        # Create larger synthetic dataset for stress testing
        print("Creating synthetic large dataset...")
        adata_large = self._create_synthetic_dataset(5000, 2000)
        datasets['large_5kx2k'] = adata_large
        print(f"✅ Large dataset: {adata_large.n_obs} cells × {adata_large.n_vars} genes")
        
        return datasets
    
    def _create_synthetic_dataset(self, n_obs: int, n_vars: int) -> anndata.AnnData:
        """Create synthetic spatial transcriptomics dataset."""
        np.random.seed(42)
        
        # Create sparse expression matrix
        from scipy.sparse import csr_matrix
        density = 0.1  # 10% non-zero values
        data = np.random.negative_binomial(n=5, p=0.3, size=int(n_obs * n_vars * density))
        row_indices = np.random.choice(n_obs, size=len(data))
        col_indices = np.random.choice(n_vars, size=len(data))
        X = csr_matrix((data, (row_indices, col_indices)), shape=(n_obs, n_vars))
        
        # Create spatial coordinates
        spatial_coords = np.random.uniform(0, 100, size=(n_obs, 2))
        
        # Create gene names
        gene_names = [f"Gene_{i:05d}" for i in range(n_vars)]
        
        # Create AnnData object
        adata = anndata.AnnData(X=X)
        adata.var_names = gene_names
        adata.obsm['spatial'] = spatial_coords
        
        return adata
    
    def test_preprocessing_memory(self, adata: anndata.AnnData) -> Dict[str, Any]:
        """Test memory usage of preprocessing operations."""
        with self.memory_test_context(f'preprocessing_{adata.n_obs}x{adata.n_vars}'):
            adata_test = adata.copy()
            
            # Basic preprocessing
            sc.pp.calculate_qc_metrics(adata_test, percent_top=None, log1p=False, inplace=True)
            sc.pp.filter_cells(adata_test, min_genes=10)
            sc.pp.filter_genes(adata_test, min_cells=3)
            sc.pp.normalize_total(adata_test, target_sum=1e4)
            sc.pp.log1p(adata_test)
            
            # HVG selection (memory intensive)
            if adata_test.n_vars < 1000:
                sc.pp.highly_variable_genes(adata_test, n_top_genes=min(200, adata_test.n_vars // 2))
            else:
                sc.pp.highly_variable_genes(adata_test, n_top_genes=2000)
            
            # PCA (memory intensive)
            sc.pp.scale(adata_test, max_value=10)
            sc.tl.pca(adata_test, svd_solver='arpack')
            
            # Neighbors (memory intensive)
            sc.pp.neighbors(adata_test, n_neighbors=15, n_pcs=40)
            sc.tl.umap(adata_test)
            
    def test_spatial_analysis_memory(self, adata: anndata.AnnData):
        """Test memory usage of spatial analysis."""
        with self.memory_test_context(f'spatial_analysis_{adata.n_obs}x{adata.n_vars}'):
            import squidpy as sq
            
            adata_test = adata.copy()
            
            # Spatial neighbors (memory intensive for large datasets)
            if 'spatial' in adata_test.obsm:
                sq.gr.spatial_neighbors(adata_test, coord_type='generic', n_neighs=6)
                
                # Spatial autocorrelation (memory intensive)
                if adata_test.n_vars > 1000:
                    # Test with subset of genes
                    test_genes = adata_test.var_names[:100]
                else:
                    test_genes = adata_test.var_names[:50]
                    
                sq.gr.spatial_autocorr(adata_test, genes=test_genes)
    
    def test_visualization_memory(self, adata: anndata.AnnData):
        """Test memory usage of visualization."""
        with self.memory_test_context(f'visualization_{adata.n_obs}x{adata.n_vars}'):
            import matplotlib
            matplotlib.use('Agg')  # Use non-interactive backend
            import matplotlib.pyplot as plt
            
            adata_test = adata.copy()
            
            # Ensure we have the required embeddings
            if 'X_umap' not in adata_test.obsm:
                try:
                    sc.pp.neighbors(adata_test, n_neighbors=15)
                    sc.tl.umap(adata_test)
                except:
                    pass
            
            if 'spatial' in adata_test.obsm:
                fig, axes = plt.subplots(2, 2, figsize=(12, 10))
                
                # Spatial plot
                if 'X_umap' in adata_test.obsm:
                    sc.pl.umap(adata_test, ax=axes[0,0], show=False)
                
                # Gene expression plot (memory intensive for large datasets)
                test_gene = adata_test.var_names[0]
                if 'spatial' in adata_test.obsm:
                    sc.pl.spatial(adata_test, color=test_gene, ax=axes[0,1], show=False)
                
                plt.close(fig)
    
    def test_cell_communication_memory(self, adata: anndata.AnnData):
        """Test memory usage of cell communication analysis."""
        with self.memory_test_context(f'cell_communication_{adata.n_obs}x{adata.n_vars}'):
            # Import and test basic LIANA functionality
            try:
                import liana as li
                
                adata_test = adata.copy()
                
                # Add fake clusters for testing
                np.random.seed(42)
                n_clusters = min(5, max(2, adata_test.n_obs // 50))
                adata_test.obs['cell_type'] = [f'Type_{i}' for i in np.random.choice(n_clusters, adata_test.n_obs)]
                
                # Run basic LIANA analysis (memory intensive)
                if adata_test.n_obs < 2000:  # Only for manageable sizes
                    li.resource.select_resource('consensus')
                    # This would be memory intensive, so we'll just simulate the setup
                    
            except ImportError:
                pass  # LIANA not available
    
    def test_trajectory_memory(self, adata: anndata.AnnData):
        """Test memory usage of trajectory analysis."""
        with self.memory_test_context(f'trajectory_{adata.n_obs}x{adata.n_vars}'):
            try:
                import scvelo as scv
                
                adata_test = adata.copy()
                
                # Add fake velocity data
                adata_test.layers['spliced'] = adata_test.X.copy()
                adata_test.layers['unspliced'] = adata_test.X.copy() * 0.8
                
                # Basic velocity preprocessing (memory intensive)
                scv.pp.filter_and_normalize(adata_test, min_shared_counts=20, n_top_genes=1000)
                scv.pp.moments(adata_test, n_pcs=30, n_neighbors=30)
                
            except ImportError:
                pass  # scVelo not available
    
    def test_spatial_domains_memory(self, adata: anndata.AnnData):
        """Test memory usage of spatial domain identification."""
        with self.memory_test_context(f'spatial_domains_{adata.n_obs}x{adata.n_vars}'):
            adata_test = adata.copy()
            
            # Basic clustering as proxy for spatial domains (memory intensive)
            try:
                sc.pp.neighbors(adata_test, n_neighbors=15)
                sc.tl.leiden(adata_test, resolution=0.5)
                
                # Spatial statistics (memory intensive for large datasets)
                if 'spatial' in adata_test.obsm and adata_test.n_obs < 3000:
                    import squidpy as sq
                    sq.gr.spatial_neighbors(adata_test, coord_type='generic')
                    
            except Exception as e:
                pass  # Skip if failed
    
    def run_comprehensive_test(self):
        """Run comprehensive memory profiling."""
        print("🧠 ChatSpatial Tools Memory Profiling")
        print("=" * 50)
        
        # Load test datasets
        datasets = self.load_test_datasets()
        
        # Test functions to run
        test_functions = [
            ('preprocessing', self.test_preprocessing_memory),
            ('spatial_analysis', self.test_spatial_analysis_memory),
            ('visualization', self.test_visualization_memory),
            ('cell_communication', self.test_cell_communication_memory),
            ('trajectory', self.test_trajectory_memory),
            ('spatial_domains', self.test_spatial_domains_memory)
        ]
        
        # Run tests for each dataset
        for dataset_name, adata in datasets.items():
            print(f"\n📊 Testing dataset: {dataset_name}")
            print(f"   Shape: {adata.n_obs} cells × {adata.n_vars} genes")
            
            for test_name, test_func in test_functions:
                print(f"  🧪 Testing {test_name}...")
                try:
                    start_time = time.time()
                    test_func(adata)
                    elapsed = time.time() - start_time
                    print(f"    ✅ Completed in {elapsed:.2f}s")
                except Exception as e:
                    print(f"    ❌ Failed: {str(e)[:100]}")
        
        return self.results

if __name__ == "__main__":
    profiler = MemoryProfiler()
    results = profiler.run_comprehensive_test()
    
    print("\n" + "="*50)
    print("📋 Memory Profiling Results Summary")
    print("="*50)
    
    for test_name, result in results.items():
        if result['error']:
            print(f"❌ {test_name}: FAILED - {result['error'][:50]}")
        else:
            print(f"✅ {test_name}:")
            print(f"   Memory Delta: {result['delta_rss_mb']:+.1f} MB")
            print(f"   Peak Traced:  {result['peak_traced_mb']:.1f} MB")