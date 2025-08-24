#!/usr/bin/env python3
"""
Quick ChatSpatial Tools Memory Analysis
Linus风格：快速测试，关注实际问题
"""

import os
import sys
import tracemalloc
import gc
import time
import psutil
from pathlib import Path
from typing import Dict, List, Any
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

sys.path.append('/Users/apple/Research/SpatialTrans_MCP/chatspatial')

def get_memory_mb():
    """Get current memory usage in MB."""
    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024

def memory_test(test_name, func, *args, **kwargs):
    """Simple memory test wrapper."""
    print(f"  🧪 {test_name}...")
    gc.collect()
    
    start_mem = get_memory_mb()
    tracemalloc.start()
    start_time = time.time()
    
    try:
        result = func(*args, **kwargs)
        success = True
        error = None
    except Exception as e:
        result = None
        success = False
        error = str(e)[:100]
    
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    end_mem = get_memory_mb()
    elapsed = time.time() - start_time
    
    delta_mem = end_mem - start_mem
    peak_mem = peak / 1024 / 1024
    
    if success:
        print(f"    ✅ Success: {delta_mem:+.1f}MB delta, {peak_mem:.1f}MB peak, {elapsed:.1f}s")
    else:
        print(f"    ❌ Failed: {error}")
    
    gc.collect()
    return {
        'success': success,
        'delta_mb': delta_mem,
        'peak_mb': peak_mem,
        'time_s': elapsed,
        'error': error
    }

def test_preprocessing_memory(adata):
    """Test preprocessing memory footprint."""
    adata_test = adata.copy()
    
    # Basic preprocessing pipeline
    sc.pp.calculate_qc_metrics(adata_test, inplace=True)
    sc.pp.filter_cells(adata_test, min_genes=10)
    sc.pp.filter_genes(adata_test, min_cells=3)
    sc.pp.normalize_total(adata_test, target_sum=1e4)
    sc.pp.log1p(adata_test)
    
    # HVG selection
    if adata_test.n_vars < 1000:
        sc.pp.highly_variable_genes(adata_test, n_top_genes=min(200, adata_test.n_vars // 2))
    else:
        sc.pp.highly_variable_genes(adata_test, n_top_genes=min(2000, adata_test.n_vars))
    
    # PCA - memory intensive step
    sc.pp.scale(adata_test, max_value=10)
    n_comps = min(50, adata_test.n_vars - 1, adata_test.n_obs - 1)
    sc.tl.pca(adata_test, n_comps=n_comps, svd_solver='arpack')
    
    return f"PCA with {n_comps} components"

def test_spatial_neighbors(adata):
    """Test spatial neighbors calculation."""
    if 'spatial' not in adata.obsm:
        raise ValueError("No spatial coordinates found")
    
    import squidpy as sq
    adata_test = adata.copy()
    
    # Spatial neighbors - memory scales with O(n^2)
    sq.gr.spatial_neighbors(adata_test, coord_type='generic', n_neighs=6)
    
    return f"Spatial neighbors for {adata_test.n_obs} spots"

def test_spatial_autocorr(adata):
    """Test spatial autocorrelation."""
    if 'spatial' not in adata.obsm:
        raise ValueError("No spatial coordinates found")
    
    import squidpy as sq
    adata_test = adata.copy()
    
    # Ensure spatial neighbors exist
    if 'spatial_connectivities' not in adata_test.obsp:
        sq.gr.spatial_neighbors(adata_test, coord_type='generic', n_neighs=6)
    
    # Test with limited genes to control memory
    n_test_genes = min(50, adata_test.n_vars)
    test_genes = adata_test.var_names[:n_test_genes]
    
    sq.gr.spatial_autocorr(adata_test, genes=test_genes)
    
    return f"Spatial autocorrelation for {n_test_genes} genes"

def test_clustering_memory(adata):
    """Test clustering algorithms."""
    adata_test = adata.copy()
    
    # Ensure we have PCA
    if 'X_pca' not in adata_test.obsm:
        n_comps = min(30, adata_test.n_vars - 1, adata_test.n_obs - 1)
        sc.tl.pca(adata_test, n_comps=n_comps)
    
    # Neighbors graph - memory intensive
    n_neighbors = min(15, adata_test.n_obs - 1)
    sc.pp.neighbors(adata_test, n_neighbors=n_neighbors, n_pcs=30)
    
    # Clustering
    sc.tl.leiden(adata_test, resolution=0.5)
    
    return f"Clustering with {n_neighbors} neighbors"

def test_umap_memory(adata):
    """Test UMAP calculation."""
    adata_test = adata.copy()
    
    # Ensure neighbors exist
    if 'connectivities' not in adata_test.obsp:
        if 'X_pca' not in adata_test.obsm:
            n_comps = min(30, adata_test.n_vars - 1, adata_test.n_obs - 1)
            sc.tl.pca(adata_test, n_comps=n_comps)
        sc.pp.neighbors(adata_test, n_neighbors=15)
    
    # UMAP - memory intensive
    sc.tl.umap(adata_test)
    
    return f"UMAP for {adata_test.n_obs} cells"

def test_visualization_memory(adata):
    """Test visualization memory."""
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    
    adata_test = adata.copy()
    
    # Create a simple spatial plot
    if 'spatial' in adata_test.obsm:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        
        # Plot first gene
        sc.pl.spatial(adata_test, color=adata_test.var_names[0], ax=ax, show=False)
        plt.close(fig)
        
        return f"Spatial plot for {adata_test.n_obs} spots"
    else:
        raise ValueError("No spatial coordinates")

def main():
    """Main memory profiling function."""
    print("🧠 ChatSpatial Quick Memory Analysis")
    print("=" * 50)
    
    data_dir = Path('/Users/apple/Research/SpatialTrans_MCP/chatspatial/data')
    
    # Test datasets
    datasets = {}
    
    # Small dataset
    small_path = data_dir / 'spatial_datasets/seqfish_developmental.h5ad'
    if small_path.exists():
        adata = sc.read_h5ad(small_path)
        datasets['seqfish_1200x400'] = adata
    
    # Medium dataset
    medium_path = data_dir / 'benchmark_datasets/benchmark_500x1k.h5ad'
    if medium_path.exists():
        adata = sc.read_h5ad(medium_path)
        datasets['benchmark_500x1k'] = adata
    
    # Create synthetic larger dataset
    print("Creating synthetic 3K dataset...")
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, size=(3000, 1500))
    adata_large = sc.AnnData(X=X.astype(np.float32))
    adata_large.var_names = [f"Gene_{i:04d}" for i in range(1500)]
    adata_large.obsm['spatial'] = np.random.uniform(0, 100, size=(3000, 2))
    datasets['synthetic_3kx1.5k'] = adata_large
    
    # Memory test functions
    tests = [
        ('Preprocessing', test_preprocessing_memory),
        ('Spatial Neighbors', test_spatial_neighbors), 
        ('Spatial Autocorr', test_spatial_autocorr),
        ('Clustering', test_clustering_memory),
        ('UMAP', test_umap_memory),
        ('Visualization', test_visualization_memory)
    ]
    
    all_results = {}
    
    # Run tests
    for dataset_name, adata in datasets.items():
        print(f"\n📊 Dataset: {dataset_name}")
        print(f"   Shape: {adata.n_obs} cells × {adata.n_vars} genes")
        
        dataset_results = {}
        
        for test_name, test_func in tests:
            result = memory_test(test_name, test_func, adata)
            dataset_results[test_name] = result
            
        all_results[dataset_name] = dataset_results
    
    # Summary
    print("\n" + "="*60)
    print("📊 MEMORY ANALYSIS SUMMARY")
    print("="*60)
    
    for dataset_name, results in all_results.items():
        print(f"\n🗂️  {dataset_name}:")
        
        # Sort by memory usage
        sorted_tests = sorted(results.items(), key=lambda x: x[1]['delta_mb'], reverse=True)
        
        for test_name, result in sorted_tests:
            if result['success']:
                status = "✅"
                memory_info = f"{result['delta_mb']:+.1f}MB ({result['peak_mb']:.1f}MB peak)"
                time_info = f"{result['time_s']:.1f}s"
            else:
                status = "❌"
                memory_info = "FAILED"
                time_info = "--"
            
            print(f"  {status} {test_name:<18} {memory_info:<20} {time_info}")
    
    # Memory intensity ranking
    print(f"\n🏆 MEMORY INTENSITY RANKING (across all datasets):")
    print("-" * 50)
    
    # Calculate average memory usage per test
    test_avg_memory = {}
    for test_name, _ in tests:
        total_memory = 0
        success_count = 0
        for dataset_results in all_results.values():
            if dataset_results[test_name]['success']:
                total_memory += dataset_results[test_name]['delta_mb']
                success_count += 1
        if success_count > 0:
            test_avg_memory[test_name] = total_memory / success_count
    
    # Sort by average memory usage
    for i, (test_name, avg_memory) in enumerate(sorted(test_avg_memory.items(), key=lambda x: x[1], reverse=True), 1):
        if avg_memory > 50:
            risk = "🔴 HIGH"
        elif avg_memory > 20:
            risk = "🟡 MEDIUM"
        else:
            risk = "🟢 LOW"
        print(f"{i}. {test_name:<18} {avg_memory:>6.1f}MB avg  {risk}")

if __name__ == "__main__":
    main()