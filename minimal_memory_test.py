#!/usr/bin/env python3
"""
Minimal Memory Test for ChatSpatial Core Tools
Linus原则：最简单的测试，最有用的结果
"""

import psutil
import tracemalloc
import gc
import time
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

def test_memory_usage():
    """Test memory usage of core operations."""
    process = psutil.Process()
    
    def get_mem(): 
        return process.memory_info().rss / 1024 / 1024
    
    def test_op(name, func):
        gc.collect()
        start_mem = get_mem()
        tracemalloc.start()
        
        try:
            start_time = time.time()
            result = func()
            elapsed = time.time() - start_time
            success = True
            error = None
        except Exception as e:
            elapsed = time.time() - start_time
            success = False
            error = str(e)[:80]
            result = None
            
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        end_mem = get_mem()
        
        return {
            'success': success,
            'delta_mb': end_mem - start_mem,
            'peak_mb': peak / 1024 / 1024,
            'time_s': elapsed,
            'error': error
        }
    
    # Create test datasets
    print("Creating test datasets...")
    
    # Small dataset (real data simulation)
    X_small = np.random.negative_binomial(5, 0.3, size=(1000, 500)).astype(np.float32)
    adata_small = sc.AnnData(X_small)
    adata_small.var_names = [f"Gene_{i:03d}" for i in range(500)]
    adata_small.obsm['spatial'] = np.random.uniform(0, 100, size=(1000, 2))
    
    # Medium dataset  
    X_medium = np.random.negative_binomial(5, 0.3, size=(3000, 1000)).astype(np.float32)
    adata_medium = sc.AnnData(X_medium)
    adata_medium.var_names = [f"Gene_{i:04d}" for i in range(1000)]
    adata_medium.obsm['spatial'] = np.random.uniform(0, 100, size=(3000, 2))
    
    # Large dataset
    X_large = np.random.negative_binomial(5, 0.3, size=(8000, 2000)).astype(np.float32)
    adata_large = sc.AnnData(X_large)
    adata_large.var_names = [f"Gene_{i:04d}" for i in range(2000)]
    adata_large.obsm['spatial'] = np.random.uniform(0, 100, size=(8000, 2))
    
    datasets = {
        'Small_1Kx500': adata_small,
        'Medium_3Kx1K': adata_medium,
        'Large_8Kx2K': adata_large
    }
    
    print(f"\n{'='*60}")
    print("🧠 MEMORY ANALYSIS RESULTS")
    print(f"{'='*60}")
    
    all_results = {}
    
    for name, adata in datasets.items():
        print(f"\n📊 Dataset: {name} ({adata.n_obs} cells × {adata.n_vars} genes)")
        
        results = {}
        
        # Test 1: Basic preprocessing
        def test_preprocessing():
            adata_test = adata.copy()
            sc.pp.normalize_total(adata_test, target_sum=1e4)
            sc.pp.log1p(adata_test)
            sc.pp.highly_variable_genes(adata_test, n_top_genes=min(1000, adata_test.n_vars//2))
            sc.pp.scale(adata_test, max_value=10)
            return "preprocessing"
        
        results['Preprocessing'] = test_op('Preprocessing', test_preprocessing)
        
        # Test 2: PCA
        def test_pca():
            adata_test = adata.copy()
            sc.pp.normalize_total(adata_test, target_sum=1e4)
            sc.pp.log1p(adata_test)
            sc.pp.scale(adata_test, max_value=10)
            n_comps = min(50, adata_test.n_vars-1, adata_test.n_obs-1)
            sc.tl.pca(adata_test, n_comps=n_comps)
            return f"PCA {n_comps} comps"
        
        results['PCA'] = test_op('PCA', test_pca)
        
        # Test 3: Neighbors
        def test_neighbors():
            adata_test = adata.copy()
            sc.pp.normalize_total(adata_test, target_sum=1e4)
            sc.pp.log1p(adata_test)
            sc.tl.pca(adata_test, n_comps=30)
            sc.pp.neighbors(adata_test, n_neighbors=15, n_pcs=30)
            return "neighbors"
        
        results['Neighbors'] = test_op('Neighbors', test_neighbors)
        
        # Test 4: UMAP  
        def test_umap():
            adata_test = adata.copy()
            sc.pp.normalize_total(adata_test, target_sum=1e4)
            sc.pp.log1p(adata_test)
            sc.tl.pca(adata_test, n_comps=30)
            sc.pp.neighbors(adata_test, n_neighbors=15)
            sc.tl.umap(adata_test)
            return "UMAP"
        
        results['UMAP'] = test_op('UMAP', test_umap)
        
        # Test 5: Spatial neighbors
        def test_spatial():
            import squidpy as sq
            adata_test = adata.copy()
            sq.gr.spatial_neighbors(adata_test, coord_type='generic', n_neighs=6)
            return "spatial neighbors"
        
        results['Spatial_Neighbors'] = test_op('Spatial_Neighbors', test_spatial)
        
        all_results[name] = results
        
        # Print results for this dataset
        for test_name, result in results.items():
            if result['success']:
                print(f"  ✅ {test_name:<18} {result['delta_mb']:+6.1f}MB  {result['peak_mb']:6.1f}MB peak  {result['time_s']:5.1f}s")
            else:
                print(f"  ❌ {test_name:<18} FAILED: {result['error']}")
    
    # Overall analysis
    print(f"\n{'='*60}")
    print("🏆 MEMORY INTENSITY RANKING")
    print(f"{'='*60}")
    
    # Calculate average memory per operation
    op_stats = {}
    for dataset_name, results in all_results.items():
        for op_name, result in results.items():
            if result['success']:
                if op_name not in op_stats:
                    op_stats[op_name] = {'deltas': [], 'peaks': [], 'times': []}
                op_stats[op_name]['deltas'].append(result['delta_mb'])
                op_stats[op_name]['peaks'].append(result['peak_mb'])
                op_stats[op_name]['times'].append(result['time_s'])
    
    # Calculate averages and sort by peak memory
    op_averages = []
    for op_name, stats in op_stats.items():
        if stats['deltas']:
            avg_delta = np.mean(stats['deltas'])
            avg_peak = np.mean(stats['peaks'])
            avg_time = np.mean(stats['times'])
            op_averages.append((op_name, avg_delta, avg_peak, avg_time))
    
    # Sort by peak memory usage
    op_averages.sort(key=lambda x: x[2], reverse=True)
    
    for i, (op_name, avg_delta, avg_peak, avg_time) in enumerate(op_averages, 1):
        if avg_peak > 100:
            risk = "🔴 HIGH"
        elif avg_peak > 50:
            risk = "🟡 MEDIUM"
        else:
            risk = "🟢 LOW"
        
        print(f"{i}. {op_name:<18} avg: {avg_delta:+6.1f}MB  peak: {avg_peak:6.1f}MB  {avg_time:5.1f}s  {risk}")
    
    # Scalability analysis
    print(f"\n{'='*60}")
    print("📈 SCALABILITY ANALYSIS")
    print(f"{'='*60}")
    
    print("Memory scaling from Small → Medium → Large datasets:")
    
    for op_name in ['Preprocessing', 'PCA', 'Neighbors', 'UMAP', 'Spatial_Neighbors']:
        if op_name in all_results['Small_1Kx500'] and op_name in all_results['Large_8Kx2K']:
            small_peak = all_results['Small_1Kx500'][op_name]['peak_mb']
            large_peak = all_results['Large_8Kx2K'][op_name]['peak_mb']
            
            if small_peak > 0 and large_peak > 0:
                scaling_factor = large_peak / small_peak
                if scaling_factor > 10:
                    scaling_risk = "🔴 POOR"
                elif scaling_factor > 5:
                    scaling_risk = "🟡 MODERATE" 
                else:
                    scaling_risk = "🟢 GOOD"
                
                print(f"{op_name:<18} {small_peak:6.1f}MB → {large_peak:6.1f}MB  ({scaling_factor:4.1f}x)  {scaling_risk}")

if __name__ == "__main__":
    test_memory_usage()