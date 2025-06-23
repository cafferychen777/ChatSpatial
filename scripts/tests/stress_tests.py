#!/usr/bin/env python
"""
Stress tests for spatial transcriptomics methods.
Tests edge cases, error handling, and robustness.
"""

import numpy as np
import pandas as pd
import anndata as ad
import warnings
warnings.filterwarnings('ignore')

def test_paste_edge_cases():
    """Test PASTE with edge cases."""
    from chatspatial.tools.spatial_registration import register_spatial_slices
    
    print("\n🔬 PASTE EDGE CASE TESTS")
    print("="*50)
    
    # Test 1: Empty gene overlap
    print("\n1. Testing with no common genes:")
    try:
        slice1 = ad.AnnData(X=np.random.rand(50, 30))
        slice2 = ad.AnnData(X=np.random.rand(50, 30))
        slice1.var_names = [f'gene_A_{i}' for i in range(30)]
        slice2.var_names = [f'gene_B_{i}' for i in range(30)]  # No overlap
        slice1.obsm['spatial'] = np.random.rand(50, 2)
        slice2.obsm['spatial'] = np.random.rand(50, 2)
        
        result = register_spatial_slices([slice1, slice2], method='paste')
        print("   ❌ Should have failed but didn't")
    except Exception as e:
        print(f"   ✅ Correctly failed: {type(e).__name__}")
    
    # Test 2: Single spot
    print("\n2. Testing with single spot:")
    try:
        slice1 = ad.AnnData(X=np.random.rand(1, 30))
        slice2 = ad.AnnData(X=np.random.rand(1, 30))
        slice1.var_names = [f'gene_{i}' for i in range(30)]
        slice2.var_names = [f'gene_{i}' for i in range(30)]
        slice1.obsm['spatial'] = np.random.rand(1, 2)
        slice2.obsm['spatial'] = np.random.rand(1, 2)
        
        result = register_spatial_slices([slice1, slice2], method='paste')
        print("   ✅ Handled single spot case")
    except Exception as e:
        print(f"   ⚠️  Failed with: {type(e).__name__}")
    
    # Test 3: Identical slices
    print("\n3. Testing with identical slices:")
    try:
        slice1 = ad.AnnData(X=np.random.rand(50, 30))
        slice1.var_names = [f'gene_{i}' for i in range(30)]
        slice1.obsm['spatial'] = np.random.rand(50, 2)
        slice2 = slice1.copy()
        
        result = register_spatial_slices([slice1, slice2], method='paste')
        print("   ✅ Handled identical slices")
    except Exception as e:
        print(f"   ❌ Failed with: {type(e).__name__}")
    
    # Test 4: Very sparse data
    print("\n4. Testing with extremely sparse data (99% zeros):")
    try:
        slice1 = ad.AnnData(X=np.random.choice([0, 1], size=(100, 50), p=[0.99, 0.01]).astype(float))
        slice2 = ad.AnnData(X=np.random.choice([0, 1], size=(100, 50), p=[0.99, 0.01]).astype(float))
        slice1.var_names = [f'gene_{i}' for i in range(50)]
        slice2.var_names = [f'gene_{i}' for i in range(50)]
        slice1.obsm['spatial'] = np.random.rand(100, 2)
        slice2.obsm['spatial'] = np.random.rand(100, 2)
        
        result = register_spatial_slices([slice1, slice2], method='paste')
        print("   ✅ Handled sparse data")
    except Exception as e:
        print(f"   ⚠️  Failed with: {type(e).__name__}")


def test_spotlight_edge_cases():
    """Test SPOTlight with edge cases."""
    from chatspatial.tools.deconvolution import deconvolve_spotlight
    
    print("\n🔬 SPOTLIGHT EDGE CASE TESTS")
    print("="*50)
    
    # Test 1: Single cell type
    print("\n1. Testing with single cell type:")
    try:
        sc_data = ad.AnnData(X=np.random.poisson(5, size=(100, 50)))
        sc_data.obs['cell_type'] = 'TypeA'  # All same type
        spatial_data = ad.AnnData(X=np.random.poisson(5, size=(20, 50)))
        
        props, _ = deconvolve_spotlight(spatial_data, sc_data)
        print(f"   ✅ Handled single cell type, shape: {props.shape}")
    except Exception as e:
        print(f"   ⚠️  Failed with: {type(e).__name__}")
    
    # Test 2: No common genes
    print("\n2. Testing with no common genes:")
    try:
        sc_data = ad.AnnData(X=np.random.poisson(5, size=(100, 50)))
        sc_data.obs['cell_type'] = np.random.choice(['A', 'B'], 100)
        sc_data.var_names = [f'gene_sc_{i}' for i in range(50)]
        
        spatial_data = ad.AnnData(X=np.random.poisson(5, size=(20, 50)))
        spatial_data.var_names = [f'gene_sp_{i}' for i in range(50)]
        
        props, _ = deconvolve_spotlight(spatial_data, sc_data)
        print("   ❌ Should have failed but didn't")
    except Exception as e:
        print(f"   ✅ Correctly failed: {type(e).__name__}")
    
    # Test 3: Very few cells per type
    print("\n3. Testing with very few cells per type (n=2):")
    try:
        sc_data = ad.AnnData(X=np.random.poisson(5, size=(10, 50)))
        sc_data.obs['cell_type'] = ['A', 'A', 'B', 'B', 'C', 'C', 'D', 'D', 'E', 'E']
        spatial_data = ad.AnnData(X=np.random.poisson(5, size=(20, 50)))
        
        props, _ = deconvolve_spotlight(spatial_data, sc_data)
        print(f"   ✅ Handled few cells per type, shape: {props.shape}")
    except Exception as e:
        print(f"   ⚠️  Failed with: {type(e).__name__}")
    
    # Test 4: All zeros in spatial data
    print("\n4. Testing with all-zero spatial data:")
    try:
        sc_data = ad.AnnData(X=np.random.poisson(5, size=(100, 50)))
        sc_data.obs['cell_type'] = np.random.choice(['A', 'B', 'C'], 100)
        spatial_data = ad.AnnData(X=np.zeros((20, 50)))
        
        props, _ = deconvolve_spotlight(spatial_data, sc_data)
        print(f"   ✅ Handled all-zero data, shape: {props.shape}")
    except Exception as e:
        print(f"   ⚠️  Failed with: {type(e).__name__}")


def test_spatialDE_edge_cases():
    """Test SpatialDE with edge cases."""
    from chatspatial.tools.spatial_statistics import find_spatial_variable_genes
    
    print("\n🔬 SPATIALDES EDGE CASE TESTS")
    print("="*50)
    
    # Test 1: Constant expression
    print("\n1. Testing with constant gene expression:")
    try:
        adata = ad.AnnData(X=np.ones((50, 10)) * 5)  # All genes have constant expression
        adata.obsm['spatial'] = np.random.rand(50, 2)
        adata.var_names = [f'gene_{i}' for i in range(10)]
        
        results = find_spatial_variable_genes(adata, method='spatialDE', n_genes=5)
        print(f"   ✅ Handled constant expression, found {len(results)} genes")
    except Exception as e:
        print(f"   ⚠️  Failed with: {type(e).__name__}")
    
    # Test 2: Colinear coordinates
    print("\n2. Testing with colinear spatial coordinates:")
    try:
        adata = ad.AnnData(X=np.random.poisson(5, size=(50, 10)).astype(float))
        # All spots on a line
        adata.obsm['spatial'] = np.column_stack([np.linspace(0, 10, 50), np.zeros(50)])
        adata.var_names = [f'gene_{i}' for i in range(10)]
        
        results = find_spatial_variable_genes(adata, method='spatialDE', n_genes=5)
        print(f"   ✅ Handled colinear coordinates, found {len(results)} genes")
    except Exception as e:
        print(f"   ⚠️  Failed with: {type(e).__name__}")
    
    # Test 3: Duplicate coordinates
    print("\n3. Testing with duplicate spatial coordinates:")
    try:
        adata = ad.AnnData(X=np.random.poisson(5, size=(50, 10)).astype(float))
        # Many spots at same location
        coords = np.random.rand(10, 2)
        adata.obsm['spatial'] = np.repeat(coords, 5, axis=0)
        adata.var_names = [f'gene_{i}' for i in range(10)]
        
        results = find_spatial_variable_genes(adata, method='spatialDE', n_genes=5)
        print(f"   ✅ Handled duplicate coordinates, found {len(results)} genes")
    except Exception as e:
        print(f"   ⚠️  Failed with: {type(e).__name__}")
    
    # Test 4: Single gene
    print("\n4. Testing with single gene:")
    try:
        adata = ad.AnnData(X=np.random.poisson(5, size=(50, 1)).astype(float))
        adata.obsm['spatial'] = np.random.rand(50, 2)
        adata.var_names = ['single_gene']
        
        results = find_spatial_variable_genes(adata, method='spatialDE')
        print(f"   ✅ Handled single gene, found {len(results)} genes")
    except Exception as e:
        print(f"   ⚠️  Failed with: {type(e).__name__}")


def test_memory_efficiency():
    """Test memory efficiency with larger datasets."""
    print("\n🔬 MEMORY EFFICIENCY TESTS")
    print("="*50)
    
    import psutil
    import os
    
    process = psutil.Process(os.getpid())
    initial_memory = process.memory_info().rss / 1024 / 1024  # MB
    
    # Test with progressively larger datasets
    sizes = [(1000, 500), (2000, 1000), (5000, 2000)]
    
    for n_spots, n_genes in sizes:
        print(f"\nTesting with {n_spots} spots, {n_genes} genes:")
        
        # Create data
        adata = ad.AnnData(X=np.random.poisson(5, size=(n_spots, n_genes)).astype(float))
        adata.obsm['spatial'] = np.random.rand(n_spots, 2)
        adata.var_names = [f'gene_{i}' for i in range(n_genes)]
        
        # Add some spatial patterns
        for i in range(min(10, n_genes)):
            pattern = np.sin(adata.obsm['spatial'][:, 0] * (i + 1))
            adata.X[:, i] += pattern * 10
        
        current_memory = process.memory_info().rss / 1024 / 1024
        print(f"   Memory before: {current_memory:.1f} MB")
        
        try:
            from chatspatial.tools.spatial_statistics import find_spatial_variable_genes
            results = find_spatial_variable_genes(adata, method='spatialDE', n_genes=20)
            
            after_memory = process.memory_info().rss / 1024 / 1024
            memory_increase = after_memory - current_memory
            
            print(f"   ✅ Success! Memory increase: {memory_increase:.1f} MB")
            print(f"   Found {len(results)} spatial genes")
            
        except Exception as e:
            print(f"   ❌ Failed: {type(e).__name__}: {str(e)[:50]}")
        
        # Cleanup
        del adata
        import gc
        gc.collect()


def test_concurrent_usage():
    """Test concurrent usage of methods."""
    print("\n🔬 CONCURRENT USAGE TEST")
    print("="*50)
    
    import concurrent.futures
    import time
    
    def run_spatialDE(idx):
        from chatspatial.tools.spatial_statistics import find_spatial_variable_genes
        
        np.random.seed(idx)
        adata = ad.AnnData(X=np.random.poisson(5, size=(100, 50)).astype(float))
        adata.obsm['spatial'] = np.random.rand(100, 2)
        adata.var_names = [f'gene_{i}' for i in range(50)]
        
        # Add spatial pattern
        adata.X[:, 0] += np.sin(adata.obsm['spatial'][:, 0] * 2) * 10
        
        start = time.time()
        results = find_spatial_variable_genes(adata, method='spatialDE', n_genes=10)
        elapsed = time.time() - start
        
        return idx, len(results), elapsed
    
    print("Running 5 SpatialDE analyses concurrently...")
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(run_spatialDE, i) for i in range(5)]
        
        for future in concurrent.futures.as_completed(futures):
            idx, n_results, elapsed = future.result()
            print(f"   Thread {idx}: Found {n_results} genes in {elapsed:.2f}s")
    
    print("✅ Concurrent execution completed successfully")


def main():
    """Run all stress tests."""
    print("🔥 RUNNING STRESS TESTS FOR SPATIAL METHODS")
    print("Testing edge cases, error handling, and robustness")
    
    # Edge case tests
    test_paste_edge_cases()
    test_spotlight_edge_cases()
    test_spatialDE_edge_cases()
    
    # Performance tests
    test_memory_efficiency()
    test_concurrent_usage()
    
    print("\n" + "="*50)
    print("✅ STRESS TESTS COMPLETED")
    print("="*50)


if __name__ == "__main__":
    main()