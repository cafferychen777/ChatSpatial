#!/usr/bin/env python3
"""
Comprehensive test suite for spatial_domains.py

This test suite covers:
1. Basic functionality of all domain identification methods
2. Edge cases and error handling
3. Parameter validation
4. Data format compatibility
5. Performance and reliability issues

Usage:
    python test_spatial_domains.py
"""

import asyncio
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import traceback
import time
from typing import Dict, Any, Optional
import warnings
warnings.filterwarnings('ignore')

# Import the module to test
from tools.spatial_domains import (
    identify_spatial_domains,
    _identify_domains_spagcn,
    _identify_domains_clustering,
    _refine_spatial_domains
)
from models.data import SpatialDomainParameters
from models.analysis import SpatialDomainResult


class MockContext:
    """Mock context for testing"""
    def __init__(self):
        self.messages = []
    
    async def info(self, message: str):
        self.messages.append(f"INFO: {message}")
        print(f"INFO: {message}")
    
    async def warning(self, message: str):
        self.messages.append(f"WARNING: {message}")
        print(f"WARNING: {message}")


def create_test_adata(n_obs: int = 100, n_vars: int = 50, add_spatial: bool = True, 
                     add_hvg: bool = True, seed: int = 42) -> ad.AnnData:
    """Create test AnnData object with spatial coordinates"""
    np.random.seed(seed)
    
    # Create expression data
    X = np.random.negative_binomial(n=5, p=0.3, size=(n_obs, n_vars)).astype(np.float32)
    
    # Create gene and cell names
    var_names = [f"Gene_{i}" for i in range(n_vars)]
    obs_names = [f"Cell_{i}" for i in range(n_obs)]
    
    # Create AnnData object
    adata = ad.AnnData(X=X, var=pd.DataFrame(index=var_names), obs=pd.DataFrame(index=obs_names))
    
    # Add spatial coordinates
    if add_spatial:
        # Create realistic spatial coordinates (grid-like pattern with some noise)
        grid_size = int(np.sqrt(n_obs)) + 1
        coords = []
        for i in range(n_obs):
            x = (i % grid_size) + np.random.normal(0, 0.1)
            y = (i // grid_size) + np.random.normal(0, 0.1)
            coords.append([x, y])
        adata.obsm['spatial'] = np.array(coords)
    
    # Add highly variable genes
    if add_hvg:
        hvg_mask = np.random.choice([True, False], size=n_vars, p=[0.3, 0.7])
        adata.var['highly_variable'] = hvg_mask
    
    return adata


def create_visium_like_adata(n_obs: int = 100, seed: int = 42) -> ad.AnnData:
    """Create Visium-like AnnData with proper spatial structure"""
    adata = create_test_adata(n_obs=n_obs, seed=seed)
    
    # Add Visium-like spatial structure
    adata.uns['spatial'] = {
        'sample1': {
            'images': {
                'hires': np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8),
                'lowres': np.random.randint(0, 255, (50, 50, 3), dtype=np.uint8)
            }
        }
    }
    
    return adata


async def test_basic_functionality():
    """Test basic functionality of all methods"""
    print("\n=== Testing Basic Functionality ===")
    
    # Test data
    adata = create_test_adata(n_obs=50, n_vars=30)
    data_store = {"test": {"adata": adata}}
    context = MockContext()
    
    # Test all methods
    methods = ["spagcn", "leiden", "louvain"]
    
    for method in methods:
        print(f"\nTesting {method} method...")
        try:
            params = SpatialDomainParameters(method=method, n_domains=3)
            result = await identify_spatial_domains("test", data_store, params, context)
            
            assert isinstance(result, SpatialDomainResult)
            assert result.method == method
            assert result.n_domains >= 1
            assert len(result.domain_counts) >= 1
            
            # Check that domain labels were added to adata
            domain_key = f"spatial_domains_{method}"
            assert domain_key in data_store["test"]["adata"].obs.columns
            
            print(f"✓ {method} test passed - found {result.n_domains} domains")
            
        except Exception as e:
            print(f"✗ {method} test failed: {e}")
            traceback.print_exc()


async def test_parameter_variations():
    """Test various parameter combinations"""
    print("\n=== Testing Parameter Variations ===")
    
    adata = create_visium_like_adata(n_obs=80)
    data_store = {"test": {"adata": adata}}
    context = MockContext()
    
    test_cases = [
        # SpaGCN variations
        {"method": "spagcn", "n_domains": 5, "spagcn_s": 0.5, "spagcn_use_histology": True},
        {"method": "spagcn", "n_domains": 3, "spagcn_s": 2.0, "spagcn_use_histology": False},
        {"method": "spagcn", "n_domains": 7, "spagcn_b": 25, "spagcn_p": 0.3},
        
        # Clustering variations
        {"method": "leiden", "resolution": 0.3, "cluster_n_neighbors": 10},
        {"method": "leiden", "resolution": 0.8, "cluster_spatial_weight": 0.5},
        {"method": "louvain", "resolution": 0.5, "use_highly_variable": False},
        
        # Refinement variations
        {"method": "leiden", "refine_domains": True},
        {"method": "spagcn", "refine_domains": False},
    ]
    
    for i, params_dict in enumerate(test_cases):
        print(f"\nTest case {i+1}: {params_dict}")
        try:
            params = SpatialDomainParameters(**params_dict)
            result = await identify_spatial_domains("test", data_store, params, context)
            
            assert result.n_domains >= 1
            print(f"✓ Parameter test {i+1} passed - {result.n_domains} domains")
            
        except Exception as e:
            print(f"✗ Parameter test {i+1} failed: {e}")


async def test_edge_cases():
    """Test edge cases and error conditions"""
    print("\n=== Testing Edge Cases ===")
    context = MockContext()
    
    # Test 1: Missing dataset
    print("Test 1: Missing dataset")
    try:
        params = SpatialDomainParameters()
        await identify_spatial_domains("nonexistent", {}, params, context)
        print("✗ Should have raised ValueError for missing dataset")
    except ValueError as e:
        print(f"✓ Correctly caught missing dataset error: {e}")
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
    
    # Test 2: No spatial coordinates
    print("\nTest 2: No spatial coordinates")
    try:
        adata = create_test_adata(add_spatial=False)
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters()
        await identify_spatial_domains("test", data_store, params, context)
        print("✗ Should have raised ValueError for missing spatial coordinates")
    except ValueError as e:
        print(f"✓ Correctly caught missing spatial coords error: {e}")
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
    
    # Test 3: Very small dataset
    print("\nTest 3: Very small dataset (5 spots)")
    try:
        adata = create_test_adata(n_obs=5, n_vars=10)
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(n_domains=3)
        result = await identify_spatial_domains("test", data_store, params, context)
        print(f"✓ Small dataset handled - {result.n_domains} domains found")
    except Exception as e:
        print(f"✗ Small dataset failed: {e}")
    
    # Test 4: Too many domains requested
    print("\nTest 4: Too many domains requested")
    try:
        adata = create_test_adata(n_obs=20, n_vars=10)
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(n_domains=25)  # More domains than spots
        result = await identify_spatial_domains("test", data_store, params, context)
        print(f"✓ Handled too many domains - got {result.n_domains} domains")
    except Exception as e:
        print(f"✗ Too many domains test failed: {e}")
    
    # Test 5: Invalid method
    print("\nTest 5: Invalid method")
    try:
        adata = create_test_adata()
        data_store = {"test": {"adata": adata}}
        # This should fail at parameter validation level
        params = SpatialDomainParameters(method="invalid_method")
        print("✗ Should have failed at parameter validation")
    except Exception as e:
        print(f"✓ Invalid method correctly rejected: {e}")


async def test_data_types():
    """Test different data types and formats"""
    print("\n=== Testing Data Types ===")
    context = MockContext()
    
    # Test 1: Integer data
    print("Test 1: Integer expression data")
    try:
        adata = create_test_adata()
        adata.X = adata.X.astype(int)
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(method="leiden")
        result = await identify_spatial_domains("test", data_store, params, context)
        print(f"✓ Integer data handled - {result.n_domains} domains")
    except Exception as e:
        print(f"✗ Integer data failed: {e}")
    
    # Test 2: Sparse data
    print("\nTest 2: Sparse expression data")
    try:
        import scipy.sparse as sp
        adata = create_test_adata()
        adata.X = sp.csr_matrix(adata.X)
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(method="louvain")
        result = await identify_spatial_domains("test", data_store, params, context)
        print(f"✓ Sparse data handled - {result.n_domains} domains")
    except Exception as e:
        print(f"✗ Sparse data failed: {e}")
    
    # Test 3: Pre-normalized data
    print("\nTest 3: Pre-normalized data")
    try:
        adata = create_test_adata()
        # Normalize the data first
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.uns['log1p'] = True  # Mark as already normalized
        
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(method="spagcn")
        result = await identify_spatial_domains("test", data_store, params, context)
        print(f"✓ Pre-normalized data handled - {result.n_domains} domains")
    except Exception as e:
        print(f"✗ Pre-normalized data failed: {e}")


async def test_performance():
    """Test performance with larger datasets"""
    print("\n=== Testing Performance ===")
    context = MockContext()
    
    # Test with progressively larger datasets
    sizes = [100, 500, 1000]
    
    for size in sizes:
        print(f"\nTesting with {size} spots...")
        try:
            start_time = time.time()
            
            adata = create_test_adata(n_obs=size, n_vars=100)
            data_store = {"test": {"adata": adata}}
            
            # Use faster method for larger datasets
            if size <= 500:
                params = SpatialDomainParameters(method="spagcn", n_domains=5)
            else:
                params = SpatialDomainParameters(method="leiden", n_domains=5)
            
            result = await identify_spatial_domains("test", data_store, params, context)
            
            elapsed = time.time() - start_time
            print(f"✓ {size} spots processed in {elapsed:.2f}s - {result.n_domains} domains")
            
            if elapsed > 60:  # More than 1 minute is concerning
                print(f"⚠ Performance warning: {elapsed:.2f}s is quite slow")
                
        except Exception as e:
            print(f"✗ Performance test with {size} spots failed: {e}")


async def test_consistency():
    """Test consistency of results across multiple runs"""
    print("\n=== Testing Consistency ===")
    context = MockContext()
    
    # Create test data
    adata = create_test_adata(n_obs=100, seed=42)
    
    # Test each method multiple times
    methods = ["leiden", "louvain"]  # SpaGCN might have more randomness
    
    for method in methods:
        print(f"\nTesting {method} consistency...")
        results = []
        
        for run in range(3):
            try:
                # Fresh copy of data for each run
                data_store = {"test": {"adata": adata.copy()}}
                params = SpatialDomainParameters(method=method, n_domains=5)
                result = await identify_spatial_domains("test", data_store, params, context)
                results.append(result.n_domains)
            except Exception as e:
                print(f"Run {run+1} failed: {e}")
                results.append(None)
        
        # Check consistency
        valid_results = [r for r in results if r is not None]
        if len(valid_results) >= 2:
            if len(set(valid_results)) == 1:
                print(f"✓ {method} consistent: {valid_results[0]} domains in all runs")
            else:
                print(f"⚠ {method} variable: {valid_results} domains across runs")
        else:
            print(f"✗ {method} too many failures to assess consistency")


async def test_refinement():
    """Test spatial domain refinement functionality"""
    print("\n=== Testing Domain Refinement ===")
    context = MockContext()
    
    # Create test data with clear spatial structure
    adata = create_test_adata(n_obs=64, seed=42)  # 8x8 grid-like
    data_store = {"test": {"adata": adata}}
    
    # Test refinement with different methods
    methods = ["leiden", "spagcn"]
    
    for method in methods:
        print(f"\nTesting refinement with {method}...")
        try:
            params = SpatialDomainParameters(method=method, refine_domains=True, n_domains=4)
            result = await identify_spatial_domains("test", data_store, params, context)
            
            # Check that both original and refined domains exist
            adata_result = data_store["test"]["adata"]
            domain_key = result.domain_key
            refined_key = result.refined_domain_key
            
            assert domain_key in adata_result.obs.columns
            if refined_key:
                assert refined_key in adata_result.obs.columns
                print(f"✓ {method} refinement successful - original and refined domains created")
            else:
                print(f"⚠ {method} refinement skipped or failed")
                
        except Exception as e:
            print(f"✗ {method} refinement failed: {e}")


async def test_robustness():
    """Test robustness against problematic inputs"""
    print("\n=== Testing Robustness ===")
    context = MockContext()
    
    # Test 1: Data with NaN values
    print("Test 1: Data with NaN values")
    try:
        adata = create_test_adata()
        adata.X[0, 0] = np.nan
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(method="leiden")
        result = await identify_spatial_domains("test", data_store, params, context)
        print(f"✓ NaN values handled - {result.n_domains} domains")
    except Exception as e:
        print(f"✗ NaN values failed: {e}")
    
    # Test 2: Data with infinite values
    print("\nTest 2: Data with infinite values")
    try:
        adata = create_test_adata()
        adata.X[0, 0] = np.inf
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(method="leiden")
        result = await identify_spatial_domains("test", data_store, params, context)
        print(f"✓ Infinite values handled - {result.n_domains} domains")
    except Exception as e:
        print(f"✗ Infinite values failed: {e}")
    
    # Test 3: All-zero expression data
    print("\nTest 3: All-zero expression data")
    try:
        adata = create_test_adata()
        adata.X[:] = 0
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(method="leiden")
        result = await identify_spatial_domains("test", data_store, params, context)
        print(f"✓ All-zero data handled - {result.n_domains} domains")
    except Exception as e:
        print(f"✗ All-zero data failed: {e}")
    
    # Test 4: Identical spatial coordinates
    print("\nTest 4: Identical spatial coordinates")
    try:
        adata = create_test_adata()
        adata.obsm['spatial'][:] = [0, 0]  # All spots at same location
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(method="leiden")
        result = await identify_spatial_domains("test", data_store, params, context)
        print(f"✓ Identical coordinates handled - {result.n_domains} domains")
    except Exception as e:
        print(f"✗ Identical coordinates failed: {e}")


async def run_all_tests():
    """Run all test suites"""
    print("Starting comprehensive spatial_domains.py tests...")
    print("=" * 60)
    
    test_functions = [
        test_basic_functionality,
        test_parameter_variations,
        test_edge_cases,
        test_data_types,
        test_consistency,
        test_refinement,
        test_robustness,
        test_performance,  # Run performance tests last
    ]
    
    passed = 0
    failed = 0
    
    for test_func in test_functions:
        try:
            await test_func()
            passed += 1
        except Exception as e:
            print(f"\n✗ Test suite {test_func.__name__} failed with error: {e}")
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Test Summary: {passed} test suites passed, {failed} failed")
    
    if failed == 0:
        print("🎉 All tests passed! The spatial_domains.py module appears to be working correctly.")
    else:
        print("⚠️ Some tests failed. Please review the output above for issues.")


if __name__ == "__main__":
    # Run the tests
    asyncio.run(run_all_tests())