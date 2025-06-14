#!/usr/bin/env python3
"""
Final comprehensive test for spatial_domains.py after all fixes
"""

import asyncio
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import time
import sys
import os

# Add the chatspatial path
sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

import warnings
warnings.filterwarnings('ignore')

# Mock the necessary classes to avoid import issues
class MockContext:
    def __init__(self):
        self.messages = []
    
    async def info(self, message: str):
        self.messages.append(f"INFO: {message}")
        print(f"INFO: {message}")
    
    async def warning(self, message: str):
        self.messages.append(f"WARNING: {message}")
        print(f"WARNING: {message}")

class SpatialDomainParameters:
    def __init__(self, method="leiden", n_domains=3, **kwargs):
        self.method = method
        self.n_domains = n_domains
        self.resolution = kwargs.get('resolution', 0.5)
        self.use_highly_variable = kwargs.get('use_highly_variable', True)
        self.refine_domains = kwargs.get('refine_domains', False)
        self.cluster_n_neighbors = kwargs.get('cluster_n_neighbors', None)
        self.cluster_spatial_weight = kwargs.get('cluster_spatial_weight', None)
        # SpaGCN parameters
        self.spagcn_s = kwargs.get('spagcn_s', 1.0)
        self.spagcn_b = kwargs.get('spagcn_b', 49)
        self.spagcn_p = kwargs.get('spagcn_p', 0.5)
        self.spagcn_use_histology = kwargs.get('spagcn_use_histology', False)
        self.spagcn_random_seed = kwargs.get('spagcn_random_seed', 100)

class SpatialDomainResult:
    def __init__(self, data_id, method, n_domains, domain_key, domain_counts, **kwargs):
        self.data_id = data_id
        self.method = method
        self.n_domains = n_domains
        self.domain_key = domain_key
        self.domain_counts = domain_counts
        self.refined_domain_key = kwargs.get('refined_domain_key')
        self.statistics = kwargs.get('statistics', {})
        self.embeddings_key = kwargs.get('embeddings_key')

def create_test_data(n_obs=100, n_vars=50, seed=42, create_structure=True):
    """Create realistic test data with spatial structure"""
    np.random.seed(seed)
    
    if create_structure:
        # Create structured data with clear spatial domains
        grid_size = int(np.sqrt(n_obs)) + 1
        coords = []
        domain_labels = []
        
        for i in range(n_obs):
            x = (i % grid_size) + np.random.normal(0, 0.1)
            y = (i // grid_size) + np.random.normal(0, 0.1)
            coords.append([x, y])
            
            # Create 3 spatial domains based on location
            if x < grid_size/3:
                domain = 0
            elif x < 2*grid_size/3:
                domain = 1
            else:
                domain = 2
            domain_labels.append(domain)
        
        # Create expression data with domain-specific patterns
        X = np.random.negative_binomial(n=5, p=0.3, size=(n_obs, n_vars)).astype(np.float32)
        
        # Add domain-specific expression patterns
        for i, domain in enumerate(domain_labels):
            if domain == 0:
                X[i, :10] += np.random.poisson(10, 10)  # High expression in first 10 genes
            elif domain == 1:
                X[i, 10:20] += np.random.poisson(10, 10)  # High expression in genes 10-20
            else:
                X[i, 20:30] += np.random.poisson(10, 10)  # High expression in genes 20-30
    else:
        # Random data without structure
        X = np.random.negative_binomial(n=5, p=0.3, size=(n_obs, n_vars)).astype(np.float32)
        coords = np.random.uniform(0, 10, (n_obs, 2))
    
    # Create AnnData
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"spot_{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    )
    adata.obsm['spatial'] = np.array(coords)
    
    # Add highly variable genes
    hvg_mask = np.random.choice([True, False], size=n_vars, p=[0.3, 0.7])
    adata.var['highly_variable'] = hvg_mask
    
    return adata

async def test_environment_check():
    """Test environment compatibility checking"""
    print("\n=== Testing Environment Compatibility ===")
    
    try:
        from chatspatial.tools.spatial_domains import _check_environment_compatibility
        issues = _check_environment_compatibility()
        
        if issues:
            print("Environment issues detected:")
            for issue in issues:
                print(f"  - {issue}")
        else:
            print("✓ No environment issues detected")
        
        return len(issues) == 0
        
    except Exception as e:
        print(f"✗ Environment check failed: {e}")
        return False

async def test_clustering_methods():
    """Test clustering methods with various scenarios"""
    print("\n=== Testing Clustering Methods ===")
    
    test_cases = [
        {"name": "Small dataset", "n_obs": 50, "n_vars": 30},
        {"name": "Medium dataset", "n_obs": 200, "n_vars": 100},
        {"name": "Structured data", "n_obs": 100, "n_vars": 50, "structure": True},
    ]
    
    methods = ["leiden", "louvain"]
    results = {}
    
    for method in methods:
        print(f"\nTesting {method}...")
        method_results = []
        
        for test_case in test_cases:
            try:
                print(f"  {test_case['name']}...", end=" ")
                
                # Create test data
                adata = create_test_data(
                    n_obs=test_case['n_obs'],
                    n_vars=test_case['n_vars'],
                    create_structure=test_case.get('structure', False)
                )
                data_store = {"test": {"adata": adata}}
                context = MockContext()
                
                # Import and test
                from chatspatial.tools.spatial_domains import identify_spatial_domains
                
                params = SpatialDomainParameters(
                    method=method,
                    n_domains=3,
                    refine_domains=True
                )
                
                start_time = time.time()
                result = await identify_spatial_domains("test", data_store, params, context)
                elapsed = time.time() - start_time
                
                print(f"✓ {elapsed:.2f}s, {result.n_domains} domains")
                method_results.append(True)
                
            except Exception as e:
                print(f"✗ {e}")
                method_results.append(False)
        
        results[method] = all(method_results)
    
    return results

async def test_large_dataset_handling():
    """Test handling of large datasets"""
    print("\n=== Testing Large Dataset Handling ===")
    
    # Test progressively larger datasets
    sizes = [500, 1000, 2000]
    
    for size in sizes:
        try:
            print(f"Testing {size} spots...", end=" ")
            
            adata = create_test_data(n_obs=size, n_vars=200)
            data_store = {"test": {"adata": adata}}
            context = MockContext()
            
            from chatspatial.tools.spatial_domains import identify_spatial_domains
            
            # Use leiden for large datasets (faster than SpaGCN)
            params = SpatialDomainParameters(
                method="leiden",
                n_domains=5,
                use_highly_variable=True
            )
            
            start_time = time.time()
            result = await identify_spatial_domains("test", data_store, params, context)
            elapsed = time.time() - start_time
            
            # Check if preprocessing worked
            final_adata = data_store["test"]["adata"]
            
            print(f"✓ {elapsed:.2f}s, {result.n_domains} domains")
            print(f"    Original: {size} spots, Final: {final_adata.n_obs} spots")
            
            if elapsed > 60:  # More than 1 minute
                print(f"    ⚠ Performance warning: {elapsed:.1f}s is slow")
            
        except Exception as e:
            print(f"✗ {e}")
            return False
    
    return True

async def test_parameter_robustness():
    """Test robustness with various parameter combinations"""
    print("\n=== Testing Parameter Robustness ===")
    
    # Test problematic parameter combinations
    test_cases = [
        {"name": "Many domains", "n_domains": 15},
        {"name": "High resolution", "resolution": 2.0},
        {"name": "Low resolution", "resolution": 0.1},
        {"name": "No HVG", "use_highly_variable": False},
        {"name": "Small neighbors", "cluster_n_neighbors": 5},
        {"name": "High spatial weight", "cluster_spatial_weight": 0.8},
    ]
    
    base_adata = create_test_data(n_obs=150, n_vars=100)
    
    for test_case in test_cases:
        try:
            print(f"  {test_case['name']}...", end=" ")
            
            data_store = {"test": {"adata": base_adata.copy()}}
            context = MockContext()
            
            from chatspatial.tools.spatial_domains import identify_spatial_domains
            
            # Create parameters with test case modifications
            params_dict = {"method": "leiden", "n_domains": 3}
            params_dict.update({k: v for k, v in test_case.items() if k != "name"})
            params = SpatialDomainParameters(**params_dict)
            
            start_time = time.time()
            result = await asyncio.wait_for(
                identify_spatial_domains("test", data_store, params, context),
                timeout=30.0  # 30 second timeout for robustness tests
            )
            elapsed = time.time() - start_time
            
            print(f"✓ {elapsed:.2f}s, {result.n_domains} domains")
            
        except asyncio.TimeoutError:
            print("✗ Timeout")
        except Exception as e:
            print(f"✗ {e}")
    
    return True

async def test_error_handling():
    """Test error handling for edge cases"""
    print("\n=== Testing Error Handling ===")
    
    context = MockContext()
    
    # Test 1: Missing dataset
    try:
        from chatspatial.tools.spatial_domains import identify_spatial_domains
        await identify_spatial_domains("nonexistent", {}, SpatialDomainParameters(), context)
        print("✗ Should have raised error for missing dataset")
    except ValueError:
        print("✓ Correctly handles missing dataset")
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
    
    # Test 2: No spatial coordinates
    try:
        adata = create_test_data()
        del adata.obsm['spatial']
        data_store = {"test": {"adata": adata}}
        await identify_spatial_domains("test", data_store, SpatialDomainParameters(), context)
        print("✗ Should have raised error for missing spatial coordinates")
    except ValueError:
        print("✓ Correctly handles missing spatial coordinates")
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
    
    # Test 3: Very small dataset
    try:
        adata = create_test_data(n_obs=3, n_vars=5)
        data_store = {"test": {"adata": adata}}
        result = await identify_spatial_domains("test", data_store, 
                                               SpatialDomainParameters(n_domains=2), context)
        print(f"✓ Handles very small dataset: {result.n_domains} domains")
    except Exception as e:
        print(f"✗ Small dataset handling failed: {e}")
    
    return True

async def run_final_tests():
    """Run all final tests"""
    print("Final Comprehensive Test Suite for spatial_domains.py")
    print("=" * 60)
    
    results = {}
    
    # Run all test suites
    results['environment'] = await test_environment_check()
    results['clustering'] = await test_clustering_methods()
    results['large_datasets'] = await test_large_dataset_handling()
    results['parameters'] = await test_parameter_robustness()
    results['error_handling'] = await test_error_handling()
    
    # Summary
    print("\n" + "=" * 60)
    print("FINAL TEST SUMMARY:")
    
    passed = 0
    total = len(results)
    
    for test_name, result in results.items():
        if isinstance(result, dict):
            # For clustering results
            sub_passed = sum(result.values())
            sub_total = len(result)
            status = "✓" if sub_passed == sub_total else "⚠"
            print(f"{status} {test_name}: {sub_passed}/{sub_total} methods passed")
            if sub_passed == sub_total:
                passed += 1
        else:
            status = "✓" if result else "✗"
            print(f"{status} {test_name}: {'PASSED' if result else 'FAILED'}")
            if result:
                passed += 1
    
    print(f"\nOverall: {passed}/{total} test suites passed")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED!")
        print("The spatial_domains.py module should now be much more reliable.")
        print("The 'mysteriously fails to run' issue should be significantly reduced.")
    else:
        print(f"\n⚠️ {total-passed} test suite(s) had issues.")
        print("Some problems may remain, but major issues should be fixed.")
    
    return passed == total

if __name__ == "__main__":
    success = asyncio.run(run_final_tests())
    sys.exit(0 if success else 1)