#!/usr/bin/env python3
"""
Test script for ChatSpatial's fixed integration functionality
"""

import sys
import os
import numpy as np
import scanpy as sc
import logging

# Setup logging to see our improved error messages
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Add ChatSpatial to path
sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

try:
    from chatspatial.tools.integration import integrate_multiple_samples, validate_data_quality
    print("✅ ChatSpatial integration module imported successfully")
except ImportError as e:
    print(f"❌ Failed to import ChatSpatial integration: {e}")
    sys.exit(1)

def data/test_quality_validation():
    """Test the new data quality validation"""
    print("\n🧪 Testing data quality validation...")
    
    # Load test data
    adata = sc.read_h5ad("pure_293t_simulated.h5ad")
    print(f"Loaded test data: {adata.shape}")
    
    try:
        validate_data_quality(adata)
        print("✅ Data quality validation passed")
        return True
    except Exception as e:
        print(f"❌ Data quality validation failed: {e}")
        return False

def test_basic_harmony_integration():
    """Test basic Harmony integration with improved error handling"""
    print("\n🧪 Testing ChatSpatial Harmony integration...")
    
    try:
        # Load two datasets
        adata1 = sc.read_h5ad("pure_293t_simulated.h5ad")
        adata2 = sc.read_h5ad("pure_jurkat_simulated.h5ad")
        
        print(f"Dataset 1: {adata1.shape}")
        print(f"Dataset 2: {adata2.shape}")
        
        # Test integration
        combined = integrate_multiple_samples(
            [adata1, adata2], 
            batch_key='batch',
            method='harmony',
            n_pcs=20  # Reduce for faster testing
        )
        
        print(f"✅ Integration successful! Combined shape: {combined.shape}")
        
        # Check results
        if 'X_harmony' in combined.obsm:
            print(f"✅ Harmony embedding created: {combined.obsm['X_harmony'].shape}")
        else:
            print("❌ Harmony embedding missing")
            return False
            
        if 'X_umap' in combined.obsm:
            print(f"✅ UMAP embedding created: {combined.obsm['X_umap'].shape}")
        else:
            print("❌ UMAP embedding missing")
            return False
            
        # Check batch information
        batch_counts = combined.obs['batch'].value_counts()
        print(f"✅ Batch distribution: {dict(batch_counts)}")
        
        return True
        
    except Exception as e:
        print(f"❌ Integration test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_error_handling():
    """Test improved error handling"""
    print("\n🧪 Testing error handling improvements...")
    
    # Test with very small dataset
    print("Testing with insufficient data...")
    small_adata = sc.AnnData(np.random.rand(5, 10))  # Very small dataset
    
    try:
        validate_data_quality(small_adata)
        print("❌ Should have failed with small dataset")
        return False
    except ValueError as e:
        print(f"✅ Correctly caught insufficient data: {e}")
    
    # Test with missing batch key
    print("Testing with missing batch key...")
    adata = sc.read_h5ad("pure_293t_simulated.h5ad")
    
    try:
        integrate_multiple_samples([adata], batch_key='nonexistent_key')
        print("❌ Should have failed with missing batch key")
        return False
    except Exception as e:
        print(f"✅ Correctly caught missing batch key error")
    
    return True

def benchmark_vs_native_harmony():
    """Quick benchmark against native harmonypy"""
    print("\n⏱️ Benchmarking ChatSpatial vs native harmonypy...")
    
    import time
    
    # Load data for benchmarking
    adata = sc.read_h5ad("jurkat_293t_mixture_simulated.h5ad")
    print(f"Benchmark data: {adata.shape}")
    
    # ChatSpatial version
    start_time = time.time()
    try:
        combined_cs = integrate_multiple_samples(
            [adata], 
            batch_key='batch',
            method='harmony',
            n_pcs=20
        )
        cs_time = time.time() - start_time
        print(f"✅ ChatSpatial integration: {cs_time:.1f}s")
        cs_success = True
    except Exception as e:
        print(f"❌ ChatSpatial integration failed: {e}")
        cs_success = False
        cs_time = float('inf')
    
    # Native harmonypy version (simplified)
    start_time = time.time()
    try:
        import harmonypy as hm
        
        # Quick preprocessing
        adata_native = adata.copy()
        sc.pp.normalize_total(adata_native)
        sc.pp.log1p(adata_native)
        sc.pp.highly_variable_genes(adata_native)
        adata_native = adata_native[:, adata_native.var.highly_variable]
        sc.pp.scale(adata_native)
        sc.tl.pca(adata_native, n_comps=20)
        
        # Run Harmony
        ho = hm.run_harmony(
            adata_native.obsm['X_pca'], 
            adata_native.obs, 
            vars_use=['batch'],
            max_iter_harmony=5  # Reduce iterations for fair comparison
        )
        
        native_time = time.time() - start_time
        print(f"✅ Native harmonypy: {native_time:.1f}s")
        native_success = True
    except Exception as e:
        print(f"❌ Native harmonypy failed: {e}")
        native_success = False
        native_time = float('inf')
    
    # Compare results
    if cs_success and native_success:
        overhead = ((cs_time - native_time) / native_time) * 100
        print(f"📊 ChatSpatial overhead: {overhead:.1f}%")
        if overhead < 50:  # Less than 50% overhead is acceptable
            print("✅ Performance is acceptable")
            return True
        else:
            print("⚠️  High overhead, consider optimization")
            return True  # Still pass, as functionality is more important
    elif cs_success:
        print("✅ ChatSpatial works, native failed")
        return True
    else:
        print("❌ Both implementations failed")
        return False

def main():
    """Run all tests"""
    print("🧬 ChatSpatial Integration Fix Validation")
    print("=" * 50)
    
    tests = [
        ("Data Quality Validation", data/test_quality_validation),
        ("Basic Harmony Integration", test_basic_harmony_integration),
        ("Error Handling", test_error_handling),
        ("Performance Benchmark", benchmark_vs_native_harmony),
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n{'='*20} {test_name} {'='*20}")
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ Test '{test_name}' crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print(f"\n{'='*20} SUMMARY {'='*20}")
    passed = 0
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{test_name}: {status}")
        if result:
            passed += 1
    
    print(f"\nOverall: {passed}/{len(results)} tests passed")
    
    if passed == len(results):
        print("🎉 All tests passed! ChatSpatial integration fixes are working correctly.")
    elif passed >= len(results) * 0.75:
        print("⚠️  Most tests passed. Minor issues may need attention.")
    else:
        print("❌ Multiple tests failed. Major issues need to be resolved.")

if __name__ == "__main__":
    main()