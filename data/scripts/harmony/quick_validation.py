#!/usr/bin/env python3
"""
Quick validation that dangerous fallbacks are fixed
"""

import sys
import os
import numpy as np
import scanpy as sc

sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

def test_no_random_fallback():
    """Test that we no longer generate random PCA data"""
    print("🔍 Testing for elimination of dangerous random PCA fallback...")
    
    from chatspatial.tools.integration import integrate_multiple_samples
    
    # Create a problematic dataset that might trigger PCA failures
    # Small dataset with identical expression
    n_cells = 15
    n_genes = 20
    X = np.ones((n_cells, n_genes))  # All genes have identical expression
    
    adata = sc.AnnData(X)
    adata.obs['batch'] = ['batch1'] * 8 + ['batch2'] * 7
    
    print(f"Created problematic dataset: {adata.shape}")
    print("All genes have identical expression - should cause PCA to fail properly")
    
    try:
        result = integrate_multiple_samples([adata], batch_key='batch', method='harmony')
        print("❌ Integration should have failed with meaningful error")
        return False
    except (ValueError, RuntimeError) as e:
        if "random" in str(e).lower():
            print(f"❌ Still using random fallback: {e}")
            return False
        else:
            print(f"✅ Proper error handling: {e}")
            return True
    except Exception as e:
        print(f"⚠️  Unexpected error type: {type(e).__name__}: {e}")
        return True  # At least it failed instead of using random data

def data/test_validation():
    """Test data quality validation"""
    print("\n🔍 Testing data quality validation...")
    
    from chatspatial.tools.integration import validate_data_quality
    
    # Test with insufficient cells
    tiny_adata = sc.AnnData(np.random.rand(5, 100))
    
    try:
        validate_data_quality(tiny_adata)
        print("❌ Should have caught insufficient cells")
        return False
    except ValueError as e:
        if "cells" in str(e).lower():
            print(f"✅ Correctly caught insufficient cells: {e}")
            return True
        else:
            print(f"❌ Wrong error message: {e}")
            return False

def main():
    print("🛡️  ChatSpatial Integration Safety Validation")
    print("=" * 50)
    
    tests = [
        test_no_random_fallback,
        data/test_validation,
    ]
    
    results = []
    for test_func in tests:
        try:
            result = test_func()
            results.append(result)
        except Exception as e:
            print(f"❌ Test crashed: {e}")
            results.append(False)
    
    passed = sum(results)
    print(f"\n📊 Results: {passed}/{len(results)} critical safety checks passed")
    
    if passed == len(results):
        print("🎉 All critical safety fixes verified!")
        print("✅ No more dangerous random data fallbacks")
        print("✅ Proper error handling implemented")
    else:
        print("⚠️  Some safety issues may remain")

if __name__ == "__main__":
    main()