#!/usr/bin/env python3
"""
Test the clean integration API that expects preprocessed data
"""

import sys
import time
import scanpy as sc
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Add ChatSpatial to path
sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

def test_with_preprocessed_data():
    """Test integration with properly preprocessed data"""
    print("🧪 Testing Clean Integration with Preprocessed Data")
    print("=" * 55)
    
    from chatspatial.tools.integration import integrate_multiple_samples
    
    # Load raw data
    print("📂 Loading raw datasets...")
    adata1 = sc.read_h5ad("quick_demo_293t.h5ad")
    adata2 = sc.read_h5ad("quick_demo_jurkat.h5ad")
    
    print(f"Dataset 1: {adata1.shape}")
    print(f"Dataset 2: {adata2.shape}")
    
    # Proper preprocessing (should be done by LLM/preprocessing.py)
    print("\n🔄 Preprocessing datasets (normally done by LLM/preprocessing.py)...")
    
    for i, adata in enumerate([adata1, adata2], 1):
        print(f"  Processing dataset {i}...")
        
        # Normalize and log transform
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Find highly variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=min(1000, adata.n_vars))
        
        # Scale data
        adata.raw = adata  # Save raw data
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata, max_value=10)
        
        # Update the reference
        if i == 1:
            adata1 = adata
        else:
            adata2 = adata
    
    print(f"  Preprocessed dataset 1: {adata1.shape}")
    print(f"  Preprocessed dataset 2: {adata2.shape}")
    
    # Test clean integration
    print("\n⚡ Running clean integration (expects preprocessed data)...")
    start_time = time.time()
    
    try:
        combined = integrate_multiple_samples(
            [adata1, adata2],
            batch_key='batch',
            method='harmony',
            n_pcs=15
        )
        
        integration_time = time.time() - start_time
        print(f"\n✅ Clean integration successful!")
        print(f"📊 Results:")
        print(f"  Combined shape: {combined.shape}")
        print(f"  Integration time: {integration_time:.1f}s")
        
        # Verify results
        if 'X_harmony' in combined.obsm:
            print(f"  ✅ Harmony embedding: {combined.obsm['X_harmony'].shape}")
        if 'X_umap' in combined.obsm:
            print(f"  ✅ UMAP embedding: {combined.obsm['X_umap'].shape}")
        
        return combined
        
    except Exception as e:
        print(f"❌ Clean integration failed: {e}")
        import traceback
        traceback.print_exc()
        return None

def test_with_unprocessed_data():
    """Test that integration properly rejects unprocessed data"""
    print("\n🧪 Testing Rejection of Unprocessed Data")
    print("=" * 45)
    
    from chatspatial.tools.integration import integrate_multiple_samples
    
    # Load raw data without preprocessing
    print("📂 Loading raw (unprocessed) datasets...")
    adata1 = sc.read_h5ad("quick_demo_293t.h5ad")
    adata2 = sc.read_h5ad("quick_demo_jurkat.h5ad")
    
    # Don't preprocess - test raw data
    print("⚠️  Attempting integration with raw data (should fail gracefully)...")
    
    try:
        combined = integrate_multiple_samples(
            [adata1, adata2],
            batch_key='batch',
            method='harmony'
        )
        print("❌ Integration should have failed with unprocessed data")
        return False
        
    except ValueError as e:
        if "negative values" in str(e) or "normalize" in str(e).lower():
            print(f"✅ Correctly rejected unprocessed data: {e}")
            return True
        else:
            print(f"❓ Unexpected error: {e}")
            return False
    except Exception as e:
        print(f"❓ Unexpected error type: {type(e).__name__}: {e}")
        return False

def test_integration_workflow():
    """Test the complete workflow: preprocessing -> integration"""
    print("\n🧪 Testing Complete Workflow")
    print("=" * 35)
    
    # This simulates what should happen in Claude Desktop:
    # 1. User loads data
    # 2. LLM calls preprocess_data 
    # 3. LLM calls integrate_multiple_samples with preprocessed data
    
    from chatspatial.tools.integration import integrate_multiple_samples
    
    print("1️⃣ Load raw data")
    adata1 = sc.read_h5ad("quick_demo_293t.h5ad")
    adata2 = sc.read_h5ad("quick_demo_jurkat.h5ad")
    
    print("2️⃣ Preprocess (done by LLM/preprocessing.py)")
    datasets = [adata1, adata2]
    
    for adata in datasets:
        # Standard preprocessing pipeline
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=500)
        adata.raw = adata
        adata = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata)
    
    print("3️⃣ Integration (clean, focused function)")
    start_time = time.time()
    
    combined = integrate_multiple_samples(
        datasets,
        batch_key='batch', 
        method='harmony',
        n_pcs=10
    )
    
    workflow_time = time.time() - start_time
    print(f"✅ Complete workflow: {workflow_time:.1f}s")
    print(f"Final result: {combined.shape}")
    
    return combined

if __name__ == "__main__":
    print("🏗️  Testing Clean Integration Architecture")
    print("=" * 60)
    
    # Test 1: Proper preprocessing workflow
    result1 = test_with_preprocessed_data()
    
    # Test 2: Rejection of unprocessed data  
    result2 = test_with_unprocessed_data()
    
    # Test 3: Complete workflow
    if result1 is not None and result2:
        result3 = test_integration_workflow()
        
        if result3 is not None:
            print(f"\n🎉 All tests passed!")
            print(f"✅ Clean architecture: integration.py focuses on integration only")
            print(f"✅ Proper separation: preprocessing handled separately") 
            print(f"✅ Good error handling: rejects unprocessed data")
        else:
            print(f"\n⚠️  Basic tests passed but workflow test failed")
    else:
        print(f"\n❌ Some basic tests failed")