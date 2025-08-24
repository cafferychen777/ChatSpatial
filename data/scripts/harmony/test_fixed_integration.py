#!/usr/bin/env python3
"""
Test the fixed ChatSpatial integration API
"""

import sys
import time
import scanpy as sc
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Add ChatSpatial to path
sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

def test_chatspatial_integration():
    """Test the fixed ChatSpatial integration"""
    print("🧪 Testing Fixed ChatSpatial Integration")
    print("=" * 50)
    
    from chatspatial.tools.integration import integrate_multiple_samples
    
    # Load datasets
    print("📂 Loading datasets...")
    adata1 = sc.read_h5ad("quick_demo_293t.h5ad")
    adata2 = sc.read_h5ad("quick_demo_jurkat.h5ad")
    
    print(f"Dataset 1: {adata1.shape}")
    print(f"Dataset 2: {adata2.shape}")
    print(f"Batch info 1: {adata1.obs['batch'].unique()}")
    print(f"Batch info 2: {adata2.obs['batch'].unique()}")
    
    # Test integration
    print("\n🔄 Running ChatSpatial Harmony integration...")
    start_time = time.time()
    
    try:
        combined = integrate_multiple_samples(
            [adata1, adata2],
            batch_key='batch',
            method='harmony',
            n_pcs=20
        )
        
        integration_time = time.time() - start_time
        print(f"\n✅ Integration successful!")
        print(f"📊 Results:")
        print(f"  Combined shape: {combined.shape}")
        print(f"  Integration time: {integration_time:.1f}s")
        
        # Check batch distribution
        if 'batch' in combined.obs:
            batch_counts = combined.obs['batch'].value_counts()
            print(f"  Batch distribution: {dict(batch_counts)}")
        
        # Check cell types
        if 'cell_type' in combined.obs:
            celltype_counts = combined.obs['cell_type'].value_counts()
            print(f"  Cell type distribution: {dict(celltype_counts)}")
        
        # Check embeddings
        embeddings = []
        if 'X_harmony' in combined.obsm:
            embeddings.append(f"X_harmony ({combined.obsm['X_harmony'].shape})")
        if 'X_umap' in combined.obsm:
            embeddings.append(f"X_umap ({combined.obsm['X_umap'].shape})")
        if 'X_pca' in combined.obsm:
            embeddings.append(f"X_pca ({combined.obsm['X_pca'].shape})")
        
        print(f"  Generated embeddings: {', '.join(embeddings)}")
        
        # Quick quality check
        from sklearn.metrics import silhouette_score
        try:
            if 'X_harmony' in combined.obsm:
                batch_silhouette = silhouette_score(
                    combined.obsm['X_harmony'][:, :5], 
                    combined.obs['batch']
                )
                print(f"  Batch mixing score: {batch_silhouette:.3f} (lower is better)")
        except Exception as e:
            print(f"  Could not compute quality metrics: {e}")
        
        return combined
        
    except Exception as e:
        print(f"❌ Integration failed: {e}")
        import traceback
        traceback.print_exc()
        return None

def test_with_larger_datasets():
    """Test with larger datasets"""
    print("\n🧪 Testing with Larger Datasets")
    print("=" * 50)
    
    from chatspatial.tools.integration import integrate_multiple_samples
    
    try:
        # Load larger datasets
        adata1 = sc.read_h5ad("pure_293t_simulated.h5ad")
        adata2 = sc.read_h5ad("pure_jurkat_simulated.h5ad")
        
        print(f"Large dataset 1: {adata1.shape}")
        print(f"Large dataset 2: {adata2.shape}")
        
        start_time = time.time()
        
        combined = integrate_multiple_samples(
            [adata1, adata2],
            batch_key='batch',
            method='harmony',
            n_pcs=30
        )
        
        integration_time = time.time() - start_time
        print(f"\n✅ Large dataset integration successful!")
        print(f"📊 Results:")
        print(f"  Combined shape: {combined.shape}")
        print(f"  Integration time: {integration_time:.1f}s")
        
        return combined
        
    except Exception as e:
        print(f"❌ Large dataset integration failed: {e}")
        return None

if __name__ == "__main__":
    # Test quick demo datasets
    result1 = test_chatspatial_integration()
    
    # Test larger datasets if quick one works
    if result1 is not None:
        print("\n" + "="*60)
        result2 = test_with_larger_datasets()
        
        if result2 is not None:
            print(f"\n🎉 All tests passed!")
            print(f"💡 ChatSpatial Harmony integration is working correctly!")
        else:
            print(f"\n⚠️  Quick demo works, but large datasets had issues")
    else:
        print(f"\n❌ Basic integration test failed")