#!/usr/bin/env python3
"""
ChatSpatial Harmony Quick Demo (5-minute version)
Demonstrates batch integration with 1000 cells in under 5 minutes
"""

import sys
import time
import scanpy as sc
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Add ChatSpatial to path
sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

def run_quick_demo():
    """Run the quick ChatSpatial Harmony demo"""
    print("🚀 ChatSpatial Harmony Quick Demo")
    print("=" * 40)
    print("Target: Complete integration in under 5 minutes")
    
    start_time = time.time()
    
    try:
        from chatspatial.tools.integration import integrate_multiple_samples
        
        # Load quick demo datasets
        print("\n📂 Loading quick demo datasets...")
        adata1 = sc.read_h5ad("quick_demo_293t.h5ad")
        adata2 = sc.read_h5ad("quick_demo_jurkat.h5ad")
        
        print(f"Dataset 1 (293T): {adata1.shape}")
        print(f"Dataset 2 (Jurkat): {adata2.shape}")
        
        load_time = time.time() - start_time
        print(f"⏱️  Loading time: {load_time:.1f}s")
        
        # Run integration
        print("\n🔄 Running Harmony integration...")
        integration_start = time.time()
        
        combined = integrate_multiple_samples(
            [adata1, adata2], 
            batch_key='batch',
            method='harmony',
            n_pcs=15  # Reduced for speed
        )
        
        integration_time = time.time() - integration_start
        total_time = time.time() - start_time
        
        print(f"\n✅ Integration completed!")
        print(f"📊 Results:")
        print(f"  Combined dataset: {combined.shape}")
        print(f"  Batch distribution: {dict(combined.obs['batch'].value_counts())}")
        print(f"  Cell type distribution: {dict(combined.obs['cell_type'].value_counts())}")
        
        print(f"\n⏱️  Performance:")
        print(f"  Integration time: {integration_time:.1f}s") 
        print(f"  Total time: {total_time:.1f}s")
        
        if total_time <= 300:  # 5 minutes
            print("🎉 SUCCESS: Completed in under 5 minutes!")
        else:
            print("⚠️  Took longer than 5 minutes, consider further optimization")
        
        # Quick quality check
        if 'X_harmony' in combined.obsm and 'X_umap' in combined.obsm:
            print("✅ Generated embeddings: X_harmony, X_umap")
        else:
            print("❌ Missing expected embeddings")
            
        return combined
        
    except Exception as e:
        print(f"❌ Demo failed: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    run_quick_demo()
