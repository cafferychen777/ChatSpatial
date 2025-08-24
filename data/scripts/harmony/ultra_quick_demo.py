#!/usr/bin/env python3
"""
ChatSpatial Harmony Ultra-Quick Demo (3-minute version)
Demonstrates batch integration with quick_mode enabled
"""

import sys
import time
import scanpy as sc
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Add ChatSpatial to path
sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

def run_ultra_quick_demo():
    """Run the ultra-quick ChatSpatial Harmony demo with quick_mode"""
    print("⚡ ChatSpatial Harmony Ultra-Quick Demo")
    print("=" * 45)
    print("Target: Complete integration in under 3 minutes with quick_mode")
    
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
        
        # Run integration with quick_mode
        print("\n⚡ Running Harmony integration with quick_mode...")
        integration_start = time.time()
        
        combined = integrate_multiple_samples(
            [adata1, adata2], 
            batch_key='batch',
            method='harmony',
            n_pcs=10,  # Further reduced for ultra speed
            quick_mode=True  # Enable quick mode!
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
        
        if total_time <= 180:  # 3 minutes
            print("🚀 SUCCESS: Completed in under 3 minutes!")
        elif total_time <= 300:  # 5 minutes
            print("✅ Good: Completed in under 5 minutes")
        else:
            print("⚠️  Took longer than expected, consider optimization")
        
        # Quick quality check
        if 'X_harmony' in combined.obsm and 'X_umap' in combined.obsm:
            print("✅ Generated embeddings: X_harmony, X_umap")
            
            # Check integration quality quickly
            from sklearn.metrics import silhouette_score
            try:
                # Silhouette score for batch mixing (lower is better for batch correction)
                batch_silhouette = silhouette_score(combined.obsm['X_harmony'][:, :5], combined.obs['batch'])
                # Silhouette score for cell type separation (higher is better)
                celltype_silhouette = silhouette_score(combined.obsm['X_harmony'][:, :5], combined.obs['cell_type'])
                
                print(f"📈 Quick Quality Metrics:")
                print(f"  Batch mixing score (lower=better): {batch_silhouette:.3f}")
                print(f"  Cell type separation (higher=better): {celltype_silhouette:.3f}")
                
                if batch_silhouette < 0.3 and celltype_silhouette > 0.3:
                    print("🎯 Good integration quality!")
                else:
                    print("📊 Integration completed (quality varies)")
                    
            except Exception as e:
                print(f"⚠️  Could not compute quality metrics: {e}")
        else:
            print("❌ Missing expected embeddings")
            return None
            
        return combined
        
    except Exception as e:
        print(f"❌ Demo failed: {e}")
        import traceback
        traceback.print_exc()
        return None

def compare_modes():
    """Compare normal vs quick mode performance"""
    print("\n🏁 Comparing Normal vs Quick Mode")
    print("=" * 40)
    
    from chatspatial.tools.integration import integrate_multiple_samples
    
    # Load data once
    adata1 = sc.read_h5ad("quick_demo_293t.h5ad")
    adata2 = sc.read_h5ad("quick_demo_jurkat.h5ad")
    
    results = {}
    
    for mode_name, quick_mode in [("Normal", False), ("Quick", True)]:
        print(f"\n🧪 Testing {mode_name} Mode...")
        start_time = time.time()
        
        try:
            combined = integrate_multiple_samples(
                [adata1.copy(), adata2.copy()], 
                batch_key='batch',
                method='harmony',
                n_pcs=10,
                quick_mode=quick_mode
            )
            
            elapsed = time.time() - start_time
            results[mode_name] = {
                'time': elapsed,
                'success': True,
                'shape': combined.shape
            }
            print(f"  ✅ {mode_name} mode: {elapsed:.1f}s")
            
        except Exception as e:
            results[mode_name] = {
                'time': float('inf'),
                'success': False,
                'error': str(e)
            }
            print(f"  ❌ {mode_name} mode failed: {e}")
    
    # Summary
    if results['Normal']['success'] and results['Quick']['success']:
        speedup = results['Normal']['time'] / results['Quick']['time']
        print(f"\n📊 Performance Comparison:")
        print(f"  Normal mode: {results['Normal']['time']:.1f}s")
        print(f"  Quick mode:  {results['Quick']['time']:.1f}s")
        print(f"  Speedup:     {speedup:.1f}x faster")
        
        if speedup > 1.5:
            print("🚀 Quick mode provides significant speedup!")
        else:
            print("📈 Quick mode provides moderate improvement")
    
    return results

if __name__ == "__main__":
    # Run ultra-quick demo
    combined = run_ultra_quick_demo()
    
    if combined is not None:
        print("\n" + "="*50)
        # Optionally run comparison (comment out to save time)
        # compare_modes()
        
        print(f"\n🎉 Ultra-Quick Demo Complete!")
        print(f"💡 For Claude Desktop usage:")
        print(f"   integrate_multiple_samples([adata1, adata2], quick_mode=True)")
        print(f"   This reduces Harmony iterations from 10 to 5 and increases sigma for faster convergence.")