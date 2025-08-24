#!/usr/bin/env python3
"""
Quick Harmony Integration Tool for Demonstrations
Separate from production integration.py, optimized for speed and demos
"""

import sys
import time
import scanpy as sc
import numpy as np
import logging
import pandas as pd

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Add ChatSpatial to path for using base functionality
sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

def quick_harmony_demo(adatas, batch_key='batch', max_iter=5, n_pcs=15):
    """
    Quick Harmony integration optimized for demonstrations
    
    NOT FOR PRODUCTION USE - This is a demo-optimized version
    """
    print("⚡ Quick Harmony Integration (Demo Version)")
    print("-" * 40)
    
    import scanpy as sc
    
    # Quick validation
    if isinstance(adatas, list):
        for i, adata in enumerate(adatas):
            if batch_key not in adata.obs:
                adata.obs[batch_key] = f"batch_{i}"
        combined = adatas[0].concatenate(adatas[1:], batch_key=batch_key, join='outer')
    else:
        combined = adatas
    
    print(f"📊 Dataset: {combined.shape}")
    
    # Quick preprocessing
    print("🔄 Quick preprocessing...")
    start_time = time.time()
    
    # Fast normalization
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    
    # Quick HVG selection
    sc.pp.highly_variable_genes(combined, n_top_genes=min(1000, combined.n_vars))
    combined = combined[:, combined.var.highly_variable]
    print(f"  Filtered to {combined.n_vars} HVG genes")
    
    # Fast scaling
    sc.pp.scale(combined, max_value=10)
    
    # Quick PCA
    n_pcs_actual = min(n_pcs, combined.n_vars, combined.n_obs - 1)
    sc.tl.pca(combined, n_comps=n_pcs_actual, svd_solver='randomized')
    
    prep_time = time.time() - start_time
    print(f"  Preprocessing: {prep_time:.1f}s")
    
    # Quick Harmony
    print("🎯 Quick Harmony integration...")
    harmony_start = time.time()
    
    try:
        import harmonypy as hm
        
        meta_data = pd.DataFrame({batch_key: combined.obs[batch_key]})
        
        # Optimized for speed
        harmony_out = hm.run_harmony(
            data_mat=combined.obsm['X_pca'],
            meta_data=meta_data,
            vars_use=[batch_key],
            sigma=0.2,  # Larger sigma for faster convergence
            max_iter_harmony=max_iter,  # Reduced iterations
            verbose=True
        )
        
        combined.obsm['X_harmony'] = harmony_out.Z_corr.T
        
        harmony_time = time.time() - harmony_start
        print(f"  Harmony: {harmony_time:.1f}s")
        
        # Quick UMAP
        print("📈 Quick UMAP...")
        umap_start = time.time()
        
        sc.pp.neighbors(combined, use_rep='X_harmony', n_neighbors=15)
        sc.tl.umap(combined, min_dist=0.3, spread=1.0)
        
        umap_time = time.time() - umap_start
        print(f"  UMAP: {umap_time:.1f}s")
        
        total_time = time.time() - start_time
        print(f"✅ Total time: {total_time:.1f}s")
        
        return combined
        
    except ImportError:
        raise ImportError("harmonypy required: pip install harmonypy")
    except Exception as e:
        raise RuntimeError(f"Quick Harmony failed: {e}")

def run_claude_desktop_demo():
    """Demo optimized for Claude Desktop usage examples"""
    print("🖥️  Claude Desktop Harmony Demo")
    print("=" * 50)
    
    # Load demo datasets
    try:
        adata1 = sc.read_h5ad("quick_demo_293t.h5ad")
        adata2 = sc.read_h5ad("quick_demo_jurkat.h5ad")
    except FileNotFoundError:
        print("❌ Quick demo datasets not found. Run create_quick_demo.py first.")
        return None
    
    print("📋 Example for Claude Desktop:")
    print("```python")
    print("# Load your datasets")
    print("adata1 = sc.read_h5ad('dataset1.h5ad')")
    print("adata2 = sc.read_h5ad('dataset2.h5ad')")
    print("")
    print("# Use ChatSpatial integration (production version)")
    print("from chatspatial.tools.integration import integrate_multiple_samples")
    print("combined = integrate_multiple_samples([adata1, adata2], method='harmony')")
    print("```")
    print("")
    
    # Run actual demo
    start_time = time.time()
    combined = quick_harmony_demo([adata1, adata2], batch_key='batch')
    total_time = time.time() - start_time
    
    print(f"\n📊 Results Summary:")
    print(f"  Final shape: {combined.shape}")
    print(f"  Batches: {dict(combined.obs['batch'].value_counts())}")
    print(f"  Cell types: {dict(combined.obs['cell_type'].value_counts())}")
    print(f"  Total time: {total_time:.1f}s")
    
    if total_time <= 120:  # 2 minutes
        print("🚀 Excellent: Under 2 minutes!")
    elif total_time <= 300:  # 5 minutes
        print("✅ Good: Under 5 minutes")
    else:
        print("⚠️  Consider data size reduction")
    
    return combined

if __name__ == "__main__":
    # Run the Claude Desktop optimized demo
    combined = run_claude_desktop_demo()
    
    if combined is not None:
        print(f"\n💡 For your Claude Desktop sessions:")
        print(f"Use the production integrate_multiple_samples() function.")
        print(f"This quick demo shows what's possible in {time.time():.0f} seconds!")