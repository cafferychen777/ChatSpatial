#!/usr/bin/env python3
"""
Create quick demo datasets for ChatSpatial Harmony integration testing
Target: 1000 cells total, ~500-1000 genes, 2-5 minute integration time
"""

import scanpy as sc
import numpy as np
import pandas as pd
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def create_quick_demo_datasets():
    """Create subsampled datasets for quick demonstration"""
    print("🚀 Creating Quick Demo Datasets for ChatSpatial Harmony")
    print("=" * 60)
    
    # Load the two pure datasets (smallest and most distinct)
    print("📂 Loading source datasets...")
    adata_293t = sc.read_h5ad("pure_293t_simulated.h5ad")
    adata_jurkat = sc.read_h5ad("pure_jurkat_simulated.h5ad")
    
    print(f"293T dataset: {adata_293t.shape}")
    print(f"Jurkat dataset: {adata_jurkat.shape}")
    
    # Subsample cells for quick demo
    np.random.seed(42)  # Reproducible subsampling
    
    # Take 500 cells from each dataset
    n_cells_per_type = 500
    
    print(f"\n🎯 Subsampling to {n_cells_per_type} cells per dataset...")
    
    # Subsample 293T
    idx_293t = np.random.choice(adata_293t.n_obs, size=n_cells_per_type, replace=False)
    adata_293t_sub = adata_293t[idx_293t, :].copy()
    
    # Subsample Jurkat  
    idx_jurkat = np.random.choice(adata_jurkat.n_obs, size=n_cells_per_type, replace=False)
    adata_jurkat_sub = adata_jurkat[idx_jurkat, :].copy()
    
    print(f"✅ Subsampled 293T: {adata_293t_sub.shape}")
    print(f"✅ Subsampled Jurkat: {adata_jurkat_sub.shape}")
    
    # Find highly variable genes using the combined dataset for better selection
    print(f"\n🧬 Selecting highly variable genes...")
    
    # Combine datasets temporarily for HVG selection
    combined_temp = sc.concat([adata_293t_sub, adata_jurkat_sub], axis=0)
    
    # Basic preprocessing for HVG selection
    sc.pp.normalize_total(combined_temp, target_sum=1e4)
    sc.pp.log1p(combined_temp)
    sc.pp.highly_variable_genes(combined_temp, min_mean=0.01, max_mean=5, min_disp=0.5, n_top_genes=800)
    
    # Get HVG genes
    hvg_genes = combined_temp.var[combined_temp.var.highly_variable].index
    n_hvg = len(hvg_genes)
    print(f"✅ Selected {n_hvg} highly variable genes")
    
    # Apply gene filtering to original subsampled datasets
    adata_293t_quick = adata_293t_sub[:, hvg_genes].copy()
    adata_jurkat_quick = adata_jurkat_sub[:, hvg_genes].copy()
    
    # Update batch information to be more descriptive for demo
    adata_293t_quick.obs['batch'] = 'demo_293t'
    adata_293t_quick.obs['cell_type'] = '293T'
    adata_293t_quick.obs['dataset'] = 'quick_demo'
    
    adata_jurkat_quick.obs['batch'] = 'demo_jurkat'  
    adata_jurkat_quick.obs['cell_type'] = 'Jurkat'
    adata_jurkat_quick.obs['dataset'] = 'quick_demo'
    
    # Save quick demo datasets
    print(f"\n💾 Saving quick demo datasets...")
    
    adata_293t_quick.write_h5ad("quick_demo_293t.h5ad")
    adata_jurkat_quick.write_h5ad("quick_demo_jurkat.h5ad")
    
    # Also create a pre-combined version for single-dataset testing
    combined_quick = sc.concat([adata_293t_quick, adata_jurkat_quick], axis=0)
    combined_quick.write_h5ad("quick_demo_combined.h5ad")
    
    # Calculate expected sizes
    size_293t = adata_293t_quick.X.nbytes / (1024*1024)
    size_jurkat = adata_jurkat_quick.X.nbytes / (1024*1024)
    size_combined = combined_quick.X.nbytes / (1024*1024)
    
    print(f"📊 Quick Demo Dataset Summary:")
    print(f"  quick_demo_293t.h5ad: {adata_293t_quick.shape} - {size_293t:.1f} MB")
    print(f"  quick_demo_jurkat.h5ad: {adata_jurkat_quick.shape} - {size_jurkat:.1f} MB") 
    print(f"  quick_demo_combined.h5ad: {combined_quick.shape} - {size_combined:.1f} MB")
    
    print(f"\n🎯 Expected performance:")
    print(f"  Total cells: {combined_quick.n_obs}")
    print(f"  Total genes: {combined_quick.n_vars}")
    print(f"  Estimated integration time: 2-5 minutes")
    print(f"  Memory usage: ~{size_combined:.0f} MB")
    
    return adata_293t_quick, adata_jurkat_quick, combined_quick

def create_quick_demo_script():
    """Create a quick demo script"""
    demo_script = '''#!/usr/bin/env python3
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
        print("\\n📂 Loading quick demo datasets...")
        adata1 = sc.read_h5ad("quick_demo_293t.h5ad")
        adata2 = sc.read_h5ad("quick_demo_jurkat.h5ad")
        
        print(f"Dataset 1 (293T): {adata1.shape}")
        print(f"Dataset 2 (Jurkat): {adata2.shape}")
        
        load_time = time.time() - start_time
        print(f"⏱️  Loading time: {load_time:.1f}s")
        
        # Run integration
        print("\\n🔄 Running Harmony integration...")
        integration_start = time.time()
        
        combined = integrate_multiple_samples(
            [adata1, adata2], 
            batch_key='batch',
            method='harmony',
            n_pcs=15  # Reduced for speed
        )
        
        integration_time = time.time() - integration_start
        total_time = time.time() - start_time
        
        print(f"\\n✅ Integration completed!")
        print(f"📊 Results:")
        print(f"  Combined dataset: {combined.shape}")
        print(f"  Batch distribution: {dict(combined.obs['batch'].value_counts())}")
        print(f"  Cell type distribution: {dict(combined.obs['cell_type'].value_counts())}")
        
        print(f"\\n⏱️  Performance:")
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
'''
    
    with open("quick_demo_harmony.py", "w") as f:
        f.write(demo_script)
    
    print("✅ Created quick_demo_harmony.py")

def main():
    """Main function"""
    try:
        # Create datasets
        adata_293t, adata_jurkat, combined = create_quick_demo_datasets()
        
        # Create demo script
        create_quick_demo_script()
        
        print(f"\n🎉 Quick Demo Setup Complete!")
        print(f"📁 Files created:")
        print(f"  - quick_demo_293t.h5ad")
        print(f"  - quick_demo_jurkat.h5ad") 
        print(f"  - quick_demo_combined.h5ad")
        print(f"  - quick_demo_harmony.py")
        
        print(f"\n🚀 To run the quick demo:")
        print(f"  python quick_demo_harmony.py")
        
    except Exception as e:
        print(f"❌ Failed to create quick demo: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()