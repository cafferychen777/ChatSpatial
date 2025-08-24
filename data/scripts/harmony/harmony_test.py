#!/usr/bin/env python3
"""
Minimal Harmony functionality test
"""

import scanpy as sc
import numpy as np
import pandas as pd

# Load data
print("Loading data...")
adata = sc.read_h5ad("jurkat_293t_mixture_simulated.h5ad")
print(f"Data shape: {adata.shape}")

# Basic preprocessing 
print("Preprocessing...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Keep only top 1000 HVGs for speed
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
adata = adata[:, :1000]  # Further reduce for speed

print(f"Filtered data shape: {adata.shape}")

# PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', n_comps=20)  # Reduce components

print("Testing Harmony integration...")
try:
    import harmonypy as hm
    
    # Run Harmony with minimal iterations for testing
    ho = hm.run_harmony(
        adata.obsm['X_pca'], 
        adata.obs, 
        vars_use=['batch'],
        max_iter_harmony=5  # Reduce iterations for speed
    )
    adata.obsm['X_harmony'] = ho.Z_corr.T
    
    print("✅ Harmony integration successful!")
    print(f"Harmony embedding shape: {adata.obsm['X_harmony'].shape}")
    
    # Quick UMAP
    sc.pp.neighbors(adata, use_rep='X_harmony', n_neighbors=10)
    sc.tl.umap(adata)
    
    print("✅ Full pipeline completed successfully!")
    print("📊 Final UMAP embedding shape:", adata.obsm['X_umap'].shape)
    
except Exception as e:
    print(f"❌ Error in Harmony integration: {e}")
    import traceback
    traceback.print_exc()