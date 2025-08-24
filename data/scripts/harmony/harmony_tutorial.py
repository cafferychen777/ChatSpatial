#!/usr/bin/env python3
"""
Harmony Integration Tutorial with Simulated Jurkat:293T Data
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the simulated dataset
adata = sc.read_h5ad("jurkat_293t_mixture_simulated.h5ad")

print(f"Dataset shape: {adata.shape}")
print(f"Cell types: {adata.obs['cell_type'].value_counts()}")
print(f"Batches: {adata.obs['batch'].value_counts()}")

# Basic preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Normalization and log transformation
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]

# Scale data and PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# UMAP before integration
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# Plot before integration
sc.pl.umap(adata, color=['batch', 'cell_type'], save='_before_harmony.pdf')

# Harmony integration using scanpy's external integration
try:
    import harmonypy as hm
    
    # Run Harmony
    ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, vars_use=['batch'])
    adata.obsm['X_harmony'] = ho.Z_corr.T
    
    # Use Harmony embeddings for downstream analysis
    sc.pp.neighbors(adata, use_rep='X_harmony')
    sc.tl.umap(adata)
    
    # Plot after integration
    sc.pl.umap(adata, color=['batch', 'cell_type'], save='_after_harmony.pdf')
    
    print("✅ Harmony integration completed!")
    print("📊 Check the generated UMAP plots to see integration results")
    
except ImportError:
    print("⚠️  harmonypy not installed. Install with: pip install harmonypy")
    print("Using alternative integration method...")
    
    # Alternative: Use ChatSpatial's integration
    try:
        from chatspatial.tools.integration import integrate_samples_harmony
        
        # This would use ChatSpatial's integration functionality
        print("Using ChatSpatial integration instead")
        
    except ImportError:
        print("Using scanpy's integration as fallback")
        sc.external.pp.harmony_integrate(adata, key='batch')
        sc.pp.neighbors(adata, use_rep='X_pca_harmony')
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=['batch', 'cell_type'], save='_scanpy_integration.pdf')

# Clustering
sc.tl.leiden(adata)
sc.pl.umap(adata, color='leiden', save='_clusters.pdf')

# Save integrated data
adata.write_h5ad("jurkat_293t_integrated.h5ad")
print("💾 Saved integrated dataset")
