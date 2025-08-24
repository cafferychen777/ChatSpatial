#!/usr/bin/env python3
"""
Verbose test to see exactly what happens in harmony tutorial
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

print("=== HARMONY TUTORIAL VERBOSE TEST ===")
print(f"Python version: {sys.version}")

# Load the simulated dataset
print("1. Loading dataset...")
adata = sc.read_h5ad("jurkat_293t_mixture_simulated.h5ad")

print(f"Dataset shape: {adata.shape}")
print(f"Cell types: {adata.obs['cell_type'].value_counts()}")
print(f"Batches: {adata.obs['batch'].value_counts()}")

# Basic preprocessing
print("2. Basic preprocessing...")
print(f"  Original cells: {adata.n_obs}")
sc.pp.filter_cells(adata, min_genes=200)
print(f"  After cell filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=3)
print(f"  After gene filter: {adata.n_vars}")

# Calculate QC metrics
print("3. QC metrics...")
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Normalization and log transformation
print("4. Normalization...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find highly variable genes
print("5. Finding HVGs...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
print(f"  HVGs selected: {adata.n_vars}")

# Scale data and PCA
print("6. Scaling and PCA...")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
print(f"  PCA shape: {adata.obsm['X_pca'].shape}")

# UMAP before integration
print("7. UMAP before integration...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
print(f"  UMAP shape: {adata.obsm['X_umap'].shape}")

# Plot before integration
print("8. Plotting before integration...")
sc.pl.umap(adata, color=['batch', 'cell_type'], save='_before_harmony_verbose.pdf')
print("  ✅ Before-integration plot saved")

# Harmony integration
print("9. Testing Harmony import...")
try:
    import harmonypy as hm
    print("  ✅ harmonypy imported successfully")
    
    print("10. Running Harmony...")
    ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, vars_use=['batch'])
    adata.obsm['X_harmony'] = ho.Z_corr.T
    print(f"  ✅ Harmony completed, shape: {adata.obsm['X_harmony'].shape}")
    
    # Use Harmony embeddings for downstream analysis
    print("11. Post-harmony analysis...")
    sc.pp.neighbors(adata, use_rep='X_harmony')
    sc.tl.umap(adata)
    print(f"  ✅ Post-harmony UMAP shape: {adata.obsm['X_umap'].shape}")
    
    # Plot after integration
    print("12. Plotting after integration...")
    sc.pl.umap(adata, color=['batch', 'cell_type'], save='_after_harmony_verbose.pdf')
    print("  ✅ After-integration plot saved")
    
    print("13. ✅ Harmony integration completed!")
    
except ImportError:
    print("  ❌ harmonypy not available")
    print("  Using scanpy's harmony as fallback...")
    sc.external.pp.harmony_integrate(adata, key='batch')
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['batch', 'cell_type'], save='_scanpy_integration_verbose.pdf')

# Clustering
print("14. Clustering...")
sc.tl.leiden(adata)
print(f"  Leiden clusters: {len(adata.obs['leiden'].unique())}")
sc.pl.umap(adata, color='leiden', save='_clusters_verbose.pdf')
print("  ✅ Clustering plot saved")

# Save integrated data
print("15. Saving results...")
adata.write_h5ad("jurkat_293t_integrated_verbose.h5ad")
print("💾 ✅ Saved integrated dataset")

print("\n=== SUMMARY ===")
print(f"Final data shape: {adata.shape}")
print(f"Available obsm keys: {list(adata.obsm.keys())}")
print("🎉 ALL TESTS COMPLETED SUCCESSFULLY!")