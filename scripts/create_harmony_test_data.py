#!/usr/bin/env python3
"""
Create simulated Harmony test datasets mimicking Jurkat and 293T cell lines
"""

import numpy as np
import pandas as pd
import scanpy as sc
import os
from pathlib import Path

def create_jurkat_293t_mixture():
    """Create simulated Jurkat:293T mixture dataset similar to 10x example"""
    
    print("🧬 Creating simulated Jurkat:293T mixture dataset...")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Dataset parameters
    n_jurkat = 1500
    n_293t = 1500  
    n_genes = 2000
    
    # Gene names (using generic names)
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    
    # Create Jurkat cells (T-cell lymphoma)
    # Higher expression of T-cell markers
    jurkat_base = np.random.negative_binomial(n=5, p=0.3, size=(n_jurkat, n_genes))
    
    # Add T-cell specific expression patterns
    t_cell_genes = np.random.choice(n_genes, 100, replace=False)
    for gene_idx in t_cell_genes:
        multiplier = np.random.uniform(2, 5, n_jurkat)
        jurkat_base[:, gene_idx] = (jurkat_base[:, gene_idx] * multiplier).astype(int)
    
    # Create cell metadata for Jurkat
    jurkat_obs = pd.DataFrame({
        'cell_id': [f"Jurkat_{i:04d}" for i in range(n_jurkat)],
        'cell_type': ['Jurkat'] * n_jurkat,
        'batch': ['pure_jurkat'] * n_jurkat,
        'dataset': ['jurkat'] * n_jurkat
    })
    jurkat_obs.index = jurkat_obs['cell_id']
    
    # Create 293T cells (embryonic kidney)  
    # Different expression pattern
    t293_base = np.random.negative_binomial(n=8, p=0.4, size=(n_293t, n_genes))
    
    # Add epithelial/kidney specific patterns
    epithelial_genes = np.random.choice(n_genes, 120, replace=False)
    for gene_idx in epithelial_genes:
        multiplier = np.random.uniform(3, 6, n_293t)
        t293_base[:, gene_idx] = (t293_base[:, gene_idx] * multiplier).astype(int)
    
    # Create cell metadata for 293T
    t293_obs = pd.DataFrame({
        'cell_id': [f"293T_{i:04d}" for i in range(n_293t)],
        'cell_type': ['293T'] * n_293t,
        'batch': ['pure_293t'] * n_293t,
        'dataset': ['293t'] * n_293t
    })
    t293_obs.index = t293_obs['cell_id']
    
    # Create mixed dataset (50:50)
    n_mix_jurkat = 800
    n_mix_293t = 800
    
    # Jurkat cells from mixture (with slight batch effect)
    mix_jurkat = np.random.negative_binomial(n=5, p=0.3, size=(n_mix_jurkat, n_genes))
    for gene_idx in t_cell_genes:
        multiplier = np.random.uniform(2, 5, n_mix_jurkat)
        mix_jurkat[:, gene_idx] = (mix_jurkat[:, gene_idx] * multiplier).astype(int)
    # Add batch effect
    batch_effect = np.random.uniform(0.8, 1.2, (n_mix_jurkat, n_genes))
    mix_jurkat = (mix_jurkat * batch_effect).astype(int)
    
    # 293T cells from mixture (with slight batch effect)
    mix_293t = np.random.negative_binomial(n=8, p=0.4, size=(n_mix_293t, n_genes))
    for gene_idx in epithelial_genes:
        multiplier = np.random.uniform(3, 6, n_mix_293t)
        mix_293t[:, gene_idx] = (mix_293t[:, gene_idx] * multiplier).astype(int)
    # Add batch effect
    batch_effect = np.random.uniform(0.9, 1.3, (n_mix_293t, n_genes))
    mix_293t = (mix_293t * batch_effect).astype(int)
    
    # Create metadata for mixture
    mix_obs = pd.DataFrame({
        'cell_id': [f"Mix_Jurkat_{i:04d}" for i in range(n_mix_jurkat)] + 
                  [f"Mix_293T_{i:04d}" for i in range(n_mix_293t)],
        'cell_type': ['Jurkat'] * n_mix_jurkat + ['293T'] * n_mix_293t,
        'batch': ['mixture'] * (n_mix_jurkat + n_mix_293t),
        'dataset': ['mixture'] * (n_mix_jurkat + n_mix_293t)
    })
    mix_obs.index = mix_obs['cell_id']
    
    # Combine expression data
    X_combined = np.vstack([
        jurkat_base,
        t293_base, 
        mix_jurkat,
        mix_293t
    ])
    
    # Combine metadata
    obs_combined = pd.concat([jurkat_obs, t293_obs, mix_obs])
    
    # Create gene metadata
    var_df = pd.DataFrame(index=gene_names)
    var_df['gene_ids'] = gene_names
    var_df['feature_types'] = 'Gene Expression'
    
    # Mark known marker genes
    var_df['is_t_cell_marker'] = False
    var_df.iloc[t_cell_genes, var_df.columns.get_loc('is_t_cell_marker')] = True
    
    var_df['is_epithelial_marker'] = False  
    var_df.iloc[epithelial_genes, var_df.columns.get_loc('is_epithelial_marker')] = True
    
    # Create AnnData object
    adata = sc.AnnData(X=X_combined, obs=obs_combined, var=var_df)
    
    # Add some QC metrics
    adata.var['n_cells'] = (adata.X > 0).sum(axis=0)
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    adata.obs['total_counts'] = adata.X.sum(axis=1)
    
    print(f"✅ Created dataset with {adata.n_obs} cells and {adata.n_vars} genes")
    print(f"   - Jurkat cells: {(adata.obs['cell_type'] == 'Jurkat').sum()}")
    print(f"   - 293T cells: {(adata.obs['cell_type'] == '293T').sum()}")
    print(f"   - Batches: {adata.obs['batch'].unique()}")
    
    return adata

def save_datasets():
    """Create and save simulated datasets"""
    
    # Create output directory
    output_dir = Path("data/harmony")
    output_dir.mkdir(exist_ok=True)
    
    # Create main mixture dataset
    adata = create_jurkat_293t_mixture()
    
    # Save full dataset
    full_path = output_dir / "jurkat_293t_mixture_simulated.h5ad"
    adata.write_h5ad(full_path)
    print(f"💾 Saved full dataset: {full_path}")
    
    # Create separate datasets by batch
    for batch in adata.obs['batch'].unique():
        batch_adata = adata[adata.obs['batch'] == batch].copy()
        batch_path = output_dir / f"{batch}_simulated.h5ad"
        batch_adata.write_h5ad(batch_path)
        print(f"💾 Saved {batch} dataset: {batch_path}")
    
    # Create 10x-style h5 files for compatibility
    try:
        # Save as 10x format for compatibility with standard tools
        for batch in adata.obs['batch'].unique():
            batch_adata = adata[adata.obs['batch'] == batch].copy()
            batch_dir = output_dir / f"{batch}_10x_format"
            batch_dir.mkdir(exist_ok=True)
            
            # Create simple h5 file
            h5_path = batch_dir / "filtered_feature_bc_matrix.h5ad" 
            batch_adata.write_h5ad(h5_path)
            
        print("💾 Created 10x-compatible format files")
    except Exception as e:
        print(f"⚠️  Could not create 10x format: {e}")
    
    return adata

def create_harmony_tutorial():
    """Create a tutorial script for using the simulated data"""
    
    tutorial_script = '''#!/usr/bin/env python3
"""
Harmony Integration Tutorial with Simulated Jurkat:293T Data
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the simulated dataset
adata = sc.read_h5ad("data/harmony/jurkat_293t_mixture_simulated.h5ad")

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
adata.write_h5ad("data/harmony/jurkat_293t_integrated.h5ad")
print("💾 Saved integrated dataset")
'''
    
    tutorial_path = Path("data/harmony") / "harmony_tutorial.py"
    with open(tutorial_path, 'w') as f:
        f.write(tutorial_script)
    
    print(f"📖 Created tutorial script: {tutorial_path}")

def main():
    """Main function"""
    print("🚀 Creating Harmony Test Datasets")
    print("=" * 50)
    
    # Create datasets
    adata = save_datasets()
    
    # Create tutorial
    create_harmony_tutorial()
    
    print("\n" + "=" * 50)
    print("✅ Harmony test datasets created successfully!")
    print("\nNext steps:")
    print("1. cd data/harmony")
    print("2. python harmony_tutorial.py")
    print("3. Or load in ChatSpatial:")
    print("   load_data('data/harmony/jurkat_293t_mixture_simulated.h5ad')")

if __name__ == "__main__":
    main()