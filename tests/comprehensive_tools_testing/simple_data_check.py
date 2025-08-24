#!/usr/bin/env python3
"""
Simple script to check existing data and prepare for deconvolution testing.
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

def main():
    """Check existing data and prepare simple datasets."""
    print("=== Checking Existing Data ===")
    
    base_path = Path(__file__).parent.parent.parent / "data"
    
    # Check single cell data
    sc_path = base_path / "single_cell_datasets" / "single_cell.h5ad"
    print("Single cell data path:", sc_path)
    print("Single cell data exists:", sc_path.exists())
    
    if sc_path.exists():
        try:
            adata = sc.read_h5ad(sc_path)
            print("Single cell data shape:", adata.shape)
            print("Single cell obs columns:", list(adata.obs.columns))
            print("First few obs index:", list(adata.obs.index[:5]))
        except Exception as e:
            print("Error reading single cell data:", e)
    
    # Check spatial datasets
    spatial_paths = [
        base_path / "real_datasets" / "st_mouse_brain.h5ad",
        base_path / "spatial_datasets" / "squidpy_visium.h5ad",
        base_path / "demo_datasets" / "visium_demo.h5ad",
    ]
    
    for path in spatial_paths:
        print("\nSpatial data path:", path)
        print("Exists:", path.exists())
        if path.exists():
            try:
                adata = sc.read_h5ad(path)
                print("Shape:", adata.shape)
                print("obsm keys:", list(adata.obsm.keys()) if hasattr(adata, 'obsm') else "No obsm")
            except Exception as e:
                print("Error reading:", e)
    
    # Create simple synthetic datasets
    print("\n=== Creating Simple Test Datasets ===")
    
    output_dir = Path(__file__).parent / "deconvolution_test_data"
    output_dir.mkdir(exist_ok=True)
    
    # Create reference data
    n_cells = 500
    n_genes = 200
    cell_types = ['T_cell', 'B_cell', 'Macrophage']
    
    print("Creating reference data: {} cells, {} genes".format(n_cells, n_genes))
    
    # Simple count matrix
    X_ref = np.random.poisson(5, size=(n_cells, n_genes))
    
    # Create cell type labels
    cell_type_labels = []
    cells_per_type = n_cells // len(cell_types)
    for i, ct in enumerate(cell_types):
        start_idx = i * cells_per_type
        end_idx = (i + 1) * cells_per_type if i < len(cell_types) - 1 else n_cells
        cell_type_labels.extend([ct] * (end_idx - start_idx))
    
    # Create reference AnnData
    ref_adata = ad.AnnData(X=X_ref.astype(np.int32))
    ref_adata.obs_names = ["cell_{}".format(i) for i in range(n_cells)]
    ref_adata.var_names = ["Gene_{}".format(i) for i in range(n_genes)]
    ref_adata.obs['cell_type'] = cell_type_labels
    
    ref_path = str(output_dir / "simple_reference.h5ad")
    ref_adata.write_h5ad(ref_path)
    print("Saved reference to:", ref_path)
    print("Reference cell types:", ref_adata.obs['cell_type'].value_counts().to_dict())
    
    # Create spatial data
    n_spots = 100
    print("Creating spatial data: {} spots, {} genes".format(n_spots, n_genes))
    
    # Simple count matrix
    X_spatial = np.random.poisson(3, size=(n_spots, n_genes))
    
    # Create spatial AnnData
    spatial_adata = ad.AnnData(X=X_spatial.astype(np.int32))
    spatial_adata.obs_names = ["spot_{}".format(i) for i in range(n_spots)]
    spatial_adata.var_names = ["Gene_{}".format(i) for i in range(n_genes)]  # Same as reference
    
    # Add spatial coordinates
    coords = []
    grid_size = int(np.ceil(np.sqrt(n_spots)))
    for i in range(n_spots):
        row = i // grid_size
        col = i % grid_size
        coords.append([col * 10, row * 10])
    spatial_adata.obsm['spatial'] = np.array(coords)
    
    spatial_path = str(output_dir / "simple_spatial.h5ad")
    spatial_adata.write_h5ad(spatial_path)
    print("Saved spatial to:", spatial_path)
    
    # Verify common genes
    common_genes = list(set(ref_adata.var_names) & set(spatial_adata.var_names))
    print("Common genes:", len(common_genes))
    
    print("\n=== Test Data Ready ===")
    print("Reference:", ref_adata.shape)
    print("Spatial:", spatial_adata.shape)
    print("Cell types:", list(ref_adata.obs['cell_type'].unique()))
    print("Output directory:", output_dir)
    
    return {
        'reference_path': ref_path,
        'spatial_path': spatial_path,
        'reference_adata': ref_adata,
        'spatial_adata': spatial_adata,
        'output_dir': str(output_dir)
    }

if __name__ == "__main__":
    main()