#!/usr/bin/env python3
"""
Prepare and validate datasets for deconvolution testing.

This script downloads/creates suitable spatial and reference datasets for 
comprehensive deconvolution testing.
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path
from typing import Tuple, Dict, Any
import warnings

# Add the parent directory to the path to import chatspatial
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

def create_synthetic_reference_data(
    n_cells=2000,
    n_genes=1000,
    cell_types=None,
    save_path=None
):
    """Create synthetic single-cell reference data for deconvolution testing.
    
    Args:
        n_cells: Number of cells
        n_genes: Number of genes  
        cell_types: List of cell type names
        save_path: Path to save the data
        
    Returns:
        AnnData object with synthetic reference data
    """
    if cell_types is None:
        cell_types = ['T_cell', 'B_cell', 'Macrophage', 'NK_cell', 'Dendritic_cell', 'Fibroblast']
    
    print("Creating synthetic reference data: {} cells, {} genes, {} cell types".format(n_cells, n_genes, len(cell_types)))
    
    # Create expression matrix with cell type specific patterns
    np.random.seed(42)
    X = np.zeros((n_cells, n_genes))
    
    # Assign cell types
    n_cells_per_type = n_cells // len(cell_types)
    cell_type_labels = []
    
    for i, cell_type in enumerate(cell_types):
        start_idx = i * n_cells_per_type
        end_idx = (i + 1) * n_cells_per_type if i < len(cell_types) - 1 else n_cells
        
        # Create cell type specific expression patterns
        base_expression = np.random.poisson(5, size=(end_idx - start_idx, n_genes))
        
        # Add cell type specific marker genes (first 50 genes per cell type)
        marker_start = i * 50
        marker_end = marker_start + 50
        if marker_end <= n_genes:
            base_expression[:, marker_start:marker_end] += np.random.poisson(10, size=(end_idx - start_idx, 50))
        
        X[start_idx:end_idx] = base_expression
        cell_type_labels.extend([cell_type] * (end_idx - start_idx))
    
    # Create AnnData object
    adata = ad.AnnData(X=X.astype(np.int32))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs['cell_type'] = cell_type_labels
    
    # Add some metadata
    adata.obs['total_counts'] = np.array(X.sum(axis=1))
    adata.obs['n_genes'] = np.array((X > 0).sum(axis=1))
    adata.var['total_counts'] = np.array(X.sum(axis=0))
    adata.var['n_cells'] = np.array((X > 0).sum(axis=0))
    
    print(f"Created reference data with cell types: {adata.obs['cell_type'].value_counts().to_dict()}")
    
    if save_path:
        adata.write_h5ad(save_path)
        print(f"Saved reference data to {save_path}")
    
    return adata


def create_synthetic_spatial_data(
    n_spots=500,
    n_genes=1000,
    reference_adata=None,
    save_path=None
):
    """Create synthetic spatial data that matches reference data genes.
    
    Args:
        n_spots: Number of spatial spots
        n_genes: Number of genes (should match reference)
        reference_adata: Reference data to match genes
        save_path: Path to save the data
        
    Returns:
        AnnData object with synthetic spatial data
    """
    print(f"Creating synthetic spatial data: {n_spots} spots, {n_genes} genes")
    
    np.random.seed(123)
    
    # Create spatial coordinates in a grid pattern
    grid_size = int(np.ceil(np.sqrt(n_spots)))
    coords = []
    for i in range(n_spots):
        row = i // grid_size
        col = i % grid_size
        coords.append([col * 100, row * 100])
    spatial_coords = np.array(coords)
    
    # Create expression matrix - mixture of cell types with spatial patterns
    X = np.random.poisson(3, size=(n_spots, n_genes))
    
    # Add spatial patterns (e.g., some genes are higher in certain regions)
    for i in range(n_spots):
        x, y = spatial_coords[i]
        # Create gradients across space
        gradient_factor = 1 + 0.5 * np.sin(x / 100) * np.cos(y / 100)
        X[i] = X[i] * gradient_factor
        
    # Ensure non-negative integers
    X = np.maximum(X, 0).astype(np.int32)
    
    # Create AnnData object
    adata = ad.AnnData(X=X)
    adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
    
    # Use same gene names as reference if provided
    if reference_adata is not None:
        adata.var_names = reference_adata.var_names[:n_genes]
    else:
        adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    
    # Add spatial coordinates
    adata.obsm['spatial'] = spatial_coords
    
    # Add metadata
    adata.obs['total_counts'] = np.array(X.sum(axis=1))
    adata.obs['n_genes'] = np.array((X > 0).sum(axis=1))
    adata.var['total_counts'] = np.array(X.sum(axis=0))
    adata.var['n_cells'] = np.array((X > 0).sum(axis=0))
    
    print(f"Created spatial data with {n_spots} spots and {n_genes} genes")
    print(f"Mean counts per spot: {adata.obs['total_counts'].mean():.1f}")
    print(f"Spatial coordinates range: X({spatial_coords[:, 0].min():.1f}, {spatial_coords[:, 0].max():.1f}), "
          f"Y({spatial_coords[:, 1].min():.1f}, {spatial_coords[:, 1].max():.1f})")
    
    if save_path:
        adata.write_h5ad(save_path)
        print(f"Saved spatial data to {save_path}")
    
    return adata


def load_and_validate_existing_data():
    """Load and validate existing datasets for deconvolution testing.
    
    Returns:
        Tuple of (spatial_adata, reference_adata)
    """
    base_path = Path(__file__).parent.parent.parent / "data"
    
    # Try to load existing single cell data
    sc_data_path = base_path / "single_cell_datasets" / "single_cell.h5ad"
    reference_adata = None
    
    if sc_data_path.exists():
        try:
            print(f"Loading existing reference data from {sc_data_path}")
            reference_adata = sc.read_h5ad(sc_data_path)
            print(f"Loaded reference data: {reference_adata.shape}")
            
            # Check if cell_type column exists
            if 'cell_type' not in reference_adata.obs.columns:
                print("Adding synthetic cell_type labels to reference data...")
                n_types = 6
                cell_types = ['T_cell', 'B_cell', 'Macrophage', 'NK_cell', 'Dendritic_cell', 'Fibroblast']
                n_cells_per_type = reference_adata.n_obs // n_types
                
                cell_type_labels = []
                for i, cell_type in enumerate(cell_types):
                    start_idx = i * n_cells_per_type
                    end_idx = (i + 1) * n_cells_per_type if i < n_types - 1 else reference_adata.n_obs
                    cell_type_labels.extend([cell_type] * (end_idx - start_idx))
                
                reference_adata.obs['cell_type'] = cell_type_labels
                
            print(f"Reference cell types: {reference_adata.obs['cell_type'].value_counts().to_dict()}")
            
        except Exception as e:
            print(f"Failed to load existing reference data: {e}")
            reference_adata = None
    
    # Try to load existing spatial data
    spatial_data_paths = [
        base_path / "real_datasets" / "st_mouse_brain.h5ad",
        base_path / "spatial_datasets" / "squidpy_visium.h5ad",
        base_path / "demo_datasets" / "visium_demo.h5ad",
    ]
    
    spatial_adata = None
    for spatial_path in spatial_data_paths:
        if spatial_path.exists():
            try:
                print(f"Loading existing spatial data from {spatial_path}")
                spatial_adata = sc.read_h5ad(spatial_path)
                print(f"Loaded spatial data: {spatial_adata.shape}")
                
                # Ensure spatial coordinates exist
                if 'spatial' not in spatial_adata.obsm and 'X_spatial' in spatial_adata.obsm:
                    spatial_adata.obsm['spatial'] = spatial_adata.obsm['X_spatial']
                elif 'spatial' not in spatial_adata.obsm:
                    print("Adding synthetic spatial coordinates...")
                    n_spots = spatial_adata.n_obs
                    grid_size = int(np.ceil(np.sqrt(n_spots)))
                    coords = []
                    for i in range(n_spots):
                        row = i // grid_size
                        col = i % grid_size
                        coords.append([col * 100, row * 100])
                    spatial_adata.obsm['spatial'] = np.array(coords)
                
                break
                
            except Exception as e:
                print(f"Failed to load {spatial_path}: {e}")
                continue
    
    return spatial_adata, reference_adata


def prepare_datasets_for_deconvolution():
    """Prepare comprehensive datasets for deconvolution testing.
    
    Returns:
        Dictionary with prepared datasets and metadata
    """
    print("=== Preparing Datasets for Deconvolution Testing ===")
    
    # Create output directory
    output_dir = Path(__file__).parent / "deconvolution_test_data"
    output_dir.mkdir(exist_ok=True)
    
    results = {
        'datasets': {},
        'metadata': {},
        'output_dir': str(output_dir)
    }
    
    # 1. Try to load existing data
    spatial_adata, reference_adata = load_and_validate_existing_data()
    
    # 2. Create synthetic data if needed
    if reference_adata is None:
        print("\nCreating synthetic reference data...")
        reference_adata = create_synthetic_reference_data(
            n_cells=1500,
            n_genes=800,
            save_path=str(output_dir / "synthetic_reference.h5ad")
        )
    
    if spatial_adata is None:
        print("\nCreating synthetic spatial data...")
        spatial_adata = create_synthetic_spatial_data(
            n_spots=400,
            n_genes=800,
            reference_adata=reference_adata,
            save_path=str(output_dir / "synthetic_spatial.h5ad")
        )
    
    # 3. Ensure compatibility between datasets
    print("\nEnsuring dataset compatibility...")
    
    # Find common genes
    common_genes = list(set(spatial_adata.var_names) & set(reference_adata.var_names))
    print(f"Found {len(common_genes)} common genes between datasets")
    
    if len(common_genes) < 100:
        print("WARNING: Too few common genes. Creating compatible datasets...")
        
        # Create new datasets with shared genes
        n_shared_genes = 800
        shared_gene_names = [f"Gene_{i}" for i in range(n_shared_genes)]
        
        # Update reference data
        reference_adata_new = create_synthetic_reference_data(
            n_cells=1500,
            n_genes=n_shared_genes,
            save_path=str(output_dir / "compatible_reference.h5ad")
        )
        reference_adata_new.var_names = shared_gene_names
        
        # Update spatial data  
        spatial_adata_new = create_synthetic_spatial_data(
            n_spots=400,
            n_genes=n_shared_genes,
            reference_adata=reference_adata_new,
            save_path=str(output_dir / "compatible_spatial.h5ad")
        )
        spatial_adata_new.var_names = shared_gene_names
        
        reference_adata = reference_adata_new
        spatial_adata = spatial_adata_new
        common_genes = shared_gene_names
    
    # 4. Create smaller test datasets for quick testing
    print("\nCreating smaller test datasets...")
    
    # Small reference data (for quick testing)
    n_cells_small = min(500, reference_adata.n_obs)
    small_ref = reference_adata[:n_cells_small, :500].copy()
    small_ref.write_h5ad(str(output_dir / "small_reference.h5ad"))
    
    # Small spatial data
    n_spots_small = min(100, spatial_adata.n_obs)
    small_spatial = spatial_adata[:n_spots_small, :500].copy()
    small_spatial.write_h5ad(str(output_dir / "small_spatial.h5ad"))
    
    # 5. Store datasets and metadata
    results['datasets'] = {
        'reference': reference_adata,
        'spatial': spatial_adata,
        'small_reference': small_ref,
        'small_spatial': small_spatial
    }
    
    results['metadata'] = {
        'n_common_genes': len(common_genes),
        'reference_shape': reference_adata.shape,
        'spatial_shape': spatial_adata.shape,
        'cell_types': list(reference_adata.obs['cell_type'].unique()),
        'n_cell_types': len(reference_adata.obs['cell_type'].unique()),
        'spatial_coord_range': {
            'x_min': float(spatial_adata.obsm['spatial'][:, 0].min()),
            'x_max': float(spatial_adata.obsm['spatial'][:, 0].max()),
            'y_min': float(spatial_adata.obsm['spatial'][:, 1].min()),
            'y_max': float(spatial_adata.obsm['spatial'][:, 1].max()),
        }
    }
    
    print(f"\n=== Dataset Preparation Complete ===")
    print(f"Reference data: {reference_adata.shape}")
    print(f"Spatial data: {spatial_adata.shape}")
    print(f"Common genes: {len(common_genes)}")
    print(f"Cell types: {results['metadata']['cell_types']}")
    print(f"Output directory: {output_dir}")
    
    return results


def main():
    """Main function to prepare datasets."""
    try:
        results = prepare_datasets_for_deconvolution()
        
        # Print summary
        print("\n=== Summary ===")
        for name, adata in results['datasets'].items():
            print(f"{name}: {adata.shape}")
            if hasattr(adata, 'obs') and 'cell_type' in adata.obs:
                print(f"  Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")
        
        return results
        
    except Exception as e:
        print(f"Error preparing datasets: {e}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == "__main__":
    main()