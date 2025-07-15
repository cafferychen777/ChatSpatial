#!/usr/bin/env python3
"""
Load pancreas velocity dataset and add spatial coordinates for testing
"""

import scanpy as sc
import numpy as np
import pandas as pd
import sys
import os

# Add chatspatial to path
sys.path.append('/Users/apple/Research/SpatialTrans_MCP/chatspatial')

def load_pancreas_velocity_data():
    """Load and prepare pancreas velocity dataset with spatial coordinates"""
    
    # Load the dataset
    adata = sc.read_h5ad('/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/pancreas_velocity.h5ad')
    
    print(f"Loaded pancreas dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Use UMAP coordinates as "spatial" coordinates for testing
    # This allows us to test trajectory analysis on this velocity dataset
    if 'X_umap' in adata.obsm:
        spatial_coords = adata.obsm['X_umap'].copy()
        # Scale coordinates to make them more "spatial-like"
        spatial_coords = spatial_coords * 100  # Scale up
        # Add some noise to make it more realistic
        np.random.seed(42)
        spatial_coords += np.random.normal(0, 5, spatial_coords.shape)
        
        # Store as spatial coordinates
        adata.obsm['spatial'] = spatial_coords
        print("✓ Added spatial coordinates based on UMAP")
    else:
        print("✗ No UMAP coordinates found, creating random spatial coordinates")
        np.random.seed(42)
        spatial_coords = np.random.randn(adata.n_obs, 2) * 50
        adata.obsm['spatial'] = spatial_coords
    
    # Print dataset info
    print("\nDataset Information:")
    print(f"Shape: {adata.shape}")
    print(f"Layers: {list(adata.layers.keys())}")
    print(f"Obsm keys: {list(adata.obsm.keys())}")
    print(f"Has spliced/unspliced: {'spliced' in adata.layers and 'unspliced' in adata.layers}")
    print(f"Has spatial coordinates: {'spatial' in adata.obsm}")
    
    # Show some cluster info
    if 'clusters' in adata.obs:
        print(f"\nCell type clusters: {adata.obs['clusters'].value_counts().head()}")
    
    return adata

if __name__ == "__main__":
    adata = load_pancreas_velocity_data()
    
    # Test that we can import and use the trajectory analysis
    try:
        from chatspatial.tools.trajectory import analyze_rna_velocity, analyze_trajectory
        print("\n✓ Successfully imported trajectory analysis functions")
    except ImportError as e:
        print(f"\n✗ Failed to import trajectory functions: {e}")
    
    print("\n" + "="*50)
    print("Dataset ready for RNA velocity and trajectory analysis!")
    print("="*50)