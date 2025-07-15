#!/usr/bin/env python3
"""
Setup pancreas velocity data for MCP testing
"""

import sys
import os
import asyncio
import numpy as np
import scanpy as sc

# Add chatspatial to path
sys.path.append('/Users/apple/Research/SpatialTrans_MCP/chatspatial')

from chatspatial.models.data import SpatialDataset
from chatspatial.server import load_spatial_data

async def setup_pancreas_data():
    """Setup pancreas velocity data for MCP testing"""
    
    print("="*60)
    print("SETTING UP PANCREAS VELOCITY DATA FOR MCP")
    print("="*60)
    
    # Load the dataset
    adata = sc.read_h5ad('/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/pancreas_velocity.h5ad')
    
    # Add spatial coordinates (using UMAP)
    if 'X_umap' in adata.obsm:
        spatial_coords = adata.obsm['X_umap'].copy() * 100
        np.random.seed(42)
        spatial_coords += np.random.normal(0, 5, spatial_coords.shape)
        adata.obsm['spatial'] = spatial_coords
    
    print(f"Dataset prepared: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"Has velocity layers: {'spliced' in adata.layers and 'unspliced' in adata.layers}")
    print(f"Has spatial coordinates: {'spatial' in adata.obsm}")
    
    # Save the prepared dataset
    output_path = '/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/pancreas_velocity_with_spatial.h5ad'
    adata.write(output_path)
    print(f"Saved prepared dataset to: {output_path}")
    
    # Test loading through MCP interface
    try:
        result = await load_spatial_data(
            data_id='pancreas_velocity_spatial',
            file_path=output_path,
            data_type='velocity'
        )
        
        print(f"✓ Successfully loaded data through MCP interface!")
        print(f"  Result: {result}")
        
    except Exception as e:
        print(f"✗ Failed to load data through MCP interface: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*60)
    print("DATA SETUP COMPLETED")
    print("="*60)

if __name__ == "__main__":
    asyncio.run(setup_pancreas_data())