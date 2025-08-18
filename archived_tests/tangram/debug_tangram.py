#!/usr/bin/env python3
"""
Debug script to understand Tangram annotation issues
"""

import asyncio
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import tangram as tg
from unittest.mock import AsyncMock
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Add project root to path
import sys
sys.path.insert(0, str(Path(__file__).parent))

from chatspatial.tools.annotation import _annotate_with_tangram
from chatspatial.models.data import AnnotationParameters

def create_simple_test_data():
    """Create simple test data for debugging"""
    print("Creating simple test data...")
    
    # Create spatial data
    np.random.seed(42)
    n_spatial = 100
    n_genes = 200
    
    # Create genes including some marker genes
    genes = ['CD3D', 'CD3E', 'CD4', 'CD19', 'MS4A1', 'CD79A'] + [f'GENE_{i}' for i in range(6, n_genes)]
    
    # Spatial expression with some structure
    X_spatial = np.random.negative_binomial(5, 0.3, (n_spatial, n_genes)).astype(float)
    
    # Make T cell markers high in first half
    X_spatial[:50, 0:3] *= 3
    # Make B cell markers high in second half  
    X_spatial[50:, 3:6] *= 3
    
    adata_spatial = ad.AnnData(X=X_spatial)
    adata_spatial.obs_names = [f'spatial_{i}' for i in range(n_spatial)]
    adata_spatial.var_names = genes
    adata_spatial.obsm['spatial'] = np.random.uniform(0, 10, (n_spatial, 2))
    
    # Create reference data
    n_ref = 200
    X_ref = np.random.negative_binomial(8, 0.4, (n_ref, n_genes)).astype(float)
    
    # Make clear cell types in reference
    # T cells: first 100 cells
    X_ref[:100, 0:3] *= 5  # CD3D, CD3E, CD4
    # B cells: next 100 cells
    X_ref[100:, 3:6] *= 5  # CD19, MS4A1, CD79A
    
    adata_ref = ad.AnnData(X=X_ref)
    adata_ref.obs_names = [f'ref_{i}' for i in range(n_ref)]
    adata_ref.var_names = genes
    
    # Add cell type annotation to reference
    cell_types = ['T_cell'] * 100 + ['B_cell'] * 100
    adata_ref.obs['cell_type'] = pd.Categorical(cell_types)
    
    print(f"Spatial data: {adata_spatial.shape}")
    print(f"Reference data: {adata_ref.shape}")
    print(f"Reference cell types: {adata_ref.obs['cell_type'].value_counts()}")
    
    return adata_spatial, adata_ref

async def debug_tangram_annotation():
    """Debug the Tangram annotation process"""
    print("=" * 60)
    print("🧬 Debugging Tangram Annotation")
    print("=" * 60)
    
    # Create test data
    adata_spatial, adata_ref = create_simple_test_data()
    
    # Create mock context
    context = AsyncMock()
    
    # Create data store
    data_store = {
        'reference': {'adata': adata_ref, 'name': 'Reference Data'}
    }
    
    # Test parameters
    params = AnnotationParameters(
        method='tangram',
        reference_data_id='reference',
        mode='cells',
        num_epochs=50
    )
    
    print("\n1. Testing basic Tangram annotation...")
    try:
        cell_types, counts, confidence_scores, mapping_score = await _annotate_with_tangram(
            adata_spatial, params, data_store, context
        )
        
        print(f"✅ Tangram completed successfully!")
        print(f"   - Cell types: {cell_types}")
        print(f"   - Counts: {counts}")  
        print(f"   - Confidence scores: {confidence_scores}")
        print(f"   - Mapping score: {mapping_score}")
        print(f"   - adata.obs columns: {list(adata_spatial.obs.columns)}")
        print(f"   - adata.obsm keys: {list(adata_spatial.obsm.keys())}")
        
        # Check if cell_type was added
        if 'cell_type' in adata_spatial.obs.columns:
            print(f"   - Cell type assignments: {adata_spatial.obs['cell_type'].value_counts()}")
        
        # Check obsm contents
        if 'tangram_ct_pred' in adata_spatial.obsm:
            pred_df = adata_spatial.obsm['tangram_ct_pred']
            print(f"   - Prediction matrix shape: {pred_df.shape}")
            print(f"   - Predicted cell types: {list(pred_df.columns)}")
            print(f"   - Sample predictions (first 5 cells):")
            print(pred_df.head())
        else:
            print("   - No tangram_ct_pred found in obsm")
            
    except Exception as e:
        print(f"❌ Tangram failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n2. Testing clusters mode...")
    params_clusters = AnnotationParameters(
        method='tangram',
        reference_data_id='reference',
        mode='clusters',
        cluster_label='cell_type',
        num_epochs=50
    )
    
    try:
        adata_spatial_2 = adata_spatial.copy()
        cell_types2, counts2, confidence_scores2, mapping_score2 = await _annotate_with_tangram(
            adata_spatial_2, params_clusters, data_store, context
        )
        
        print(f"✅ Clusters mode completed!")
        print(f"   - Cell types: {cell_types2}")
        print(f"   - Counts: {counts2}")
        print(f"   - Confidence scores: {confidence_scores2}")
        print(f"   - Mapping score: {mapping_score2}")
        
    except Exception as e:
        print(f"❌ Clusters mode failed: {str(e)}")
        import traceback
        traceback.print_exc()
    
    return True

if __name__ == "__main__":
    success = asyncio.run(debug_tangram_annotation())
    print("\n" + "=" * 60)
    if success:
        print("🎉 Debug completed successfully!")
    else:
        print("❌ Debug failed!")
    print("=" * 60)