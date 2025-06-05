"""
Data loading utilities for spatial transcriptomics data.
"""

from typing import Dict, Optional, Any, Literal
import os
import numpy as np
import scanpy as sc
import squidpy as sq
import anndata as ad


async def load_spatial_data(
    data_path: str,
    data_type: Literal["10x_visium", "slide_seq", "merfish", "seqfish", "other", "auto"] = "auto",
    name: Optional[str] = None
) -> Dict[str, Any]:
    """Load spatial transcriptomics data

    Args:
        data_path: Path to the data file or directory
        data_type: Type of spatial data. If 'auto', will try to determine the type from the file extension or directory structure.
        name: Optional name for the dataset

    Returns:
        Dictionary with dataset information and AnnData object
    """
    # Validate path
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data path not found: {data_path}")

    # Auto-detect data type if set to 'auto'
    if data_type == "auto":
        if os.path.isfile(data_path) and data_path.endswith('.h5ad'):
            # It's an h5ad file
            data_type = "h5ad"
        elif os.path.isdir(data_path):
            # Check if it has the structure of a 10x Visium dataset
            if (os.path.exists(os.path.join(data_path, 'filtered_feature_bc_matrix')) or
                os.path.exists(os.path.join(data_path, 'filtered_feature_bc_matrix.h5'))):
                data_type = "10x_visium"
            else:
                # Default to other if we can't determine
                data_type = "other"
        else:
            # Default to other for unknown file types
            data_type = "other"

    # Convert h5ad to other for backward compatibility
    if data_type == "h5ad":
        data_type = "other"

    # Load data based on data_type
    if data_type == "10x_visium":
        # For 10x Visium, we need to provide the path to the directory containing the data
        try:
            # Check if it's a directory or an h5ad file
            if os.path.isdir(data_path):
                # Check if the directory has the expected structure
                if os.path.exists(os.path.join(data_path, 'filtered_feature_bc_matrix')):
                    # Standard 10x Visium directory structure
                    adata = sc.read_visium(data_path)
                elif os.path.exists(os.path.join(data_path, 'filtered_feature_bc_matrix.h5')):
                    # H5 file based 10x Visium directory structure
                    adata = sc.read_visium(data_path)
                elif os.path.exists(os.path.join(data_path, 'filtered_feature_bc_matrix', 'matrix.mtx.gz')):
                    # Matrix files based 10x Visium directory structure
                    # Use scanpy's read_10x_mtx function
                    adata = sc.read_10x_mtx(
                        os.path.join(data_path, 'filtered_feature_bc_matrix'),
                        var_names='gene_symbols',
                        cache=True
                    )
                    # Try to load spatial coordinates if available
                    spatial_dir = os.path.join(data_path, 'spatial')
                    if os.path.exists(spatial_dir):
                        try:
                            # Add spatial information manually
                            import json
                            import pandas as pd

                            # Load tissue positions
                            positions_path = os.path.join(spatial_dir, 'tissue_positions_list.csv')
                            if os.path.exists(positions_path):
                                positions = pd.read_csv(positions_path, header=None)
                                positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
                                positions.set_index('barcode', inplace=True)

                                # Filter for spots in tissue
                                positions = positions[positions['in_tissue'] == 1]

                                # Add spatial coordinates to adata
                                adata.obsm['spatial'] = positions.loc[adata.obs_names, ['pxl_col_in_fullres', 'pxl_row_in_fullres']].values

                                # Load scalefactors
                                scalefactors_path = os.path.join(spatial_dir, 'scalefactors_json.json')
                                if os.path.exists(scalefactors_path):
                                    with open(scalefactors_path, 'r') as f:
                                        scalefactors = json.load(f)

                                    # Add scalefactors to adata
                                    adata.uns['spatial'] = {
                                        'scalefactors': scalefactors
                                    }
                        except Exception as e:
                            print(f"Warning: Could not load spatial information: {str(e)}")
                else:
                    raise ValueError(f"Directory {data_path} does not have the expected 10x Visium structure")
            elif os.path.isfile(data_path) and data_path.endswith('.h5ad'):
                # If it's an h5ad file but marked as 10x_visium, read it as h5ad
                adata = sc.read_h5ad(data_path)
                # Check if it has the necessary spatial information
                if 'spatial' not in adata.uns and not any('spatial' in key for key in adata.obsm.keys()):
                    raise ValueError("The h5ad file does not contain spatial information required for 10x Visium data")
            else:
                raise ValueError(f"Unsupported file format for 10x_visium: {data_path}")

            # Add spatial neighborhood graph if not already present
            if 'spatial_connectivities' not in adata.obsp and 'spatial' in adata.obsm:
                try:
                    sq.gr.spatial_neighbors(adata)
                except Exception as e:
                    print(f"Warning: Could not compute spatial neighbors: {str(e)}")
        except Exception as e:
            raise ValueError(f"Error loading 10x Visium data: {str(e)}")
    elif data_type == "h5ad" or data_type in ["slide_seq", "merfish", "seqfish", "other"]:
        # For h5ad files or other data types
        try:
            adata = sc.read_h5ad(data_path)
        except Exception as e:
            raise ValueError(f"Error loading {data_type} data: {str(e)}")
    else:
        raise ValueError(f"Unsupported data type: {data_type}")

    # Set dataset name
    dataset_name = name or os.path.basename(data_path).split('.')[0]

    # Calculate basic statistics
    n_cells = adata.n_obs
    n_genes = adata.n_vars

    # Check if spatial coordinates are available
    spatial_coordinates_available = 'spatial' in adata.uns or hasattr(adata, 'obsm') and any('spatial' in key for key in adata.obsm.keys())

    # Check if tissue image is available (for Visium data)
    tissue_image_available = data_type == "10x_visium" and 'spatial' in adata.uns and 'images' in adata.uns['spatial']

    # Make variable names unique to avoid reindexing issues
    if hasattr(adata, 'var_names_make_unique'):
        adata.var_names_make_unique()

    # Return dataset info and AnnData object
    return {
        "name": dataset_name,
        "type": data_type,
        "path": data_path,
        "adata": adata,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "spatial_coordinates_available": spatial_coordinates_available,
        "tissue_image_available": tissue_image_available
    }