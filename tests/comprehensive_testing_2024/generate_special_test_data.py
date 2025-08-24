"""
Generate Special Test Datasets for Error Handling and Compatibility Testing

Linus principle: "Generate test data that breaks things in realistic ways."

This script creates additional test datasets that current synthetic data might not cover:
- Corrupted data with specific patterns
- Edge case formats 
- Version-specific compatibility issues
- Platform-specific data encoding problems
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from typing import Dict, Any, List
import warnings

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../'))

from chatspatial.utils.error_handling import suppress_output


def create_corrupted_h5ad_datasets(output_dir: str) -> Dict[str, str]:
    """
    Create datasets that test specific corruption/error scenarios.
    These test real failure modes we've seen in production.
    """
    datasets = {}
    
    with suppress_output():
        # 1. Dataset with NaN coordinates
        print("Creating dataset with NaN spatial coordinates...")
        adata_nan_coords = ad.AnnData(
            X=np.random.poisson(5, (100, 500)).astype(np.float32),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(100)]),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(500)])
        )
        # Add spatial coordinates with NaN values
        coords = np.random.uniform(0, 100, (100, 2))
        coords[10:15, :] = np.nan  # Inject NaN values
        adata_nan_coords.obsm['spatial'] = coords
        
        nan_coords_path = os.path.join(output_dir, "nan_coordinates.h5ad")
        adata_nan_coords.write(nan_coords_path)
        datasets["nan_coordinates"] = nan_coords_path
        
        # 2. Dataset with identical coordinates (no spatial variation)
        print("Creating dataset with identical coordinates...")
        adata_identical = ad.AnnData(
            X=np.random.poisson(3, (50, 200)).astype(np.float32),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(50)]),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(200)])
        )
        # All cells at same location
        adata_identical.obsm['spatial'] = np.tile([50.0, 50.0], (50, 1))
        
        identical_coords_path = os.path.join(output_dir, "identical_coordinates.h5ad")
        adata_identical.write(identical_coords_path)
        datasets["identical_coordinates"] = identical_coords_path
        
        # 3. Dataset that will test dimension validation
        print("Creating dataset for dimension validation testing...")
        adata_dim_test = ad.AnnData(
            X=np.random.poisson(2, (75, 300)).astype(np.float32),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(75)]),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(300)])
        )
        # Correct dimensions for writing, but we'll use this to test validation
        adata_dim_test.obsm['spatial'] = np.random.uniform(0, 100, (75, 2))
        
        dim_test_path = os.path.join(output_dir, "dimension_validation_test.h5ad")
        adata_dim_test.write(dim_test_path)
        datasets["dimension_validation_test"] = dim_test_path
        
        # 4. Dataset with extreme coordinate values
        print("Creating dataset with extreme coordinate values...")
        adata_extreme = ad.AnnData(
            X=np.random.poisson(4, (60, 250)).astype(np.float32),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(60)]),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(250)])
        )
        # Mix of normal and extreme coordinates
        coords_extreme = np.random.uniform(0, 100, (60, 2))
        coords_extreme[:5, :] = [np.inf, np.inf]  # Infinite coordinates
        coords_extreme[5:10, :] = [1e10, 1e10]    # Very large coordinates
        coords_extreme[10:15, :] = [-1e6, -1e6]   # Very negative coordinates
        adata_extreme.obsm['spatial'] = coords_extreme
        
        extreme_coords_path = os.path.join(output_dir, "extreme_coordinates.h5ad")
        adata_extreme.write(extreme_coords_path)
        datasets["extreme_coordinates"] = extreme_coords_path
        
        # 5. Dataset with wrong spatial dimension (1D instead of 2D)
        print("Creating dataset with 1D coordinates...")
        adata_1d = ad.AnnData(
            X=np.random.poisson(3, (40, 150)).astype(np.float32),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(40)]),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(150)])
        )
        # 1D coordinates instead of 2D
        adata_1d.obsm['spatial'] = np.random.uniform(0, 100, (40, 1))
        
        coords_1d_path = os.path.join(output_dir, "coordinates_1d.h5ad")
        adata_1d.write(coords_1d_path)
        datasets["coordinates_1d"] = coords_1d_path
        
        # 6. Dataset with missing required layers for velocity analysis
        print("Creating dataset missing velocity layers...")
        adata_no_velocity = ad.AnnData(
            X=np.random.poisson(5, (80, 400)).astype(np.float32),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(80)]),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(400)])
        )
        adata_no_velocity.obsm['spatial'] = np.random.uniform(0, 100, (80, 2))
        # Deliberately missing 'spliced' and 'unspliced' layers
        
        no_velocity_path = os.path.join(output_dir, "missing_velocity_layers.h5ad")
        adata_no_velocity.write(no_velocity_path)
        datasets["missing_velocity_layers"] = no_velocity_path
        
        # 7. Dataset with empty/all-zero layers
        print("Creating dataset with empty layers...")
        adata_empty_layers = ad.AnnData(
            X=np.random.poisson(3, (70, 350)).astype(np.float32),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(70)]),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(350)])
        )
        adata_empty_layers.obsm['spatial'] = np.random.uniform(0, 100, (70, 2))
        # Add empty layers
        adata_empty_layers.layers['spliced'] = np.zeros((70, 350))
        adata_empty_layers.layers['unspliced'] = np.zeros((70, 350))
        
        empty_layers_path = os.path.join(output_dir, "empty_velocity_layers.h5ad")
        adata_empty_layers.write(empty_layers_path)
        datasets["empty_velocity_layers"] = empty_layers_path
    
    return datasets


def create_format_compatibility_datasets(output_dir: str) -> Dict[str, str]:
    """
    Create datasets for testing format compatibility issues.
    """
    datasets = {}
    
    with suppress_output():
        # Base dataset
        base_adata = ad.AnnData(
            X=np.random.poisson(4, (100, 500)).astype(np.float32),
            obs=pd.DataFrame({
                'cell_type': np.random.choice(['TypeA', 'TypeB', 'TypeC'], 100),
                'batch': np.random.choice([1, 2, 3], 100)
            }, index=[f"cell_{i}" for i in range(100)]),
            var=pd.DataFrame({
                'highly_variable': np.random.choice([True, False], 500),
                'gene_length': np.random.randint(1000, 5000, 500)
            }, index=[f"gene_{i}" for i in range(500)])
        )
        base_adata.obsm['spatial'] = np.random.uniform(0, 100, (100, 2))
        
        # 1. Dataset with string indices (compatibility issue)
        print("Creating dataset with string indices...")
        adata_strings = base_adata.copy()
        adata_strings.obs.index = [f"cell_string_{i}" for i in range(100)]
        adata_strings.var.index = [f"GENE_STRING_{i}" for i in range(500)]
        
        string_indices_path = os.path.join(output_dir, "string_indices.h5ad")
        adata_strings.write(string_indices_path)
        datasets["string_indices"] = string_indices_path
        
        # 2. Dataset with mixed data types (all strings to avoid type conflicts)
        print("Creating dataset with mixed data types...")
        adata_mixed = base_adata.copy()
        adata_mixed.obs['mixed_column'] = [f"str_{i}" if i % 2 == 0 else f"num_{i}" for i in range(100)]
        adata_mixed.var['mixed_var'] = [f"float_{i}" if i % 3 == 0 else f"str_{i}" for i in range(500)]
        
        mixed_types_path = os.path.join(output_dir, "mixed_data_types.h5ad")
        adata_mixed.write(mixed_types_path)
        datasets["mixed_data_types"] = mixed_types_path
        
        # 3. Dataset with unicode/special characters
        print("Creating dataset with unicode characters...")
        adata_unicode = base_adata.copy()
        adata_unicode.obs['unicode_col'] = [f"测试_{i}α" for i in range(100)]
        adata_unicode.var.index = [f"基因_{i}β" for i in range(500)]
        
        unicode_path = os.path.join(output_dir, "unicode_characters.h5ad")
        adata_unicode.write(unicode_path)
        datasets["unicode_characters"] = unicode_path
        
        # 4. Very sparse dataset (>99% zeros)
        print("Creating extremely sparse dataset...")
        sparse_X = np.zeros((100, 500))
        # Only add a few non-zero values
        for _ in range(50):  # Only 50 non-zero entries out of 50,000
            i, j = np.random.randint(0, 100), np.random.randint(0, 500)
            sparse_X[i, j] = np.random.randint(1, 10)
        
        adata_very_sparse = ad.AnnData(
            X=sparse_X,
            obs=base_adata.obs.copy(),
            var=base_adata.var.copy()
        )
        adata_very_sparse.obsm['spatial'] = base_adata.obsm['spatial'].copy()
        
        very_sparse_path = os.path.join(output_dir, "extremely_sparse.h5ad")
        adata_very_sparse.write(very_sparse_path)
        datasets["extremely_sparse"] = very_sparse_path
    
    return datasets


def create_legacy_parameter_test_data(output_dir: str) -> Dict[str, Any]:
    """
    Create test cases for legacy parameter compatibility.
    This tests old parameter names/formats that should still work.
    """
    test_cases = {}
    
    # Legacy parameter formats that should be supported
    test_cases["legacy_clustering_params"] = {
        # Old parameter name vs new parameter name
        "clustering_key": "cluster_key",  # Old scanpy style
        "cluster_method": "clustering_method",
        "n_clusters": "k_clusters",
    }
    
    test_cases["legacy_spatial_params"] = {
        "coord_key": "spatial_key",
        "coordinate_name": "spatial_key",
        "spatial_coordinates": "spatial_key"
    }
    
    test_cases["legacy_plotting_params"] = {
        "color_by": "color",
        "size_by": "size", 
        "group_by": "groupby"
    }
    
    # Save parameter mapping test data
    import json
    param_test_path = os.path.join(output_dir, "legacy_parameter_mappings.json")
    with open(param_test_path, 'w') as f:
        json.dump(test_cases, f, indent=2)
    
    return {"legacy_parameter_mappings": param_test_path}


def generate_all_special_datasets():
    """
    Main function to generate all special test datasets.
    """
    output_dir = os.path.join(os.path.dirname(__file__), "datasets")
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 50)
    print("Generating Special Test Datasets for Error Handling")
    print("=" * 50)
    
    # Track all generated datasets
    all_datasets = {}
    
    try:
        # Generate corrupted datasets
        print("\n1. Creating corrupted/problematic datasets...")
        corrupted = create_corrupted_h5ad_datasets(output_dir)
        all_datasets.update(corrupted)
        
        # Generate format compatibility datasets
        print("\n2. Creating format compatibility datasets...")
        format_compat = create_format_compatibility_datasets(output_dir)
        all_datasets.update(format_compat)
        
        # Generate legacy parameter test data
        print("\n3. Creating legacy parameter test data...")
        legacy_params = create_legacy_parameter_test_data(output_dir)
        all_datasets.update(legacy_params)
        
        # Update summary
        print("\n4. Updating dataset summary...")
        update_datasets_summary(output_dir, all_datasets)
        
        print(f"\n✓ Generated {len(all_datasets)} special test datasets")
        print(f"✓ Datasets saved to: {output_dir}")
        
        return all_datasets
        
    except Exception as e:
        print(f"✗ Error generating datasets: {e}")
        return {}


def update_datasets_summary(output_dir: str, new_datasets: Dict[str, str]):
    """Update the datasets summary CSV with new datasets"""
    summary_path = os.path.join(output_dir, "datasets_summary.csv")
    
    # Read existing summary if it exists
    if os.path.exists(summary_path):
        existing_df = pd.read_csv(summary_path)
    else:
        existing_df = pd.DataFrame(columns=['dataset', 'n_cells', 'n_genes', 'has_spatial', 'file_size_mb', 'sparsity'])
    
    # Analyze new datasets
    new_rows = []
    for name, path in new_datasets.items():
        if path.endswith('.h5ad'):
            try:
                with suppress_output():
                    adata = sc.read(path)
                    
                file_size_mb = os.path.getsize(path) / (1024 * 1024)
                has_spatial = 'spatial' in adata.obsm
                
                # Calculate sparsity
                if hasattr(adata.X, 'nnz'):
                    sparsity = 1 - (adata.X.nnz / (adata.n_obs * adata.n_vars))
                else:
                    sparsity = 1 - np.count_nonzero(adata.X) / (adata.n_obs * adata.n_vars)
                
                new_rows.append({
                    'dataset': os.path.basename(path),
                    'n_cells': adata.n_obs,
                    'n_genes': adata.n_vars,
                    'has_spatial': has_spatial,
                    'file_size_mb': file_size_mb,
                    'sparsity': sparsity
                })
            except Exception as e:
                # For datasets that are intentionally corrupted
                new_rows.append({
                    'dataset': os.path.basename(path),
                    'n_cells': 'corrupted',
                    'n_genes': 'corrupted',
                    'has_spatial': 'unknown',
                    'file_size_mb': os.path.getsize(path) / (1024 * 1024) if os.path.exists(path) else 0,
                    'sparsity': 'unknown'
                })
    
    # Append new rows if any
    if new_rows:
        new_df = pd.DataFrame(new_rows)
        updated_df = pd.concat([existing_df, new_df], ignore_index=True)
        updated_df.to_csv(summary_path, index=False)
        print(f"✓ Updated datasets summary with {len(new_rows)} new entries")


if __name__ == "__main__":
    generated = generate_all_special_datasets()
    
    print("\n" + "=" * 50)
    print("Generation Summary:")
    print("=" * 50)
    for name, path in generated.items():
        print(f"✓ {name}: {os.path.basename(path)}")
    print("\nSpecial datasets ready for error handling and compatibility testing!")