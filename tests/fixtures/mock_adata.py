"""
Mock AnnData object generation utilities

This module provides tools for creating test AnnData objects without relying on real datasets.
"""
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from typing import Optional, List


def create_mock_adata(
    n_obs: int = 100,
    n_vars: int = 200,
    add_spatial: bool = False,
    add_clusters: bool = False,
    add_cell_types: bool = False,
    sparse_matrix: bool = True,
    random_seed: int = 42,
) -> AnnData:
    """
    Create mock AnnData object for testing

    Args:
        n_obs: Number of cells/spots
        n_vars: Number of genes
        add_spatial: Whether to add spatial coordinates
        add_clusters: Whether to add cluster labels (leiden)
        add_cell_types: Whether to add cell type labels
        sparse_matrix: Whether to use sparse matrix (recommended)
        random_seed: Random seed for reproducibility

    Returns:
        Mock AnnData object in ChatSpatial standard format
    """
    np.random.seed(random_seed)

    # Generate expression matrix
    if sparse_matrix:
        # Sparse matrix (90% zeros, mimics real scRNA-seq data)
        X = sparse.random(n_obs, n_vars, density=0.1, format="csr", random_state=random_seed)
        X.data = np.abs(X.data) * 10  # Positive expression values
    else:
        X = np.random.rand(n_obs, n_vars) * 10

    # Create obs (cell metadata)
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])

    # Create var (gene metadata)
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])

    adata = AnnData(X=X, obs=obs, var=var)

    # Add spatial coordinates
    if add_spatial:
        spatial_coords = np.random.rand(n_obs, 2) * 100
        adata.obsm["spatial"] = spatial_coords

    # Add cluster labels
    if add_clusters:
        n_clusters = min(5, max(2, n_obs // 20))  # 2-5 clusters
        adata.obs["leiden"] = [str(i % n_clusters) for i in range(n_obs)]
        adata.obs["leiden"] = adata.obs["leiden"].astype("category")

    # Add cell type labels
    if add_cell_types:
        cell_types = ["TypeA", "TypeB", "TypeC", "TypeD"]
        n_types = min(len(cell_types), max(2, n_obs // 25))
        adata.obs["cell_type"] = np.random.choice(
            cell_types[:n_types], size=n_obs
        )
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    return adata


def create_mock_reference_adata(
    n_obs: int = 500,
    n_vars: int = 200,
    cell_types: Optional[List[str]] = None,
    random_seed: int = 42,
) -> AnnData:
    """
    Create reference dataset for cell annotation

    Args:
        n_obs: Number of cells (reference datasets are usually larger)
        n_vars: Number of genes
        cell_types: List of cell types
        random_seed: Random seed

    Returns:
        Reference AnnData object with cell_type labels
    """
    if cell_types is None:
        cell_types = ["T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial"]

    np.random.seed(random_seed)

    adata = create_mock_adata(
        n_obs=n_obs,
        n_vars=n_vars,
        sparse_matrix=True,
        random_seed=random_seed,
    )

    # Add cell type labels (required for reference datasets)
    adata.obs["cell_type"] = np.random.choice(cell_types, size=n_obs)
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    return adata


def create_mock_velocity_adata(
    n_obs: int = 100, n_vars: int = 200, random_seed: int = 42
) -> AnnData:
    """
    Create AnnData with RNA velocity layers

    Args:
        n_obs: Number of cells
        n_vars: Number of genes
        random_seed: Random seed

    Returns:
        AnnData with spliced/unspliced layers
    """
    adata = create_mock_adata(
        n_obs=n_obs,
        n_vars=n_vars,
        sparse_matrix=True,
        random_seed=random_seed,
    )

    np.random.seed(random_seed)

    # Add spliced and unspliced layers
    adata.layers["spliced"] = sparse.random(
        n_obs, n_vars, density=0.1, format="csr", random_state=random_seed
    )
    adata.layers["spliced"].data = np.abs(adata.layers["spliced"].data) * 10

    adata.layers["unspliced"] = sparse.random(
        n_obs, n_vars, density=0.05, format="csr", random_state=random_seed + 1
    )
    adata.layers["unspliced"].data = np.abs(adata.layers["unspliced"].data) * 5

    return adata


def generate_sample_datasets(output_dir: str = "tests/fixtures/sample_data"):
    """
    Generate all sample datasets for testing (run once for CI/CD)

    This function creates small test dataset files to avoid regenerating them for each test run.

    Args:
        output_dir: Output directory path
    """
    import os

    os.makedirs(output_dir, exist_ok=True)

    print(f"📦 Generating test datasets to {output_dir}/")

    # 1. Visium sample (100 spots × 100 genes)
    print("  → visium_sample.h5ad (100 spots × 100 genes)")
    visium = create_mock_adata(
        n_obs=100, n_vars=100, add_spatial=True, add_clusters=True, random_seed=42
    )
    visium.write_h5ad(f"{output_dir}/visium_sample.h5ad")

    # 2. Reference dataset (500 cells × 100 genes)
    print("  → reference.h5ad (500 cells × 100 genes, 5 cell types)")
    reference = create_mock_reference_adata(n_obs=500, n_vars=100, random_seed=42)
    reference.write_h5ad(f"{output_dir}/reference.h5ad")

    # 3. Velocity dataset (100 cells × 100 genes with spliced/unspliced)
    print("  → velocity_sample.h5ad (100 cells × 100 genes, spliced/unspliced)")
    velocity = create_mock_velocity_adata(n_obs=100, n_vars=100, random_seed=44)
    velocity.write_h5ad(f"{output_dir}/velocity_sample.h5ad")

    print(f"✅ Done! Generated 3 test datasets")
    total_size = sum(os.path.getsize(f'{output_dir}/{f}') for f in os.listdir(output_dir)) / 1024 / 1024
    print(f"   Total size: ~{total_size:.2f} MB")


if __name__ == "__main__":
    # Generate sample datasets when run directly
    generate_sample_datasets()
