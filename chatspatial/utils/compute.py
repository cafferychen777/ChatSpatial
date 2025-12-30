"""
Compute utilities for ChatSpatial.

This module provides lazy computation functions that ensure required
computations are available before analysis. These functions follow the
"ensure" pattern: check if computation exists, compute if missing.

Design Principles:
1. Single Responsibility: Each function ensures one computation
2. Idempotent: Safe to call multiple times
3. Transparent: Returns whether computation was performed
4. Composable: Functions can depend on each other

Usage:
    # In analysis tools, use these to ensure prerequisites
    computed = ensure_pca(adata)

    # Or use the async version with context
    await ensure_pca_async(adata, ctx)
"""

from typing import TYPE_CHECKING, Literal, Optional

import scanpy as sc

from .adata_utils import ensure_categorical
from .exceptions import DataNotFoundError

if TYPE_CHECKING:
    import anndata as ad


# =============================================================================
# Core Computation Functions
# =============================================================================


def ensure_pca(
    adata: "ad.AnnData",
    n_comps: int = 30,
    use_highly_variable: bool = True,
    random_state: int = 0,
) -> bool:
    """
    Ensure PCA is computed on the dataset.

    Args:
        adata: AnnData object (modified in-place)
        n_comps: Number of principal components
        use_highly_variable: Use only HVG if available
        random_state: Random seed for reproducibility

    Returns:
        True if PCA was computed, False if already existed
    """
    if "X_pca" in adata.obsm:
        return False

    # Adjust n_comps if necessary
    max_comps = min(adata.n_obs, adata.n_vars) - 1
    n_comps = min(n_comps, max_comps)

    sc.tl.pca(
        adata,
        n_comps=n_comps,
        use_highly_variable=use_highly_variable and "highly_variable" in adata.var,
        random_state=random_state,
    )
    return True


def ensure_neighbors(
    adata: "ad.AnnData",
    n_neighbors: int = 15,
    n_pcs: Optional[int] = None,
    use_rep: str = "X_pca",
    random_state: int = 0,
) -> bool:
    """
    Ensure neighborhood graph is computed.

    Automatically ensures PCA is available first.

    Args:
        adata: AnnData object (modified in-place)
        n_neighbors: Number of neighbors for k-NN graph
        n_pcs: Number of PCs to use (None = auto)
        use_rep: Representation to use (default: X_pca)
        random_state: Random seed

    Returns:
        True if neighbors was computed, False if already existed
    """
    if "neighbors" in adata.uns and "connectivities" in adata.obsp:
        return False

    # Ensure PCA exists if using X_pca
    if use_rep == "X_pca":
        ensure_pca(adata)

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep,
        random_state=random_state,
    )
    return True


def ensure_umap(
    adata: "ad.AnnData",
    min_dist: float = 0.5,
    spread: float = 1.0,
    random_state: int = 0,
) -> bool:
    """
    Ensure UMAP embedding is computed.

    Automatically ensures neighbors are available first.

    Args:
        adata: AnnData object (modified in-place)
        min_dist: Minimum distance parameter for UMAP
        spread: Spread parameter for UMAP
        random_state: Random seed

    Returns:
        True if UMAP was computed, False if already existed
    """
    if "X_umap" in adata.obsm:
        return False

    ensure_neighbors(adata)

    sc.tl.umap(
        adata,
        min_dist=min_dist,
        spread=spread,
        random_state=random_state,
    )
    return True


def ensure_leiden(
    adata: "ad.AnnData",
    resolution: float = 1.0,
    key_added: str = "leiden",
    random_state: int = 0,
) -> bool:
    """
    Ensure Leiden clustering is computed.

    Automatically ensures neighbors are available first.

    Args:
        adata: AnnData object (modified in-place)
        resolution: Clustering resolution (higher = more clusters)
        key_added: Key for storing results in adata.obs
        random_state: Random seed

    Returns:
        True if clustering was computed, False if already existed
    """
    if key_added in adata.obs:
        return False

    ensure_neighbors(adata)

    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=key_added,
        random_state=random_state,
    )

    ensure_categorical(adata, key_added)
    return True


def ensure_louvain(
    adata: "ad.AnnData",
    resolution: float = 1.0,
    key_added: str = "louvain",
    random_state: int = 0,
) -> bool:
    """
    Ensure Louvain clustering is computed.

    Automatically ensures neighbors are available first.

    Args:
        adata: AnnData object (modified in-place)
        resolution: Clustering resolution
        key_added: Key for storing results in adata.obs
        random_state: Random seed

    Returns:
        True if clustering was computed, False if already existed
    """
    if key_added in adata.obs:
        return False

    ensure_neighbors(adata)

    sc.tl.louvain(
        adata,
        resolution=resolution,
        key_added=key_added,
        random_state=random_state,
    )

    ensure_categorical(adata, key_added)
    return True


def ensure_diffmap(
    adata: "ad.AnnData",
    n_comps: int = 15,
) -> bool:
    """
    Ensure diffusion map is computed (for trajectory analysis).

    Automatically ensures neighbors are available first.

    Args:
        adata: AnnData object (modified in-place)
        n_comps: Number of diffusion components

    Returns:
        True if diffmap was computed, False if already existed
    """
    if "X_diffmap" in adata.obsm:
        return False

    ensure_neighbors(adata)

    sc.tl.diffmap(adata, n_comps=n_comps)
    return True


def ensure_spatial_neighbors(
    adata: "ad.AnnData",
    coord_type: Literal["grid", "generic"] = "generic",
    n_neighs: int = 6,
    n_rings: int = 1,
    spatial_key: str = "spatial",
) -> bool:
    """
    Ensure spatial neighborhood graph is computed.

    Args:
        adata: AnnData object (modified in-place)
        coord_type: Type of coordinate system ('grid' for Visium, 'generic' for others)
        n_neighs: Number of neighbors (for generic coord_type)
        n_rings: Number of rings (for grid coord_type)
        spatial_key: Key for spatial coordinates in obsm

    Returns:
        True if spatial neighbors was computed, False if already existed
    """
    if "spatial_connectivities" in adata.obsp:
        return False

    if spatial_key not in adata.obsm:
        raise DataNotFoundError(
            f"Spatial coordinates not found in adata.obsm['{spatial_key}']"
        )

    import squidpy as sq

    if coord_type == "grid":
        sq.gr.spatial_neighbors(adata, coord_type="grid", n_rings=n_rings)
    else:
        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=n_neighs)

    return True


# =============================================================================
# Async Wrapper for Spatial Neighbors (used by spatial_statistics.py)
# =============================================================================


async def ensure_spatial_neighbors_async(
    adata: "ad.AnnData",
    ctx,  # ToolContext  # noqa: ARG001
    coord_type: Literal["grid", "generic"] = "generic",
    n_neighs: int = 6,
) -> bool:
    """Async version of ensure_spatial_neighbors with context logging."""
    return ensure_spatial_neighbors(adata, coord_type, n_neighs)


# =============================================================================
# Validation Functions (Check-only, no computation)
# =============================================================================


def has_pca(adata: "ad.AnnData") -> bool:
    """Check if PCA is available."""
    return "X_pca" in adata.obsm


def has_neighbors(adata: "ad.AnnData") -> bool:
    """Check if neighborhood graph is available."""
    return "neighbors" in adata.uns and "connectivities" in adata.obsp


def has_umap(adata: "ad.AnnData") -> bool:
    """Check if UMAP embedding is available."""
    return "X_umap" in adata.obsm


def has_clustering(adata: "ad.AnnData", key: str = "leiden") -> bool:
    """Check if clustering results are available."""
    return key in adata.obs


def has_spatial_neighbors(adata: "ad.AnnData") -> bool:
    """Check if spatial neighborhood graph is available."""
    return "spatial_connectivities" in adata.obsp


def has_hvg(adata: "ad.AnnData") -> bool:
    """Check if highly variable genes are marked."""
    return "highly_variable" in adata.var and adata.var["highly_variable"].any()
