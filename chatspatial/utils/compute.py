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

import inspect
from typing import TYPE_CHECKING, Literal, Optional

import numpy as np
import scanpy as sc
import scipy.sparse

from .adata_utils import ensure_categorical
from .dependency_manager import require
from .exceptions import DataError, DataNotFoundError

if TYPE_CHECKING:
    import anndata as ad


# =============================================================================
# Core Computation Functions
# =============================================================================


def select_hvgs_by_variance(adata: "ad.AnnData", n_hvgs: int) -> None:
    """Mark the genes with the highest finite variance as highly variable.

    Variance is defined for any finite matrix, which makes this usable when
    scanpy's mean-binned dispersion is not.
    """
    if scipy.sparse.issparse(adata.X):
        means = np.asarray(adata.X.mean(axis=0)).ravel()
        sq_means = np.asarray(adata.X.power(2).mean(axis=0)).ravel()
        var = sq_means - np.square(means)
    else:
        var = np.var(np.asarray(adata.X), axis=0)
    var = np.nan_to_num(var, nan=-np.inf, posinf=-np.inf, neginf=-np.inf)
    n_select = min(max(int(n_hvgs), 1), adata.n_vars)
    top_idx = np.argpartition(var, -n_select)[-n_select:]
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[top_idx] = True
    adata.var["highly_variable"] = mask


def ensure_highly_variable_genes(
    adata: "ad.AnnData",
    n_top_genes: int,
    **scanpy_kwargs: object,
) -> bool:
    """Mark highly variable genes, falling back to variance ranking.

    Scanpy bins genes by mean expression to estimate dispersion. Degenerate
    inputs leave those bins unusable — genes that are all zero after merging
    samples with different gene sets, or bins that collapse to one value — and
    scanpy then fails in whichever way it reaches first: a ``KeyError`` for NaN
    bin edges, a ``ValueError`` for infinite ones, an ``IndexError`` when every
    normalized dispersion came out NaN and none is left to rank. All mean the
    same thing, so the fallback is chosen on the failure itself rather than on
    its wording.

    ``flavor="pearson_residuals"`` is routed to scanpy's experimental
    implementation, the only one that accepts it. Dispatching on the flavor
    before the call matters: the stable implementation rejects the name with a
    ``ValueError``, which the fallback below would otherwise swallow and
    silently answer with plain variance ranking.

    Args:
        adata: Dataset to annotate in place.
        n_top_genes: Number of genes to mark.
        **scanpy_kwargs: Extra arguments forwarded to scanpy.

    Returns:
        True when the variance fallback was used instead of scanpy's binning.
    """
    select_genes = (
        sc.experimental.pp.highly_variable_genes
        if scanpy_kwargs.get("flavor") == "pearson_residuals"
        else sc.pp.highly_variable_genes
    )
    try:
        select_genes(adata, n_top_genes=n_top_genes, **scanpy_kwargs)
        return False
    except (KeyError, ValueError, IndexError):
        select_hvgs_by_variance(adata, n_top_genes)
        # The fallback must leave the dataset in the shape scanpy would have,
        # or callers get a marked-but-unfiltered matrix on degenerate input.
        if scanpy_kwargs.get("subset"):
            adata._inplace_subset_var(adata.var["highly_variable"].to_numpy())
        return True


def top_n_desc_indices(
    values: np.ndarray,
    n_top: int,
    *,
    sanitize_nonfinite: bool = False,
) -> np.ndarray:
    """Return indices of top-n values in descending order.

    Args:
        values: 1D score array.
        n_top: Number of top indices to return.
        sanitize_nonfinite: Replace non-finite values with -inf before ranking.

    Returns:
        Indices sorted from highest score to lowest score.
    """
    if n_top <= 0 or values.size == 0:
        return np.array([], dtype=int)

    rank_values = np.asarray(values, dtype=float)
    if sanitize_nonfinite:
        rank_values = np.where(np.isfinite(rank_values), rank_values, -np.inf)

    n = min(n_top, rank_values.size)
    top_idx = np.argpartition(rank_values, -n)[-n:]
    return top_idx[np.argsort(rank_values[top_idx])[::-1]]


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

    if n_comps < 1:
        raise DataError(
            f"Cannot compute PCA: dataset has {adata.n_obs} cells and "
            f"{adata.n_vars} genes (need at least 2 of each)."
        )

    pca_kwargs: dict[str, object] = {
        "n_comps": n_comps,
        "random_state": random_state,
    }
    use_hvg_mask = use_highly_variable and "highly_variable" in adata.var
    if "mask_var" in inspect.signature(sc.tl.pca).parameters:
        pca_kwargs["mask_var"] = "highly_variable" if use_hvg_mask else None
    else:
        pca_kwargs["use_highly_variable"] = use_hvg_mask
    sc.tl.pca(adata, **pca_kwargs)
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
        ensure_pca(
            adata,
            n_comps=n_pcs if n_pcs is not None else 30,
            random_state=random_state,
        )

    # PCA may contain fewer dimensions than requested for small datasets.
    # Clamp to the representation that was actually computed instead of
    # passing an impossible value to Scanpy's representation selector.
    if n_pcs is not None and use_rep in adata.obsm:
        n_pcs = min(n_pcs, adata.obsm[use_rep].shape[1])

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

    # scanpy imports the louvain package lazily, so without this the caller
    # sees a bare "No module named 'louvain'" with no way to act on it.
    require("louvain", feature="Louvain clustering")

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
        sq.gr.spatial_neighbors(
            adata, coord_type="grid", n_rings=n_rings, spatial_key=spatial_key
        )
    else:
        sq.gr.spatial_neighbors(
            adata, coord_type="generic", n_neighs=n_neighs, spatial_key=spatial_key
        )

    return True


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


# =============================================================================
# Gaussian Mixture Model Clustering (replaces R mclust dependency)
# =============================================================================


def gmm_clustering(
    data: "np.ndarray",
    n_clusters: int,
    covariance_type: str = "tied",
    random_state: int = 42,
    n_init: int = 10,
    max_iter: int = 300,
) -> "np.ndarray":
    """
    Gaussian Mixture Model clustering using sklearn.

    This is a pure Python implementation that replaces R's mclust package.
    The 'tied' covariance_type is equivalent to mclust's EEE model
    (Equal volume, Equal shape, Equal orientation).

    Covariance type mapping (mclust -> sklearn):
        - EEE -> 'tied' (all clusters share same covariance)
        - VVV -> 'full' (each cluster has its own covariance)
        - EII -> 'spherical' (spherical, same across clusters)
        - VII -> 'diag' (diagonal, different across clusters)

    Why this replaces mclust:
        1. Eliminates R dependency (no rpy2 required)
        2. Faster execution (no R interop overhead)
        3. Better error handling (Python native exceptions)
        4. Tested equivalent: ARI = 1.0 with mclust EEE

    Args:
        data: Input data matrix (n_samples, n_features)
        n_clusters: Target number of clusters (equivalent to mclust's G parameter)
        covariance_type: GMM covariance type
            - 'tied': Same covariance for all clusters (EEE model, default)
            - 'full': Each cluster has its own general covariance
            - 'diag': Each cluster has diagonal covariance
            - 'spherical': Each cluster has single variance
        random_state: Random seed for reproducibility
        n_init: Number of initializations (EM is sensitive to initialization)
        max_iter: Maximum EM iterations per initialization

    Returns:
        Cluster labels as integer array (1-indexed like mclust for compatibility)

    Example:
        >>> labels = gmm_clustering(embeddings, n_clusters=7)
        >>> adata.obs['domain'] = labels
    """
    import numpy as np
    from sklearn.mixture import GaussianMixture

    # Validate input
    if data.ndim != 2:
        raise ValueError(f"Expected 2D array, got {data.ndim}D")

    if n_clusters < 1:
        raise ValueError(f"n_clusters must be >= 1, got {n_clusters}")

    if n_clusters > data.shape[0]:
        raise ValueError(
            f"n_clusters ({n_clusters}) cannot exceed n_samples ({data.shape[0]})"
        )

    # Ensure float64 for numerical stability (same as mclust requirement)
    data = np.asarray(data, dtype=np.float64)

    # Initialize and fit GMM
    gmm = GaussianMixture(
        n_components=n_clusters,
        covariance_type=covariance_type,
        random_state=random_state,
        n_init=n_init,
        max_iter=max_iter,
        init_params="k-means++",  # Better initialization than random
    )

    # Fit and predict
    labels = gmm.fit_predict(data)

    if not gmm.converged_:
        import warnings

        warnings.warn(
            f"GMM did not converge after {gmm.n_iter_} iterations. "
            "Results may be unreliable. Consider increasing max_iter.",
            UserWarning,
            stacklevel=2,
        )

    # Convert to 1-indexed labels (mclust compatibility)
    # mclust returns labels starting from 1, sklearn from 0
    labels = labels + 1

    return labels
