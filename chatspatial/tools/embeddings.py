"""
Embedding computation tools for spatial transcriptomics data.

This module provides explicit control over dimensionality reduction and clustering
computations. While analysis tools compute these lazily using ensure_* functions,
users can use this tool to control computation parameters directly.
"""

from typing import Literal

from pydantic import Field

from ..models.analysis import EmbeddingResult
from ..models.base import StrictParameters
from ..spatial_mcp_adapter import ToolContext
from ..utils.adata_utils import get_spatial_key, store_analysis_metadata
from ..utils.compute import (
    ensure_diffmap,
    ensure_leiden,
    ensure_louvain,
    ensure_neighbors,
    ensure_pca,
    ensure_spatial_neighbors,
    ensure_umap,
)
from ..utils.exceptions import DataNotFoundError
from ..utils.results_export import export_analysis_result


class EmbeddingParameters(StrictParameters):
    """Parameters for embedding computation."""

    # What to compute
    compute_pca: bool = Field(
        default=True,
        description="Compute PCA dimensionality reduction",
    )
    compute_neighbors: bool = Field(
        default=True,
        description="Compute k-NN neighbor graph (requires PCA)",
    )
    compute_umap: bool = Field(
        default=True,
        description="Compute UMAP embedding (requires neighbors)",
    )
    compute_clustering: bool = Field(
        default=True,
        description="Compute Leiden clustering (requires neighbors)",
    )
    compute_diffmap: bool = Field(
        default=False,
        description="Compute diffusion map for trajectory analysis (requires neighbors)",
    )
    compute_spatial_neighbors: bool = Field(
        default=True,
        description="Compute spatial neighborhood graph for spatial analysis",
    )

    # PCA parameters
    n_pcs: int = Field(
        default=30,
        ge=2,
        le=100,
        description="Number of principal components",
    )
    use_highly_variable: bool = Field(
        default=True,
        description="Use only highly variable genes for PCA",
    )

    # Neighbor graph parameters
    n_neighbors: int = Field(
        default=15,
        ge=2,
        le=100,
        description="Number of neighbors for k-NN graph",
    )

    # UMAP parameters
    umap_min_dist: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="UMAP minimum distance parameter",
    )

    # Clustering parameters
    clustering_method: Literal["leiden", "louvain"] = Field(
        default="leiden",
        description="Clustering algorithm",
    )
    clustering_resolution: float = Field(
        default=1.0,
        ge=0.1,
        le=2.0,
        description="Clustering resolution (higher = more clusters)",
    )
    clustering_key: str = Field(
        default="leiden",
        description="Key to store clustering results in adata.obs",
    )

    # Diffusion map parameters
    diffmap_n_comps: int = Field(
        default=15,
        ge=2,
        le=50,
        description="Number of diffusion components",
    )

    # Spatial neighbor parameters
    spatial_coord_type: Literal["grid", "generic"] = Field(
        default="generic",
        description="Coordinate type: 'grid' for Visium hexagonal, 'generic' for others",
    )
    spatial_n_neighs: int = Field(
        default=6,
        ge=1,
        le=30,
        description="Number of spatial neighbors (for generic coord_type)",
    )
    spatial_n_rings: int = Field(
        default=1,
        ge=1,
        le=3,
        description="Number of rings (for grid coord_type)",
    )

    # Force recomputation
    force: bool = Field(
        default=False,
        description="Force recomputation even if results already exist",
    )

    # Random seed
    random_state: int = Field(
        default=0,
        description="Random seed for reproducibility",
    )


_SPATIAL_OBSP_KEYS = ("spatial_connectivities", "spatial_distances")


def _clear_requested_embeddings(
    adata,
    params: EmbeddingParameters,
    *,
    needs_neighbors: bool,
) -> None:
    """Remove only outputs selected for forced recomputation."""
    if (params.compute_pca or needs_neighbors) and "X_pca" in adata.obsm:
        del adata.obsm["X_pca"]
        adata.uns.pop("pca", None)
    if needs_neighbors:
        adata.uns.pop("neighbors", None)
        adata.obsp.pop("connectivities", None)
        adata.obsp.pop("distances", None)
    if params.compute_umap:
        adata.obsm.pop("X_umap", None)
    if params.compute_clustering and params.clustering_key in adata.obs:
        adata.obs.pop(params.clustering_key)
    if params.compute_diffmap:
        adata.obsm.pop("X_diffmap", None)
    if params.compute_spatial_neighbors:
        for key in _SPATIAL_OBSP_KEYS:
            adata.obsp.pop(key, None)
        adata.uns.pop("spatial_neighbors", None)


def _compute_clustering(
    adata,
    params: EmbeddingParameters,
    computed: list[str],
    skipped: list[str],
) -> int | None:
    """Compute the requested cluster labels and return their count."""
    ensure_cluster = (
        ensure_leiden if params.clustering_method == "leiden" else ensure_louvain
    )
    created = ensure_cluster(
        adata,
        resolution=params.clustering_resolution,
        key_added=params.clustering_key,
        random_state=params.random_state,
    )
    if params.clustering_key not in adata.obs:
        skipped.append(f"{params.clustering_key} (missing; clustering not computed)")
        return None

    n_clusters = int(adata.obs[params.clustering_key].nunique())
    if created:
        method_name = params.clustering_method.capitalize()
        computed.append(f"{method_name} clustering ({n_clusters} clusters)")
    else:
        skipped.append(f"{params.clustering_key} (already exists)")
    return n_clusters


async def _compute_spatial_neighbors(
    adata,
    params: EmbeddingParameters,
    ctx: ToolContext,
    computed: list[str],
    skipped: list[str],
    previous_obsp: dict,
    previous_metadata,
    had_metadata: bool,
) -> None:
    """Compute the optional spatial graph transactionally."""
    try:
        detected_key = get_spatial_key(adata) or "spatial"
        created = ensure_spatial_neighbors(
            adata,
            coord_type=params.spatial_coord_type,
            n_neighs=params.spatial_n_neighs,
            n_rings=params.spatial_n_rings,
            spatial_key=detected_key,
        )
    except (DataNotFoundError, ValueError) as exc:
        for key in _SPATIAL_OBSP_KEYS:
            adata.obsp.pop(key, None)
        adata.obsp.update(previous_obsp)
        if had_metadata:
            adata.uns["spatial_neighbors"] = previous_metadata
        else:
            adata.uns.pop("spatial_neighbors", None)
        await ctx.warning(f"Could not compute spatial neighbors: {exc}")
        skipped.append(f"spatial neighbors (error: {exc})")
        return

    target = computed if created else skipped
    target.append(
        "spatial neighbors" if created else "spatial neighbors (already exists)"
    )


def _store_embedding_metadata(
    adata,
    data_id: str,
    params: EmbeddingParameters,
    *,
    computed: list[str],
    skipped: list[str],
    n_clusters: int | None,
    pca_variance_ratio: float | None,
) -> None:
    """Store and export clustering provenance when labels are available."""
    if not params.compute_clustering or params.clustering_key not in adata.obs:
        return

    actual_n_pcs = params.n_pcs
    actual_n_neighbors = params.n_neighbors
    pca_params = adata.uns.get("pca", {}).get("params", {})
    if hasattr(pca_params, "get"):
        actual_n_pcs = pca_params.get("n_comps", actual_n_pcs)
    neighbor_params = adata.uns.get("neighbors", {}).get("params", {})
    if hasattr(neighbor_params, "get"):
        actual_n_neighbors = neighbor_params.get("n_neighbors", actual_n_neighbors)

    analysis_name = f"embeddings_{params.clustering_method}"
    store_analysis_metadata(
        adata,
        analysis_name=analysis_name,
        method=params.clustering_method,
        parameters={
            "n_pcs": actual_n_pcs,
            "n_neighbors": actual_n_neighbors,
            "clustering_resolution": params.clustering_resolution,
            "clustering_key": params.clustering_key,
        },
        results_keys={"obs": [params.clustering_key]},
        statistics={
            "n_clusters": n_clusters,
            "pca_variance_ratio": pca_variance_ratio,
            "computed": computed,
            "skipped": skipped,
        },
    )
    export_analysis_result(adata, data_id, analysis_name)


async def compute_embeddings(
    data_id: str,
    ctx: ToolContext,
    params: EmbeddingParameters | None = None,
) -> EmbeddingResult:
    """Compute dimensionality reduction, clustering, and neighbor graphs.

    This tool provides explicit control over embedding computations.
    Analysis tools compute these lazily, but you can use this tool to:
    - Control computation parameters (n_pcs, n_neighbors, resolution)
    - Force recomputation with different parameters
    - Compute specific embeddings without running full preprocessing

    Args:
        data_id: Dataset ID
        ctx: Tool context
        params: Embedding computation parameters

    Returns:
        Summary of computed embeddings
    """
    if params is None:
        params = EmbeddingParameters()

    source_adata = await ctx.get_adata(data_id)
    adata = source_adata.copy()
    computed = []
    skipped = []

    # Spatial-neighbor computation is intentionally optional. Preserve its
    # previous complete state so a backend that writes some keys before raising
    # cannot turn a non-fatal warning into committed graph corruption.
    previous_spatial_obsp = {
        key: adata.obsp[key] for key in _SPATIAL_OBSP_KEYS if key in adata.obsp
    }
    had_spatial_neighbors_metadata = "spatial_neighbors" in adata.uns
    previous_spatial_neighbors_metadata = adata.uns.get("spatial_neighbors")

    # UMAP, clustering, and diffusion maps all depend on the neighbor graph,
    # which in turn depends on PCA.
    needs_neighbors = any(
        (
            params.compute_neighbors,
            params.compute_umap,
            params.compute_clustering,
            params.compute_diffmap,
        )
    )

    # Handle force recomputation by removing existing results
    if params.force:
        _clear_requested_embeddings(adata, params, needs_neighbors=needs_neighbors)

    # 1. PCA
    if params.compute_pca:
        if ensure_pca(
            adata,
            n_comps=params.n_pcs,
            use_highly_variable=params.use_highly_variable,
            random_state=params.random_state,
        ):
            computed.append("PCA")
        else:
            skipped.append("PCA (already exists)")

    # 2. Neighbors (required by UMAP, clustering, and diffusion maps)
    if needs_neighbors:
        had_pca = "X_pca" in adata.obsm
        neighbors_computed = ensure_neighbors(
            adata,
            n_neighbors=params.n_neighbors,
            n_pcs=params.n_pcs,
            random_state=params.random_state,
        )
        if not params.compute_pca and not had_pca and "X_pca" in adata.obsm:
            computed.append("PCA (dependency)")
        if neighbors_computed:
            computed.append(
                "neighbors" if params.compute_neighbors else "neighbors (dependency)"
            )
        elif params.compute_neighbors:
            skipped.append("neighbors (already exists)")

    # 3. UMAP (requires neighbors)
    if params.compute_umap:
        if ensure_umap(
            adata,
            min_dist=params.umap_min_dist,
            random_state=params.random_state,
        ):
            computed.append("UMAP")
        else:
            skipped.append("UMAP (already exists)")

    # 4. Clustering (requires neighbors)
    n_clusters = None
    if params.compute_clustering:
        n_clusters = _compute_clustering(adata, params, computed, skipped)

    # 5. Diffusion map (requires neighbors)
    if params.compute_diffmap:
        if ensure_diffmap(adata, n_comps=params.diffmap_n_comps):
            computed.append("diffusion map")
        else:
            skipped.append("diffusion map (already exists)")

    # 6. Spatial neighbors
    if params.compute_spatial_neighbors:
        await _compute_spatial_neighbors(
            adata,
            params,
            ctx,
            computed,
            skipped,
            previous_spatial_obsp,
            previous_spatial_neighbors_metadata,
            had_spatial_neighbors_metadata,
        )

    # Get PCA variance ratio if available
    pca_variance_ratio = None
    if "pca" in adata.uns and "variance_ratio" in adata.uns["pca"]:
        pca_variance_ratio = float(adata.uns["pca"]["variance_ratio"].sum())

    # Store metadata and export results (only if clustering was computed)
    # Note: Only clustering results are exported as CSV - PCA/UMAP coordinates
    # are too large for CSV export and are better accessed via adata directly
    _store_embedding_metadata(
        adata,
        data_id,
        params,
        computed=computed,
        skipped=skipped,
        n_clusters=n_clusters,
        pca_variance_ratio=pca_variance_ratio,
    )

    result = EmbeddingResult(
        data_id=data_id,
        computed=computed,
        skipped=skipped,
        n_clusters=n_clusters,
        pca_variance_ratio=pca_variance_ratio,
    )
    await ctx.set_adata(data_id, adata)
    return result
