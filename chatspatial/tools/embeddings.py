"""
Embedding computation tools for spatial transcriptomics data.

This module provides explicit control over dimensionality reduction and clustering
computations. While analysis tools compute these lazily using ensure_* functions,
users can use this tool to control computation parameters directly.
"""

from typing import List, Literal, Optional

from pydantic import BaseModel, Field

from ..spatial_mcp_adapter import ToolContext
from ..utils.compute import (
    ensure_diffmap,
    ensure_leiden,
    ensure_louvain,
    ensure_neighbors,
    ensure_pca,
    ensure_spatial_neighbors,
    ensure_umap,
)
from ..utils.mcp_utils import mcp_tool_error_handler


class EmbeddingParameters(BaseModel):
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


class EmbeddingResult(BaseModel):
    """Result of embedding computation."""

    data_id: str
    computed: List[str]
    skipped: List[str]
    n_clusters: Optional[int] = None
    pca_variance_ratio: Optional[float] = None


@mcp_tool_error_handler()
async def compute_embeddings(
    data_id: str,
    ctx: ToolContext,
    params: EmbeddingParameters = EmbeddingParameters(),
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
    adata = await ctx.get_adata(data_id)
    computed = []
    skipped = []

    # Handle force recomputation by removing existing results
    if params.force:
        if params.compute_pca and "X_pca" in adata.obsm:
            del adata.obsm["X_pca"]
            if "pca" in adata.uns:
                del adata.uns["pca"]
        if params.compute_neighbors:
            if "neighbors" in adata.uns:
                del adata.uns["neighbors"]
            if "connectivities" in adata.obsp:
                del adata.obsp["connectivities"]
            if "distances" in adata.obsp:
                del adata.obsp["distances"]
        if params.compute_umap and "X_umap" in adata.obsm:
            del adata.obsm["X_umap"]
        if params.compute_clustering and params.clustering_key in adata.obs:
            del adata.obs[params.clustering_key]
        if params.compute_diffmap and "X_diffmap" in adata.obsm:
            del adata.obsm["X_diffmap"]
        if params.compute_spatial_neighbors and "spatial_connectivities" in adata.obsp:
            del adata.obsp["spatial_connectivities"]
            if "spatial_distances" in adata.obsp:
                del adata.obsp["spatial_distances"]

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

    # 2. Neighbors (requires PCA)
    if params.compute_neighbors:
        if ensure_neighbors(
            adata,
            n_neighbors=params.n_neighbors,
            n_pcs=params.n_pcs,
            random_state=params.random_state,
        ):
            computed.append("neighbors")
        else:
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
        if params.clustering_method == "leiden":
            if ensure_leiden(
                adata,
                resolution=params.clustering_resolution,
                key_added=params.clustering_key,
                random_state=params.random_state,
            ):
                n_clusters = adata.obs[params.clustering_key].nunique()
                computed.append(f"Leiden clustering ({n_clusters} clusters)")
            else:
                skipped.append(f"{params.clustering_key} (already exists)")
                n_clusters = adata.obs[params.clustering_key].nunique()
        else:
            if ensure_louvain(
                adata,
                resolution=params.clustering_resolution,
                key_added=params.clustering_key,
                random_state=params.random_state,
            ):
                n_clusters = adata.obs[params.clustering_key].nunique()
                computed.append(f"Louvain clustering ({n_clusters} clusters)")
            else:
                skipped.append(f"{params.clustering_key} (already exists)")
                n_clusters = adata.obs[params.clustering_key].nunique()

    # 5. Diffusion map (requires neighbors)
    if params.compute_diffmap:
        if ensure_diffmap(adata, n_comps=params.diffmap_n_comps):
            computed.append("diffusion map")
        else:
            skipped.append("diffusion map (already exists)")

    # 6. Spatial neighbors
    if params.compute_spatial_neighbors:
        try:
            if ensure_spatial_neighbors(
                adata,
                coord_type=params.spatial_coord_type,
                n_neighs=params.spatial_n_neighs,
                n_rings=params.spatial_n_rings,
            ):
                computed.append("spatial neighbors")
            else:
                skipped.append("spatial neighbors (already exists)")
        except ValueError as e:
            await ctx.warning(f"Could not compute spatial neighbors: {e}")
            skipped.append(f"spatial neighbors (error: {e})")

    # Get PCA variance ratio if available
    pca_variance_ratio = None
    if "pca" in adata.uns and "variance_ratio" in adata.uns["pca"]:
        pca_variance_ratio = float(adata.uns["pca"]["variance_ratio"].sum())

    # Store updated data
    await ctx.set_adata(data_id, adata)

    return EmbeddingResult(
        data_id=data_id,
        computed=computed,
        skipped=skipped,
        n_clusters=n_clusters,
        pca_variance_ratio=pca_variance_ratio,
    )
