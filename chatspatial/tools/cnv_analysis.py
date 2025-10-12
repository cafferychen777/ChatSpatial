"""
Copy Number Variation (CNV) analysis tools for spatial transcriptomics data.
"""

from typing import Any, Dict, Optional

import numpy as np
import scanpy as sc
from mcp.server.fastmcp import Context

from ..models.analysis import CNVResult
from ..models.data import CNVParameters

# Dependency checking for infercnvpy
try:
    import infercnvpy as cnv

    INFERCNVPY_AVAILABLE = True
except ImportError:
    INFERCNVPY_AVAILABLE = False


async def infer_cnv(
    data_id: str,
    data_store: Dict[str, Any],
    params: CNVParameters,
    context: Optional[Context] = None,
) -> CNVResult:
    """Infer copy number variations using infercnvpy

    This function performs CNV inference on spatial transcriptomics data using
    infercnvpy, which detects chromosomal copy number alterations by comparing
    gene expression patterns across chromosomes between tumor and normal cells.

    Args:
        data_id: Dataset identifier
        data_store: Dictionary storing loaded datasets
        params: CNV analysis parameters (reference_key, reference_categories, etc.)
        context: MCP context for logging

    Returns:
        CNVResult containing CNV analysis results and statistics

    Raises:
        ValueError: If dataset not found or parameters are invalid
        RuntimeError: If infercnvpy is not available or analysis fails
    """
    # Check if infercnvpy is available
    if not INFERCNVPY_AVAILABLE:
        raise RuntimeError(
            "infercnvpy is not installed. Please install it with:\n"
            "  pip install 'chatspatial[cnv]'\n"
            "or:\n"
            "  pip install infercnvpy>=0.4.0"
        )

    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise ValueError(f"Dataset '{data_id}' not found in data store")

    adata = data_store[data_id]["adata"]

    # Validate reference_key exists
    if params.reference_key not in adata.obs.columns:
        available_keys = ", ".join(adata.obs.columns[:10])
        raise ValueError(
            f"Reference key '{params.reference_key}' not found in adata.obs.\n"
            f"Available columns: {available_keys}..."
        )

    # Validate reference_categories exist in the reference_key column
    available_categories = set(adata.obs[params.reference_key].unique())
    missing_categories = set(params.reference_categories) - available_categories
    if missing_categories:
        raise ValueError(
            f"Reference categories {missing_categories} not found in "
            f"adata.obs['{params.reference_key}'].\n"
            f"Available categories: {sorted(available_categories)}"
        )

    if context:
        await context.info(
            f"Running CNV inference with {len(params.reference_categories)} "
            f"reference cell types: {', '.join(params.reference_categories)}"
        )

    # Create a copy of adata for CNV analysis
    adata_cnv = adata.copy()

    # Check if gene position information is available
    if "chromosome" not in adata_cnv.var.columns:
        if context:
            await context.warning(
                "No chromosome information found in adata.var. "
                "Attempting to infer from gene names..."
            )
        try:
            # Try to infer gene positions from infercnvpy's built-in database
            cnv.tl.infercnv(
                adata_cnv,
                reference_key=params.reference_key,
                reference_cat=params.reference_categories,
                window_size=params.window_size,
                step=params.step,
                dynamic_threshold=params.dynamic_threshold,
            )
        except Exception as e:
            raise RuntimeError(
                "Failed to run CNV inference. Gene position information is "
                "required but not found in adata.var.\n"
                "Please ensure your data includes chromosome and position "
                "information, or use infercnvpy's built-in gene annotation.\n"
                f"Error: {str(e)}"
            )
    else:
        # Gene positions are available, run CNV inference
        if context:
            await context.info("Running infercnvpy CNV inference...")

        # Exclude chromosomes if specified
        if params.exclude_chromosomes:
            genes_to_keep = ~adata_cnv.var["chromosome"].isin(
                params.exclude_chromosomes
            )
            adata_cnv = adata_cnv[:, genes_to_keep].copy()
            if context:
                await context.info(
                    f"Excluded chromosomes: {', '.join(params.exclude_chromosomes)}"
                )

        # Run infercnvpy
        cnv.tl.infercnv(
            adata_cnv,
            reference_key=params.reference_key,
            reference_cat=params.reference_categories,
            window_size=params.window_size,
            step=params.step,
            dynamic_threshold=params.dynamic_threshold,
        )

    # Optional: Cluster cells by CNV pattern
    if params.cluster_cells:
        if context:
            await context.info("Clustering cells by CNV pattern...")
        try:
            sc.pp.neighbors(adata_cnv, use_rep="X_cnv", n_neighbors=15)
            sc.tl.leiden(adata_cnv, key_added="cnv_clusters")
            if context:
                n_clusters = len(adata_cnv.obs["cnv_clusters"].unique())
                await context.info(f"Identified {n_clusters} CNV-based clusters")
        except Exception as e:
            if context:
                await context.warning(f"Failed to cluster cells by CNV: {str(e)}")

    # Optional: Compute dendrogram
    if params.dendrogram and params.cluster_cells:
        if context:
            await context.info("Computing hierarchical clustering dendrogram...")
        try:
            sc.tl.dendrogram(adata_cnv, groupby="cnv_clusters")
        except Exception as e:
            if context:
                await context.warning(f"Failed to compute dendrogram: {str(e)}")

    # Extract CNV statistics
    if context:
        await context.info("Extracting CNV statistics...")

    # Check what data is available
    cnv_score_key = None
    if "X_cnv" in adata_cnv.obsm:
        cnv_score_key = "X_cnv"
    elif "cnv" in adata_cnv.layers:
        cnv_score_key = "cnv"

    # Calculate statistics
    statistics = {}
    if cnv_score_key and cnv_score_key in adata_cnv.obsm:
        cnv_matrix = adata_cnv.obsm[cnv_score_key]

        # Convert sparse matrix to dense if needed
        import scipy.sparse
        if scipy.sparse.issparse(cnv_matrix):
            cnv_matrix = cnv_matrix.toarray()

        statistics["mean_cnv"] = float(np.mean(cnv_matrix))
        statistics["std_cnv"] = float(np.std(cnv_matrix))
        statistics["median_cnv"] = float(np.median(cnv_matrix))

        # Calculate per-cell CNV scores
        cell_cnv_scores = np.mean(np.abs(cnv_matrix), axis=1)
        statistics["mean_cell_cnv_score"] = float(np.mean(cell_cnv_scores))
        statistics["max_cell_cnv_score"] = float(np.max(cell_cnv_scores))

    # Count reference vs non-reference cells
    is_reference = adata_cnv.obs[params.reference_key].isin(params.reference_categories)
    statistics["n_reference_cells"] = int(is_reference.sum())
    statistics["n_non_reference_cells"] = int((~is_reference).sum())

    # Get chromosome information
    if "chromosome" in adata_cnv.var.columns:
        n_chromosomes = len(adata_cnv.var["chromosome"].unique())
    else:
        n_chromosomes = 0  # Unknown

    n_genes_analyzed = adata_cnv.n_vars

    # Store CNV results back in the original adata object
    if cnv_score_key and cnv_score_key in adata_cnv.obsm:
        adata.obsm[cnv_score_key] = adata_cnv.obsm[cnv_score_key]

    # Store CNV metadata (required for infercnvpy plotting functions)
    if "cnv" in adata_cnv.uns:
        adata.uns["cnv"] = adata_cnv.uns["cnv"]

    if params.cluster_cells and "cnv_clusters" in adata_cnv.obs:
        adata.obs["cnv_clusters"] = adata_cnv.obs["cnv_clusters"]

    if params.dendrogram and "dendrogram_cnv_clusters" in adata_cnv.uns:
        adata.uns["dendrogram_cnv_clusters"] = adata_cnv.uns["dendrogram_cnv_clusters"]

    # Store CNV analysis parameters in adata.uns for reference
    adata.uns["cnv_analysis"] = {
        "reference_key": params.reference_key,
        "reference_categories": list(params.reference_categories),  # Convert to list
        "window_size": params.window_size,
        "step": params.step,
        "cnv_score_key": cnv_score_key,
    }

    if context:
        await context.info(
            f"CNV analysis complete. Analyzed {n_genes_analyzed} genes "
            f"across {n_chromosomes} chromosomes."
        )

    return CNVResult(
        data_id=data_id,
        method="infercnvpy",
        reference_key=params.reference_key,
        reference_categories=list(params.reference_categories),  # Convert to list
        n_chromosomes=n_chromosomes,
        n_genes_analyzed=n_genes_analyzed,
        cnv_score_key=cnv_score_key,
        statistics=statistics,
        visualization_available=cnv_score_key is not None,
    )
