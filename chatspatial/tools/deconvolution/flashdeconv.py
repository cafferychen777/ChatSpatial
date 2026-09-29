"""Spatial deconvolution with deterministic gene weighting and spatial regularization."""

from typing import Any, Literal

import pandas as pd

from ...utils.dependency_manager import require
from ...utils.exceptions import ChatSpatialError, ProcessingError
from .base import PreparedDeconvolutionData, create_deconvolution_stats


def deconvolve(
    data: PreparedDeconvolutionData,
    sketch_dim: int = 512,
    lambda_spatial: float | Literal["auto"] = "auto",
    n_hvg: int = 2000,
    n_markers_per_type: int = 50,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Deconvolve spatial data using FlashDeconv.

    Args:
        data: Prepared deconvolution data (immutable)
        sketch_dim: Bucket count used to compute deterministic gene weights
        lambda_spatial: Spatial regularization; "auto" adapts to the data scale, 0 disables it
        n_hvg: Number of highly variable genes to use (default: 2000)
        n_markers_per_type: Number of marker genes per cell type (default: 50)

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    fd = require("flashdeconv", feature="FlashDeconv deconvolution")

    try:
        # Data already copied in prepare_deconvolution
        adata_st = data.spatial
        reference = data.reference

        # Run FlashDeconv
        fd.tl.deconvolve(
            adata_st,
            reference,
            cell_type_key=data.cell_type_key,
            sketch_dim=sketch_dim,
            lambda_spatial=lambda_spatial,
            n_hvg=n_hvg,
            n_markers_per_type=n_markers_per_type,
        )

        # Extract proportions
        if "flashdeconv" not in adata_st.obsm:
            raise ProcessingError(
                "FlashDeconv did not produce output in adata.obsm['flashdeconv']"
            )

        proportions = adata_st.obsm["flashdeconv"].copy()

        params = adata_st.uns["flashdeconv_params"]

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            data.common_genes,
            method="FlashDeconv",
            device="CPU",
            backend_version=fd.__version__,
            genes_used=params["n_genes_used"],
            **{
                key: params[key]
                for key in (
                    "sketch_dim",
                    "lambda_spatial",
                    "n_hvg",
                    "n_markers_per_type",
                    "gene_weighting",
                    "max_iter",
                    "tol",
                    "preprocess",
                    "rho_sparsity",
                    "spatial_method",
                    "k_neighbors",
                    "converged",
                    "n_iterations",
                )
            },
        )

        return proportions, stats

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"FlashDeconv deconvolution failed: {e}") from e
