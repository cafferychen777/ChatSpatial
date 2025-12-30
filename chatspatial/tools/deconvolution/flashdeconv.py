"""
FlashDeconv deconvolution method.

FlashDeconv is an ultra-fast spatial transcriptomics deconvolution method
that uses random sketching for O(N) time complexity.
"""

from typing import TYPE_CHECKING, Any, Dict, Tuple

import pandas as pd

if TYPE_CHECKING:
    pass

from ...utils.dependency_manager import is_available
from ...utils.exceptions import DependencyError, ProcessingError
from .base import DeconvolutionContext, create_deconvolution_stats


async def deconvolve(
    deconv_ctx: DeconvolutionContext,
    sketch_dim: int = 512,
    lambda_spatial: float = 5000.0,
    n_hvg: int = 2000,
    n_markers_per_type: int = 50,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using FlashDeconv.

    FlashDeconv is an ultra-fast deconvolution method with:
    - O(N) time complexity via random sketching
    - Processes 1M spots in ~3 minutes on CPU
    - No GPU required
    - Automatic marker gene selection
    - Spatial regularization for smooth proportions

    Args:
        deconv_ctx: Prepared DeconvolutionContext with validated data
        sketch_dim: Dimension for random sketching (default: 512)
        lambda_spatial: Spatial regularization strength (default: 5000.0)
        n_hvg: Number of highly variable genes to use (default: 2000)
        n_markers_per_type: Number of marker genes per cell type (default: 50)

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)

    Raises:
        DependencyError: If flashdeconv is not installed
        ProcessingError: If deconvolution fails
    """
    if not is_available("flashdeconv"):
        raise DependencyError(
            "FlashDeconv is not available. Install with: pip install flashdeconv"
        )

    try:
        import flashdeconv as fd

        # Create a copy for FlashDeconv (it modifies in place)
        adata_st = deconv_ctx.spatial_prepared.copy()
        reference = deconv_ctx.reference_prepared

        # Run FlashDeconv
        fd.tl.deconvolve(
            adata_st,
            reference,
            cell_type_key=deconv_ctx.cell_type_key,
            sketch_dim=sketch_dim,
            lambda_spatial=lambda_spatial,
            n_hvg=n_hvg,
            n_markers_per_type=n_markers_per_type,
        )

        # Extract proportions
        if "flashdeconv" not in adata_st.obsm:
            raise ProcessingError(
                "FlashDeconv did not produce expected output in adata.obsm['flashdeconv']"
            )

        proportions = adata_st.obsm["flashdeconv"].copy()

        # Ensure DataFrame format
        if not isinstance(proportions, pd.DataFrame):
            proportions = pd.DataFrame(
                proportions,
                index=deconv_ctx.spatial_prepared.obs_names,
                columns=deconv_ctx.cell_types,
            )
        else:
            proportions.index = deconv_ctx.spatial_prepared.obs_names

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            deconv_ctx.common_genes,
            method="FlashDeconv",
            device="CPU",
            sketch_dim=sketch_dim,
            lambda_spatial=lambda_spatial,
            n_hvg=n_hvg,
            n_markers_per_type=n_markers_per_type,
        )

        return proportions, stats

    except Exception as e:
        if isinstance(e, (DependencyError, ProcessingError)):
            raise
        raise ProcessingError(f"FlashDeconv deconvolution failed: {str(e)}") from e
