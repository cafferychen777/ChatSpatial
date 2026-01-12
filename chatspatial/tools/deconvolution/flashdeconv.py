"""
FlashDeconv deconvolution method.

FlashDeconv is an ultra-fast spatial transcriptomics deconvolution method
that uses random sketching for O(N) time complexity.
"""

from typing import Any, Dict, Tuple

import pandas as pd

from ...utils.dependency_manager import is_available
from ...utils.exceptions import DependencyError, ProcessingError
from .base import PreparedDeconvolutionData, create_deconvolution_stats


async def deconvolve(
    data: PreparedDeconvolutionData,
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
        data: Prepared deconvolution data (immutable)
        sketch_dim: Dimension for random sketching (default: 512)
        lambda_spatial: Spatial regularization strength (default: 5000.0)
        n_hvg: Number of highly variable genes to use (default: 2000)
        n_markers_per_type: Number of marker genes per cell type (default: 50)

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    if not is_available("flashdeconv"):
        raise DependencyError(
            "FlashDeconv is not available. Install with: pip install flashdeconv"
        )

    try:
        import flashdeconv as fd

        # Create copies (FlashDeconv modifies in place)
        adata_st = data.spatial.copy()
        reference = data.reference.copy()

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

        # Ensure DataFrame format
        if not isinstance(proportions, pd.DataFrame):
            proportions = pd.DataFrame(
                proportions,
                index=data.spatial.obs_names,
                columns=data.cell_types,
            )
        else:
            proportions.index = data.spatial.obs_names

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            data.common_genes,
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
