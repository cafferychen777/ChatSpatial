"""
Tangram deconvolution method.

Tangram maps single-cell RNA-seq data to spatial transcriptomics
data using the native tangram-sc library.
"""

import gc
from typing import TYPE_CHECKING, Any, Dict, Tuple

import pandas as pd

if TYPE_CHECKING:
    pass

from ...utils.exceptions import DependencyError, ProcessingError
from .base import DeconvolutionContext, create_deconvolution_stats


async def deconvolve(
    deconv_ctx: DeconvolutionContext,
    n_epochs: int = 1000,
    mode: str = "cells",
    learning_rate: float = 0.1,
    density_prior: str = "rna_count_based",
    use_gpu: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using native Tangram library.

    Args:
        deconv_ctx: Prepared DeconvolutionContext
        n_epochs: Number of training epochs
        mode: Mapping mode - 'cells' or 'clusters'
        learning_rate: Optimizer learning rate
        density_prior: Spatial density prior - 'rna_count_based' or 'uniform'
        use_gpu: Whether to use GPU acceleration

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    # Check for tangram package (installed as tangram-sc, imported as tangram)
    try:
        import tangram as tg
    except ImportError as e:
        raise DependencyError(
            "tangram-sc is required for Tangram. Install with: pip install tangram-sc"
        ) from e

    cell_type_key = deconv_ctx.cell_type_key

    try:
        # Get subset data
        spatial_data, ref_data = deconv_ctx.get_subset_data()
        common_genes = deconv_ctx.common_genes

        # Tangram requires 'cell_type' column for cluster mode
        if "cell_type" not in ref_data.obs.columns:
            ref_data.obs["cell_type"] = ref_data.obs[cell_type_key]

        # Select marker genes for training (tangram recommendation: 100-1000 genes)
        # Use up to 500 genes from common genes for efficiency
        n_training_genes = min(500, len(common_genes))
        training_genes = common_genes[:n_training_genes]

        # Preprocess with tangram (this sets up required annotations)
        tg.pp_adatas(ref_data, spatial_data, genes=training_genes)

        # Set device
        device = "cuda:0" if use_gpu else "cpu"

        # Map cells to space
        if mode == "clusters":
            # Cluster mode: aggregate by cell type before mapping
            ad_map = tg.map_cells_to_space(
                ref_data,
                spatial_data,
                mode="clusters",
                cluster_label=cell_type_key,
                density_prior=density_prior,
                num_epochs=n_epochs,
                learning_rate=learning_rate,
                device=device,
            )
        else:
            # Default cells mode
            ad_map = tg.map_cells_to_space(
                ref_data,
                spatial_data,
                mode="cells",
                density_prior=density_prior,
                num_epochs=n_epochs,
                learning_rate=learning_rate,
                device=device,
            )

        # Get mapping matrix (cells x spots)
        mapping_matrix = ad_map.X

        # Calculate cell type proportions from mapping
        cell_types = ref_data.obs[cell_type_key].unique()

        # Create cell type indicator matrix
        cell_type_series = ref_data.obs[cell_type_key]
        type_indicators = pd.get_dummies(cell_type_series)
        type_indicators = type_indicators.reindex(columns=cell_types, fill_value=0)

        # Matrix multiply: (n_types x n_cells) @ (n_cells x n_spots) = (n_types x n_spots)
        proportions_array = type_indicators.values.T @ mapping_matrix

        # Create DataFrame
        proportions = pd.DataFrame(
            proportions_array.T, index=spatial_data.obs_names, columns=cell_types
        )

        # Normalize to proportions
        row_sums = proportions.sum(axis=1)
        row_sums = row_sums.replace(0, 1)  # Avoid division by zero
        proportions = proportions.div(row_sums, axis=0)

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            common_genes,
            "Tangram",
            device=device.upper(),
            n_epochs=n_epochs,
            mode=mode,
            density_prior=density_prior,
            n_training_genes=n_training_genes,
        )

        # Memory cleanup
        del ad_map, mapping_matrix, type_indicators
        del spatial_data, ref_data
        gc.collect()

        return proportions, stats

    except Exception as e:
        if isinstance(e, (DependencyError, ProcessingError)):
            raise
        raise ProcessingError(f"Tangram deconvolution failed: {str(e)}") from e
