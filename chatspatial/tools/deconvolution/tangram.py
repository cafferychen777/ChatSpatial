"""
Tangram deconvolution method.

Tangram maps single-cell RNA-seq data to spatial transcriptomics
data using the scvi.external.Tangram wrapper from scvi-tools.
"""

import gc
from typing import TYPE_CHECKING, Any, Dict, Tuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    pass

from ...utils.dependency_manager import is_available, require
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
    """Deconvolve spatial data using Tangram via scvi-tools.

    Args:
        deconv_ctx: Prepared DeconvolutionContext
        n_epochs: Number of training epochs
        mode: Mapping mode - 'cells', 'clusters', or 'constrained'
        learning_rate: Optimizer learning rate
        density_prior: Spatial density prior - 'rna_count_based' or 'uniform'
        use_gpu: Whether to use GPU acceleration

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    if not is_available("scvi-tools"):
        raise DependencyError(
            "scvi-tools is required for Tangram. Install with: pip install scvi-tools"
        )

    require("mudata")
    import mudata as md
    from scvi.external import Tangram

    cell_type_key = deconv_ctx.cell_type_key

    try:
        # Get subset data
        spatial_data, ref_data = deconv_ctx.get_subset_data()
        common_genes = deconv_ctx.common_genes

        # Create density prior
        # Memory optimization: ravel() returns view when possible, flatten() always copies
        if density_prior == "rna_count_based":
            density_values = np.asarray(spatial_data.X.sum(axis=1)).ravel()
        else:
            density_values = np.ones(spatial_data.n_obs)

        density_values = density_values / density_values.sum()
        spatial_data.obs["density_prior"] = density_values

        # Create MuData object
        mdata = md.MuData({"sc_train": ref_data, "sp_train": spatial_data})

        # Setup MuData for Tangram
        Tangram.setup_mudata(
            mdata,
            density_prior_key="density_prior",
            modalities={
                "density_prior_key": "sp_train",
                "sc_layer": "sc_train",
                "sp_layer": "sp_train",
            },
        )

        # Create model
        if mode == "constrained":
            target_count = max(1, int(spatial_data.n_obs * 0.1))
            tangram_model = Tangram(mdata, constrained=True, target_count=target_count)
        else:
            tangram_model = Tangram(mdata, constrained=False)

        # Train
        train_kwargs = {"max_epochs": n_epochs, "lr": learning_rate}
        if use_gpu:
            train_kwargs["accelerator"] = "gpu"

        tangram_model.train(**train_kwargs)

        # Get mapping matrix and calculate cell type proportions
        mapping_matrix = tangram_model.get_mapper_matrix()
        cell_types = ref_data.obs[cell_type_key].unique()

        # Memory-efficient vectorized aggregation using pandas groupby logic
        # Instead of loop with repeated array creation, use matrix multiplication
        # Create cell type indicator matrix: (n_cells,) categorical -> one-hot (n_cells x n_types)
        cell_type_series = ref_data.obs[cell_type_key]
        type_indicators = pd.get_dummies(cell_type_series)

        # Reorder columns to match cell_types order
        type_indicators = type_indicators.reindex(columns=cell_types, fill_value=0)

        # Matrix multiply: (n_types x n_cells) @ (n_cells x n_spots) = (n_types x n_spots)
        # This computes sum of mapping weights for each cell type in one operation
        proportions_array = type_indicators.values.T @ mapping_matrix

        # Create DataFrame (proportions_array is n_types x n_spots, need to transpose)
        proportions = pd.DataFrame(
            proportions_array.T, index=spatial_data.obs_names, columns=cell_types
        )

        # Normalize to proportions (Tangram returns equivalent cell counts)
        row_sums = proportions.sum(axis=1)
        row_sums = row_sums.replace(0, 1)  # Avoid division by zero
        proportions = proportions.div(row_sums, axis=0)

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            common_genes,
            "Tangram",
            device="GPU" if use_gpu else "CPU",
            n_epochs=n_epochs,
            mode=mode,
            density_prior=density_prior,
        )

        # Memory cleanup: release model and intermediate data
        del tangram_model, mdata, mapping_matrix, type_indicators
        del spatial_data, ref_data
        gc.collect()

        return proportions, stats

    except Exception as e:
        if isinstance(e, (DependencyError, ProcessingError)):
            raise
        raise ProcessingError(f"Tangram deconvolution failed: {str(e)}") from e
