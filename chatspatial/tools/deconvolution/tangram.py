"""
Tangram deconvolution method.

Tangram maps single-cell RNA-seq data to spatial transcriptomics
data using the native tangram-sc library.
"""

import gc
from typing import Any

import pandas as pd

from ...utils.dependency_manager import require
from ...utils.device_utils import get_device
from ...utils.exceptions import ChatSpatialError, ProcessingError
from .base import PreparedDeconvolutionData, create_deconvolution_stats


def deconvolve(
    data: PreparedDeconvolutionData,
    n_epochs: int = 1000,
    mode: str = "cells",
    learning_rate: float = 0.1,
    density_prior: str = "rna_count_based",
    use_gpu: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Deconvolve spatial data using native Tangram library.

    Args:
        data: Prepared deconvolution data (immutable)
        n_epochs: Number of training epochs
        mode: Mapping mode - 'cells' or 'clusters'
        learning_rate: Optimizer learning rate
        density_prior: Spatial density prior - 'rna_count_based' or 'uniform'
        use_gpu: Whether to use GPU acceleration

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    tg = require("tangram", feature="Tangram deconvolution")

    try:
        # Data already copied in prepare_deconvolution
        spatial_data = data.spatial
        ref_data = data.reference

        # Tangram requires 'cell_type' column for cluster mode
        if "cell_type" not in ref_data.obs.columns:
            ref_data.obs["cell_type"] = ref_data.obs[data.cell_type_key]

        # Select training genes (tangram recommendation: 100-1000 genes).
        # Prefer highly variable genes from reference for biological relevance;
        # fall back to all common genes if HVG annotation is unavailable.
        n_max = min(500, len(data.common_genes))
        if "highly_variable" in ref_data.var.columns:
            hvg_set = set(ref_data.var_names[ref_data.var["highly_variable"]])
            hvg_common = [g for g in data.common_genes if g in hvg_set]
            training_genes = (
                hvg_common[:n_max] if hvg_common else data.common_genes[:n_max]
            )
        else:
            training_genes = data.common_genes[:n_max]

        # Preprocess with tangram (this sets up required annotations)
        tg.pp_adatas(ref_data, spatial_data, genes=training_genes)

        # Set device (supports CUDA, MPS, and CPU)
        device = get_device(prefer_gpu=use_gpu)

        # Map cells to space
        if mode == "clusters":
            ad_map = tg.map_cells_to_space(
                ref_data,
                spatial_data,
                mode="clusters",
                cluster_label=data.cell_type_key,
                density_prior=density_prior,
                num_epochs=n_epochs,
                learning_rate=learning_rate,
                device=device,
            )
        else:
            ad_map = tg.map_cells_to_space(
                ref_data,
                spatial_data,
                mode="cells",
                density_prior=density_prior,
                num_epochs=n_epochs,
                learning_rate=learning_rate,
                device=device,
            )

        reference_cell_types = list(ref_data.obs[data.cell_type_key].unique())
        mapping_matrix = ad_map.X

        if mode == "clusters":
            if mapping_matrix.shape != (
                len(reference_cell_types),
                spatial_data.n_obs,
            ):
                raise ProcessingError(
                    "Unexpected Tangram cluster mapping shape: "
                    f"{mapping_matrix.shape}, expected "
                    f"({len(reference_cell_types)}, {spatial_data.n_obs})"
                )
            if (
                not hasattr(ad_map, "obs")
                or data.cell_type_key not in ad_map.obs.columns
            ):
                raise ProcessingError(
                    "Tangram cluster mapping is missing its cell-type row labels."
                )
            cell_types = ad_map.obs[data.cell_type_key].map(str).tolist()
            if len(cell_types) != len(set(cell_types)):
                raise ProcessingError(
                    "Tangram cluster mapping returned duplicate cell-type labels."
                )
            proportions = pd.DataFrame(
                mapping_matrix.T, index=spatial_data.obs_names, columns=cell_types
            )
        else:
            if mapping_matrix.shape != (ref_data.n_obs, spatial_data.n_obs):
                raise ProcessingError(
                    "Unexpected Tangram cell mapping shape: "
                    f"{mapping_matrix.shape}, expected "
                    f"({ref_data.n_obs}, {spatial_data.n_obs})"
                )
            cell_type_series = ref_data.obs[data.cell_type_key]
            type_indicators = pd.get_dummies(cell_type_series)
            type_indicators = type_indicators.reindex(
                columns=reference_cell_types, fill_value=0
            )
            proportions_array = type_indicators.values.T @ mapping_matrix
            proportions = pd.DataFrame(
                proportions_array.T,
                index=spatial_data.obs_names,
                columns=[str(cell_type) for cell_type in reference_cell_types],
            )

        # Normalize to proportions
        row_sums = proportions.sum(axis=1)
        row_sums = row_sums.replace(0, 1)  # Avoid division by zero
        proportions = proportions.div(row_sums, axis=0)

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            data.common_genes,
            method="Tangram",
            device=device.upper(),
            n_epochs=n_epochs,
            mode=mode,
            density_prior=density_prior,
            n_training_genes=len(training_genes),
        )

        # Memory cleanup
        del ad_map, mapping_matrix
        del spatial_data, ref_data
        gc.collect()

        return proportions, stats

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"Tangram deconvolution failed: {e}") from e
