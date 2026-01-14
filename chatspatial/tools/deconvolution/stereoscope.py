"""
Stereoscope deconvolution method.

Stereoscope uses a two-stage training workflow:
1. Train RNAStereoscope model on reference data
2. Train SpatialStereoscope model on spatial data using RNA model
"""

import gc
from typing import Any

import pandas as pd

from ...utils.adata_utils import ensure_categorical
from ...utils.exceptions import ProcessingError
from .base import PreparedDeconvolutionData, create_deconvolution_stats


async def deconvolve(
    data: PreparedDeconvolutionData,
    n_epochs: int = 150000,
    learning_rate: float = 0.01,
    batch_size: int = 128,
    use_gpu: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Deconvolve spatial data using Stereoscope from scvi-tools.

    Args:
        data: Prepared deconvolution data (immutable)
        n_epochs: Total epochs (default: 150000, split 75K+75K)
        learning_rate: Learning rate (default: 0.01)
        batch_size: Minibatch size (default: 128)
        use_gpu: Use GPU acceleration

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    from scvi.external import RNAStereoscope, SpatialStereoscope

    try:
        # Create working copies (scvi-tools may modify in place)
        spatial_data = data.spatial.copy()
        ref_data = data.reference.copy()

        # Ensure categorical cell type
        ensure_categorical(ref_data, data.cell_type_key)

        cell_types = list(ref_data.obs[data.cell_type_key].cat.categories)

        # Calculate epoch split
        if n_epochs == 150000:
            rna_epochs, spatial_epochs = 75000, 75000
        else:
            rna_epochs = n_epochs // 2
            spatial_epochs = n_epochs - rna_epochs

        plan_kwargs = {"lr": learning_rate}
        accelerator = "gpu" if use_gpu else "cpu"

        # ===== Stage 1: Train RNAStereoscope =====
        RNAStereoscope.setup_anndata(ref_data, labels_key=data.cell_type_key)
        rna_model = RNAStereoscope(ref_data)

        train_kwargs = {
            "max_epochs": rna_epochs,
            "batch_size": batch_size,
            "plan_kwargs": plan_kwargs,
        }
        if use_gpu:
            train_kwargs["accelerator"] = accelerator
        rna_model.train(**train_kwargs)

        # ===== Stage 2: Train SpatialStereoscope =====
        SpatialStereoscope.setup_anndata(spatial_data)
        spatial_model = SpatialStereoscope.from_rna_model(spatial_data, rna_model)

        train_kwargs["max_epochs"] = spatial_epochs
        spatial_model.train(**train_kwargs)

        # Extract proportions
        proportions = pd.DataFrame(
            spatial_model.get_proportions(),
            index=spatial_data.obs_names,
            columns=cell_types,
        )

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            data.common_genes,
            method="Stereoscope",
            device="gpu" if use_gpu else "cpu",
            n_epochs=n_epochs,
            rna_epochs=rna_epochs,
            spatial_epochs=spatial_epochs,
            learning_rate=learning_rate,
        )

        # Memory cleanup
        del spatial_model, rna_model
        del spatial_data, ref_data
        gc.collect()

        return proportions, stats

    except Exception as e:
        if isinstance(e, ProcessingError):
            raise
        raise ProcessingError(f"Stereoscope deconvolution failed: {e}") from e
