"""
Stereoscope deconvolution method.

Stereoscope uses a two-stage training workflow:
1. Train RNAStereoscope model on reference data
2. Train SpatialStereoscope model on spatial data using RNA model
"""

import gc
from typing import TYPE_CHECKING, Any, Dict, Tuple

import pandas as pd

if TYPE_CHECKING:
    pass

from ...utils.adata_utils import ensure_categorical
from ...utils.exceptions import ProcessingError
from .base import DeconvolutionContext, create_deconvolution_stats


async def deconvolve(
    deconv_ctx: DeconvolutionContext,
    n_epochs: int = 150000,
    learning_rate: float = 0.01,
    batch_size: int = 128,
    use_gpu: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using Stereoscope from scvi-tools.

    Args:
        deconv_ctx: Prepared DeconvolutionContext
        n_epochs: Total epochs (default: 150000, split 75K+75K)
        learning_rate: Learning rate (default: 0.01)
        batch_size: Minibatch size (default: 128)
        use_gpu: Use GPU acceleration

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    from scvi.external import RNAStereoscope, SpatialStereoscope

    try:
        spatial_data, ref_data = deconv_ctx.get_subset_data()
        cell_type_key = deconv_ctx.cell_type_key
        common_genes = deconv_ctx.common_genes

        # Ensure categorical cell type
        ensure_categorical(ref_data, cell_type_key)

        cell_types = list(ref_data.obs[cell_type_key].cat.categories)

        # Calculate epoch split
        if n_epochs == 150000:
            rna_epochs, spatial_epochs = 75000, 75000
        else:
            rna_epochs = n_epochs // 2
            spatial_epochs = n_epochs - rna_epochs

        plan_kwargs = {"lr": learning_rate}
        accelerator = "gpu" if use_gpu else "cpu"

        # ===== Stage 1: Train RNAStereoscope =====
        RNAStereoscope.setup_anndata(ref_data, labels_key=cell_type_key)
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
            common_genes,
            "Stereoscope",
            device="gpu" if use_gpu else "cpu",
            n_epochs=n_epochs,
            rna_epochs=rna_epochs,
            spatial_epochs=spatial_epochs,
            learning_rate=learning_rate,
        )

        # Memory cleanup: release models and intermediate data
        del spatial_model, rna_model
        del spatial_data, ref_data
        gc.collect()

        return proportions, stats

    except Exception as e:
        if isinstance(e, ProcessingError):
            raise
        raise ProcessingError(f"Stereoscope deconvolution failed: {str(e)}") from e
