"""
DestVI deconvolution method.

DestVI performs multi-resolution deconvolution by first training a CondSCVI
model on reference data, then using it to initialize a DestVI model.
"""

import gc
from typing import TYPE_CHECKING, Any, Dict, Tuple

import pandas as pd

if TYPE_CHECKING:
    pass

from ...utils.dependency_manager import is_available
from ...utils.exceptions import DataError, DependencyError, ProcessingError
from .base import DeconvolutionContext, create_deconvolution_stats


async def deconvolve(
    deconv_ctx: DeconvolutionContext,
    n_epochs: int = 10000,
    n_hidden: int = 128,
    n_latent: int = 10,
    n_layers: int = 1,
    dropout_rate: float = 0.1,
    learning_rate: float = 1e-3,
    train_size: float = 0.9,
    vamp_prior_p: int = 15,
    l1_reg: float = 10.0,
    use_gpu: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using DestVI from scvi-tools.

    Args:
        deconv_ctx: Prepared DeconvolutionContext
        n_epochs: Total epochs (split between CondSCVI and DestVI)
        n_hidden: Hidden units in neural networks
        n_latent: Latent space dimensionality
        n_layers: Number of layers
        dropout_rate: Dropout rate
        learning_rate: Learning rate
        train_size: Fraction for training (default: 0.9)
        vamp_prior_p: VampPrior components (default: 15)
        l1_reg: L1 regularization (default: 10.0)
        use_gpu: Use GPU acceleration

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    if not is_available("scvi-tools"):
        raise DependencyError(
            "scvi-tools is required for DestVI. Install with: pip install scvi-tools"
        )

    import scvi

    try:
        # Get subset data from context
        spatial_data, ref_data = deconv_ctx.get_subset_data()
        cell_type_key = deconv_ctx.cell_type_key
        common_genes = deconv_ctx.common_genes

        # Validate cell types
        if len(deconv_ctx.cell_types) < 2:
            raise DataError(
                f"Reference needs at least 2 cell types, found {len(deconv_ctx.cell_types)}"
            )

        # Calculate epoch distribution
        condscvi_epochs = max(400, n_epochs // 5)
        destvi_epochs = max(200, n_epochs // 10)

        # Device setting
        accelerator = "gpu" if use_gpu else "cpu"
        plan_kwargs = {"lr": learning_rate}

        # ===== Stage 1: Train CondSCVI on reference =====
        scvi.model.CondSCVI.setup_anndata(
            ref_data,
            labels_key=cell_type_key,
            batch_key=None,
        )

        condscvi_model = scvi.model.CondSCVI(
            ref_data,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
        )

        condscvi_model.train(
            max_epochs=condscvi_epochs,
            accelerator=accelerator,
            train_size=train_size,
            plan_kwargs=plan_kwargs,
        )

        # ===== Stage 2: Train DestVI on spatial =====
        scvi.model.DestVI.setup_anndata(spatial_data)

        destvi_model = scvi.model.DestVI.from_rna_model(
            spatial_data,
            condscvi_model,
            vamp_prior_p=vamp_prior_p,
            l1_reg=l1_reg,
        )

        destvi_model.train(
            max_epochs=destvi_epochs,
            accelerator=accelerator,
            train_size=train_size,
            plan_kwargs=plan_kwargs,
        )

        # Get proportions
        proportions = destvi_model.get_proportions()
        proportions.index = spatial_data.obs_names

        if proportions.empty or len(proportions) != spatial_data.n_obs:
            raise ProcessingError("Failed to extract valid proportions from DestVI")

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            common_genes,
            "DestVI",
            device="gpu" if use_gpu else "cpu",
            n_epochs=n_epochs,
            condscvi_epochs=condscvi_epochs,
            destvi_epochs=destvi_epochs,
            n_hidden=n_hidden,
            n_latent=n_latent,
        )

        # Memory cleanup: release models and intermediate data
        del destvi_model, condscvi_model
        del spatial_data, ref_data
        gc.collect()

        return proportions, stats

    except Exception as e:
        if isinstance(e, (DependencyError, DataError, ProcessingError)):
            raise
        raise ProcessingError(f"DestVI deconvolution failed: {str(e)}") from e
