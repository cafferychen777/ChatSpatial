"""
DestVI deconvolution method.

DestVI performs multi-resolution deconvolution by first training a CondSCVI
model on reference data, then using it to initialize a DestVI model.
"""

import gc
from typing import Any

import numpy as np
import pandas as pd
import torch
from scipy import sparse

from ...utils.dependency_manager import is_available
from ...utils.device_utils import get_accelerator
from ...utils.exceptions import DataError, DependencyError, ProcessingError
from .base import PreparedDeconvolutionData, create_deconvolution_stats


def _ensure_float32_x(adata: Any) -> None:
    if adata.X.dtype != np.float32:
        adata.X = adata.X.astype(np.float32, copy=False)
    if sparse.issparse(adata.X):
        adata.X.data = adata.X.data.astype(np.float32, copy=False)


def _patch_condscvi_vamp_prior_tensors(condscvi_model: Any) -> None:
    original_get_vamp_prior = getattr(condscvi_model, "get_vamp_prior", None)
    if original_get_vamp_prior is None:
        return

    def _get_vamp_prior_tensors(*args: Any, **kwargs: Any) -> dict[str, torch.Tensor]:
        return {
            key: torch.as_tensor(value, dtype=torch.float32)
            for key, value in original_get_vamp_prior(*args, **kwargs).items()
        }

    condscvi_model.get_vamp_prior = _get_vamp_prior_tensors


def deconvolve(
    data: PreparedDeconvolutionData,
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
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Deconvolve spatial data using DestVI from scvi-tools.

    Args:
        data: Prepared deconvolution data (immutable)
        n_epochs: Base epoch count. CondSCVI trains for max(400, n_epochs//5),
            DestVI trains for max(200, n_epochs//10). Default 2500 yields
            500 CondSCVI + 250 DestVI epochs.
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
        # Data already copied in prepare_deconvolution
        spatial_data = data.spatial
        ref_data = data.reference
        _ensure_float32_x(spatial_data)
        _ensure_float32_x(ref_data)

        # Validate cell types
        if data.n_cell_types < 2:
            raise DataError(
                f"Reference needs at least 2 cell types, found {data.n_cell_types}"
            )

        # Calculate epoch distribution
        condscvi_epochs = max(400, n_epochs // 5)
        destvi_epochs = max(200, n_epochs // 10)

        # Device setting
        accelerator = get_accelerator(prefer_gpu=use_gpu)
        condscvi_plan_kwargs = {"lr": learning_rate}
        destvi_plan_kwargs = {"lr": learning_rate, "ct_sparsity_weight": l1_reg}

        # ===== Stage 1: Train CondSCVI on reference =====
        scvi.model.CondSCVI.setup_anndata(
            ref_data,
            labels_key=data.cell_type_key,
            batch_key=None,
        )

        condscvi_model = scvi.model.CondSCVI(
            ref_data,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            prior="normal",
        )

        condscvi_model.train(
            max_epochs=condscvi_epochs,
            accelerator=accelerator,
            train_size=train_size,
            plan_kwargs=condscvi_plan_kwargs,
        )

        # ===== Stage 2: Train DestVI on spatial =====
        scvi.model.DestVI.setup_anndata(spatial_data)
        _patch_condscvi_vamp_prior_tensors(condscvi_model)

        destvi_model = scvi.model.DestVI.from_rna_model(
            spatial_data,
            condscvi_model,
            vamp_prior_p=vamp_prior_p,
        )

        destvi_model.train(
            max_epochs=destvi_epochs,
            accelerator=accelerator,
            train_size=train_size,
            plan_kwargs=destvi_plan_kwargs,
        )

        # Get proportions
        proportions = destvi_model.get_proportions()
        proportions.index = spatial_data.obs_names

        if proportions.empty or len(proportions) != spatial_data.n_obs:
            raise ProcessingError("Failed to extract valid proportions from DestVI")

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            data.common_genes,
            method="DestVI",
            device=accelerator,
            n_epochs=n_epochs,
            condscvi_epochs=condscvi_epochs,
            destvi_epochs=destvi_epochs,
            n_hidden=n_hidden,
            n_latent=n_latent,
        )

        # Memory cleanup
        del destvi_model, condscvi_model
        del spatial_data, ref_data
        gc.collect()

        return proportions, stats

    except Exception as e:
        if isinstance(e, (DependencyError, DataError, ProcessingError)):
            raise
        raise ProcessingError(f"DestVI deconvolution failed: {e}") from e
