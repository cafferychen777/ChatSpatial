"""
Cell2location deconvolution method.

Cell2location uses a two-stage training process:
1. Reference model (NB regression) learns cell type gene expression signatures
2. Cell2location model performs spatial mapping using these signatures
"""

import gc
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...utils.adata_utils import find_common_genes, validate_gene_overlap
from ...utils.dependency_manager import is_available, require
from ...utils.exceptions import DataError, ProcessingError
from ...utils.mcp_utils import suppress_output
from .base import (
    DeconvolutionContext,
    check_model_convergence,
    create_deconvolution_stats,
    prepare_anndata_for_counts,
)


async def _apply_gene_filtering(
    adata: "ad.AnnData",
    ctx: "ToolContext",
    cell_count_cutoff: int = 5,
    cell_percentage_cutoff2: float = 0.03,
    nonz_mean_cutoff: float = 1.12,
) -> "ad.AnnData":
    """Apply cell2location's official gene filtering.

    Reference: cell2location tutorial - "very permissive gene selection"

    Note: The original filter_genes function creates a matplotlib figure.
    We suppress this by using Agg backend and closing the figure immediately.
    """
    if not is_available("cell2location"):
        await ctx.warning(
            "cell2location.utils.filtering not available. "
            "Skipping gene filtering (may degrade results)."
        )
        return adata.copy()

    # Suppress matplotlib figure from filter_genes
    import matplotlib
    import matplotlib.pyplot as plt

    original_backend = matplotlib.get_backend()
    matplotlib.use("Agg")
    plt.ioff()

    try:
        from cell2location.utils.filtering import filter_genes

        selected = filter_genes(
            adata,
            cell_count_cutoff=cell_count_cutoff,
            cell_percentage_cutoff2=cell_percentage_cutoff2,
            nonz_mean_cutoff=nonz_mean_cutoff,
        )
        # Close any figures created by filter_genes
        plt.close("all")

        return adata[:, selected].copy()
    finally:
        # Restore original backend if needed
        if original_backend != "Agg":
            try:
                matplotlib.use(original_backend)
            except Exception:
                pass  # Backend might not be switchable


async def deconvolve(
    deconv_ctx: DeconvolutionContext,
    ref_model_epochs: int = 250,
    n_epochs: int = 30000,
    n_cells_per_spot: int = 30,
    detection_alpha: float = 20.0,
    use_gpu: bool = False,
    batch_key: Optional[str] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    apply_gene_filtering: bool = True,
    gene_filter_cell_count_cutoff: int = 5,
    gene_filter_cell_percentage_cutoff2: float = 0.03,
    gene_filter_nonz_mean_cutoff: float = 1.12,
    ref_model_lr: float = 0.002,
    cell2location_lr: float = 0.005,
    ref_model_train_size: float = 1.0,
    cell2location_train_size: float = 1.0,
    early_stopping: bool = False,
    early_stopping_patience: int = 45,
    early_stopping_threshold: float = 0.0,
    use_aggressive_training: bool = False,
    validation_size: float = 0.1,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using Cell2location.

    Args:
        deconv_ctx: Prepared DeconvolutionContext (note: we re-prepare with gene filtering)
        ref_model_epochs: Epochs for reference model (default: 250)
        n_epochs: Epochs for Cell2location model (default: 30000)
        n_cells_per_spot: Expected cells per location (default: 30)
        detection_alpha: RNA detection sensitivity (default: 20, NEW 2024)
        use_gpu: Use GPU acceleration
        batch_key: Column for batch correction
        categorical_covariate_keys: Technical covariates
        apply_gene_filtering: Apply official gene filtering (default: True)
        gene_filter_*: Gene filtering parameters
        ref_model_lr: Reference model learning rate (default: 0.002)
        cell2location_lr: Cell2location learning rate (default: 0.005)
        *_train_size: Training data fractions
        early_stopping*: Early stopping parameters
        use_aggressive_training: Use train_aggressive() method
        validation_size: Validation set size

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    require("cell2location")
    from cell2location.models import Cell2location, RegressionModel

    ctx = deconv_ctx.ctx
    cell_type_key = deconv_ctx.cell_type_key

    try:
        # Device selection
        device = "cpu"
        if use_gpu and is_available("torch"):
            import torch

            if torch.cuda.is_available():
                device = "cuda"

        # Cell2location needs special preparation with gene filtering
        # So we re-prepare from original data instead of using deconv_ctx.prepared
        ref = await prepare_anndata_for_counts(
            deconv_ctx.reference_adata, "Reference", ctx, require_int_dtype=False
        )
        sp = await prepare_anndata_for_counts(
            deconv_ctx.spatial_adata, "Spatial", ctx, require_int_dtype=False
        )

        # Apply gene filtering (cell2location-specific preprocessing)
        if apply_gene_filtering:
            ref = await _apply_gene_filtering(
                ref,
                ctx,
                cell_count_cutoff=gene_filter_cell_count_cutoff,
                cell_percentage_cutoff2=gene_filter_cell_percentage_cutoff2,
                nonz_mean_cutoff=gene_filter_nonz_mean_cutoff,
            )
            sp = await _apply_gene_filtering(
                sp,
                ctx,
                cell_count_cutoff=gene_filter_cell_count_cutoff,
                cell_percentage_cutoff2=gene_filter_cell_percentage_cutoff2,
                nonz_mean_cutoff=gene_filter_nonz_mean_cutoff,
            )

        # Find common genes after filtering
        common_genes = find_common_genes(ref.var_names, sp.var_names)
        validate_gene_overlap(
            common_genes,
            sp.n_vars,
            ref.n_vars,
            min_genes=100,
            source_name="spatial",
            target_name="reference",
        )

        # Subset to common genes
        ref = ref[:, common_genes]
        sp = sp[:, common_genes]

        # Ensure float32 for scvi-tools compatibility
        # Note: scipy sparse matrices' astype() preserves sparsity
        # Dense arrays are converted in-place without extra memory allocation
        if ref.X.dtype != np.float32:
            ref.X = ref.X.astype(np.float32)
        if sp.X.dtype != np.float32:
            sp.X = sp.X.astype(np.float32)

        # Handle NaN in cell types
        if ref.obs[cell_type_key].isna().any():
            await ctx.warning(f"Reference has NaN in {cell_type_key}. Excluding.")
            ref = ref[~ref.obs[cell_type_key].isna()].copy()

        # ===== Stage 1: Train Reference Model =====
        RegressionModel.setup_anndata(
            adata=ref,
            labels_key=cell_type_key,
            batch_key=batch_key,
            categorical_covariate_keys=categorical_covariate_keys,
        )

        ref_model = RegressionModel(ref)
        with suppress_output():
            train_kwargs = _build_train_kwargs(
                epochs=ref_model_epochs,
                lr=ref_model_lr,
                train_size=ref_model_train_size,
                device=device,
                early_stopping=early_stopping,
                early_stopping_patience=early_stopping_patience,
                validation_size=validation_size,
                use_aggressive=use_aggressive_training,
            )
            ref_model.train(**train_kwargs)

        # Check convergence
        converged, warning = check_model_convergence(ref_model, "ReferenceModel")
        if not converged and warning:
            await ctx.warning(warning)

        # Export reference signatures
        ref = ref_model.export_posterior(
            ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
        )
        ref_signatures = _extract_reference_signatures(ref)

        # ===== Stage 2: Train Cell2location Model =====
        Cell2location.setup_anndata(
            adata=sp,
            batch_key=batch_key,
            categorical_covariate_keys=categorical_covariate_keys,
        )

        cell2loc_model = Cell2location(
            sp,
            cell_state_df=ref_signatures,
            N_cells_per_location=n_cells_per_spot,
            detection_alpha=detection_alpha,
        )

        with suppress_output():
            train_kwargs = _build_train_kwargs(
                epochs=n_epochs,
                lr=cell2location_lr,
                train_size=cell2location_train_size,
                device=device,
                early_stopping=early_stopping,
                early_stopping_patience=early_stopping_patience,
                validation_size=validation_size,
                use_aggressive=use_aggressive_training,
            )
            cell2loc_model.train(**train_kwargs)

        # Check convergence
        converged, warning = check_model_convergence(cell2loc_model, "Cell2location")
        if not converged and warning:
            await ctx.warning(warning)

        # Export results
        sp = cell2loc_model.export_posterior(
            sp, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
        )

        # Extract cell abundance
        cell_abundance = _extract_cell_abundance(sp)

        # Create proportions DataFrame
        proportions = pd.DataFrame(
            cell_abundance,
            index=sp.obs_names,
            columns=ref_signatures.columns,
        )

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            common_genes,
            "Cell2location",
            device,
            n_epochs=n_epochs,
            n_cells_per_spot=n_cells_per_spot,
            detection_alpha=detection_alpha,
        )

        # Add model performance metrics
        if hasattr(cell2loc_model, "history") and cell2loc_model.history is not None:
            history = cell2loc_model.history
            if "elbo_train" in history and not history["elbo_train"].empty:
                stats["final_elbo"] = float(history["elbo_train"].iloc[-1])

        # Memory cleanup: release models and intermediate data
        # Cell2location is particularly memory-intensive with two-stage training
        del cell2loc_model, ref_model
        del ref, sp, ref_signatures
        gc.collect()

        return proportions, stats

    except Exception as e:
        if isinstance(e, (ProcessingError, DataError)):
            raise
        raise ProcessingError(f"Cell2location deconvolution failed: {str(e)}") from e


def _build_train_kwargs(
    epochs: int,
    lr: float,
    train_size: float,
    device: str,
    early_stopping: bool,
    early_stopping_patience: int,
    validation_size: float,
    use_aggressive: bool,
) -> Dict[str, Any]:
    """Build training kwargs for scvi-tools models."""
    if use_aggressive:
        kwargs = {"max_epochs": epochs, "lr": lr}
        if device == "cuda":
            kwargs["accelerator"] = "gpu"
        if early_stopping:
            kwargs["early_stopping"] = True
            kwargs["early_stopping_patience"] = early_stopping_patience
            kwargs["check_val_every_n_epoch"] = 1
            kwargs["train_size"] = 1.0 - validation_size
        else:
            kwargs["train_size"] = train_size
    else:
        kwargs = {
            "max_epochs": epochs,
            "batch_size": 2500,
            "lr": lr,
            "train_size": train_size,
        }
        if device == "cuda":
            kwargs["accelerator"] = "gpu"
    return kwargs


def _extract_reference_signatures(ref: "ad.AnnData") -> pd.DataFrame:
    """Extract reference signatures from trained RegressionModel."""
    factor_names = ref.uns["mod"]["factor_names"]
    cols = [f"means_per_cluster_mu_fg_{i}" for i in factor_names]

    if "means_per_cluster_mu_fg" in ref.varm:
        signatures = ref.varm["means_per_cluster_mu_fg"][cols].copy()
    else:
        signatures = ref.var[cols].copy()

    signatures.columns = factor_names
    return signatures


def _extract_cell_abundance(sp: "ad.AnnData"):
    """Extract cell abundance from Cell2location results.

    Cell2location stores results as DataFrames with prefixed column names like
    'q05cell_abundance_w_sf_CellType'. We need to extract the values and
    return them as a numpy array for consistent downstream processing.
    """
    possible_keys = [
        "q05_cell_abundance_w_sf",
        "means_cell_abundance_w_sf",
        "q50_cell_abundance_w_sf",
    ]

    for key in possible_keys:
        if key in sp.obsm:
            result = sp.obsm[key]
            # Cell2location returns DataFrame with prefixed column names
            # Extract values as numpy array for consistent processing
            if hasattr(result, "values"):
                return result.values
            return result

    raise ProcessingError(
        f"Cell2location did not produce expected output. "
        f"Available keys: {list(sp.obsm.keys())}"
    )
