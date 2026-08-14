"""
Integration tools for spatial transcriptomics data.
"""

import logging
from collections.abc import Sequence
from typing import TYPE_CHECKING, Optional, TypedDict, Union

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from ..models.analysis import IntegrationResult
from ..models.data import IntegrationParameters
from ..utils.adata_utils import check_is_integer_counts
from ..utils.async_utils import run_sync
from ..utils.compute import ensure_highly_variable_genes
from ..utils.dependency_manager import require
from ..utils.device_utils import get_accelerator
from ..utils.exceptions import (
    ChatSpatialError,
    DataError,
    ParameterError,
    ProcessingError,
)

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

from ..utils.adata_utils import (
    get_spatial_key,
    require_spatial_coords,
    store_analysis_metadata,
    validate_adata_basics,
)
from ..utils.results_export import export_analysis_result

logger = logging.getLogger(__name__)

_ORIGINAL_OBS_NAME_KEY = "_chatspatial_original_obs_name"


class _ScviSettings(TypedDict):
    """Typed scVI arguments shared by execution and provenance metadata."""

    n_hidden: int
    n_latent: int
    n_layers: int
    dropout_rate: float
    gene_likelihood: str
    n_epochs: int | None
    use_gpu: bool


def _available_obs_key(columns: pd.Index, preferred: str) -> str:
    """Return a deterministic observation-column name without overwriting data."""
    if preferred not in columns:
        return preferred

    suffix = 1
    while f"{preferred}_{suffix}" in columns:
        suffix += 1
    return f"{preferred}_{suffix}"


def _normalize_sample_keys(
    adatas: Sequence[ad.AnnData], sample_keys: Sequence[str] | None
) -> list[str]:
    """Validate or create stable per-input keys for observation identities."""
    keys = (
        [f"sample_{index}" for index in range(len(adatas))]
        if sample_keys is None
        else [str(key) for key in sample_keys]
    )
    if len(keys) != len(adatas):
        raise ParameterError(
            "sample_keys must contain exactly one key for each input dataset."
        )
    if any(not key for key in keys):
        raise ParameterError("sample_keys cannot contain empty values.")
    if len(set(keys)) != len(keys):
        raise ParameterError("sample_keys must be unique.")
    return keys


def _concatenate_samples(
    adatas: Sequence[ad.AnnData],
    *,
    batch_key: str,
    sample_keys: Sequence[str] | None,
) -> ad.AnnData:
    """Concatenate inputs without mutation and create globally unique cell IDs."""
    keys = _normalize_sample_keys(adatas, sample_keys)
    for key, adata in zip(keys, adatas, strict=True):
        if not adata.obs_names.is_unique:
            raise DataError(
                f"Input dataset '{key}' contains duplicate observation names; "
                "integration requires unique cell IDs within each sample."
            )

    original_obs_names = np.concatenate(
        [adata.obs_names.astype(str).to_numpy() for adata in adatas]
    )
    source_samples = np.concatenate(
        [np.repeat(key, adata.n_obs) for key, adata in zip(keys, adatas, strict=True)]
    )

    combined = ad.concat(
        adatas,
        keys=keys,
        join="outer",
        merge="first",
        uns_merge="first",
        index_unique=":",
    )
    if not combined.obs_names.is_unique:
        raise DataError(
            "Integration could not construct globally unique observation names."
        )

    original_name_key = _available_obs_key(combined.obs.columns, _ORIGINAL_OBS_NAME_KEY)
    combined.obs[original_name_key] = original_obs_names

    if batch_key not in combined.obs:
        combined.obs[batch_key] = source_samples
    elif combined.obs[batch_key].isna().any():
        # Concatenating a mixture of labeled and unlabeled inputs creates missing
        # values. Fill only those cells, preserving labels supplied by callers.
        batch_values = combined.obs[batch_key].astype(object)
        missing_batch = batch_values.isna().to_numpy()
        batch_values.iloc[missing_batch] = source_samples[missing_batch]
        combined.obs[batch_key] = batch_values

    return combined


def _orient_harmony_embedding(embedding: object, n_obs: int) -> np.ndarray:
    """Return a finite Harmony embedding with observations on the first axis."""
    try:
        normalized = np.asarray(embedding, dtype=np.float32)
    except (TypeError, ValueError) as exc:
        raise ProcessingError("Harmony returned a non-numeric embedding.") from exc

    if normalized.ndim != 2:
        raise ProcessingError("Harmony embedding must be two-dimensional.")

    if normalized.shape[0] == n_obs:
        oriented = normalized
    elif normalized.shape[1] == n_obs:
        oriented = normalized.T
    else:
        raise ProcessingError(
            "Harmony embedding does not match the dataset: "
            f"expected {n_obs} observations, received shape {normalized.shape}."
        )

    if oriented.shape[1] == 0:
        raise ProcessingError("Harmony returned zero embedding components.")
    if not np.isfinite(oriented).all():
        raise ProcessingError("Harmony returned non-finite embedding values.")
    return oriented


def _reassemble_scanorama_embeddings(
    batch_indices: list[np.ndarray],
    embeddings: list[np.ndarray],
    n_obs: int,
) -> np.ndarray:
    """Restore batch-wise Scanorama embeddings to the original observation order."""
    if not embeddings:
        raise ProcessingError("Scanorama returned no batch embeddings.")
    if len(batch_indices) != len(embeddings):
        raise ProcessingError(
            "Scanorama returned a different number of embeddings than input batches: "
            f"expected {len(batch_indices)}, received {len(embeddings)}."
        )

    n_components: int | None = None
    reassembled: np.ndarray | None = None
    assigned = np.zeros(n_obs, dtype=bool)

    for batch_number, (indices, embedding) in enumerate(
        zip(batch_indices, embeddings, strict=True)
    ):
        normalized_indices = np.asarray(indices, dtype=int)
        try:
            normalized_embedding = np.asarray(embedding, dtype=np.float32)
        except (TypeError, ValueError) as exc:
            raise ProcessingError(
                f"Scanorama embedding {batch_number} must be numeric."
            ) from exc
        if normalized_indices.ndim != 1:
            raise ProcessingError(
                f"Scanorama batch indices {batch_number} must be one-dimensional."
            )
        if normalized_embedding.ndim != 2:
            raise ProcessingError(
                f"Scanorama embedding {batch_number} must be two-dimensional."
            )
        if normalized_embedding.shape[0] != len(normalized_indices):
            raise ProcessingError(
                f"Scanorama embedding {batch_number} has "
                f"{normalized_embedding.shape[0]} rows for "
                f"{len(normalized_indices)} observations."
            )

        if n_components is None:
            n_components = normalized_embedding.shape[1]
            if n_components == 0:
                raise ProcessingError("Scanorama returned zero embedding components.")
            reassembled = np.empty((n_obs, n_components), dtype=np.float32)
        elif normalized_embedding.shape[1] != n_components:
            raise ProcessingError(
                "Scanorama returned inconsistent embedding dimensions across batches."
            )
        if not np.isfinite(normalized_embedding).all():
            raise ProcessingError(
                f"Scanorama embedding {batch_number} contains non-finite values."
            )

        if (
            np.any(normalized_indices < 0)
            or np.any(normalized_indices >= n_obs)
            or len(np.unique(normalized_indices)) != len(normalized_indices)
            or assigned[normalized_indices].any()
        ):
            raise ProcessingError(
                "Scanorama batch indices are duplicated or outside the dataset bounds."
            )

        assert reassembled is not None
        reassembled[normalized_indices] = normalized_embedding
        assigned[normalized_indices] = True

    if not assigned.all():
        raise ProcessingError(
            f"Scanorama output is missing embeddings for {int((~assigned).sum())} "
            "observations."
        )

    assert reassembled is not None
    return reassembled


def _clean_concatenated_metadata(combined: ad.AnnData) -> None:
    """Normalize outer-join metadata and remove incomplete graph artifacts."""
    for column in combined.var.columns:
        values = combined.var[column]
        is_label_column = pd.api.types.is_object_dtype(
            values.dtype
        ) or pd.api.types.is_string_dtype(values.dtype)
        if not is_label_column or not values.isna().any():
            continue

        unique_values = values.dropna().unique()
        if set(unique_values).issubset({True, False, "True", "False"}):
            mapped = values.map(
                {True: True, False: False, "True": True, "False": False}
            )
            combined.var[column] = mapped.astype("boolean").fillna(False).astype(bool)
        else:
            combined.var[column] = values.fillna("").astype(str)

    combined.obsm.pop("X_diffmap", None)
    combined.uns.pop("diffmap_evals", None)


def _prepare_combined_input(
    adatas: Union[list[ad.AnnData], ad.AnnData],
    *,
    batch_key: str,
    sample_keys: Sequence[str] | None,
) -> ad.AnnData:
    """Build and validate the integration workspace."""
    if isinstance(adatas, list):
        if len(adatas) < 2:
            raise ParameterError(
                f"Integration requires at least 2 datasets, got {len(adatas)}. "
                "Use preprocess_data for single dataset processing."
            )
        combined = _concatenate_samples(
            adatas,
            batch_key=batch_key,
            sample_keys=sample_keys,
        )
        _clean_concatenated_metadata(combined)
    else:
        combined = adatas
        if batch_key not in combined.obs:
            raise ParameterError(
                f"Merged dataset is missing batch information key '{batch_key}'"
            )

    is_integer, _, _ = check_is_integer_counts(combined.X)
    if is_integer:
        raise DataError("Data appears to be raw counts. Run preprocessing first.")
    validate_adata_basics(combined, min_obs=10, min_vars=10, check_empty_ratio=True)
    return combined


def _ensure_integration_hvgs(combined: ad.AnnData, batch_key: str) -> None:
    """Ensure a useful highly-variable-gene mask exists after concatenation."""
    n_hvg = int(combined.var.get("highly_variable", pd.Series(dtype=bool)).sum())
    if "highly_variable" not in combined.var:
        reason = "No highly variable genes marked after merge"
    elif n_hvg == 0:
        reason = "No genes marked as highly variable after merge"
    elif n_hvg < 50:
        reason = f"Very few HVGs ({n_hvg})"
    else:
        return

    logger.warning("%s, recalculating with batch correction", reason)
    # Merging samples with different gene sets leaves genes that are all zero,
    # which makes scanpy's mean-binned dispersion undefined.
    if ensure_highly_variable_genes(
        combined,
        n_top_genes=2000,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        batch_key=batch_key,
    ):
        logger.warning(
            "Scanpy HVG binning was unusable after the merge; "
            "selected HVGs by finite variance instead."
        )


def _run_scvi_integration(
    combined: ad.AnnData,
    *,
    batch_key: str,
    params: IntegrationParameters | None,
) -> ad.AnnData:
    """Run scVI integration and record its method-specific provenance."""
    settings: _ScviSettings = {
        "n_hidden": params.scvi_n_hidden if params else 128,
        "n_latent": params.scvi_n_latent if params else 10,
        "n_layers": params.scvi_n_layers if params else 1,
        "dropout_rate": params.scvi_dropout_rate if params else 0.1,
        "gene_likelihood": params.scvi_gene_likelihood if params else "zinb",
        "n_epochs": params.n_epochs if params else None,
        "use_gpu": params.use_gpu if params else False,
    }
    try:
        combined = integrate_with_scvi(combined, batch_key=batch_key, **settings)
    except ChatSpatialError:
        raise
    except Exception as exc:
        raise ProcessingError(
            f"scVI integration failed: {exc}. "
            "Ensure data is preprocessed and has ≥2 batches."
        ) from exc

    sc.tl.umap(combined)
    _store_integration_metadata(
        combined,
        method="scvi",
        batch_key=batch_key,
        parameters={"batch_key": batch_key, **settings},
        results_keys={"obsm": ["X_scvi"]},
    )
    return combined


def _nonzero_variance_mask(matrix) -> np.ndarray:
    """Return a dense boolean mask without densifying sparse expression data."""
    if sparse.issparse(matrix):
        means = np.asarray(matrix.mean(axis=0)).ravel()
        mean_squares = np.asarray(matrix.power(2).mean(axis=0)).ravel()
        variances = mean_squares - means**2
    else:
        variances = np.var(matrix, axis=0)
    return np.asarray(variances > 0).ravel()


def _scale_integration_data(combined: ad.AnnData) -> ad.AnnData:
    """Scale on isolated candidates so a failed strategy leaves no partial state."""
    errors: list[Exception] = []
    for zero_center in (True, False):
        candidate = combined.copy()
        try:
            sc.pp.scale(candidate, zero_center=zero_center, max_value=10)
            values = candidate.X.data if sparse.issparse(candidate.X) else candidate.X
            if not np.isfinite(values).all():
                raise ProcessingError("Scaling produced non-finite expression values")
            return candidate
        except Exception as exc:
            errors.append(exc)
            strategy = "Zero-centered" if zero_center else "Non-zero-centered"
            logger.warning("%s scaling failed: %s", strategy, exc)

    raise ProcessingError(
        "Data scaling failed completely. "
        f"Zero-center error: {errors[0]}. "
        f"Non-zero-center error: {errors[1]}. "
        "This usually indicates extreme outliers or invalid values. "
        "Consider additional quality control or outlier removal."
    ) from errors[-1]


def _compute_integration_pca(combined: ad.AnnData, n_pcs: int) -> None:
    """Validate the scaled matrix and compute PCA with ordered fallbacks."""
    max_components = min(n_pcs, combined.n_vars, combined.n_obs - 1)
    if max_components < 2:
        raise DataError(
            f"Cannot perform PCA: only {max_components} components possible. "
            f"Dataset has {combined.n_obs} cells and {combined.n_vars} genes. "
            "Minimum 2 components required for downstream analysis."
        )

    values = combined.X.data if sparse.issparse(combined.X) else combined.X
    if np.isnan(values).any():
        raise DataError("Data contains NaN values after scaling")
    if np.isinf(values).any():
        raise DataError("Data contains infinite values after scaling")

    for solver, component_cap in (("arpack", 50), ("randomized", 50), ("full", 20)):
        try:
            sc.tl.pca(
                combined,
                n_comps=min(max_components, component_cap),
                svd_solver=solver,
                zero_center=False,
            )
            return
        except Exception as exc:
            logger.warning("PCA with %s solver failed: %s", solver, exc)

    raise ProcessingError(
        f"PCA failed for {combined.n_obs}×{combined.n_vars} data. Check data quality."
    )


def _prepare_classical_integration(combined: ad.AnnData, n_pcs: int) -> ad.AnnData:
    """Select HVGs, remove constant genes, scale, and compute PCA."""
    if "highly_variable" in combined.var:
        hvg_mask = combined.var["highly_variable"].astype(bool).to_numpy()
        if not hvg_mask.any():
            raise DataError(
                "No highly variable genes found. Check HVG selection parameters."
            )
        combined = combined[:, hvg_mask]

    nonzero_variance = _nonzero_variance_mask(combined.X)
    if not nonzero_variance.all():
        logger.warning(
            "Removing %d genes with zero variance before scaling",
            int((~nonzero_variance).sum()),
        )
        combined = combined[:, nonzero_variance]

    combined = _scale_integration_data(combined)
    _compute_integration_pca(combined, n_pcs)
    return combined


def _integrate_harmony(combined: ad.AnnData, batch_key: str) -> None:
    """Apply Harmony and construct the corrected neighbor graph."""
    harmonypy = require("harmonypy", feature="Harmony integration")
    try:
        harmony_out = harmonypy.run_harmony(
            data_mat=combined.obsm["X_pca"],
            meta_data=pd.DataFrame({batch_key: combined.obs[batch_key].values}),
            vars_use=[batch_key],
            max_iter_harmony=10,
            verbose=True,
        )
        combined.obsm["X_pca_harmony"] = _orient_harmony_embedding(
            harmony_out.Z_corr, combined.n_obs
        )
        sc.pp.neighbors(combined, use_rep="X_pca_harmony")
    except Exception as exc:
        raise ProcessingError(
            f"Harmony integration failed: {exc}. "
            f"Check batch_key '{batch_key}' has ≥2 valid batches."
        ) from exc


def _integrate_scanorama(combined: ad.AnnData, batch_key: str, n_pcs: int) -> None:
    """Apply the Scanpy wrapper or its row-order-preserving raw fallback."""
    scanorama = require("scanorama", feature="Scanorama integration")
    try:
        import scanpy.external as sce

        if hasattr(sce.pp, "scanorama_integrate"):
            sce.pp.scanorama_integrate(
                combined,
                key=batch_key,
                basis="X_pca",
                adjusted_basis="X_scanorama",
            )
        else:
            batch_indices = [
                np.flatnonzero((combined.obs[batch_key] == batch).to_numpy())
                for batch in combined.obs[batch_key].unique()
            ]
            datasets = [combined[index].X for index in batch_indices]
            genes = [combined[index].var_names.tolist() for index in batch_indices]
            integrated, _ = scanorama.integrate(datasets, genes, dimred=n_pcs)
            combined.obsm["X_scanorama"] = _reassemble_scanorama_embeddings(
                batch_indices, integrated, combined.n_obs
            )
        sc.pp.neighbors(combined, use_rep="X_scanorama")
    except Exception as exc:
        raise ProcessingError(
            f"Scanorama integration failed: {exc}. Check gene overlap between batches."
        ) from exc


def _apply_classical_batch_correction(
    combined: ad.AnnData,
    *,
    method: str,
    batch_key: str,
    n_pcs: int,
) -> None:
    """Dispatch one classical integration backend."""
    if method == "harmony":
        _integrate_harmony(combined, batch_key)
    elif method == "bbknn":
        bbknn = require("bbknn", feature="BBKNN integration")
        bbknn.bbknn(combined, batch_key=batch_key, neighbors_within_batch=3)
    elif method == "scanorama":
        _integrate_scanorama(combined, batch_key, n_pcs)
    else:
        logger.warning(
            "Integration method '%s' not recognized. Using uncorrected PCA embedding.",
            method,
        )
        sc.pp.neighbors(combined)


def _integration_results_keys(
    combined: ad.AnnData, method: str
) -> dict[str, list[str]]:
    """Describe the representation produced by an integration method."""
    if method == "harmony":
        key = "X_pca_harmony" if "X_pca_harmony" in combined.obsm else "X_harmony"
        return {"obsm": [key]}
    if method == "bbknn":
        return {}
    if method == "scanorama":
        return {"obsm": ["X_scanorama"]}
    return {"obsm": ["X_pca"]}


def _store_integration_metadata(
    combined: ad.AnnData,
    *,
    method: str,
    batch_key: str,
    parameters: dict,
    results_keys: dict[str, list[str]],
) -> None:
    """Store H5AD-safe integration provenance."""
    batch_sizes = {
        str(key): int(value)
        for key, value in combined.obs[batch_key].value_counts().to_dict().items()
    }
    store_analysis_metadata(
        combined,
        analysis_name=f"integration_{method}",
        method=method,
        parameters=parameters,
        results_keys=results_keys,
        statistics={
            "n_batches": int(combined.obs[batch_key].nunique()),
            "batch_sizes": batch_sizes,
            "n_cells_total": int(combined.n_obs),
            "n_genes": int(combined.n_vars),
        },
    )


def integrate_multiple_samples(
    adatas: Union[list[ad.AnnData], ad.AnnData],
    batch_key: str = "batch",
    method: str = "harmony",
    n_pcs: int = 30,
    params: Optional[IntegrationParameters] = None,
    sample_keys: Sequence[str] | None = None,
) -> ad.AnnData:
    """Integrate multiple spatial transcriptomics samples.

    This function expects preprocessed data (normalized, log-transformed, with HVGs marked).
    Use preprocessing.py or preprocess_data() before calling this function.

    Args:
        adatas: List of preprocessed AnnData objects or a single combined AnnData.
        batch_key: Batch information key in obs.
        method: Integration method ('harmony', 'bbknn', 'scanorama', 'scvi').
        n_pcs: Number of principal components for integration.
        params: Optional IntegrationParameters for method-specific settings.
        sample_keys: Stable source identifiers used to make cell IDs globally unique.

    Returns:
        Integrated AnnData with batch correction applied.

    Raises:
        ParameterError: If fewer than 2 datasets provided for integration.
        DataError: If data is not properly preprocessed.
    """

    combined = _prepare_combined_input(
        adatas,
        batch_key=batch_key,
        sample_keys=sample_keys,
    )
    _ensure_integration_hvgs(combined, batch_key)

    # NOTE: Do NOT set combined.raw here if it is None.
    # Input data is already normalized+log-transformed (see docstring).
    # Storing it in .raw would violate the contract that .raw holds raw counts,
    # poisoning downstream tools (differential, deconvolution, annotation) that
    # rely on get_raw_data_source() treating .raw as the highest-priority count source.
    # If .raw or layers["counts"] was set during preprocessing, it's already present;
    # if not, integration should not fabricate one from normalized data.
    if combined.raw is None and "counts" not in combined.layers:
        logger.warning(
            "No raw counts found (adata.raw or layers['counts']). "
            "Downstream analyses requiring raw counts may be limited. "
            "Ensure preprocess_data() was run before integration."
        )

    # ========================================================================
    # EARLY BRANCH FOR scVI-TOOLS METHODS
    # scVI requires normalized+log data WITHOUT scaling/PCA
    # It generates its own latent representation
    # NOTE: scVI-tools methods work better with ALL genes, not just HVGs
    # ========================================================================
    if method == "scvi":
        return _run_scvi_integration(combined, batch_key=batch_key, params=params)

    # ========================================================================
    # CLASSICAL METHODS: Continue with scale → PCA → integration
    # ========================================================================

    combined = _prepare_classical_integration(combined, n_pcs)

    _apply_classical_batch_correction(
        combined,
        method=method,
        batch_key=batch_key,
        n_pcs=n_pcs,
    )
    sc.tl.umap(combined)
    _store_integration_metadata(
        combined,
        method=method,
        batch_key=batch_key,
        parameters={
            "batch_key": batch_key,
            "n_pcs": n_pcs,
            "n_batches": int(combined.obs[batch_key].nunique()),
        },
        results_keys=_integration_results_keys(combined, method),
    )

    return combined


def rescale_spatial_coordinates(
    combined_adata: ad.AnnData,
    batch_key: str = "batch",
    reference_batch: Optional[str] = None,
) -> ad.AnnData:
    """Rescale spatial coordinates to a common scale across batches.

    Applies z-score standardization using the reference batch's statistics,
    so all batches share comparable coordinate ranges. This is NOT geometric
    registration (rotation/translation) -- use spatial_registration tools
    for that.

    Args:
        combined_adata: Combined AnnData containing multiple samples.
        batch_key: Batch information key in obs.
        reference_batch: Reference batch for rescaling. String values are also
            matched against non-string batch labels. Uses the first batch if None.

    Returns:
        AnnData with rescaled spatial coordinates in obsm['spatial_aligned'].

    Raises:
        DataError: If the dataset, batch labels, or spatial coordinates are invalid.
        ParameterError: If the batch column or reference batch is not found.
    """
    from sklearn.preprocessing import StandardScaler

    if batch_key not in combined_adata.obs:
        raise ParameterError(
            f"Batch key '{batch_key}' not found in adata.obs. "
            f"Available columns: {list(combined_adata.obs.columns)}"
        )
    if combined_adata.n_obs == 0:
        raise DataError("Dataset is empty, cannot perform spatial registration")

    batch_values = combined_adata.obs[batch_key]
    if batch_values.isna().any():
        raise DataError(f"Batch column '{batch_key}' contains missing values")

    spatial_coords = require_spatial_coords(combined_adata)

    # Get batch information
    batches = list(batch_values.unique())
    batches_by_label = {str(batch): batch for batch in batches}

    # If reference batch not specified, use the first batch
    if reference_batch is None:
        reference_batch_value = batches[0]
    elif reference_batch in batches:
        reference_batch_value = reference_batch
    elif reference_batch in batches_by_label:
        reference_batch_value = batches_by_label[reference_batch]
    else:
        raise ParameterError(f"Reference batch '{reference_batch}' not found in data")
    reference_batch_label = str(reference_batch_value)

    # Get reference batch spatial coordinates
    reference_mask = (batch_values == reference_batch_value).to_numpy()
    ref_coords = spatial_coords[reference_mask]

    # Standardize reference coordinates
    scaler = StandardScaler()
    scaler.fit(ref_coords)

    # Transform and restore every batch directly in original observation order.
    aligned_coords = np.empty(spatial_coords.shape, dtype=float)
    for batch in batches:
        batch_mask = (batch_values == batch).to_numpy()
        aligned_coords[batch_mask] = scaler.transform(spatial_coords[batch_mask])

    combined_adata.obsm["spatial_aligned"] = aligned_coords

    # Store metadata for scientific provenance tracking
    n_batches = len(batches)
    # Convert keys to strings for H5AD compatibility (mirrors integration metadata)
    batch_sizes = {
        str(batch): int(np.sum(combined_adata.obs[batch_key] == batch))
        for batch in batches
    }

    store_analysis_metadata(
        combined_adata,
        analysis_name="spatial_alignment",
        method="standardization",
        parameters={
            "batch_key": batch_key,
            "reference_batch": reference_batch_label,
        },
        results_keys={"obsm": ["spatial_aligned"]},
        statistics={
            "n_batches": n_batches,
            "batch_sizes": batch_sizes,
            "reference_batch": reference_batch_label,
        },
    )

    return combined_adata


def integrate_with_scvi(
    combined: sc.AnnData,
    batch_key: str = "batch",
    n_hidden: int = 128,
    n_latent: int = 10,
    n_layers: int = 1,
    dropout_rate: float = 0.1,
    gene_likelihood: str = "zinb",
    n_epochs: Optional[int] = None,
    use_gpu: bool = False,
) -> sc.AnnData:
    """Integrate data using scVI for batch correction

    scVI is a deep generative model for single-cell RNA-seq that can perform
    batch correction by learning a low-dimensional latent representation.

    Args:
        combined: Combined AnnData object with multiple batches
        batch_key: Column name in obs containing batch labels
        n_hidden: Number of nodes per hidden layer (default: 128)
        n_latent: Dimensionality of the latent space (default: 10)
        n_layers: Number of hidden layers (default: 1)
        dropout_rate: Dropout rate for neural networks (default: 0.1)
        gene_likelihood: Distribution for gene expression (default: "zinb")
        n_epochs: Number of training epochs (None = auto-determine)
        use_gpu: Whether to use GPU acceleration (default: False)

    Returns:
        AnnData object with scVI latent representation in obsm['X_scvi']

    Raises:
        DependencyError: If scvi-tools is not installed or cannot be imported
        ValueError: If data is not preprocessed or invalid

    Reference:
        Lopez et al. (2018) "Deep generative modeling for single-cell transcriptomics"
        Nature Methods 15, 1053–1058
    """
    scvi = require("scvi", feature="scVI integration")

    # Validate data is preprocessed (HVG selection uses normalized X)
    max_val = combined.X.max() if hasattr(combined.X, "max") else np.max(combined.X)
    if max_val > 50:
        raise DataError(
            f"scVI requires preprocessed data. Max value {max_val:.1f} too high."
        )

    # Validate batch key
    if batch_key not in combined.obs:
        raise ParameterError(
            f"Batch key '{batch_key}' not found in adata.obs. "
            f"Available columns: {list(combined.obs.columns)}"
        )

    # Check for batch diversity
    n_batches = combined.obs[batch_key].nunique()
    if n_batches < 2:
        raise DataError(
            f"scVI requires at least 2 batches, found {n_batches}. "
            "Check your batch labels."
        )

    # scVI's generative model requires raw counts, not log-normalized data.
    # Use layers["counts"] when available; fall back to adata.X only if
    # it appears to contain integer counts.
    layer_for_scvi: str | None = None
    if "counts" in combined.layers:
        is_int_counts, has_neg_counts, _ = check_is_integer_counts(
            combined.layers["counts"]
        )
        if is_int_counts and not has_neg_counts:
            layer_for_scvi = "counts"
        else:
            raise DataError(
                "layers['counts'] exists but does not contain valid "
                "integer counts (found normalized or negative values). "
                "scVI requires raw integer counts."
            )
    else:
        is_int, has_neg, _ = check_is_integer_counts(combined.X)
        if is_int and not has_neg:
            layer_for_scvi = None  # X is already counts
        else:
            # Try to salvage counts from .raw before giving up
            from ..utils.adata_utils import ensure_counts_layer
            from ..utils.exceptions import DataNotFoundError as _DNF

            try:
                ensure_counts_layer(combined)
                layer_for_scvi = "counts"
            except _DNF:
                raise DataError(
                    "scVI requires raw count data but only normalized "
                    "values are available (no integer counts in X, "
                    "layers['counts'], or .raw).\n\n"
                    "Solutions:\n"
                    "1. Load data with raw counts before integration\n"
                    "2. Use method='harmony' or method='scanorama' "
                    "which work with normalized data"
                ) from None

    # Setup AnnData for scVI
    scvi.model.SCVI.setup_anndata(
        combined,
        batch_key=batch_key,
        layer=layer_for_scvi,
    )

    # Initialize scVI model
    model = scvi.model.SCVI(
        combined,
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers,
        dropout_rate=dropout_rate,
        gene_likelihood=gene_likelihood,
    )

    # Auto-determine epochs based on dataset size if not specified
    if n_epochs is None:
        n_cells = combined.n_obs
        if n_cells < 1000:
            n_epochs = 400
        elif n_cells < 10000:
            n_epochs = 200
        else:
            n_epochs = 100

    # Train model
    accelerator = get_accelerator(prefer_gpu=use_gpu)
    model.train(max_epochs=n_epochs, early_stopping=True, accelerator=accelerator)

    # Get latent representation
    combined.obsm["X_scvi"] = model.get_latent_representation()

    # Compute neighbors using scVI embedding
    sc.pp.neighbors(combined, use_rep="X_scvi")

    return combined


async def integrate_samples(
    data_ids: list[str],
    ctx: "ToolContext",
    params: IntegrationParameters | None = None,
) -> IntegrationResult:
    """Integrate multiple spatial transcriptomics samples and perform batch correction

    Args:
        data_ids: List of dataset IDs to integrate
        ctx: ToolContext for unified data access and logging
        params: Integration parameters

    Returns:
        Integration result
    """
    if params is None:
        params = IntegrationParameters()
    if len(data_ids) < 2:
        raise ParameterError("Integration requires at least 2 dataset IDs")
    if len(set(data_ids)) != len(data_ids):
        raise ParameterError("Integration dataset IDs must be distinct")

    # Collect all AnnData objects
    # Memory optimization: concatenate() creates new object without modifying sources
    # Verified by comprehensive testing: all operations preserve original datasets
    # Users can still access A, B, C after integration via ctx references
    adatas = []
    for data_id in data_ids:
        adata = await ctx.get_adata(data_id)
        adatas.append(adata)

    # Integrate samples (pass full params for method-specific settings like scVI)
    combined_adata = await run_sync(
        integrate_multiple_samples,
        adatas,
        batch_key=params.batch_key,
        method=params.method,
        n_pcs=params.n_pcs,
        params=params,
        sample_keys=data_ids,
    )

    # Rescale spatial coordinates if requested and available
    # Note: Spatial rescaling is optional - BBKNN, Harmony, MNN, Scanorama
    # work on gene expression/PCA space without spatial coordinates
    if params.align_spatial and get_spatial_key(combined_adata):
        combined_adata = await run_sync(
            rescale_spatial_coordinates,
            combined_adata,
            batch_key=params.batch_key,
            reference_batch=params.reference_batch,
        )

    # Reserve an ID for reproducibility exports, but do not publish the dataset
    # until every required export has completed. A failed export must not leave
    # an orphaned dataset whose ID was never returned to the client.
    integrated_id = ctx.reserve_dataset_id(prefix="integrated")

    # Export results for reproducibility
    # Note: Metadata was stored in helper functions; export uses the appropriate analysis names
    if params.method == "scvi":
        await run_sync(
            export_analysis_result,
            combined_adata,
            integrated_id,
            "integration_scvi",
        )
    else:
        await run_sync(
            export_analysis_result,
            combined_adata,
            integrated_id,
            f"integration_{params.method}",
        )

    if params.align_spatial and "spatial_aligned" in combined_adata.obsm:
        await run_sync(
            export_analysis_result,
            combined_adata,
            integrated_id,
            "spatial_alignment",
        )

    # The embedding key is what downstream tools need to consume the result,
    # so report it instead of leaving the caller to guess the method's naming.
    embedding_keys = _integration_results_keys(combined_adata, params.method).get(
        "obsm", []
    )
    result = IntegrationResult(
        data_id=integrated_id,
        n_samples=len(data_ids),
        integration_method=params.method,
        n_cells=int(combined_adata.n_obs),
        batch_key=params.batch_key,
        embedding_key=embedding_keys[0] if embedding_keys else None,
    )
    await ctx.add_dataset(
        combined_adata,
        prefix="integrated",
        data_id=integrated_id,
    )
    return result
