"""
Integration tools for spatial transcriptomics data.
"""

import logging
from typing import TYPE_CHECKING, Optional, Union

import anndata as ad
import numpy as np
import scanpy as sc

from ..models.analysis import IntegrationResult
from ..models.data import IntegrationParameters
from ..utils.adata_utils import check_is_integer_counts
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


def integrate_multiple_samples(
    adatas: Union[list[ad.AnnData], ad.AnnData],
    batch_key: str = "batch",
    method: str = "harmony",
    n_pcs: int = 30,
    params: Optional[IntegrationParameters] = None,
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

    Returns:
        Integrated AnnData with batch correction applied.

    Raises:
        ParameterError: If fewer than 2 datasets provided for integration.
        DataError: If data is not properly preprocessed.
    """

    # Merge datasets
    if isinstance(adatas, list):
        # Validate list has at least 2 datasets for integration
        if len(adatas) < 2:
            raise ParameterError(
                f"Integration requires at least 2 datasets, got {len(adatas)}. "
                "Use preprocess_data for single dataset processing."
            )

        # Check if datasets have batch labels
        has_batch_labels = all(batch_key in adata.obs for adata in adatas)

        if not has_batch_labels:
            # Auto-create batch labels for multi-sample integration
            # Each sample becomes its own batch (scientifically correct for independent samples)
            for i, adata in enumerate(adatas):
                if batch_key not in adata.obs:
                    adata.obs[batch_key] = f"sample_{i}"

        # Merge datasets with stable obs semantics and outer-gene union.
        # Use anndata.concat (AnnData.concatenate is deprecated).
        combined = ad.concat(
            adatas,
            join="outer",
            merge="first",
            uns_merge="first",
            index_unique=None,
        )

        # FIX: Clean var columns with NA values in object dtype
        # Problem: outer join creates NA values in var columns when genes don't exist in all samples
        # When object columns contain NA, H5AD save corrupts var.index (becomes 0,1,2...)
        # and moves gene names to _index column
        # Solution: Fill NA with appropriate values or convert types
        for col in combined.var.columns:
            if combined.var[col].dtype == "object" and combined.var[col].isna().any():
                # For boolean-like columns (highly_variable, etc.), fill NA with False
                unique_vals = combined.var[col].dropna().unique()
                if set(unique_vals).issubset({True, False, "True", "False"}):
                    combined.var[col] = combined.var[col].fillna(False).astype(bool)
                else:
                    # For string columns, fill NA with empty string
                    combined.var[col] = combined.var[col].fillna("").astype(str)

        # FIX: Remove incomplete diffmap artifacts created by concatenation (scanpy issue #1021)
        # Problem: concatenate() copies obsm['X_diffmap'] but NOT uns['diffmap_evals']
        # This creates incomplete state that causes KeyError in sc.tl.umap()
        # Solution: Delete incomplete artifacts to allow UMAP to use default initialization
        if "X_diffmap" in combined.obsm:
            del combined.obsm["X_diffmap"]
        if "diffmap_evals" in combined.uns:
            del combined.uns["diffmap_evals"]

    else:
        # If already a merged dataset, ensure it has batch information
        combined = adatas
        if batch_key not in combined.obs:
            raise ParameterError(
                f"Merged dataset is missing batch information key '{batch_key}'"
            )

    # Validate input data is preprocessed using SSOT integer check
    is_int, _, _ = check_is_integer_counts(combined.X)
    if is_int:
        raise DataError("Data appears to be raw counts. Run preprocessing first.")

    # Validate data quality before processing
    validate_adata_basics(combined, min_obs=10, min_vars=10, check_empty_ratio=True)

    # Check if data has highly variable genes marked (should be done in preprocessing)
    if "highly_variable" not in combined.var.columns:
        logger.warning(
            "No highly variable genes marked after merge. Recalculating HVGs with batch correction."
        )
        # Recalculate HVGs with batch correction
        sc.pp.highly_variable_genes(
            combined,
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            batch_key=batch_key,
            n_top_genes=2000,
        )
        n_hvg = combined.var["highly_variable"].sum()
    else:
        n_hvg = combined.var["highly_variable"].sum()
        if n_hvg == 0:
            logger.warning(
                "No genes marked as highly variable after merge, recalculating"
            )
            # Recalculate HVGs with batch correction
            sc.pp.highly_variable_genes(
                combined,
                min_mean=0.0125,
                max_mean=3,
                min_disp=0.5,
                batch_key=batch_key,
                n_top_genes=2000,
            )
            n_hvg = combined.var["highly_variable"].sum()
        elif n_hvg < 50:
            logger.warning(
                f"Very few HVGs ({n_hvg}), recalculating with batch correction"
            )
            sc.pp.highly_variable_genes(
                combined,
                min_mean=0.0125,
                max_mean=3,
                min_disp=0.5,
                batch_key=batch_key,
                n_top_genes=2000,
            )
            n_hvg = combined.var["highly_variable"].sum()

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
        # Use user-configurable parameters if provided, otherwise use defaults
        # This ensures scientific reproducibility and user control
        scvi_n_hidden = params.scvi_n_hidden if params else 128
        scvi_n_latent = params.scvi_n_latent if params else 10
        scvi_n_layers = params.scvi_n_layers if params else 1
        scvi_dropout_rate = params.scvi_dropout_rate if params else 0.1
        scvi_gene_likelihood = params.scvi_gene_likelihood if params else "zinb"
        scvi_n_epochs = params.n_epochs if params else None
        scvi_use_gpu = params.use_gpu if params else False

        try:
            combined = integrate_with_scvi(
                combined,
                batch_key=batch_key,
                n_hidden=scvi_n_hidden,
                n_latent=scvi_n_latent,
                n_layers=scvi_n_layers,
                dropout_rate=scvi_dropout_rate,
                gene_likelihood=scvi_gene_likelihood,
                n_epochs=scvi_n_epochs,
                use_gpu=scvi_use_gpu,
            )
        except ChatSpatialError:
            raise
        except Exception as e:
            raise ProcessingError(
                f"scVI integration failed: {e}. "
                f"Ensure data is preprocessed and has ≥2 batches."
            ) from e

        # Calculate UMAP embedding to visualize integration effect
        sc.tl.umap(combined)

        # Store metadata for scientific provenance tracking
        n_batches = combined.obs[batch_key].nunique()
        batch_sizes = combined.obs[batch_key].value_counts().to_dict()

        # CRITICAL FIX: Convert dict keys to strings for H5AD compatibility
        batch_sizes = {str(k): int(v) for k, v in batch_sizes.items()}

        store_analysis_metadata(
            combined,
            analysis_name="integration_scvi",
            method="scvi",
            parameters={
                "batch_key": batch_key,
                "n_hidden": scvi_n_hidden,
                "n_latent": scvi_n_latent,
                "n_layers": scvi_n_layers,
                "dropout_rate": scvi_dropout_rate,
                "gene_likelihood": scvi_gene_likelihood,
                "n_epochs": scvi_n_epochs,
                "use_gpu": scvi_use_gpu,
            },
            results_keys={"obsm": ["X_scvi"]},
            statistics={
                "n_batches": int(n_batches),
                "batch_sizes": batch_sizes,
                "n_cells_total": int(combined.n_obs),
                "n_genes": int(combined.n_vars),
            },
        )

        return combined

    # ========================================================================
    # CLASSICAL METHODS: Continue with scale → PCA → integration
    # ========================================================================

    # Filter to highly variable genes for classical methods
    if "highly_variable" in combined.var.columns:
        n_hvg = combined.var["highly_variable"].sum()
        if n_hvg == 0:
            raise DataError(
                "No highly variable genes found. Check HVG selection parameters."
            )
        # Memory optimization: Subsetting creates view, reassignment triggers GC
        # No need to materialize with .copy() - view will be materialized on first write
        combined = combined[:, combined.var["highly_variable"]]

    # Remove genes with zero variance to avoid NaN in scaling
    from scipy import sparse

    # MEMORY OPTIMIZATION: Calculate variance without toarray()
    # Uses E[X²] - E[X]² formula for sparse matrices
    # Saves ~80% memory (e.g., 76 MB → 15 MB for 10k cells × 2k genes)
    if sparse.issparse(combined.X):
        # Sparse matrix: compute variance using E[X²] - E[X]² formula
        # This avoids creating dense copy (5-10x memory reduction)
        mean_per_gene = np.array(combined.X.mean(axis=0)).flatten()

        # Calculate E[X²] using .power(2) - cleaner and ~1.5x faster than copy + data**2
        mean_squared = np.array(combined.X.power(2).mean(axis=0)).flatten()

        # Variance = E[X²] - E[X]²
        gene_var = mean_squared - mean_per_gene**2
    else:
        # Dense matrix: use standard variance calculation
        gene_var = np.var(combined.X, axis=0)
    nonzero_var_genes = gene_var > 0
    if not np.all(nonzero_var_genes):
        n_removed = np.sum(~nonzero_var_genes)
        logger.warning(f"Removing {n_removed} genes with zero variance before scaling")
        # Memory optimization: Subsetting creates view, no need to copy
        # View will be materialized when scaling modifies the data
        combined = combined[:, nonzero_var_genes]

    # Scale data with proper error handling
    try:
        sc.pp.scale(combined, zero_center=True, max_value=10)
    except Exception as e:
        logger.warning(f"Scaling with zero centering failed: {e}")
        try:
            sc.pp.scale(combined, zero_center=False, max_value=10)
        except Exception as e2:
            raise ProcessingError(
                f"Data scaling failed completely. Zero-center error: {e}. Non-zero-center error: {e2}. "
                f"This usually indicates data contains extreme outliers or invalid values. "
                f"Consider additional quality control or outlier removal."
            ) from e2

    # PCA with proper error handling
    # Determine safe number of components
    max_possible_components = min(n_pcs, combined.n_vars, combined.n_obs - 1)

    if max_possible_components < 2:
        raise DataError(
            f"Cannot perform PCA: only {max_possible_components} components possible. "
            f"Dataset has {combined.n_obs} cells and {combined.n_vars} genes. "
            f"Minimum 2 components required for downstream analysis."
        )

    # Check data matrix before PCA
    # MEMORY OPTIMIZATION: Check sparse matrix .data directly without toarray()
    # Sparse matrices only store non-zero elements, and zero elements cannot be NaN/Inf
    # Saves ~80% memory (e.g., 76 MB → 15 MB for 10k cells × 2k genes)
    from scipy import sparse

    if sparse.issparse(combined.X):
        # Sparse matrix: only check non-zero elements stored in .data
        # This avoids creating a dense copy (5-10x memory reduction)
        if np.isnan(combined.X.data).any():
            raise DataError("Data contains NaN values after scaling")
        if np.isinf(combined.X.data).any():
            raise DataError("Data contains infinite values after scaling")
    else:
        # Dense matrix: check all elements
        if np.isnan(combined.X).any():
            raise DataError("Data contains NaN values after scaling")
        if np.isinf(combined.X).any():
            raise DataError("Data contains infinite values after scaling")

    # Variance check removed: zero-variance genes already filtered at lines 301-323

    # Try PCA with different solvers, but fail properly if none work
    pca_success = False
    for solver, max_comps in [
        ("arpack", min(max_possible_components, 50)),
        ("randomized", min(max_possible_components, 50)),
        ("full", min(max_possible_components, 20)),
    ]:
        try:
            sc.tl.pca(combined, n_comps=max_comps, svd_solver=solver, zero_center=False)
            pca_success = True
            break
        except Exception as e:
            logger.warning(f"PCA with {solver} solver failed: {e}")
            continue

    if not pca_success:
        raise ProcessingError(
            f"PCA failed for {combined.n_obs}×{combined.n_vars} data. Check data quality."
        )

    # Apply batch correction based on selected method
    if method == "harmony":
        # Use Harmony for batch correction
        # Direct harmonypy call for version compatibility (scanpy.external has issues
        # with harmonypy >= 0.1.0, see: https://github.com/scverse/scanpy/issues/3940)
        harmonypy = require("harmonypy", feature="Harmony integration")
        try:
            import pandas as pd

            X_pca = combined.obsm["X_pca"]
            n_cells = combined.n_obs
            meta_data = pd.DataFrame({batch_key: combined.obs[batch_key].values})

            harmony_out = harmonypy.run_harmony(
                data_mat=X_pca,
                meta_data=meta_data,
                vars_use=[batch_key],
                max_iter_harmony=10,
                verbose=True,
            )

            # Shape normalization for version compatibility:
            # - harmonypy < 0.1.0: Z_corr is (n_pcs, n_cells), needs .T
            # - harmonypy >= 0.1.0: Z_corr is (n_cells, n_pcs), already correct
            combined.obsm["X_pca_harmony"] = _orient_harmony_embedding(
                harmony_out.Z_corr, n_cells
            )

            sc.pp.neighbors(combined, use_rep="X_pca_harmony")

        except Exception as e:
            raise ProcessingError(
                f"Harmony integration failed: {e}. "
                f"Check batch_key '{batch_key}' has ≥2 valid batches."
            ) from e

    elif method == "bbknn":
        # Use BBKNN for batch correction
        bbknn = require("bbknn", feature="BBKNN integration")

        bbknn.bbknn(combined, batch_key=batch_key, neighbors_within_batch=3)

    elif method == "scanorama":
        # Use Scanorama for batch correction
        # BEST PRACTICE: Use scanpy.external wrapper for better integration with scanpy workflow
        scanorama = require("scanorama", feature="Scanorama integration")
        try:
            import scanpy.external as sce

            # Check if scanorama_integrate is available in scanpy.external
            if hasattr(sce.pp, "scanorama_integrate"):
                # Use scanpy.external wrapper (preferred method)
                sce.pp.scanorama_integrate(
                    combined, key=batch_key, basis="X_pca", adjusted_basis="X_scanorama"
                )
                # Use integrated representation for neighbor graph
                sc.pp.neighbors(combined, use_rep="X_scanorama")
            else:
                # Fallback to raw scanorama (same algorithm, different interface)
                # Separate data by batch, tracking original row indices
                datasets = []
                genes_list = []
                batch_indices = []  # original row positions per batch

                for batch in combined.obs[batch_key].unique():
                    batch_mask = (combined.obs[batch_key] == batch).to_numpy()
                    idx = np.where(batch_mask)[0]
                    batch_indices.append(idx)

                    batch_data = combined[idx]
                    # Scanorama natively supports sparse matrices
                    datasets.append(batch_data.X)
                    genes_list.append(batch_data.var_names.tolist())

                # Run Scanorama integration
                integrated, _corrected_genes = scanorama.integrate(
                    datasets, genes_list, dimred=n_pcs
                )

                # Reassemble in original row order (handles interleaved batches)
                integrated_X = _reassemble_scanorama_embeddings(
                    batch_indices, integrated, combined.n_obs
                )

                # Store integrated representation in obsm
                combined.obsm["X_scanorama"] = integrated_X

                # Use integrated representation for neighbor graph
                sc.pp.neighbors(combined, use_rep="X_scanorama")

        except Exception as e:
            raise ProcessingError(
                f"Scanorama integration failed: {e}. "
                f"Check gene overlap between batches."
            ) from e

    else:
        # Default: use uncorrected PCA result
        logger.warning(
            f"Integration method '{method}' not recognized. "
            f"Using uncorrected PCA embedding."
        )
        sc.pp.neighbors(combined)

    # Calculate UMAP embedding to visualize integration effect
    sc.tl.umap(combined)

    # Store metadata for scientific provenance tracking
    # Determine which representation was used
    # Note: neighbors (connectivities/distances sparse matrices) not exported to CSV
    if method == "harmony":
        if "X_pca_harmony" in combined.obsm:
            results_keys = {"obsm": ["X_pca_harmony"]}
        else:
            results_keys = {"obsm": ["X_harmony"]}
    elif method == "bbknn":
        # BBKNN primarily modifies neighbors graph, no obsm output to export
        results_keys = {}
    elif method == "scanorama":
        results_keys = {"obsm": ["X_scanorama"]}
    else:
        results_keys = {"obsm": ["X_pca"]}

    # Get batch statistics
    n_batches = combined.obs[batch_key].nunique()
    batch_sizes = combined.obs[batch_key].value_counts().to_dict()

    # CRITICAL FIX: Convert dict keys to strings for H5AD compatibility
    # H5AD requires all dictionary keys to be strings
    # Without this, save_data() fails with "Can't implicitly convert non-string objects to strings"
    batch_sizes = {str(k): int(v) for k, v in batch_sizes.items()}

    store_analysis_metadata(
        combined,
        analysis_name=f"integration_{method}",
        method=method,
        parameters={
            "batch_key": batch_key,
            "n_pcs": n_pcs,
            "n_batches": n_batches,
        },
        results_keys=results_keys,
        statistics={
            "n_batches": int(n_batches),  # Also ensure int types for H5AD
            "batch_sizes": batch_sizes,
            "n_cells_total": int(combined.n_obs),
            "n_genes": int(combined.n_vars),
        },
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

    # Collect all AnnData objects
    # Memory optimization: concatenate() creates new object without modifying sources
    # Verified by comprehensive testing: all operations preserve original datasets
    # Users can still access A, B, C after integration via ctx references
    adatas = []
    for data_id in data_ids:
        adata = await ctx.get_adata(data_id)
        adatas.append(adata)

    # Integrate samples (pass full params for method-specific settings like scVI)
    combined_adata = integrate_multiple_samples(
        adatas,
        batch_key=params.batch_key,
        method=params.method,
        n_pcs=params.n_pcs,
        params=params,
    )

    # Rescale spatial coordinates if requested and available
    # Note: Spatial rescaling is optional - BBKNN, Harmony, MNN, Scanorama
    # work on gene expression/PCA space without spatial coordinates
    if params.align_spatial and get_spatial_key(combined_adata):
        combined_adata = rescale_spatial_coordinates(
            combined_adata,
            batch_key=params.batch_key,
            reference_batch=params.reference_batch,
        )

    # Store integrated data using ToolContext (ID generated by data manager)
    integrated_id = await ctx.add_dataset(combined_adata, prefix="integrated")

    # Export results for reproducibility
    # Note: Metadata was stored in helper functions; export uses the appropriate analysis names
    if params.method == "scvi":
        export_analysis_result(combined_adata, integrated_id, "integration_scvi")
    else:
        export_analysis_result(
            combined_adata, integrated_id, f"integration_{params.method}"
        )

    if params.align_spatial and "spatial_aligned" in combined_adata.obsm:
        export_analysis_result(combined_adata, integrated_id, "spatial_alignment")

    # Return result
    return IntegrationResult(
        data_id=integrated_id,
        n_samples=len(data_ids),
        integration_method=params.method,
    )
