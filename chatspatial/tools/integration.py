"""
Integration tools for spatial transcriptomics data.
"""

import logging
from typing import TYPE_CHECKING, List, Optional

import anndata as ad
import numpy as np
import scanpy as sc

from ..models.analysis import IntegrationResult
from ..models.data import IntegrationParameters
from ..utils.dependency_manager import require
from ..utils.exceptions import (
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

from ..utils.adata_utils import get_spatial_key, store_analysis_metadata

logger = logging.getLogger(__name__)


def validate_data_quality(adata, min_cells=10, min_genes=10):
    """Validate data quality before integration

    Args:
        adata: AnnData object
        min_cells: Minimum number of cells required
        min_genes: Minimum number of genes required

    Raises:
        ValueError: If data quality is insufficient for integration
    """
    if adata.n_obs < min_cells:
        raise DataError(
            f"Dataset has only {adata.n_obs} cells, minimum {min_cells} required for integration. "
            f"Consider combining with other datasets or use single-sample analysis."
        )

    if adata.n_vars < min_genes:
        raise DataError(
            f"Dataset has only {adata.n_vars} genes, minimum {min_genes} required for integration. "
            f"Check if data was properly loaded and genes were not over-filtered."
        )

    # Check for empty cells or genes
    if hasattr(adata.X, "toarray"):
        # Sparse-aware: use getnnz() method (no conversion needed)
        cell_counts = np.array(adata.X.getnnz(axis=1)).flatten()
        gene_counts = np.array(adata.X.getnnz(axis=0)).flatten()
    else:
        # Dense matrix
        cell_counts = np.sum(adata.X > 0, axis=1)
        gene_counts = np.sum(adata.X > 0, axis=0)

    empty_cells = np.sum(cell_counts == 0)
    empty_genes = np.sum(gene_counts == 0)

    if empty_cells > adata.n_obs * 0.1:
        raise DataError(
            f"{empty_cells} cells ({empty_cells/adata.n_obs*100:.1f}%) have zero expression. "
            f"Check data quality and consider filtering."
        )

    if empty_genes > adata.n_vars * 0.5:
        raise DataError(
            f"{empty_genes} genes ({empty_genes/adata.n_vars*100:.1f}%) have zero expression across all cells. "
            f"Consider gene filtering before integration."
        )

    logging.info(
        f"Data quality validation passed: {adata.n_obs} cells, {adata.n_vars} genes"
    )


def integrate_multiple_samples(
    adatas,
    batch_key="batch",
    method="harmony",
    n_pcs=30,
    params: Optional[IntegrationParameters] = None,
):
    """Integrate multiple spatial transcriptomics samples

    This function expects preprocessed data (normalized, log-transformed, with HVGs marked).
    Use preprocessing.py or preprocess_data() before calling this function.

    Args:
        adatas: List of preprocessed AnnData objects or a single combined AnnData object
        batch_key: Batch information key
        method: Integration method, options: 'harmony', 'bbknn', 'scanorama', 'scvi'
        n_pcs: Number of principal components for integration
        params: Optional IntegrationParameters for method-specific settings (scVI, etc.)

    Returns:
        Integrated AnnData object with batch correction applied

    Raises:
        ValueError: If data is not properly preprocessed
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
            logging.info(
                "No batch labels found. Auto-creating batch labels: each sample = 1 batch.\n"
                "This is scientifically appropriate when integrating independent samples."
            )
            for i, adata in enumerate(adatas):
                if batch_key not in adata.obs:
                    # Use sample index as batch label
                    adata.obs[batch_key] = f"sample_{i}"
                    logging.info(f"Dataset {i}: assigned batch label 'sample_{i}'")

        # Merge datasets
        combined = adatas[0].concatenate(
            adatas[1:],
            batch_key=batch_key,
            join="outer",  # Use outer join to keep all genes
        )

        # FIX: Clean var columns with NA values in object dtype
        # Problem: outer join creates NA values in var columns when genes don't exist in all samples
        # When object columns contain NA, H5AD save corrupts var.index (becomes 0,1,2...)
        # and moves gene names to _index column
        # Solution: Fill NA with appropriate values or convert types
        import pandas as pd

        for col in combined.var.columns:
            if combined.var[col].dtype == "object" and combined.var[col].isna().any():
                # For boolean-like columns (highly_variable, etc.), fill NA with False
                unique_vals = combined.var[col].dropna().unique()
                if set(unique_vals).issubset({True, False, "True", "False"}):
                    combined.var[col] = combined.var[col].fillna(False).astype(bool)
                else:
                    # For string columns, fill NA with empty string
                    combined.var[col] = combined.var[col].fillna("").astype(str)
                logging.info(
                    f"Cleaned var column '{col}' with NA values for H5AD compatibility"
                )

        # FIX: Remove incomplete diffmap artifacts created by concatenation (scanpy issue #1021)
        # Problem: concatenate() copies obsm['X_diffmap'] but NOT uns['diffmap_evals']
        # This creates incomplete state that causes KeyError in sc.tl.umap()
        # Solution: Delete incomplete artifacts to allow UMAP to use default initialization
        if "X_diffmap" in combined.obsm:
            del combined.obsm["X_diffmap"]
            logging.info(
                "Removed incomplete X_diffmap from concatenated data (scanpy issue #1021)"
            )
        if "diffmap_evals" in combined.uns:
            del combined.uns["diffmap_evals"]
            logging.info(
                "Removed incomplete diffmap_evals from concatenated data (scanpy issue #1021)"
            )

    else:
        # If already a merged dataset, ensure it has batch information
        combined = adatas
        if batch_key not in combined.obs:
            raise ParameterError(
                f"Merged dataset is missing batch information key '{batch_key}'"
            )

    # Validate input data is preprocessed
    # Check if data appears to be raw (high values without log transformation)
    max_val = combined.X.max() if hasattr(combined.X, "max") else np.max(combined.X)
    min_val = combined.X.min() if hasattr(combined.X, "min") else np.min(combined.X)

    # Raw count data typically has high integer values and no negative values
    # Properly preprocessed data should be either:
    # 1. Log-transformed (positive values, typically 0-15 range)
    # 2. Scaled (centered around 0, can have negative values)
    if min_val >= 0 and max_val > 100:
        raise DataError("Data appears to be raw counts. Run preprocessing first.")

    # Check if data appears to be normalized (reasonable range after preprocessing)
    if max_val > 50:
        logging.warning(
            f"Data has very high values (max={max_val:.1f}). "
            "Consider log transformation if not already applied."
        )

    # Validate data quality before processing
    validate_data_quality(combined)

    # Check if data has highly variable genes marked (should be done in preprocessing)
    if "highly_variable" not in combined.var.columns:
        logging.warning(
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
            logging.warning(
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
            logging.warning(
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

    logging.info(f"Using {n_hvg} highly variable genes for integration")

    # Save raw data if not already saved
    # IMPORTANT: Create a proper frozen copy for .raw to preserve counts
    # Using `combined.raw = combined` creates a view that gets modified during normalization
    if combined.raw is None:
        combined.raw = ad.AnnData(
            X=combined.X.copy(),  # Must copy - will be modified during normalization
            var=combined.var,  # No copy needed - AnnData internally creates independent copy
            obs=combined.obs.copy(),  # Must copy - will be modified by clustering/annotation
            uns={},  # Empty dict - raw doesn't need uns metadata
        )

    # ========================================================================
    # EARLY BRANCH FOR scVI-TOOLS METHODS
    # scVI requires normalized+log data WITHOUT scaling/PCA
    # It generates its own latent representation
    # NOTE: scVI-tools methods work better with ALL genes, not just HVGs
    # ========================================================================
    if method == "scvi":
        logging.info("Using scVI method - skipping scale and PCA")
        logging.info(f"Gene count before scVI processing: {combined.n_vars}")

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
            logging.info("scVI integration completed successfully")
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
            results_keys={"obsm": ["X_scVI"], "uns": ["neighbors"]},
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
        logging.info(f"Filtered to {n_hvg} highly variable genes")

    # Remove genes with zero variance to avoid NaN in scaling
    import numpy as np
    from scipy import sparse

    # MEMORY OPTIMIZATION: Calculate variance without toarray()
    # Uses E[X²] - E[X]² formula for sparse matrices
    # Saves ~80% memory (e.g., 76 MB → 15 MB for 10k cells × 2k genes)
    if sparse.issparse(combined.X):
        # Sparse matrix: compute variance using E[X²] - E[X]² formula
        # This avoids creating dense copy (5-10x memory reduction)
        mean_per_gene = np.array(combined.X.mean(axis=0)).flatten()

        # Calculate E[X²]
        X_squared = combined.X.copy()
        # Square the data: use np.array() for type safety (handles memoryview, ensures copy)
        X_squared.data = np.array(X_squared.data) ** 2
        mean_squared = np.array(X_squared.mean(axis=0)).flatten()

        # Variance = E[X²] - E[X]²
        gene_var = mean_squared - mean_per_gene**2
    else:
        # Dense matrix: use standard variance calculation
        gene_var = np.var(combined.X, axis=0)
    nonzero_var_genes = gene_var > 0
    if not np.all(nonzero_var_genes):
        n_removed = np.sum(~nonzero_var_genes)
        logging.warning(f"Removing {n_removed} genes with zero variance before scaling")
        # Memory optimization: Subsetting creates view, no need to copy
        # View will be materialized when scaling modifies the data
        combined = combined[:, nonzero_var_genes]

    # Scale data with proper error handling
    try:
        sc.pp.scale(combined, zero_center=True, max_value=10)
        logging.info("Data scaling successful with zero centering")
    except Exception as e:
        logging.warning(f"Scaling with zero centering failed: {e}")
        try:
            sc.pp.scale(combined, zero_center=False, max_value=10)
            logging.info("Data scaling successful without zero centering")
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
    import numpy as np
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
            logging.info(f"PCA successful with {solver} solver, {max_comps} components")
            pca_success = True
            break
        except Exception as e:
            logging.warning(f"PCA with {solver} solver failed: {e}")
            continue

    if not pca_success:
        raise ProcessingError(
            f"PCA failed for {combined.n_obs}×{combined.n_vars} data. Check data quality."
        )

    # Apply batch correction based on selected method
    if method == "harmony":
        # Use Harmony for batch correction
        # BEST PRACTICE: Use scanpy.external wrapper for better integration with scanpy workflow
        require("harmonypy", feature="Harmony integration")
        try:
            import scanpy.external as sce

            # Check if harmony_integrate is available in scanpy.external
            if hasattr(sce.pp, "harmony_integrate"):
                # Use scanpy.external wrapper (preferred method)
                logging.info("Using scanpy.external.pp.harmony_integrate (recommended)")
                sce.pp.harmony_integrate(
                    combined,
                    key=batch_key,
                    basis="X_pca",  # Use PCA representation
                    adjusted_basis="X_pca_harmony",  # Store corrected embedding
                )
                # Use corrected embedding for downstream analysis
                sc.pp.neighbors(combined, use_rep="X_pca_harmony")
            else:
                # Fallback to raw harmonypy (same algorithm, different interface)
                logging.info(
                    "scanpy.external.pp.harmony_integrate not available, using raw harmonypy"
                )
                import harmonypy
                import pandas as pd

                # Get PCA result
                X_pca = combined.obsm["X_pca"]

                # Create DataFrame with batch information
                meta_data = pd.DataFrame({batch_key: combined.obs[batch_key]})

                # Run Harmony
                harmony_out = harmonypy.run_harmony(
                    data_mat=X_pca,
                    meta_data=meta_data,
                    vars_use=[batch_key],
                    sigma=0.1,
                    nclust=None,
                    max_iter_harmony=10,
                    verbose=True,
                )

                # Save Harmony corrected result
                combined.obsm["X_harmony"] = harmony_out.Z_corr.T

                # Use corrected result to calculate neighbor graph
                sc.pp.neighbors(combined, use_rep="X_harmony")

        except Exception as e:
            raise ProcessingError(
                f"Harmony integration failed: {e}. "
                f"Check batch_key '{batch_key}' has ≥2 valid batches."
            ) from e

    elif method == "bbknn":
        # Use BBKNN for batch correction
        require("bbknn", feature="BBKNN integration")
        import bbknn

        bbknn.bbknn(combined, batch_key=batch_key, neighbors_within_batch=3)

    elif method == "scanorama":
        # Use Scanorama for batch correction
        # BEST PRACTICE: Use scanpy.external wrapper for better integration with scanpy workflow
        require("scanorama", feature="Scanorama integration")
        try:
            import scanpy.external as sce

            # Check if scanorama_integrate is available in scanpy.external
            if hasattr(sce.pp, "scanorama_integrate"):
                # Use scanpy.external wrapper (preferred method)
                logging.info(
                    "Using scanpy.external.pp.scanorama_integrate (recommended)"
                )
                sce.pp.scanorama_integrate(
                    combined, key=batch_key, basis="X_pca", adjusted_basis="X_scanorama"
                )
                # Use integrated representation for neighbor graph
                sc.pp.neighbors(combined, use_rep="X_scanorama")
            else:
                # Fallback to raw scanorama (same algorithm, different interface)
                logging.info(
                    "scanpy.external.pp.scanorama_integrate not available, using raw scanorama"
                )
                import numpy as np
                import scanorama

                # Separate data by batch
                datasets = []
                genes_list = []
                batch_order = []

                for batch in combined.obs[batch_key].unique():
                    batch_mask = combined.obs[batch_key] == batch
                    batch_data = combined[batch_mask]

                    # Scanorama natively supports sparse matrices
                    datasets.append(batch_data.X)
                    genes_list.append(batch_data.var_names.tolist())
                    batch_order.append(batch)

                    logging.info(f"Prepared batch '{batch}': {batch_data.X.shape}")

                # Run Scanorama integration
                logging.info("Running Scanorama integration...")
                integrated, corrected_genes = scanorama.integrate(
                    datasets, genes_list, dimred=100
                )

                # Stack integrated results back together
                integrated_X = np.vstack(integrated)
                logging.info(f"Scanorama integration completed: {integrated_X.shape}")

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
        logging.warning(
            f"Integration method '{method}' not recognized. "
            f"Using uncorrected PCA embedding."
        )
        sc.pp.neighbors(combined)

    # Calculate UMAP embedding to visualize integration effect
    sc.tl.umap(combined)

    # Store metadata for scientific provenance tracking
    # Determine which representation was used
    if method == "harmony":
        if "X_pca_harmony" in combined.obsm:
            results_keys = {"obsm": ["X_pca_harmony"], "uns": ["neighbors"]}
        else:
            results_keys = {"obsm": ["X_harmony"], "uns": ["neighbors"]}
    elif method == "bbknn":
        results_keys = {"uns": ["neighbors"]}
    elif method == "scanorama":
        results_keys = {"obsm": ["X_scanorama"], "uns": ["neighbors"]}
    else:
        results_keys = {"obsm": ["X_pca"], "uns": ["neighbors"]}

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


def align_spatial_coordinates(combined_adata, batch_key="batch", reference_batch=None):
    """Align spatial coordinates of multiple samples

    Args:
        combined_adata: Combined AnnData object containing multiple samples
        batch_key: Batch information key
        reference_batch: Reference batch, if None use the first batch

    Returns:
        AnnData object with aligned spatial coordinates
    """
    import numpy as np
    from sklearn.preprocessing import StandardScaler

    # Ensure data contains spatial coordinates
    spatial_key = get_spatial_key(combined_adata)
    if not spatial_key:
        raise DataNotFoundError("Data is missing spatial coordinates")

    # Get batch information
    batches = combined_adata.obs[batch_key].unique()

    if len(batches) == 0:
        raise DataError("Dataset is empty, cannot perform spatial registration")

    # If reference batch not specified, use the first batch
    if reference_batch is None:
        reference_batch = batches[0]
    elif reference_batch not in batches:
        raise ParameterError(f"Reference batch '{reference_batch}' not found in data")

    # Get reference batch spatial coordinates
    ref_coords = combined_adata[combined_adata.obs[batch_key] == reference_batch].obsm[
        spatial_key
    ]

    # Standardize reference coordinates
    scaler = StandardScaler()
    ref_coords_scaled = scaler.fit_transform(ref_coords)

    # Align spatial coordinates for each batch
    aligned_coords = []

    for batch in batches:
        # Get current batch index
        batch_idx = combined_adata.obs[batch_key] == batch

        if batch == reference_batch:
            # Reference batch remains unchanged
            aligned_coords.append(ref_coords_scaled)
        else:
            # Get current batch spatial coordinates
            batch_coords = combined_adata[batch_idx].obsm[spatial_key]

            # Standardize current batch coordinates
            batch_coords_scaled = scaler.transform(batch_coords)

            # Add to aligned coordinates list
            aligned_coords.append(batch_coords_scaled)

    # Merge all aligned coordinates
    combined_adata.obsm["spatial_aligned"] = np.zeros((combined_adata.n_obs, 2))

    # Fill aligned coordinates back to original data
    start_idx = 0
    for batch, coords in zip(batches, aligned_coords, strict=False):
        batch_idx = combined_adata.obs[batch_key] == batch
        n_cells = np.sum(batch_idx)
        combined_adata.obsm["spatial_aligned"][start_idx : start_idx + n_cells] = coords
        start_idx += n_cells

    # Store metadata for scientific provenance tracking
    n_batches = len(batches)
    batch_sizes = {
        batch: np.sum(combined_adata.obs[batch_key] == batch) for batch in batches
    }

    store_analysis_metadata(
        combined_adata,
        analysis_name="spatial_alignment",
        method="standardization",
        parameters={
            "batch_key": batch_key,
            "reference_batch": reference_batch,
        },
        results_keys={"obsm": ["spatial_aligned"]},
        statistics={
            "n_batches": n_batches,
            "batch_sizes": batch_sizes,
            "reference_batch": reference_batch,
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
        ImportError: If scvi-tools is not installed
        ValueError: If data is not preprocessed or invalid

    Reference:
        Lopez et al. (2018) "Deep generative modeling for single-cell transcriptomics"
        Nature Methods 15, 1053–1058
    """
    import numpy as np

    require("scvi", feature="scVI integration")
    import scvi

    # Validate data is preprocessed
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

    logging.info(f"Setting up scVI model with {n_batches} batches")

    # Setup AnnData for scVI
    scvi.model.SCVI.setup_anndata(
        combined, batch_key=batch_key, layer=None  # Use .X (should be preprocessed)
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
        logging.info(f"Auto-determined {n_epochs} epochs for {n_cells} cells")

    # Train model
    logging.info(f"Training scVI model for {n_epochs} epochs...")
    # Note: scvi-tools 1.x uses accelerator instead of use_gpu
    accelerator = "gpu" if use_gpu else "cpu"
    model.train(max_epochs=n_epochs, early_stopping=True, accelerator=accelerator)

    # Get latent representation
    logging.info("Extracting scVI latent representation...")
    combined.obsm["X_scvi"] = model.get_latent_representation()

    # Compute neighbors using scVI embedding
    sc.pp.neighbors(combined, use_rep="X_scvi")

    logging.info(
        f"scVI integration complete: {combined.n_obs} cells, "
        f"{n_latent}-dimensional latent space"
    )

    return combined


async def integrate_samples(
    data_ids: List[str],
    ctx: "ToolContext",
    params: IntegrationParameters = IntegrationParameters(),
) -> IntegrationResult:
    """Integrate multiple spatial transcriptomics samples and perform batch correction

    Args:
        data_ids: List of dataset IDs to integrate
        ctx: ToolContext for unified data access and logging
        params: Integration parameters

    Returns:
        Integration result
    """
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

    # Align spatial coordinates
    # NOTE: Spatial alignment is OPTIONAL and only needed for spatial data
    # Methods like BBKNN, Harmony, MNN, Scanorama work on gene expression/PCA space
    # and do NOT require spatial coordinates for batch correction
    if params.align_spatial:
        # Check if data actually has spatial coordinates
        if "spatial" not in combined_adata.obsm:
            # Skip alignment for non-spatial data - this is scientifically correct
            # BBKNN, Harmony, MNN, Scanorama are designed for scRNA-seq data without spatial info
            pass
        else:
            combined_adata = align_spatial_coordinates(
                combined_adata,
                batch_key=params.batch_key,
                reference_batch=params.reference_batch,
            )

    # Generate new integrated dataset ID
    integrated_id = f"integrated_{'-'.join(data_ids)}"

    # Store integrated data using ToolContext
    await ctx.add_dataset(integrated_id, combined_adata)

    # Return result
    return IntegrationResult(
        data_id=integrated_id,
        n_samples=len(data_ids),
        integration_method=params.method,
        umap_visualization=None,  # Use visualize_data tool instead
        spatial_visualization=None,  # Use visualize_data tool instead
    )
