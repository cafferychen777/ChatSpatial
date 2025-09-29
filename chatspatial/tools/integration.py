"""
Integration tools for spatial transcriptomics data.
"""

import logging
from typing import Any, Dict, List, Optional

import numpy as np
import scanpy as sc
from mcp.server.fastmcp import Context

from ..models.analysis import IntegrationResult
from ..models.data import IntegrationParameters


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
        raise ValueError(
            f"Dataset has only {adata.n_obs} cells, minimum {min_cells} required for integration. "
            f"Consider combining with other datasets or use single-sample analysis."
        )

    if adata.n_vars < min_genes:
        raise ValueError(
            f"Dataset has only {adata.n_vars} genes, minimum {min_genes} required for integration. "
            f"Check if data was properly loaded and genes were not over-filtered."
        )

    # Check for empty cells or genes
    if hasattr(adata.X, "toarray"):
        X_dense = adata.X.toarray()
    else:
        X_dense = adata.X

    cell_counts = np.sum(X_dense > 0, axis=1)
    gene_counts = np.sum(X_dense > 0, axis=0)

    empty_cells = np.sum(cell_counts == 0)
    empty_genes = np.sum(gene_counts == 0)

    if empty_cells > adata.n_obs * 0.1:
        raise ValueError(
            f"{empty_cells} cells ({empty_cells/adata.n_obs*100:.1f}%) have zero expression. "
            f"Check data quality and consider filtering."
        )

    if empty_genes > adata.n_vars * 0.5:
        raise ValueError(
            f"{empty_genes} genes ({empty_genes/adata.n_vars*100:.1f}%) have zero expression across all cells. "
            f"Consider gene filtering before integration."
        )

    logging.info(
        f"Data quality validation passed: {adata.n_obs} cells, {adata.n_vars} genes"
    )


def integrate_multiple_samples(adatas, batch_key="batch", method="harmony", n_pcs=30):
    """Integrate multiple spatial transcriptomics samples

    This function expects preprocessed data (normalized, log-transformed, with HVGs marked).
    Use preprocessing.py or preprocess_data() before calling this function.

    Args:
        adatas: List of preprocessed AnnData objects or a single combined AnnData object
        batch_key: Batch information key
        method: Integration method, options: 'harmony', 'bbknn', 'scanorama', 'mnn'
        n_pcs: Number of principal components for integration

    Returns:
        Integrated AnnData object with batch correction applied

    Raises:
        ValueError: If data is not properly preprocessed
    """

    # Merge datasets
    if isinstance(adatas, list):
        # Validate that each dataset has batch labels - NO AUTOMATIC CREATION
        for i, adata in enumerate(adatas):
            if batch_key not in adata.obs:
                raise ValueError(
                    f"Dataset {i} missing batch information in column '{batch_key}'.\n\n"
                    f"🔬 SAMPLE INTEGRATION REQUIRES BATCH LABELS\n"
                    f"Each dataset must have batch information for proper integration.\n\n"
                    f"🛠️ HOW TO ADD BATCH LABELS:\n"
                    f"• Same batch: adata.obs['{batch_key}'] = 'experiment_1'\n"
                    f"• Different batches: adata.obs['{batch_key}'] = 'batch_A', 'batch_B', etc.\n\n"
                    f"⚠️ IMPORTANT: Only use real batch information from your experiment.\n"
                    f"Don't create fake batch labels - this violates scientific integrity."
                )

        # Merge datasets
        combined = adatas[0].concatenate(
            adatas[1:],
            batch_key=batch_key,
            join="outer",  # Use outer join to keep all genes
        )
    else:
        # If already a merged dataset, ensure it has batch information
        combined = adatas
        if batch_key not in combined.obs:
            raise ValueError(
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
        raise ValueError(
            "Data appears to be raw counts (high positive values). "
            "Please normalize and log-transform data before integration. "
            "Use preprocessing.py or run: sc.pp.normalize_total(adata); sc.pp.log1p(adata)"
        )

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
    if combined.raw is None:
        combined.raw = combined

    # Filter to highly variable genes
    if "highly_variable" in combined.var.columns:
        n_hvg = combined.var["highly_variable"].sum()
        if n_hvg == 0:
            raise ValueError(
                "No highly variable genes found. Check HVG selection parameters."
            )
        combined = combined[:, combined.var["highly_variable"]].copy()
        logging.info(f"Filtered to {n_hvg} highly variable genes")

    # Remove genes with zero variance to avoid NaN in scaling
    import numpy as np

    if hasattr(combined.X, "todense"):
        X_check = np.asarray(combined.X.todense())
    else:
        X_check = np.asarray(combined.X)

    gene_var = np.var(X_check, axis=0)
    nonzero_var_genes = gene_var > 0
    if not np.all(nonzero_var_genes):
        n_removed = np.sum(~nonzero_var_genes)
        logging.warning(f"Removing {n_removed} genes with zero variance before scaling")
        combined = combined[:, nonzero_var_genes].copy()

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
            raise RuntimeError(
                f"Data scaling failed completely. Zero-center error: {e}. Non-zero-center error: {e2}. "
                f"This usually indicates data contains extreme outliers or invalid values. "
                f"Consider additional quality control or outlier removal."
            )

    # PCA with proper error handling
    # Determine safe number of components
    max_possible_components = min(n_pcs, combined.n_vars, combined.n_obs - 1)

    if max_possible_components < 2:
        raise ValueError(
            f"Cannot perform PCA: only {max_possible_components} components possible. "
            f"Dataset has {combined.n_obs} cells and {combined.n_vars} genes. "
            f"Minimum 2 components required for downstream analysis."
        )

    # Check data matrix before PCA
    import numpy as np

    if hasattr(combined.X, "todense"):
        X_check = np.asarray(combined.X.todense())
    else:
        X_check = np.asarray(combined.X)

    # Check for NaN or Inf
    if np.isnan(X_check).any():
        raise ValueError("Data contains NaN values after scaling")
    if np.isinf(X_check).any():
        raise ValueError("Data contains infinite values after scaling")

    # Check variance
    var_per_gene = np.var(X_check, axis=0)
    zero_var_genes = np.sum(var_per_gene == 0)
    if zero_var_genes > 0:
        logging.warning(
            f"Found {zero_var_genes} genes with zero variance after scaling"
        )

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
        raise RuntimeError(
            f"All PCA methods failed for dataset with {combined.n_obs} cells and {combined.n_vars} genes. "
            f"This usually indicates: \n"
            f"1. Data contains NaN or infinite values\n"
            f"2. All genes have identical expression\n"
            f"3. Data matrix is rank-deficient\n"
            f"4. Insufficient memory for computation\n"
            f"Please check data quality and preprocessing steps."
        )

    # Apply batch correction based on selected method
    if method == "harmony":
        # Use Harmony for batch correction
        try:
            import harmonypy

            # Get PCA result
            X_pca = combined.obsm["X_pca"]

            # Run Harmony - need to pass the DataFrame with batch info, not just the labels
            # Create a temporary DataFrame with batch information
            import pandas as pd

            meta_data = pd.DataFrame({batch_key: combined.obs[batch_key]})

            # Run Harmony with proper parameters
            harmony_out = harmonypy.run_harmony(
                data_mat=X_pca,  # PCA matrix
                meta_data=meta_data,  # DataFrame with batch info
                vars_use=[batch_key],  # Column name in meta_data
                sigma=0.1,  # Default parameter
                nclust=None,  # Let harmony determine the number of clusters
                max_iter_harmony=10,  # Default number of iterations
                verbose=True,  # Show progress
            )

            # Save Harmony corrected result
            combined.obsm["X_harmony"] = harmony_out.Z_corr.T

            # Use corrected result to calculate neighbor graph
            sc.pp.neighbors(combined, use_rep="X_harmony")
        except ImportError:
            raise ImportError(
                "harmonypy package is required for harmony integration. Install with 'pip install harmonypy'"
            )
        except Exception as e:
            # Provide clear error message instead of silent fallback
            logging.error(f"Harmony integration failed: {e}")

            # Check if it's an import error (harmonypy not installed)
            if "harmonypy" in str(e).lower():
                raise ImportError(
                    "Harmony integration failed due to missing harmonypy package. "
                    "Please install with: pip install harmonypy"
                )
            else:
                raise RuntimeError(
                    f"Harmony integration failed with error: {e}. "
                    f"This may be due to: \n"
                    f"1. Insufficient batch diversity (need at least 2 different batches)\n"
                    f"2. Batch key '{batch_key}' contains invalid values\n"
                    f"3. PCA result is corrupted or has NaN values\n"
                    f"Consider using method='mnn' or checking your batch labels."
                )

    elif method == "bbknn":
        # Use BBKNN for batch correction
        try:
            import bbknn

            bbknn.bbknn(combined, batch_key=batch_key, neighbors_within_batch=3)
        except ImportError:
            raise ImportError(
                "bbknn package is required for BBKNN integration. Install with 'pip install bbknn'"
            )

    elif method == "scanorama":
        # Use Scanorama for batch correction
        try:
            import numpy as np
            import scanorama

            # Separate data by batch
            datasets = []
            genes_list = []
            batch_order = []

            for batch in combined.obs[batch_key].unique():
                batch_mask = combined.obs[batch_key] == batch
                batch_data = combined[batch_mask]

                # Convert to dense array if sparse
                if hasattr(batch_data.X, "toarray"):
                    X_batch = batch_data.X.toarray()
                else:
                    X_batch = batch_data.X

                datasets.append(X_batch)
                genes_list.append(batch_data.var_names.tolist())
                batch_order.append(batch)

                logging.info(f"Prepared batch '{batch}': {X_batch.shape}")

            # Run Scanorama integration (returns low-dimensional embeddings)
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

        except ImportError:
            raise ImportError(
                "scanorama package is required for Scanorama integration. Install with 'pip install scanorama'"
            )
        except Exception as e:
            # Scanorama failed - do not fallback to inferior methods
            error_msg = (
                f"Scanorama integration failed: {e}. "
                f"Scanorama requires compatible data and proper preprocessing. "
                f"Please check: 1) Data quality and batch structure, 2) Gene overlap between batches, "
                f"or 3) Consider using method='harmony' or 'mnn' for more robust integration."
            )
            logging.error(error_msg)
            raise RuntimeError(error_msg)

    elif method == "mnn":
        # Use MNN for batch correction
        sc.pp.combat(combined, key=batch_key)
        sc.pp.neighbors(combined)

    else:
        # Default: use uncorrected PCA result
        sc.pp.neighbors(combined)

    # Calculate UMAP embedding to visualize integration effect
    sc.tl.umap(combined)

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
    if "spatial" not in combined_adata.obsm:
        raise ValueError("Data is missing spatial coordinates")

    # Get batch information
    batches = combined_adata.obs[batch_key].unique()

    # If reference batch not specified, use the first batch
    if reference_batch is None:
        reference_batch = batches[0]
    elif reference_batch not in batches:
        raise ValueError(f"Reference batch '{reference_batch}' not found in data")

    # Get reference batch spatial coordinates
    ref_coords = combined_adata[combined_adata.obs[batch_key] == reference_batch].obsm[
        "spatial"
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
            batch_coords = combined_adata[batch_idx].obsm["spatial"]

            # Standardize current batch coordinates
            batch_coords_scaled = scaler.transform(batch_coords)

            # Add to aligned coordinates list
            aligned_coords.append(batch_coords_scaled)

    # Merge all aligned coordinates
    combined_adata.obsm["spatial_aligned"] = np.zeros((combined_adata.n_obs, 2))

    # Fill aligned coordinates back to original data
    start_idx = 0
    for batch, coords in zip(batches, aligned_coords):
        batch_idx = combined_adata.obs[batch_key] == batch
        n_cells = np.sum(batch_idx)
        combined_adata.obsm["spatial_aligned"][start_idx : start_idx + n_cells] = coords
        start_idx += n_cells

    return combined_adata




# Import scvi-tools for advanced integration methods
try:
    import scvi
    from scvi.external import ContrastiveVI
except ImportError:
    scvi = None
    ContrastiveVI = None


async def integrate_with_contrastive_vi(
    adata,
    batch_key: str = "batch",
    condition_key: Optional[str] = None,
    n_epochs: int = 400,
    n_hidden: int = 128,
    n_background_latent: int = 10,
    n_salient_latent: int = 10,
    use_gpu: bool = False,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """Integrate data using ContrastiveVI for identifying condition-specific variations

    ContrastiveVI is particularly useful for identifying variations that are
    specific to certain conditions (e.g., disease vs. healthy) while accounting
    for batch effects.

    Args:
        adata: AnnData object with batch information
        batch_key: Key in adata.obs for batch information
        condition_key: Key in adata.obs for condition information (e.g., 'disease_state')
        n_epochs: Number of epochs for training
        n_hidden: Number of hidden units in neural networks
        n_background_latent: Dimensionality of background (shared) latent space
        n_salient_latent: Dimensionality of salient (condition-specific) latent space
        use_gpu: Whether to use GPU for training
        context: MCP context for logging

    Returns:
        Dictionary containing integration results

    Raises:
        ImportError: If scvi-tools package is not available
        ValueError: If required keys are missing
        RuntimeError: If integration fails
    """
    try:
        if scvi is None or ContrastiveVI is None:
            raise ImportError(
                "scvi-tools package with ContrastiveVI is required. Install with 'pip install scvi-tools'"
            )

        if context:
            await context.info("Starting ContrastiveVI integration...")

        # Validate batch key
        if batch_key not in adata.obs:
            raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")

        # If no condition key provided, use batch key as condition
        if condition_key is None:
            condition_key = batch_key
            if context:
                await context.info(
                    f"No condition key provided, using batch key '{batch_key}' as condition"
                )
        elif condition_key not in adata.obs:
            raise ValueError(f"Condition key '{condition_key}' not found in adata.obs")

        if context:
            await context.info(
                f"Integrating {adata.n_obs} cells with {len(adata.obs[batch_key].unique())} batches"
            )
            if condition_key != batch_key:
                await context.info(
                    f"Conditions: {len(adata.obs[condition_key].unique())} unique values"
                )

        # Setup ContrastiveVI
        ContrastiveVI.setup_anndata(
            adata, batch_key=batch_key, labels_key=condition_key
        )

        # Create ContrastiveVI model
        model = ContrastiveVI(
            adata,
            n_hidden=n_hidden,
            n_background_latent=n_background_latent,
            n_salient_latent=n_salient_latent,
        )

        if context:
            await context.info("Training ContrastiveVI model...")

        # ContrastiveVI requires background and target indices
        # Background indices: cells from a reference condition (e.g., healthy)
        # Target indices: cells from condition of interest (e.g., disease)

        # Get unique conditions
        conditions = adata.obs[condition_key].unique()
        if len(conditions) < 2:
            raise ValueError(
                f"ContrastiveVI requires at least 2 conditions, found {len(conditions)}"
            )

        # Use first condition as background, others as target
        background_condition = conditions[0]
        background_indices = np.where(adata.obs[condition_key] == background_condition)[
            0
        ]
        target_indices = np.where(adata.obs[condition_key] != background_condition)[0]

        if context:
            await context.info(
                f"Using '{background_condition}' as background ({len(background_indices)} cells)"
            )
            await context.info(
                f"Using other conditions as target ({len(target_indices)} cells)"
            )

        # Train model
        if use_gpu:
            model.train(
                background_indices=background_indices,
                target_indices=target_indices,
                max_epochs=n_epochs,
                accelerator="gpu",
            )
        else:
            model.train(
                background_indices=background_indices,
                target_indices=target_indices,
                max_epochs=n_epochs,
            )

        if context:
            await context.info("ContrastiveVI training completed")

        # Get results
        if context:
            await context.info("Extracting integrated representations...")

        # Get background (shared) latent representation
        background_latent = model.get_latent_representation(
            adata, representation_kind="background"
        )

        # Get salient (condition-specific) latent representation
        salient_latent = model.get_latent_representation(
            adata, representation_kind="salient"
        )

        # Store results in adata
        adata.obsm["X_contrastive_background"] = background_latent
        adata.obsm["X_contrastive_salient"] = salient_latent

        # Use background representation for standard analyses (UMAP, clustering)
        adata.obsm["X_integrated"] = background_latent

        # Calculate integration metrics
        # Compute silhouette score for batch mixing in background space
        from sklearn.metrics import silhouette_score

        try:
            background_silhouette = silhouette_score(
                background_latent, adata.obs[batch_key]
            )

            # For salient space, we expect separation by condition
            if condition_key != batch_key:
                salient_silhouette = silhouette_score(
                    salient_latent, adata.obs[condition_key]
                )
            else:
                salient_silhouette = background_silhouette
        except (ValueError, KeyError) as e:
            # Handle expected data-related errors (e.g., insufficient unique labels, missing keys)
            logger.warning(
                f"Could not calculate silhouette scores due to data issue: {e}"
            )
            background_silhouette = None  # Explicitly indicate unavailable metric
            salient_silhouette = None
        except Exception as e:
            # Handle unexpected errors but don't mask them
            error_msg = f"Integration quality metrics calculation failed: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        # Calculate summary statistics
        results = {
            "method": "ContrastiveVI",
            "n_background_latent": n_background_latent,
            "n_salient_latent": n_salient_latent,
            "n_epochs": n_epochs,
            "n_batches": len(adata.obs[batch_key].unique()),
            "n_conditions": len(adata.obs[condition_key].unique()),
            "background_mixing_score": (
                float(background_silhouette)
                if background_silhouette is not None
                else None
            ),
            "salient_separation_score": (
                float(salient_silhouette) if salient_silhouette is not None else None
            ),
            "background_latent_shape": background_latent.shape,
            "salient_latent_shape": salient_latent.shape,
            "integration_completed": True,
            "device": "GPU" if use_gpu else "CPU",
        }

        if context:
            await context.info("ContrastiveVI integration completed successfully")
            await context.info(
                "Background representation captures batch-corrected shared variation"
            )
            await context.info(
                "Salient representation captures condition-specific variation"
            )
            await context.info(
                "Stored in adata.obsm['X_contrastive_background'] and adata.obsm['X_contrastive_salient']"
            )

        return results

    except Exception as e:
        error_msg = f"ContrastiveVI integration failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg) from e


async def integrate_samples(
    data_ids: List[str],
    data_store: Dict[str, Any],
    params: IntegrationParameters = IntegrationParameters(),
    context: Optional[Context] = None,
) -> IntegrationResult:
    """Integrate multiple spatial transcriptomics samples and perform batch correction

    Args:
        data_ids: List of dataset IDs to integrate
        data_store: Dictionary storing datasets
        params: Integration parameters
        context: MCP context

    Returns:
        Integration result
    """
    if context:
        await context.info(
            f"Integrating {len(data_ids)} samples using {params.method} method"
        )

    # Collect all AnnData objects
    adatas = []
    for data_id in data_ids:
        if data_id not in data_store:
            raise ValueError(f"Dataset {data_id} not found in data store")
        adatas.append(data_store[data_id]["adata"].copy())

    # Integrate samples
    combined_adata = integrate_multiple_samples(
        adatas, batch_key=params.batch_key, method=params.method, n_pcs=params.n_pcs
    )

    if context:
        await context.info(
            f"Integration complete. Combined dataset has {combined_adata.n_obs} cells and {combined_adata.n_vars} genes"
        )

    # Align spatial coordinates
    if params.align_spatial:
        if context:
            await context.info("Aligning spatial coordinates")
        combined_adata = align_spatial_coordinates(
            combined_adata,
            batch_key=params.batch_key,
            reference_batch=params.reference_batch,
        )

    # Note: Visualizations should be created using the separate visualize_data tool
    # This maintains clean separation between analysis and visualization
    if context:
        await context.info(
            "Integration analysis complete. Use visualize_data tool with plot_type='integration_umap' or 'integration_spatial' to visualize results"
        )

    # Generate new integrated dataset ID
    integrated_id = f"integrated_{'-'.join(data_ids)}"

    # Store integrated data
    data_store[integrated_id] = {"adata": combined_adata}

    if context:
        await context.info(f"Integration complete. New dataset ID: {integrated_id}")

    # Return result
    return IntegrationResult(
        data_id=integrated_id,
        n_samples=len(data_ids),
        integration_method=params.method,
        umap_visualization=None,  # Use visualize_data tool instead
        spatial_visualization=None,  # Use visualize_data tool instead
    )
