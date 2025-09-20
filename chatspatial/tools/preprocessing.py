"""
Preprocessing tools for spatial transcriptomics data.
"""

import traceback
from typing import Any, Dict, Optional

import numpy as np
import scanpy as sc
import squidpy as sq
from mcp.server.fastmcp import Context

from ..models.analysis import PreprocessingResult
from ..models.data import AnalysisParameters
from ..utils.data_adapter import standardize_adata
from ..utils.tool_error_handling import mcp_tool_error_handler

# Import scvi-tools for advanced preprocessing
try:
    import scvi
    from scvi.external import RESOLVI
except ImportError:
    scvi = None
    RESOLVI = None

# Import scvelo for RNA velocity preprocessing
try:
    import scvelo as scv
except ImportError:
    scv = None

# Import SpaGCN for spatial domain-specific preprocessing
try:
    import SpaGCN as spg
except ImportError:
    spg = None

# Constants for preprocessing
DEFAULT_TARGET_SUM = 1e4
MAX_SCALE_VALUE = 10
MERFISH_GENE_THRESHOLD = 200
MIN_NEIGHBORS = 3
MAX_NEIGHBORS_RATIO = 0.1
MIN_KMEANS_CLUSTERS = 2
MAX_TSNE_PCA_COMPONENTS = 50

# Clustering resolution by dataset size
CLUSTERING_RESOLUTIONS = {
    "small": 0.4,  # < 100 cells
    "medium": 0.6,  # 100-500 cells
    "large": 0.8,  # > 500 cells
}


def _detect_data_type(adata) -> str:
    """Detect the type of spatial transcriptomics data"""
    if adata.n_vars < MERFISH_GENE_THRESHOLD:
        return "merfish"
    elif adata.n_vars > 10000:
        return "visium"
    else:
        return "other"


def _should_use_all_genes_for_hvg(adata) -> bool:
    """Check if we should use all genes for HVG selection (for small gene sets)"""
    # Only use all genes if we have a very small gene set (like MERFISH)
    return adata.n_vars < 100


def _safe_matrix_operation(adata, operation: str):
    """Safely perform matrix operations on sparse or dense matrices"""
    try:
        if hasattr(adata.X, "toarray"):
            # Sparse matrix - avoid converting to dense if possible
            if operation == "variance":
                # Use sparse-compatible variance calculation
                mean = np.array(adata.X.mean(axis=0)).flatten()
                var = np.array(adata.X.power(2).mean(axis=0)).flatten() - mean**2
                return var
            elif operation == "sum_axis1":
                return np.array(adata.X.sum(axis=1)).flatten()
            elif operation == "count_nonzero_axis1":
                return np.array((adata.X > 0).sum(axis=1)).flatten()
        else:
            # Dense matrix
            if operation == "variance":
                return np.var(adata.X, axis=0)
            elif operation == "sum_axis1":
                return np.sum(adata.X, axis=1)
            elif operation == "count_nonzero_axis1":
                return np.sum(adata.X > 0, axis=1)
    except Exception:
        # Matrix operation failed - returning None
        return None


@mcp_tool_error_handler()
async def preprocess_data(
    data_id: str,
    data_store: Dict[str, Any],
    params: AnalysisParameters = AnalysisParameters(),
    context: Optional[Context] = None,
) -> PreprocessingResult:
    """Preprocess spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Preprocessing parameters
        context: MCP context

    Returns:
        Preprocessing result summary
    """
    try:
        if context:
            await context.info(
                f"Preprocessing dataset {data_id} with {params.normalization} normalization"
            )

        # Retrieve the AnnData object from data store
        if data_id not in data_store:
            raise ValueError(f"Dataset {data_id} not found in data store")

        # Make a copy of the AnnData object to avoid modifying the original
        adata = data_store[data_id]["adata"].copy()

        # LINUS FIX: Standardize data format at the entry point
        # This eliminates all downstream special cases for data format handling
        if context:
            await context.info("Standardizing data structure to ChatSpatial format...")
        try:
            adata = standardize_adata(
                adata, copy=False, strict=False, preserve_original=True
            )
            if context:
                await context.info("✓ Data structure standardized successfully")
        except Exception as e:
            if context:
                await context.warning(
                    f"Data standardization failed: {e}. Proceeding with original data."
                )
            # Continue with original data if standardization fails

        # Validate input data
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError(
                f"Dataset {data_id} is empty: {adata.n_obs} cells, {adata.n_vars} genes"
            )

        # Detect data type for informational purposes
        data_type = _detect_data_type(adata)

        if context:
            await context.info(
                f"Detected data type: {data_type} ({adata.n_obs} cells, {adata.n_vars} genes)"
            )

        # 1. Calculate QC metrics
        if context:
            await context.info("Calculating QC metrics...")
        try:
            sc.pp.calculate_qc_metrics(adata, inplace=True)
        except Exception as e:
            if context:
                await context.warning(
                    f"Could not calculate QC metrics: {str(e)}. Computing from raw data."
                )

            # Create realistic QC metrics based on actual data using safe operations
            gene_counts = _safe_matrix_operation(adata, "sum_axis1")
            n_genes = _safe_matrix_operation(adata, "count_nonzero_axis1")

            # Ensure we have valid data
            if (
                gene_counts is None
                or n_genes is None
                or len(gene_counts) == 0
                or np.all(gene_counts == 0)
            ):
                gene_counts = (
                    np.ones(adata.n_obs) * DEFAULT_TARGET_SUM // 10
                )  # Realistic fallback
                n_genes = np.ones(adata.n_obs) * min(100, adata.n_vars)  # Fallback

            adata.obs["total_counts"] = gene_counts
            adata.obs["n_genes_by_counts"] = n_genes
            # Set mitochondrial percentage to 0 for spatial data (usually not available)
            adata.obs["pct_counts_mt"] = np.zeros(adata.n_obs)

        # Store original QC metrics before filtering
        qc_metrics = {
            "n_cells_before_filtering": int(adata.n_obs),
            "n_genes_before_filtering": int(adata.n_vars),
            "median_genes_per_cell": float(np.median(adata.obs.n_genes_by_counts)),
            "median_umi_per_cell": float(np.median(adata.obs.total_counts)),
        }

        # 2. Apply user-controlled data filtering and subsampling
        if context:
            await context.info("Applying data filtering and subsampling...")

        # Apply gene filtering using LLM-controlled parameters
        min_cells = params.filter_genes_min_cells
        if min_cells is not None and min_cells > 0:
            if context:
                await context.info(f"Filtering genes: min_cells={min_cells}")
            sc.pp.filter_genes(adata, min_cells=min_cells)

        # Apply cell filtering using LLM-controlled parameters
        min_genes = params.filter_cells_min_genes
        if min_genes is not None and min_genes > 0:
            if context:
                await context.info(f"Filtering cells: min_genes={min_genes}")
            sc.pp.filter_cells(adata, min_genes=min_genes)

        # Apply spot subsampling if requested
        if params.subsample_spots is not None and params.subsample_spots < adata.n_obs:
            if context:
                await context.info(
                    f"Subsampling from {adata.n_obs} to {params.subsample_spots} spots"
                )
            sc.pp.subsample(
                adata,
                n_obs=params.subsample_spots,
                random_state=params.subsample_random_seed,
            )

        # Apply gene subsampling if requested (after HVG selection)
        gene_subsample_requested = params.subsample_genes is not None

        # Save raw data before normalization (required for some analysis methods)
        if context:
            await context.info("Saving raw data for downstream analysis...")
        adata.raw = adata

        # Update QC metrics after filtering
        qc_metrics.update(
            {
                "n_cells_after_filtering": int(adata.n_obs),
                "n_genes_after_filtering": int(adata.n_vars),
            }
        )

        # 3. Normalize data
        if context:
            await context.info(
                f"Normalizing data using {params.normalization} method..."
            )

        if params.normalization == "log":
            # Standard log normalization
            sc.pp.normalize_total(adata, target_sum=DEFAULT_TARGET_SUM)
            sc.pp.log1p(adata)
        elif params.normalization == "sct":
            # SCTransform-like normalization
            if context:
                await context.info("Using simplified SCTransform normalization...")
            sc.pp.normalize_total(adata, target_sum=DEFAULT_TARGET_SUM)
            sc.pp.log1p(adata)
            # Note: Full SCTransform requires additional variance stabilization
            # For production use, consider integrating sctransform package
        elif params.normalization == "pearson_residuals":
            # Modern Pearson residuals normalization (recommended for UMI data)
            if context:
                await context.info("Using modern Pearson residuals normalization...")
            try:
                # Try experimental version first (more recent)
                if hasattr(sc.experimental.pp, "normalize_pearson_residuals"):
                    sc.experimental.pp.normalize_pearson_residuals(
                        adata, n_top_genes=min(params.n_hvgs, adata.n_vars)
                    )
                else:
                    # Fallback to regular version if experimental not available
                    sc.pp.normalize_total(adata, target_sum=DEFAULT_TARGET_SUM)
                    sc.pp.log1p(adata)
                    if context:
                        await context.warning(
                            "Pearson residuals not available, using log normalization fallback"
                        )
            except Exception as e:
                if context:
                    await context.warning(
                        f"Pearson residuals normalization failed: {e}. Using log normalization fallback."
                    )
                sc.pp.normalize_total(adata, target_sum=DEFAULT_TARGET_SUM)
                sc.pp.log1p(adata)

        # 4. Find highly variable genes and apply gene subsampling
        if context:
            await context.info("Finding highly variable genes...")

        # Determine number of HVGs to select
        if gene_subsample_requested:
            # User wants to subsample genes
            n_hvgs = min(params.subsample_genes, adata.n_vars - 1, params.n_hvgs)
            if context:
                await context.info(
                    f"User requested {params.subsample_genes} genes, selecting {n_hvgs} highly variable genes"
                )
        else:
            # Use standard HVG selection
            n_hvgs = min(params.n_hvgs, adata.n_vars - 1)
            if context:
                await context.info(f"Using {n_hvgs} highly variable genes...")

        # Check if we should use all genes (for very small gene sets)
        if _should_use_all_genes_for_hvg(adata):
            if context:
                await context.info(
                    f"Small gene set detected ({adata.n_vars} genes), using all genes for analysis"
                )
            adata.var["highly_variable"] = True
        else:
            try:
                sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs)
            except Exception as e:
                if context:
                    await context.warning(
                        f"HVG selection failed: {e}. Using all genes."
                    )
                adata.var["highly_variable"] = True

        # Apply gene subsampling if requested
        if gene_subsample_requested and params.subsample_genes < adata.n_vars:
            if "highly_variable" in adata.var and adata.var["highly_variable"].any():
                # Keep only highly variable genes
                adata = adata[:, adata.var["highly_variable"]].copy()
                if context:
                    await context.info(
                        f"Subsampled to {adata.n_vars} highly variable genes"
                    )
            else:
                # Fallback: keep top variable genes by variance using safe operation
                gene_var = _safe_matrix_operation(adata, "variance")
                if gene_var is not None:
                    top_genes_idx = np.argsort(gene_var)[-params.subsample_genes :]
                    adata = adata[:, top_genes_idx].copy()
                    if context:
                        await context.info(
                            f"Subsampled to {adata.n_vars} top variable genes"
                        )
                else:
                    if context:
                        await context.warning(
                            "Could not compute gene variance, keeping all genes"
                        )

        # 5. Batch effect correction (if applicable)
        if (
            params.batch_key in adata.obs
            and len(adata.obs[params.batch_key].unique()) > 1
        ):
            if context:
                await context.info(
                    "Detected batch information. Applying batch effect correction..."
                )
                # Warn about Combat limitations for large sparse matrices
                if (
                    hasattr(adata.X, "toarray") and adata.X.nnz > 1000000
                ):  # >1M non-zero elements
                    await context.warning(
                        "ComBat requires dense matrix conversion, may impact memory usage for large datasets. Consider using scVI or Harmony for better performance."
                    )
            try:
                sc.pp.combat(adata, key=params.batch_key)
                if context:
                    await context.info("Batch effect correction completed using ComBat")
            except Exception as e:
                if context:
                    await context.warning(
                        f"Batch effect correction failed: {e}. Continuing without correction. For complex batch effects, consider using integration tools like scVI or Harmony."
                    )

        # 6. Scale data (if requested)
        if params.scale:
            if context:
                await context.info("Scaling data...")
            try:
                # Trust scanpy's internal zero-variance handling and sparse matrix optimization
                sc.pp.scale(adata, max_value=MAX_SCALE_VALUE)

                # Clean up any NaN/Inf values that might remain (sparse-matrix safe)
                if hasattr(adata.X, "data"):
                    # Sparse matrix - only modify the data array
                    adata.X.data = np.nan_to_num(
                        adata.X.data,
                        nan=0.0,
                        posinf=MAX_SCALE_VALUE,
                        neginf=-MAX_SCALE_VALUE,
                    )
                else:
                    # Dense matrix
                    adata.X = np.nan_to_num(
                        adata.X,
                        nan=0.0,
                        posinf=MAX_SCALE_VALUE,
                        neginf=-MAX_SCALE_VALUE,
                    )

            except Exception as e:
                if context:
                    await context.warning(
                        f"Scaling failed: {e}. Continuing without scaling."
                    )

        # 7. Run PCA
        if context:
            await context.info("Running PCA...")

        # Adjust n_pcs based on dataset size
        n_pcs = min(
            params.n_pcs, adata.n_vars - 1, adata.n_obs - 1
        )  # Ensure we don't use more PCs than possible

        if context:
            await context.info(f"Using {n_pcs} principal components...")

        try:
            sc.tl.pca(adata, n_comps=n_pcs)
        except Exception as e:
            if context:
                await context.warning(f"PCA failed: {e}. Using reduced components.")
            # Fallback with fewer components
            n_pcs_fallback = min(10, adata.n_vars - 1, adata.n_obs - 1)
            try:
                sc.tl.pca(adata, n_comps=n_pcs_fallback)
                n_pcs = n_pcs_fallback  # Update n_pcs for downstream use
            except Exception as e2:
                # Both PCA methods failed - this is a critical failure
                raise RuntimeError(
                    f"All PCA methods failed. Original error: {str(e)}. Fallback error: {str(e2)}. "
                    "Cannot perform dimensionality reduction. Please check data quality or try different preprocessing parameters."
                )

        # 8. Compute neighbors graph
        if context:
            await context.info("Computing neighbors graph...")

        # Determine n_neighbors: user-specified or adaptive
        if params.n_neighbors is not None:
            n_neighbors = params.n_neighbors
            if context:
                await context.info(f"Using user-specified n_neighbors: {n_neighbors}")
        else:
            # Adaptive n_neighbors based on dataset size
            n_neighbors = min(
                10, int(adata.n_obs * MAX_NEIGHBORS_RATIO)
            )  # Use at most 10% of cells as neighbors
            n_neighbors = max(
                n_neighbors, MIN_NEIGHBORS
            )  # But at least MIN_NEIGHBORS neighbors
            if context:
                await context.info(
                    f"Auto-detected n_neighbors: {n_neighbors} (based on {adata.n_obs} cells)"
                )

        if context:
            await context.info(
                f"Using {n_neighbors} neighbors and {n_pcs} PCs for graph construction..."
            )

        try:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

            # 9. Run UMAP for visualization
            if context:
                await context.info("Running UMAP...")
            sc.tl.umap(adata)

            # 10. Run clustering
            if context:
                await context.info("Running Leiden clustering...")

            # Determine clustering resolution: user-specified or adaptive
            if params.clustering_resolution is not None:
                resolution = params.clustering_resolution
                if context:
                    await context.info(
                        f"Using user-specified clustering resolution: {resolution}"
                    )
            else:
                # Adaptive resolution based on dataset size
                if adata.n_obs < 100:
                    resolution = CLUSTERING_RESOLUTIONS["small"]
                    dataset_size = "small"
                elif adata.n_obs < 500:
                    resolution = CLUSTERING_RESOLUTIONS["medium"]
                    dataset_size = "medium"
                else:
                    resolution = CLUSTERING_RESOLUTIONS["large"]
                    dataset_size = "large"
                if context:
                    await context.info(
                        f"Auto-detected clustering resolution: {resolution} (for {dataset_size} dataset with {adata.n_obs} cells)"
                    )

            if context:
                await context.info(
                    f"Using Leiden clustering with resolution {resolution}..."
                )

            sc.tl.leiden(adata, resolution=resolution, key_added=params.cluster_key)

            # Count clusters
            n_clusters = len(adata.obs[params.cluster_key].unique())

            # Compute diffusion map for trajectory analysis
            if context:
                await context.info("Computing diffusion map for trajectory analysis...")
            sc.tl.diffmap(adata)
        except Exception as e:
            if context:
                await context.error(f"Clustering failed: {str(e)}")
            # Don't silently create fallback clustering - let the LLM handle the failure
            raise RuntimeError(
                f"Clustering failed: {str(e)}. "
                "Please check data quality or try different preprocessing parameters."
            )

        # 11. Add RNA velocity preprocessing if requested
        if getattr(params, "enable_rna_velocity", False):
            if context:
                await context.info("Adding RNA velocity preprocessing...")
            try:
                # Check for velocity data (spliced/unspliced layers)
                has_velocity_data = (
                    "spliced" in adata.layers and "unspliced" in adata.layers
                )

                if has_velocity_data:
                    # Use unified velocity computation from trajectory module
                    from ..models.data import RNAVelocityParameters
                    from .trajectory import compute_rna_velocity

                    # Create RNAVelocityParameters from AnalysisParameters
                    # Use embedded velocity_params if available, otherwise create from individual fields
                    if (
                        hasattr(params, "velocity_params")
                        and params.velocity_params is not None
                    ):
                        velocity_params = params.velocity_params
                    else:
                        # Build from individual fields for backward compatibility
                        velocity_params = RNAVelocityParameters(
                            mode=getattr(params, "velocity_mode", "stochastic"),
                            min_shared_counts=getattr(
                                params, "velocity_min_shared_counts", 30
                            ),
                            n_top_genes=getattr(params, "velocity_n_top_genes", 2000),
                            n_pcs=getattr(params, "velocity_n_pcs", 30),
                            n_neighbors=getattr(params, "velocity_n_neighbors", 30),
                        )

                    # Compute velocity with unified function
                    adata = compute_rna_velocity(adata, params=velocity_params)

                    if context:
                        await context.info(
                            f"RNA velocity preprocessing completed using {velocity_params.mode} mode"
                        )
                else:
                    if context:
                        await context.warning(
                            "RNA velocity requested but no spliced/unspliced layers found"
                        )
            except Exception as e:
                if context:
                    await context.warning(f"RNA velocity preprocessing failed: {e}")

        # 12. Add trajectory analysis preprocessing if requested
        if getattr(params, "enable_trajectory_analysis", False):
            if context:
                await context.info("Adding trajectory analysis preprocessing...")
            try:
                # Compute diffusion map for pseudotime analysis
                sc.tl.diffmap(adata)

                # Optionally compute DPT if root cell is specified
                if (
                    hasattr(params, "dpt_root_cell")
                    and params.dpt_root_cell is not None
                ):
                    if params.dpt_root_cell in adata.obs_names:
                        adata.uns["iroot"] = np.where(
                            adata.obs_names == params.dpt_root_cell
                        )[0][0]
                        sc.tl.dpt(adata)
                        if context:
                            await context.info(
                                "Diffusion pseudotime computed with specified root cell"
                            )
                    else:
                        if context:
                            await context.warning(
                                f"Root cell {params.dpt_root_cell} not found"
                            )

                if context:
                    await context.info("Trajectory analysis preprocessing completed")
            except Exception as e:
                if context:
                    await context.warning(
                        f"Trajectory analysis preprocessing failed: {e}"
                    )

        # 13. Add spatial domain-specific preprocessing if requested
        if getattr(params, "enable_spatial_domains", False) and spg is not None:
            if context:
                await context.info("Adding spatial domain-specific preprocessing...")
            try:
                # Apply SpaGCN-specific gene filtering
                spg.prefilter_genes(adata, min_cells=3)
                spg.prefilter_specialgenes(adata)
                if context:
                    await context.info("SpaGCN gene filtering completed")
            except Exception as e:
                if context:
                    await context.warning(f"Spatial domain preprocessing failed: {e}")

        # 14. For spatial data, compute spatial neighbors if not already done
        if params.spatial_key in adata.uns or any(
            params.spatial_key in key for key in adata.obsm.keys()
        ):
            if context:
                await context.info("Computing spatial neighbors...")
            try:
                # Check if spatial neighbors already exist
                if "spatial_connectivities" not in adata.obsp:
                    # For MERFISH data, we need to ensure the spatial coordinates are correctly formatted
                    if params.spatial_key in adata.obsm:
                        # Check if this is MERFISH data by looking at the shape and content
                        if adata.obsm[params.spatial_key].shape[1] == 2:
                            # This is likely MERFISH or other single-cell resolution spatial data
                            # Use delaunay method which works better for single-cell data
                            sq.gr.spatial_neighbors(
                                adata, coord_type="generic", delaunay=True
                            )
                        else:
                            # Use default parameters for other spatial data types
                            sq.gr.spatial_neighbors(adata)
                    else:
                        # Use default parameters if spatial key is in uns but not in obsm
                        sq.gr.spatial_neighbors(adata)
            except Exception as e:
                if context:
                    await context.warning(
                        f"Could not compute spatial neighbors: {str(e)}"
                    )
                    await context.info("Continuing without spatial neighbors...")

        # Store the processed AnnData object back in the data store
        data_store[data_id]["adata"] = adata

        # Return preprocessing result
        return PreprocessingResult(
            data_id=data_id,
            n_cells=adata.n_obs,
            n_genes=adata.n_vars,
            n_hvgs=(
                int(sum(adata.var.highly_variable))
                if "highly_variable" in adata.var
                else n_hvgs
            ),
            clusters=n_clusters,
            qc_metrics=qc_metrics,
        )

    except Exception as e:
        error_msg = f"Error in preprocessing: {str(e)}"
        tb = traceback.format_exc()
        if context:
            await context.warning(error_msg)
            await context.info(f"Error details: {tb}")
        raise RuntimeError(f"{error_msg}\n{tb}")


async def preprocess_with_resolvi(
    adata,
    n_epochs: int = 50,  # Reduced for testing
    n_hidden: int = 128,
    n_latent: int = 10,
    use_gpu: bool = False,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """Preprocess spatial transcriptomics data using RESOLVI

    RESOLVI addresses noise and bias in single-cell resolved spatial
    transcriptomics data through denoising and bias correction.

    Args:
        adata: Spatial transcriptomics AnnData object
        n_epochs: Number of epochs for training
        n_hidden: Number of hidden units in neural networks
        n_latent: Dimensionality of latent space
        use_gpu: Whether to use GPU for training
        context: MCP context for logging

    Returns:
        Dictionary containing RESOLVI preprocessing results

    Raises:
        ImportError: If scvi-tools package is not available
        ValueError: If input data is invalid
        RuntimeError: If RESOLVI computation fails
    """
    try:
        if scvi is None or RESOLVI is None:
            raise ImportError(
                "scvi-tools package is required for RESOLVI preprocessing. Install with 'pip install scvi-tools'"
            )

        if context:
            await context.info("Starting RESOLVI spatial data preprocessing...")

        # Validate spatial coordinates
        spatial_key = "spatial"  # Default spatial key for RESOLVI
        if spatial_key not in adata.obsm:
            if context:
                await context.warning(
                    "Spatial coordinates not found. RESOLVI works best with spatial information."
                )

        # Ensure data is in the correct format for RESOLVI
        import scipy.sparse as sp

        # Convert sparse matrices to dense if needed
        if sp.issparse(adata.X):
            adata.X = adata.X.toarray()

        # RESOLVI expects count data (integers)
        # Check if data appears to be normalized (has decimals)
        if np.any(adata.X != adata.X.astype(int)):
            if context:
                await context.warning(
                    "RESOLVI expects integer count data, but found decimals. Converting to integer counts."
                )
            # Round to nearest integer and ensure non-negative
            adata.X = np.maximum(0, np.round(adata.X)).astype(np.int32)
        else:
            # Ensure integer type
            adata.X = adata.X.astype(np.int32)

        if context:
            await context.info(
                f"Preprocessing {adata.n_obs} spots and {adata.n_vars} genes with RESOLVI"
            )

        # Setup RESOLVI
        if spatial_key in adata.obsm:
            # Ensure spatial coordinates are in the correct format for RESOLVI
            spatial_coords = adata.obsm[spatial_key]
            if spatial_coords.shape[1] == 2:
                # Add spatial coordinates as X_spatial for RESOLVI
                adata.obsm["X_spatial"] = spatial_coords
        else:
            # RESOLVI requires spatial coordinates, create dummy ones if not available
            if context:
                await context.warning(
                    "No spatial coordinates found. Creating dummy spatial coordinates for RESOLVI."
                )
            # Create a simple 2D grid layout
            n_spots = adata.n_obs
            grid_size = int(np.ceil(np.sqrt(n_spots)))
            x_coords = np.arange(n_spots) % grid_size
            y_coords = np.arange(n_spots) // grid_size
            adata.obsm["X_spatial"] = np.column_stack([x_coords, y_coords]).astype(
                np.float32
            )

        # RESOLVI setup without spatial_key parameter (not supported in current version)
        RESOLVI.setup_anndata(adata)

        # Create RESOLVI model with proper parameter types
        # Disable downsample_counts to avoid the LogNormal parameter issue
        resolvi_model = RESOLVI(
            adata,
            n_hidden=int(n_hidden),
            n_latent=int(n_latent),
            downsample_counts=False,  # Disable to avoid torch tensor type issues
        )

        if context:
            await context.info("Training RESOLVI model...")

        # Train model
        if use_gpu:
            resolvi_model.train(max_epochs=n_epochs, accelerator="gpu")
        else:
            resolvi_model.train(max_epochs=n_epochs)

        if context:
            await context.info("RESOLVI training completed")

        # Get results
        if context:
            await context.info("Extracting RESOLVI denoised data...")

        # Get denoised expression
        denoised_expression = resolvi_model.get_normalized_expression()

        # Get latent representation
        latent = resolvi_model.get_latent_representation()

        # Store results in adata
        adata.layers["resolvi_denoised"] = denoised_expression
        adata.obsm["X_resolvi_latent"] = latent

        # Calculate preprocessing statistics
        original_mean = np.mean(adata.X)
        denoised_mean = np.mean(denoised_expression)

        # Calculate noise reduction metrics
        if hasattr(adata.X, "toarray"):
            orig_data = adata.X.toarray()
        else:
            orig_data = adata.X

        if hasattr(denoised_expression, "toarray"):
            denoised_data = denoised_expression.toarray()
        else:
            denoised_data = denoised_expression

        noise_reduction = np.mean(np.abs(orig_data - denoised_data))

        # Calculate summary statistics
        results = {
            "method": "RESOLVI",
            "n_latent_dims": n_latent,
            "n_epochs": n_epochs,
            "denoising_completed": True,
            "original_mean_expression": float(original_mean),
            "denoised_mean_expression": float(denoised_mean),
            "noise_reduction_metric": float(noise_reduction),
            "latent_shape": latent.shape,
            "denoised_shape": denoised_expression.shape,
            "training_completed": True,
            "device": "GPU" if use_gpu else "CPU",
        }

        if context:
            await context.info("RESOLVI preprocessing completed successfully")
            await context.info(
                "Stored denoised data in adata.layers['resolvi_denoised']"
            )
            await context.info(
                "Stored latent representation in adata.obsm['X_resolvi_latent']"
            )
            await context.info(f"Noise reduction metric: {noise_reduction:.4f}")

        return results

    except Exception as e:
        error_msg = f"RESOLVI preprocessing failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg) from e
