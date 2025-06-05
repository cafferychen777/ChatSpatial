"""
Deconvolution tools for spatial transcriptomics data.

This module provides functions for deconvolving spatial transcriptomics data
to estimate cell type proportions in each spatial location.
"""

from typing import Dict, List, Optional, Any, Tuple, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from scipy.optimize import nnls
from sklearn.preprocessing import normalize
import traceback
import warnings
import sys
import os
from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image

# Try to import Spotiphy directly
# If it's not installed, we'll handle the ImportError when the function is called

from ..models.data import DeconvolutionParameters
from ..models.analysis import DeconvolutionResult
from ..utils.image_utils import fig_to_image, fig_to_base64, create_placeholder_image


def deconvolve_nnls(
    spatial_adata: ad.AnnData,
    reference_profiles: pd.DataFrame,
    n_top_genes: int = 2000,
    min_common_genes: int = 100
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using Non-negative Least Squares

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_profiles: Reference expression profiles (genes x cell_types)
        n_top_genes: Number of top genes to use
        min_common_genes: Minimum number of common genes required

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)

    Raises:
        ValueError: If input data is invalid or insufficient common genes
        RuntimeError: If NNLS computation fails
    """
    # Validate input data
    if spatial_adata is None:
        raise ValueError("Spatial AnnData object cannot be None")
    if reference_profiles is None:
        raise ValueError("Reference profiles cannot be None")
    if spatial_adata.n_obs == 0:
        raise ValueError("Spatial data contains no observations")
    if reference_profiles.shape[1] == 0:
        raise ValueError("Reference profiles contain no cell types")

    # Find common genes
    # Ensure we're working with unique gene names
    spatial_genes = list(set(spatial_adata.var_names))
    reference_genes = list(set(reference_profiles.index))

    common_genes = list(set(spatial_genes) & set(reference_genes))
    if len(common_genes) < min_common_genes:
        raise ValueError(
            f"Only {len(common_genes)} genes in common between spatial data and reference profiles. "
            f"Need at least {min_common_genes}. Consider using a different reference dataset or "
            f"reducing the min_common_genes parameter."
        )

    # Select top variable genes
    if len(common_genes) > n_top_genes:
        try:
            # Calculate gene variance
            X = spatial_adata[:, common_genes].X
            if isinstance(X, np.matrix) or isinstance(X, np.ndarray):
                gene_var = X.var(axis=0)
                if isinstance(gene_var, np.matrix):
                    gene_var = gene_var.A1
            else:
                # Handle sparse matrix
                gene_var = np.zeros(len(common_genes))
                X_dense = X.toarray() if hasattr(X, 'toarray') else np.array(X)
                gene_var = X_dense.var(axis=0)

            # Select top variable genes
            top_genes = [common_genes[i] for i in np.argsort(-gene_var)[:n_top_genes]]
        except Exception as e:
            warnings.warn(f"Error selecting top variable genes: {str(e)}. Using all common genes instead.")
            top_genes = common_genes
    else:
        top_genes = common_genes

    # Extract expression matrices
    try:
        X_spatial = spatial_adata[:, top_genes].X
        if isinstance(X_spatial, np.matrix):
            X_spatial = X_spatial.A
        elif hasattr(X_spatial, 'toarray'):
            X_spatial = X_spatial.toarray()
    except Exception as e:
        raise RuntimeError(f"Failed to extract spatial expression matrix: {str(e)}")

    # Extract reference profiles for the same genes
    try:
        ref_profiles = reference_profiles.loc[top_genes]

        # Check for NaN values
        if ref_profiles.isna().any().any():
            warnings.warn("Reference profiles contain NaN values. Filling with zeros.")
            ref_profiles = ref_profiles.fillna(0)

        # Normalize reference profiles
        ref_profiles = normalize(ref_profiles, axis=0, norm='l1')
    except Exception as e:
        raise RuntimeError(f"Failed to process reference profiles: {str(e)}")

    # Run NNLS for each spot
    proportions = np.zeros((spatial_adata.n_obs, ref_profiles.shape[1]))
    residuals = np.zeros(spatial_adata.n_obs)
    failed_spots = 0

    for i in range(spatial_adata.n_obs):
        try:
            spot_expression = X_spatial[i]
            if isinstance(spot_expression, np.ndarray) and spot_expression.ndim == 0:
                spot_expression = np.array([float(spot_expression)])
            elif hasattr(spot_expression, 'toarray'):
                spot_expression = spot_expression.toarray().flatten()

            # Check for NaN values
            if np.isnan(spot_expression).any():
                warnings.warn(f"Spot {i} contains NaN values. Filling with zeros.")
                spot_expression = np.nan_to_num(spot_expression)

            # Solve NNLS
            props, res = nnls(ref_profiles, spot_expression)
            proportions[i] = props
            residuals[i] = res
        except Exception as e:
            warnings.warn(f"Failed to deconvolve spot {i}: {str(e)}. Setting to zero.")
            proportions[i] = np.zeros(ref_profiles.shape[1])
            residuals[i] = np.nan
            failed_spots += 1

    # Check if too many spots failed
    if failed_spots > 0:
        warnings.warn(f"Failed to deconvolve {failed_spots} out of {spatial_adata.n_obs} spots.")
    if failed_spots == spatial_adata.n_obs:
        raise RuntimeError("Failed to deconvolve any spots. Check your input data.")

    # Create DataFrame with results
    prop_df = pd.DataFrame(
        proportions,
        index=spatial_adata.obs_names,
        columns=reference_profiles.columns
    )

    # Normalize proportions to sum to 1
    row_sums = prop_df.sum(axis=1)
    if (row_sums == 0).any():
        warnings.warn(f"Some spots have zero total proportion. These will remain as zeros.")
    prop_df = prop_df.div(row_sums, axis=0).fillna(0)

    # Calculate statistics
    stats = {
        "mean_residual": np.nanmean(residuals),
        "median_residual": np.nanmedian(residuals),
        "genes_used": len(top_genes),
        "common_genes": len(common_genes),
        "cell_types": list(reference_profiles.columns),
        "n_cell_types": len(reference_profiles.columns),
        "mean_proportions": prop_df.mean().to_dict(),
        "failed_spots": failed_spots,
        "success_rate": 1.0 - (failed_spots / spatial_adata.n_obs)
    }

    return prop_df, stats





def deconvolve_cell2location(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str = 'cell_type',
    n_epochs: int = 10000,
    n_cells_per_spot: int = 10,
    use_gpu: bool = False,
    min_common_genes: int = 100
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using Cell2location

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_epochs: Number of epochs for training
        n_cells_per_spot: Expected number of cells per spot
        use_gpu: Whether to use GPU for training
        min_common_genes: Minimum number of common genes required

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)

    Raises:
        ImportError: If cell2location package is not installed
        ValueError: If input data is invalid or insufficient common genes
        RuntimeError: If cell2location computation fails
    """
    # Validate input data
    if spatial_adata is None:
        raise ValueError("Spatial AnnData object cannot be None")
    if reference_adata is None:
        raise ValueError("Reference AnnData object cannot be None")
    if spatial_adata.n_obs == 0:
        raise ValueError("Spatial data contains no observations")
    if reference_adata.n_obs == 0:
        raise ValueError("Reference data contains no observations")
    if cell_type_key not in reference_adata.obs:
        raise ValueError(f"Cell type key '{cell_type_key}' not found in reference data")
    if len(reference_adata.obs[cell_type_key].unique()) < 2:
        raise ValueError(f"Reference data must contain at least 2 cell types, found {len(reference_adata.obs[cell_type_key].unique())}")

    # Import cell2location
    try:
        import cell2location
        from cell2location.models import RegressionModel, Cell2location
    except ImportError:
        raise ImportError(
            "Cell2location package is required but not installed. "
            "Install with 'pip install cell2location' or "
            "'pip install spatial-transcriptomics-mcp[deconvolution]'"
        )

    try:
        # Prepare reference data
        ref = reference_adata.copy()

        # Normalize reference data if not already normalized
        if 'normalized' not in ref.uns:
            sc.pp.normalize_total(ref)
            sc.pp.log1p(ref)
            ref.uns['normalized'] = True

        # Prepare spatial data
        sp = spatial_adata.copy()

        # Normalize spatial data if not already normalized
        if 'normalized' not in sp.uns:
            sc.pp.normalize_total(sp)
            sc.pp.log1p(sp)
            sp.uns['normalized'] = True

        # Find common genes
        common_genes = list(set(sp.var_names) & set(ref.var_names))
        if len(common_genes) < min_common_genes:
            raise ValueError(
                f"Only {len(common_genes)} genes in common between spatial data and reference data. "
                f"Need at least {min_common_genes}. Consider using a different reference dataset or "
                f"reducing the min_common_genes parameter."
            )

        # Subset to common genes
        ref = ref[:, common_genes]
        sp = sp[:, common_genes]

        # Log progress
        print(f"Training cell2location model with {len(common_genes)} common genes and {len(ref.obs[cell_type_key].unique())} cell types")
        print(f"Reference data shape: {ref.shape}, Spatial data shape: {sp.shape}")

        # Check if cell type key has valid values
        if ref.obs[cell_type_key].isna().any():
            warnings.warn(f"Reference data contains NaN values in {cell_type_key}. These cells will be excluded.")
            ref = ref[~ref.obs[cell_type_key].isna()].copy()

        # Train regression model to get reference cell type signatures
        try:
            RegressionModel.setup_anndata(
                adata=ref,
                labels_key=cell_type_key
            )
        except Exception as e:
            raise RuntimeError(f"Failed to setup reference data for RegressionModel: {str(e)}")

        try:
            mod = RegressionModel(ref)
            # Remove use_gpu parameter which is no longer supported in newer versions
            mod.train(max_epochs=n_epochs)
        except Exception as e:
            error_msg = str(e)
            tb = traceback.format_exc()
            raise RuntimeError(f"Failed to train RegressionModel: {error_msg}\n{tb}")

        # Export reference signatures
        try:
            # Remove use_gpu parameter from sample_kwargs
            ref_signatures = mod.export_posterior(
                adata=ref,
                sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
            )
        except Exception as e:
            raise RuntimeError(f"Failed to export reference signatures: {str(e)}")

        # Prepare spatial data for cell2location model
        try:
            Cell2location.setup_anndata(
                adata=sp,
                batch_key=None
            )
        except Exception as e:
            raise RuntimeError(f"Failed to setup spatial data for Cell2location: {str(e)}")

        # Run cell2location model
        try:
            mod = Cell2location(
                sp,
                cell_state_df=ref_signatures,
                N_cells_per_location=n_cells_per_spot,
                detection_alpha=20.0
            )
            # Remove use_gpu parameter which is no longer supported in newer versions
            mod.train(max_epochs=n_epochs)
        except Exception as e:
            error_msg = str(e)
            tb = traceback.format_exc()
            raise RuntimeError(f"Failed to train Cell2location model: {error_msg}\n{tb}")

        # Export results
        try:
            # Remove use_gpu parameter from sample_kwargs
            sp = mod.export_posterior(
                sp,
                sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
            )
        except Exception as e:
            raise RuntimeError(f"Failed to export Cell2location results: {str(e)}")

        # Get cell abundance
        if 'q05_cell_abundance_w_sf' not in sp.obsm:
            raise RuntimeError("Cell2location did not produce expected output 'q05_cell_abundance_w_sf'")

        cell_abundance = sp.obsm['q05_cell_abundance_w_sf']

        # Create DataFrame with results
        proportions = pd.DataFrame(
            cell_abundance,
            index=sp.obs_names,
            columns=ref_signatures.index
        )

        # Check for NaN values
        if proportions.isna().any().any():
            warnings.warn("Cell2location returned NaN values. Filling with zeros.")
            proportions = proportions.fillna(0)

        # Check for negative values
        if (proportions < 0).any().any():
            warnings.warn("Cell2location returned negative values. Setting to zero.")
            proportions[proportions < 0] = 0

        # Normalize proportions to sum to 1
        row_sums = proportions.sum(axis=1)
        if (row_sums == 0).any():
            warnings.warn(f"Some spots have zero total proportion. These will remain as zeros.")
        proportions = proportions.div(row_sums, axis=0).fillna(0)

        # Calculate statistics
        stats = {
            "cell_types": list(proportions.columns),
            "n_cell_types": len(proportions.columns),
            "mean_proportions": proportions.mean().to_dict(),
            "genes_used": len(common_genes),
            "common_genes": len(common_genes),
            "n_epochs": n_epochs,
            "n_cells_per_spot": n_cells_per_spot,
            "method": "Cell2location",
            "use_gpu": use_gpu
        }

        # Add model performance metrics if available
        if hasattr(mod, 'history') and mod.history is not None:
            try:
                history = mod.history
                if 'elbo_train' in history and len(history['elbo_train']) > 0:
                    stats["final_elbo"] = float(history['elbo_train'][-1])
                if 'elbo_validation' in history and len(history['elbo_validation']) > 0:
                    stats["final_elbo_validation"] = float(history['elbo_validation'][-1])
            except Exception as e:
                warnings.warn(f"Failed to extract model history: {str(e)}")

        return proportions, stats

    except Exception as e:
        if not isinstance(e, (ValueError, ImportError, RuntimeError)):
            error_msg = str(e)
            tb = traceback.format_exc()
            print(f"Cell2location deconvolution failed: {error_msg}")
            print("Falling back to NNLS method...")

            # Fall back to NNLS method
            return _fallback_to_nnls(spatial_adata, reference_adata, cell_type_key)
        else:
            # For specific errors, also try to fall back to NNLS
            print(f"Cell2location failed with error: {str(e)}")
            print("Falling back to NNLS method...")

            # Fall back to NNLS method
            return _fallback_to_nnls(spatial_adata, reference_adata, cell_type_key)


def _fallback_to_nnls(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Fallback to NNLS method when Cell2location fails

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    print("Using NNLS method as fallback for deconvolution")

    # Get cell type information
    cell_types = sorted(list(reference_adata.obs[cell_type_key].unique()))

    # Create reference profiles
    ref_profiles = {}
    for cell_type in cell_types:
        try:
            # Get cells of this type
            cells = reference_adata[reference_adata.obs[cell_type_key] == cell_type]
            if cells.n_obs == 0:
                print(f"No cells found for cell type '{cell_type}'. Skipping.")
                continue

            # Calculate mean expression
            mean_expr = np.mean(cells.X, axis=0)
            if isinstance(mean_expr, np.matrix):
                mean_expr = mean_expr.A1
            elif hasattr(mean_expr, 'toarray'):
                mean_expr = mean_expr.toarray().flatten()

            ref_profiles[cell_type] = mean_expr
        except Exception as e:
            print(f"Failed to process cell type '{cell_type}': {str(e)}. Skipping.")

    if not ref_profiles:
        raise ValueError("Failed to create any reference profiles from reference data")

    # Create reference profiles DataFrame
    reference_profiles = pd.DataFrame(ref_profiles, index=reference_adata.var_names)

    # Run NNLS deconvolution
    try:
        proportions, stats = deconvolve_nnls(
            spatial_adata,
            reference_profiles
        )
        print(f"NNLS deconvolution completed successfully. Found {stats['n_cell_types']} cell types.")

        # Update method in stats
        stats["method"] = "NNLS (fallback from Cell2location)"

        return proportions, stats
    except Exception as e:
        error_msg = str(e)
        tb = traceback.format_exc()
        print(f"NNLS deconvolution failed: {error_msg}")
        raise RuntimeError(f"Both Cell2location and fallback NNLS failed: {error_msg}\n{tb}")


# Check if Spotiphy is available
def is_spotiphy_available() -> Tuple[bool, str]:
    """Check if Spotiphy and its dependencies are available

    Returns:
        Tuple of (is_available, error_message)
    """
    try:
        # Try to import Spotiphy
        import spotiphy

        # Check for required modules
        try:
            # First check if we can import the specific modules we need
            from spotiphy import sc_reference, deconvolution
            import torch

            # Try to import stardist but don't fail if it's not available
            try:
                import stardist
            except ImportError:
                pass  # We'll handle this later if needed

            # Try to import tensorflow but don't fail if it's not available
            try:
                import tensorflow
            except ImportError:
                pass  # We'll handle this later if needed

            return True, ""
        except ImportError as e:
            if "stardist" in str(e):
                return False, "Spotiphy requires stardist which is not installed. Install with 'pip install stardist'"
            elif "tensorflow" in str(e) or "keras" in str(e):
                return False, "Spotiphy requires tensorflow which is not installed. Install with 'pip install tensorflow'"
            elif "torch" in str(e):
                return False, "Spotiphy requires PyTorch which is not installed. Install with 'pip install torch'"
            else:
                return False, f"Spotiphy dependency missing: {str(e)}"
    except ImportError:
        return False, "Spotiphy package is not installed. Install with 'pip install spotiphy'"


def deconvolve_spotiphy(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str = 'cell_type',
    n_epochs: int = 8000,
    batch_prior: float = 2.0,
    adam_params: Optional[Dict[str, Any]] = None,
    use_gpu: bool = False,
    min_common_genes: int = 100
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using Spotiphy

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_epochs: Number of epochs for training
        batch_prior: Parameter of the prior distribution of the batch effect
        adam_params: Parameters for the Adam optimizer
        use_gpu: Whether to use GPU for training
        min_common_genes: Minimum number of common genes required

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)

    Raises:
        ImportError: If Spotiphy package is not installed
        ValueError: If input data is invalid or insufficient common genes
        RuntimeError: If Spotiphy computation fails
    """
    # Validate input data
    if spatial_adata is None:
        raise ValueError("Spatial AnnData object cannot be None")
    if reference_adata is None:
        raise ValueError("Reference AnnData object cannot be None")
    if spatial_adata.n_obs == 0:
        raise ValueError("Spatial data contains no observations")
    if reference_adata.n_obs == 0:
        raise ValueError("Reference data contains no observations")
    if cell_type_key not in reference_adata.obs:
        raise ValueError(f"Cell type key '{cell_type_key}' not found in reference data")
    if len(reference_adata.obs[cell_type_key].unique()) < 2:
        raise ValueError(f"Reference data must contain at least 2 cell types, found {len(reference_adata.obs[cell_type_key].unique())}")

    # Check if Spotiphy is available
    is_available, error_message = is_spotiphy_available()
    if not is_available:
        raise ImportError(
            f"{error_message}. "
            "For complete installation, run: "
            "'pip install spotiphy stardist tensorflow torch'"
        )

    # Import required modules
    from spotiphy import sc_reference, deconvolution
    import torch

    # Set device
    device = 'cuda' if use_gpu and torch.cuda.is_available() else 'cpu'
    if use_gpu and not torch.cuda.is_available():
        warnings.warn("CUDA is not available. Using CPU instead.")

    try:
        # Prepare data
        # Find common genes
        common_genes = list(set(spatial_adata.var_names) & set(reference_adata.var_names))
        if len(common_genes) < min_common_genes:
            raise ValueError(
                f"Only {len(common_genes)} genes in common between spatial data and reference data. "
                f"Need at least {min_common_genes}. Consider using a different reference dataset or "
                f"reducing the min_common_genes parameter."
            )

        # Subset to common genes
        spatial_data = spatial_adata[:, common_genes].copy()
        reference_data = reference_adata[:, common_genes].copy()

        # Normalize data
        sc_ref_data, spatial_data = sc_reference.initialization(
            reference_data,
            spatial_data,
            filtering=False,
            verbose=1
        )

        # Get cell types
        type_list = sorted(list(sc_ref_data.obs[cell_type_key].unique()))

        # Construct reference profiles
        sc_ref = sc_reference.construct_sc_ref(sc_ref_data, cell_type_key)

        # Extract expression matrices
        X = spatial_data.X
        if isinstance(X, np.matrix):
            X = X.A
        elif hasattr(X, 'toarray'):
            X = X.toarray()

        # Set Adam parameters
        if adam_params is None:
            adam_params = {"lr": 0.003, "betas": (0.95, 0.999)}

        # Run Spotiphy deconvolution
        print(f"Running Spotiphy deconvolution with {len(common_genes)} common genes and {len(type_list)} cell types")
        print(f"Device: {device}, Epochs: {n_epochs}")

        # Estimate cell proportions
        cell_proportions = deconvolution.estimation_proportion(
            X,
            sc_ref_data,
            sc_ref,
            type_list,
            cell_type_key,
            device=device,
            n_epoch=n_epochs,
            adam_params=adam_params,
            batch_prior=batch_prior,
            plot=False
        )

        # Create DataFrame with results
        proportions = pd.DataFrame(
            cell_proportions,
            index=spatial_data.obs_names,
            columns=type_list
        )

        # Calculate statistics
        stats = {
            "cell_types": type_list,
            "n_cell_types": len(type_list),
            "mean_proportions": proportions.mean().to_dict(),
            "genes_used": len(common_genes),
            "common_genes": len(common_genes),
            "n_epochs": n_epochs,
            "batch_prior": batch_prior,
            "method": "Spotiphy",
            "device": device
        }

        return proportions, stats

    except Exception as e:
        if not isinstance(e, (ValueError, ImportError)):
            error_msg = str(e)
            tb = traceback.format_exc()
            raise RuntimeError(f"Spotiphy deconvolution failed: {error_msg}\n{tb}")
        else:
            raise


def visualize_deconvolution_results(
    spatial_adata: ad.AnnData,
    proportions: pd.DataFrame,
    n_cell_types: int = 6,
    cmap: str = 'viridis',
    size: int = 50,
    title_prefix: str = ""
) -> plt.Figure:
    """Visualize deconvolution results

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        proportions: Cell type proportions DataFrame
        n_cell_types: Number of top cell types to visualize
        cmap: Colormap to use for visualization
        size: Size of points in spatial plots
        title_prefix: Prefix to add to plot titles

    Returns:
        Matplotlib figure with visualization

    Raises:
        ValueError: If input data is invalid
        RuntimeError: If visualization fails
    """
    # Validate input data
    if spatial_adata is None:
        raise ValueError("Spatial AnnData object cannot be None")
    if proportions is None or proportions.empty:
        raise ValueError("Proportions DataFrame cannot be None or empty")
    if n_cell_types < 1:
        raise ValueError("Number of cell types to visualize must be at least 1")

    try:
        # Get top cell types by mean proportion
        n_cell_types = min(n_cell_types, proportions.shape[1])
        top_cell_types = proportions.mean().sort_values(ascending=False).index[:n_cell_types]

        # Create figure
        n_cols = min(3, n_cell_types)
        n_rows = (n_cell_types + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
        axes = axes.flatten() if n_cell_types > 1 else [axes]

        # Set figure title if provided
        if title_prefix:
            fig.suptitle(f"{title_prefix} Cell Type Proportions", fontsize=16)
            plt.subplots_adjust(top=0.9)

        # Plot each cell type
        for i, cell_type in enumerate(top_cell_types):
            if i < len(axes):
                ax = axes[i]
                try:
                    # Create temporary AnnData with cell type proportions
                    temp_adata = spatial_adata.copy()
                    temp_adata.obs['proportion'] = proportions[cell_type]

                    # Check for NaN values
                    if temp_adata.obs['proportion'].isna().any():
                        warnings.warn(f"Cell type {cell_type} contains NaN values. Filling with zeros.")
                        temp_adata.obs['proportion'] = temp_adata.obs['proportion'].fillna(0)

                    # Plot spatial distribution
                    if 'spatial' in temp_adata.obsm:
                        sc.pl.embedding(
                            temp_adata,
                            basis='spatial',
                            color='proportion',
                            title=cell_type,
                            color_map=cmap,
                            size=size,
                            ax=ax,
                            show=False
                        )
                    else:
                        # Fallback to UMAP if available
                        if 'X_umap' in temp_adata.obsm:
                            sc.pl.umap(
                                temp_adata,
                                color='proportion',
                                title=cell_type,
                                color_map=cmap,
                                size=size,
                                ax=ax,
                                show=False
                            )
                        else:
                            # Just show a bar plot
                            sorted_props = proportions[cell_type].sort_values(ascending=False)
                            ax.bar(range(len(sorted_props)), sorted_props.values)
                            ax.set_title(cell_type)
                            ax.set_xlabel('Spots (sorted)')
                            ax.set_ylabel('Proportion')

                    # Add mean proportion to title
                    mean_prop = proportions[cell_type].mean()
                    ax.set_title(f"{cell_type}\n(Mean: {mean_prop:.3f})")

                except Exception as e:
                    warnings.warn(f"Failed to visualize cell type {cell_type}: {str(e)}")
                    ax.text(0.5, 0.5, f"Failed to visualize {cell_type}:\n{str(e)}",
                            ha='center', va='center', transform=ax.transAxes, wrap=True)
                    ax.set_xticks([])
                    ax.set_yticks([])

        # Hide empty axes
        for i in range(n_cell_types, len(axes)):
            axes[i].axis('off')

        plt.tight_layout()
        return fig

    except Exception as e:
        error_msg = str(e)
        tb = traceback.format_exc()
        raise RuntimeError(f"Failed to visualize deconvolution results: {error_msg}\n{tb}")


async def deconvolve_spatial_data(
    data_id: str,
    data_store: Dict[str, Any],
    params: DeconvolutionParameters = DeconvolutionParameters(),
    context: Optional[Context] = None
) -> DeconvolutionResult:
    """Deconvolve spatial transcriptomics data to estimate cell type proportions

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing datasets
        params: Deconvolution parameters
        context: MCP context

    Returns:
        Deconvolution result with cell type proportions and visualization

    Raises:
        ValueError: If input data is invalid or required parameters are missing
        RuntimeError: If deconvolution computation fails
    """
    try:
        # Validate input parameters
        if not data_id:
            raise ValueError("Dataset ID cannot be empty")
        if not data_store:
            raise ValueError("Data store cannot be empty")

        if context:
            await context.info(f"Deconvolving spatial data using {params.method} method")
            await context.info(f"Parameters: {params.model_dump()}")

        # Get spatial data
        if data_id not in data_store:
            raise ValueError(f"Dataset {data_id} not found in data store. Available datasets: {list(data_store.keys())}")
        if "adata" not in data_store[data_id]:
            raise ValueError(f"Dataset {data_id} does not contain AnnData object")

        spatial_adata = data_store[data_id]["adata"]
        if spatial_adata.n_obs == 0:
            raise ValueError(f"Dataset {data_id} contains no observations")

        # Ensure spatial data has unique gene names
        if hasattr(spatial_adata, 'var_names_make_unique'):
            if context:
                await context.info("Ensuring spatial dataset has unique gene names")
            spatial_adata.var_names_make_unique()

        if context:
            await context.info(f"Spatial dataset shape: {spatial_adata.shape}")

        # Get reference data if provided
        reference_adata = None
        if params.reference_data_id:
            if params.reference_data_id not in data_store:
                raise ValueError(
                    f"Reference dataset {params.reference_data_id} not found in data store. "
                    f"Available datasets: {list(data_store.keys())}"
                )
            if "adata" not in data_store[params.reference_data_id]:
                raise ValueError(f"Reference dataset {params.reference_data_id} does not contain AnnData object")

            reference_adata = data_store[params.reference_data_id]["adata"]
            if reference_adata.n_obs == 0:
                raise ValueError(f"Reference dataset {params.reference_data_id} contains no observations")

            # Ensure reference data has unique gene names
            if hasattr(reference_adata, 'var_names_make_unique'):
                if context:
                    await context.info("Ensuring reference dataset has unique gene names")
                reference_adata.var_names_make_unique()

            if context:
                await context.info(f"Reference dataset shape: {reference_adata.shape}")
                if params.cell_type_key in reference_adata.obs:
                    cell_types = reference_adata.obs[params.cell_type_key].unique()
                    await context.info(f"Reference dataset contains {len(cell_types)} cell types: {list(cell_types)}")

        # Run deconvolution based on selected method
        if params.method == "nnls":
            if context:
                await context.info("Running NNLS deconvolution")

            # Check if reference profiles are provided
            if params.reference_profiles is None:
                if reference_adata is None:
                    raise ValueError(
                        "Reference data or profiles required for NNLS deconvolution. "
                        "Please provide either reference_data_id or reference_profiles."
                    )

                # Create reference profiles from reference data
                if context:
                    await context.info("Creating reference profiles from reference data")

                # Get cell type information
                if params.cell_type_key not in reference_adata.obs:
                    raise ValueError(
                        f"Cell type key '{params.cell_type_key}' not found in reference data. "
                        f"Available keys: {list(reference_adata.obs.columns)}"
                    )

                # Aggregate by cell type
                cell_types = reference_adata.obs[params.cell_type_key].unique()
                ref_profiles = {}

                for cell_type in cell_types:
                    try:
                        # Get cells of this type
                        cells = reference_adata[reference_adata.obs[params.cell_type_key] == cell_type]
                        if cells.n_obs == 0:
                            if context:
                                await context.warning(f"No cells found for cell type '{cell_type}'. Skipping.")
                            continue

                        # Calculate mean expression
                        mean_expr = cells.X.mean(axis=0)
                        if isinstance(mean_expr, np.matrix):
                            mean_expr = mean_expr.A1
                        ref_profiles[cell_type] = mean_expr
                    except Exception as e:
                        if context:
                            await context.warning(f"Failed to process cell type '{cell_type}': {str(e)}. Skipping.")

                if not ref_profiles:
                    raise ValueError("Failed to create any reference profiles from reference data")

                # Create reference profiles DataFrame
                reference_profiles = pd.DataFrame(ref_profiles, index=reference_adata.var_names)

                # Ensure reference profiles have unique gene names
                if hasattr(reference_profiles.index, 'duplicated') and reference_profiles.index.duplicated().any():
                    if context:
                        await context.info("Ensuring reference profiles have unique gene names")
                    reference_profiles = reference_profiles.loc[~reference_profiles.index.duplicated(keep='first')]

                if context:
                    await context.info(f"Created reference profiles for {len(ref_profiles)} cell types")
            else:
                # Use provided reference profiles
                try:
                    reference_profiles = pd.DataFrame(params.reference_profiles)
                    if context:
                        await context.info(f"Using provided reference profiles for {reference_profiles.shape[1]} cell types")
                except Exception as e:
                    raise ValueError(f"Failed to convert provided reference profiles to DataFrame: {str(e)}")

            # Run NNLS deconvolution
            try:
                proportions, stats = deconvolve_nnls(
                    spatial_adata,
                    reference_profiles,
                    n_top_genes=params.n_top_genes
                )
                if context:
                    await context.info(f"NNLS deconvolution completed successfully. Found {stats['n_cell_types']} cell types.")
            except Exception as e:
                error_msg = str(e)
                tb = traceback.format_exc()
                if context:
                    await context.warning(f"NNLS deconvolution failed: {error_msg}")
                raise RuntimeError(f"NNLS deconvolution failed: {error_msg}\n{tb}")



        elif params.method == "cell2location":
            if context:
                await context.info("Running Cell2location deconvolution")

            # Check if reference data is provided
            if reference_adata is None:
                raise ValueError(
                    "Reference data required for Cell2location deconvolution. "
                    "Please provide reference_data_id."
                )

            # Check cell type key
            if params.cell_type_key not in reference_adata.obs:
                raise ValueError(
                    f"Cell type key '{params.cell_type_key}' not found in reference data. "
                    f"Available keys: {list(reference_adata.obs.columns)}"
                )

            # Run Cell2location
            try:
                proportions, stats = deconvolve_cell2location(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.n_epochs,
                    n_cells_per_spot=params.n_cells_per_spot or 10,
                    use_gpu=params.use_gpu
                )
                if context:
                    await context.info(f"Cell2location completed successfully. Found {stats['n_cell_types']} cell types.")
            except Exception as e:
                error_msg = str(e)
                tb = traceback.format_exc()
                if context:
                    await context.warning(f"Cell2location failed: {error_msg}")
                    await context.info("Falling back to NNLS method for deconvolution")

                # Fall back to NNLS method
                if context:
                    await context.info("Using NNLS method as fallback")

                # Create reference profiles from reference data
                if context:
                    await context.info("Creating reference profiles from reference data")

                # Get cell type information
                cell_types = sorted(list(reference_adata.obs[params.cell_type_key].unique()))

                # Create reference profiles
                ref_profiles = {}
                for cell_type in cell_types:
                    try:
                        # Get cells of this type
                        cells = reference_adata[reference_adata.obs[params.cell_type_key] == cell_type]
                        if cells.n_obs == 0:
                            if context:
                                await context.warning(f"No cells found for cell type '{cell_type}'. Skipping.")
                            continue

                        # Calculate mean expression
                        mean_expr = np.mean(cells.X, axis=0)
                        if isinstance(mean_expr, np.matrix):
                            mean_expr = mean_expr.A1
                        elif hasattr(mean_expr, 'toarray'):
                            mean_expr = mean_expr.toarray().flatten()

                        ref_profiles[cell_type] = mean_expr
                    except Exception as e:
                        if context:
                            await context.warning(f"Failed to process cell type '{cell_type}': {str(e)}. Skipping.")

                if not ref_profiles:
                    raise ValueError("Failed to create any reference profiles from reference data")

                # Create reference profiles DataFrame
                reference_profiles = pd.DataFrame(ref_profiles, index=reference_adata.var_names)

                # Run NNLS deconvolution
                try:
                    proportions, stats = deconvolve_nnls(
                        spatial_adata,
                        reference_profiles,
                        n_top_genes=params.n_top_genes
                    )
                    if context:
                        await context.info(f"NNLS deconvolution completed successfully. Found {stats['n_cell_types']} cell types.")

                    # Update method in stats
                    stats["method"] = "NNLS (fallback from Cell2location)"
                except Exception as e:
                    error_msg = str(e)
                    tb = traceback.format_exc()
                    if context:
                        await context.warning(f"NNLS deconvolution failed: {error_msg}")
                    raise RuntimeError(f"Both Cell2location and fallback NNLS failed: {error_msg}\n{tb}")

        elif params.method == "spotiphy":
            if context:
                await context.info("Running Spotiphy deconvolution")

            # Check if reference data is provided
            if reference_adata is None:
                raise ValueError(
                    "Reference data required for Spotiphy deconvolution. "
                    "Please provide reference_data_id."
                )

            # Check cell type key
            if params.cell_type_key not in reference_adata.obs:
                raise ValueError(
                    f"Cell type key '{params.cell_type_key}' not found in reference data. "
                    f"Available keys: {list(reference_adata.obs.columns)}"
                )

            # Check if Spotiphy is available
            is_available, error_message = is_spotiphy_available()
            if not is_available:
                if context:
                    await context.warning(f"Spotiphy is not available: {error_message}")
                    await context.info("Falling back to NNLS method for deconvolution")

                # Fall back to NNLS method
                if context:
                    await context.info("Using NNLS method as fallback")

                # Create reference profiles from reference data
                if context:
                    await context.info("Creating reference profiles from reference data")

                # Get cell type information
                cell_types = sorted(list(reference_adata.obs[params.cell_type_key].unique()))

                # Create reference profiles
                ref_profiles = {}
                for cell_type in cell_types:
                    try:
                        # Get cells of this type
                        cells = reference_adata[reference_adata.obs[params.cell_type_key] == cell_type]
                        if cells.n_obs == 0:
                            if context:
                                await context.warning(f"No cells found for cell type '{cell_type}'. Skipping.")
                            continue

                        # Calculate mean expression
                        mean_expr = np.mean(cells.X, axis=0)
                        if isinstance(mean_expr, np.matrix):
                            mean_expr = mean_expr.A1
                        elif hasattr(mean_expr, 'toarray'):
                            mean_expr = mean_expr.toarray().flatten()

                        ref_profiles[cell_type] = mean_expr
                    except Exception as e:
                        if context:
                            await context.warning(f"Failed to process cell type '{cell_type}': {str(e)}. Skipping.")

                if not ref_profiles:
                    raise ValueError("Failed to create any reference profiles from reference data")

                # Create reference profiles DataFrame
                reference_profiles = pd.DataFrame(ref_profiles, index=reference_adata.var_names)

                # Run NNLS deconvolution
                try:
                    proportions, stats = deconvolve_nnls(
                        spatial_adata,
                        reference_profiles,
                        n_top_genes=params.n_top_genes
                    )
                    if context:
                        await context.info(f"NNLS deconvolution completed successfully. Found {stats['n_cell_types']} cell types.")
                except Exception as e:
                    error_msg = str(e)
                    tb = traceback.format_exc()
                    if context:
                        await context.warning(f"NNLS deconvolution failed: {error_msg}")
                    raise RuntimeError(f"NNLS deconvolution failed: {error_msg}\n{tb}")
            else:
                # Set Adam parameters
                adam_params = {
                    "lr": params.spotiphy_adam_lr,
                    "betas": params.spotiphy_adam_betas
                }

                # Run Spotiphy
                try:
                    proportions, stats = deconvolve_spotiphy(
                        spatial_adata,
                        reference_adata,
                        cell_type_key=params.cell_type_key,
                        n_epochs=params.n_epochs,
                        batch_prior=params.spotiphy_batch_prior,
                        adam_params=adam_params,
                        use_gpu=params.use_gpu
                    )
                    if context:
                        await context.info(f"Spotiphy completed successfully. Found {stats['n_cell_types']} cell types.")
                except Exception as e:
                    error_msg = str(e)
                    tb = traceback.format_exc()
                    if context:
                        await context.warning(f"Spotiphy failed: {error_msg}")
                        await context.info("Falling back to NNLS method for deconvolution")

                    # Fall back to NNLS method
                    if context:
                        await context.info("Using NNLS method as fallback")

                    # Create reference profiles from reference data
                    if context:
                        await context.info("Creating reference profiles from reference data")

                    # Get cell type information
                    cell_types = sorted(list(reference_adata.obs[params.cell_type_key].unique()))

                    # Create reference profiles
                    ref_profiles = {}
                    for cell_type in cell_types:
                        try:
                            # Get cells of this type
                            cells = reference_adata[reference_adata.obs[params.cell_type_key] == cell_type]
                            if cells.n_obs == 0:
                                if context:
                                    await context.warning(f"No cells found for cell type '{cell_type}'. Skipping.")
                                continue

                            # Calculate mean expression
                            mean_expr = np.mean(cells.X, axis=0)
                            if isinstance(mean_expr, np.matrix):
                                mean_expr = mean_expr.A1
                            elif hasattr(mean_expr, 'toarray'):
                                mean_expr = mean_expr.toarray().flatten()

                            ref_profiles[cell_type] = mean_expr
                        except Exception as e:
                            if context:
                                await context.warning(f"Failed to process cell type '{cell_type}': {str(e)}. Skipping.")

                    if not ref_profiles:
                        raise ValueError("Failed to create any reference profiles from reference data")

                    # Create reference profiles DataFrame
                    reference_profiles = pd.DataFrame(ref_profiles, index=reference_adata.var_names)

                    # Run NNLS deconvolution
                    try:
                        proportions, stats = deconvolve_nnls(
                            spatial_adata,
                            reference_profiles,
                            n_top_genes=params.n_top_genes
                        )
                        if context:
                            await context.info(f"NNLS deconvolution completed successfully. Found {stats['n_cell_types']} cell types.")
                    except Exception as e:
                        error_msg = str(e)
                        tb = traceback.format_exc()
                        if context:
                            await context.warning(f"NNLS deconvolution failed: {error_msg}")
                        raise RuntimeError(f"NNLS deconvolution failed: {error_msg}\n{tb}")

        else:
            raise ValueError(
                f"Unsupported deconvolution method: {params.method}. "
                f"Supported methods are: nnls, cell2location, spotiphy"
            )

        # Store proportions in AnnData object
        proportions_key = f"deconvolution_{params.method}"
        spatial_adata.obsm[proportions_key] = proportions.values

        # Store cell type names in uns
        cell_types_key = f"{proportions_key}_cell_types"
        spatial_adata.uns[cell_types_key] = list(proportions.columns)

        # Also store individual cell type proportions in obs for easier visualization
        for cell_type in proportions.columns:
            obs_key = f"{proportions_key}_{cell_type}"
            spatial_adata.obs[obs_key] = proportions[cell_type].values

        if context:
            await context.info(f"Stored cell type proportions in adata.obsm['{proportions_key}']")
            await context.info(f"Stored cell type names in adata.uns['{cell_types_key}']")
            await context.info(f"Stored individual cell type proportions in adata.obs with prefix '{proportions_key}_'")

        # Create visualization
        if context:
            await context.info("Creating visualization")

        try:
            fig = visualize_deconvolution_results(
                spatial_adata,
                proportions,
                n_cell_types=min(6, proportions.shape[1]),
                title_prefix=f"{params.method.upper()}"
            )

            # Convert figure to image
            img = fig_to_image(fig)
            if context:
                await context.info("Visualization created successfully")
        except Exception as e:
            if context:
                await context.warning(f"Failed to create visualization: {str(e)}. Using placeholder image.")
            img = create_placeholder_image(f"Failed to create visualization: {str(e)}")

        # Add cell type annotation based on deconvolution results
        if context:
            await context.info("Adding cell type annotation based on deconvolution results")

        # Determine the dominant cell type for each spot
        dominant_cell_types = []
        for i in range(proportions.shape[0]):
            row = proportions.iloc[i]
            max_idx = row.argmax()
            dominant_cell_types.append(proportions.columns[max_idx])

        # Add to adata.obs
        spatial_adata.obs['cell_type'] = dominant_cell_types

        # Make it categorical
        spatial_adata.obs['cell_type'] = spatial_adata.obs['cell_type'].astype('category')

        if context:
            await context.info(f"Added cell type annotation with {len(proportions.columns)} cell types")
            await context.info(f"Most common cell type: {spatial_adata.obs['cell_type'].value_counts().index[0]}")

        # Return result
        result = DeconvolutionResult(
            data_id=data_id,
            method=params.method,
            cell_types=list(proportions.columns),
            n_cell_types=proportions.shape[1],
            proportions_key=proportions_key,
            visualization=img,
            statistics=stats,
            visualization_params={
                "feature": f"{proportions_key}:{proportions.columns[0]}",
                "show_deconvolution": True,
                "n_cell_types": 4,
                "plot_type": "spatial",
                "colormap": "viridis"
            }
        )

        if context:
            await context.info(f"Deconvolution completed successfully with {result.n_cell_types} cell types")

        return result

    except Exception as e:
        if not isinstance(e, (ValueError, ImportError, RuntimeError)):
            error_msg = str(e)
            tb = traceback.format_exc()
            if context:
                await context.warning(f"Deconvolution failed with unexpected error: {error_msg}")
            raise RuntimeError(f"Deconvolution failed with unexpected error: {error_msg}\n{tb}")
        else:
            if context:
                await context.warning(f"Deconvolution failed: {str(e)}")
            raise
