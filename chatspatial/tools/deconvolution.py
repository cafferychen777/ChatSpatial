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

# Import context manager utilities
from contextlib import contextmanager
import logging
from io import StringIO
from contextlib import redirect_stdout, redirect_stderr


# Helper functions to eliminate redundancy


def _validate_deconvolution_inputs(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str,
    min_common_genes: int = 100
) -> List[str]:
    """Validate inputs and return a list of common genes.
    
    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        min_common_genes: Minimum number of common genes required
        
    Returns:
        List of common genes
        
    Raises:
        ValueError: If input data is invalid or insufficient common genes
    """
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

    # Find common genes
    common_genes = list(set(spatial_adata.var_names) & set(reference_adata.var_names))
    if len(common_genes) < min_common_genes:
        raise ValueError(
            f"Only {len(common_genes)} genes in common between spatial data and reference data. "
            f"Need at least {min_common_genes}. Consider using a different reference dataset or "
            f"reducing the min_common_genes parameter."
        )
    
    return common_genes


def _get_device(use_gpu: bool, method: str) -> str:
    """Determine the appropriate compute device.
    
    Args:
        use_gpu: Whether to use GPU for training
        method: Name of the method (for logging purposes)
        
    Returns:
        Device string ('cuda', 'mps', or 'cpu')
    """
    try:
        import torch
    except ImportError:
        warnings.warn("PyTorch not available. Using CPU.")
        return "cpu"

    if not use_gpu:
        print(f"Using CPU for {method} (GPU acceleration disabled)")
        return "cpu"
    
    if torch.cuda.is_available():
        print(f"Using CUDA GPU acceleration for {method}")
        return "cuda"
    
    # MPS support - handle differently for different methods
    if hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
        if method.lower() == 'spotiphy':
            print(f"Using Apple MPS acceleration for {method}")
            return "mps"
        else:
            # Cell2location has issues with MPS
            warnings.warn(
                f"MPS acceleration is available but disabled for {method} due to numerical instability issues "
                f"with cell2location 0.1.4. Using CPU instead. Consider upgrading cell2location "
                f"or PyTorch for better MPS support in the future."
            )
            print(f"Using CPU for {method} (MPS disabled due to compatibility issues)")
            return "cpu"
    
    warnings.warn(f"GPU requested for {method} but neither CUDA nor MPS is available. Using CPU instead.")
    print(f"Using CPU for {method}")
    return "cpu"


def _prepare_anndata_for_counts(adata: ad.AnnData, data_name: str) -> ad.AnnData:
    """Ensure AnnData object has raw integer counts in .X
    
    Args:
        adata: AnnData object to prepare
        data_name: Name of the data for logging purposes
        
    Returns:
        AnnData object with integer counts
    """
    # Check if already in correct format
    if adata.X.dtype == np.int32 or adata.X.dtype == np.int64:
        return adata
    
    warnings.warn(f"{data_name} data is not in integer format. Attempting conversion to counts.")
    
    # Convert sparse matrix to dense if needed for processing
    if hasattr(adata.X, 'toarray'):
        X_dense = adata.X.toarray()
    else:
        X_dense = adata.X.copy()

    # If data appears to be log-transformed, try to reverse it
    if X_dense.max() < 20:  # Likely log-transformed
        warnings.warn(f"{data_name} data appears to be log-transformed. Converting back to counts.")
        X_dense = np.expm1(X_dense)
    
    # Round to integers and ensure non-negative
    adata.X = np.round(np.maximum(X_dense, 0)).astype(np.int32)
    return adata


@contextmanager
def suppress_output():
    """Context manager to temporarily suppress stdout, stderr, and logging."""
    old_level = logging.getLogger().level
    logging.getLogger().setLevel(logging.ERROR)
    
    with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
        try:
            yield
        finally:
            logging.getLogger().setLevel(old_level)


def _create_deconvolution_stats(
    proportions: pd.DataFrame,
    common_genes: List[str],
    method: str,
    device: str,
    **method_specific_params
) -> Dict[str, Any]:
    """Create standardized statistics dictionary for deconvolution results.
    
    Args:
        proportions: Cell type proportions DataFrame
        common_genes: List of common genes used
        method: Deconvolution method name
        device: Device used for computation
        **method_specific_params: Additional method-specific parameters
        
    Returns:
        Statistics dictionary
    """
    stats = {
        "cell_types": list(proportions.columns),
        "n_cell_types": len(proportions.columns),
        "mean_proportions": proportions.mean().to_dict(),
        "genes_used": len(common_genes),
        "common_genes": len(common_genes),
        "method": method,
        "device": device
    }
    
    # Add method-specific parameters
    stats.update(method_specific_params)
    
    return stats


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
    # Unified validation and gene finding
    common_genes = _validate_deconvolution_inputs(spatial_adata, reference_adata, cell_type_key, min_common_genes)

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
        # Unified device selection
        device = _get_device(use_gpu, 'Cell2location')

        # Prepare data using helper functions - Cell2location expects raw count data
        ref = _prepare_anndata_for_counts(reference_adata.copy(), "Reference")
        sp = _prepare_anndata_for_counts(spatial_adata.copy(), "Spatial")

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
            # Create RegressionModel
            mod = RegressionModel(ref)

            # Use the suppress_output context manager
            with suppress_output():
                # Train with appropriate accelerator
                if device == "cuda":
                    mod.train(max_epochs=n_epochs, accelerator='gpu')
                else:
                    mod.train(max_epochs=n_epochs)

        except Exception as e:
            error_msg = str(e)
            tb = traceback.format_exc()
            raise RuntimeError(f"Failed to train RegressionModel: {error_msg}\n{tb}")

        # Export reference signatures
        try:
            # Export the estimated cell abundance (summary of the posterior distribution)
            ref = mod.export_posterior(
                ref,
                sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
            )

            # Extract estimated expression in each cluster
            if "means_per_cluster_mu_fg" in ref.varm.keys():
                ref_signatures = ref.varm["means_per_cluster_mu_fg"][
                    [f"means_per_cluster_mu_fg_{i}" for i in ref.uns["mod"]["factor_names"]]
                ].copy()
            else:
                ref_signatures = ref.var[[f"means_per_cluster_mu_fg_{i}" for i in ref.uns["mod"]["factor_names"]]].copy()
            ref_signatures.columns = ref.uns["mod"]["factor_names"]
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
            # Create Cell2location model (device specification is handled in train method)
            mod = Cell2location(
                sp,
                cell_state_df=ref_signatures,
                N_cells_per_location=n_cells_per_spot,
                detection_alpha=20.0
            )

            # Use the suppress_output context manager
            with suppress_output():
                # Train with appropriate accelerator
                if device == "cuda":
                    mod.train(max_epochs=n_epochs, batch_size=2500, accelerator='gpu')
                else:
                    mod.train(max_epochs=n_epochs, batch_size=2500)

        except Exception as e:
            error_msg = str(e)
            tb = traceback.format_exc()
            raise RuntimeError(f"Failed to train Cell2location model: {error_msg}\n{tb}")

        # Export results
        try:
            # Export the estimated cell abundance (summary of the posterior distribution)
            sp = mod.export_posterior(
                sp,
                sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
            )
        except Exception as e:
            raise RuntimeError(f"Failed to export Cell2location results: {str(e)}")

        # Get cell abundance - try different possible keys
        cell_abundance = None
        possible_keys = ['q05_cell_abundance_w_sf', 'means_cell_abundance_w_sf', 'q50_cell_abundance_w_sf']

        for key in possible_keys:
            if key in sp.obsm:
                cell_abundance = sp.obsm[key]
                print(f"Using cell abundance from key: {key}")
                break

        if cell_abundance is None:
            # List available keys for debugging
            available_keys = list(sp.obsm.keys())
            raise RuntimeError(f"Cell2location did not produce expected output. Available keys: {available_keys}")

        # Create DataFrame with results
        proportions = pd.DataFrame(
            cell_abundance,
            index=sp.obs_names,
            columns=ref_signatures.columns
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

        # Create standardized statistics using helper function
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            "Cell2location",
            device,
            n_epochs=n_epochs,
            n_cells_per_spot=n_cells_per_spot,
            use_gpu=use_gpu
        )

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

            # Re-raise the error since we don't have fallback
            raise
        else:
            # Re-raise the error since we don't have fallback
            raise





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
    # Unified validation and gene finding
    common_genes = _validate_deconvolution_inputs(spatial_adata, reference_adata, cell_type_key, min_common_genes)

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

    try:
        # Unified device selection
        device = _get_device(use_gpu, 'Spotiphy')

        # Subset to common genes
        spatial_data = spatial_adata[:, common_genes].copy()
        reference_data = reference_adata[:, common_genes].copy()

        # Ensure data is in the right format for Spotiphy
        # Convert to dense arrays if sparse
        if hasattr(spatial_data.X, 'toarray'):
            spatial_data.X = spatial_data.X.toarray()
        if hasattr(reference_data.X, 'toarray'):
            reference_data.X = reference_data.X.toarray()

        # Ensure non-negative values
        spatial_data.X = np.maximum(spatial_data.X, 0)
        reference_data.X = np.maximum(reference_data.X, 0)

        # Normalize data using Spotiphy's initialization function
        # This performs CPM normalization and filtering
        sc_ref_data, spatial_data = sc_reference.initialization(
            reference_data,
            spatial_data,
            filtering=False,  # We already filtered for common genes
            verbose=1
        )

        # Get cell types from the processed reference data
        type_list = sorted(list(sc_ref_data.obs[cell_type_key].unique()))

        print(f"Found {len(type_list)} cell types: {type_list}")

        # Construct reference profiles using Spotiphy's method
        sc_ref = sc_reference.construct_sc_ref(sc_ref_data, cell_type_key)

        print(f"Constructed reference profiles with shape: {sc_ref.shape}")

        # Extract expression matrices - ensure it's a numpy array
        X = spatial_data.X
        if isinstance(X, np.matrix):
            X = X.A
        elif hasattr(X, 'toarray'):
            X = X.toarray()

        # Ensure X is float64 for compatibility with PyTorch
        X = X.astype(np.float64)

        # Set Adam parameters
        if adam_params is None:
            adam_params = {"lr": 0.003, "betas": (0.95, 0.999)}

        # Run Spotiphy deconvolution
        print(f"Running Spotiphy deconvolution with {len(common_genes)} common genes and {len(type_list)} cell types")
        print(f"Device: {device}, Epochs: {n_epochs}")

        # Estimate cell proportions using the core deconvolute function
        print("Running Spotiphy deconvolute function...")

        # Run the core deconvolution
        pyro_params = deconvolution.deconvolute(
            X,
            sc_ref,
            device=device,
            n_epoch=n_epochs,
            adam_params=adam_params,
            batch_prior=batch_prior,
            plot=False
        )

        # Extract sigma parameter and convert to cell proportions
        sigma = pyro_params['sigma'].cpu().detach().numpy()

        # Calculate mean expression per cell type from reference data
        Y = np.array(sc_ref_data.X)
        if hasattr(sc_ref_data.X, 'toarray'):
            Y = sc_ref_data.X.toarray()

        mean_exp = np.array([np.mean(np.sum(Y[sc_ref_data.obs[cell_type_key]==cell_type], axis=1))
                             for cell_type in type_list])

        # Avoid division by zero
        mean_exp = np.maximum(mean_exp, 1e-8)

        # Convert to cell proportions
        cell_proportions = sigma / mean_exp
        cell_proportions = cell_proportions / np.sum(cell_proportions, axis=1)[:, np.newaxis]

        # Create DataFrame with results
        proportions = pd.DataFrame(
            cell_proportions,
            index=spatial_data.obs_names,
            columns=type_list
        )

        # Create standardized statistics using helper function
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            "Spotiphy",
            device,
            n_epochs=n_epochs,
            batch_prior=batch_prior
        )

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

        # Load and validate reference data ONCE for methods that need it
        reference_adata = None
        if params.method in ["cell2location", "spotiphy"]:
            if not params.reference_data_id:
                raise ValueError(f"Reference data is required for method '{params.method}'. Please provide reference_data_id.")
            
            if params.reference_data_id not in data_store:
                raise ValueError(
                    f"Reference dataset '{params.reference_data_id}' not found in data store. "
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

            # Check cell type key
            if params.cell_type_key not in reference_adata.obs:
                raise ValueError(
                    f"Cell type key '{params.cell_type_key}' not found in reference data. "
                    f"Available keys: {list(reference_adata.obs.columns)}"
                )

            if context:
                await context.info(f"Reference dataset shape: {reference_adata.shape}")
                cell_types = reference_adata.obs[params.cell_type_key].unique()
                await context.info(f"Using reference with {len(cell_types)} cell types: {list(cell_types)}")

        # Run deconvolution with cleaner calls and centralized error handling
        proportions, stats = None, None
        
        try:
            if params.method == "cell2location":
                if context:
                    await context.info("Running Cell2location deconvolution")
                proportions, stats = deconvolve_cell2location(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.n_epochs,
                    n_cells_per_spot=params.n_cells_per_spot or 10,
                    use_gpu=params.use_gpu
                )

            elif params.method == "spotiphy":
                if context:
                    await context.info("Running Spotiphy deconvolution")
                
                # Check if Spotiphy is available
                is_available, error_message = is_spotiphy_available()
                if not is_available:
                    if context:
                        await context.warning(f"Spotiphy is not available: {error_message}")
                    raise ImportError(f"Spotiphy is not available: {error_message}")
                
                # Set Adam parameters
                adam_params = {
                    "lr": params.spotiphy_adam_lr,
                    "betas": params.spotiphy_adam_betas
                }

                proportions, stats = deconvolve_spotiphy(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.n_epochs,
                    batch_prior=params.spotiphy_batch_prior,
                    adam_params=adam_params,
                    use_gpu=params.use_gpu
                )

            else:
                raise ValueError(
                    f"Unsupported deconvolution method: {params.method}. "
                    f"Supported methods are: cell2location, spotiphy"
                )

        except Exception as e:
            # Centralized error handling for deconvolution calls
            error_msg = str(e)
            tb = traceback.format_exc()
            if context:
                await context.warning(f"{params.method.capitalize()} failed: {error_msg}")
            raise RuntimeError(f"{params.method.capitalize()} deconvolution failed: {error_msg}\n{tb}") from e

        if context:
            await context.info(f"{params.method.capitalize()} completed successfully. Found {stats['n_cell_types']} cell types.")

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

        # Create visualization using unified visualization system
        if context:
            await context.info("Creating visualization using unified visualization system")

        try:
            # Import unified visualization
            from .visualization import visualize_data
            from ..models.data import VisualizationParameters

            # Create visualization parameters
            viz_params = VisualizationParameters(
                plot_type="deconvolution",
                n_cell_types=min(6, proportions.shape[1]),
                title=f"{params.method.upper()} Cell Type Proportions",
                colormap="viridis",
                figure_size=(12, 8),
                dpi=100
            )

            # Create temporary data store
            temp_data_store = {data_id: {"adata": spatial_adata}}

            # Generate visualization
            img = await visualize_data(data_id, temp_data_store, viz_params, context)

            if context:
                await context.info("Visualization created successfully using unified system")
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
                "plot_type": "deconvolution",
                "n_cell_types": min(6, proportions.shape[1]),
                "title": f"{params.method.upper()} Cell Type Proportions",
                "colormap": "viridis",
                "figure_size": [12, 8],
                "dpi": 100
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
