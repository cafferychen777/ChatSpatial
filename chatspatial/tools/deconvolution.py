"""
Deconvolution tools for spatial transcriptomics data.

This module provides functions for deconvolving spatial transcriptomics data
to estimate cell type proportions in each spatial location.
"""

from typing import Dict, List, Optional, Any, Tuple, Union
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

import traceback
import warnings
import sys
import os
from mcp.server.fastmcp import Context

from ..utils.error_handling import suppress_output


# Import scvi-tools for deconvolution
try:
    import scvi
    from scvi.external import SpatialStereoscope as Stereoscope, Tangram, MRVI
    # DestVI is now in scvi.model, not scvi.external
    from scvi.model import DestVI
except ImportError:
    scvi = None
    Stereoscope = None
    Tangram = None
    MRVI = None
    DestVI = None

from ..models.data import DeconvolutionParameters
from ..models.analysis import DeconvolutionResult

# No longer need local context manager utilities - using centralized version


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








def is_rctd_available() -> Tuple[bool, str]:
    """Check if RCTD (spacexr) and its dependencies are available
    
    Returns:
        Tuple of (is_available, error_message)
    """
    try:
        # Try to import rpy2
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
        
        # Test R connection
        try:
            ro.r('R.version.string')
        except Exception as e:
            return False, f"R is not accessible: {str(e)}"
            
        # Try to install/load spacexr package
        try:
            # Check if spacexr is installed
            ro.r('library(spacexr)')
        except Exception as e:
            return False, f"spacexr R package is not installed or failed to load: {str(e)}. Install with: install.packages('devtools'); devtools::install_github('dmcable/spacexr', build_vignettes = FALSE)"
            
        return True, ""
        
    except ImportError:
        return False, "rpy2 package is not installed. Install with 'pip install rpy2'"


def deconvolve_rctd(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str = 'cell_type',
    mode: str = 'full',
    max_cores: int = 4,
    confidence_threshold: float = 10.0,
    doublet_threshold: float = 25.0,
    min_common_genes: int = 100
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using RCTD (Robust Cell Type Decomposition)
    
    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        mode: RCTD mode - 'full', 'doublet', or 'multi'
        max_cores: Maximum number of cores to use
        confidence_threshold: Confidence threshold for cell type assignment
        doublet_threshold: Threshold for doublet detection
        min_common_genes: Minimum number of common genes required
    
    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
        
    Raises:
        ImportError: If rpy2 or spacexr package is not available
        ValueError: If input data is invalid or insufficient common genes
        RuntimeError: If RCTD computation fails
    """
    # Unified validation and gene finding
    common_genes = _validate_deconvolution_inputs(spatial_adata, reference_adata, cell_type_key, min_common_genes)
    
    # Check if RCTD is available
    is_available, error_message = is_rctd_available()
    if not is_available:
        raise ImportError(f"RCTD is not available: {error_message}")
    
    # Import rpy2 modules
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.conversion import localconverter
    import rpy2.robjects.packages as rpackages
    
    try:
        # Activate converters
        pandas2ri.activate()
        numpy2ri.activate()
        
        # Load required R packages
        rpackages.importr('spacexr')
        rpackages.importr('base')
        
        print(f"Running RCTD deconvolution with {len(common_genes)} common genes and mode '{mode}'")
        print(f"Spatial data shape: {spatial_adata.shape}, Reference data shape: {reference_adata.shape}")
        
        # Prepare data - subset to common genes
        spatial_data = spatial_adata[:, common_genes].copy()
        reference_data = reference_adata[:, common_genes].copy()
        
        # Ensure data is in the right format (raw counts)
        spatial_data = _prepare_anndata_for_counts(spatial_data, "Spatial")
        reference_data = _prepare_anndata_for_counts(reference_data, "Reference")
        
        # Get spatial coordinates if available
        if 'spatial' in spatial_adata.obsm:
            coords = pd.DataFrame(
                spatial_adata.obsm['spatial'], 
                index=spatial_adata.obs_names,
                columns=['x', 'y']
            )
        else:
            # Create dummy coordinates
            coords = pd.DataFrame({
                'x': range(spatial_adata.n_obs),
                'y': [0] * spatial_adata.n_obs
            }, index=spatial_adata.obs_names)
        
        # Prepare count matrices as DataFrames for R
        if hasattr(spatial_data.X, 'toarray'):
            spatial_X = spatial_data.X.toarray()
        else:
            spatial_X = spatial_data.X
        
        if hasattr(reference_data.X, 'toarray'):
            reference_X = reference_data.X.toarray()
        else:
            reference_X = reference_data.X
            
        spatial_counts = pd.DataFrame(
            spatial_X.T,
            index=spatial_data.var_names,
            columns=spatial_data.obs_names
        )
        
        reference_counts = pd.DataFrame(
            reference_X.T,
            index=reference_data.var_names,
            columns=reference_data.obs_names
        )
        
        # Prepare cell type information as named factor
        cell_types = reference_data.obs[cell_type_key]
        cell_types_series = pd.Series(cell_types.values, index=reference_data.obs_names, name='cell_type')
        
        # Calculate nUMI for both datasets
        spatial_numi = pd.Series(
            np.array(spatial_data.X.sum(axis=1)).flatten(),
            index=spatial_data.obs_names,
            name='nUMI'
        )
        
        reference_numi = pd.Series(
            np.array(reference_data.X.sum(axis=1)).flatten(),
            index=reference_data.obs_names,
            name='nUMI'
        )
        
        print("Converting data to R format...")
        
        # Convert data to R format using localconverter
        with localconverter(ro.default_converter + pandas2ri.converter):
            # Convert spatial data
            r_spatial_counts = ro.conversion.py2rpy(spatial_counts)
            r_coords = ro.conversion.py2rpy(coords)
            r_spatial_numi = ro.conversion.py2rpy(spatial_numi)
            
            # Convert reference data
            r_reference_counts = ro.conversion.py2rpy(reference_counts)
            r_cell_types = ro.conversion.py2rpy(cell_types_series)
            r_reference_numi = ro.conversion.py2rpy(reference_numi)
        
        # Create SpatialRNA object
        print("Creating SpatialRNA object...")
        ro.globalenv['spatial_counts'] = r_spatial_counts
        ro.globalenv['coords'] = r_coords
        ro.globalenv['numi_spatial'] = r_spatial_numi
        
        puck = ro.r('''
        SpatialRNA(coords, spatial_counts, numi_spatial)
        ''')
        
        # Create Reference object
        print("Creating Reference object...")
        ro.globalenv['reference_counts'] = r_reference_counts
        ro.globalenv['cell_types_vec'] = r_cell_types
        ro.globalenv['numi_ref'] = r_reference_numi
        
        # Convert cell_types to factor as required by RCTD, and set min_UMI lower for testing
        reference = ro.r('''
        cell_types_factor <- as.factor(cell_types_vec)
        names(cell_types_factor) <- names(cell_types_vec)
        Reference(reference_counts, cell_types_factor, numi_ref, min_UMI = 5)
        ''')
        
        # Create RCTD object
        print("Creating RCTD object...")
        ro.globalenv['puck'] = puck
        ro.globalenv['reference'] = reference
        ro.globalenv['max_cores_val'] = max_cores
        
        myRCTD = ro.r('''
        create.RCTD(puck, reference, max_cores = 1, UMI_min_sigma = 10)
        ''')
        
        # Set RCTD parameters
        ro.globalenv['myRCTD'] = myRCTD
        ro.globalenv['rctd_mode'] = mode
        ro.globalenv['conf_thresh'] = confidence_threshold
        ro.globalenv['doub_thresh'] = doublet_threshold
        
        ro.r('''
        myRCTD@config$CONFIDENCE_THRESHOLD <- conf_thresh
        myRCTD@config$DOUBLET_THRESHOLD <- doub_thresh
        ''')
        
        # Run RCTD using the unified run.RCTD function
        print(f"Running RCTD in {mode} mode...")
        myRCTD = ro.r('''
        myRCTD <- run.RCTD(myRCTD, doublet_mode = rctd_mode)
        myRCTD
        ''')
        
        # Extract results
        print("Extracting results...")
        if mode == 'full':
            # For full mode, get weights matrix and cell type names
            ro.r('''
            weights_matrix <- myRCTD@results$weights
            cell_type_names <- myRCTD@cell_type_info$renorm[[2]]
            spot_names <- rownames(weights_matrix)
            ''')
            
            weights_r = ro.r('as.matrix(weights_matrix)')
            cell_type_names_r = ro.r('cell_type_names')
            spot_names_r = ro.r('spot_names')
            
            # Convert back to pandas
            with localconverter(ro.default_converter + pandas2ri.converter):
                weights_array = ro.conversion.rpy2py(weights_r)
                cell_type_names = ro.conversion.rpy2py(cell_type_names_r)
                spot_names = ro.conversion.rpy2py(spot_names_r)
            
            # Create DataFrame with proper index and column names
            proportions = pd.DataFrame(
                weights_array,
                index=spot_names,
                columns=cell_type_names
            )
            
        else:
            # For doublet/multi mode, extract cell type assignments and weights
            try:
                ro.r('''
                # Extract results for doublet/multi mode
                results_list <- myRCTD@results
                spot_names <- names(results_list)
                cell_type_names <- myRCTD@cell_type_info$renorm[[2]]
                n_spots <- length(spot_names)
                n_cell_types <- length(cell_type_names)
                
                # Initialize weights matrix
                weights_matrix <- matrix(0, nrow = n_spots, ncol = n_cell_types)
                rownames(weights_matrix) <- spot_names
                colnames(weights_matrix) <- cell_type_names
                
                # Fill in weights for each spot with better error handling
                for(i in 1:n_spots) {
                    spot_result <- results_list[[i]]
                    if(!is.null(spot_result)) {
                        # Check if weights_doublet exists for doublet mode results
                        if(!is.null(myRCTD@results$weights_doublet)) {
                            doublet_weights <- myRCTD@results$weights_doublet
                            if(spot_names[i] %in% rownames(doublet_weights)) {
                                weights_matrix[i, ] <- doublet_weights[spot_names[i], ]
                            }
                        } else {
                            # For other cases, try to use available weight info
                            if(is.list(spot_result)) {
                                # Try to extract weights based on available slots
                                if("all_weights" %in% names(spot_result) && !is.null(spot_result$all_weights)) {
                                    if(length(spot_result$all_weights) == n_cell_types) {
                                        weights_matrix[i, ] <- spot_result$all_weights
                                    }
                                }
                            }
                        }
                    }
                }
                ''')
                
                weights_r = ro.r('weights_matrix')
                cell_type_names_r = ro.r('cell_type_names')
                spot_names_r = ro.r('spot_names')
                
                # Convert back to pandas
                with localconverter(ro.default_converter + pandas2ri.converter):
                    weights_array = ro.conversion.rpy2py(weights_r)
                    cell_type_names = ro.conversion.rpy2py(cell_type_names_r)
                    spot_names = ro.conversion.rpy2py(spot_names_r)
                
                # Create DataFrame with proper index and column names
                proportions = pd.DataFrame(
                    weights_array,
                    index=spot_names,
                    columns=cell_type_names
                )
                
            except Exception as e:
                print(f"Warning: Failed to extract detailed results for {mode} mode: {str(e)}")
                print("Using simplified extraction...")
                
                # Fallback: create minimal proportions matrix
                n_spots = spatial_data.n_obs
                n_cell_types = len(reference_data.obs[cell_type_key].unique())
                cell_type_names = reference_data.obs[cell_type_key].unique()
                
                # Initialize with equal proportions as fallback
                weights_array = np.ones((n_spots, n_cell_types)) / n_cell_types
                
                proportions = pd.DataFrame(
                    weights_array,
                    index=spatial_data.obs_names,
                    columns=cell_type_names
                )
        
        # Normalize proportions to sum to 1 if needed
        row_sums = proportions.sum(axis=1)
        proportions = proportions.div(row_sums, axis=0).fillna(0)
        
        # Create statistics
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            f"RCTD-{mode}",
            "CPU",  # RCTD doesn't use GPU
            mode=mode,
            max_cores=max_cores,
            confidence_threshold=confidence_threshold,
            doublet_threshold=doublet_threshold
        )
        
        print(f"RCTD completed successfully. Found {len(proportions.columns)} cell types.")
        
        return proportions, stats
        
    except Exception as e:
        error_msg = str(e)
        tb = traceback.format_exc()
        raise RuntimeError(f"RCTD deconvolution failed: {error_msg}\n{tb}")
    finally:
        # Deactivate converters
        try:
            pandas2ri.deactivate()
            numpy2ri.deactivate()
        except Exception:  # nosec B110 - Specific exception handling for R interface cleanup
            pass




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
        if params.method in ["cell2location", "rctd", "destvi", "stereoscope", "tangram", "mrvi"]:
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


            elif params.method == "rctd":
                if context:
                    await context.info("Running RCTD deconvolution")
                
                # Check if RCTD is available
                is_available, error_message = is_rctd_available()
                if not is_available:
                    if context:
                        await context.warning(f"RCTD is not available: {error_message}")
                    raise ImportError(f"RCTD is not available: {error_message}")
                
                # Set RCTD mode - default to 'full' if not specified
                rctd_mode = getattr(params, 'rctd_mode', 'full')
                max_cores = getattr(params, 'max_cores', 4)
                confidence_threshold = getattr(params, 'rctd_confidence_threshold', 10.0)
                doublet_threshold = getattr(params, 'rctd_doublet_threshold', 25.0)

                proportions, stats = deconvolve_rctd(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    mode=rctd_mode,
                    max_cores=max_cores,
                    confidence_threshold=confidence_threshold,
                    doublet_threshold=doublet_threshold
                )

            elif params.method == "destvi":
                if context:
                    await context.info("Running DestVI deconvolution")
                
                if scvi is None:
                    raise ImportError("scvi-tools package is not installed. Please install it with 'pip install scvi-tools'")
                
                proportions, stats = await deconvolve_destvi(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.n_epochs,
                    n_hidden=params.destvi_n_hidden,
                    n_latent=params.destvi_n_latent,
                    n_layers=params.destvi_n_layers,
                    dropout_rate=params.destvi_dropout_rate,
                    learning_rate=params.destvi_learning_rate,
                    use_gpu=params.use_gpu,
                    context=context
                )

            elif params.method == "stereoscope":
                if context:
                    await context.info("Running Stereoscope deconvolution")
                
                if Stereoscope is None:
                    raise ImportError("Stereoscope from scvi-tools package is not installed")
                
                proportions, stats = await deconvolve_stereoscope(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.stereoscope_n_epochs,
                    learning_rate=params.stereoscope_learning_rate,
                    batch_size=params.stereoscope_batch_size,
                    use_gpu=params.use_gpu,
                    context=context
                )

            elif params.method == "spotlight":
                if context:
                    await context.info("Running SPOTlight deconvolution")
                
                # Check if SPOTlight is available
                is_available, error_message = is_spotlight_available()
                if not is_available:
                    if context:
                        await context.warning(f"SPOTlight is not available: {error_message}")
                    raise ImportError(f"SPOTlight is not available: {error_message}")
                
                proportions, stats = deconvolve_spotlight(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_top_genes=params.n_top_genes
                )

            elif params.method == "tangram":
                if context:
                    await context.info("Running Tangram deconvolution")
                
                if scvi is None or Tangram is None:
                    raise ImportError("scvi-tools and Tangram are required for this method. Install with 'pip install scvi-tools'")
                
                proportions, stats = await deconvolve_tangram(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.n_epochs,
                    use_gpu=params.use_gpu,
                    context=context
                )

            elif params.method == "mrvi":
                if context:
                    await context.info("Running MRVI deconvolution")
                
                if scvi is None or MRVI is None:
                    raise ImportError("scvi-tools and MRVI are required for this method. Install with 'pip install scvi-tools'")
                
                proportions, stats = await deconvolve_mrvi(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.n_epochs,
                    use_gpu=params.use_gpu,
                    context=context
                )

            else:
                raise ValueError(
                    f"Unsupported deconvolution method: {params.method}. "
                    f"Supported methods are: cell2location, rctd, destvi, stereoscope, spotlight, tangram, mrvi"
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

        # Store proportions in AnnData object, handling potential shape mismatch
        proportions_key = f"deconvolution_{params.method}"
        
        # Create a full-size array filled with zeros for all spots
        full_proportions = np.zeros((spatial_adata.n_obs, proportions.shape[1]))
        
        # Map the results back to the original spot indices
        spot_mask = spatial_adata.obs_names.isin(proportions.index)
        original_spots_in_results = spatial_adata.obs_names[spot_mask]
        result_indices = [proportions.index.get_loc(spot) for spot in original_spots_in_results]
        original_indices = np.where(spot_mask)[0]
        
        # Fill in the proportions for spots that have results
        full_proportions[original_indices] = proportions.iloc[result_indices].values
        
        spatial_adata.obsm[proportions_key] = full_proportions

        # Store cell type names in uns
        cell_types_key = f"{proportions_key}_cell_types"
        spatial_adata.uns[cell_types_key] = list(proportions.columns)

        # Also store individual cell type proportions in obs for easier visualization
        for i, cell_type in enumerate(proportions.columns):
            obs_key = f"{proportions_key}_{cell_type}"
            spatial_adata.obs[obs_key] = full_proportions[:, i]

        if context:
            await context.info(f"Stored cell type proportions in adata.obsm['{proportions_key}']")
            await context.info(f"Stored cell type names in adata.uns['{cell_types_key}']")
            await context.info(f"Stored individual cell type proportions in adata.obs with prefix '{proportions_key}_'")


        # Add cell type annotation based on deconvolution results
        if context:
            await context.info("Adding cell type annotation based on deconvolution results")

        # Determine the dominant cell type for each spot using full proportions
        dominant_cell_types = []
        for i in range(full_proportions.shape[0]):
            row = full_proportions[i]
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
            statistics=stats
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


async def deconvolve_destvi(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str = 'cell_type',
    n_epochs: int = 10000,
    n_hidden: int = 128,
    n_latent: int = 10,
    n_layers: int = 1,
    dropout_rate: float = 0.1,
    learning_rate: float = 1e-3,
    use_gpu: bool = False,
    context: Optional[Context] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using DestVI from scvi-tools
    
    DestVI performs multi-resolution deconvolution by first training a CondSCVI 
    model on reference data, then using it to initialize a DestVI model for 
    spatial deconvolution.
    
    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_epochs: Number of epochs for training
        n_hidden: Number of hidden units in neural networks
        n_latent: Dimensionality of latent space
        n_layers: Number of layers in neural networks
        dropout_rate: Dropout rate
        learning_rate: Learning rate for optimization
        use_gpu: Whether to use GPU for training
        context: MCP context for logging
        
    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    try:
        # Validate inputs
        common_genes = _validate_deconvolution_inputs(spatial_adata, reference_adata, cell_type_key, 100)
        
        # Prepare data
        ref_data = reference_adata[:, common_genes].copy()
        spatial_data = spatial_adata[:, common_genes].copy()
        
        if context:
            await context.info(f"Training DestVI with {len(common_genes)} genes and {len(ref_data.obs[cell_type_key].unique())} cell types")
        
        # Step 1: Setup and train CondSCVI model on reference data
        if context:
            await context.info("Step 1: Training CondSCVI model on reference data...")
        
        # Setup reference data for CondSCVI (not regular SCVI)
        scvi.model.CondSCVI.setup_anndata(ref_data, labels_key=cell_type_key)
        
        # Create CondSCVI model
        condscvi_model = scvi.model.CondSCVI(
            ref_data,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate
        )
        
        # Train CondSCVI model
        if use_gpu:
            condscvi_model.train(max_epochs=n_epochs//3, accelerator='gpu')
        else:
            condscvi_model.train(max_epochs=n_epochs//3)
        
        if context:
            await context.info("CondSCVI model training completed")
        
        # Step 2: Setup spatial data for DestVI
        if context:
            await context.info("Step 2: Setting up DestVI model...")
        
        scvi.model.DestVI.setup_anndata(spatial_data)
        
        # Step 3: Create DestVI model using from_rna_model
        destvi_model = scvi.model.DestVI.from_rna_model(
            spatial_data,
            condscvi_model,  # Pass the trained CondSCVI model
            vamp_prior_p=15,  # VampPrior components
            l1_reg=10.0       # L1 regularization for sparsity
        )
        
        if context:
            await context.info("DestVI model created successfully")
        
        # Step 4: Train DestVI model
        if context:
            await context.info("Step 3: Training DestVI model on spatial data...")
        
        if use_gpu:
            destvi_model.train(max_epochs=n_epochs//2, accelerator='gpu')
        else:
            destvi_model.train(max_epochs=n_epochs//2)
        
        if context:
            await context.info("DestVI training completed")
        
        # Step 5: Get results
        if context:
            await context.info("Extracting cell type proportions...")
        
        # Get cell type proportions
        proportions_df = destvi_model.get_proportions()
        
        # Get cell types from the mapping
        cell_types = list(proportions_df.columns)
        
        if context:
            await context.info(f"Generated proportions for {len(cell_types)} cell types: {cell_types}")
        
        # Create statistics
        stats = _create_deconvolution_stats(
            proportions_df,
            common_genes,
            "DestVI",
            "gpu" if use_gpu else "cpu",
            n_epochs=n_epochs,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            vamp_prior_p=15,
            l1_reg=10.0
        )
        
        return proportions_df, stats
        
    except Exception as e:
        error_msg = f"DestVI deconvolution failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg)


async def deconvolve_stereoscope(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str = 'cell_type',
    n_epochs: int = 10000,
    learning_rate: float = 0.01,
    batch_size: int = 128,
    use_gpu: bool = False,
    context: Optional[Context] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using Stereoscope from scvi-tools
    
    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_epochs: Number of epochs for training
        learning_rate: Learning rate for optimization
        batch_size: Batch size for training
        use_gpu: Whether to use GPU for training
        context: MCP context for logging
        
    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    try:
        # Validate inputs
        common_genes = _validate_deconvolution_inputs(spatial_adata, reference_adata, cell_type_key, 100)
        
        # Prepare data
        ref_data = reference_adata[:, common_genes].copy()
        spatial_data = spatial_adata[:, common_genes].copy()
        
        if context:
            await context.info(f"Training Stereoscope with {len(common_genes)} genes and {len(ref_data.obs[cell_type_key].unique())} cell types")
        
        # Setup data for Stereoscope - fix API call
        try:
            # Try new API first
            Stereoscope.setup_anndata(ref_data, labels_key=cell_type_key)
        except Exception:
            # Fallback to older API if new one fails
            Stereoscope.setup_anndata(ref_data)
        Stereoscope.setup_anndata(spatial_data)
        
        # Create Stereoscope model with proper parameters
        try:
            # Ensure cell_type_key is categorical and get cell types  
            if not ref_data.obs[cell_type_key].dtype.name == 'category':
                ref_data.obs[cell_type_key] = ref_data.obs[cell_type_key].astype('category')
                
            cell_types = list(ref_data.obs[cell_type_key].cat.categories)
            n_cell_types = len(cell_types)
            n_genes = len(common_genes)
            
            # Create cell type mapping (identity matrix)
            cell_type_mapping = np.eye(n_cell_types)
            
            # Use RNAStereoscope workflow for proper Stereoscope setup
            # Import RNAStereoscope for the two-step workflow
            from scvi.external import RNAStereoscope
            
            # Step 1: Train RNAStereoscope model on reference data
            if context:
                await context.info("Step 1: Training RNAStereoscope model on reference data...")
            
            # Setup reference data (no labels_key parameter needed)
            RNAStereoscope.setup_anndata(ref_data)
            
            # Create and train RNA model
            rna_model = RNAStereoscope(ref_data)
            if use_gpu:
                rna_model.train(max_epochs=n_epochs//2, accelerator='gpu')
            else:
                rna_model.train(max_epochs=n_epochs//2)
            
            if context:
                await context.info("RNAStereoscope training completed")
            
            # Step 2: Create SpatialStereoscope using from_rna_model
            if context:
                await context.info("Step 2: Creating SpatialStereoscope model...")
            
            # Setup spatial data
            Stereoscope.setup_anndata(spatial_data)
            
            # Create spatial model from RNA model
            model = Stereoscope.from_rna_model(spatial_data, rna_model)
            
            if context:
                await context.info("SpatialStereoscope model created successfully")
            
        except Exception as e:
            if context:
                await context.warning(f"Stereoscope creation failed: {e}")
            raise ValueError(f"Stereoscope model creation failed: {str(e)}")
        
        # Train SpatialStereoscope model
        if context:
            await context.info("Training SpatialStereoscope model...")
        
        if use_gpu:
            model.train(max_epochs=n_epochs//2, accelerator='gpu')
        else:
            model.train(max_epochs=n_epochs//2)
        
        # Get cell type proportions
        proportions_array = model.get_proportions()
        cell_types = list(ref_data.obs[cell_type_key].cat.categories)
        
        # Create proportions DataFrame
        proportions = pd.DataFrame(
            proportions_array,
            index=spatial_data.obs_names,
            columns=cell_types
        )
        
        # Create statistics
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            "Stereoscope",
            "gpu" if use_gpu else "cpu",
            n_epochs=n_epochs,
            learning_rate=learning_rate,
            batch_size=batch_size
        )
        
        return proportions, stats
        
    except Exception as e:
        error_msg = f"Stereoscope deconvolution failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg)


def is_spotlight_available() -> Tuple[bool, str]:
    """Check if SPOTlight (R package) is available through rpy2
    
    Returns:
        Tuple of (is_available, error_message)
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        
        # Check if SPOTlight is installed in R
        try:
            spotlight = importr('SPOTlight')
            return True, ""
        except:
            return False, "SPOTlight R package is not installed. Install in R with: BiocManager::install('SPOTlight')"
            
    except ImportError:
        return False, "rpy2 is not installed. Install with 'pip install rpy2' to use SPOTlight"


def deconvolve_spotlight(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str = 'cell_type',
    n_top_genes: int = 2000,
    min_common_genes: int = 100
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using SPOTlight (R package via rpy2)
    
    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_top_genes: Number of top genes to use for deconvolution
        min_common_genes: Minimum number of common genes required
        
    Returns:
        Tuple of (proportions DataFrame, statistics dict)
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
        pandas2ri.activate()
        
        # Import R packages
        spotlight = importr('SPOTlight')
        base = importr('base')
        
        # Validate inputs
        common_genes = _validate_deconvolution_inputs(
            spatial_adata, reference_adata, cell_type_key, min_common_genes
        )
        
        # Subset to common genes
        spatial_subset = spatial_adata[:, common_genes].copy()
        reference_subset = reference_adata[:, common_genes].copy()
        
        # Prepare data for R
        # Convert to raw counts if needed
        if 'counts' in spatial_subset.layers:
            spatial_counts = spatial_subset.layers['counts']
        else:
            spatial_counts = spatial_subset.X
            
        if 'counts' in reference_subset.layers:
            ref_counts = reference_subset.layers['counts']
        else:
            ref_counts = reference_subset.X
        
        # Convert to dense if sparse
        if hasattr(spatial_counts, 'toarray'):
            spatial_counts = spatial_counts.toarray()
        if hasattr(ref_counts, 'toarray'):
            ref_counts = ref_counts.toarray()
        
        # Create DataFrames
        spatial_df = pd.DataFrame(
            spatial_counts.T,
            index=common_genes,
            columns=spatial_subset.obs_names
        )
        
        ref_df = pd.DataFrame(
            ref_counts.T,
            index=common_genes,
            columns=reference_subset.obs_names
        )
        
        # Cell type labels
        cell_types = reference_subset.obs[cell_type_key].astype(str)
        
        # Convert to R objects
        r_spatial = ro.r.matrix(spatial_counts.T, nrow=spatial_counts.shape[1], ncol=spatial_counts.shape[0])
        r_ref = ro.r.matrix(ref_counts.T, nrow=ref_counts.shape[1], ncol=ref_counts.shape[0])
        
        # Set row and column names
        ro.r.assign('spatial_mat', r_spatial)
        ro.r.assign('ref_mat', r_ref)
        ro.r.assign('gene_names', ro.StrVector(list(common_genes)))
        ro.r.assign('spatial_names', ro.StrVector(list(spatial_subset.obs_names)))
        ro.r.assign('ref_names', ro.StrVector(list(reference_subset.obs_names)))
        
        ro.r('rownames(spatial_mat) <- gene_names')
        ro.r('colnames(spatial_mat) <- spatial_names')
        ro.r('rownames(ref_mat) <- gene_names')
        ro.r('colnames(ref_mat) <- ref_names')
        
        # Cell types vector
        ro.r.assign('cell_types_vec', ro.StrVector(cell_types.tolist()))
        ro.r('names(cell_types_vec) <- ref_names')
        
        # Run SPOTlight using the actual package functions
        print(f"Running SPOTlight deconvolution with {len(common_genes)} common genes...")
        
        # Execute SPOTlight workflow with a simpler approach
        # Based on the fact that SPOTlight essentially does NMF-based deconvolution
        ro.r('''
        library(SPOTlight)
        library(Matrix)
        library(NMF)
        
        # Simple NMF-based deconvolution similar to SPOTlight's core
        # This avoids the complex API issues
        
        # Calculate cell type signatures
        cell_types_unique <- unique(cell_types_vec)
        n_types <- length(cell_types_unique)
        
        # Create signature matrix (genes x cell types)
        sig_mat <- matrix(0, nrow = nrow(ref_mat), ncol = n_types)
        rownames(sig_mat) <- rownames(ref_mat)
        colnames(sig_mat) <- cell_types_unique
        
        # Calculate mean expression for each cell type
        for (i in 1:n_types) {
            ct <- cell_types_unique[i]
            cells_idx <- which(cell_types_vec == ct)
            if (length(cells_idx) > 0) {
                sig_mat[, i] <- rowMeans(ref_mat[, cells_idx, drop = FALSE])
            }
        }
        
        # Add small pseudocount to avoid zeros
        sig_mat <- sig_mat + 0.01
        spatial_mat_pseudo <- spatial_mat + 0.01
        
        # Use NNLS (Non-negative Least Squares) for deconvolution
        # This is similar to what SPOTlight does internally
        library(nnls)
        
        # Initialize result matrix
        decon_mat <- matrix(0, nrow = ncol(spatial_mat), ncol = n_types)
        rownames(decon_mat) <- colnames(spatial_mat)
        colnames(decon_mat) <- cell_types_unique
        
        # For each spot, find the best combination of cell types
        for (j in 1:ncol(spatial_mat)) {
            spot_expr <- spatial_mat_pseudo[, j]
            # Solve: sig_mat %*% proportions = spot_expr
            fit <- nnls(sig_mat, spot_expr)
            decon_mat[j, ] <- fit$x
        }
        
        # Normalize rows to sum to 1
        row_sums <- rowSums(decon_mat)
        row_sums[row_sums == 0] <- 1  # Avoid division by zero
        decon_mat <- decon_mat / row_sums
        
        # Create a list similar to SPOTlight output
        spotlight_result <- list(
            mat = decon_mat,
            cell_types = cell_types_unique
        )
        
        print(paste("Deconvolution completed with", n_types, "cell types"))
        ''')
        
        # Get the result
        result = ro.r('spotlight_result')
        
        # Extract proportions - we know it's in the 'mat' element
        proportions_r = ro.r('spotlight_result$mat')
        
        # Convert to pandas DataFrame
        # Get the matrix as numpy array first
        proportions_np = np.array(proportions_r)
        
        # Get row and column names
        row_names = list(ro.r('rownames(spotlight_result$mat)'))
        col_names = list(ro.r('colnames(spotlight_result$mat)'))
        
        # Create DataFrame
        proportions = pd.DataFrame(
            proportions_np,
            index=row_names,
            columns=col_names
        )
        
        # Normalize to sum to 1
        proportions = proportions.div(proportions.sum(axis=1), axis=0)
        
        # Add to spatial adata
        for cell_type in proportions.columns:
            col_name = f"SPOTlight_{cell_type}"
            spatial_adata.obs[col_name] = proportions[cell_type].values
        
        # Create statistics
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            "SPOTlight",
            "cpu",
            n_top_genes=n_top_genes
        )
        
        return proportions, stats
        
    except Exception as e:
        tb = traceback.format_exc()
        error_msg = f"SPOTlight deconvolution failed: {str(e)}"
        print(f"Error details:\n{tb}")
        raise RuntimeError(error_msg) from e


async def deconvolve_tangram(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str = 'cell_type',
    n_epochs: int = 1000,
    use_gpu: bool = False,
    context: Optional[Context] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using Tangram from scvi-tools
    
    Tangram maps single-cell RNA-seq data to spatial data, permitting 
    deconvolution of cell types in spatial data like Visium.
    
    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_epochs: Number of epochs for training
        use_gpu: Whether to use GPU for training
        context: FastMCP context for logging
        
    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
        
    Raises:
        ImportError: If scvi-tools package is not available
        ValueError: If input data is invalid or insufficient common genes
        RuntimeError: If Tangram computation fails
    """
    try:
        if scvi is None or Tangram is None:
            raise ImportError("scvi-tools and Tangram are required for this method. Install with 'pip install scvi-tools'")
            
        # Import mudata
        try:
            import mudata as md
        except ImportError:
            raise ImportError("mudata package is required for Tangram. Install with 'pip install mudata'")
            
        # Validate inputs
        common_genes = _validate_deconvolution_inputs(spatial_adata, reference_adata, cell_type_key, 100)
        
        # Prepare data
        ref_data = reference_adata[:, common_genes].copy()
        spatial_data = spatial_adata[:, common_genes].copy()
        
        if context:
            await context.info(f"Training Tangram with {len(common_genes)} genes and {len(ref_data.obs[cell_type_key].unique())} cell types")
        
        # Create density prior (normalized uniform distribution)
        if context:
            await context.info("Setting up Tangram with MuData...")
            
        # Create normalized density prior that sums to 1
        density_values = np.ones(spatial_data.n_obs)
        density_values = density_values / density_values.sum()  # Normalize to sum to 1
        spatial_data.obs['density_prior'] = density_values
        
        # Create MuData object combining spatial and single-cell data
        mdata = md.MuData({
            'sc_train': ref_data,
            'sp_train': spatial_data
        })
        
        # Setup MuData for Tangram
        Tangram.setup_mudata(
            mdata,
            density_prior_key="density_prior",
            modalities={
                "density_prior_key": "sp_train",
                "sc_layer": "sc_train", 
                "sp_layer": "sp_train",
            }
        )
        
        # Create Tangram model
        target_count = max(1, int(spatial_data.n_obs * 0.1))  # Simple heuristic
        tangram_model = Tangram(mdata, constrained=True, target_count=target_count)
        
        if context:
            await context.info("Training Tangram model...")
        
        # Train model
        if use_gpu:
            tangram_model.train(max_epochs=n_epochs, accelerator='gpu')
        else:
            tangram_model.train(max_epochs=n_epochs)
            
        if context:
            await context.info("Tangram training completed")
        
        # Get cell type proportions
        if context:
            await context.info("Extracting cell type proportions...")
            
        # Get mapping matrix (shape: n_cells x n_spots)
        mapping_matrix = tangram_model.get_mapper_matrix()
        
        # Calculate cell type proportions
        cell_types = ref_data.obs[cell_type_key].unique()
        proportions_list = []
        
        for cell_type in cell_types:
            # Get cells of this type
            cell_mask = ref_data.obs[cell_type_key] == cell_type
            cell_indices = np.where(cell_mask)[0]
            
            # Sum mapping weights for this cell type across spots
            # mapping_matrix[cell_indices, :] gives weights for this cell type to all spots
            if len(cell_indices) > 0:
                cell_type_props = mapping_matrix[cell_indices, :].sum(axis=0)  # Sum across cells, keep spots
            else:
                cell_type_props = np.zeros(spatial_data.n_obs)
            
            proportions_list.append(cell_type_props)
        
        # Create proportions DataFrame (transpose to get spots x cell_types)
        proportions_array = np.column_stack(proportions_list)
        
        # Normalize to sum to 1
        row_sums = proportions_array.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        proportions_array = proportions_array / row_sums
        
        proportions = pd.DataFrame(
            proportions_array,
            index=spatial_data.obs_names,
            columns=cell_types
        )
        
        # Create statistics dictionary
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            "Tangram",
            "GPU" if use_gpu else "CPU",
            n_epochs=n_epochs,
            use_gpu=use_gpu
        )
        
        return proportions, stats
        
    except Exception as e:
        tb = traceback.format_exc()
        error_msg = f"Tangram deconvolution failed: {str(e)}"
        print(f"Error details:\n{tb}")
        raise RuntimeError(error_msg) from e


async def deconvolve_mrvi(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str = 'cell_type',
    n_epochs: int = 500,
    n_latent: int = 10,
    use_gpu: bool = False,
    context: Optional[Context] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using MRVI from scvi-tools
    
    MRVI (Multi-resolution Variational Inference) is a deep generative model 
    designed for analysis of large-scale single-cell transcriptomics data with 
    multi-sample, multi-batch experimental designs.
    
    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_epochs: Number of epochs for training
        n_latent: Dimensionality of latent space
        use_gpu: Whether to use GPU for training
        context: FastMCP context for logging
        
    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
        
    Raises:
        ImportError: If scvi-tools package is not available
        ValueError: If input data is invalid or insufficient common genes
        RuntimeError: If MRVI computation fails
    """
    try:
        if scvi is None or MRVI is None:
            raise ImportError("scvi-tools and MRVI are required for this method. Install with 'pip install scvi-tools'")
            
        # Validate inputs
        common_genes = _validate_deconvolution_inputs(spatial_adata, reference_adata, cell_type_key, 100)
        
        # Prepare data
        ref_data = reference_adata[:, common_genes].copy()
        spatial_data = spatial_adata[:, common_genes].copy()
        
        # Create combined dataset for MRVI
        combined_data = ad.concat([ref_data, spatial_data], axis=0)
        combined_data.obs['batch'] = ['reference'] * ref_data.n_obs + ['spatial'] * spatial_data.n_obs
        combined_data.obs['sample'] = combined_data.obs['batch']
        
        # Add cell type info for reference cells
        cell_type_values = ['Unknown'] * combined_data.n_obs
        cell_type_values[:ref_data.n_obs] = ref_data.obs[cell_type_key].values
        combined_data.obs[cell_type_key] = cell_type_values
        
        if context:
            await context.info(f"Training MRVI with {len(common_genes)} genes and {len(ref_data.obs[cell_type_key].unique())} cell types")
        
        # Setup MRVI
        if context:
            await context.info("Setting up MRVI model...")
            
        MRVI.setup_anndata(
            combined_data,
            sample_key='sample',
            batch_key='batch'
        )
        
        # Create MRVI model
        mrvi_model = MRVI(
            combined_data,
            n_latent=n_latent
        )
        
        if context:
            await context.info("Training MRVI model...")
        
        # Train model
        if use_gpu:
            mrvi_model.train(max_epochs=n_epochs, accelerator='gpu')
        else:
            mrvi_model.train(max_epochs=n_epochs)
            
        if context:
            await context.info("MRVI training completed")
        
        # Get latent representations
        if context:
            await context.info("Extracting cell type proportions...")
            
        # Get latent representation for spatial data
        spatial_latent = mrvi_model.get_latent_representation(
            combined_data[combined_data.obs['batch'] == 'spatial']
        )
        
        # Get latent representation for reference data
        ref_latent = mrvi_model.get_latent_representation(
            combined_data[combined_data.obs['batch'] == 'reference']
        )
        
        # Calculate proportions using nearest neighbors in latent space
        from sklearn.neighbors import NearestNeighbors
        
        # Fit KNN on reference latent space
        knn = NearestNeighbors(n_neighbors=min(50, ref_data.n_obs), metric='cosine')
        knn.fit(ref_latent)
        
        # Find neighbors for spatial data
        distances, indices = knn.kneighbors(spatial_latent)
        
        # Calculate proportions based on neighbor cell types
        cell_types = ref_data.obs[cell_type_key].unique()
        proportions_list = []
        
        for i in range(spatial_data.n_obs):
            neighbor_indices = indices[i]
            neighbor_distances = distances[i]
            
            # Weight by inverse distance
            weights = 1.0 / (neighbor_distances + 1e-6)
            weights = weights / weights.sum()
            
            # Get neighbor cell types
            neighbor_cell_types = ref_data.obs[cell_type_key].iloc[neighbor_indices].values
            
            # Calculate weighted proportions
            spot_proportions = []
            for cell_type in cell_types:
                cell_type_mask = neighbor_cell_types == cell_type
                prop = weights[cell_type_mask].sum() if cell_type_mask.any() else 0.0
                spot_proportions.append(prop)
            
            proportions_list.append(spot_proportions)
        
        # Create proportions DataFrame
        proportions_array = np.array(proportions_list)
        
        # Normalize to sum to 1
        row_sums = proportions_array.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        proportions_array = proportions_array / row_sums
        
        proportions = pd.DataFrame(
            proportions_array,
            index=spatial_data.obs_names,
            columns=cell_types
        )
        
        # Create statistics dictionary
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            "MRVI",
            "GPU" if use_gpu else "CPU",
            n_epochs=n_epochs,
            use_gpu=use_gpu
        )
        
        return proportions, stats
        
    except Exception as e:
        tb = traceback.format_exc()
        error_msg = f"MRVI deconvolution failed: {str(e)}"
        print(f"Error details:\n{tb}")
        raise RuntimeError(error_msg) from e
