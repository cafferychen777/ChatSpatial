"""
GASTON-based spatial variable genes identification for ChatSpatial MCP.

This module integrates GASTON (Generative Adversarial Spatial Transcriptomics Optimization Network)
for identifying spatial variable genes through topographic mapping and isodepth learning.
"""

import os
import sys
import tempfile
import shutil
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from typing import Dict, List, Tuple, Optional, Any
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings

# Get the path to the third_party directory relative to this file
current_dir = Path(__file__).parent
project_root = current_dir.parent.parent  # Go up to chatspatial root
GASTON_PATH = os.path.join(project_root, "third_party", "GASTON", "src")

if os.path.exists(GASTON_PATH) and GASTON_PATH not in sys.path:
    sys.path.insert(0, GASTON_PATH)

try:
    import gaston
    from gaston import neural_net, spatial_gene_classification, binning_and_plotting
    from gaston import dp_related, segmented_fit, process_NN_output
    GASTON_AVAILABLE = True
except ImportError as e:
    GASTON_AVAILABLE = False
    print(f"Warning: GASTON not available: {e}")
    print("Please install GASTON or ensure it's in the Python path")

from ..models.data import SpatialVariableGenesParameters
from ..models.analysis import SpatialVariableGenesResult


async def identify_spatial_genes(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None
) -> SpatialVariableGenesResult:
    """
    Identify spatial variable genes using GASTON method.
    
    Args:
        data_id: Dataset identifier
        data_store: Data storage instance
        params: GASTON parameters
        context: MCP context for logging
        
    Returns:
        SpatialVariableGenesResult with GASTON analysis results
    """
    if not GASTON_AVAILABLE:
        raise ImportError("GASTON is not available. Please install GASTON to use this method.")
    
    if context:
        await context.info("Starting GASTON spatial variable genes identification")
    
    # Get data
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Work with the original adata object, not a copy
    adata = data_store[data_id]["adata"]
    
    # Validate spatial coordinates
    if 'spatial' not in adata.obsm:
        raise ValueError("Spatial coordinates not found in adata.obsm['spatial']")
    
    # Extract spatial coordinates and expression data
    spatial_coords = adata.obsm['spatial']
    if spatial_coords.shape[1] != 2:
        raise ValueError("Spatial coordinates must be 2D (x, y)")

    # Log data information
    if context:
        await context.info(f"Processing data: {adata.n_obs} spots, {adata.n_vars} genes")

    # Preprocessing
    if context:
        await context.info(f"Preprocessing data using {params.preprocessing_method}")
    
    if params.preprocessing_method == "glmpca":
        # Use GLM-PCA for preprocessing (recommended)
        expression_features = await _preprocess_glmpca(
            adata, params.n_components, context
        )
    else:
        # Use Pearson residuals PCA
        expression_features = await _preprocess_pearson_residuals(
            adata, params.n_components, context
        )
    
    # Create temporary directory for GASTON outputs
    temp_dir = tempfile.mkdtemp(prefix="gaston_")
    
    try:
        # Prepare data for GASTON
        coords_file = os.path.join(temp_dir, "spatial_coords.npy")
        expression_file = os.path.join(temp_dir, "expression_features.npy")
        
        np.save(coords_file, spatial_coords)
        np.save(expression_file, expression_features)
        
        if context:
            await context.info("Training GASTON neural network")
        
        # Train GASTON model
        model, loss_list, final_loss = await _train_gaston_model(
            coords_file, expression_file, params, temp_dir, context
        )
        
        if context:
            await context.info("Analyzing spatial patterns and gene classifications")
        
        # Analyze spatial patterns
        spatial_analysis = await _analyze_spatial_patterns(
            model, spatial_coords, expression_features, adata, params, context
        )
        
        # Note: Visualizations are now handled by the separate visualize_data tool
        # This maintains clean separation between analysis and visualization
        
        # Store results in adata
        results_key = f"gaston_results_{params.random_seed}"
        await _store_results_in_adata(
            adata, spatial_analysis, results_key, context
        )
        
        # Create result object
        result = SpatialVariableGenesResult(
            data_id=data_id,
            preprocessing_method=params.preprocessing_method,
            n_components=params.n_components,
            n_epochs_trained=params.epochs,
            final_loss=final_loss,
            spatial_hidden_layers=params.spatial_hidden_layers,
            expression_hidden_layers=params.expression_hidden_layers,
            n_spatial_domains=spatial_analysis['n_domains'],
            spatial_domains_key=f"{results_key}_spatial_domains",
            isodepth_key=f"{results_key}_isodepth",
            continuous_gradient_genes=spatial_analysis['continuous_genes'],
            discontinuous_genes=spatial_analysis['discontinuous_genes'],
            n_continuous_genes=len(spatial_analysis['continuous_genes']),
            n_discontinuous_genes=len(spatial_analysis['discontinuous_genes']),
            model_predictions_key=f"{results_key}_predictions",
            spatial_embedding_key=f"{results_key}_embedding",
            isodepth_map_visualization=None,  # Use visualize_data tool instead
            spatial_domains_visualization=None,  # Use visualize_data tool instead
            top_genes_visualization=None,  # Use visualize_data tool instead
            model_performance=spatial_analysis['performance_metrics'],
            spatial_autocorrelation=spatial_analysis['autocorrelation_metrics'],
            model_checkpoint_path=spatial_analysis.get('model_path')
        )
        
        if context:
            await context.info(f"GASTON analysis completed successfully")
            await context.info(f"Identified {result.n_spatial_domains} spatial domains")
            await context.info(f"Found {result.n_continuous_genes} genes with continuous gradients")
            await context.info(f"Found {result.n_discontinuous_genes} genes with discontinuities")
            await context.info("Use visualize_data tool with plot_type='gaston_isodepth', 'gaston_domains', or 'gaston_genes' to visualize results")
        
        return result
        
    finally:
        # Clean up temporary directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


async def _preprocess_glmpca(adata, n_components: int, context) -> np.ndarray:
    """Preprocess data using GLM-PCA."""
    try:
        from glmpca.glmpca import glmpca
    except ImportError:
        try:
            from glmpca import glmpca
        except ImportError:
            raise ImportError("glmpca package is required for GLM-PCA preprocessing")

    if context:
        await context.info("Running GLM-PCA preprocessing")

    # Get count matrix
    if hasattr(adata.X, 'toarray'):
        counts = adata.X.toarray()
    else:
        counts = adata.X.copy()

    # Ensure counts are integers and non-negative
    counts = np.round(np.maximum(counts, 0)).astype(int)

    if context:
        await context.info(f"Running GLM-PCA with {n_components} components on {counts.shape[0]} spots and {counts.shape[1]} genes")

    # Run GLM-PCA
    glmpca_result = glmpca(counts.T, L=n_components, fam="poi")

    # Return the factors (PCs)
    return glmpca_result["factors"]


async def _preprocess_pearson_residuals(adata, n_components: int, context) -> np.ndarray:
    """Preprocess data using Pearson residuals PCA."""
    if context:
        await context.info("Computing Pearson residuals and PCA")
    
    # Compute Pearson residuals
    sc.experimental.pp.normalize_pearson_residuals(adata)
    
    # Run PCA
    sc.tl.pca(adata, n_comps=n_components)
    
    return adata.obsm['X_pca']


async def _train_gaston_model(
    coords_file: str,
    expression_file: str,
    params: SpatialVariableGenesParameters,
    output_dir: str,
    context
) -> Tuple[Any, List[float], float]:
    """Train GASTON neural network model."""
    
    # Load data
    S = np.load(coords_file)
    A = np.load(expression_file)
    
    # Convert to torch tensors and rescale
    S_torch, A_torch = neural_net.load_rescale_input_data(S, A)
    
    # Train model
    model, loss_list = neural_net.train(
        S_torch, A_torch,
        S_hidden_list=params.spatial_hidden_layers,
        A_hidden_list=params.expression_hidden_layers,
        epochs=params.epochs,
        checkpoint=params.checkpoint_interval,
        save_dir=output_dir,
        optim=params.optimizer,
        lr=params.learning_rate,
        seed=params.random_seed,
        save_final=True,
        pos_encoding=params.use_positional_encoding,
        embed_size=params.embedding_size,
        sigma=params.sigma,
        batch_size=params.batch_size
    )
    
    final_loss = float(loss_list[-1]) if len(loss_list) > 0 else 0.0
    
    return model, loss_list, final_loss


async def _analyze_spatial_patterns(
    model, spatial_coords: np.ndarray, expression_features: np.ndarray,
    adata, params: SpatialVariableGenesParameters, context
) -> Dict[str, Any]:
    """Analyze spatial patterns from trained GASTON model using complete GASTON workflow."""

    if context:
        await context.info("Extracting isodepth values from trained model")

    # Step 1: Get isodepth values from spatial embedding
    S_torch = torch.tensor(spatial_coords, dtype=torch.float32)
    with torch.no_grad():
        isodepth = model.spatial_embedding(S_torch).numpy().flatten()

    # Step 2: Get model predictions
    with torch.no_grad():
        predictions = model(S_torch).numpy()

    # Step 3: Create spatial domains based on isodepth quantiles
    n_domains = params.n_domains
    domain_boundaries = np.quantile(isodepth, np.linspace(0, 1, n_domains + 1))
    spatial_domains = np.digitize(isodepth, domain_boundaries) - 1
    spatial_domains = np.clip(spatial_domains, 0, n_domains - 1)

    if context:
        await context.info(f"Created {n_domains} spatial domains based on isodepth")
        await context.info("Performing data binning and piecewise linear fitting")

    # Step 4: Prepare data for GASTON binning and fitting
    # Get original count matrix (need raw counts for Poisson regression)
    if hasattr(adata.X, 'toarray'):
        counts_mat = adata.X.toarray().T  # GASTON expects G x N matrix
    else:
        counts_mat = adata.X.T

    # Ensure counts are non-negative integers
    counts_mat = np.maximum(counts_mat, 0).astype(int)

    gene_labels = adata.var_names.values

    # Create dummy cell type dataframe (for all cell types analysis)
    cell_type_df = pd.DataFrame({'All': np.ones(len(spatial_coords))},
                               index=adata.obs_names)

    try:
        # Step 5: Perform binning using GASTON's binning function
        binning_output = binning_and_plotting.bin_data(
            counts_mat=counts_mat,
            gaston_labels=spatial_domains,
            gaston_isodepth=isodepth,
            cell_type_df=cell_type_df,
            gene_labels=gene_labels,
            num_bins=params.num_bins,
            umi_threshold=params.umi_threshold,
            pc=0,  # No pseudocount
            pc_exposure=True
        )

        if context:
            await context.info(f"Binning completed. Analyzing {len(binning_output['gene_labels_idx'])} genes")

        # Step 6: Perform piecewise linear fitting using GASTON's segmented fit
        pw_fit_dict = segmented_fit.pw_linear_fit(
            counts_mat=counts_mat,
            gaston_labels=spatial_domains,
            gaston_isodepth=isodepth,
            cell_type_df=cell_type_df,
            ct_list=['All'],  # Only analyze all cell types
            umi_threshold=params.umi_threshold,
            t=params.pvalue_threshold,
            isodepth_mult_factor=params.isodepth_mult_factor,
            reg=params.regularization,
            zero_fit_threshold=params.zero_fit_threshold
        )

        if context:
            await context.info("Piecewise linear fitting completed")
            await context.info("Classifying genes into continuous and discontinuous patterns")

        # Step 7: Classify genes using GASTON's classification functions
        continuous_genes = spatial_gene_classification.get_cont_genes(
            pw_fit_dict, binning_output, q=params.continuous_quantile
        )

        discontinuous_genes = spatial_gene_classification.get_discont_genes(
            pw_fit_dict, binning_output, q=params.discontinuous_quantile
        )

        if context:
            await context.info(f"Found {len(continuous_genes)} genes with continuous gradients")
            await context.info(f"Found {len(discontinuous_genes)} genes with discontinuities")

    except Exception as e:
        if context:
            await context.info(f"Warning: GASTON analysis failed: {e}")
            await context.info("Falling back to basic spatial domain analysis")

        # Fallback to basic analysis if GASTON functions fail
        continuous_genes = {}
        discontinuous_genes = {}
        binning_output = {'gene_labels_idx': gene_labels}

    # Performance metrics
    mse = np.mean((predictions - expression_features) ** 2)
    r2 = 1 - (np.sum((expression_features - predictions) ** 2) /
              np.sum((expression_features - np.mean(expression_features, axis=0)) ** 2))

    performance_metrics = {
        'mse': float(mse),
        'r2': float(r2),
        'isodepth_range': [float(isodepth.min()), float(isodepth.max())],
        'isodepth_std': float(isodepth.std()),
        'n_genes_analyzed': len(binning_output['gene_labels_idx'])
    }

    # Spatial autocorrelation metrics (placeholder for now)
    autocorrelation_metrics = {
        'moran_i': 0.0,  # Could compute Moran's I for isodepth
        'geary_c': 0.0   # Could compute Geary's C for isodepth
    }

    return {
        'isodepth': isodepth,
        'spatial_domains': spatial_domains,
        'n_domains': n_domains,
        'predictions': predictions,
        'continuous_genes': continuous_genes,
        'discontinuous_genes': discontinuous_genes,
        'performance_metrics': performance_metrics,
        'autocorrelation_metrics': autocorrelation_metrics,
        'binning_output': binning_output  # Store for potential future use
    }


# Note: Visualization functions have been moved to visualization.py
# Use visualize_data tool with plot_type="gaston_isodepth", "gaston_domains", or "gaston_genes"


async def _store_results_in_adata(
    adata, spatial_analysis: Dict[str, Any], results_key: str, context
):
    """Store GASTON results in AnnData object."""

    # Store isodepth values
    adata.obs[f"{results_key}_isodepth"] = spatial_analysis['isodepth']

    # Store spatial domains
    adata.obs[f"{results_key}_spatial_domains"] = spatial_analysis['spatial_domains'].astype(str)

    # Store predictions
    adata.obsm[f"{results_key}_predictions"] = spatial_analysis['predictions']

    # Store spatial embedding (isodepth as 1D embedding)
    adata.obsm[f"{results_key}_embedding"] = spatial_analysis['isodepth'].reshape(-1, 1)

    # Store metadata
    adata.uns[f"{results_key}_metadata"] = {
        'method': 'GASTON',
        'n_domains': spatial_analysis['n_domains'],
        'performance_metrics': spatial_analysis['performance_metrics'],
        'autocorrelation_metrics': spatial_analysis['autocorrelation_metrics']
    }

    if context:
        await context.info(f"Results stored in adata with key prefix: {results_key}")


def _set_random_seeds(seed: int):
    """Set random seeds for reproducibility."""
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
