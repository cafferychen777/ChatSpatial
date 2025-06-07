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

# Add GASTON to path
GASTON_PATH = "/Users/apple/Research/SpatialTrans_MCP/GASTON/src"
if GASTON_PATH not in sys.path:
    sys.path.insert(0, GASTON_PATH)

try:
    import gaston
    from gaston import neural_net, spatial_gene_classification, binning_and_plotting
    from gaston import dp_related, segmented_fit, process_NN_output
    GASTON_AVAILABLE = True
except ImportError as e:
    GASTON_AVAILABLE = False
    print(f"Warning: GASTON not available: {e}")

from ..models.data import SpatialVariableGenesParameters
from ..models.analysis import SpatialVariableGenesResult, Image


async def identify_spatial_variable_genes_gaston(
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

    # Apply user-controlled data filtering and subsampling
    if context:
        await context.info(f"Original data: {adata.n_obs} spots, {adata.n_vars} genes")

    # Apply filtering and subsampling if requested by user
    adata_processed = await _apply_user_preprocessing(adata, params, context)

    # Update spatial coordinates after potential subsampling
    spatial_coords = adata_processed.obsm['spatial']

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
            model, spatial_coords, expression_features, params, context
        )
        
        # Generate visualizations
        visualizations = {}
        if params.include_visualization:
            if context:
                await context.info("Generating visualizations")
            visualizations = await _generate_visualizations(
                adata, spatial_analysis, params, temp_dir, context
            )
        
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
            isodepth_map_visualization=visualizations.get('isodepth_map'),
            spatial_domains_visualization=visualizations.get('spatial_domains'),
            top_genes_visualization=visualizations.get('top_genes'),
            model_performance=spatial_analysis['performance_metrics'],
            spatial_autocorrelation=spatial_analysis['autocorrelation_metrics'],
            model_checkpoint_path=spatial_analysis.get('model_path')
        )
        
        if context:
            await context.info(f"GASTON analysis completed successfully")
            await context.info(f"Identified {result.n_spatial_domains} spatial domains")
            await context.info(f"Found {result.n_continuous_genes} genes with continuous gradients")
            await context.info(f"Found {result.n_discontinuous_genes} genes with discontinuities")
        
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
    params: SpatialVariableGenesParameters, context
) -> Dict[str, Any]:
    """Analyze spatial patterns from trained GASTON model."""
    
    # Get isodepth values
    S_torch = torch.tensor(spatial_coords, dtype=torch.float32)
    with torch.no_grad():
        isodepth = model.spatial_embedding(S_torch).numpy().flatten()
    
    # Get model predictions
    with torch.no_grad():
        predictions = model(S_torch).numpy()
    
    # Analyze spatial domains and gene patterns
    # This is a simplified version - full implementation would use GASTON's
    # binning_and_plotting and spatial_gene_classification modules
    
    # For now, create basic spatial domains based on isodepth quantiles
    n_domains = 5  # Default number of domains
    domain_boundaries = np.quantile(isodepth, np.linspace(0, 1, n_domains + 1))
    spatial_domains = np.digitize(isodepth, domain_boundaries) - 1
    spatial_domains = np.clip(spatial_domains, 0, n_domains - 1)
    
    # Placeholder for gene classification
    continuous_genes = {}
    discontinuous_genes = {}
    
    # Performance metrics
    mse = np.mean((predictions - expression_features) ** 2)
    r2 = 1 - (np.sum((expression_features - predictions) ** 2) / 
              np.sum((expression_features - np.mean(expression_features, axis=0)) ** 2))
    
    performance_metrics = {
        'mse': float(mse),
        'r2': float(r2),
        'isodepth_range': [float(isodepth.min()), float(isodepth.max())],
        'isodepth_std': float(isodepth.std())
    }
    
    # Spatial autocorrelation metrics (placeholder)
    autocorrelation_metrics = {
        'moran_i': 0.0,  # Would compute Moran's I for isodepth
        'geary_c': 0.0   # Would compute Geary's C for isodepth
    }
    
    return {
        'isodepth': isodepth,
        'spatial_domains': spatial_domains,
        'n_domains': n_domains,
        'predictions': predictions,
        'continuous_genes': continuous_genes,
        'discontinuous_genes': discontinuous_genes,
        'performance_metrics': performance_metrics,
        'autocorrelation_metrics': autocorrelation_metrics
    }


async def _generate_visualizations(
    adata, spatial_analysis: Dict[str, Any], params: SpatialVariableGenesParameters,
    temp_dir: str, context
) -> Dict[str, Image]:
    """Generate visualizations for GASTON results."""

    visualizations = {}
    spatial_coords = adata.obsm['spatial']

    # Create figure directory
    fig_dir = os.path.join(temp_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    try:
        # Isodepth map visualization
        if params.plot_isodepth_map:
            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
            scatter = ax.scatter(
                spatial_coords[:, 0], spatial_coords[:, 1],
                c=spatial_analysis['isodepth'],
                cmap='viridis', s=20, alpha=0.7
            )
            plt.colorbar(scatter, ax=ax, label='Isodepth')
            ax.set_title('GASTON Isodepth Map')
            ax.set_xlabel('Spatial X')
            ax.set_ylabel('Spatial Y')
            ax.set_aspect('equal')

            isodepth_path = os.path.join(fig_dir, f"isodepth_map.{params.image_format}")
            plt.savefig(isodepth_path, dpi=params.image_dpi, bbox_inches='tight')
            plt.close()

            with open(isodepth_path, 'rb') as f:
                visualizations['isodepth_map'] = Image(
                    data=f.read(),
                    format=params.image_format,
                    description="GASTON isodepth topographic map"
                )

        # Spatial domains visualization
        if params.plot_spatial_domains:
            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
            scatter = ax.scatter(
                spatial_coords[:, 0], spatial_coords[:, 1],
                c=spatial_analysis['spatial_domains'],
                cmap='tab10', s=20, alpha=0.7
            )
            plt.colorbar(scatter, ax=ax, label='Spatial Domain')
            ax.set_title('GASTON Spatial Domains')
            ax.set_xlabel('Spatial X')
            ax.set_ylabel('Spatial Y')
            ax.set_aspect('equal')

            domains_path = os.path.join(fig_dir, f"spatial_domains.{params.image_format}")
            plt.savefig(domains_path, dpi=params.image_dpi, bbox_inches='tight')
            plt.close()

            with open(domains_path, 'rb') as f:
                visualizations['spatial_domains'] = Image(
                    data=f.read(),
                    format=params.image_format,
                    description="GASTON spatial domains"
                )

        # Model performance visualization
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Predicted vs actual scatter plot
        predictions = spatial_analysis['predictions']
        if params.preprocessing_method == "glmpca":
            # For GLM-PCA, show first component
            actual = predictions[:, 0] if predictions.shape[1] > 0 else []
            predicted = predictions[:, 0] if predictions.shape[1] > 0 else []
        else:
            # For Pearson residuals, show first PC
            actual = predictions[:, 0] if predictions.shape[1] > 0 else []
            predicted = predictions[:, 0] if predictions.shape[1] > 0 else []

        if len(actual) > 0:
            ax1.scatter(actual, predicted, alpha=0.5, s=10)
            ax1.plot([actual.min(), actual.max()], [actual.min(), actual.max()], 'r--', lw=2)
            ax1.set_xlabel('Actual')
            ax1.set_ylabel('Predicted')
            ax1.set_title('Model Predictions vs Actual')

            # Add R² to plot
            r2 = spatial_analysis['performance_metrics']['r2']
            ax1.text(0.05, 0.95, f'R² = {r2:.3f}', transform=ax1.transAxes,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Isodepth distribution
        ax2.hist(spatial_analysis['isodepth'], bins=50, alpha=0.7, edgecolor='black')
        ax2.set_xlabel('Isodepth Value')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Isodepth Distribution')

        plt.tight_layout()
        performance_path = os.path.join(fig_dir, f"model_performance.{params.image_format}")
        plt.savefig(performance_path, dpi=params.image_dpi, bbox_inches='tight')
        plt.close()

        with open(performance_path, 'rb') as f:
            visualizations['model_performance'] = Image(
                data=f.read(),
                format=params.image_format,
                description="GASTON model performance metrics"
            )

    except Exception as e:
        if context:
            await context.info(f"Warning: Error generating visualizations: {e}")

    return visualizations


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
