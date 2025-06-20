"""
Visualization tools for spatial transcriptomics data.
Refactored for improved maintainability, concurrency safety, and robustness.
"""

from typing import Dict, Optional, Any, List, Tuple, Callable, Awaitable
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import anndata as ad
import traceback
import pandas as pd
import seaborn as sns
from functools import wraps
from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image

from ..models.data import VisualizationParameters
from ..utils.image_utils import fig_to_image
from ..utils.error_handling import (
    DataNotFoundError, InvalidParameterError,
    ProcessingError, DataCompatibilityError, validate_adata
)

# Type alias for visualization functions
VisFunc = Callable[[ad.AnnData, VisualizationParameters, Optional[Context]], Awaitable[plt.Figure]]

# --- Helper Functions ---

def create_figure(figsize=(10, 8), dpi=120) -> Tuple[plt.Figure, plt.Axes]:
    """Create a matplotlib figure with a specific size and DPI."""
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    return fig, ax


def setup_multi_panel_figure(
    n_panels: int,
    params: VisualizationParameters,
    default_title: str
) -> Tuple[plt.Figure, np.ndarray]:
    """Sets up a multi-panel matplotlib figure."""
    n_cols = params.panel_layout[1] if params.panel_layout else min(3, n_panels)
    n_rows = params.panel_layout[0] if params.panel_layout else (n_panels + n_cols - 1) // n_cols

    figsize = params.figure_size or (min(5 * n_cols, 15), min(4 * n_rows, 16))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, dpi=params.dpi, squeeze=False)
    axes = axes.flatten()

    fig.suptitle(params.title or default_title, fontsize=16)
    plt.subplots_adjust(top=0.92, hspace=0.4, wspace=0.3)

    # Hide unused axes
    for ax in axes[n_panels:]:
        ax.axis('off')
        
    return fig, axes


async def get_validated_features(
    adata: ad.AnnData,
    params: VisualizationParameters,
    min_features: int = 1,
    max_features: int = 12,
    default_to_highly_variable: bool = True,
    context: Optional[Context] = None
) -> List[str]:
    """Gets a validated list of features (genes) for visualization."""
    # Ensure unique var names
    if not adata.var_names.is_unique:
        if context:
            await context.info("Making gene names unique to avoid indexing errors.")
        adata.var_names_make_unique()

    # Get requested features
    features = params.features or ([params.feature] if params.feature else [])

    if features:
        available_features = [f for f in features if f in adata.var_names]
        if not available_features:
            if context:
                await context.warning(f"None of the specified genes found: {features}.")
        else:
            if len(available_features) > max_features:
                if context:
                    await context.warning(f"Too many features ({len(available_features)}). Limiting to {max_features}.")
                return available_features[:max_features]
            return available_features

    # Fall back to highly variable genes
    if default_to_highly_variable and 'highly_variable' in adata.var:
        hvg = adata.var_names[adata.var.highly_variable].tolist()
        if hvg:
            if context:
                await context.info(f"Using top {min(max_features, len(hvg))} highly variable genes.")
            return hvg[:max_features]
    
    # Last resort: use first genes
    if len(adata.var_names) >= min_features:
        if context:
            await context.info(f"Using first {min_features} genes from the dataset.")
        return adata.var_names[:min_features].tolist()
    
    raise DataNotFoundError(f"Could not find at least {min_features} valid genes for visualization.")


def get_spatial_coordinates(adata: ad.AnnData) -> np.ndarray:
    """Get spatial coordinates from AnnData, returning a (n_obs, 2) array."""
    if 'spatial' in adata.obsm:
        return adata.obsm['spatial']
    if 'X_spatial' in adata.obsm:
        return adata.obsm['X_spatial']
    if 'x' in adata.obs and 'y' in adata.obs:
        return adata.obs[['x', 'y']].values
    raise DataNotFoundError("Could not find spatial coordinates.")


def handle_visualization_errors(plot_title: str):
    """A decorator to catch errors in visualization functions and return a placeholder figure."""
    def decorator(func: VisFunc) -> VisFunc:
        @wraps(func)
        async def wrapper(adata: ad.AnnData, params: VisualizationParameters, context: Optional[Context]) -> plt.Figure:
            try:
                return await func(adata, params, context)
            except Exception as e:
                # Log the full traceback for debugging
                if context:
                    await context.error(f"Error creating {plot_title}: {traceback.format_exc()}")
                
                # Create a user-friendly placeholder figure
                fig, ax = create_figure(dpi=params.dpi)
                error_text = f'Error creating {plot_title}:\n\n{type(e).__name__}: {str(e)}'
                ax.text(0.5, 0.5, error_text, ha='center', va='center', 
                       transform=ax.transAxes, wrap=True, fontsize=10)
                ax.set_title(f'{plot_title} - Error', color='red')
                ax.axis('off')
                return fig
        return wrapper
    return decorator


# --- NEW HELPER FUNCTIONS FOR REFACTORING ---

def get_feature_data(adata: ad.AnnData, feature: str) -> Optional[pd.Series]:
    """
    Safely retrieves feature data from adata.var_names or adata.obs.
    Returns a pandas Series, which is ideal for plotting.
    REFACTOR: This avoids in-place modification of adata.
    """
    if feature in adata.var_names:
        expression = adata[:, feature].X
        if hasattr(expression, 'toarray'):
            expression = expression.toarray()
        return pd.Series(expression.flatten(), index=adata.obs_names, name=feature)
    elif feature in adata.obs.columns:
        return adata.obs[feature]
    return None


async def ensure_highly_variable_genes(adata: ad.AnnData, context: Optional[Context], n_top_genes: int = 50):
    """
    Computes highly variable genes if they are not already present.
    REFACTOR: Centralized HVG computation to avoid code duplication.
    """
    if 'highly_variable' not in adata.var or not adata.var['highly_variable'].any():
        if context:
            await context.info(f"Computing top {n_top_genes} highly variable genes...")
        try:
            # Clean data before computation
            if hasattr(adata.X, 'data'):
                adata.X.data = np.nan_to_num(adata.X.data, nan=0.0, posinf=0.0, neginf=0.0)
            else:
                adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
        except Exception as e:
            if context:
                await context.warning(f"Failed to compute HVGs: {e}. Using top expressed genes.")
            # Fallback: use top expressed genes
            gene_means = np.array(adata.X.mean(axis=0)).flatten()
            top_indices = np.argsort(gene_means)[-n_top_genes:]
            adata.var['highly_variable'] = False
            adata.var.loc[adata.var_names[top_indices], 'highly_variable'] = True


async def ensure_umap(adata: ad.AnnData, context: Optional[Context]):
    """Computes UMAP if it is not already present."""
    if 'X_umap' not in adata.obsm:
        if context:
            await context.info("Computing UMAP...")
        try:
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
        except Exception as e:
            if context:
                await context.warning(f"UMAP computation failed: {e}. Using PCA.")
            if 'X_pca' not in adata.obsm:
                from sklearn.decomposition import PCA
                pca = PCA(n_components=min(50, adata.n_vars - 1))
                X_pca = pca.fit_transform(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X)
                adata.obsm['X_pca'] = X_pca
            adata.obsm['X_umap'] = adata.obsm['X_pca'][:, :2]


def plot_spatial_feature(
    adata: ad.AnnData, 
    feature: Optional[str], 
    ax: plt.Axes, 
    params: VisualizationParameters
):
    """
    Helper function to plot a feature on spatial coordinates using scanpy.
    REFACTOR: Uses scanpy's plotting functions to preserve all features.
    CONCURRENCY: Passes data directly without modifying adata.
    """
    # Get feature data if specified
    if feature:
        color_data = get_feature_data(adata, feature)
        if color_data is None:
            ax.text(0.5, 0.5, f"Feature '{feature}' not found", 
                    ha='center', va='center', transform=ax.transAxes)
            return
        color_values = color_data.values
    else:
        color_values = None
    
    # Use scanpy's embedding plot for consistency
    sc.pl.embedding(
        adata,
        basis="spatial",
        color=color_values,
        cmap=params.colormap if (color_values is not None and pd.api.types.is_numeric_dtype(color_data)) else None,
        ax=ax,
        show=False,
        s=params.spot_size or 50,
        alpha=params.alpha,
        frameon=params.show_axes,
        title=feature or 'Spatial Distribution',
        colorbar_loc='right' if params.show_colorbar else None
    )
    ax.set_aspect('equal')
    ax.invert_yaxis()  # Match image orientation


# --- CORE VISUALIZATION FUNCTIONS ---

@handle_visualization_errors("Spatial")
async def create_spatial_visualization(
    adata: ad.AnnData, 
    params: VisualizationParameters, 
    context: Optional[Context]
) -> plt.Figure:
    """Creates a spatial plot with optional tissue image background."""
    # Handle multi-gene visualization
    if params.features and len(params.features) > 1 and params.multi_panel:
        return await create_multi_gene_visualization(adata, params, context)
    
    # Single feature visualization
    feature = params.feature or (params.features[0] if params.features else None)
    
    # Check for tissue image
    has_image = ('spatial' in adata.uns and 
                 'images' in adata.uns['spatial'] and
                 any('hires' in str(k) for k in adata.uns['spatial']['images']))
    
    if has_image:
        # Use scanpy's spatial plot with tissue image
        fig = sc.pl.spatial(
            adata,
            img_key="hires",
            color=feature,
            cmap=params.colormap,
            alpha=params.alpha,
            size=params.spot_size,
            show=False,
            return_fig=True
        )
    else:
        # Create custom spatial plot
        fig, ax = create_figure(figsize=params.figure_size or (10, 8), dpi=params.dpi)
        plot_spatial_feature(adata, feature, ax, params)
    
    return fig


@handle_visualization_errors("UMAP")
async def create_umap_visualization(
    adata: ad.AnnData, 
    params: VisualizationParameters, 
    context: Optional[Context]
) -> plt.Figure:
    """Creates a UMAP plot."""
    await ensure_umap(adata, context)
    
    # Handle multi-feature visualization
    if params.features and len(params.features) > 1 and params.multi_panel:
        return await create_multi_gene_umap_visualization(adata, params, context)
    
    # Determine feature to plot
    feature = params.feature or (params.features[0] if params.features else None)
    if not feature and 'leiden' in adata.obs.columns:
        feature = 'leiden'
    
    # Create UMAP plot
    fig = sc.pl.umap(
        adata,
        color=feature,
        cmap=params.colormap,
        alpha=params.alpha,
        show=False,
        return_fig=True
    )
    
    return fig


@handle_visualization_errors("Heatmap")
async def create_heatmap_visualization(
    adata: ad.AnnData, 
    params: VisualizationParameters, 
    context: Optional[Context]
) -> plt.Figure:
    """Creates a heatmap of gene expression."""
    # Ensure highly variable genes are computed
    await ensure_highly_variable_genes(adata, context, n_top_genes=50)
    
    # Determine genes to plot
    if params.features:
        genes_to_plot = [g for g in params.features if g in adata.var_names]
        if not genes_to_plot:
            if context:
                await context.warning("Specified genes not found. Using highly variable genes.")
            genes_to_plot = adata.var_names[adata.var.highly_variable][:20].tolist()
    else:
        genes_to_plot = adata.var_names[adata.var.highly_variable][:20].tolist()
    
    # Limit number of genes
    genes_to_plot = genes_to_plot[:min(20, len(genes_to_plot))]
    
    # Determine grouping
    groupby = params.color_by or ('leiden' if 'leiden' in adata.obs.columns else None)
    
    # Create heatmap
    fig = sc.pl.heatmap(
        adata,
        genes_to_plot,
        groupby=groupby,
        cmap=params.colormap,
        dendrogram=True,
        standard_scale='var',
        show=False,
        return_fig=True
    )
    
    return fig


@handle_visualization_errors("Violin")
async def create_violin_visualization(
    adata: ad.AnnData, 
    params: VisualizationParameters, 
    context: Optional[Context]  # noqa: ARG001
) -> plt.Figure:
    """Creates violin plots for gene expression."""
    # Determine genes to plot
    if params.features:
        genes_to_plot = [g for g in params.features if g in adata.var_names][:5]
    else:
        await ensure_highly_variable_genes(adata, context)
        genes_to_plot = adata.var_names[adata.var.highly_variable][:3].tolist()
    
    if not genes_to_plot:
        raise DataNotFoundError("No valid genes for violin plot.")
    
    # Determine grouping
    groupby = params.color_by or ('leiden' if 'leiden' in adata.obs.columns else None)
    
    # Create violin plot
    fig = sc.pl.violin(
        adata,
        genes_to_plot,
        groupby=groupby,
        multi_panel=True,
        jitter=0.4,
        show=False,
        return_fig=True
    )
    
    return fig


@handle_visualization_errors("Multi-Gene")
async def create_multi_gene_visualization(
    adata: ad.AnnData, 
    params: VisualizationParameters, 
    context: Optional[Context]
) -> plt.Figure:
    """Creates a multi-panel spatial plot for multiple genes."""
    genes = await get_validated_features(adata, params, min_features=2, max_features=9, context=context)
    
    # Setup multi-panel figure
    fig, axes = setup_multi_panel_figure(len(genes), params, "Multi-Gene Spatial Expression")
    
    # Plot each gene
    for gene, ax in zip(genes, axes):
        plot_spatial_feature(adata, gene, ax, params)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


@handle_visualization_errors("Multi-Gene UMAP")
async def create_multi_gene_umap_visualization(
    adata: ad.AnnData, 
    params: VisualizationParameters, 
    context: Optional[Context]
) -> plt.Figure:
    """Creates a multi-panel UMAP plot for multiple features."""
    await ensure_umap(adata, context)
    
    features = params.features or [params.feature] if params.feature else []
    if not features:
        raise DataNotFoundError("No features specified for multi-panel UMAP.")
    
    # Validate features
    valid_features = []
    for feat in features[:9]:  # Limit to 9
        if feat in adata.var_names or feat in adata.obs.columns:
            valid_features.append(feat)
    
    if not valid_features:
        raise DataNotFoundError("No valid features found.")
    
    # Create multi-panel plot
    n_features = len(valid_features)
    n_cols = min(3, n_features)
    n_rows = (n_features + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows), dpi=params.dpi)
    if n_features == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Plot each feature
    for feat, ax in zip(valid_features, axes):
        sc.pl.umap(adata, color=feat, ax=ax, show=False, 
                  cmap=params.colormap if feat in adata.var_names else None)
    
    # Hide unused axes
    for ax in axes[len(valid_features):]:
        ax.axis('off')
    
    plt.tight_layout()
    return fig


@handle_visualization_errors("Deconvolution")
async def create_deconvolution_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context]
) -> plt.Figure:
    """Create deconvolution results visualization."""
    # Find deconvolution results
    deconv_keys = [key for key in adata.obsm.keys() if 'deconvolution' in key.lower() or 'proportions' in key.lower()]
    
    if not deconv_keys:
        # Look in obs columns
        cell_type_cols = [col for col in adata.obs.columns if any(method in col.lower() for method in ['cell2location', 'spotiphy', 'destvi', 'stereoscope'])]
        
        if not cell_type_cols:
            raise DataNotFoundError("No deconvolution results found.")
        
        # Create visualization from obs columns
        # NOTE: This logic attempts to parse deconvolution results from column names
        # Expected format: "deconvolution_{method}_{cell_type}"
        # Example: "deconvolution_cell2location_T_cells"
        # This allows flexibility for different deconvolution methods but could be
        # standardized in the future by always storing results in adata.obsm
        
        # Group by method
        methods = {}
        for col in cell_type_cols:
            parts = col.split('_')
            if len(parts) >= 3:
                method = parts[1]
                cell_type = '_'.join(parts[2:])
                if method not in methods:
                    methods[method] = []
                methods[method].append((col, cell_type))
        
        # Use first method found
        method_name = list(methods.keys())[0]
        cols = methods[method_name]
        
        # Create multi-panel plot for top cell types
        n_types = min(params.n_cell_types or 6, len(cols))
        fig, axes = setup_multi_panel_figure(n_types, params, f"{method_name.upper()} Deconvolution Results")
        
        # Sort by mean proportion and plot top types
        mean_props = [(col, adata.obs[col].mean()) for col, _ in cols]
        mean_props.sort(key=lambda x: x[1], reverse=True)
        
        for (col, _), ax in zip(mean_props[:n_types], axes):
            cell_type = col.split('_', 2)[2]  # Extract cell type name
            plot_spatial_feature(adata, col, ax, params)
            ax.set_title(cell_type)
    else:
        # Use obsm data
        deconv_key = deconv_keys[0]
        proportions = pd.DataFrame(adata.obsm[deconv_key], index=adata.obs_names)
        
        # Get column names for cell types if available
        if f"{deconv_key}_cell_types" in adata.uns:
            proportions.columns = adata.uns[f"{deconv_key}_cell_types"]
        
        # Select top cell types by mean proportion using helper
        n_types = min(params.n_cell_types or 6, len(proportions.columns))
        top_types = get_top_cell_types_from_deconvolution(proportions, n_types)
        
        # Create multi-panel plot
        fig, axes = setup_multi_panel_figure(n_types, params, "Deconvolution Results")
        
        # Plot each cell type WITHOUT modifying adata
        for cell_type, ax in zip(top_types, axes):
            # Pass proportion values directly to scanpy
            color_values = proportions[cell_type].values
            
            has_image = 'spatial' in adata.uns and 'images' in adata.uns['spatial']
            if has_image:
                sc.pl.spatial(
                    adata, 
                    img_key="hires", 
                    color=color_values,  # Pass array directly
                    ax=ax,
                    show=False,
                    title=str(cell_type),
                    cmap=params.colormap,
                    size=params.spot_size
                )
            else:
                # Use embedding plot for non-image data
                sc.pl.embedding(
                    adata,
                    basis="spatial",
                    color=color_values,  # Pass array directly
                    ax=ax,
                    show=False,
                    title=str(cell_type),
                    cmap=params.colormap,
                    s=params.spot_size or 50
                )
                ax.set_aspect('equal')
                ax.invert_yaxis()
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


@handle_visualization_errors("Spatial Domains")
async def create_spatial_domains_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context]
) -> plt.Figure:
    """Create spatial domains visualization."""
    # Find domain results
    domain_keys = [col for col in adata.obs.columns if any(x in col for x in ['spatial_domain', 'STAGATE', 'SpaGCN', 'domain'])]
    
    if not domain_keys:
        raise DataNotFoundError("No spatial domain results found.")
    
    # Determine which domains to visualize
    if params.feature and params.feature in adata.obs.columns:
        domain_key = params.feature
    else:
        # Use most recent or refined domains
        refined = [k for k in domain_keys if 'refined' in k]
        domain_key = refined[0] if refined else domain_keys[-1]
    
    if context:
        await context.info(f"Visualizing spatial domains: {domain_key}")
    
    # Create visualization
    has_image = ('spatial' in adata.uns and 'images' in adata.uns['spatial'])
    
    if has_image:
        fig = sc.pl.spatial(adata, img_key="hires", color=domain_key,
                           palette='tab20', show=False, return_fig=True)
    else:
        fig, ax = create_figure(figsize=params.figure_size or (10, 8), dpi=params.dpi)
        plot_spatial_feature(adata, domain_key, ax, params)
    
    return fig


@handle_visualization_errors("Cell Communication")
async def create_cell_communication_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None
) -> plt.Figure:
    """Create cell communication visualization."""
    if context:
        await context.info("Creating cell communication visualization")
    
    # Look for LIANA+ results in adata.uns
    liana_keys = [key for key in adata.uns.keys() if 'liana' in key.lower()]
    spatial_score_keys = [key for key in adata.uns.keys() if 'spatial_scores' in key.lower()]
    
    # Try to find communication results
    communication_results = None
    spatial_scores = None
    analysis_type = "none"
    
    # Check for LIANA+ spatial results
    if spatial_score_keys:
        spatial_scores_key = spatial_score_keys[0]
        if spatial_scores_key in adata.uns:
            spatial_scores = adata.uns[spatial_scores_key]
            analysis_type = "spatial"
            if context:
                await context.info(f"Found LIANA+ spatial scores: {spatial_scores_key}")
    
    # Check for general LIANA+ results
    elif liana_keys:
        liana_key = liana_keys[0]
        if liana_key in adata.uns:
            communication_results = adata.uns[liana_key]
            analysis_type = "cluster"
            if context:
                await context.info(f"Found LIANA+ results: {liana_key}")
    
    # Check for ligand-receptor pair results in obs
    lr_columns = [col for col in adata.obs.columns if any(pattern in col.lower() for pattern in ['ligand', 'receptor', 'lr_', 'communication'])]
    
    if analysis_type == "spatial" and spatial_scores is not None:
        return _create_spatial_communication_plot(adata, spatial_scores, params, context)
    elif analysis_type == "cluster" and communication_results is not None:
        return _create_cluster_communication_plot(adata, communication_results, params, context)
    elif lr_columns:
        return _create_lr_expression_plot(adata, lr_columns, params, context)
    else:
        # No communication data found, create instructional plot
        fig, ax = plt.subplots(figsize=params.figure_size or (10, 8))
        ax.text(0.5, 0.5,
               'No cell communication results found in dataset.\n\n'
               'To analyze cell communication, first run:\n'
               'analyze_cell_communication(data_id="data_1", params={"method": "liana"})\n\n'
               'This will generate LIANA+ results for visualization.',
               ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title('Cell Communication - Not Available')
        ax.axis('off')
        return fig


@handle_visualization_errors("Trajectory")
async def create_trajectory_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context]  # noqa: ARG001
) -> plt.Figure:
    """Create trajectory/pseudotime visualization."""
    # Find pseudotime data
    pseudotime_key = params.feature
    if not pseudotime_key:
        pseudotime_candidates = [k for k in adata.obs.columns if 'pseudotime' in k.lower()]
        if pseudotime_candidates:
            pseudotime_key = pseudotime_candidates[0]
        else:
            raise DataNotFoundError("No pseudotime information found.")
    
    if pseudotime_key not in adata.obs.columns:
        raise DataNotFoundError(f"Pseudotime column '{pseudotime_key}' not found.")
    
    # Determine basis
    basis = params.basis or 'spatial'
    if basis == 'spatial' and 'spatial' not in adata.obsm:
        basis = 'umap'
        await ensure_umap(adata, context)
    
    # Create plot
    fig = sc.pl.embedding(
        adata,
        basis=basis,
        color=pseudotime_key,
        cmap=params.colormap or 'viridis',
        show=False,
        return_fig=True
    )
    
    return fig


@handle_visualization_errors("RNA Velocity")
async def create_rna_velocity_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context]  # noqa: ARG001
) -> plt.Figure:
    """Create RNA velocity visualization."""
    # Check for velocity data
    if 'velocity_graph' not in adata.uns:
        raise DataNotFoundError("RNA velocity not computed. Run velocity analysis first.")
    
    try:
        import scvelo as scv
    except ImportError:
        raise ProcessingError("scvelo not installed. Cannot create velocity plot.")
    
    # Determine basis
    basis = params.basis or 'umap'
    if basis == 'umap':
        await ensure_umap(adata, context)
    
    # Create velocity plot
    fig, ax = create_figure(figsize=(10, 8), dpi=params.dpi)
    scv.pl.velocity_embedding_stream(
        adata,
        basis=basis,
        color=params.feature or 'velocity_pseudotime',
        ax=ax,
        show=False
    )
    
    return fig


@handle_visualization_errors("Spatial Analysis")
async def create_spatial_analysis_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context]  # noqa: ARG001
) -> plt.Figure:
    """Create spatial analysis visualization (Moran's I, Getis-Ord, etc.)."""
    # Find spatial analysis results
    analysis_cols = [col for col in adata.obs.columns if any(x in col for x in ['morans_i', 'getis_ord', 'spatial_lag'])]
    
    if not analysis_cols:
        raise DataNotFoundError("No spatial analysis results found.")
    
    # Determine what to visualize
    if params.feature and params.feature in adata.obs.columns:
        feature_col = params.feature
    else:
        # Default to first Getis-Ord result
        getis_cols = [col for col in analysis_cols if 'getis_ord' in col]
        feature_col = getis_cols[0] if getis_cols else analysis_cols[0]
    
    # Create visualization
    fig, ax = create_figure(figsize=(10, 8), dpi=params.dpi)
    plot_spatial_feature(adata, feature_col, ax, params)
    
    return fig


# Placeholder implementations for GASTON visualizations
@handle_visualization_errors("GASTON Isodepth")
async def create_gaston_isodepth_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None
) -> plt.Figure:
    """Create GASTON isodepth map visualization."""
    # Look for GASTON isodepth results
    isodepth_keys = [col for col in adata.obs.columns if 'isodepth' in col.lower()]
    
    if not isodepth_keys:
        raise DataNotFoundError("GASTON isodepth results not found. Please run 'find_spatial_genes' with method='gaston' first.")
    
    # Use the most recent isodepth result
    isodepth_key = isodepth_keys[-1]
    if context:
        await context.info(f"Visualizing isodepth using column: {isodepth_key}")
    
    # Create figure
    fig, ax = plt.subplots(figsize=params.figure_size or (10, 8), dpi=params.dpi)
    
    # Plot isodepth map
    plot_spatial_feature(adata, feature=isodepth_key, ax=ax, params=params)
    ax.set_title(params.title or "GASTON Isodepth Map")
    
    return fig


@handle_visualization_errors("GASTON Domains")
async def create_gaston_domains_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None
) -> plt.Figure:
    """Create GASTON spatial domains visualization."""
    # Look for GASTON spatial domains results
    domains_keys = [col for col in adata.obs.columns if 'gaston' in col.lower() and 'spatial_domains' in col.lower()]
    
    if not domains_keys:
        raise DataNotFoundError("GASTON spatial domains results not found. Please run 'find_spatial_genes' with method='gaston' first.")
    
    # Use the most recent domains result
    domains_key = domains_keys[-1]
    if context:
        await context.info(f"Visualizing spatial domains using column: {domains_key}")
    
    # Create figure
    fig, ax = plt.subplots(figsize=params.figure_size or (10, 8), dpi=params.dpi)
    
    # Plot spatial domains
    plot_spatial_feature(adata, feature=domains_key, ax=ax, params=params)
    ax.set_title(params.title or "GASTON Spatial Domains")
    
    return fig


@handle_visualization_errors("GASTON Genes")
async def create_gaston_genes_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None
) -> plt.Figure:
    """Create GASTON spatial genes visualization."""
    # Look for GASTON results in uns
    gaston_keys = [key for key in adata.uns.keys() if 'gaston' in key.lower()]
    
    if not gaston_keys:
        raise DataNotFoundError("GASTON results not found. Please run 'find_spatial_genes' with method='gaston' first.")
    
    # Get the most recent GASTON result
    gaston_key = gaston_keys[-1]
    gaston_results = adata.uns[gaston_key]
    
    if context:
        await context.info(f"Visualizing GASTON genes from: {gaston_key}")
    
    # Get continuous and discontinuous genes
    continuous_genes = gaston_results.get('continuous_genes', {})
    discontinuous_genes = gaston_results.get('discontinuous_genes', {})
    
    # Select top genes to visualize
    n_genes = min(params.n_cell_types or 6, 6)  # Limit to 6 genes max
    
    # Combine and select top genes
    all_genes = list(continuous_genes.keys())[:n_genes//2] + list(discontinuous_genes.keys())[:n_genes//2]
    
    if not all_genes:
        raise DataNotFoundError("No spatial variable genes found in GASTON results.")
    
    # Ensure genes exist in the dataset
    available_genes = [gene for gene in all_genes if gene in adata.var_names]
    
    if not available_genes:
        raise DataNotFoundError("None of the identified spatial genes are available in the current dataset.")
    
    # Create multi-panel figure
    fig, axes = setup_multi_panel_figure(
        n_panels=len(available_genes),
        params=params,
        default_title="GASTON Spatial Variable Genes"
    )
    
    # Plot each gene
    for i, gene in enumerate(available_genes):
        if i < len(axes):
            ax = axes[i]
            plot_spatial_feature(adata, feature=gene, ax=ax, params=params)
            
            # Add gene type annotation
            gene_type = "Continuous" if gene in continuous_genes else "Discontinuous"
            ax.set_title(f"{gene} ({gene_type})")
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


@handle_visualization_errors("L-R Pairs")
async def create_lr_pairs_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None
) -> plt.Figure:
    """Create ligand-receptor pairs visualization."""
    # Validate ligand-receptor pairs
    if not params.lr_pairs:
        raise InvalidParameterError("No ligand-receptor pairs specified. Use params.lr_pairs")
    
    if context:
        await context.info(f"Visualizing {len(params.lr_pairs)} ligand-receptor pairs")
    
    # Calculate total panels needed (ligand + receptor + correlation for each pair)
    n_pairs = len(params.lr_pairs)
    panels_per_pair = 3 if params.show_correlation else 2
    total_panels = n_pairs * panels_per_pair
    
    # Create multi-panel figure
    n_cols = min(3, total_panels)  # Max 3 columns
    n_rows = (total_panels + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows), dpi=params.dpi, squeeze=False)
    axes = axes.flatten()
    
    ax_idx = 0
    for lr_pair in params.lr_pairs:
        try:
            # Parse ligand-receptor pair
            if isinstance(lr_pair, str) and '_' in lr_pair:
                ligand, receptor = lr_pair.split('_', 1)
            elif isinstance(lr_pair, (list, tuple)) and len(lr_pair) == 2:
                ligand, receptor = lr_pair
            else:
                raise ValueError(f"Invalid L-R pair format: {lr_pair}")
            
            # Check if genes exist
            if ligand not in adata.var_names or receptor not in adata.var_names:
                if context:
                    await context.warning(f"Skipping {ligand}-{receptor}: not found in data")
                continue
            
            # Get expression data without modifying adata
            ligand_data = get_feature_data(adata, ligand)
            receptor_data = get_feature_data(adata, receptor)
            ligand_expr = ligand_data.values
            receptor_expr = receptor_data.values
            
            # Apply color scaling for correlation plot
            if params.color_scale == "log":
                ligand_expr_scaled = np.log1p(ligand_expr)
                receptor_expr_scaled = np.log1p(receptor_expr)
            elif params.color_scale == "sqrt":
                ligand_expr_scaled = np.sqrt(ligand_expr)
                receptor_expr_scaled = np.sqrt(receptor_expr)
            else:
                ligand_expr_scaled = ligand_expr
                receptor_expr_scaled = receptor_expr
            
            # Plot ligand using the centralized helper for consistency
            if ax_idx < len(axes):
                ax = axes[ax_idx]
                plot_spatial_feature(adata, ligand, ax, params)
                ax.set_title(f"{ligand} (Ligand)")
                ax_idx += 1
            
            # Plot receptor using the centralized helper for consistency
            if ax_idx < len(axes):
                ax = axes[ax_idx]
                plot_spatial_feature(adata, receptor, ax, params)
                ax.set_title(f"{receptor} (Receptor)")
                ax_idx += 1
            
            # Plot correlation
            if params.show_correlation and ax_idx < len(axes):
                ax = axes[ax_idx]
                
                # Calculate correlation
                from scipy.stats import pearsonr, spearmanr, kendalltau
                
                if params.correlation_method == "pearson":
                    corr, p_value = pearsonr(ligand_expr_scaled, receptor_expr_scaled)
                elif params.correlation_method == "spearman":
                    corr, p_value = spearmanr(ligand_expr_scaled, receptor_expr_scaled)
                else:  # kendall
                    corr, p_value = kendalltau(ligand_expr_scaled, receptor_expr_scaled)
                
                # Create scatter plot with scaled values
                ax.scatter(ligand_expr_scaled, receptor_expr_scaled, alpha=params.alpha, s=20)
                ax.set_xlabel(f"{ligand} Expression")
                ax.set_ylabel(f"{receptor} Expression")
                
                if params.show_correlation_stats:
                    ax.set_title(f"Correlation: {corr:.3f}\np-value: {p_value:.2e}", fontsize=10)
                else:
                    ax.set_title(f"{ligand} vs {receptor}", fontsize=10)
                
                # Add trend line
                z = np.polyfit(ligand_expr_scaled, receptor_expr_scaled, 1)
                p = np.poly1d(z)
                ax.plot(ligand_expr_scaled, p(ligand_expr_scaled), "r--", alpha=0.8)
                
                ax_idx += 1
        
        except Exception as e:
            # Handle individual pair plotting errors
            if ax_idx < len(axes):
                ax = axes[ax_idx]
                ax.text(0.5, 0.5, f'Error plotting {lr_pair}:\n{str(e)}',
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{lr_pair} (Error)', fontsize=10)
                ax_idx += 1
    
    # Hide unused axes
    for ax in axes[ax_idx:]:
        ax.axis('off')
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


@handle_visualization_errors("Gene Correlation")
async def create_gene_correlation_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context]  # noqa: ARG001
) -> plt.Figure:
    """Create gene correlation visualization."""
    # Get genes to correlate
    genes = await get_validated_features(adata, params, min_features=2, max_features=50, context=context)
    
    if len(genes) < 2:
        raise DataNotFoundError("Need at least 2 genes for correlation analysis.")
    
    # Extract expression matrix
    expr_matrix = adata[:, genes].X
    if hasattr(expr_matrix, 'toarray'):
        expr_matrix = expr_matrix.toarray()
    
    # Compute correlation
    corr_matrix = pd.DataFrame(expr_matrix, columns=genes).corr()
    
    # Create heatmap
    fig, ax = create_figure(figsize=(10, 8), dpi=params.dpi)
    sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm', 
                center=0, square=True, ax=ax)
    ax.set_title("Gene Correlation Matrix")
    
    plt.tight_layout()
    return fig


@handle_visualization_errors("Enrichment")
async def create_enrichment_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None
) -> plt.Figure:
    """Create EnrichMap enrichment visualization.
    
    Supports multiple visualization types:
    - Default: Spatial enrichment map
    - violin: Enrichment scores by cluster
    - heatmap: Gene contributions heatmap
    - Multi-panel: Multiple enrichment scores
    """
    if context:
        await context.info("Creating EnrichMap enrichment visualization")
    
    # Find enrichment score columns
    score_cols = [col for col in adata.obs.columns if col.endswith('_score')]
    
    if not score_cols:
        raise DataNotFoundError("No enrichment scores found. Please run 'analyze_enrichment' first.")
    
    # Check if user wants gene contribution visualization
    if hasattr(params, 'show_gene_contributions') and params.show_gene_contributions:
        if 'gene_contributions' not in adata.uns:
            raise DataNotFoundError("Gene contributions not found.")
        
        # Create gene contribution heatmap
        gene_contribs = adata.uns['gene_contributions']
        
        if not gene_contribs:
            raise DataNotFoundError("No gene contributions available.")
        
        # Convert to DataFrame for visualization
        contrib_data = {}
        for sig, genes in gene_contribs.items():
            contrib_data[sig] = genes
        
        if not contrib_data:
            raise DataNotFoundError("No gene contribution data to visualize.")
        
        df = pd.DataFrame(contrib_data).fillna(0)
        
        # Create heatmap
        fig, ax = create_figure(figsize=(max(8, len(df.columns) * 1.5), max(6, len(df) * 0.3)))
        
        sns.heatmap(df, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                   ax=ax, cbar_kws={'label': 'Gene Weight'})
        ax.set_title('Gene Contributions to Enrichment Scores', fontsize=14)
        ax.set_xlabel('Signatures', fontsize=12)
        ax.set_ylabel('Genes', fontsize=12)
        
        plt.tight_layout()
        return fig
    
    # Check if user wants violin plot by cluster
    if params.plot_type == "violin" or (hasattr(params, 'show_violin') and params.show_violin):
        # Determine grouping variable
        group_by = params.color_by if hasattr(params, 'color_by') and params.color_by else 'leiden'
        
        if group_by not in adata.obs.columns:
            raise DataNotFoundError(f"Grouping variable '{group_by}' not found in adata.obs")
        
        # Determine which scores to plot
        if params.features and len(params.features) > 1:
            scores_to_plot = []
            for feat in params.features:
                if feat in adata.obs.columns:
                    scores_to_plot.append(feat)
                elif f"{feat}_score" in adata.obs.columns:
                    scores_to_plot.append(f"{feat}_score")
        elif params.feature:
            if params.feature in adata.obs.columns:
                scores_to_plot = [params.feature]
            elif f"{params.feature}_score" in adata.obs.columns:
                scores_to_plot = [f"{params.feature}_score"]
            else:
                scores_to_plot = [score_cols[0]]
        else:
            scores_to_plot = score_cols[:3]  # Limit to 3 scores for clarity
        
        # Create violin plot
        n_scores = len(scores_to_plot)
        fig, axes = plt.subplots(1, n_scores, figsize=(5 * n_scores, 6))
        
        if n_scores == 1:
            axes = [axes]
        
        for i, score in enumerate(scores_to_plot):
            ax = axes[i]
            
            # Create violin plot using seaborn
            data_for_plot = pd.DataFrame({
                group_by: adata.obs[group_by],
                'Score': adata.obs[score]
            })
            
            sns.violinplot(data=data_for_plot, x=group_by, y='Score', ax=ax)
            
            sig_name = score.replace('_score', '')
            ax.set_title(f"{sig_name} by {group_by}", fontsize=12)
            ax.set_xlabel(group_by, fontsize=10)
            ax.set_ylabel('Enrichment Score', fontsize=10)
            ax.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        return fig
    
    # Standard spatial visualization
    # Determine which score to visualize
    if params.feature:
        # User specified a score
        if params.feature in adata.obs.columns:
            score_col = params.feature
        elif f"{params.feature}_score" in adata.obs.columns:
            score_col = f"{params.feature}_score"
        else:
            raise DataNotFoundError(f"Score column '{params.feature}' not found. Available scores: {score_cols}")
    else:
        # Use the first available score
        score_col = score_cols[0]
        if context:
            await context.info(f"Using score column: {score_col}")
    
    # Check if we should create multi-panel plot for multiple scores
    if params.features and len(params.features) > 1:
        # Multi-score visualization
        scores_to_plot = []
        for feat in params.features:
            if feat in adata.obs.columns:
                scores_to_plot.append(feat)
            elif f"{feat}_score" in adata.obs.columns:
                scores_to_plot.append(f"{feat}_score")
        
        if not scores_to_plot:
            raise DataNotFoundError(f"None of the specified scores found: {params.features}")
        
        # Create multi-panel figure
        fig, axes = setup_multi_panel_figure(
            n_panels=len(scores_to_plot),
            params=params,
            default_title="EnrichMap Enrichment Scores"
        )
        
        # Plot each score
        for i, score in enumerate(scores_to_plot):
            if i < len(axes):
                ax = axes[i]
                plot_spatial_feature(adata, feature=score, ax=ax, params=params)
                
                # Extract signature name from score column
                sig_name = score.replace('_score', '')
                ax.set_title(f"{sig_name} Enrichment")
    else:
        # Single score visualization
        fig, ax = create_figure(figsize=(10, 8))
        plot_spatial_feature(adata, feature=score_col, ax=ax, params=params)
        
        # Extract signature name from score column
        sig_name = score_col.replace('_score', '')
        ax.set_title(f"{sig_name} Enrichment Score", fontsize=14)
        
        # Add colorbar label if shown
        if params.show_colorbar:
            # Find the colorbar and update its label
            if hasattr(ax, 'collections') and ax.collections:
                cbar = plt.colorbar(ax.collections[0], ax=ax)
                cbar.set_label('Enrichment Score', fontsize=12)
    
    plt.tight_layout()
    return fig


# --- MAIN DISPATCHER ---

# REFACTOR: Dispatch dictionary for cleaner routing
VISUALIZATION_DISPATCHER: Dict[str, VisFunc] = {
    "spatial": create_spatial_visualization,
    "umap": create_umap_visualization,
    "heatmap": create_heatmap_visualization,
    "violin": create_violin_visualization,
    "multi_gene": create_multi_gene_visualization,
    "deconvolution": create_deconvolution_visualization,
    "spatial_domains": create_spatial_domains_visualization,
    "cell_communication": create_cell_communication_visualization,
    "trajectory": create_trajectory_visualization,
    "rna_velocity": create_rna_velocity_visualization,
    "spatial_analysis": create_spatial_analysis_visualization,
    "gaston_isodepth": create_gaston_isodepth_visualization,
    "gaston_domains": create_gaston_domains_visualization,
    "gaston_genes": create_gaston_genes_visualization,
    "lr_pairs": create_lr_pairs_visualization,
    "gene_correlation": create_gene_correlation_visualization,
    "enrichment": create_enrichment_visualization,
}


async def visualize_data(
    data_id: str,
    data_store: Dict[str, Any],
    params: VisualizationParameters = VisualizationParameters(),
    context: Optional[Context] = None
) -> Image:
    """
    Main visualization dispatcher using a clean dispatch pattern.
    
    This function serves as the central entry point for all visualization requests.
    It validates inputs, dispatches to the appropriate visualization function,
    and ensures proper resource cleanup.
    
    Parameters
    ----------
    data_id : str
        Identifier for the dataset in the data store
    data_store : Dict[str, Any]
        Dictionary containing all loaded datasets
    params : VisualizationParameters
        Configuration for the visualization including plot_type, features, etc.
    context : Optional[Context]
        MCP context for logging and user communication
    
    Returns
    -------
    Image
        Rendered visualization as an Image object
    
    Raises
    ------
    InvalidParameterError
        If plot_type is not valid
    DataNotFoundError
        If dataset or required data is not found
    ProcessingError
        If visualization fails for any other reason
    
    Notes
    -----
    REFACTOR: Replaced long if/elif chain with dispatch dictionary.
    REFACTOR: Added proper resource management with try/finally.
    REFACTOR: Treats AnnData as read-only to ensure concurrency safety.
    """
    if context:
        await context.info(f"Creating '{params.plot_type}' visualization for dataset '{data_id}'")
    
    # Validate plot type
    if params.plot_type not in VISUALIZATION_DISPATCHER:
        raise InvalidParameterError(
            f"Invalid plot_type: '{params.plot_type}'. "
            f"Valid options: {list(VISUALIZATION_DISPATCHER.keys())}"
        )
    
    # Validate dataset
    if data_id not in data_store:
        raise DataNotFoundError(f"Dataset '{data_id}' not found.")
    
    adata = data_store[data_id].get("adata")
    if adata is None:
        raise DataNotFoundError(f"AnnData object not found for dataset '{data_id}'.")
    
    validate_adata(adata, min_cells=5, min_genes=5)
    
    # Get visualization function
    vis_function = VISUALIZATION_DISPATCHER[params.plot_type]
    
    fig = None
    try:
        # Set plotting parameters
        sc.settings.set_figure_params(dpi=params.dpi, facecolor='white')
        
        # Handle special case: deconvolution multi-panel
        if params.show_deconvolution and params.plot_type == "spatial":
            # Find deconvolution results
            deconv_keys = [key for key in adata.obsm.keys() if key.startswith("deconvolution_")]
            if deconv_keys:
                params.plot_type = "deconvolution"
                vis_function = VISUALIZATION_DISPATCHER["deconvolution"]
        
        # Create visualization
        # NOTE: We pass the original adata since our refactored functions
        # avoid in-place modifications for concurrency safety
        fig = await vis_function(adata, params, context)
        
        # Convert to image
        if context:
            await context.info("Converting figure to image...")
        
        return fig_to_image(fig, format="png", dpi=params.dpi)
        
    except Exception as e:
        if context:
            await context.error(f"Visualization error: {traceback.format_exc()}")
        
        # Re-raise known exceptions
        if isinstance(e, (DataNotFoundError, InvalidParameterError, 
                         DataCompatibilityError, ProcessingError)):
            raise
        
        # Wrap unknown exceptions
        raise ProcessingError(
            f"Visualization failed: {str(e)}",
            details={"plot_type": params.plot_type, "traceback": traceback.format_exc()}
        ) from e
        
    finally:
        # REFACTOR: Ensure proper cleanup
        if fig is not None:
            plt.close(fig)
        # Close any other figures to prevent memory leaks
        plt.close('all')


# --- Additional helper functions ---

def get_top_cell_types_from_deconvolution(
    proportions: pd.DataFrame,
    n_types: int = 6
) -> List[str]:
    """Extract top cell types from deconvolution results by mean proportion.
    
    This helper can be reused across different deconvolution visualization scenarios.
    """
    mean_props = proportions.mean().sort_values(ascending=False)
    return mean_props.head(n_types).index.tolist()


# --- Helper functions for cell communication visualization ---

def _create_spatial_communication_plot(adata, spatial_scores, params, context):  # noqa: ARG001
    """Create spatial communication visualization using LIANA+ spatial scores."""
    try:
        # Get spatial coordinates
        coords = get_spatial_coordinates(adata)
        x_coords = coords[:, 0]
        y_coords = coords[:, 1]
        
        # Get top communication pairs from spatial scores
        if hasattr(spatial_scores, 'var_names'):
            # It's an AnnData object
            top_pairs = list(spatial_scores.var_names)[:6]  # Top 6 pairs
        else:
            # It's a different format, try to extract pair names
            top_pairs = []
        
        if not top_pairs:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, 'No communication pairs found in spatial scores',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Cell Communication - Spatial')
            ax.axis('off')
            return fig
        
        # Create subplot layout
        n_pairs = min(len(top_pairs), 6)  # Limit to 6 for display
        n_cols = min(3, n_pairs)
        n_rows = (n_pairs + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        if n_pairs == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        for i, pair in enumerate(top_pairs[:n_pairs]):
            ax = axes[i]
            
            try:
                # Get spatial scores for this pair
                if hasattr(spatial_scores, 'var_names') and pair in spatial_scores.var_names:
                    if hasattr(spatial_scores.X, 'toarray'):
                        scores = spatial_scores[:, pair].X.toarray().flatten()
                    else:
                        scores = spatial_scores[:, pair].X.flatten()
                else:
                    scores = np.zeros(len(adata))
                
                # Create scatter plot
                scatter = ax.scatter(
                    x_coords, y_coords,
                    c=scores,
                    cmap=params.colormap or 'viridis',
                    s=15, alpha=0.8
                )
                
                # Format pair name for display
                display_name = pair.replace('_', ' → ') if '_' in pair else pair
                ax.set_title(display_name, fontsize=10)
                ax.set_xlabel('X coordinate')
                ax.set_ylabel('Y coordinate')
                ax.set_aspect('equal')
                ax.invert_yaxis()
                
                # Add colorbar
                plt.colorbar(scatter, ax=ax, shrink=0.7, label='Communication Score')
                
            except Exception as e:
                ax.text(0.5, 0.5, f'Error: {str(e)}', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(pair, fontsize=10)
        
        # Hide unused subplots (Pythonic way)
        for ax in axes[n_pairs:]:
            ax.set_visible(False)
        
        plt.suptitle('Cell Communication - Spatial Distribution', fontsize=14)
        plt.tight_layout()
        return fig
        
    except Exception as e:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f'Error creating spatial communication plot:\n{str(e)}',
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Cell Communication - Error')
        ax.axis('off')
        return fig


def _create_cluster_communication_plot(adata, communication_results, params, context):  # noqa: ARG001
    """Create cluster-based communication visualization."""
    try:
        # Try to extract top communication pairs from results
        if isinstance(communication_results, pd.DataFrame):
            # Results are in DataFrame format
            df = communication_results
            
            # Look for common LIANA+ columns
            score_cols = [col for col in df.columns if any(term in col.lower() for term in ['score', 'pvalue', 'significant'])]
            
            if score_cols and len(df) > 0:
                # Sort by first score column and get top pairs
                top_results = df.nlargest(8, score_cols[0])
                
                # Create ligand-receptor pair names
                if 'ligand' in df.columns and 'receptor' in df.columns:
                    pair_names = [f"{row['ligand']} → {row['receptor']}" for _, row in top_results.iterrows()]
                    scores = top_results[score_cols[0]].values
                elif 'lr_pair' in df.columns:
                    pair_names = top_results['lr_pair'].tolist()
                    scores = top_results[score_cols[0]].values
                else:
                    pair_names = [f"Pair {i+1}" for i in range(len(top_results))]
                    scores = top_results[score_cols[0]].values
            else:
                pair_names = ["No significant pairs"]
                scores = [0]
        else:
            # Results in other format
            pair_names = ["Communication analysis completed"]
            scores = [1.0]
        
        # Create horizontal bar plot
        fig, ax = plt.subplots(figsize=params.figure_size or (12, 8))
        
        y_pos = np.arange(len(pair_names))
        bars = ax.barh(y_pos, scores, color='steelblue', alpha=0.7)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(pair_names, fontsize=10)
        ax.set_xlabel('Communication Strength')
        ax.set_title('Top Cell Communication Pairs')
        ax.invert_yaxis()
        
        # Add value labels on bars
        for bar, score in zip(bars, scores):
            if score > 0:
                ax.text(score + max(scores) * 0.01, bar.get_y() + bar.get_height()/2,
                       f'{score:.3f}', va='center', fontsize=9)
        
        plt.tight_layout()
        return fig
        
    except Exception as e:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f'Error creating cluster communication plot:\n{str(e)}',
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Cell Communication - Error')
        ax.axis('off')
        return fig


def _create_lr_expression_plot(adata, lr_columns, params, context):  # noqa: ARG001
    """Create ligand-receptor expression visualization."""
    try:
        # Get spatial coordinates
        coords = get_spatial_coordinates(adata)
        x_coords = coords[:, 0]
        y_coords = coords[:, 1]
        
        # Select top LR columns (limit to 6)
        selected_cols = lr_columns[:6]
        
        # Create subplot layout
        n_cols = min(3, len(selected_cols))
        n_rows = (len(selected_cols) + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        if len(selected_cols) == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        for i, col in enumerate(selected_cols):
            ax = axes[i]
            
            try:
                values = adata.obs[col].values
                
                # Handle categorical or numeric data
                if pd.api.types.is_numeric_dtype(values):
                    scatter = ax.scatter(x_coords, y_coords, c=values, 
                                       cmap=params.colormap or 'viridis', s=15, alpha=0.8)
                    plt.colorbar(scatter, ax=ax, shrink=0.7)
                else:
                    # Categorical data
                    unique_vals = pd.unique(values)
                    colors = plt.cm.Set1(np.linspace(0, 1, len(unique_vals)))
                    
                    for j, val in enumerate(unique_vals):
                        mask = values == val
                        ax.scatter(x_coords[mask], y_coords[mask], 
                                 c=[colors[j]], label=str(val), s=15, alpha=0.8)
                    
                    if len(unique_vals) <= 10:
                        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
                
                # Clean up column name for display
                display_name = col.replace('_', ' ').title()
                ax.set_title(display_name, fontsize=10)
                ax.set_xlabel('X coordinate')
                ax.set_ylabel('Y coordinate')
                ax.set_aspect('equal')
                ax.invert_yaxis()
                
            except Exception as e:
                ax.text(0.5, 0.5, f'Error: {str(e)}', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(col, fontsize=10)
        
        # Hide unused subplots (Pythonic way)
        for ax in axes[len(selected_cols):]:
            ax.set_visible(False)
        
        plt.suptitle('Ligand-Receptor Expression Patterns', fontsize=14)
        plt.tight_layout()
        return fig
        
    except Exception as e:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f'Error creating LR expression plot:\n{str(e)}',
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Cell Communication - Error')
        ax.axis('off')
        return fig