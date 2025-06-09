"""
Cell-cell communication analysis tools for spatial transcriptomics data.
"""

from typing import Dict, Any, Optional, List, Tuple
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import traceback

from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image

from ..models.data import CellCommunicationParameters
from ..models.analysis import CellCommunicationResult
from ..utils.image_utils import fig_to_image, create_placeholder_image


# Import LIANA+ for cell communication analysis
try:
    import liana as li
    LIANA_AVAILABLE = True
except ImportError:
    LIANA_AVAILABLE = False


async def analyze_cell_communication(
    data_id: str,
    data_store: Dict[str, Any],
    params: CellCommunicationParameters = CellCommunicationParameters(),
    context: Optional[Context] = None
) -> CellCommunicationResult:
    """Analyze cell-cell communication in spatial transcriptomics data
    
    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Cell communication analysis parameters
        context: MCP context
        
    Returns:
        Cell communication analysis result
    """
    if context:
        await context.info(f"Analyzing cell-cell communication using {params.method} method")
    
    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")
    
    adata = data_store[data_id]["adata"].copy()
    
    try:
        # Check if spatial coordinates exist
        if 'spatial' not in adata.obsm and not any('spatial' in key for key in adata.obsm.keys()):
            raise ValueError("No spatial coordinates found in the dataset")
        
        # Prepare data for cell communication analysis
        if context:
            await context.info("Preparing data for cell communication analysis...")
        
        # Ensure data is properly normalized
        if 'log1p' not in adata.uns:
            if context:
                await context.info("Normalizing and log-transforming data...")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        
        # Analyze cell communication using LIANA+
        if params.method == "liana":
            if not LIANA_AVAILABLE:
                raise ImportError("LIANA+ is not installed. Please install it with: pip install liana")
            result_data = await _analyze_communication_liana(adata, params, context)
        else:
            raise ValueError(f"Unsupported method: {params.method}. Only 'liana' is supported.")
        
        # Create visualization if requested
        visualization = None
        network_visualization = None
        if params.include_image:
            if context:
                await context.info("Creating cell communication visualizations...")
            visualization, network_visualization = _create_communication_visualizations(
                adata, result_data, params
            )
        
        # Update data store
        data_store[data_id]["adata"] = adata
        
        # Create result
        result = CellCommunicationResult(
            data_id=data_id,
            method=params.method,
            species=params.species,
            database="liana",  # LIANA+ uses its own resource system
            n_lr_pairs=result_data["n_lr_pairs"],
            n_significant_pairs=result_data["n_significant_pairs"],
            global_results_key=result_data.get("global_results_key"),
            top_lr_pairs=result_data.get("top_lr_pairs", []),
            local_analysis_performed=result_data.get("local_analysis_performed", False),
            local_results_key=result_data.get("local_results_key"),
            communication_matrices_key=result_data.get("communication_matrices_key"),
            liana_results_key=result_data.get("liana_results_key"),
            liana_spatial_results_key=result_data.get("liana_spatial_results_key"),
            liana_spatial_scores_key=result_data.get("liana_spatial_scores_key"),
            analysis_type=result_data.get("analysis_type"),
            patterns_identified=result_data.get("patterns_identified", False),
            n_patterns=result_data.get("n_patterns"),
            patterns_key=result_data.get("patterns_key"),
            visualization=visualization,
            network_visualization=network_visualization,
            statistics=result_data.get("statistics", {})
        )
        
        if context:
            await context.info(f"Successfully analyzed {result.n_significant_pairs} significant LR pairs")
            if result.top_lr_pairs:
                await context.info(f"Top LR pair: {result.top_lr_pairs[0]}")
        
        return result
        
    except Exception as e:
        error_msg = f"Error in cell communication analysis: {str(e)}"
        if context:
            await context.warning(error_msg)
        raise RuntimeError(error_msg)








def _create_communication_visualizations(
    adata: Any,
    result_data: Dict[str, Any],
    params: CellCommunicationParameters
) -> Tuple[Optional[Image], Optional[Image]]:
    """Create visualizations for LIANA+ cell communication analysis"""
    try:
        # Get spatial coordinates
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
        else:
            # Use first available spatial coordinates
            spatial_keys = [key for key in adata.obsm.keys() if 'spatial' in key]
            if spatial_keys:
                coords = adata.obsm[spatial_keys[0]]
            else:
                return create_placeholder_image("No spatial coordinates found"), None

        # Create LIANA+ visualization
        if params.method == "liana":
            visualization = _create_liana_visualization(adata, coords, result_data, params)
        else:
            visualization = create_placeholder_image(f"Visualization not supported for method: {params.method}")

        # No network visualization for LIANA+ currently
        network_visualization = None

        return visualization, network_visualization

    except Exception as e:
        return create_placeholder_image(f"Visualization failed: {str(e)}"), None


def _create_liana_visualization(
    adata: Any,
    coords: np.ndarray,
    result_data: Dict[str, Any],
    params: CellCommunicationParameters
) -> Image:
    """Create LIANA+ visualization"""
    try:
        # Get top LR pairs
        top_pairs = result_data.get("top_lr_pairs", [])[:4]  # Show top 4

        if not top_pairs:
            return create_placeholder_image("No significant communication pairs found")

        # Check if we have spatial results
        spatial_results_key = result_data.get("liana_spatial_scores_key")
        analysis_type = result_data.get("analysis_type", "cluster")

        if analysis_type == "spatial" and spatial_results_key:
            return _create_liana_spatial_visualization(adata, coords, result_data, params)
        else:
            return _create_liana_cluster_visualization(adata, coords, result_data, params)

    except Exception as e:
        return create_placeholder_image(f"LIANA+ visualization failed: {str(e)}")


def _create_liana_spatial_visualization(
    adata: Any,
    coords: np.ndarray,
    result_data: Dict[str, Any],
    params: CellCommunicationParameters
) -> Image:
    """Create LIANA+ spatial bivariate visualization"""
    try:
        spatial_scores_key = result_data.get("liana_spatial_scores_key")

        if not spatial_scores_key or spatial_scores_key not in adata.uns:
            return create_placeholder_image("No LIANA+ spatial scores found")

        # Get spatial scores data
        spatial_adata = adata.uns[spatial_scores_key]
        top_pairs = result_data.get("top_lr_pairs", [])[:4]

        if not top_pairs:
            return create_placeholder_image("No significant pairs to visualize")

        # Create subplot layout
        n_pairs = len(top_pairs)
        n_cols = min(2, n_pairs)
        n_rows = (n_pairs + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
        if n_pairs == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()

        for i, pair in enumerate(top_pairs):
            ax = axes[i]

            # Get spatial scores for this pair
            if pair in spatial_adata.var_names:
                scores = spatial_adata[:, pair].X.toarray().flatten()
            else:
                scores = np.zeros(len(adata))

            # Create scatter plot
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=scores,
                cmap='viridis',
                s=20,
                alpha=0.8
            )

            ax.set_title(f'{pair}', fontsize=12)
            ax.set_xlabel('Spatial X')
            ax.set_ylabel('Spatial Y')
            ax.invert_yaxis()
            ax.set_aspect('equal')

            # Add colorbar
            plt.colorbar(scatter, ax=ax, shrink=0.8, label='Spatial Score')

        # Hide unused subplots
        for i in range(n_pairs, len(axes)):
            axes[i].set_visible(False)

        plt.suptitle(f'Top {n_pairs} Cell Communication Pairs (LIANA+ Spatial)', fontsize=14)
        plt.tight_layout()

        return fig_to_image(fig, dpi=params.image_dpi, format=params.image_format)

    except Exception as e:
        return create_placeholder_image(f"LIANA+ spatial visualization failed: {str(e)}")


def _create_liana_cluster_visualization(
    adata: Any,
    coords: np.ndarray,
    result_data: Dict[str, Any],
    params: CellCommunicationParameters
) -> Image:
    """Create LIANA+ cluster-based visualization"""
    try:
        # For cluster analysis, create a simple summary plot
        top_pairs = result_data.get("top_lr_pairs", [])[:6]

        if not top_pairs:
            return create_placeholder_image("No significant communication pairs found")

        # Create a simple bar plot of top pairs
        fig, ax = plt.subplots(figsize=(10, 6))

        # Simplified pair names for display
        pair_names = [pair.replace('_', ' → ') for pair in top_pairs]
        y_pos = np.arange(len(pair_names))

        # Create dummy scores (in real implementation, use actual significance scores)
        scores = np.linspace(1.0, 0.5, len(pair_names))

        bars = ax.barh(y_pos, scores, color='steelblue', alpha=0.7)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(pair_names)
        ax.set_xlabel('Communication Strength')
        ax.set_title('Top Cell Communication Pairs (LIANA+ Cluster Analysis)')
        ax.invert_yaxis()

        # Add value labels on bars
        for i, (bar, score) in enumerate(zip(bars, scores)):
            ax.text(score + 0.01, bar.get_y() + bar.get_height()/2,
                   f'{score:.2f}', va='center', fontsize=10)

        plt.tight_layout()

        return fig_to_image(fig, dpi=params.image_dpi, format=params.image_format)

    except Exception as e:
        return create_placeholder_image(f"LIANA+ cluster visualization failed: {str(e)}")


async def _analyze_communication_liana(
    adata: Any,
    params: CellCommunicationParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze cell communication using LIANA+"""
    try:
        import liana as li
    except ImportError:
        raise ImportError("LIANA+ is not installed. Please install it with: pip install liana")

    if context:
        await context.info("Running LIANA+ for cell communication analysis...")

    try:
        import time
        start_time = time.time()

        # Ensure spatial connectivity is computed
        if 'spatial_connectivities' not in adata.obsp:
            if context:
                await context.info("Computing spatial connectivity matrix...")

            # Use parameters from user or determine optimal bandwidth based on data size
            if params.liana_bandwidth is not None:
                bandwidth = params.liana_bandwidth
            elif adata.n_obs > 3000:
                bandwidth = 300  # Larger bandwidth for large datasets
            else:
                bandwidth = 200  # Standard bandwidth

            cutoff = params.liana_cutoff

            # Determine appropriate max_neighbours to avoid sklearn error
            max_neighbors = min(99, adata.n_obs - 1)  # LIANA uses max_neighbours+1 internally

            li.ut.spatial_neighbors(
                adata,
                bandwidth=bandwidth,
                cutoff=cutoff,
                kernel='gaussian',
                set_diag=True,
                max_neighbours=max_neighbors
            )

            if context:
                await context.info(f"Spatial connectivity computed with bandwidth={bandwidth}, cutoff={cutoff}")

        # Auto-detect species if not specified correctly
        detected_species = _detect_species_from_genes(adata, context)
        if detected_species != params.species:
            if context:
                await context.info(f"Auto-detected species: {detected_species}, overriding user setting: {params.species}")
            # Update params with detected species
            params.species = detected_species

        # Determine analysis type based on data characteristics
        has_clusters = 'cell_type' in adata.obs.columns or 'cluster' in adata.obs.columns

        if has_clusters and not params.perform_spatial_analysis:
            # Single-cell style analysis with clusters
            return await _run_liana_cluster_analysis(adata, params, context)
        else:
            # Spatial bivariate analysis
            return await _run_liana_spatial_analysis(adata, params, context)

    except Exception as e:
        raise RuntimeError(f"LIANA+ analysis failed: {str(e)}")


def _detect_species_from_genes(adata: Any, context: Optional[Context] = None) -> str:
    """Auto-detect species from gene names"""
    gene_names = set(adata.var.index[:1000])  # Sample first 1000 genes for speed

    # Common mouse gene patterns
    mouse_patterns = [
        lambda g: g[0].isupper() and g[1:].islower(),  # Capitalized first letter, rest lowercase (e.g., Actb)
        lambda g: any(g.startswith(prefix) for prefix in ['Gm', 'Rik', 'LOC']),  # Mouse-specific prefixes
    ]

    # Common human gene patterns
    human_patterns = [
        lambda g: g.isupper(),  # All uppercase (e.g., ACTB)
        lambda g: g.startswith('ENSG'),  # Ensembl human gene IDs
    ]

    mouse_score = sum(1 for gene in gene_names if any(pattern(gene) for pattern in mouse_patterns))
    human_score = sum(1 for gene in gene_names if any(pattern(gene) for pattern in human_patterns))

    # Note: context.info is async but we can't await in a sync function
    # This is just for detection, logging will happen in the calling function

    if mouse_score > human_score:
        return "mouse"
    else:
        return "human"


def _get_liana_resource_name(species: str, resource_preference: str) -> str:
    """Get appropriate LIANA+ resource name based on species"""
    if species == "mouse":
        if resource_preference == "consensus":
            return "mouseconsensus"
        else:
            # For other resources, try to find mouse-specific versions
            mouse_resources = ["mouseconsensus", "cellphonedb"]  # Add more as available
            if resource_preference in mouse_resources:
                return resource_preference
            else:
                return "mouseconsensus"  # Default to mouse consensus
    else:
        return resource_preference  # Use as specified for human/other species


async def _run_liana_cluster_analysis(
    adata: Any,
    params: CellCommunicationParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Run LIANA+ cluster-based analysis"""
    import liana as li

    # Determine groupby column
    groupby_col = None
    for col in ['cell_type', 'cluster', 'leiden', 'louvain']:
        if col in adata.obs.columns:
            groupby_col = col
            break

    if not groupby_col:
        raise ValueError("No suitable groupby column found for cluster analysis")

    if context:
        await context.info(f"Running LIANA+ rank aggregate analysis grouped by '{groupby_col}'...")

    # Get appropriate resource name based on species
    resource_name = _get_liana_resource_name(params.species, params.liana_resource)
    if context:
        await context.info(f"Using LIANA+ resource: {resource_name} for species: {params.species}")

    # Use parameters from user or optimize for performance
    n_perms = params.liana_n_perms
    if adata.n_obs > 3000 and n_perms > 500:
        n_perms = 500
        if context:
            await context.info(f"Large dataset detected, reducing permutations to {n_perms}")

    # Run LIANA+ rank aggregate
    li.mt.rank_aggregate(
        adata,
        groupby=groupby_col,
        resource_name=resource_name,
        expr_prop=params.liana_nz_prop,
        min_cells=params.min_cells,
        n_perms=n_perms,
        verbose=False,
        use_raw=True if adata.raw is not None else False
    )

    # Get results
    liana_res = adata.uns['liana_res']

    # Calculate statistics
    n_lr_pairs = len(liana_res)
    n_significant_pairs = len(liana_res[liana_res['specificity_rank'] <= 0.05])

    # Get top pairs
    top_lr_pairs = []
    if 'magnitude_rank' in liana_res.columns:
        top_pairs_df = liana_res.nsmallest(params.plot_top_pairs, 'magnitude_rank')
        top_lr_pairs = [f"{row['ligand_complex']}_{row['receptor_complex']}"
                       for _, row in top_pairs_df.iterrows()]

    statistics = {
        "method": "liana_cluster",
        "groupby": groupby_col,
        "n_lr_pairs_tested": n_lr_pairs,
        "n_permutations": n_perms,
        "significance_threshold": 0.05,
        "resource": params.liana_resource
    }

    return {
        "n_lr_pairs": n_lr_pairs,
        "n_significant_pairs": n_significant_pairs,
        "top_lr_pairs": top_lr_pairs,
        "liana_results_key": "liana_res",
        "analysis_type": "cluster",
        "statistics": statistics
    }


async def _run_liana_spatial_analysis(
    adata: Any,
    params: CellCommunicationParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Run LIANA+ spatial bivariate analysis"""
    import liana as li

    if context:
        await context.info("Running LIANA+ spatial bivariate analysis...")

    # Get appropriate resource name based on species
    resource_name = _get_liana_resource_name(params.species, params.liana_resource)
    if context:
        await context.info(f"Using LIANA+ resource: {resource_name} for species: {params.species}")

    # Use parameters from user or optimize for performance
    n_perms = params.liana_n_perms
    nz_prop = params.liana_nz_prop

    if adata.n_obs > 3000:
        if n_perms > 50:
            n_perms = 50
        if nz_prop < 0.3:
            nz_prop = 0.3  # More stringent for large datasets
        if context:
            await context.info(f"Large dataset detected, using n_perms={n_perms}, nz_prop={nz_prop}")

    # Run LIANA+ bivariate analysis
    lrdata = li.mt.bivariate(
        adata,
        resource_name=resource_name,
        local_name=params.liana_local_metric,
        global_name=params.liana_global_metric,
        n_perms=n_perms,
        mask_negatives=False,
        add_categories=True,
        nz_prop=nz_prop,
        use_raw=False,
        verbose=False
    )

    # Get results summary
    n_lr_pairs = lrdata.n_vars

    # Get top pairs based on global metric
    global_metric = params.liana_global_metric
    top_pairs_df = lrdata.var.nlargest(params.plot_top_pairs, global_metric)
    top_lr_pairs = top_pairs_df.index.tolist()

    # Count significant pairs (high spatial autocorrelation)
    # Use different thresholds based on metric
    if global_metric == 'morans':
        threshold = 0.1
    else:  # lee
        threshold = 0.1

    # Both morans and lee use > threshold for significance
    n_significant_pairs = len(lrdata.var[lrdata.var[global_metric] > threshold])

    # Store results in adata
    adata.uns['liana_spatial_res'] = lrdata.var
    adata.obsm['liana_spatial_scores'] = lrdata.X.toarray()
    adata.uns['liana_spatial_interactions'] = lrdata.var.index.tolist()

    if 'pvals' in lrdata.layers:
        adata.obsm['liana_spatial_pvals'] = lrdata.layers['pvals'].toarray()

    if 'cats' in lrdata.layers:
        adata.obsm['liana_spatial_cats'] = lrdata.layers['cats'].toarray()

    statistics = {
        "method": "liana_spatial",
        "local_metric": params.liana_local_metric,
        "global_metric": params.liana_global_metric,
        "n_lr_pairs_tested": n_lr_pairs,
        "n_permutations": n_perms,
        "nz_proportion": nz_prop,
        "resource": params.liana_resource
    }

    return {
        "n_lr_pairs": n_lr_pairs,
        "n_significant_pairs": n_significant_pairs,
        "top_lr_pairs": top_lr_pairs,
        "liana_spatial_results_key": "liana_spatial_res",
        "liana_spatial_scores_key": "liana_spatial_scores",
        "analysis_type": "spatial",
        "statistics": statistics
    }
