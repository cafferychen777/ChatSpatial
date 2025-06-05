"""
Spatial variable genes identification tools for spatial transcriptomics data.
"""

from typing import Dict, Any, Optional, List, Tuple
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import traceback

from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image

from ..models.data import SpatialVariableGenesParameters
from ..models.analysis import SpatialVariableGenesResult
from ..utils.image_utils import fig_to_image, create_placeholder_image


# Import optional dependencies with fallbacks
try:
    import SpatialDE as spd
    import NaiveDE
    SPATIALDE_AVAILABLE = True
except ImportError:
    SPATIALDE_AVAILABLE = False

try:
    import SPARK
    SPARK_AVAILABLE = True
except ImportError:
    SPARK_AVAILABLE = False


async def identify_spatial_variable_genes(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters = SpatialVariableGenesParameters(),
    context: Optional[Context] = None
) -> SpatialVariableGenesResult:
    """Identify spatial variable genes in spatial transcriptomics data
    
    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Spatial variable genes identification parameters
        context: MCP context
        
    Returns:
        Spatial variable genes identification result
    """
    if context:
        await context.info(f"Identifying spatial variable genes using {params.method} method")
    
    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")
    
    adata = data_store[data_id]["adata"].copy()
    
    try:
        # Check if spatial coordinates exist
        if 'spatial' not in adata.obsm and not any('spatial' in key for key in adata.obsm.keys()):
            raise ValueError("No spatial coordinates found in the dataset")
        
        # Prepare data for spatial variable genes analysis
        if context:
            await context.info("Preparing data for spatial variable genes analysis...")
        
        # Get spatial coordinates
        if 'spatial' in adata.obsm:
            spatial_coords = adata.obsm['spatial']
        else:
            # Use first available spatial coordinates
            spatial_key = [key for key in adata.obsm.keys() if 'spatial' in key][0]
            spatial_coords = adata.obsm[spatial_key]
        
        # Create coordinate DataFrame
        coord_df = pd.DataFrame(spatial_coords, columns=['x', 'y'], index=adata.obs.index)
        
        # Identify spatial variable genes based on method
        if params.method == "spatialde":
            results, aeh_results, patterns = await _identify_svg_spatialde(
                adata, coord_df, params, context
            )
        elif params.method == "spark":
            results, aeh_results, patterns = await _identify_svg_spark(
                adata, coord_df, params, context
            )
        elif params.method == "trendsceek":
            results, aeh_results, patterns = await _identify_svg_trendsceek(
                adata, coord_df, params, context
            )
        else:
            raise ValueError(f"Unsupported method: {params.method}")
        
        # Store results in adata
        results_key = f"spatial_variable_genes_{params.method}"
        
        # Add results to adata.var
        for col in results.columns:
            if col != 'g':  # 'g' is the gene name column
                adata.var[f"{results_key}_{col}"] = results.set_index('g')[col].reindex(adata.var.index)
        
        # Get significant genes
        significant_genes = results[results['qval'] < params.significance_threshold]['g'].tolist()
        top_genes = results.nsmallest(params.n_top_genes, 'qval')['g'].tolist()
        
        # Store AEH results if performed
        aeh_patterns_key = None
        aeh_membership_key = None
        n_patterns = None
        patterns_visualization = None
        
        if params.perform_aeh and aeh_results is not None:
            aeh_patterns_key = f"{results_key}_aeh_patterns"
            aeh_membership_key = f"{results_key}_aeh_membership"
            
            # Store patterns in obsm
            adata.obsm[aeh_patterns_key] = patterns
            
            # Store membership in var
            for col in aeh_results.columns:
                if col != 'g':
                    adata.var[f"{aeh_membership_key}_{col}"] = aeh_results.set_index('g')[col].reindex(adata.var.index)
            
            n_patterns = patterns.shape[1] if patterns is not None else 0
            
            # Create patterns visualization
            if params.include_image and patterns is not None:
                patterns_visualization = _create_patterns_visualization(
                    coord_df, patterns, aeh_results, params
                )
        
        # Create visualization if requested
        visualization = None
        if params.include_image:
            if context:
                await context.info("Creating spatial variable genes visualization...")
            visualization = _create_svg_visualization(
                adata, coord_df, top_genes[:params.plot_top_genes], params
            )
        
        # Update data store
        data_store[data_id]["adata"] = adata
        
        # Create statistics
        statistics = {
            "method": params.method,
            "n_tested_genes": len(results),
            "n_significant_genes": len(significant_genes),
            "significance_threshold": params.significance_threshold,
            "mean_qval_significant": float(results[results['qval'] < params.significance_threshold]['qval'].mean()) if significant_genes else None,
            "top_gene": top_genes[0] if top_genes else None,
            "min_qval": float(results['qval'].min()) if len(results) > 0 else None
        }
        
        # Create result
        result = SpatialVariableGenesResult(
            data_id=data_id,
            method=params.method,
            n_significant_genes=len(significant_genes),
            n_tested_genes=len(results),
            significance_threshold=params.significance_threshold,
            top_genes=top_genes,
            results_key=results_key,
            spatialde_results=results.to_dict('records') if params.method == "spatialde" else None,
            aeh_performed=params.perform_aeh and aeh_results is not None,
            aeh_patterns_key=aeh_patterns_key,
            aeh_membership_key=aeh_membership_key,
            n_patterns=n_patterns,
            visualization=visualization,
            patterns_visualization=patterns_visualization,
            statistics=statistics
        )
        
        if context:
            await context.info(f"Successfully identified {len(significant_genes)} significant spatial variable genes")
            if top_genes:
                await context.info(f"Top spatial variable gene: {top_genes[0]}")
        
        return result
        
    except Exception as e:
        error_msg = f"Error in spatial variable genes identification: {str(e)}"
        if context:
            await context.warning(error_msg)
        raise RuntimeError(error_msg)


async def _identify_svg_spatialde(
    adata: Any,
    coord_df: pd.DataFrame,
    params: SpatialVariableGenesParameters,
    context: Optional[Context] = None
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[np.ndarray]]:
    """Identify spatial variable genes using SpatialDE"""
    if not SPATIALDE_AVAILABLE:
        raise ImportError("SpatialDE is not installed. Please install it with: pip install spatialde")
    
    if context:
        await context.info("Running SpatialDE for spatial variable genes identification...")
    
    try:
        # Prepare expression data
        if context:
            await context.info("Preparing expression data...")
        
        # Get expression matrix (spots x genes)
        expr_df = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
            index=adata.obs.index,
            columns=adata.var.index
        )
        
        # Normalize and stabilize if requested
        if params.spatialde_normalize:
            if context:
                await context.info("Normalizing and stabilizing expression data...")
            
            # Stabilize variance (Anscombe transformation)
            norm_expr = NaiveDE.stabilize(expr_df.T).T
            
            # Regress out total counts if requested
            if params.spatialde_regress_out_total_counts:
                # Add total counts to coordinate dataframe
                coord_df_with_counts = coord_df.copy()
                coord_df_with_counts['total_counts'] = expr_df.sum(axis=1)
                
                # Regress out total counts
                resid_expr = NaiveDE.regress_out(
                    coord_df_with_counts, 
                    norm_expr.T, 
                    'np.log(total_counts)'
                ).T
            else:
                resid_expr = norm_expr
        else:
            resid_expr = expr_df
        
        # Run SpatialDE
        if context:
            await context.info(f"Running SpatialDE on {resid_expr.shape[1]} genes...")
        
        # Limit number of genes for performance if needed
        if resid_expr.shape[1] > 5000:
            if context:
                await context.info(f"Large number of genes detected ({resid_expr.shape[1]}), using top 5000 most variable genes")
            
            # Calculate gene variance and select top genes
            gene_var = resid_expr.var(axis=0)
            top_var_genes = gene_var.nlargest(5000).index
            resid_expr = resid_expr[top_var_genes]
        
        # Run SpatialDE
        X = coord_df[['x', 'y']]
        results = spd.run(X, resid_expr)
        
        # Sort by significance
        results = results.sort_values('qval')
        
        # Perform AEH if requested
        aeh_results = None
        patterns = None
        
        if params.perform_aeh:
            if context:
                await context.info("Performing Automatic Expression Histology (AEH)...")
            
            # Filter significant genes for AEH
            sign_results = results[results['qval'] < params.significance_threshold]
            
            if len(sign_results) > 0:
                # Determine length scale
                if params.aeh_length_scale is None:
                    # Use median length scale from significant genes
                    length_scale = sign_results['l'].median()
                    if context:
                        await context.info(f"Using median length scale: {length_scale:.3f}")
                else:
                    length_scale = params.aeh_length_scale
                
                try:
                    aeh_results, patterns = spd.aeh.spatial_patterns(
                        X, resid_expr, sign_results, 
                        C=params.aeh_n_patterns, 
                        l=length_scale,
                        verbosity=0
                    )
                    
                    if context:
                        await context.info(f"AEH identified {params.aeh_n_patterns} spatial patterns")
                        
                except Exception as e:
                    if context:
                        await context.warning(f"AEH failed: {str(e)}")
                    aeh_results = None
                    patterns = None
            else:
                if context:
                    await context.warning("No significant genes found for AEH")
        
        return results, aeh_results, patterns
        
    except Exception as e:
        raise RuntimeError(f"SpatialDE failed: {str(e)}")


async def _identify_svg_spark(
    adata: Any,
    coord_df: pd.DataFrame,
    params: SpatialVariableGenesParameters,
    context: Optional[Context] = None
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[np.ndarray]]:
    """Identify spatial variable genes using SPARK"""
    if not SPARK_AVAILABLE:
        raise ImportError("SPARK is not installed. Please install it with: pip install pyspark")
    
    if context:
        await context.info("SPARK method is not yet implemented")
    
    # Placeholder implementation
    # In a full implementation, you would integrate SPARK here
    raise NotImplementedError("SPARK method is not yet implemented")


async def _identify_svg_trendsceek(
    adata: Any,
    coord_df: pd.DataFrame,
    params: SpatialVariableGenesParameters,
    context: Optional[Context] = None
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[np.ndarray]]:
    """Identify spatial variable genes using Trendsceek"""
    if context:
        await context.info("Trendsceek method is not yet implemented")
    
    # Placeholder implementation
    # In a full implementation, you would integrate Trendsceek here
    raise NotImplementedError("Trendsceek method is not yet implemented")


def _create_svg_visualization(
    adata: Any,
    coord_df: pd.DataFrame,
    top_genes: List[str],
    params: SpatialVariableGenesParameters
) -> Image:
    """Create visualization of top spatial variable genes"""
    try:
        n_genes = len(top_genes)
        if n_genes == 0:
            return create_placeholder_image("No spatial variable genes to visualize")
        
        # Determine subplot layout
        n_cols = min(3, n_genes)
        n_rows = (n_genes + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
        if n_genes == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        # Plot each gene
        for i, gene in enumerate(top_genes):
            ax = axes[i]
            
            if gene in adata.var.index:
                # Get expression values
                if hasattr(adata.X, 'toarray'):
                    expr_values = adata.X[:, adata.var.index == gene].toarray().flatten()
                else:
                    expr_values = adata.X[:, adata.var.index == gene].flatten()
                
                # Create scatter plot
                scatter = ax.scatter(
                    coord_df['x'], 
                    coord_df['y'], 
                    c=expr_values,
                    cmap='viridis',
                    s=20,
                    alpha=0.8
                )
                
                ax.set_title(f'{gene}', fontsize=12)
                ax.set_xlabel('Spatial X')
                ax.set_ylabel('Spatial Y')
                ax.invert_yaxis()
                ax.set_aspect('equal')
                
                # Add colorbar
                plt.colorbar(scatter, ax=ax, shrink=0.8)
            else:
                ax.text(0.5, 0.5, f'Gene {gene}\nnot found', 
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{gene} (not found)')
        
        # Hide unused subplots
        for i in range(n_genes, len(axes)):
            axes[i].set_visible(False)
        
        plt.suptitle(f'Top {n_genes} Spatial Variable Genes', fontsize=14)
        plt.tight_layout()
        
        return fig_to_image(fig, dpi=params.image_dpi, format=params.image_format)
        
    except Exception as e:
        return create_placeholder_image(f"Visualization failed: {str(e)}")


def _create_patterns_visualization(
    coord_df: pd.DataFrame,
    patterns: np.ndarray,
    aeh_results: pd.DataFrame,
    params: SpatialVariableGenesParameters
) -> Image:
    """Create visualization of spatial patterns from AEH"""
    try:
        n_patterns = patterns.shape[1]
        
        # Create subplot layout
        n_cols = min(3, n_patterns)
        n_rows = (n_patterns + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
        if n_patterns == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        # Plot each pattern
        for i in range(n_patterns):
            ax = axes[i]
            
            # Get pattern expression
            pattern_expr = patterns[:, i]
            
            # Create scatter plot
            scatter = ax.scatter(
                coord_df['x'], 
                coord_df['y'], 
                c=pattern_expr,
                cmap='RdBu_r',
                s=20,
                alpha=0.8
            )
            
            # Count genes in this pattern
            n_genes_in_pattern = len(aeh_results[aeh_results['pattern'] == i])
            
            ax.set_title(f'Pattern {i} ({n_genes_in_pattern} genes)', fontsize=12)
            ax.set_xlabel('Spatial X')
            ax.set_ylabel('Spatial Y')
            ax.invert_yaxis()
            ax.set_aspect('equal')
            
            # Add colorbar
            plt.colorbar(scatter, ax=ax, shrink=0.8)
        
        # Hide unused subplots
        for i in range(n_patterns, len(axes)):
            axes[i].set_visible(False)
        
        plt.suptitle(f'Spatial Expression Patterns (AEH)', fontsize=14)
        plt.tight_layout()
        
        return fig_to_image(fig, dpi=params.image_dpi, format=params.image_format)
        
    except Exception as e:
        return create_placeholder_image(f"Patterns visualization failed: {str(e)}")
