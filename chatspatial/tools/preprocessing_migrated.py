"""
Example of migrated preprocessing tool with proper MCP error handling.
This shows how to convert the existing preprocessing.py to use the new error handling.
"""

from typing import Dict, Any, Optional
import numpy as np
import scanpy as sc
import squidpy as sq
import traceback
from mcp.server.fastmcp import Context

from ..models.data import AnalysisParameters
from ..models.analysis import PreprocessingResult
from ..utils.tool_error_handling import (
    mcp_tool_error_handler,
    create_success_result,
    create_error_result,
    dataset_not_found_error,
    invalid_parameter_error,
    analysis_failed_error
)

# Import constants from original file
from .preprocessing import (
    DEFAULT_TARGET_SUM,
    MAX_SCALE_VALUE,
    MERFISH_GENE_THRESHOLD,
    MIN_NEIGHBORS,
    MAX_NEIGHBORS_RATIO,
    MIN_KMEANS_CLUSTERS,
    MAX_TSNE_PCA_COMPONENTS,
    CLUSTERING_RESOLUTIONS,
    _detect_data_type,
    _get_adaptive_parameters,
    _safe_matrix_operation
)


# Method 1: Using decorator (simplest approach)
@mcp_tool_error_handler(include_traceback=True)
async def preprocess_data_with_decorator(
    data_id: str,
    data_store: Dict[str, Any],
    params: AnalysisParameters = AnalysisParameters(),
    context: Optional[Context] = None
) -> PreprocessingResult:
    """
    Preprocess spatial transcriptomics data using decorator for error handling.
    
    The decorator automatically catches any exceptions and returns them
    in the proper MCP format with isError: true.
    """
    # Validate dataset exists
    if data_id not in data_store:
        raise ValueError(f"Dataset '{data_id}' not found")
    
    # Get data
    dataset_info = data_store[data_id]
    adata = dataset_info["adata"]
    
    if context:
        await context.info(f"Starting preprocessing for dataset: {data_id}")
        await context.info(f"Dataset shape: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # The rest of the logic remains the same as original
    # Any exception will be automatically caught and formatted
    # ... (preprocessing logic) ...
    
    return PreprocessingResult(
        data_id=data_id,
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        n_hvgs=2000,
        clusters=10,
        qc_metrics={}
    )


# Method 2: Manual error handling (more control)
async def preprocess_data_manual(
    data_id: str,
    data_store: Dict[str, Any],
    params: AnalysisParameters = AnalysisParameters(),
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Preprocess spatial transcriptomics data with manual error handling.
    
    This approach gives you more control over specific error messages
    and when to include tracebacks.
    """
    # Step 1: Validate dataset exists
    if data_id not in data_store:
        return dataset_not_found_error(data_id).to_dict()
    
    # Step 2: Validate parameters
    if params.n_pcs < 1:
        return invalid_parameter_error(
            "n_pcs", 
            "positive integer", 
            params.n_pcs
        ).to_dict()
    
    if params.n_hvgs < 100:
        return invalid_parameter_error(
            "n_hvgs",
            "integer >= 100",
            params.n_hvgs
        ).to_dict()
    
    try:
        # Get data
        dataset_info = data_store[data_id]
        adata = dataset_info["adata"]
        
        if context:
            await context.info(f"Starting preprocessing for dataset: {data_id}")
            await context.info(f"Dataset shape: {adata.n_obs} cells x {adata.n_vars} genes")
        
        # Detect data type
        data_type = _detect_data_type(adata)
        adaptive_params = _get_adaptive_parameters(adata, data_type)
        
        # Quality control
        if params.filter_cells or params.filter_genes:
            if context:
                await context.info("Performing quality control...")
            
            try:
                # Calculate QC metrics
                sc.pp.calculate_qc_metrics(adata, inplace=True)
                
                # Filter cells and genes
                if params.filter_cells:
                    sc.pp.filter_cells(adata, min_genes=adaptive_params['min_genes_per_cell'])
                if params.filter_genes:
                    sc.pp.filter_genes(adata, min_cells=adaptive_params['min_cells_per_gene'])
                    
            except Exception as e:
                # Specific error for QC failure
                return analysis_failed_error(
                    "Quality control",
                    f"Failed to filter cells/genes: {str(e)}"
                ).to_dict()
        
        # Normalization
        if params.normalization != "none":
            if context:
                await context.info(f"Normalizing data using {params.normalization} method...")
            
            try:
                if params.normalization == "log":
                    sc.pp.normalize_total(adata, target_sum=DEFAULT_TARGET_SUM)
                    sc.pp.log1p(adata)
                elif params.normalization == "sct":
                    # Simplified SCTransform
                    sc.pp.normalize_total(adata, target_sum=DEFAULT_TARGET_SUM)
                    sc.pp.log1p(adata)
                else:
                    return invalid_parameter_error(
                        "normalization",
                        "one of: log, sct, none",
                        params.normalization
                    ).to_dict()
                    
            except Exception as e:
                return analysis_failed_error(
                    "Normalization",
                    str(e)
                ).to_dict()
        
        # Find highly variable genes
        if context:
            await context.info("Finding highly variable genes...")
        
        try:
            n_hvgs = min(params.n_hvgs, adata.n_vars - 1)
            sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs)
        except Exception as e:
            # Non-fatal error - continue with all genes
            if context:
                await context.warning(f"HVG selection failed: {e}. Using all genes.")
            adata.var['highly_variable'] = True
        
        # PCA
        if context:
            await context.info("Running PCA...")
        
        try:
            n_pcs = min(params.n_pcs, adata.n_vars - 1, adata.n_obs - 1)
            sc.tl.pca(adata, n_comps=n_pcs)
        except Exception as e:
            return analysis_failed_error(
                "PCA",
                f"Failed to compute principal components: {str(e)}"
            ).to_dict()
        
        # Clustering
        if context:
            await context.info("Computing neighborhood graph and clustering...")
        
        try:
            # Neighborhood graph
            n_neighbors = min(
                params.n_neighbors,
                max(MIN_NEIGHBORS, int(adata.n_obs * MAX_NEIGHBORS_RATIO))
            )
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
            
            # Clustering
            resolution = CLUSTERING_RESOLUTIONS.get(
                'large' if adata.n_obs > 500 else 'medium' if adata.n_obs > 100 else 'small',
                0.6
            )
            sc.tl.leiden(adata, resolution=resolution)
            n_clusters = len(adata.obs['leiden'].unique())
            
        except Exception as e:
            return analysis_failed_error(
                "Clustering",
                str(e)
            ).to_dict()
        
        # UMAP embedding
        if context:
            await context.info("Computing UMAP embedding...")
        
        try:
            sc.tl.umap(adata)
        except Exception as e:
            # Non-fatal - preprocessing can continue without UMAP
            if context:
                await context.warning(f"UMAP computation failed: {e}")
        
        # Spatial neighbors (if spatial data exists)
        if 'spatial' in adata.obsm:
            if context:
                await context.info("Computing spatial neighbors...")
            
            try:
                sq.gr.spatial_neighbors(adata, coord_type="generic")
            except Exception as e:
                if context:
                    await context.warning(f"Spatial neighbor computation failed: {e}")
        
        # Update dataset info
        dataset_info["preprocessing_done"] = True
        dataset_info["data_type"] = data_type
        
        # Prepare QC metrics
        qc_metrics = {}
        if 'n_genes_by_counts' in adata.obs.columns:
            qc_metrics = {
                "mean_genes_per_cell": float(adata.obs['n_genes_by_counts'].mean()),
                "mean_counts_per_cell": float(adata.obs['total_counts'].mean()) if 'total_counts' in adata.obs else 0
            }
        
        # Create result
        result = PreprocessingResult(
            data_id=data_id,
            n_cells=adata.n_obs,
            n_genes=adata.n_vars,
            n_hvgs=int(sum(adata.var.highly_variable)) if 'highly_variable' in adata.var else n_hvgs,
            clusters=n_clusters,
            qc_metrics=qc_metrics
        )
        
        # Return success result
        return create_success_result(result).to_dict()
        
    except Exception as e:
        # Catch-all for unexpected errors
        if context:
            await context.warning(f"Preprocessing failed: {str(e)}")
        
        # Return error with full traceback for debugging
        return create_error_result(e, include_traceback=True).to_dict()


# Example of a hybrid approach for complex tools
async def preprocess_data_hybrid(
    data_id: str,
    data_store: Dict[str, Any],
    params: AnalysisParameters = AnalysisParameters(),
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Hybrid approach: manual validation + decorator for main logic.
    
    This is useful when you want custom error messages for validation
    but automatic handling for the main processing logic.
    """
    # Manual validation with custom error messages
    if data_id not in data_store:
        return dataset_not_found_error(data_id).to_dict()
    
    if params.n_pcs < 1 or params.n_pcs > 1000:
        return invalid_parameter_error(
            "n_pcs",
            "integer between 1 and 1000",
            params.n_pcs
        ).to_dict()
    
    # Use inner function with decorator for main logic
    @mcp_tool_error_handler()
    async def _do_preprocessing():
        dataset_info = data_store[data_id]
        adata = dataset_info["adata"]
        
        # Main preprocessing logic...
        # Any exception here will be automatically handled
        
        return PreprocessingResult(
            data_id=data_id,
            n_cells=adata.n_obs,
            n_genes=adata.n_vars,
            n_hvgs=2000,
            clusters=10,
            qc_metrics={}
        )
    
    # Call the inner function
    return await _do_preprocessing()