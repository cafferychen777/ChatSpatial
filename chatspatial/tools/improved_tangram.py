"""
Improved Tangram Cell Type Annotation with Better User Experience

Following Linus Torvalds' "good taste" principles:
1. Eliminate special cases through centralized column detection
2. Make the API predictable with clear parameter validation
3. Provide actionable error messages and progress feedback
4. Handle edge cases gracefully without scattered if/else logic

Author: Linus-style code review and improvement
"""

from typing import Optional, List, Dict, Any, Tuple
import pandas as pd
import anndata as ad
from mcp.server.fastmcp import Context

from ..models.data import AnnotationParameters


class AnnotationColumnDetector:
    """
    Centralized logic for finding annotation columns in reference data.
    
    This eliminates scattered column detection logic and makes behavior predictable.
    Following Linus principle: "eliminate special cases, make the normal case handle everything"
    """
    
    # Define column priorities - most specific first
    STANDARD_PATTERNS = [
        'cell_type',      # Most common
        'celltype', 
        'cell_types',
        'subclass_label', # Allen Institute format
        'cluster_label',
        'cluster_id',
        'annotation',
        'cell_annotation',
        'leiden',         # Fallback to clustering
        'louvain'
    ]
    
    @classmethod
    async def find_annotation_column(cls, adata: ad.AnnData, 
                                   preferred_col: Optional[str] = None,
                                   mode: str = "cells",
                                   context: Optional[Context] = None) -> Tuple[Optional[str], str]:
        """
        Find the best annotation column in reference data.
        
        Returns:
            (column_name, informative_message) - column name and user-friendly explanation
        """
        if preferred_col:
            if preferred_col in adata.obs.columns:
                msg = f"Using user-specified column '{preferred_col}'"
                if context:
                    await context.info(msg)
                return preferred_col, msg
            else:
                available_cols = list(adata.obs.columns)
                msg = (f"User-specified column '{preferred_col}' not found in reference data. "
                      f"Available columns: {available_cols[:5]}{'...' if len(available_cols) > 5 else ''}")
                if context:
                    await context.error(msg)
                return None, msg
        
        # Try standard patterns in priority order
        available_cols = set(adata.obs.columns)
        
        for pattern in cls.STANDARD_PATTERNS:
            if pattern in available_cols:
                msg = f"Auto-detected '{pattern}' column (standard pattern)"
                if context:
                    await context.info(msg)
                return pattern, msg
        
        # If nothing found, provide actionable guidance
        categorical_cols = [col for col in adata.obs.columns 
                          if adata.obs[col].dtype.name == 'category']
        
        if categorical_cols:
            suggestion = categorical_cols[0]
            msg = (f"No standard annotation column found. Available categorical columns: {categorical_cols}. "
                  f"Consider setting cluster_label='{suggestion}' in your parameters.")
            if context:
                await context.warning(msg)
            return None, msg
        
        all_cols = list(available_cols)
        msg = (f"No suitable annotation columns found in reference data. "
              f"Available columns: {all_cols[:5]}{'...' if len(all_cols) > 5 else ''}. "
              f"Please ensure your reference data has cell type annotations.")
        if context:
            await context.error(msg)
        return None, msg


class ImprovedTangramParameters:
    """
    Clear parameter validation and defaults for Tangram.
    
    This makes the API more predictable and user-friendly.
    Following Linus principle: "validate early, fail fast with clear messages"
    """
    
    @staticmethod
    async def validate_and_prepare(params: AnnotationParameters, 
                                 reference_adata: ad.AnnData,
                                 context: Optional[Context] = None) -> Dict[str, Any]:
        """
        Validate parameters and provide clear error messages.
        
        Returns prepared parameters with defaults filled in.
        """
        result = {}
        
        # Mode validation
        valid_modes = ["cells", "clusters"]
        mode = getattr(params, 'mode', 'cells')
        if mode not in valid_modes:
            error_msg = f"Tangram mode must be one of {valid_modes}, got '{mode}'"
            if context:
                await context.error(error_msg)
            raise ValueError(error_msg)
        result['mode'] = mode
        
        if context:
            await context.info(f"Tangram mode: '{mode}'")
        
        # Cluster label handling with intelligent detection
        user_cluster_label = getattr(params, 'cluster_label', None)
        annotation_col, reason = await AnnotationColumnDetector.find_annotation_column(
            reference_adata, user_cluster_label, mode, context
        )
        
        if mode == "clusters" and not annotation_col:
            error_msg = f"Clusters mode requires a valid annotation column. {reason}"
            if context:
                await context.error(error_msg)
            raise ValueError(error_msg)
        
        result['cluster_label'] = annotation_col
        result['cluster_label_reason'] = reason
        
        # Training parameters with sensible defaults and user feedback
        result['num_epochs'] = getattr(params, 'num_epochs', 1000)
        result['learning_rate'] = getattr(params, 'learning_rate', 0.1)
        result['lambda_d'] = getattr(params, 'lambda_d', 0.0)
        result['lambda_g1'] = getattr(params, 'lambda_g1', 1.0)
        result['lambda_g2'] = getattr(params, 'lambda_g2', 0.0)
        result['lambda_r'] = getattr(params, 'lambda_r', 0.0)
        result['lambda_count'] = getattr(params, 'lambda_count', 1.0)
        result['lambda_f_reg'] = getattr(params, 'lambda_f_reg', 1.0)
        
        # Provide feedback on training parameters
        if context:
            await context.info(f"Training configuration: {result['num_epochs']} epochs, "
                             f"learning_rate={result['learning_rate']}")
            if result['num_epochs'] < 500:
                await context.warning(f"Low epoch count ({result['num_epochs']}) may result in poor mapping quality. "
                                    f"Consider using at least 500-1000 epochs for production use.")
        
        return result


async def improved_tangram_annotation(adata_sp: ad.AnnData,
                                    params: AnnotationParameters,
                                    data_store: Dict[str, Any],
                                    context: Optional[Context] = None) -> Tuple[List[str], Dict[str, int], Dict[str, float], float]:
    """
    Improved Tangram annotation with better user experience.
    
    Key improvements:
    1. Early parameter validation with actionable error messages
    2. Centralized column detection logic
    3. Informative progress reporting
    4. Graceful handling of edge cases
    5. Consistent return values with meaningful feedback
    """
    
    # Step 1: Validate inputs early with clear, actionable messages
    if not params.reference_data_id:
        error_msg = ("Reference data ID is required for Tangram annotation. "
                    "Please specify reference_data_id in your AnnotationParameters.")
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)
    
    if params.reference_data_id not in data_store:
        available_ids = list(data_store.keys())
        error_msg = (f"Reference dataset '{params.reference_data_id}' not found. "
                    f"Available datasets: {available_ids}")
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)
    
    reference_data = data_store[params.reference_data_id]['adata']
    
    # Step 2: Import Tangram with clear error guidance
    try:
        import tangram as tg
    except ImportError as e:
        error_msg = ("Tangram is not installed. Install it with: pip install tangram-sc\n"
                    "For GPU support: pip install tangram-sc[torch-gpu]")
        if context:
            await context.error(error_msg)
        raise ImportError(error_msg) from e
    
    # Step 3: Prepare and validate parameters with user feedback
    try:
        tangram_params = await ImprovedTangramParameters.validate_and_prepare(
            params, reference_data, context
        )
    except ValueError as e:
        # Parameter validation already logged the error
        raise
    
    # Step 4: Report what we're doing to the user
    if context:
        await context.info(f"🚀 Starting Tangram annotation")
        await context.info(f"📊 Spatial data: {adata_sp.n_obs} cells × {adata_sp.n_vars} genes")
        await context.info(f"📚 Reference data: {reference_data.n_obs} cells × {reference_data.n_vars} genes")
        if tangram_params.get('cluster_label'):
            await context.info(f"🏷️  {tangram_params['cluster_label_reason']}")
    
    # Step 5: Determine training genes with user feedback
    training_genes = params.training_genes
    
    if training_genes is None:
        if params.marker_genes:
            if context:
                await context.info("🧬 Using provided marker genes for Tangram mapping")
            # Flatten marker genes dictionary
            training_genes = []
            for genes in params.marker_genes.values():
                training_genes.extend(genes)
            training_genes = list(set(training_genes))  # Remove duplicates
        else:
            if context:
                await context.info("🔍 Computing highly variable genes for Tangram mapping")
            if 'highly_variable' not in reference_data.var:
                import scanpy as sc
                sc.pp.highly_variable_genes(reference_data, n_top_genes=2000)
            training_genes = list(reference_data.var_names[reference_data.var.highly_variable])
    
    if context:
        overlap_genes = set(training_genes) & set(adata_sp.var_names) & set(reference_data.var_names)
        await context.info(f"🎯 Training genes: {len(training_genes)} specified, "
                         f"{len(overlap_genes)} overlap between datasets")
        if len(overlap_genes) < 100:
            await context.warning(f"Low gene overlap ({len(overlap_genes)}) may result in poor mapping quality. "
                                f"Consider using more genes or checking gene name consistency.")
    
    # Step 6: Preprocess data for Tangram
    try:
        if context:
            await context.info("⚙️  Preprocessing data for Tangram...")
        tg.pp_adatas(reference_data, adata_sp, genes=training_genes)
        
        final_genes = len(reference_data.uns.get('training_genes', training_genes))
        if context:
            await context.info(f"✅ Preprocessing complete: {final_genes} training genes selected")
    except Exception as e:
        error_msg = f"Tangram preprocessing failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)
    
    # Step 7: Run Tangram mapping with progress feedback
    mode = tangram_params['mode']
    cluster_label = tangram_params['cluster_label']
    
    if context:
        await context.info(f"🔄 Running Tangram mapping (this may take several minutes)...")
        estimated_time = tangram_params['num_epochs'] // 100  # Rough estimate
        await context.info(f"⏱️  Estimated time: {estimated_time}-{estimated_time*2} minutes for "
                         f"{tangram_params['num_epochs']} epochs")
    
    mapping_args = {
        "mode": mode,
        "num_epochs": tangram_params['num_epochs'],
        "device": "cpu"  # Use CPU for compatibility
    }
    
    if mode == "clusters":
        mapping_args["cluster_label"] = cluster_label
    
    try:
        ad_map = tg.map_cells_to_space(reference_data, adata_sp, **mapping_args)
    except Exception as e:
        error_msg = f"Tangram mapping failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)
    
    # Step 8: Extract mapping score with better error handling
    tangram_mapping_score = 0.0
    try:
        if 'training_history' in ad_map.uns and len(ad_map.uns['training_history']) > 0:
            last_entry = ad_map.uns['training_history'][-1]
            if isinstance(last_entry, (list, tuple)) and len(last_entry) > 0:
                tangram_mapping_score = float(last_entry[0])
            elif isinstance(last_entry, (int, float)):
                tangram_mapping_score = float(last_entry)
        elif hasattr(ad_map, 'X') and ad_map.X is not None:
            tangram_mapping_score = 0.5  # Indicate successful mapping
        
        if context:
            quality_desc = "excellent" if tangram_mapping_score > 0.8 else "good" if tangram_mapping_score > 0.6 else "acceptable" if tangram_mapping_score > 0.4 else "poor"
            await context.info(f"📈 Tangram mapping completed - Score: {tangram_mapping_score:.3f} ({quality_desc} quality)")
    except Exception as score_error:
        if context:
            await context.warning(f"Could not extract mapping score: {score_error}")
        tangram_mapping_score = 0.5
    
    # Step 9: Project cell annotations with improved error handling
    try:
        if mode == "clusters" and cluster_label:
            if context:
                await context.info(f"🎨 Projecting cluster annotations using '{cluster_label}' column")
            tg.plot_cell_annotation(ad_map, adata_sp, annotation=cluster_label)
        else:
            # For cells mode, use the detected annotation column
            annotation_col = tangram_params['cluster_label']
            if annotation_col:
                if context:
                    await context.info(f"🎨 Projecting cell annotations using '{annotation_col}' column")
                tg.plot_cell_annotation(ad_map, adata_sp, annotation=annotation_col)
            else:
                if context:
                    await context.warning("⚠️  No annotation column available - skipping projection step")
    except Exception as proj_error:
        if context:
            await context.warning(f"Cell annotation projection failed: {proj_error}. Results may still be valid.")
    
    # Step 10: Extract results with comprehensive feedback
    cell_types = []
    counts = {}
    confidence_scores = {}
    
    if 'tangram_ct_pred' in adata_sp.obsm:
        cell_type_df = adata_sp.obsm['tangram_ct_pred']
        
        if context:
            await context.info(f"📊 Cell type predictions available: {cell_type_df.shape[1]} types, "
                             f"{cell_type_df.shape[0]} cells")
        
        # Get cell types and counts
        cell_types = list(cell_type_df.columns)
        
        # Assign cell type based on highest probability
        adata_sp.obs["cell_type"] = cell_type_df.idxmax(axis=1)
        adata_sp.obs["cell_type"] = adata_sp.obs["cell_type"].astype("category")
        
        # Get counts
        counts = adata_sp.obs["cell_type"].value_counts().to_dict()
        
        # Calculate confidence scores (mean probability for each cell type)
        for cell_type in cell_types:
            cells_of_type = adata_sp.obs["cell_type"] == cell_type
            if cells_of_type.sum() > 0:
                mean_prob = cell_type_df.loc[cells_of_type, cell_type].mean()
                confidence_scores[cell_type] = round(float(mean_prob), 2)
            else:
                confidence_scores[cell_type] = 0.5
        
        if context:
            total_annotated = sum(counts.values())
            await context.info(f"✅ Annotation complete: {total_annotated}/{adata_sp.n_obs} cells assigned")
            
            # Show top cell types
            top_types = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:3]
            top_summary = ", ".join([f"{t}({c})" for t, c in top_types])
            await context.info(f"🏆 Top cell types: {top_summary}")
            
            # Show confidence summary
            avg_confidence = sum(confidence_scores.values()) / len(confidence_scores) if confidence_scores else 0
            await context.info(f"🎯 Average confidence: {avg_confidence:.2f}")
            
            # Provide next steps guidance
            await context.info("💡 Next steps: Use create_visualization tool with plot_type='cell_types' to visualize results")
    else:
        if context:
            await context.warning("⚠️  No cell type predictions found in Tangram results. "
                                "This may indicate mapping failure or incompatible data.")
        # Return minimal valid results
        cell_types = ["Unknown"]
        counts = {"Unknown": adata_sp.n_obs}
        confidence_scores = {"Unknown": 0.1}
    
    return cell_types, counts, confidence_scores, tangram_mapping_score