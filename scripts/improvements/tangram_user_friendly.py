#!/usr/bin/env python3
"""
Improved Tangram Annotation with Better User Experience

Following Linus Torvalds' "good taste" principles:
1. Eliminate special cases
2. Make the API predictable  
3. Provide clear feedback to users
4. Handle edge cases gracefully

Author: Linus-style code review
"""

from typing import Optional, List, Dict, Any
import pandas as pd
import anndata as ad


class AnnotationColumnDetector:
    """
    Centralized logic for finding annotation columns in reference data.
    
    This eliminates scattered column detection logic and makes behavior predictable.
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
    def find_annotation_column(cls, adata: ad.AnnData, 
                             preferred_col: Optional[str] = None,
                             mode: str = "cells") -> tuple[Optional[str], str]:
        """
        Find the best annotation column in reference data.
        
        Returns:
            (column_name, reason) - column name and explanation
        """
        if preferred_col:
            if preferred_col in adata.obs.columns:
                return preferred_col, f"User specified '{preferred_col}'"
            else:
                return None, f"User specified '{preferred_col}' not found in reference data"
        
        # Try standard patterns in order
        available_cols = set(adata.obs.columns)
        
        for pattern in cls.STANDARD_PATTERNS:
            if pattern in available_cols:
                return pattern, f"Auto-detected '{pattern}' (standard pattern)"
        
        # If nothing found, list what's available for user guidance
        categorical_cols = [col for col in adata.obs.columns 
                          if adata.obs[col].dtype.name == 'category']
        
        if categorical_cols:
            suggestion = categorical_cols[0]
            return None, f"No standard annotation column found. Available categorical columns: {categorical_cols}. Consider using '{suggestion}'"
        
        return None, f"No suitable annotation columns found in reference data. Available columns: {list(available_cols)[:5]}..."


class ImprovedTangramParameters:
    """
    Clear parameter validation and defaults for Tangram.
    
    This makes the API more predictable and user-friendly.
    """
    
    @staticmethod
    def validate_and_prepare(params: Any, reference_adata: ad.AnnData) -> Dict[str, Any]:
        """
        Validate parameters and provide clear error messages.
        
        Returns prepared parameters with defaults filled in.
        """
        result = {}
        
        # Mode validation
        valid_modes = ["cells", "clusters"]
        mode = getattr(params, 'mode', 'cells')
        if mode not in valid_modes:
            raise ValueError(f"Mode must be one of {valid_modes}, got '{mode}'")
        result['mode'] = mode
        
        # Cluster label handling
        cluster_label = getattr(params, 'cluster_label', None)
        annotation_col, reason = AnnotationColumnDetector.find_annotation_column(
            reference_adata, cluster_label, mode
        )
        
        if mode == "clusters" and not annotation_col:
            raise ValueError(f"Clusters mode requires annotation column. {reason}")
        
        result['cluster_label'] = annotation_col
        result['cluster_label_reason'] = reason
        
        # Training parameters with sensible defaults
        result['num_epochs'] = getattr(params, 'num_epochs', 1000)  # Default for production
        result['learning_rate'] = getattr(params, 'learning_rate', 0.1)
        result['lambda_d'] = getattr(params, 'lambda_d', 0.0)
        result['lambda_g1'] = getattr(params, 'lambda_g1', 1.0)
        result['lambda_g2'] = getattr(params, 'lambda_g2', 0.0)
        result['lambda_r'] = getattr(params, 'lambda_r', 0.0)
        result['lambda_count'] = getattr(params, 'lambda_count', 1.0)
        result['lambda_f_reg'] = getattr(params, 'lambda_f_reg', 1.0)
        
        return result


async def improved_tangram_annotation(adata_sp: ad.AnnData,
                                    params: Any,
                                    data_store: Dict[str, Any],
                                    context: Optional[Any] = None) -> tuple:
    """
    Improved Tangram annotation with better user experience.
    
    Key improvements:
    1. Clear parameter validation with helpful error messages
    2. Centralized column detection logic
    3. Informative progress reporting
    4. Graceful handling of edge cases
    5. Consistent return values
    """
    
    # Step 1: Validate inputs early with clear messages
    if not params.reference_data_id:
        raise ValueError("Reference data ID is required for Tangram annotation. Please specify reference_data_id in parameters.")
    
    if params.reference_data_id not in data_store:
        available_ids = list(data_store.keys())
        raise ValueError(f"Reference dataset '{params.reference_data_id}' not found. Available datasets: {available_ids}")
    
    reference_data = data_store[params.reference_data_id]['adata']
    
    # Step 2: Prepare and validate parameters
    try:
        tangram_params = ImprovedTangramParameters.validate_and_prepare(params, reference_data)
    except ValueError as e:
        if context:
            await context.error(f"Parameter validation failed: {str(e)}")
        raise
    
    # Step 3: Report what we're doing to the user
    if context:
        await context.info(f"Starting Tangram annotation in '{tangram_params['mode']}' mode")
        if tangram_params.get('cluster_label'):
            await context.info(f"Using annotation column: '{tangram_params['cluster_label']}' - {tangram_params['cluster_label_reason']}")
        await context.info(f"Training parameters: {tangram_params['num_epochs']} epochs, learning_rate={tangram_params['learning_rate']}")
    
    # Step 4: Import and run Tangram
    try:
        import tangram as tg
    except ImportError as e:
        error_msg = "Tangram is not installed. Install it with: pip install tangram-sc"
        if context:
            await context.error(error_msg)
        raise ImportError(error_msg) from e
    
    # ... rest of implementation would follow similar patterns ...
    
    # The key point is: every step provides clear feedback to users
    # and handles edge cases gracefully without special case logic
    
    return [], {}, {}, 0.0  # Placeholder return


# Usage example showing improved user experience:
"""
Before (confusing):
> No suitable annotation column found for cells mode projection

After (helpful):
> No standard annotation column found. Available categorical columns: ['cluster_id', 'sample_type']. Consider using 'cluster_id'
> To use a specific column, set cluster_label='cluster_id' in your parameters.

Before (silent failure):
> Tangram completed (but no cell types were actually assigned)

After (clear feedback):  
> Tangram annotation completed successfully
> - Found 5 cell types using 'cell_type' column (auto-detected standard pattern)
> - Annotated 200/200 spatial cells
> - Mapping score: 0.73 (good quality)
> - Results stored in adata.obs['cell_type'] and adata.obsm['tangram_ct_pred']
"""