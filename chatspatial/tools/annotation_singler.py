"""
SingleR-based cell type annotation implementation for spatial transcriptomics.

This module provides SingleR integration as a replacement for the simple marker_genes method.
SingleR uses reference-based annotation with correlation analysis for more accurate results.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from typing import Dict, List, Optional, Any, Tuple
import logging

try:
    from mcp.server.fastmcp import Context
except ImportError:
    Context = None

logger = logging.getLogger(__name__)

# Try to import SingleR
SINGLER_AVAILABLE = False
SINGLER_IMPORT_ERROR = None

try:
    import singler
    import singlecellexperiment as sce
    SINGLER_AVAILABLE = True
except ImportError as e:
    SINGLER_IMPORT_ERROR = str(e)
    logger.warning(f"SingleR not available: {e}")

# Try to import celldex for reference datasets
CELLDEX_AVAILABLE = False
try:
    import celldex
    CELLDEX_AVAILABLE = True
except ImportError:
    logger.warning("celldex not available for reference datasets")


class SingleRAnnotator:
    """SingleR-based cell type annotation"""
    
    def __init__(self):
        self.reference_cache = {}
        
    async def annotate_with_singler(
        self, 
        adata,
        reference_data=None,
        reference_labels=None,
        reference_name=None,
        marker_genes=None,
        use_integrated=False,
        num_threads=4,
        context: Optional[Context] = None
    ) -> Tuple[List[str], Dict[str, int], Dict[str, float], Optional[float]]:
        """
        Annotate cells using SingleR method.
        
        Args:
            adata: AnnData object to annotate
            reference_data: Reference dataset (SingleCellExperiment or AnnData)
            reference_labels: Labels for reference data
            reference_name: Name of pre-built reference (e.g., "blueprint_encode", "dice")
            marker_genes: Optional marker genes to use
            use_integrated: Whether to use integrated annotation with multiple references
            num_threads: Number of threads for parallel processing
            context: MCP context for logging
            
        Returns:
            Tuple of (cell_types, counts, confidence_scores, mapping_score)
        """
        
        if not SINGLER_AVAILABLE:
            error_msg = f"SingleR is not installed. {SINGLER_IMPORT_ERROR}\nInstall with: pip install singler singlecellexperiment"
            if context:
                await context.error(error_msg)
            raise ImportError(error_msg)
        
        try:
            if context:
                await context.info("🔬 Starting SingleR annotation...")
                await context.info(f"   Cells: {adata.n_obs}, Genes: {adata.n_vars}")
            
            # Step 1: Convert AnnData to SingleCellExperiment
            if context:
                await context.info("Converting data format...")
            
            # Create SingleCellExperiment from AnnData
            test_sce = sce.SingleCellExperiment.from_anndata(adata)
            
            # Get expression matrix and features
            # Use normalized data if available, otherwise raw
            if "X_normalized" in adata.layers:
                test_mat = adata.layers["X_normalized"]
            elif adata.X is not None:
                test_mat = adata.X
            else:
                test_mat = adata.raw.X if adata.raw else adata.X
            
            # Convert sparse to dense if needed
            if hasattr(test_mat, 'toarray'):
                test_mat = test_mat.toarray()
            
            # Transpose for SingleR (genes x cells)
            test_mat = test_mat.T
            
            # Get feature names
            test_features = list(adata.var_names)
            
            # Step 2: Get or prepare reference
            ref_data, ref_labels = await self._prepare_reference(
                reference_data, reference_labels, reference_name, 
                marker_genes, context
            )
            
            if context:
                await context.info(f"Reference prepared with {len(set(ref_labels))} cell types")
            
            # Step 3: Run SingleR annotation
            if context:
                await context.info("Running SingleR classification...")
            
            if use_integrated and isinstance(ref_data, list):
                # Integrated annotation with multiple references
                single_results, integrated = singler.annotate_integrated(
                    test_mat,
                    ref_data=ref_data,
                    ref_labels=ref_labels,
                    test_features=test_features,
                    num_threads=num_threads
                )
                
                # Use integrated results
                best_labels = integrated.column("best_label")
                scores = integrated.column("scores")
                
            else:
                # Single reference annotation
                results = singler.annotate_single(
                    test_data=test_mat,
                    test_features=test_features,
                    ref_data=ref_data,
                    ref_labels=ref_labels,
                    num_threads=num_threads
                )
                
                best_labels = results.column("best")
                scores = results.column("scores")
            
            # Step 4: Process results
            cell_types = list(best_labels)
            
            # Calculate statistics
            unique_types = list(set(cell_types))
            counts = pd.Series(cell_types).value_counts().to_dict()
            
            # Calculate confidence scores
            confidence_scores = {}
            
            if scores is not None:
                # scores is a BiocFrame where columns are cell types
                # Convert BiocFrame to pandas DataFrame
                try:
                    scores_df = pd.DataFrame(scores.to_dict())
                except AttributeError:
                    # If to_dict doesn't work, try converting to numpy array
                    scores_df = pd.DataFrame(scores.to_numpy() if hasattr(scores, 'to_numpy') else scores)
                
                # For each cell type, calculate average confidence
                for cell_type in unique_types:
                    # Get cells assigned to this type
                    mask = [ct == cell_type for ct in cell_types]
                    
                    if cell_type in scores_df.columns and any(mask):
                        # Average score for cells assigned to this type
                        type_scores = scores_df.loc[mask, cell_type]
                        avg_score = type_scores.mean()
                        
                        # Normalize to 0-1 range
                        # SingleR scores are typically correlations (-1 to 1)
                        # Convert to 0-1 confidence
                        confidence = (avg_score + 1) / 2
                        confidence_scores[cell_type] = float(confidence)
                    else:
                        confidence_scores[cell_type] = 0.5
            else:
                # Default confidence if scores not available
                for cell_type in unique_types:
                    confidence_scores[cell_type] = 0.8
            
            # Step 5: Add to AnnData
            adata.obs['singler_celltype'] = cell_types
            adata.obs['singler_celltype'] = adata.obs['singler_celltype'].astype('category')
            
            # Add confidence scores
            confidence_array = [confidence_scores.get(ct, 0.5) for ct in cell_types]
            adata.obs['singler_confidence'] = confidence_array
            
            if context:
                await context.info(f"✅ SingleR annotation completed!")
                await context.info(f"   Found {len(unique_types)} cell types")
                top_types = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:3]
                await context.info(f"   Top types: {', '.join([f'{t}({c})' for t, c in top_types])}")
            
            return unique_types, counts, confidence_scores, None
            
        except Exception as e:
            error_msg = f"SingleR annotation failed: {str(e)}"
            if context:
                await context.error(error_msg)
            raise ValueError(error_msg)
    
    async def _prepare_reference(
        self, 
        reference_data, 
        reference_labels, 
        reference_name,
        marker_genes,
        context
    ):
        """Prepare reference dataset for SingleR"""
        
        # Priority: reference_name > reference_data > marker_genes
        if reference_name and CELLDEX_AVAILABLE:
            # Use pre-built reference from celldex
            if context:
                await context.info(f"Loading reference: {reference_name}")
            
            # Check cache
            if reference_name in self.reference_cache:
                ref = self.reference_cache[reference_name]
            else:
                # Fetch reference
                ref = celldex.fetch_reference(
                    reference_name, 
                    "2024-02-26", 
                    realize_assays=True
                )
                self.reference_cache[reference_name] = ref
            
            # Get labels
            if reference_labels is None:
                # Try common label columns
                for label_col in ["label.main", "label.fine", "cell_type"]:
                    try:
                        reference_labels = ref.get_column_data().column(label_col)
                        break
                    except:
                        continue
                
                if reference_labels is None:
                    raise ValueError(f"Could not find labels in reference {reference_name}")
            
            return ref, reference_labels
        
        elif reference_data is not None:
            # Use provided reference data
            if context:
                await context.info("Using provided reference data")
            
            # Convert AnnData to SingleCellExperiment if needed
            if hasattr(reference_data, 'X'):  # AnnData
                ref_sce = sce.SingleCellExperiment.from_anndata(reference_data)
                
                # Get labels from AnnData
                if reference_labels is None:
                    # Try to find cell type column
                    for col in ['cell_type', 'celltype', 'CellType', 'leiden']:
                        if col in reference_data.obs.columns:
                            reference_labels = list(reference_data.obs[col])
                            break
                    
                    if reference_labels is None:
                        raise ValueError("No cell type labels found in reference data")
                
                return ref_sce, reference_labels
            else:
                # Already SingleCellExperiment
                return reference_data, reference_labels
        
        elif marker_genes:
            # Note: marker genes alone cannot create a proper reference for SingleR
            # SingleR requires actual expression profiles from reference cells
            if context:
                await context.info("Warning: SingleR requires reference expression data, not just marker genes")
                await context.info("Falling back to default BlueprintEncode reference")
            
            # Fall back to default reference when only marker genes are provided
            if CELLDEX_AVAILABLE:
                ref = celldex.fetch_reference(
                    "blueprint_encode", 
                    "2024-02-26", 
                    realize_assays=True
                )
                ref_labels = ref.get_column_data().column("label.main")
                return ref, ref_labels
            else:
                raise ValueError(
                    "SingleR requires reference expression data. Please either:\n"
                    "1. Provide a reference_name (e.g., 'blueprint_encode')\n"
                    "2. Provide reference_data with expression profiles\n"
                    "3. Install celldex for access to pre-built references"
                )
        
        else:
            # Try to use default references
            if CELLDEX_AVAILABLE:
                if context:
                    await context.info("Using default BlueprintEncode reference")
                
                ref = celldex.fetch_reference(
                    "blueprint_encode", 
                    "2024-02-26", 
                    realize_assays=True
                )
                ref_labels = ref.get_column_data().column("label.main")
                
                return ref, ref_labels
            else:
                raise ValueError(
                    "No reference data provided. Please provide either:\n"
                    "1. reference_data and reference_labels\n"
                    "2. reference_name (requires celldex)\n"
                    "3. marker_genes for synthetic reference"
                )
    


# Global annotator instance
_singler_annotator = SingleRAnnotator()


async def annotate_with_singler(
    adata,
    params,
    data_store=None,
    context: Optional[Context] = None
) -> Tuple[List[str], Dict[str, int], Dict[str, float], Optional[float]]:
    """
    Main entry point for SingleR annotation.
    
    This function replaces the old marker_genes implementation.
    """
    
    # Extract parameters
    reference_data = None
    reference_labels = None
    reference_name = params.singler_reference if hasattr(params, 'singler_reference') else None
    
    # Check for reference data in data_store
    if hasattr(params, 'reference_data_id') and params.reference_data_id and data_store:
        if params.reference_data_id in data_store:
            reference_data = data_store[params.reference_data_id]["adata"]
    
    # Use marker genes as fallback
    marker_genes = params.marker_genes if hasattr(params, 'marker_genes') else None
    
    # Run SingleR
    return await _singler_annotator.annotate_with_singler(
        adata=adata,
        reference_data=reference_data,
        reference_labels=reference_labels,
        reference_name=reference_name,
        marker_genes=marker_genes,
        use_integrated=getattr(params, 'singler_integrated', False),
        num_threads=getattr(params, 'num_threads', 4),
        context=context
    )