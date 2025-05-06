"""
Cell type annotation tools for spatial transcriptomics data.
"""

from typing import Dict, List, Optional, Any
import numpy as np
import pandas as pd
import scanpy as sc
from mcp.server.fastmcp import Context

from ..models.data import AnnotationParameters
from ..models.analysis import AnnotationResult


# Default marker genes for common cell types
DEFAULT_MARKER_GENES = {
    "T cells": ["CD3D", "CD3E", "CD3G", "CD8A", "CD4", "IL7R", "CCR7", "LCK", "PTPRC"],
    "B cells": ["CD19", "MS4A1", "CD79A", "CD79B", "BANK1", "CD22", "IGHM", "IGHD"],
    "Macrophages": ["CD68", "CD14", "CSF1R", "FCGR3A", "FCGR1A", "ITGAM", "MARCO", "MSR1"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "ACTA2", "FAP", "THY1"],
    "Epithelial cells": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "CLDN3", "CLDN4", "MUC1"],
    "Neurons": ["RBFOX3", "MAP2", "TUBB3", "SYP", "SNAP25", "GRIN1", "GRIN2A", "GRIA1"],
    "Oligodendrocytes": ["MBP", "MOG", "MAG", "MOBP", "PLP1", "OLIG1", "OLIG2", "SOX10"],
    "Astrocytes": ["GFAP", "AQP4", "SLC1A3", "SLC1A2", "ALDH1L1", "GJA1", "SOX9", "ALDOC"],
    "Endothelial cells": ["PECAM1", "CDH5", "VWF", "CLDN5", "FLT1", "KDR", "TEK", "EMCN"]
}


async def annotate_cell_types(
    data_id: str,
    data_store: Dict[str, Any],
    params: AnnotationParameters = AnnotationParameters(),
    context: Optional[Context] = None
) -> AnnotationResult:
    """Annotate cell types in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Annotation parameters
        context: MCP context

    Returns:
        Annotation result
    """
    if context:
        await context.info(f"Annotating cell types using {params.method} method")

    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    adata = data_store[data_id]["adata"]

    # Different methods would be implemented:
    if params.method == "marker_genes":
        if context:
            await context.info("Using marker genes method for annotation")

        # Use provided marker genes or default ones
        marker_genes = params.marker_genes if params.marker_genes else DEFAULT_MARKER_GENES

        # Check if marker genes exist in the dataset
        all_genes = set(adata.var_names)
        valid_cell_types = []
        valid_marker_genes = {}

        for cell_type, genes in marker_genes.items():
            # Filter genes that exist in the dataset
            existing_genes = [gene for gene in genes if gene in all_genes]
            if existing_genes:
                valid_cell_types.append(cell_type)
                valid_marker_genes[cell_type] = existing_genes
                if context:
                    await context.info(f"Found {len(existing_genes)} marker genes for {cell_type}")
            else:
                if context:
                    await context.warning(f"No marker genes found for {cell_type}")

        if not valid_cell_types:
            raise ValueError("No valid marker genes found for any cell type")

        # Calculate scores for each cell type
        cell_type_scores = {}
        for cell_type, genes in valid_marker_genes.items():
            score_name = f"{cell_type}_score"

            # Use scanpy's score_genes function
            sc.tl.score_genes(adata, gene_list=genes, score_name=score_name, use_raw=True)
            cell_type_scores[cell_type] = score_name

        # Create a DataFrame with all scores
        scores_df = pd.DataFrame(index=adata.obs_names)
        for cell_type, score_name in cell_type_scores.items():
            scores_df[cell_type] = adata.obs[score_name]

        # Assign cell type based on highest score
        adata.obs["cell_type"] = scores_df.idxmax(axis=1)
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

        # Calculate confidence scores (normalized scores)
        confidence_scores = {}
        for cell_type in valid_cell_types:
            # Get cells of this type
            cells_of_type = adata.obs["cell_type"] == cell_type
            if np.sum(cells_of_type) > 0:
                # Calculate mean score for this cell type
                mean_score = np.mean(adata.obs[cell_type_scores[cell_type]][cells_of_type])
                # Normalize to 0-1 range (assuming scores are roughly between -1 and 1)
                confidence = (mean_score + 1) / 2
                confidence_scores[cell_type] = round(min(max(confidence, 0.5), 0.99), 2)
            else:
                confidence_scores[cell_type] = 0.5

    elif params.method == "correlation":
        # Based on correlation with reference
        if context:
            await context.info("Using correlation method for annotation")

        if not params.reference_data:
            raise ValueError("Reference data is required for correlation method")

        # This would be implemented with a reference dataset
        # For now, we'll use the marker genes method as a fallback
        if context:
            await context.warning("Reference data not provided, falling back to marker genes method")

        # Use provided marker genes or default ones
        marker_genes = params.marker_genes if params.marker_genes else DEFAULT_MARKER_GENES

        # Check if marker genes exist in the dataset
        all_genes = set(adata.var_names)
        valid_cell_types = []
        valid_marker_genes = {}

        for cell_type, genes in marker_genes.items():
            # Filter genes that exist in the dataset
            existing_genes = [gene for gene in genes if gene in all_genes]
            if existing_genes:
                valid_cell_types.append(cell_type)
                valid_marker_genes[cell_type] = existing_genes

        # Calculate scores for each cell type
        cell_type_scores = {}
        for cell_type, genes in valid_marker_genes.items():
            score_name = f"{cell_type}_score"
            sc.tl.score_genes(adata, gene_list=genes, score_name=score_name, use_raw=True)
            cell_type_scores[cell_type] = score_name

        # Create a DataFrame with all scores
        scores_df = pd.DataFrame(index=adata.obs_names)
        for cell_type, score_name in cell_type_scores.items():
            scores_df[cell_type] = adata.obs[score_name]

        # Assign cell type based on highest score
        adata.obs["cell_type"] = scores_df.idxmax(axis=1)
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

        # Calculate confidence scores
        confidence_scores = {}
        for cell_type in valid_cell_types:
            cells_of_type = adata.obs["cell_type"] == cell_type
            if np.sum(cells_of_type) > 0:
                mean_score = np.mean(adata.obs[cell_type_scores[cell_type]][cells_of_type])
                confidence = (mean_score + 1) / 2
                confidence_scores[cell_type] = round(min(max(confidence, 0.5), 0.99), 2)
            else:
                confidence_scores[cell_type] = 0.5

    elif params.method in ["supervised", "popv", "gptcelltype", "scrgcl"]:
        # These would be implemented with more advanced methods
        if context:
            await context.warning(f"{params.method} method not fully implemented, falling back to marker genes method")

        # Use provided marker genes or default ones
        marker_genes = params.marker_genes if params.marker_genes else DEFAULT_MARKER_GENES

        # Check if marker genes exist in the dataset
        all_genes = set(adata.var_names)
        valid_cell_types = []
        valid_marker_genes = {}

        for cell_type, genes in marker_genes.items():
            # Filter genes that exist in the dataset
            existing_genes = [gene for gene in genes if gene in all_genes]
            if existing_genes:
                valid_cell_types.append(cell_type)
                valid_marker_genes[cell_type] = existing_genes

        # Calculate scores for each cell type
        cell_type_scores = {}
        for cell_type, genes in valid_marker_genes.items():
            score_name = f"{cell_type}_score"
            sc.tl.score_genes(adata, gene_list=genes, score_name=score_name, use_raw=True)
            cell_type_scores[cell_type] = score_name

        # Create a DataFrame with all scores
        scores_df = pd.DataFrame(index=adata.obs_names)
        for cell_type, score_name in cell_type_scores.items():
            scores_df[cell_type] = adata.obs[score_name]

        # Assign cell type based on highest score
        adata.obs["cell_type"] = scores_df.idxmax(axis=1)
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

        # Calculate confidence scores
        confidence_scores = {}
        for cell_type in valid_cell_types:
            cells_of_type = adata.obs["cell_type"] == cell_type
            if np.sum(cells_of_type) > 0:
                mean_score = np.mean(adata.obs[cell_type_scores[cell_type]][cells_of_type])
                confidence = (mean_score + 1) / 2
                confidence_scores[cell_type] = round(min(max(confidence, 0.5), 0.99), 2)
            else:
                confidence_scores[cell_type] = 0.5

    # Count cells per cell type
    cell_types = list(adata.obs["cell_type"].unique())
    counts = adata.obs["cell_type"].value_counts().to_dict()

    # Update the AnnData object in the data store
    data_store[data_id]["adata"] = adata

    # Return result
    return AnnotationResult(
        data_id=data_id,
        method=params.method,
        cell_types=cell_types,
        counts=counts,
        confidence_scores=confidence_scores
    )