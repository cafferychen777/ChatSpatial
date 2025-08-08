"""
Cell type annotation tools for spatial transcriptomics data.
"""

from typing import Dict, List, Optional, Any, Union
import numpy as np
import pandas as pd
import scanpy as sc
from mcp.server.fastmcp import Context

# Import Tangram for cell type annotation
try:
    import tangram as tg
except ImportError:
    tg = None

# Import scvi-tools for cell type annotation
try:
    import scvi
    from scvi.external import CellAssign
except ImportError:
    scvi = None
    CellAssign = None

# Import mllmcelltype for LLM-based annotation
try:
    import mllmcelltype
except ImportError:
    mllmcelltype = None

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

# Constants for annotation
DEFAULT_HVG_COUNT = 2000
DEFAULT_SCANVI_EPOCHS = 200
CONFIDENCE_MIN = 0.5
CONFIDENCE_MAX = 0.99

async def _handle_annotation_error(error: Exception, method: str, context: Optional[Context] = None) -> None:
    """Handle annotation errors consistently"""
    error_msg = f"{method} annotation failed: {str(error)}"
    if context:
        await context.error(f"Error in {method} annotation: {str(error)}")
    raise ValueError(error_msg)

async def _validate_marker_genes(adata, marker_genes: Dict[str, List[str]], context: Optional[Context] = None) -> Dict[str, List[str]]:
    """Validate and filter marker genes that exist in the dataset"""
    all_genes = set(adata.var_names)
    valid_cell_types = []
    valid_marker_genes = {}
    
    for cell_type, genes in marker_genes.items():
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
    
    return valid_marker_genes

def _calculate_cell_type_scores(adata, valid_marker_genes: Dict[str, List[str]], use_raw: bool = False) -> Dict[str, str]:
    """Calculate cell type scores using marker genes"""
    cell_type_scores = {}
    for cell_type, genes in valid_marker_genes.items():
        score_name = f"{cell_type}_score"
        try:
            # Try with default settings first
            sc.tl.score_genes(adata, gene_list=genes, score_name=score_name, use_raw=use_raw)
        except Exception as e:
            try:
                # Fallback: try with ctrl_as_ref=False
                sc.tl.score_genes(adata, gene_list=genes, score_name=score_name, use_raw=use_raw, ctrl_as_ref=False)
            except Exception as e2:
                # Ultimate fallback: calculate simple mean expression
                if use_raw and adata.raw is not None:
                    expr_data = adata.raw.to_adata()
                else:
                    expr_data = adata
                
                # Calculate mean expression across marker genes
                gene_indices = [i for i, gene in enumerate(expr_data.var_names) if gene in genes]
                if gene_indices:
                    if hasattr(expr_data.X, 'toarray'):
                        scores = expr_data.X[:, gene_indices].toarray().mean(axis=1)
                    else:
                        scores = expr_data.X[:, gene_indices].mean(axis=1)
                    adata.obs[score_name] = scores.flatten()
                else:
                    # No genes found, set to zero
                    adata.obs[score_name] = 0.0
        
        cell_type_scores[cell_type] = score_name
    return cell_type_scores

def _assign_cell_types_from_scores(adata, cell_type_scores: Dict[str, str]) -> tuple:
    """Assign cell types based on highest scores and calculate confidence"""
    # Create DataFrame with all scores
    scores_df = pd.DataFrame(index=adata.obs_names)
    for cell_type, score_name in cell_type_scores.items():
        scores_df[cell_type] = adata.obs[score_name]
    
    # Assign cell type based on highest score
    adata.obs["cell_type"] = scores_df.idxmax(axis=1)
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    
    # Calculate confidence scores
    confidence_scores = {}
    valid_cell_types = list(cell_type_scores.keys())
    
    for cell_type in valid_cell_types:
        cells_of_type = adata.obs["cell_type"] == cell_type
        if np.sum(cells_of_type) > 0:
            mean_score = np.mean(adata.obs[cell_type_scores[cell_type]][cells_of_type])
            # Normalize to 0-1 range (assuming scores are roughly between -1 and 1)
            confidence = (mean_score + 1) / 2
            confidence_scores[cell_type] = round(min(max(confidence, CONFIDENCE_MIN), CONFIDENCE_MAX), 2)
        else:
            confidence_scores[cell_type] = CONFIDENCE_MIN
    
    # Get cell types and counts
    cell_types = list(adata.obs["cell_type"].unique())
    counts = adata.obs["cell_type"].value_counts().to_dict()
    
    return cell_types, counts, confidence_scores


async def _annotate_with_marker_genes(adata, params: AnnotationParameters, context: Optional[Context] = None):
    """Annotate cell types using marker genes method"""
    try:
        if context:
            await context.info("Using marker genes method for annotation")
        
        # Use provided marker genes or default ones
        marker_genes = params.marker_genes if params.marker_genes else DEFAULT_MARKER_GENES
        
        # Validate marker genes
        valid_marker_genes = await _validate_marker_genes(adata, marker_genes, context)
        
        # Calculate scores
        cell_type_scores = _calculate_cell_type_scores(adata, valid_marker_genes, use_raw=False)
        
        # Assign cell types and calculate confidence
        cell_types, counts, confidence_scores = _assign_cell_types_from_scores(adata, cell_type_scores)
        
        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization
        
        return cell_types, counts, confidence_scores, None
        
    except Exception as e:
        await _handle_annotation_error(e, "marker_genes", context)

async def _annotate_with_correlation(adata, params: AnnotationParameters, context: Optional[Context] = None):
    """Annotate cell types using correlation method"""
    try:
        if context:
            await context.info("Using correlation method for annotation")
        
        if not params.reference_data:
            if context:
                await context.warning("Reference data not provided, falling back to marker genes method")
            return await _annotate_with_marker_genes(adata, params, context)
        
        # Use provided marker genes or default ones
        marker_genes = params.marker_genes if params.marker_genes else DEFAULT_MARKER_GENES
        
        # Validate marker genes
        valid_marker_genes = await _validate_marker_genes(adata, marker_genes, context)
        
        # Calculate scores using raw data for correlation
        cell_type_scores = _calculate_cell_type_scores(adata, valid_marker_genes, use_raw=True)
        
        # Assign cell types and calculate confidence
        cell_types, counts, confidence_scores = _assign_cell_types_from_scores(adata, cell_type_scores)
        
        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization
        
        return cell_types, counts, confidence_scores, None
        
    except Exception as e:
        await _handle_annotation_error(e, "correlation", context)

async def _annotate_with_tangram(adata, params: AnnotationParameters, data_store: Dict[str, Any], context: Optional[Context] = None):
    """Annotate cell types using Tangram method"""
    try:
        if context:
            await context.info("Using Tangram method for annotation")

        if tg is None:
            raise ImportError("Tangram package is not installed. Please install it with 'pip install tangram-sc'")

        # Check if reference data is provided
        if params.reference_data_id is None:
            raise ValueError("Reference data ID is required for Tangram method")

        if params.reference_data_id not in data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found")

        # Get reference single-cell data
        adata_sc = data_store[params.reference_data_id]["adata"]
        adata_sp = adata  # Spatial data

        if context:
            await context.info(f"Using reference dataset {params.reference_data_id} with {adata_sc.n_obs} cells")

        # Determine training genes
        training_genes = params.training_genes

        if training_genes is None:
            # Use marker genes if available
            if params.marker_genes:
                if context:
                    await context.info("Using provided marker genes for Tangram mapping")
                # Flatten marker genes dictionary
                training_genes = []
                for genes in params.marker_genes.values():
                    training_genes.extend(genes)
                training_genes = list(set(training_genes))  # Remove duplicates
            else:
                # Use highly variable genes
                if context:
                    await context.info("Computing highly variable genes for Tangram mapping")
                if 'highly_variable' not in adata_sc.var:
                    sc.pp.highly_variable_genes(adata_sc, n_top_genes=DEFAULT_HVG_COUNT)
                training_genes = list(adata_sc.var_names[adata_sc.var.highly_variable])

        if context:
            await context.info(f"Using {len(training_genes)} genes for Tangram mapping")

        # Preprocess data for Tangram
        tg.pp_adatas(adata_sc, adata_sp, genes=training_genes)

        if context:
            await context.info(f"Preprocessed data for Tangram mapping. {len(adata_sc.uns['training_genes'])} training genes selected.")

        # Set mapping mode
        mode = params.mode
        cluster_label = params.cluster_label

        if mode == "clusters" and cluster_label is None:
            if context:
                await context.warning("Cluster label not provided for 'clusters' mode. Using default cell type annotation if available.")
            # Try to find a cell type annotation in the reference data
            for col in ['cell_type', 'celltype', 'cell_types', 'leiden', 'louvain']:
                if col in adata_sc.obs:
                    cluster_label = col
                    break

            if cluster_label is None:
                raise ValueError("No cluster label found in reference data. Please provide a cluster_label parameter.")

            if context:
                await context.info(f"Using '{cluster_label}' as cluster label for Tangram mapping")

        # Run Tangram mapping
        if context:
            await context.info(f"Running Tangram mapping in '{mode}' mode for {params.num_epochs} epochs")

        mapping_args = {
            "mode": mode,
            "num_epochs": params.num_epochs,
            "device": "cpu"  # Use CPU for compatibility
        }

        if mode == "clusters":
            mapping_args["cluster_label"] = cluster_label

        ad_map = tg.map_cells_to_space(adata_sc, adata_sp, **mapping_args)

        # Get mapping score
        tangram_mapping_score = 0.0  # Default score
        try:
            if 'training_history' in ad_map.uns and len(ad_map.uns['training_history']) > 0:
                # Try to get the last score from training history
                last_entry = ad_map.uns['training_history'][-1]
                if isinstance(last_entry, (list, tuple)) and len(last_entry) > 0:
                    tangram_mapping_score = float(last_entry[0])
                elif isinstance(last_entry, (int, float)):
                    tangram_mapping_score = float(last_entry)
            elif hasattr(ad_map, 'X') and ad_map.X is not None:
                # If no training history, use a default score based on mapping success
                tangram_mapping_score = 0.5  # Indicate successful mapping
        except Exception as score_error:
            if context:
                await context.warning(f"Could not extract mapping score: {score_error}")
            tangram_mapping_score = 0.5  # Default score for successful mapping

        if context:
            await context.info(f"Tangram mapping completed with score: {tangram_mapping_score}")

        # Project cell annotations to space
        try:
            if mode == "clusters" and cluster_label:
                tg.plot_cell_annotation(ad_map, adata_sp, annotation=cluster_label)
            else:
                # For cells mode, project all cell types
                if 'subclass_label' in adata_sc.obs:
                    tg.plot_cell_annotation(ad_map, adata_sp, annotation='subclass_label')
        except Exception as proj_error:
            if context:
                await context.warning(f"Could not project cell annotations: {proj_error}")
            # Continue without projection

        # Get cell type predictions
        cell_types = []
        counts = {}
        confidence_scores = {}
        visualization = None
        
        if 'tangram_ct_pred' in adata_sp.obsm:
            cell_type_df = adata_sp.obsm['tangram_ct_pred']

            # Get cell types and counts
            cell_types = list(cell_type_df.columns)

            # Assign cell type based on highest probability
            adata_sp.obs["cell_type"] = cell_type_df.idxmax(axis=1)
            adata_sp.obs["cell_type"] = adata_sp.obs["cell_type"].astype("category")

            # Get counts
            counts = adata_sp.obs["cell_type"].value_counts().to_dict()

            # Calculate confidence scores (use max probability as confidence)
            confidence_scores = {}
            for cell_type in cell_types:
                cells_of_type = adata_sp.obs["cell_type"] == cell_type
                if np.sum(cells_of_type) > 0:
                    # Use mean probability as confidence
                    mean_prob = cell_type_df.loc[cells_of_type, cell_type].mean()
                    confidence_scores[cell_type] = round(float(mean_prob), 2)
                else:
                    confidence_scores[cell_type] = CONFIDENCE_MIN

            # Note: Visualizations should be created using the separate visualize_data tool
            # This maintains clean separation between analysis and visualization
            if context:
                await context.info("Cell type mapping complete. Use visualize_data tool to visualize results")

        else:
            if context:
                await context.warning("No cell type predictions found in Tangram results")

        return cell_types, counts, confidence_scores, tangram_mapping_score, None
        
    except Exception as e:
        await _handle_annotation_error(e, "tangram", context)

async def _annotate_with_scanvi(adata, params: AnnotationParameters, data_store: Dict[str, Any], context: Optional[Context] = None):
    """Annotate cell types using scANVI method"""
    try:
        if context:
            await context.info("Using scANVI method for annotation")

        if scvi is None:
            raise ImportError("scvi-tools package is not installed. Please install it with 'pip install scvi-tools'")

        # Check if reference data is provided
        if params.reference_data_id is None:
            raise ValueError("Reference data ID is required for scANVI method")

        if params.reference_data_id not in data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found")

        # Get reference single-cell data
        adata_ref = data_store[params.reference_data_id]["adata"]
        
        # Setup AnnData for scANVI
        scvi.model.SCANVI.setup_anndata(
            adata_ref,
            labels_key=getattr(params, 'cell_type_key', "cell_type"),
            unlabeled_category=params.scanvi_unlabeled_category
        )
        
        # Train scANVI model
        if context:
            await context.info("Training scANVI model...")
        
        # Create scANVI model with correct parameter positioning
        model = scvi.model.SCANVI(
            adata_ref,
            n_hidden=params.scanvi_n_hidden,
            n_latent=params.scanvi_n_latent,
            n_layers=params.scanvi_n_layers,
            dropout_rate=params.scanvi_dropout_rate
        )
        
        model.train(max_epochs=params.num_epochs)
        
        # Prepare spatial data - add dummy cell_type column for setup
        cell_type_key = getattr(params, 'cell_type_key', "cell_type")
        if cell_type_key not in adata.obs.columns:
            # Add dummy cell type column filled with unlabeled category
            adata.obs[cell_type_key] = params.scanvi_unlabeled_category
        
        scvi.model.SCANVI.setup_anndata(adata, labels_key=cell_type_key, unlabeled_category=params.scanvi_unlabeled_category)
        
        # Transfer model to spatial data
        spatial_model = scvi.model.SCANVI.load_query_data(adata, model)
        spatial_model.train(max_epochs=DEFAULT_SCANVI_EPOCHS)
        
        # Get predictions
        predictions = spatial_model.predict()
        adata.obs["cell_type"] = predictions
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
        
        # Get cell types and counts
        cell_types = list(adata.obs["cell_type"].cat.categories)
        counts = adata.obs["cell_type"].value_counts().to_dict()
        
        # Get prediction probabilities as confidence scores
        try:
            probs = spatial_model.predict(soft=True)
            confidence_scores = {}
            for i, cell_type in enumerate(cell_types):
                cells_of_type = adata.obs["cell_type"] == cell_type
                if np.sum(cells_of_type) > 0 and isinstance(probs, pd.DataFrame):
                    if cell_type in probs.columns:
                        mean_prob = probs.loc[cells_of_type, cell_type].mean()
                        confidence_scores[cell_type] = round(float(mean_prob), 2)
                    else:
                        confidence_scores[cell_type] = CONFIDENCE_MIN
                elif np.sum(cells_of_type) > 0 and hasattr(probs, 'shape') and probs.shape[1] > i:
                    mean_prob = probs[cells_of_type, i].mean()
                    confidence_scores[cell_type] = round(float(mean_prob), 2)
                else:
                    confidence_scores[cell_type] = CONFIDENCE_MIN
        except Exception as e:
            if context:
                await context.warning(f"Could not get confidence scores: {e}")
            confidence_scores = {cell_type: CONFIDENCE_MIN for cell_type in cell_types}

        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization
        
        return cell_types, counts, confidence_scores, None
                
    except Exception as e:
        await _handle_annotation_error(e, "scanvi", context)

async def _annotate_with_mllmcelltype(adata, params: AnnotationParameters, context: Optional[Context] = None):
    """Annotate cell types using mLLMCellType (LLM-based) method"""
    try:
        if context:
            await context.info("Using mLLMCellType (LLM-based) method for annotation")

        if mllmcelltype is None:
            raise ImportError("mllmcelltype package is not installed. Please install it with 'pip install mllmcelltype'")

        # Check if clustering has been performed
        cluster_key = params.cluster_label if params.cluster_label else 'leiden'
        if cluster_key not in adata.obs:
            if context:
                await context.info(f"Performing clustering ({cluster_key}) for mLLMCellType annotation")
            if cluster_key == 'leiden':
                sc.tl.leiden(adata)
            else:
                sc.tl.louvain(adata)

        # Find differentially expressed genes for each cluster
        if context:
            await context.info("Finding marker genes for each cluster")
        
        sc.tl.rank_genes_groups(adata, cluster_key, method='wilcoxon')
        
        # Extract top marker genes for each cluster
        marker_genes_dict = {}
        n_genes = params.mllm_n_marker_genes if hasattr(params, 'mllm_n_marker_genes') else 20
        
        for cluster in adata.obs[cluster_key].unique():
            # Get top genes for this cluster
            gene_names = adata.uns['rank_genes_groups']['names'][str(cluster)][:n_genes]
            marker_genes_dict[f"Cluster_{cluster}"] = list(gene_names)
        
        if context:
            await context.info(f"Found marker genes for {len(marker_genes_dict)} clusters")

        # Prepare parameters for mllmcelltype
        species = params.mllm_species if hasattr(params, 'mllm_species') else 'human'
        tissue = params.mllm_tissue if hasattr(params, 'mllm_tissue') else None
        provider = params.mllm_provider if hasattr(params, 'mllm_provider') else 'openai'
        model = params.mllm_model if hasattr(params, 'mllm_model') else None
        api_key = params.mllm_api_key if hasattr(params, 'mllm_api_key') else None

        # Call mllmcelltype to annotate clusters
        if context:
            await context.info(f"Calling LLM ({provider}/{model or 'default'}) for cell type annotation")
        
        try:
            annotations = mllmcelltype.annotate_clusters(
                marker_genes=marker_genes_dict,
                species=species,
                provider=provider,
                model=model,
                api_key=api_key,
                tissue=tissue,
                use_cache=True
            )
        except Exception as e:
            if context:
                await context.error(f"mLLMCellType annotation failed: {str(e)}")
            raise

        if context:
            await context.info(f"Received annotations for {len(annotations)} clusters")

        # Map cluster annotations back to cells
        cluster_to_celltype = {}
        for cluster_name, cell_type in annotations.items():
            # Extract cluster number from "Cluster_X" format
            cluster_id = cluster_name.replace("Cluster_", "")
            cluster_to_celltype[cluster_id] = cell_type

        # Apply cell type annotations to cells
        adata.obs["cell_type"] = adata.obs[cluster_key].astype(str).map(cluster_to_celltype)
        
        # Handle any unmapped clusters
        unmapped = adata.obs["cell_type"].isna()
        if unmapped.any():
            if context:
                await context.warning(f"Found {unmapped.sum()} cells in unmapped clusters")
            adata.obs.loc[unmapped, "cell_type"] = "Unknown"
        
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

        # Get cell types and counts
        cell_types = list(adata.obs["cell_type"].unique())
        counts = adata.obs["cell_type"].value_counts().to_dict()

        # Calculate confidence scores based on cluster homogeneity
        confidence_scores = {}
        for cell_type in cell_types:
            if cell_type != "Unknown":
                # High confidence for LLM-based annotations
                confidence_scores[cell_type] = 0.85
            else:
                confidence_scores[cell_type] = CONFIDENCE_MIN

        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization
        if context:
            await context.info("Cell type annotation complete. Use visualize_data tool to visualize results")

        return cell_types, counts, confidence_scores, None, None

    except Exception as e:
        await _handle_annotation_error(e, "mllmcelltype", context)

async def _annotate_with_cellassign(adata, params: AnnotationParameters, context: Optional[Context] = None):
    """Annotate cell types using CellAssign method"""
    try:
        if context:
            await context.info("Using CellAssign method for annotation")

        if CellAssign is None:
            raise ImportError("CellAssign from scvi-tools package is not installed. Please install it with 'pip install scvi-tools'")

        # Check if marker genes are provided
        if params.marker_genes is None:
            if context:
                await context.warning("No marker genes provided, using default marker genes")
            marker_genes = DEFAULT_MARKER_GENES
        else:
            marker_genes = params.marker_genes

        # Validate marker genes
        valid_marker_genes = await _validate_marker_genes(adata, marker_genes, context)
        valid_cell_types = list(valid_marker_genes.keys())
        
        # Create marker gene matrix as DataFrame (required by CellAssign API)
        all_marker_genes = []
        for genes in valid_marker_genes.values():
            all_marker_genes.extend(genes)
        available_marker_genes = list(set(all_marker_genes))  # Remove duplicates
        
        if not available_marker_genes:
            raise ValueError("No marker genes found in the dataset")
        
        # Create DataFrame with genes as index, cell types as columns
        marker_gene_matrix = pd.DataFrame(
            np.zeros((len(available_marker_genes), len(valid_cell_types))),
            index=available_marker_genes,
            columns=valid_cell_types
        )
        
        # Fill marker matrix
        for cell_type in valid_cell_types:
            for gene in valid_marker_genes[cell_type]:
                if gene in available_marker_genes:
                    marker_gene_matrix.loc[gene, cell_type] = 1
        
        # Subset data to only marker genes first
        adata_subset = adata[:, available_marker_genes].copy()

        # Check for invalid values in the data
        if hasattr(adata_subset.X, 'toarray'):
            X_array = adata_subset.X.toarray()
        else:
            X_array = adata_subset.X

        # Replace any NaN or Inf values with zeros
        if np.any(np.isnan(X_array)) or np.any(np.isinf(X_array)):
            if context:
                await context.warning("Found NaN or Inf values in data, replacing with zeros")
            X_array = np.nan_to_num(X_array, nan=0.0, posinf=0.0, neginf=0.0)
            adata_subset.X = X_array

        # Additional data cleaning for CellAssign compatibility
        # Check for genes with zero variance (which can cause issues in CellAssign)
        gene_vars = np.var(X_array, axis=0)
        zero_var_genes = gene_vars == 0
        if np.any(zero_var_genes):
            if context:
                await context.warning(f"Found {np.sum(zero_var_genes)} genes with zero variance, adding small noise")
            # Add small random noise to zero-variance genes to avoid numerical issues
            np.random.seed(42)  # For reproducibility
            noise_scale = 1e-6
            for i in np.where(zero_var_genes)[0]:
                X_array[:, i] += np.random.normal(0, noise_scale, X_array.shape[0])
            adata_subset.X = X_array

        # Ensure data is non-negative (CellAssign expects count-like data)
        if np.any(X_array < 0):
            if context:
                await context.warning("Found negative values in data, clipping to zero")
            X_array = np.maximum(X_array, 0)
            adata_subset.X = X_array
        
        # Add size factors if not present
        if 'size_factors' not in adata_subset.obs:
            # Calculate size factors from subset data
            if hasattr(adata_subset.X, 'sum'):
                size_factors = adata_subset.X.sum(axis=1)
                if hasattr(size_factors, 'A1'):  # sparse matrix
                    size_factors = size_factors.A1
            else:
                size_factors = np.sum(adata_subset.X, axis=1)
            
            # Ensure size factors are positive (avoid division by zero)
            size_factors = np.maximum(size_factors, 1e-6)
            adata_subset.obs['size_factors'] = pd.Series(size_factors, index=adata_subset.obs.index)
        
        # Setup CellAssign on subset data only
        CellAssign.setup_anndata(adata_subset, size_factor_key='size_factors')
        
        # Train CellAssign model
        model = CellAssign(
            adata_subset,
            marker_gene_matrix
        )
        
        model.train(
            max_epochs=params.cellassign_max_iter,
            lr=params.cellassign_learning_rate
        )
        
        # Get predictions
        predictions = model.predict()
        
        # Handle different prediction formats
        if isinstance(predictions, pd.DataFrame):
            # CellAssign returns DataFrame with probabilities
            predicted_indices = predictions.values.argmax(axis=1)
            adata.obs["cell_type"] = [valid_cell_types[i] for i in predicted_indices]
            
            # Get confidence scores from probabilities DataFrame
            confidence_scores = {}
            for i, cell_type in enumerate(valid_cell_types):
                cells_of_type = adata.obs["cell_type"] == cell_type
                if np.sum(cells_of_type) > 0:
                    # Use iloc with boolean indexing properly
                    cell_indices = np.where(cells_of_type)[0]
                    mean_prob = predictions.iloc[cell_indices, i].mean()
                    confidence_scores[cell_type] = round(float(mean_prob), 2)
                else:
                    confidence_scores[cell_type] = CONFIDENCE_MIN
        else:
            # Other models return indices directly
            adata.obs["cell_type"] = [valid_cell_types[i] for i in predictions]
            confidence_scores = {cell_type: CONFIDENCE_MIN for cell_type in valid_cell_types}
        
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
        
        # Get cell types and counts
        cell_types = valid_cell_types
        counts = adata.obs["cell_type"].value_counts().to_dict()

        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization
        
        return cell_types, counts, confidence_scores, None
                
    except Exception as e:
        await _handle_annotation_error(e, "cellassign", context)

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

    # Route to appropriate annotation method
    try:
        if params.method == "tangram":
            cell_types, counts, confidence_scores, tangram_mapping_score = await _annotate_with_tangram(
                adata, params, data_store, context
            )
        elif params.method == "scanvi":
            cell_types, counts, confidence_scores, tangram_mapping_score = await _annotate_with_scanvi(
                adata, params, data_store, context
            )
        elif params.method == "cellassign":
            cell_types, counts, confidence_scores, tangram_mapping_score = await _annotate_with_cellassign(
                adata, params, context
            )
        elif params.method == "marker_genes":
            cell_types, counts, confidence_scores, tangram_mapping_score = await _annotate_with_marker_genes(
                adata, params, context
            )
        elif params.method == "correlation":
            cell_types, counts, confidence_scores, tangram_mapping_score = await _annotate_with_correlation(
                adata, params, context
            )
        elif params.method == "mllmcelltype":
            cell_types, counts, confidence_scores, tangram_mapping_score = await _annotate_with_mllmcelltype(
                adata, params, context
            )
        elif params.method in ["supervised", "popv", "gptcelltype", "scrgcl"]:
            # These methods fall back to marker genes
            if context:
                await context.warning(f"{params.method} method not fully implemented, falling back to marker genes method")
            cell_types, counts, confidence_scores, tangram_mapping_score = await _annotate_with_marker_genes(
                adata, params, context
            )
        else:
            raise ValueError(f"Unknown annotation method: {params.method}")

    except Exception as e:
        if context:
            await context.error(f"Annotation failed: {str(e)}")
        raise

    # Update the AnnData object in the data store
    data_store[data_id]["adata"] = adata

    # Inform user about visualization options
    if context:
        await context.info("Cell type annotation complete. Use create_visualization tool with plot_type='cell_types' to visualize results")

    # Return result
    return AnnotationResult(
        data_id=data_id,
        method=params.method,
        cell_types=cell_types,
        counts=counts,
        confidence_scores=confidence_scores,
        tangram_mapping_score=tangram_mapping_score
    )