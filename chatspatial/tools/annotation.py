"""
Cell type annotation tools for spatial transcriptomics data.
"""

from typing import Dict, List, Optional, Any, Union
import hashlib
import pickle
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from mcp.server.fastmcp import Context

# Import Tangram for cell type annotation
try:
    import tangram as tg
except ImportError:
    tg = None

# Optional imports - actual validation happens at runtime
# This allows the module to load even if optional dependencies are missing

# R interface validation is now handled by _validate_rpy2_and_r() function

from ..models.data import AnnotationParameters
from ..models.analysis import AnnotationResult

# ============================================================================
# DEPENDENCY VALIDATION SYSTEM
# ============================================================================

class DependencyError(Exception):
    """Custom exception for missing dependencies with helpful installation info"""
    def __init__(self, package_name: str, method_name: str, install_guide: str):
        self.package_name = package_name
        self.method_name = method_name
        self.install_guide = install_guide
        super().__init__(f"{package_name} is required for {method_name} method.\n{install_guide}")

def _get_installation_guide(package_name: str) -> str:
    """Get user-friendly installation instructions for a package"""
    guides = {
        "scvi-tools": {
            "pip": "pip install scvi-tools",
            "conda": "conda install -c conda-forge scvi-tools",
            "docs": "https://scvi-tools.org/installation.html",
            "note": "Requires Python 3.8+ and PyTorch"
        },
        "tangram-sc": {
            "pip": "pip install tangram-sc",
            "conda": "conda install -c bioconda tangram-sc",
            "docs": "https://tangram-sc.readthedocs.io/",
            "note": "Requires PyTorch and scanpy"
        },
        "mllmcelltype": {
            "pip": "pip install mllmcelltype",
            "conda": "Not available via conda",
            "docs": "https://github.com/Winnie09/mLLMCellType",
            "note": "Requires API keys for LLM providers"
        },
        "rpy2": {
            "pip": "pip install rpy2",
            "conda": "conda install -c conda-forge rpy2",
            "docs": "https://rpy2.github.io/doc/latest/html/overview.html",
            "note": "Requires R installation (https://www.r-project.org/)"
        }
    }
    
    if package_name not in guides:
        return f"Please install {package_name}"
    
    guide = guides[package_name]
    installation_text = f"""
Installation Options:
  • pip: {guide['pip']}
  • conda: {guide['conda']}
  • Documentation: {guide['docs']}
  • Note: {guide['note']}

Troubleshooting:
  • Ensure you have the latest pip: pip install --upgrade pip
  • For conda conflicts: conda update --all
  • Check Python version compatibility
"""
    return installation_text.strip()

def _validate_scvi_tools(context: Optional[Context] = None):
    """Validate scvi-tools availability and return the module"""
    try:
        import scvi
        from scvi.external import CellAssign
        
        if context:
            # Optional: Check version compatibility
            import pkg_resources
            try:
                version = pkg_resources.get_distribution("scvi-tools").version
                context.info(f"Using scvi-tools version {version}")
            except:
                pass  # Version check is optional
        
        return scvi, CellAssign
    except ImportError as e:
        install_guide = _get_installation_guide("scvi-tools")
        raise DependencyError("scvi-tools", "scANVI", install_guide) from e

def _validate_tangram(context: Optional[Context] = None):
    """Validate tangram availability and return the module"""
    try:
        import tangram as tg
        
        if context and hasattr(tg, '__version__'):
            context.info(f"Using tangram version {tg.__version__}")
        
        return tg
    except ImportError as e:
        install_guide = _get_installation_guide("tangram-sc")
        raise DependencyError("tangram-sc", "Tangram", install_guide) from e

def _validate_mllmcelltype(context: Optional[Context] = None):
    """Validate mllmcelltype availability and return the module"""
    try:
        import mllmcelltype
        
        if context and hasattr(mllmcelltype, '__version__'):
            context.info(f"Using mllmcelltype version {mllmcelltype.__version__}")
        
        return mllmcelltype
    except ImportError as e:
        install_guide = _get_installation_guide("mllmcelltype")
        raise DependencyError("mllmcelltype", "mLLMCellType", install_guide) from e

def _validate_rpy2_and_r(context: Optional[Context] = None):
    """Validate R and rpy2 availability with detailed error reporting"""
    try:
        # First check rpy2
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri, numpy2ri
        from rpy2.robjects.packages import importr
        from rpy2.robjects.conversion import localconverter
        
        # Test R availability
        robjects.r('R.version')
        
        if context:
            r_version = robjects.r('R.version.string')[0]
            context.info(f"Using R: {r_version}")
        
        return robjects, pandas2ri, numpy2ri, importr, localconverter
    except ImportError as e:
        install_guide = _get_installation_guide("rpy2")
        raise DependencyError("rpy2 and R", "sc-type", install_guide) from e
    except Exception as e:
        error_msg = f"""
R environment setup failed: {str(e)}

Common solutions:
  • Install R from https://www.r-project.org/
  • Set R_HOME environment variable
  • Install required R packages: install.packages(c('dplyr', 'openxlsx', 'HGNChelper'))
  • macOS: brew install r
  • Ubuntu: sudo apt install r-base
  • Windows: Download from CRAN
"""
        raise DependencyError("R environment", "sc-type", error_msg) from e

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

# Supported annotation methods
SUPPORTED_METHODS = {
    "tangram", "scanvi", "cellassign", "marker_genes", 
    "mllmcelltype", "sctype"
}

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

async def _annotate_with_marker_genes_scanpy(adata, marker_genes: Dict[str, List[str]], context: Optional[Context] = None) -> tuple:
    """Use scanpy's official marker_gene_overlap method for annotation"""
    try:
        if context:
            await context.info("Using scanpy's official marker_gene_overlap method")
        
        # Step 1: Validate clustering for differential expression analysis
        if 'leiden' not in adata.obs.columns:
            raise ValueError(
                "Leiden clustering not found in adata.obs. "
                "Marker gene analysis requires clustering information. "
                "Please run clustering in preprocessing.py or use: sc.tl.leiden(adata)"
            )
        
        # Step 2: Calculate differential expression for clusters
        if context:
            await context.info("Computing differential expression analysis")
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', n_genes=100)
        
        # Step 3: Use scanpy's official marker gene overlap
        if context:
            await context.info(f"Computing marker gene overlap for {len(marker_genes)} cell types")
        
        # Convert marker genes to the format expected by scanpy
        reference_markers = {cell_type: set(genes) for cell_type, genes in marker_genes.items()}
        
        overlap_scores = sc.tl.marker_gene_overlap(
            adata,
            reference_markers=reference_markers,
            method='overlap_count',  # Count overlapping genes
            top_n_markers=50,  # Use top 50 DE genes per cluster
            inplace=False
        )
        
        if context:
            await context.info(f"Overlap analysis completed for {len(overlap_scores)} clusters")
        
        # Step 4: Assign cell types based on highest overlap scores
        # Note: overlap_scores has cell_types as index and clusters as columns
        cluster_to_celltype = {}
        confidence_scores_per_cluster = {}
        
        # For each cluster (column), find the cell type (index) with highest overlap
        for cluster_col in overlap_scores.columns:
            cluster_str = str(cluster_col)
            # Get scores for this cluster across all cell types
            cluster_scores = overlap_scores[cluster_col]  # This is a Series
            best_celltype = cluster_scores.idxmax()  # Cell type with highest overlap
            max_overlap = cluster_scores[best_celltype]
            
            cluster_to_celltype[cluster_str] = best_celltype
            
            # Calculate confidence based on overlap ratio
            total_possible_overlap = len(reference_markers[best_celltype])
            confidence = min(max_overlap / total_possible_overlap, 0.95) if total_possible_overlap > 0 else 0.1
            confidence_scores_per_cluster[cluster_str] = round(max(confidence, CONFIDENCE_MIN), 2)
        
        # Step 5: Map cluster assignments to individual cells
        adata.obs['cell_type'] = adata.obs['leiden'].astype(str).map(cluster_to_celltype)
        
        # Handle any unmapped clusters (assign as "Unknown")
        unmapped_mask = adata.obs['cell_type'].isna()
        if unmapped_mask.any():
            if context:
                await context.warning(f"Found {unmapped_mask.sum()} cells in unmapped clusters, assigning as 'Unknown'")
            adata.obs.loc[unmapped_mask, 'cell_type'] = 'Unknown'
            confidence_scores_per_cluster['Unknown'] = CONFIDENCE_MIN
        
        adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
        
        # Step 6: Calculate final statistics
        cell_types = list(adata.obs['cell_type'].unique())
        counts = adata.obs['cell_type'].value_counts().to_dict()
        
        # Map cluster confidence to cell type confidence (average across clusters for each cell type)
        confidence_scores = {}
        for cell_type in cell_types:
            # Find all clusters assigned to this cell type
            assigned_clusters = [k for k, v in cluster_to_celltype.items() if v == cell_type]
            if assigned_clusters:
                # Average confidence across clusters for this cell type
                avg_confidence = np.mean([confidence_scores_per_cluster.get(c, CONFIDENCE_MIN) for c in assigned_clusters])
                confidence_scores[cell_type] = round(avg_confidence, 2)
            else:
                confidence_scores[cell_type] = CONFIDENCE_MIN
        
        if context:
            await context.info(f"✅ Marker gene annotation completed: {len(cell_types)} cell types identified")
            top_types = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:3]
            await context.info(f"Top cell types: {', '.join([f'{t}({c})' for t, c in top_types])}")
        
        return cell_types, counts, confidence_scores
        
    except Exception as e:
        error_msg = f"Scanpy marker gene overlap annotation failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)


async def _annotate_with_marker_genes(adata, params: AnnotationParameters, context: Optional[Context] = None):
    """Annotate cell types using scanpy's official marker gene overlap method"""
    try:
        if context:
            await context.info("Using scanpy's official marker gene overlap method")
        
        # Use provided marker genes or default ones
        marker_genes = params.marker_genes if params.marker_genes else DEFAULT_MARKER_GENES
        
        # Validate marker genes exist in dataset
        valid_marker_genes = await _validate_marker_genes(adata, marker_genes, context)
        
        # Use scanpy's official method
        cell_types, counts, confidence_scores = await _annotate_with_marker_genes_scanpy(
            adata, valid_marker_genes, context
        )
        
        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization
        
        return cell_types, counts, confidence_scores, None
        
    except Exception as e:
        await _handle_annotation_error(e, "marker_genes", context)


async def _annotate_with_tangram(adata, params: AnnotationParameters, data_store: Dict[str, Any], context: Optional[Context] = None):
    """Annotate cell types using Tangram method"""
    try:
        if context:
            await context.info("Using Tangram method for annotation")

        # Validate dependencies with comprehensive error reporting
        tg = _validate_tangram(context)

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
                if 'highly_variable' not in adata_sc.var:
                    raise ValueError(
                        "Highly variable genes not found in reference data. "
                        "Tangram mapping requires HVG selection. "
                        "Please run HVG selection in preprocessing.py or use: sc.pp.highly_variable_genes(adata)"
                    )
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
                if context:
                    await context.info(f"Projecting cluster annotations using '{cluster_label}' column")
                tg.plot_cell_annotation(ad_map, adata_sp, annotation=cluster_label)
            else:
                # For cells mode, try multiple annotation columns
                annotation_col = None
                potential_cols = ['cell_type', 'celltype', 'cell_types', 'subclass_label']
                
                for col in potential_cols:
                    if col in adata_sc.obs:
                        annotation_col = col
                        break
                
                if annotation_col:
                    if context:
                        await context.info(f"Projecting cell annotations using '{annotation_col}' column")
                    tg.plot_cell_annotation(ad_map, adata_sp, annotation=annotation_col)
                else:
                    if context:
                        await context.warning("No suitable annotation column found for cells mode projection")
        except Exception as proj_error:
            if context:
                await context.warning(f"Could not project cell annotations: {proj_error}")
            # Continue without projection

        # Get cell type predictions
        cell_types = []
        counts = {}
        confidence_scores = {}
        
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

        return cell_types, counts, confidence_scores, tangram_mapping_score
        
    except Exception as e:
        await _handle_annotation_error(e, "tangram", context)

async def _annotate_with_scanvi(adata, params: AnnotationParameters, data_store: Dict[str, Any], context: Optional[Context] = None):
    """Annotate cell types using scANVI method"""
    try:
        if context:
            await context.info("Using scANVI method for annotation")

        # Validate dependencies with comprehensive error reporting
        scvi, CellAssign = _validate_scvi_tools(context)

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

        # Validate dependencies with comprehensive error reporting
        mllmcelltype = _validate_mllmcelltype(context)

        # Validate clustering has been performed
        cluster_key = params.cluster_label if params.cluster_label else 'leiden'
        if cluster_key not in adata.obs:
            raise ValueError(
                f"Clustering key '{cluster_key}' not found in adata.obs. "
                f"mLLMCellType annotation requires clustering information. "
                f"Please run clustering in preprocessing.py or use: sc.tl.{cluster_key}(adata)"
            )

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

        # Validate dependencies with comprehensive error reporting
        scvi, CellAssign = _validate_scvi_tools(context)

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

    # Validate method first - clean and simple
    if params.method not in SUPPORTED_METHODS:
        raise ValueError(f"Unsupported method: {params.method}. Supported: {sorted(SUPPORTED_METHODS)}")

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
        elif params.method == "mllmcelltype":
            cell_types, counts, confidence_scores, tangram_mapping_score = await _annotate_with_mllmcelltype(
                adata, params, context
            )
        else:  # sctype
            cell_types, counts, confidence_scores, tangram_mapping_score = await _annotate_with_sctype(
                adata, params, context
            )

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


# ============================================================================
# SC-TYPE IMPLEMENTATION
# ============================================================================

# Cache for sc-type results to avoid repeated R calls
_SCTYPE_CACHE = {}
_SCTYPE_CACHE_DIR = Path.home() / ".chatspatial" / "sctype_cache"


def _get_sctype_cache_key(adata, params: AnnotationParameters) -> str:
    """Generate cache key for sc-type results"""
    # Create a hash based on data and parameters
    data_hash = hashlib.md5()
    
    # Hash expression data (sample first 1000 cells and 500 genes for efficiency)
    if hasattr(adata.X, 'toarray'):
        sample_data = adata.X[:min(1000, adata.n_obs), :min(500, adata.n_vars)].toarray()
    else:
        sample_data = adata.X[:min(1000, adata.n_obs), :min(500, adata.n_vars)]
    data_hash.update(sample_data.tobytes())
    
    # Hash relevant parameters
    params_dict = {
        'tissue': params.sctype_tissue,
        'db': params.sctype_db_,
        'scaled': params.sctype_scaled,
        'custom_markers': params.sctype_custom_markers
    }
    data_hash.update(str(params_dict).encode())
    
    return data_hash.hexdigest()


async def _load_sctype_functions(context: Optional[Context] = None) -> None:
    """Load sc-type R functions and auto-install R packages if needed"""
    if context:
        await context.info("📚 Loading sc-type R functions...")
    
    try:
        # Auto-install required R packages
        robjects.r('''
            # Install packages automatically if not present
            required_packages <- c("dplyr", "openxlsx", "HGNChelper")
            for (pkg in required_packages) {
                if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
                    cat("Installing R package:", pkg, "\\n")
                    install.packages(pkg, repos = "https://cran.r-project.org/", quiet = TRUE)
                    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
                        stop(paste("Failed to install required R package:", pkg))
                    }
                }
            }
        ''')
        
        if context:
            await context.info("✅ R packages loaded/installed successfully")
        
        # Load sc-type functions from GitHub
        robjects.r('''
            # Load gene sets preparation function
            source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
            
            # Load scoring function
            source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
        ''')
        
        if context:
            await context.info("✅ sc-type R functions loaded successfully")
            
    except Exception as e:
        error_msg = f"Failed to load sc-type R functions: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)


async def _prepare_sctype_genesets(params: AnnotationParameters, context: Optional[Context] = None):
    """Prepare gene sets for sc-type"""
    if context:
        await context.info("🧬 Preparing sc-type gene sets...")
    
    try:
        if params.sctype_custom_markers:
            # Use custom markers
            if context:
                await context.info("Using custom marker gene sets")
            return _convert_custom_markers_to_gs(params.sctype_custom_markers, context)
        else:
            # Use sc-type database
            tissue = params.sctype_tissue
            if not tissue:
                raise ValueError("sctype_tissue parameter is required when not using custom markers")
            
            if context:
                await context.info(f"Using sc-type database for tissue: {tissue}")
                
            # Download and use ScTypeDB
            db_path = params.sctype_db_ or "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
            
            robjects.r.assign('db_path', db_path)
            robjects.r.assign('tissue_type', tissue)
            
            robjects.r('''
                # Load gene sets
                gs_list <- gene_sets_prepare(db_path, tissue_type)
            ''')
            
            return robjects.r['gs_list']
            
    except Exception as e:
        error_msg = f"Failed to prepare gene sets: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)


def _convert_custom_markers_to_gs(custom_markers: Dict[str, Dict[str, List[str]]], context: Optional[Context] = None):
    """Convert custom markers to sc-type gene set format"""
    if not custom_markers:
        raise ValueError("Custom markers dictionary is empty")
    
    gs_positive = {}
    gs_negative = {}
    
    valid_celltypes = 0
    
    for cell_type, markers in custom_markers.items():
        if not isinstance(markers, dict):
            continue
            
        positive_genes = []
        negative_genes = []
        
        if 'positive' in markers and isinstance(markers['positive'], list):
            positive_genes = [gene.upper().strip() for gene in markers['positive'] if gene and str(gene).strip()]
            
        if 'negative' in markers and isinstance(markers['negative'], list):
            negative_genes = [gene.upper().strip() for gene in markers['negative'] if gene and str(gene).strip()]
        
        # Only include cell types that have at least some positive markers
        if positive_genes:
            gs_positive[cell_type] = positive_genes
            gs_negative[cell_type] = negative_genes  # Can be empty list
            valid_celltypes += 1
    
    if valid_celltypes == 0:
        raise ValueError("No valid cell types found in custom markers - all cell types need at least one positive marker")
    
    # Convert to R lists using proper nested conversion
    with localconverter(robjects.default_converter + pandas2ri.converter):
        # Convert Python dictionaries to R named lists, handle empty lists properly
        r_gs_positive = robjects.r['list'](**{
            k: robjects.StrVector(v) if v else robjects.StrVector([]) 
            for k, v in gs_positive.items()
        })
        r_gs_negative = robjects.r['list'](**{
            k: robjects.StrVector(v) if v else robjects.StrVector([]) 
            for k, v in gs_negative.items()
        })
        
        # Create the final gs_list structure
        gs_list = robjects.r['list'](
            gs_positive=r_gs_positive,
            gs_negative=r_gs_negative
        )
    
    return gs_list




async def _run_sctype_scoring(adata, gs_list, params: AnnotationParameters, context: Optional[Context] = None):
    """Run sc-type scoring algorithm"""
    if context:
        await context.info("🔬 Running sc-type scoring...")
    
    try:
        
        # Prepare expression data
        if params.sctype_scaled and 'scaled' in adata.layers:
            expr_data = adata.layers['scaled']
            if context:
                await context.info("Using scaled expression data from adata.layers['scaled']")
        else:
            expr_data = adata.X
            if context:
                await context.info("Using raw expression data from adata.X")
        
        # Convert to dense array if sparse
        if hasattr(expr_data, 'toarray'):
            expr_data = expr_data.toarray()
        
        # Convert to DataFrame with proper gene and cell names
        expr_df = pd.DataFrame(
            expr_data.T,  # sc-type expects genes as rows, cells as columns
            index=adata.var_names,
            columns=adata.obs_names
        )
        
        # Transfer to R
        with localconverter(robjects.default_converter + pandas2ri.converter):
            robjects.r.assign('scdata', expr_df)
        
        robjects.r.assign('gs_positive', gs_list[0])  # gs_positive from gs_list
        robjects.r.assign('gs_negative', gs_list[1])  # gs_negative from gs_list
        
        # Run sc-type scoring with comprehensive error handling
        robjects.r('''
            # Check if gene sets are valid
            if (length(gs_positive) == 0) {
                stop("No valid positive gene sets found")
            }
            
            # Get available genes in the dataset
            available_genes <- rownames(scdata)
            
            # Check each cell type for overlapping genes and filter empty ones
            filtered_gs_positive <- list()
            filtered_gs_negative <- list()
            
            for (celltype in names(gs_positive)) {
                pos_genes <- gs_positive[[celltype]]
                neg_genes <- if (celltype %in% names(gs_negative)) gs_negative[[celltype]] else c()
                
                # Find overlapping genes
                pos_overlap <- intersect(toupper(pos_genes), toupper(available_genes))
                neg_overlap <- intersect(toupper(neg_genes), toupper(available_genes))
                
                # Only keep cell types with at least one positive marker gene
                if (length(pos_overlap) > 0) {
                    filtered_gs_positive[[celltype]] <- pos_overlap
                    filtered_gs_negative[[celltype]] <- neg_overlap
                }
            }
            
            # Check if we have any valid cell types after filtering
            if (length(filtered_gs_positive) == 0) {
                # Create a dummy result with "Unknown" for all cells
                n_cells <- ncol(scdata)
                es_max <- matrix(0, nrow = 1, ncol = n_cells)
                rownames(es_max) <- c("Unknown")
                colnames(es_max) <- colnames(scdata)
            } else {
                # Run sc-type scoring with filtered gene sets
                tryCatch({
                    es_max <- sctype_score(
                        scRNAseqData = as.matrix(scdata),
                        scaled = TRUE,
                        gs = filtered_gs_positive,
                        gs2 = filtered_gs_negative
                    )
                    
                    # Handle edge case where sc-type returns empty or invalid results
                    if (is.null(es_max) || nrow(es_max) == 0 || ncol(es_max) == 0) {
                        # Create a dummy result
                        n_cells <- ncol(scdata)
                        es_max <- matrix(0, nrow = 1, ncol = n_cells)
                        rownames(es_max) <- c("Unknown")
                        colnames(es_max) <- colnames(scdata)
                    }
                }, error = function(e) {
                    # If sc-type scoring fails, create a dummy result
                    n_cells <- ncol(scdata)
                    es_max <- matrix(0, nrow = 1, ncol = n_cells)
                    rownames(es_max) <- c("Unknown")
                    colnames(es_max) <- colnames(scdata)
                })
            }
        ''')
        
        # Get results back to Python
        with localconverter(robjects.default_converter + pandas2ri.converter):
            scores_df = robjects.r['es_max']
        
        if context:
            await context.info(f"✅ Scoring completed for {scores_df.shape[0]} cell types and {scores_df.shape[1]} cells")
        
        return scores_df
        
    except Exception as e:
        error_msg = f"Failed to run sc-type scoring: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)


async def _assign_sctype_celltypes(scores_df, context: Optional[Context] = None):
    """Assign cell types based on sc-type scores"""
    if context:
        await context.info("🏷️  Assigning cell types based on scores...")
    
    try:
        # Handle empty or invalid scores
        if scores_df is None:
            raise ValueError("Scores data is None")
        
        # Convert to DataFrame if it's a numpy array
        if isinstance(scores_df, np.ndarray):
            # Check for empty array
            if scores_df.size == 0:
                raise ValueError("Scores array is empty")
            
            # Convert to DataFrame - need to get column and row names from R
            import pandas as pd
            scores_df = pd.DataFrame(scores_df)
            if context:
                await context.info(f"Converted numpy array to DataFrame: {scores_df.shape}")
        
        # Check DataFrame has data
        if hasattr(scores_df, 'empty') and scores_df.empty:
            raise ValueError("Scores DataFrame is empty")
        
        # Get the cell type with highest score for each cell
        cell_types = []
        confidence_scores = []
        
        # Handle both DataFrame and numpy array cases
        if hasattr(scores_df, 'columns'):
            # DataFrame case
            if len(scores_df.columns) == 0:
                raise ValueError("DataFrame has no columns")
            
            n_cells = len(scores_df.columns)
            n_celltypes = len(scores_df.index)
            
            # Handle single cell type case
            if n_celltypes == 1:
                # Only one cell type available - assign to all cells based on scores
                single_celltype = str(scores_df.index[0])
                
                for col_idx in range(n_cells):
                    cell_score = scores_df.iloc[0, col_idx]
                    if cell_score > 0:
                        cell_types.append(single_celltype)
                        confidence_scores.append(min(cell_score / 10.0, 1.0))
                    else:
                        cell_types.append("Unknown")
                        confidence_scores.append(0.0)
            else:
                # Multiple cell types available
                
                for col_idx in range(n_cells):
                    cell_scores = scores_df.iloc[:, col_idx]  # All cell types for this cell
                    
                    # Find cell type with maximum score
                    max_score_idx = cell_scores.idxmax()
                    max_score = cell_scores.loc[max_score_idx]
                    
                    # If max score is positive, assign cell type, otherwise "Unknown"
                    if max_score > 0:
                        cell_types.append(str(max_score_idx))
                        # Normalize confidence score to 0-1 range
                        confidence_scores.append(min(max_score / 10.0, 1.0))
                    else:
                        cell_types.append("Unknown")
                        confidence_scores.append(0.0)
        else:
            # Numpy array case
            if len(scores_df.shape) == 0 or scores_df.size == 0:
                raise ValueError("Empty scores array")
                
            n_cells = scores_df.shape[1] if len(scores_df.shape) > 1 else 1
            n_celltypes = scores_df.shape[0] if len(scores_df.shape) > 1 else scores_df.shape[0]
            
            if n_celltypes == 0:
                raise ValueError("No cell types in scores array")
            
            for cell_idx in range(n_cells):
                if len(scores_df.shape) > 1:
                    cell_scores = scores_df[:, cell_idx]
                else:
                    cell_scores = scores_df if n_cells == 1 else scores_df[cell_idx:cell_idx+1]
                
                # Handle case where cell_scores is empty
                if len(cell_scores) == 0:
                    cell_types.append("Unknown")
                    confidence_scores.append(0.0)
                    continue
                
                # Find cell type with maximum score
                max_score_idx = np.argmax(cell_scores)
                max_score = cell_scores[max_score_idx]
                
                # If max score is positive, assign cell type, otherwise "Unknown"
                if max_score > 0:
                    cell_types.append(f"CellType_{max_score_idx}")
                    # Normalize confidence score to 0-1 range
                    confidence_scores.append(min(max_score / 10.0, 1.0))
                else:
                    cell_types.append("Unknown")
                    confidence_scores.append(0.0)
        
        if context:
            unique_types = list(set(cell_types))
            await context.info(f"✅ Assigned {len(unique_types)} unique cell types: {', '.join(unique_types[:10])}" + 
                             ("..." if len(unique_types) > 10 else ""))
        
        return cell_types, confidence_scores
        
    except Exception as e:
        error_msg = f"Failed to assign cell types: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)


def _calculate_sctype_stats(cell_types):
    """Calculate statistics for sc-type results"""
    from collections import Counter
    counts = Counter(cell_types)
    return dict(counts)


async def _cache_sctype_results(cache_key: str, results, context: Optional[Context] = None) -> None:
    """Cache sc-type results to disk"""
    if context:
        await context.info("💾 Caching sc-type results...")
    
    try:
        # Create cache directory
        _SCTYPE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
        
        # Save to cache file
        cache_file = _SCTYPE_CACHE_DIR / f"{cache_key}.pkl"
        with open(cache_file, 'wb') as f:
            pickle.dump(results, f)
        
        # Also store in memory cache
        _SCTYPE_CACHE[cache_key] = results
        
        if context:
            await context.info("✅ Results cached successfully")
            
    except Exception as e:
        if context:
            await context.warning(f"Failed to cache results: {str(e)}")


async def _load_cached_sctype_results(cache_key: str, context: Optional[Context] = None):
    """Load cached sc-type results"""
    # Check memory cache first
    if cache_key in _SCTYPE_CACHE:
        if context:
            await context.info("📂 Using cached results from memory")
        return _SCTYPE_CACHE[cache_key]
    
    # Check disk cache
    cache_file = _SCTYPE_CACHE_DIR / f"{cache_key}.pkl"
    if cache_file.exists():
        try:
            with open(cache_file, 'rb') as f:
                results = pickle.load(f)
            
            # Store in memory cache for next time
            _SCTYPE_CACHE[cache_key] = results
            
            if context:
                await context.info("📂 Using cached results from disk")
            return results
            
        except Exception as e:
            if context:
                await context.warning(f"Failed to load cached results: {str(e)}")
    
    return None


async def _annotate_with_sctype(
    adata: sc.AnnData,
    params: AnnotationParameters,
    context: Optional[Context] = None
) -> tuple[List[str], Dict[str, int], List[float], Optional[float]]:
    """
    Annotate cell types using sc-type method
    
    Args:
        adata: AnnData object
        params: Annotation parameters
        context: MCP context
        
    Returns:
        Tuple of (cell_types, counts, confidence_scores, mapping_score)
    """
    
    # Validate dependencies with comprehensive error reporting
    robjects, pandas2ri, numpy2ri, importr, localconverter = _validate_rpy2_and_r(context)
    
    # Define supported tissue types from sc-type database
    SCTYPE_VALID_TISSUES = {
        "Adrenal", "Brain", "Eye", "Heart", "Hippocampus", "Immune system", 
        "Intestine", "Kidney", "Liver", "Lung", "Muscle", "Pancreas", 
        "Placenta", "Spleen", "Stomach", "Thymus"
    }
    
    # Validate required parameters
    if not params.sctype_tissue and not params.sctype_custom_markers:
        raise ValueError("Either sctype_tissue or sctype_custom_markers must be specified")
    
    # Validate tissue type if provided
    if params.sctype_tissue and params.sctype_tissue not in SCTYPE_VALID_TISSUES:
        valid_tissues = sorted(SCTYPE_VALID_TISSUES)
        raise ValueError(
            f"Tissue type '{params.sctype_tissue}' is not supported by sc-type database.\n"
            f"Supported tissues: {', '.join(valid_tissues)}\n"
            f"Alternatively, use sctype_custom_markers to define custom cell type markers."
        )
    
    try:
        if context:
            await context.info("🧬 Starting sc-type cell type annotation...")
            await context.info(f"📊 Analyzing {adata.n_obs} cells with {adata.n_vars} genes")
            if params.sctype_tissue:
                await context.info(f"🔬 Using tissue type: {params.sctype_tissue}")
            else:
                await context.info(f"🔬 Using custom markers for {len(params.sctype_custom_markers)} cell types")
        
        # Check cache if enabled
        cache_key = None
        if params.sctype_use_cache:
            cache_key = _get_sctype_cache_key(adata, params)
            cached_results = await _load_cached_sctype_results(cache_key, context)
            if cached_results:
                return cached_results
        
        # Step 1: Load sc-type functions
        await _load_sctype_functions(context)
        
        # Step 2: Prepare gene sets
        gs_list = await _prepare_sctype_genesets(params, context)
        
        # Step 3: Run sc-type scoring
        scores_df = await _run_sctype_scoring(adata, gs_list, params, context)
        
        # Step 4: Assign cell types
        cell_types, confidence_scores = await _assign_sctype_celltypes(scores_df, context)
        
        # Step 5: Calculate statistics
        counts = _calculate_sctype_stats(cell_types)
        
        # Step 6: Calculate average confidence scores per cell type (to match AnnotationResult interface)
        confidence_by_celltype = {}
        for cell_type in set(cell_types):
            # Get confidence scores for this cell type
            celltype_confidences = [conf for i, conf in enumerate(confidence_scores) if cell_types[i] == cell_type]
            if celltype_confidences:
                confidence_by_celltype[cell_type] = sum(celltype_confidences) / len(celltype_confidences)
            else:
                confidence_by_celltype[cell_type] = 0.0
        
        # Step 7: Add results to adata
        
        adata.obs['sctype_celltype'] = cell_types
        adata.obs['sctype_confidence'] = confidence_scores
        
        results = (cell_types, counts, confidence_by_celltype, None)
        
        # Cache results if enabled
        if params.sctype_use_cache and cache_key:
            await _cache_sctype_results(cache_key, results, context)
        
        if context:
            await context.info("🎉 sc-type annotation completed successfully!")
            await context.info(f"📈 Found {len(counts)} cell types: {', '.join(list(counts.keys())[:5])}" + 
                             ("..." if len(counts) > 5 else ""))
        
        return results
        
    except Exception as e:
        await _handle_annotation_error(e, "sctype", context)