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


# Import optional dependencies with fallbacks
try:
    import commot as ct
    COMMOT_AVAILABLE = True
except ImportError:
    COMMOT_AVAILABLE = False

try:
    import spatialdm as sdm
    SPATIALDM_AVAILABLE = True
except ImportError:
    SPATIALDM_AVAILABLE = False


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
        
        # Analyze cell communication based on method
        if params.method == "commot":
            result_data = await _analyze_communication_commot(adata, params, context)
        elif params.method == "spatialdm":
            result_data = await _analyze_communication_spatialdm(adata, params, context)
        elif params.method == "cellphonedb":
            result_data = await _analyze_communication_cellphonedb(adata, params, context)
        else:
            raise ValueError(f"Unsupported method: {params.method}")
        
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
            database=params.database,
            n_lr_pairs=result_data["n_lr_pairs"],
            n_significant_pairs=result_data["n_significant_pairs"],
            global_results_key=result_data.get("global_results_key"),
            top_lr_pairs=result_data.get("top_lr_pairs", []),
            local_analysis_performed=result_data.get("local_analysis_performed", False),
            local_results_key=result_data.get("local_results_key"),
            communication_matrices_key=result_data.get("communication_matrices_key"),
            commot_sender_key=result_data.get("commot_sender_key"),
            commot_receiver_key=result_data.get("commot_receiver_key"),
            spatialdm_selected_spots_key=result_data.get("spatialdm_selected_spots_key"),
            spatialdm_weight_matrix_key=result_data.get("spatialdm_weight_matrix_key"),
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


async def _analyze_communication_commot(
    adata: Any,
    params: CellCommunicationParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze cell communication using COMMOT"""
    if not COMMOT_AVAILABLE:
        raise ImportError("COMMOT is not installed. Please install it with: pip install commot")
    
    if context:
        await context.info("Running COMMOT for cell communication analysis...")
    
    try:
        # Get ligand-receptor database
        if context:
            await context.info(f"Loading {params.database} ligand-receptor database for {params.species}...")
        
        if params.database == "user" and params.custom_lr_pairs:
            # Create custom LR dataframe
            lr_data = []
            for i, (ligand, receptor) in enumerate(params.custom_lr_pairs):
                lr_data.append([ligand, receptor, f"custom_pathway_{i}"])
            df_ligrec = pd.DataFrame(lr_data, columns=['ligand', 'receptor', 'pathway'])
        else:
            # Use built-in database
            database_name = "CellChat" if params.database == "cellchat" else "CellPhoneDB"
            df_ligrec = ct.pp.ligand_receptor_database(
                database=database_name, 
                species=params.species
            )
        
        if context:
            await context.info(f"Found {len(df_ligrec)} ligand-receptor pairs in database")
        
        # Run COMMOT spatial communication analysis
        if context:
            await context.info("Computing spatial communication networks...")
        
        ct.tl.spatial_communication(
            adata,
            database_name=params.database,
            df_ligrec=df_ligrec,
            dis_thr=params.commot_dis_thr,
            heteromeric=params.commot_heteromeric
        )
        
        # Get communication results
        sender_key = f'commot-{params.database}-sum-sender'
        receiver_key = f'commot-{params.database}-sum-receiver'
        
        # Calculate significance if requested
        n_significant_pairs = 0
        top_lr_pairs = []
        
        if sender_key in adata.obsm and receiver_key in adata.obsm:
            sender_df = pd.DataFrame(adata.obsm[sender_key], index=adata.obs.index)
            receiver_df = pd.DataFrame(adata.obsm[receiver_key], index=adata.obs.index)
            
            # Get top pairs based on total communication strength
            total_communication = sender_df.sum() + receiver_df.sum()
            top_lr_pairs = total_communication.nlargest(params.plot_top_pairs).index.tolist()
            n_significant_pairs = len([x for x in total_communication if x > 0])
        
        # Identify communication patterns if requested
        patterns_identified = False
        n_patterns = None
        patterns_key = None
        
        if params.identify_communication_patterns:
            if context:
                await context.info("Identifying communication patterns...")
            try:
                # Use clustering on communication profiles
                from sklearn.cluster import KMeans
                
                if sender_key in adata.obsm:
                    comm_data = adata.obsm[sender_key] + adata.obsm[receiver_key]
                    # Filter out spots with no communication
                    active_spots = comm_data.sum(axis=1) > 0
                    
                    if active_spots.sum() > 10:  # Need enough spots for clustering
                        n_patterns = min(5, active_spots.sum() // 10)  # Adaptive number of patterns
                        kmeans = KMeans(n_clusters=n_patterns, random_state=42)
                        patterns = np.full(len(adata), -1)  # -1 for non-communicating spots
                        patterns[active_spots] = kmeans.fit_predict(comm_data[active_spots])
                        
                        patterns_key = f"commot_{params.database}_patterns"
                        adata.obs[patterns_key] = patterns.astype(str)
                        adata.obs[patterns_key] = adata.obs[patterns_key].astype('category')
                        patterns_identified = True
                        
                        if context:
                            await context.info(f"Identified {n_patterns} communication patterns")
            except Exception as e:
                if context:
                    await context.warning(f"Pattern identification failed: {str(e)}")
        
        statistics = {
            "method": "commot",
            "database": params.database,
            "species": params.species,
            "distance_threshold": params.commot_dis_thr,
            "heteromeric": params.commot_heteromeric,
            "n_lr_pairs_tested": len(df_ligrec),
            "n_communicating_spots": int((adata.obsm[sender_key].sum(axis=1) > 0).sum()) if sender_key in adata.obsm else 0
        }
        
        return {
            "n_lr_pairs": len(df_ligrec),
            "n_significant_pairs": n_significant_pairs,
            "top_lr_pairs": top_lr_pairs,
            "commot_sender_key": sender_key,
            "commot_receiver_key": receiver_key,
            "patterns_identified": patterns_identified,
            "n_patterns": n_patterns,
            "patterns_key": patterns_key,
            "statistics": statistics
        }
        
    except Exception as e:
        raise RuntimeError(f"COMMOT failed: {str(e)}")


async def _analyze_communication_spatialdm(
    adata: Any,
    params: CellCommunicationParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze cell communication using SpatialDM"""
    if not SPATIALDM_AVAILABLE:
        raise ImportError("SpatialDM is not installed. Please install it with: pip install SpatialDM")
    
    if context:
        await context.info("Running SpatialDM for cell communication analysis...")
    
    try:
        # Compute weight matrix
        if context:
            await context.info("Computing spatial weight matrix...")
        
        sdm.weight_matrix(
            adata,
            l=params.spatialdm_l,
            cutoff=params.spatialdm_cutoff,
            n_neighbors=params.spatialdm_n_neighbors,
            single_cell=False
        )
        
        # Extract ligand-receptor pairs
        if context:
            await context.info(f"Extracting ligand-receptor pairs for {params.species}...")
        
        sdm.extract_lr(adata, params.species, min_cell=params.min_cells)
        
        n_lr_pairs = adata.uns.get('num_pairs', 0)
        if context:
            await context.info(f"Found {n_lr_pairs} valid ligand-receptor pairs")
        
        if n_lr_pairs == 0:
            raise ValueError("No valid ligand-receptor pairs found")
        
        # Global analysis
        global_results_key = None
        if params.perform_global_analysis:
            if context:
                await context.info("Performing global LR pair selection...")
            
            sdm.spatialdm_global(
                adata,
                n_perm=params.spatialdm_n_permutations,
                method=params.spatialdm_method
            )
            
            sdm.sig_pairs(
                adata,
                method=params.spatialdm_method.split('_')[0] if '_' in params.spatialdm_method else params.spatialdm_method,
                fdr=params.spatialdm_fdr,
                threshold=params.spatialdm_threshold
            )
            
            global_results_key = "global_res"
        
        # Local analysis
        local_analysis_performed = False
        local_results_key = None
        selected_spots_key = None
        
        if params.perform_local_analysis and global_results_key:
            if context:
                await context.info("Performing local spot analysis...")
            
            try:
                sdm.spatialdm_local(
                    adata,
                    n_perm=params.spatialdm_n_permutations,
                    method=params.spatialdm_method
                )
                
                sdm.sig_spots(
                    adata,
                    method=params.spatialdm_method.split('_')[0] if '_' in params.spatialdm_method else params.spatialdm_method,
                    fdr=params.spatialdm_fdr,
                    threshold=params.spatialdm_threshold
                )
                
                local_analysis_performed = True
                local_results_key = "local_stat"
                selected_spots_key = "selected_spots"
                
            except Exception as e:
                if context:
                    await context.warning(f"Local analysis failed: {str(e)}")
        
        # Get significant pairs and top pairs
        n_significant_pairs = 0
        top_lr_pairs = []
        
        if global_results_key and global_results_key in adata.uns:
            global_res = adata.uns[global_results_key]
            if 'selected' in global_res.columns:
                n_significant_pairs = int(global_res['selected'].sum())
                
                # Get top pairs based on significance
                if params.spatialdm_method in ['z-score', 'both'] and 'z_pval' in global_res.columns:
                    top_pairs_df = global_res[global_res['selected']].nsmallest(params.plot_top_pairs, 'z_pval')
                elif 'perm_pval' in global_res.columns:
                    top_pairs_df = global_res[global_res['selected']].nsmallest(params.plot_top_pairs, 'perm_pval')
                else:
                    top_pairs_df = global_res[global_res['selected']].head(params.plot_top_pairs)
                
                top_lr_pairs = top_pairs_df.index.tolist()
        
        statistics = {
            "method": "spatialdm",
            "species": params.species,
            "n_lr_pairs_tested": n_lr_pairs,
            "weight_cutoff": params.spatialdm_cutoff,
            "n_permutations": params.spatialdm_n_permutations,
            "statistical_method": params.spatialdm_method,
            "fdr_correction": params.spatialdm_fdr,
            "significance_threshold": params.spatialdm_threshold,
            "global_analysis": params.perform_global_analysis,
            "local_analysis": local_analysis_performed
        }
        
        return {
            "n_lr_pairs": n_lr_pairs,
            "n_significant_pairs": n_significant_pairs,
            "global_results_key": global_results_key,
            "top_lr_pairs": top_lr_pairs,
            "local_analysis_performed": local_analysis_performed,
            "local_results_key": local_results_key,
            "spatialdm_selected_spots_key": selected_spots_key,
            "spatialdm_weight_matrix_key": "weight",
            "statistics": statistics
        }
        
    except Exception as e:
        raise RuntimeError(f"SpatialDM failed: {str(e)}")


async def _analyze_communication_cellphonedb(
    adata: Any,
    params: CellCommunicationParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze cell communication using CellPhoneDB approach"""
    if context:
        await context.info("CellPhoneDB method is not yet implemented")
    
    # Placeholder implementation
    # In a full implementation, you would integrate CellPhoneDB here
    raise NotImplementedError("CellPhoneDB method is not yet implemented")


def _create_communication_visualizations(
    adata: Any,
    result_data: Dict[str, Any],
    params: CellCommunicationParameters
) -> Tuple[Optional[Image], Optional[Image]]:
    """Create visualizations for cell communication analysis"""
    try:
        visualization = None
        network_visualization = None
        
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
        
        # Create main visualization based on method
        if params.method == "commot":
            visualization = _create_commot_visualization(adata, coords, result_data, params)
        elif params.method == "spatialdm":
            visualization = _create_spatialdm_visualization(adata, coords, result_data, params)
        
        # Create network visualization if patterns were identified
        if result_data.get("patterns_identified") and result_data.get("patterns_key"):
            network_visualization = _create_network_visualization(adata, coords, result_data, params)
        
        return visualization, network_visualization
        
    except Exception as e:
        return create_placeholder_image(f"Visualization failed: {str(e)}"), None


def _create_commot_visualization(
    adata: Any,
    coords: np.ndarray,
    result_data: Dict[str, Any],
    params: CellCommunicationParameters
) -> Image:
    """Create COMMOT-specific visualization"""
    try:
        sender_key = result_data.get("commot_sender_key")
        receiver_key = result_data.get("commot_receiver_key")
        
        if not sender_key or sender_key not in adata.obsm:
            return create_placeholder_image("No COMMOT communication data found")
        
        # Get top LR pairs
        top_pairs = result_data.get("top_lr_pairs", [])[:4]  # Show top 4
        
        if not top_pairs:
            return create_placeholder_image("No significant communication pairs found")
        
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
        
        sender_data = adata.obsm[sender_key]
        receiver_data = adata.obsm[receiver_key]
        
        for i, pair in enumerate(top_pairs):
            ax = axes[i]
            
            # Get communication strength for this pair
            if isinstance(sender_data, pd.DataFrame):
                if pair in sender_data.columns:
                    comm_strength = sender_data[pair] + receiver_data[pair]
                else:
                    comm_strength = sender_data.iloc[:, i] + receiver_data.iloc[:, i] if i < sender_data.shape[1] else np.zeros(len(adata))
            else:
                comm_strength = sender_data[:, i] + receiver_data[:, i] if i < sender_data.shape[1] else np.zeros(len(adata))
            
            # Create scatter plot
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=comm_strength,
                cmap='Reds',
                s=20,
                alpha=0.8
            )
            
            ax.set_title(f'{pair}', fontsize=12)
            ax.set_xlabel('Spatial X')
            ax.set_ylabel('Spatial Y')
            ax.invert_yaxis()
            ax.set_aspect('equal')
            
            # Add colorbar
            plt.colorbar(scatter, ax=ax, shrink=0.8, label='Communication Strength')
        
        # Hide unused subplots
        for i in range(n_pairs, len(axes)):
            axes[i].set_visible(False)
        
        plt.suptitle(f'Top {n_pairs} Cell Communication Pairs (COMMOT)', fontsize=14)
        plt.tight_layout()
        
        return fig_to_image(fig, dpi=params.image_dpi, format=params.image_format)
        
    except Exception as e:
        return create_placeholder_image(f"COMMOT visualization failed: {str(e)}")


def _create_spatialdm_visualization(
    adata: Any,
    coords: np.ndarray,
    result_data: Dict[str, Any],
    params: CellCommunicationParameters
) -> Image:
    """Create SpatialDM-specific visualization"""
    try:
        selected_spots_key = result_data.get("spatialdm_selected_spots_key")
        
        if not selected_spots_key or selected_spots_key not in adata.uns:
            return create_placeholder_image("No SpatialDM selected spots found")
        
        # Get top LR pairs
        top_pairs = result_data.get("top_lr_pairs", [])[:4]  # Show top 4
        
        if not top_pairs:
            return create_placeholder_image("No significant communication pairs found")
        
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
        
        selected_spots = adata.uns[selected_spots_key]
        
        for i, pair in enumerate(top_pairs):
            ax = axes[i]
            
            # Get selected spots for this pair
            if pair in selected_spots.index:
                is_selected = selected_spots.loc[pair].values.astype(bool)
            else:
                is_selected = np.zeros(len(adata), dtype=bool)
            
            # Create scatter plot
            # Non-selected spots in gray
            ax.scatter(
                coords[~is_selected, 0],
                coords[~is_selected, 1],
                c='lightgray',
                s=15,
                alpha=0.5,
                label='Non-communicating'
            )
            
            # Selected spots in red
            if is_selected.sum() > 0:
                ax.scatter(
                    coords[is_selected, 0],
                    coords[is_selected, 1],
                    c='red',
                    s=25,
                    alpha=0.8,
                    label='Communicating'
                )
            
            ax.set_title(f'{pair}\n({is_selected.sum()} spots)', fontsize=12)
            ax.set_xlabel('Spatial X')
            ax.set_ylabel('Spatial Y')
            ax.invert_yaxis()
            ax.set_aspect('equal')
            ax.legend()
        
        # Hide unused subplots
        for i in range(n_pairs, len(axes)):
            axes[i].set_visible(False)
        
        plt.suptitle(f'Top {n_pairs} Cell Communication Pairs (SpatialDM)', fontsize=14)
        plt.tight_layout()
        
        return fig_to_image(fig, dpi=params.image_dpi, format=params.image_format)
        
    except Exception as e:
        return create_placeholder_image(f"SpatialDM visualization failed: {str(e)}")


def _create_network_visualization(
    adata: Any,
    coords: np.ndarray,
    result_data: Dict[str, Any],
    params: CellCommunicationParameters
) -> Image:
    """Create communication network visualization"""
    try:
        patterns_key = result_data.get("patterns_key")
        
        if not patterns_key or patterns_key not in adata.obs:
            return create_placeholder_image("No communication patterns found")
        
        patterns = adata.obs[patterns_key]
        unique_patterns = [p for p in patterns.unique() if p != '-1']  # Exclude non-communicating spots
        
        if len(unique_patterns) == 0:
            return create_placeholder_image("No communication patterns identified")
        
        # Create visualization
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot non-communicating spots in gray
        non_comm_mask = patterns == '-1'
        if non_comm_mask.sum() > 0:
            ax.scatter(
                coords[non_comm_mask, 0],
                coords[non_comm_mask, 1],
                c='lightgray',
                s=15,
                alpha=0.3,
                label='Non-communicating'
            )
        
        # Plot each communication pattern
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_patterns)))
        
        for i, pattern in enumerate(unique_patterns):
            pattern_mask = patterns == pattern
            if pattern_mask.sum() > 0:
                ax.scatter(
                    coords[pattern_mask, 0],
                    coords[pattern_mask, 1],
                    c=[colors[i]],
                    s=30,
                    alpha=0.8,
                    label=f'Pattern {pattern} ({pattern_mask.sum()} spots)'
                )
        
        ax.set_xlabel('Spatial X')
        ax.set_ylabel('Spatial Y')
        ax.set_title(f'Communication Patterns ({len(unique_patterns)} patterns)')
        ax.invert_yaxis()
        ax.set_aspect('equal')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        
        return fig_to_image(fig, dpi=params.image_dpi, format=params.image_format)
        
    except Exception as e:
        return create_placeholder_image(f"Network visualization failed: {str(e)}")
