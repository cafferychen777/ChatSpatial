"""
Cell-cell communication analysis tools for spatial transcriptomics data.
"""

from typing import Dict, Any, Optional, List, Tuple, Union
import numpy as np
import pandas as pd
import scanpy as sc
import traceback
import warnings
from scipy import sparse

from mcp.server.fastmcp import Context

from ..models.data import CellCommunicationParameters
from ..models.analysis import CellCommunicationResult


# Import LIANA+ for cell communication analysis
try:
    import liana as li
    LIANA_AVAILABLE = True
except ImportError:
    LIANA_AVAILABLE = False

# Import CellPhoneDB for cell communication analysis
try:
    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
    from cellphonedb.src.core.methods import cpdb_degs_analysis_method
    CELLPHONEDB_AVAILABLE = True
except ImportError:
    CELLPHONEDB_AVAILABLE = False


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
        # Comprehensive data validation
        if context:
            await context.info("Performing comprehensive data quality validation...")
        
        validation_result = await _comprehensive_data_validation(adata, context)
        if not validation_result["passed"]:
            error_msg = f"Data validation failed: {validation_result['error_message']}"
            if validation_result.get('suggestions'):
                error_msg += f"\n\nSuggestions to fix:\n{validation_result['suggestions']}"
            raise ValueError(error_msg)
        
        # Prepare data for cell communication analysis
        if context:
            await context.info("Data validation passed. Preparing for cell communication analysis...")
        
        # Validate input data is preprocessed according to cell communication best practices
        # Based on single-cell best practices (2024) and LIANA+ recommendations
        max_val = adata.X.max() if hasattr(adata.X, 'max') else np.max(adata.X)
        min_val = adata.X.min() if hasattr(adata.X, 'min') else np.min(adata.X)
        
        # Check for raw count data (typically high integer values without log transformation)
        if min_val >= 0 and max_val > 100:
            raise ValueError(
                "Data appears to be raw counts. Cell communication analysis requires: "
                "1) Count normalization (e.g., to 10,000 counts per cell), "
                "2) Log1p transformation to stabilize variance. "
                "According to single-cell best practices, this is essential for CCC inference. "
                "Please use preprocessing.py or run: sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)"
            )
        
        # Validate log transformation was applied (recommended range after log1p is typically 0-15)
        if min_val >= 0 and max_val > 20:
            import warnings
            warnings.warn(
                f"High maximum expression value ({max_val:.1f}) may indicate missing log transformation. "
                "Cell communication inference works best with log1p-transformed data.",
                UserWarning
            )
        
        # Analyze cell communication using selected method
        if params.method == "liana":
            if not LIANA_AVAILABLE:
                raise ImportError("LIANA+ is not installed. Please install it with: pip install liana")
            result_data = await _analyze_communication_liana(adata, params, context)
        
        elif params.method == "cellphonedb":
            if not CELLPHONEDB_AVAILABLE:
                raise ImportError("CellPhoneDB is not installed. Please install it with: pip install cellphonedb")
            result_data = await _analyze_communication_cellphonedb(adata, params, context)
            
        elif params.method == "cellchat_liana":
            if not LIANA_AVAILABLE:
                raise ImportError("LIANA+ is not installed. Please install it with: pip install liana")
            result_data = await _analyze_communication_cellchat_liana(adata, params, context)
            
        else:
            raise ValueError(f"Unsupported method: {params.method}. Supported methods: 'liana', 'cellphonedb', 'cellchat_liana'")
        
        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization
        if context:
            await context.info("Cell communication analysis complete. Use visualize_data tool with plot_type='cell_communication' to visualize results")
        
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
            visualization=None,  # Use visualize_data tool instead
            network_visualization=None,  # Use visualize_data tool instead
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
    """Auto-detect species from gene names using adaptive sampling strategy"""
    # Linus's Rule: No magic numbers, no special cases
    gene_names = _get_representative_gene_sample(adata)

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

    # Calculate confidence scores for robustness
    total_classified = mouse_score + human_score
    confidence = total_classified / len(gene_names) if gene_names else 0
    
    # Note: context.info is async but we can't await in a sync function
    # Store confidence for potential logging by caller
    result_species = "mouse" if mouse_score > human_score else "human"
    
    # Store detection metadata for debugging
    if hasattr(adata, 'uns'):
        adata.uns['species_detection'] = {
            'detected_species': result_species,
            'mouse_score': mouse_score,
            'human_score': human_score,
            'total_genes_sampled': len(gene_names),
            'classification_confidence': confidence
        }
    
    return result_species


def _get_representative_gene_sample(adata: Any) -> set:
    """
    Get representative gene sample using adaptive sampling strategy.
    
    Linus's principle: Eliminate special cases by making the algorithm adapt to data,
    not the other way around.
    """
    total_genes = len(adata.var.index)
    
    if total_genes == 0:
        return set()
    
    # Adaptive sample size based on dataset characteristics
    if total_genes <= 100:
        # Small dataset: use all genes
        sample_size = total_genes
        return set(adata.var.index)
    elif total_genes <= 1000:
        # Medium dataset: use all genes, no sampling needed
        sample_size = total_genes
        return set(adata.var.index)
    else:
        # Large dataset: intelligent sampling
        # Use square root scaling with reasonable bounds
        sample_size = max(500, min(2000, int(np.sqrt(total_genes) * 50)))
        
        # Stratified sampling to avoid ordering bias
        # Take samples from different parts of the gene list
        return _stratified_gene_sampling(adata.var.index, sample_size)


def _stratified_gene_sampling(gene_index: Any, sample_size: int) -> set:
    """
    Perform stratified sampling to get representative genes from different regions.
    
    This eliminates the ordering bias problem that plagued the original [:1000] approach.
    """
    total_genes = len(gene_index)
    
    if sample_size >= total_genes:
        return set(gene_index)
    
    # Divide gene space into strata and sample from each
    n_strata = min(10, sample_size // 50)  # At least 50 genes per stratum
    if n_strata < 2:
        # Fallback to simple random sampling for very small sample sizes
        indices = np.random.choice(total_genes, sample_size, replace=False)
        return set(gene_index[indices])
    
    genes_per_stratum = sample_size // n_strata
    remaining_genes = sample_size % n_strata
    
    sampled_genes = set()
    
    for i in range(n_strata):
        # Calculate stratum boundaries
        stratum_start = (total_genes * i) // n_strata
        stratum_end = (total_genes * (i + 1)) // n_strata
        
        # Number of genes to sample from this stratum
        stratum_sample_size = genes_per_stratum
        if i < remaining_genes:
            stratum_sample_size += 1
        
        # Sample from this stratum
        stratum_size = stratum_end - stratum_start
        if stratum_sample_size >= stratum_size:
            # Take all genes from this stratum
            sampled_genes.update(gene_index[stratum_start:stratum_end])
        else:
            # Random sample from stratum
            stratum_indices = np.random.choice(
                stratum_size, 
                stratum_sample_size, 
                replace=False
            ) + stratum_start
            sampled_genes.update(gene_index[stratum_indices])
    
    return sampled_genes


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


async def _analyze_communication_cellphonedb(
    adata: Any,
    params: CellCommunicationParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze cell communication using CellPhoneDB"""
    try:
        from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
        import tempfile
        import os
    except ImportError:
        raise ImportError("CellPhoneDB is not installed. Please install it with: pip install cellphonedb")

    if context:
        await context.info("Running CellPhoneDB statistical analysis...")

    try:
        import time
        start_time = time.time()

        # Prepare data for CellPhoneDB
        if context:
            await context.info("Preparing data for CellPhoneDB analysis...")

        # Check for required cell type information
        cell_type_col = None
        for col in ['cell_type', 'celltype', 'cluster', 'leiden', 'louvain']:
            if col in adata.obs.columns:
                cell_type_col = col
                break
        
        if not cell_type_col:
            raise ValueError("No cell type information found. Please ensure your data has one of: 'cell_type', 'celltype', 'cluster', 'leiden', 'louvain' columns")

        # Prepare counts data (genes x cells)
        if hasattr(adata.X, 'toarray'):
            counts_df = pd.DataFrame(
                adata.X.toarray().T,
                index=adata.var.index,
                columns=adata.obs.index
            )
        else:
            counts_df = pd.DataFrame(
                adata.X.T,
                index=adata.var.index,
                columns=adata.obs.index
            )

        # Prepare meta data
        meta_df = pd.DataFrame({
            'Cell': adata.obs.index,
            'cell_type': adata.obs[cell_type_col].astype(str)
        })

        if context:
            await context.info(f"Data prepared: {counts_df.shape[0]} genes, {counts_df.shape[1]} cells, {meta_df['cell_type'].nunique()} cell types")

        # Create microenvironments file if spatial data is available and requested
        microenvs_file = None
        if params.cellphonedb_use_microenvironments and 'spatial' in adata.obsm:
            if context:
                await context.info("Creating spatial microenvironments...")
            microenvs_file = await _create_microenvironments_file(adata, params, context)

        # Set random seed for reproducibility
        if params.cellphonedb_debug_seed is not None:
            np.random.seed(params.cellphonedb_debug_seed)

        # Run CellPhoneDB statistical analysis
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save data to temporary files
            counts_file = os.path.join(temp_dir, 'counts.txt')
            meta_file = os.path.join(temp_dir, 'meta.txt')
            
            counts_df.to_csv(counts_file, sep='\t')
            meta_df.to_csv(meta_file, sep='\t', index=False)
            
            if context:
                await context.info("Running CellPhoneDB statistical analysis (this may take several minutes)...")

            try:
                # Run the analysis
                deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
                    cpdb_file_path=None,  # Use default database
                    meta_file_path=meta_file,
                    counts_file_path=counts_file,
                    counts_data='ensembl',
                    threshold=params.cellphonedb_threshold,
                    result_precision=params.cellphonedb_result_precision,
                    pvalue=params.cellphonedb_pvalue,
                    iterations=params.cellphonedb_iterations,
                    debug_seed=params.cellphonedb_debug_seed,
                    output_path=temp_dir,
                    microenvs_file_path=microenvs_file
                )
            except Exception as api_error:
                # Fallback to simpler API call if the above fails
                if context:
                    await context.warning(f"CellPhoneDB API call failed, trying alternative approach: {str(api_error)}")
                
                deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
                    cpdb_file_path=None,
                    meta_file_path=meta_file,
                    counts_file_path=counts_file,
                    threshold=params.cellphonedb_threshold,
                    result_precision=params.cellphonedb_result_precision,
                    pvalue=params.cellphonedb_pvalue,
                    iterations=params.cellphonedb_iterations
                )

            # Store results in AnnData object
            adata.uns['cellphonedb_deconvoluted'] = deconvoluted
            adata.uns['cellphonedb_means'] = means
            adata.uns['cellphonedb_pvalues'] = pvalues
            adata.uns['cellphonedb_significant_means'] = significant_means

        # Calculate statistics
        n_lr_pairs = len(means) if means is not None and hasattr(means, '__len__') else 0
        n_significant_pairs = len(significant_means) if significant_means is not None and hasattr(significant_means, '__len__') else 0

        # Get top LR pairs
        top_lr_pairs = []
        if significant_means is not None and hasattr(significant_means, 'head'):
            # CellPhoneDB returns interactions in 'interacting_pair' column
            if hasattr(significant_means, 'columns') and 'interacting_pair' in significant_means.columns:
                top_pairs_df = significant_means.head(params.plot_top_pairs)
                top_lr_pairs = top_pairs_df['interacting_pair'].tolist()

        end_time = time.time()
        analysis_time = end_time - start_time

        if context:
            await context.info(f"CellPhoneDB analysis completed in {analysis_time:.2f} seconds")
            await context.info(f"Found {n_significant_pairs} significant interactions out of {n_lr_pairs} tested")

        statistics = {
            "method": "cellphonedb",
            "iterations": params.cellphonedb_iterations,
            "threshold": params.cellphonedb_threshold,
            "pvalue_threshold": params.cellphonedb_pvalue,
            "n_cell_types": meta_df['cell_type'].nunique(),
            "microenvironments_used": microenvs_file is not None,
            "analysis_time_seconds": analysis_time
        }

        return {
            "n_lr_pairs": n_lr_pairs,
            "n_significant_pairs": n_significant_pairs,
            "top_lr_pairs": top_lr_pairs,
            "cellphonedb_results_key": "cellphonedb_means",
            "cellphonedb_pvalues_key": "cellphonedb_pvalues",
            "cellphonedb_significant_key": "cellphonedb_significant_means",
            "analysis_type": "statistical",
            "statistics": statistics
        }

    except Exception as e:
        raise RuntimeError(f"CellPhoneDB analysis failed: {str(e)}")


async def _analyze_communication_cellchat_liana(
    adata: Any,
    params: CellCommunicationParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze cell communication using CellChat via LIANA"""
    try:
        import liana as li
    except ImportError:
        raise ImportError("LIANA+ is not installed. Please install it with: pip install liana")

    if context:
        await context.info("Running CellChat analysis via LIANA...")

    try:
        import time
        start_time = time.time()

        # Determine groupby column for CellChat
        groupby_col = None
        for col in ['cell_type', 'celltype', 'cluster', 'leiden', 'louvain']:
            if col in adata.obs.columns:
                groupby_col = col
                break

        if not groupby_col:
            raise ValueError("No suitable groupby column found. Please ensure your data has cell type annotations")

        if context:
            await context.info(f"Running CellChat analysis grouped by '{groupby_col}'...")

        # Get appropriate resource name based on species
        resource_name = _get_liana_resource_name(params.species, "cellchat")
        if context:
            await context.info(f"Using CellChat resource: {resource_name} for species: {params.species}")

        # Use parameters from user or optimize for performance
        n_perms = params.liana_n_perms
        if adata.n_obs > 3000 and n_perms > 500:
            n_perms = 500
            if context:
                await context.info(f"Large dataset detected, reducing permutations to {n_perms}")

        # Run LIANA with CellChat resource
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

        end_time = time.time()
        analysis_time = end_time - start_time

        if context:
            await context.info(f"CellChat via LIANA analysis completed in {analysis_time:.2f} seconds")

        statistics = {
            "method": "cellchat_liana",
            "groupby": groupby_col,
            "n_lr_pairs_tested": n_lr_pairs,
            "n_permutations": n_perms,
            "significance_threshold": 0.05,
            "resource": "cellchat",
            "cellchat_type": params.cellchat_type,
            "analysis_time_seconds": analysis_time
        }

        return {
            "n_lr_pairs": n_lr_pairs,
            "n_significant_pairs": n_significant_pairs,
            "top_lr_pairs": top_lr_pairs,
            "liana_results_key": "liana_res",
            "analysis_type": "cellchat_via_liana",
            "statistics": statistics
        }

    except Exception as e:
        raise RuntimeError(f"CellChat via LIANA analysis failed: {str(e)}")


async def _create_microenvironments_file(
    adata: Any,
    params: CellCommunicationParameters,
    context: Optional[Context] = None
) -> Optional[str]:
    """Create microenvironments file for CellPhoneDB spatial analysis"""
    try:
        import tempfile
        import os
        from sklearn.neighbors import NearestNeighbors
        
        if 'spatial' not in adata.obsm:
            return None
            
        spatial_coords = adata.obsm['spatial']
        
        # Determine spatial radius
        if params.cellphonedb_spatial_radius is not None:
            radius = params.cellphonedb_spatial_radius
        else:
            # Auto-determine radius based on data density
            # Use median distance to 5th nearest neighbor as a heuristic
            nn = NearestNeighbors(n_neighbors=6)
            nn.fit(spatial_coords)
            distances, _ = nn.kneighbors(spatial_coords)
            radius = np.median(distances[:, 5]) * 2  # 5th neighbor (0-indexed), doubled
            
        if context:
            await context.info(f"Using spatial radius: {radius:.2f}")
        
        # Find spatial neighbors for each cell
        nn = NearestNeighbors(radius=radius)
        nn.fit(spatial_coords)
        neighbor_matrix = nn.radius_neighbors_graph(spatial_coords)
        
        # Create microenvironments based on spatial connectivity
        microenvs = []
        for i in range(adata.n_obs):
            neighbors = neighbor_matrix[i].indices
            if len(neighbors) > 1:  # At least one neighbor besides itself
                microenv_cells = [adata.obs.index[j] for j in neighbors]
                microenv_name = f"microenv_{i}"
                for cell in microenv_cells:
                    microenvs.append([cell, microenv_name])
        
        # Save to temporary file
        temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_microenvironments.txt')
        temp_file.write("cell\tmicroenvironment\n")
        for cell, microenv in microenvs:
            temp_file.write(f"{cell}\t{microenv}\n")
        temp_file.close()
        
        if context:
            await context.info(f"Created microenvironments file with {len(set([m[1] for m in microenvs]))} microenvironments")
            
        return temp_file.name
        
    except Exception as e:
        if context:
            await context.warning(f"Failed to create microenvironments file: {str(e)}")
        return None


async def _comprehensive_data_validation(
    adata: Any,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Comprehensive data quality validation for AnnData objects
    
    This function implements Linus's "good taste" principle by:
    1. Checking real problems that occur in production
    2. Failing early with clear error messages
    3. Providing actionable suggestions for fixes
    4. Using a unified validation framework
    
    Args:
        adata: AnnData object to validate
        context: MCP context for logging
        
    Returns:
        Dict with validation results:
        {
            "passed": bool,
            "error_message": str,
            "warnings": List[str],
            "suggestions": str,
            "validation_details": Dict[str, Any]
        }
    """
    validation_result = {
        "passed": True,
        "error_message": "",
        "warnings": [],
        "suggestions": "",
        "validation_details": {}
    }
    
    errors = []
    warnings_list = []
    suggestions = []
    
    try:
        # 1. Basic structure validation
        structure_check = await _validate_basic_structure(adata, context)
        validation_result["validation_details"]["structure"] = structure_check
        if not structure_check["passed"]:
            errors.extend(structure_check["errors"])
            suggestions.extend(structure_check["suggestions"])
        warnings_list.extend(structure_check["warnings"])
        
        # 2. Expression matrix validation
        expression_check = await _validate_expression_matrix(adata, context)
        validation_result["validation_details"]["expression"] = expression_check
        if not expression_check["passed"]:
            errors.extend(expression_check["errors"])
            suggestions.extend(expression_check["suggestions"])
        warnings_list.extend(expression_check["warnings"])
        
        # 3. Spatial coordinates validation
        spatial_check = await _validate_spatial_coordinates(adata, context)
        validation_result["validation_details"]["spatial"] = spatial_check
        if not spatial_check["passed"]:
            errors.extend(spatial_check["errors"])
            suggestions.extend(spatial_check["suggestions"])
        warnings_list.extend(spatial_check["warnings"])
        
        # 4. Metadata validation
        metadata_check = await _validate_metadata(adata, context)
        validation_result["validation_details"]["metadata"] = metadata_check
        if not metadata_check["passed"]:
            errors.extend(metadata_check["errors"])
            suggestions.extend(metadata_check["suggestions"])
        warnings_list.extend(metadata_check["warnings"])
        
        # 5. Cell communication specific validation
        comm_check = await _validate_communication_requirements(adata, context)
        validation_result["validation_details"]["communication"] = comm_check
        if not comm_check["passed"]:
            errors.extend(comm_check["errors"])
            suggestions.extend(comm_check["suggestions"])
        warnings_list.extend(comm_check["warnings"])
        
        # Compile final results
        if errors:
            validation_result["passed"] = False
            validation_result["error_message"] = "\n".join([f"• {error}" for error in errors])
        
        validation_result["warnings"] = warnings_list
        if suggestions:
            validation_result["suggestions"] = "\n".join([f"• {suggestion}" for suggestion in suggestions])
        
        # Log summary
        if context:
            if validation_result["passed"]:
                await context.info(f"✅ Data validation passed with {len(warnings_list)} warnings")
                if warnings_list:
                    for warning in warnings_list[:3]:  # Show first 3 warnings
                        await context.info(f"⚠️  {warning}")
            else:
                await context.warning(f"❌ Data validation failed: {len(errors)} critical issues found")
        
        return validation_result
        
    except Exception as e:
        validation_result["passed"] = False
        validation_result["error_message"] = f"Validation process failed: {str(e)}"
        validation_result["suggestions"] = "Please check your data format and try again"
        return validation_result


async def _validate_basic_structure(
    adata: Any,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Validate basic AnnData structure"""
    result = {"passed": True, "errors": [], "warnings": [], "suggestions": []}
    
    try:
        # Check if it's actually an AnnData object
        if not hasattr(adata, 'X') or not hasattr(adata, 'obs') or not hasattr(adata, 'var'):
            result["errors"].append("Invalid AnnData object: missing X, obs, or var attributes")
            result["suggestions"].append("Ensure you're passing a valid AnnData object")
            result["passed"] = False
            return result
        
        # Check dimensions consistency
        if adata.X.shape[0] != len(adata.obs):
            result["errors"].append(f"Dimension mismatch: X has {adata.X.shape[0]} cells but obs has {len(adata.obs)} entries")
            result["suggestions"].append("Check data loading process - cells and observations must match")
            result["passed"] = False
        
        if adata.X.shape[1] != len(adata.var):
            result["errors"].append(f"Dimension mismatch: X has {adata.X.shape[1]} genes but var has {len(adata.var)} entries")
            result["suggestions"].append("Check data loading process - genes and variables must match")
            result["passed"] = False
        
        # Check for empty data
        if adata.n_obs == 0:
            result["errors"].append("Dataset is empty: no cells found")
            result["suggestions"].append("Load data with actual cell measurements")
            result["passed"] = False
        
        if adata.n_vars == 0:
            result["errors"].append("Dataset is empty: no genes found")
            result["suggestions"].append("Load data with actual gene measurements")
            result["passed"] = False
        
        # Data size warnings
        if adata.n_obs < 50:
            result["warnings"].append(f"Very few cells ({adata.n_obs}) - cell communication analysis may be unreliable")
        
        if adata.n_vars < 1000:
            result["warnings"].append(f"Few genes ({adata.n_vars}) - consider using more comprehensive gene set")
        
    except Exception as e:
        result["errors"].append(f"Structure validation failed: {str(e)}")
        result["suggestions"].append("Check if the input is a valid AnnData object")
        result["passed"] = False
    
    return result


async def _validate_expression_matrix(
    adata: Any,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Validate expression matrix data quality"""
    result = {"passed": True, "errors": [], "warnings": [], "suggestions": []}
    
    try:
        X = adata.X
        
        # Convert sparse to dense for NaN/Inf checking (sample if too large)
        if sparse.issparse(X):
            # For large matrices, sample to check quality
            if X.shape[0] > 5000 or X.shape[1] > 2000:
                sample_cells = min(1000, X.shape[0])
                sample_genes = min(500, X.shape[1])
                cell_idx = np.random.choice(X.shape[0], sample_cells, replace=False)
                gene_idx = np.random.choice(X.shape[1], sample_genes, replace=False)
                X_sample = X[cell_idx, :][:, gene_idx].toarray()
            else:
                X_sample = X.toarray()
        else:
            X_sample = X if X.size <= 10_000_000 else X[:1000, :500]  # Sample large dense matrices
        
        # Check for NaN values
        if np.isnan(X_sample).any():
            nan_count = np.isnan(X_sample).sum()
            nan_fraction = nan_count / X_sample.size
            if nan_fraction > 0.1:  # More than 10% NaN
                result["errors"].append(f"High proportion of NaN values ({nan_fraction:.2%}) in expression matrix")
                result["suggestions"].append("Remove or impute NaN values before analysis")
                result["passed"] = False
            else:
                result["warnings"].append(f"Found {nan_count} NaN values ({nan_fraction:.2%}) in expression matrix")
        
        # Check for infinite values
        if np.isinf(X_sample).any():
            inf_count = np.isinf(X_sample).sum()
            result["errors"].append(f"Found {inf_count} infinite values in expression matrix")
            result["suggestions"].append("Replace infinite values with finite numbers or remove affected cells/genes")
            result["passed"] = False
        
        # Check for negative values (suspicious in count data)
        if (X_sample < 0).any():
            neg_count = (X_sample < 0).sum()
            neg_fraction = neg_count / X_sample.size
            if neg_fraction > 0.01:  # More than 1% negative
                result["warnings"].append(f"Found {neg_count} negative values ({neg_fraction:.2%}) - unusual for count data")
        
        # Check for extremely large values
        max_val = X_sample.max()
        if max_val > 1e6:
            result["warnings"].append(f"Very large expression values detected (max: {max_val:.2e}) - consider normalization")
        
        # Check sparsity
        if sparse.issparse(X):
            sparsity = 1.0 - (X.nnz / (X.shape[0] * X.shape[1]))
        else:
            sparsity = (X == 0).sum() / X.size
        
        if sparsity > 0.99:
            result["warnings"].append(f"Very sparse data ({sparsity:.1%} zeros) - ensure sufficient sequencing depth")
        elif sparsity < 0.3:
            result["warnings"].append(f"Unusually dense data ({sparsity:.1%} zeros) - verify data normalization")
        
        # Check data type
        if not np.issubdtype(X.dtype, np.number):
            result["errors"].append(f"Expression matrix has non-numeric data type: {X.dtype}")
            result["suggestions"].append("Convert expression data to numeric format")
            result["passed"] = False
        
    except Exception as e:
        result["errors"].append(f"Expression matrix validation failed: {str(e)}")
        result["suggestions"].append("Check expression matrix format and values")
        result["passed"] = False
    
    return result


async def _validate_spatial_coordinates(
    adata: Any,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Validate spatial coordinates for cell communication analysis"""
    result = {"passed": True, "errors": [], "warnings": [], "suggestions": []}
    
    try:
        # Check for spatial coordinates existence
        spatial_key = None
        for key in adata.obsm.keys():
            if 'spatial' in key.lower():
                spatial_key = key
                break
        
        if spatial_key is None:
            result["errors"].append("No spatial coordinates found in adata.obsm")
            result["suggestions"].append("Add spatial coordinates to adata.obsm['spatial'] or similar key")
            result["passed"] = False
            return result
        
        spatial_coords = adata.obsm[spatial_key]
        
        # Check spatial coordinates shape
        if spatial_coords.shape[0] != adata.n_obs:
            result["errors"].append(f"Spatial coordinates shape mismatch: {spatial_coords.shape[0]} != {adata.n_obs} cells")
            result["suggestions"].append("Ensure spatial coordinates match number of cells")
            result["passed"] = False
        
        if spatial_coords.shape[1] < 2:
            result["errors"].append(f"Spatial coordinates must have at least 2 dimensions, found {spatial_coords.shape[1]}")
            result["suggestions"].append("Provide x,y coordinates (and optionally z) for spatial analysis")
            result["passed"] = False
        
        # Check for NaN/Inf in spatial coordinates
        if np.isnan(spatial_coords).any():
            nan_count = np.isnan(spatial_coords).sum()
            result["errors"].append(f"Found {nan_count} NaN values in spatial coordinates")
            result["suggestions"].append("Remove cells with missing spatial coordinates")
            result["passed"] = False
        
        if np.isinf(spatial_coords).any():
            inf_count = np.isinf(spatial_coords).sum()
            result["errors"].append(f"Found {inf_count} infinite values in spatial coordinates")
            result["suggestions"].append("Replace infinite spatial coordinates with valid values")
            result["passed"] = False
        
        # Check for duplicate coordinates (suspicious)
        unique_coords = np.unique(spatial_coords, axis=0)
        if len(unique_coords) < len(spatial_coords) * 0.95:  # More than 5% duplicates
            duplicate_fraction = 1 - (len(unique_coords) / len(spatial_coords))
            result["warnings"].append(f"High proportion of duplicate spatial coordinates ({duplicate_fraction:.2%})")
        
        # Check coordinate range reasonableness
        coord_ranges = []
        for dim in range(min(3, spatial_coords.shape[1])):  # Check up to 3 dimensions
            coord_min, coord_max = spatial_coords[:, dim].min(), spatial_coords[:, dim].max()
            coord_range = coord_max - coord_min
            coord_ranges.append(coord_range)
            
            if coord_range == 0:
                result["warnings"].append(f"All cells have identical coordinate in dimension {dim + 1}")
            elif coord_range > 1e6:
                result["warnings"].append(f"Very large coordinate range in dimension {dim + 1}: {coord_range:.2e}")
        
        # Check if coordinates are integer-like (suggesting pixel coordinates)
        if len(coord_ranges) >= 2:
            x_coords, y_coords = spatial_coords[:, 0], spatial_coords[:, 1]
            x_is_int = np.allclose(x_coords, np.round(x_coords))
            y_is_int = np.allclose(y_coords, np.round(y_coords))
            
            if x_is_int and y_is_int and coord_ranges[0] > 100 and coord_ranges[1] > 100:
                result["warnings"].append("Coordinates appear to be in pixel units - consider converting to physical units")
        
    except Exception as e:
        result["errors"].append(f"Spatial coordinates validation failed: {str(e)}")
        result["suggestions"].append("Check spatial coordinates format and values")
        result["passed"] = False
    
    return result


async def _validate_metadata(
    adata: Any,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Validate observation and variable metadata"""
    result = {"passed": True, "errors": [], "warnings": [], "suggestions": []}
    
    try:
        # Validate observation metadata (adata.obs)
        if adata.obs.index.duplicated().any():
            dup_count = adata.obs.index.duplicated().sum()
            result["errors"].append(f"Found {dup_count} duplicate cell IDs in adata.obs")
            result["suggestions"].append("Ensure all cell IDs are unique")
            result["passed"] = False
        
        # Check for empty cell IDs
        empty_ids = adata.obs.index.isna() | (adata.obs.index == '')
        if empty_ids.any():
            empty_count = empty_ids.sum()
            result["errors"].append(f"Found {empty_count} empty cell IDs")
            result["suggestions"].append("Provide valid cell IDs for all cells")
            result["passed"] = False
        
        # Validate variable metadata (adata.var)
        if adata.var.index.duplicated().any():
            dup_count = adata.var.index.duplicated().sum()
            result["errors"].append(f"Found {dup_count} duplicate gene names in adata.var")
            result["suggestions"].append("Ensure all gene names are unique or use adata.var_names_unique()")
            result["passed"] = False
        
        # Check for empty gene names
        empty_genes = adata.var.index.isna() | (adata.var.index == '')
        if empty_genes.any():
            empty_count = empty_genes.sum()
            result["errors"].append(f"Found {empty_count} empty gene names")
            result["suggestions"].append("Provide valid gene names for all variables")
            result["passed"] = False
        
        # Check for suspicious gene name patterns
        gene_names = adata.var.index.astype(str)
        numeric_genes = pd.to_numeric(gene_names, errors='coerce').notna()
        if numeric_genes.sum() > len(gene_names) * 0.5:
            result["warnings"].append(f"Many genes have numeric names ({numeric_genes.sum()}/{len(gene_names)})")
        
    except Exception as e:
        result["errors"].append(f"Metadata validation failed: {str(e)}")
        result["suggestions"].append("Check observation and variable metadata format")
        result["passed"] = False
    
    return result


async def _validate_communication_requirements(
    adata: Any,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Validate specific requirements for cell communication analysis"""
    result = {"passed": True, "errors": [], "warnings": [], "suggestions": []}
    
    try:
        # Check for minimum cell count for meaningful analysis
        min_cells_recommended = 100
        if adata.n_obs < min_cells_recommended:
            result["warnings"].append(
                f"Low cell count ({adata.n_obs}) may lead to unreliable cell communication results. "
                f"Recommended: >{min_cells_recommended} cells"
            )
        
        # Check for cell type annotations (helpful for interpretation)
        cell_type_cols = [col for col in adata.obs.columns if 'type' in col.lower() or 'cluster' in col.lower()]
        if not cell_type_cols:
            result["warnings"].append(
                "No cell type annotations found. Consider adding cell type information for better interpretation"
            )
        else:
            # Check cell type diversity
            for col in cell_type_cols[:1]:  # Check first cell type column
                n_types = adata.obs[col].nunique()
                if n_types < 2:
                    result["warnings"].append(f"Only {n_types} cell type(s) found - communication analysis needs diversity")
                elif n_types > 50:
                    result["warnings"].append(f"Very many cell types ({n_types}) - consider grouping similar types")
        
        # Check gene expression levels
        if hasattr(adata.X, 'max'):
            max_expr = adata.X.max()
            if max_expr < 1:
                result["warnings"].append(
                    f"Low maximum expression value ({max_expr:.3f}) suggests data may need normalization"
                )
            elif max_expr > 50 and 'log1p' not in adata.uns:
                result["warnings"].append(
                    f"High maximum expression value ({max_expr:.1f}) suggests data may need log transformation"
                )
        
        # Check for mitochondrial genes (quality control indicator)
        gene_names = adata.var.index.astype(str).str.upper()
        mt_genes = gene_names.str.startswith('MT-') | gene_names.str.startswith('MT.') | gene_names.str.contains('^MT-')
        n_mt_genes = mt_genes.sum()
        
        if n_mt_genes == 0:
            result["warnings"].append("No mitochondrial genes detected - consider quality control filtering")
        elif n_mt_genes > len(gene_names) * 0.1:
            result["warnings"].append(f"High proportion of mitochondrial genes ({n_mt_genes}/{len(gene_names)})")
        
    except Exception as e:
        result["errors"].append(f"Communication requirements validation failed: {str(e)}")
        result["suggestions"].append("Check data quality and preprocessing steps")
        result["passed"] = False
    
    return result
