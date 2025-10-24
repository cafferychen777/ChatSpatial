"""
Cell-cell communication analysis tools for spatial transcriptomics data.
"""

from typing import Any, Dict, Optional

import numpy as np
from mcp.server.fastmcp import Context

from ..models.analysis import CellCommunicationResult
from ..models.data import CellCommunicationParameters
from ..utils import validate_obs_column

# Import LIANA+ for cell communication analysis
try:
    import liana as li  # noqa: F401

    LIANA_AVAILABLE = True
except ImportError:
    LIANA_AVAILABLE = False

# Check CellPhoneDB availability
try:
    import importlib.util

    CELLPHONEDB_AVAILABLE = importlib.util.find_spec("cellphonedb") is not None
except ImportError:
    CELLPHONEDB_AVAILABLE = False


async def _validate_liana_requirements(
    adata: Any, params: CellCommunicationParameters, context: Optional[Context] = None
) -> None:
    """Validate LIANA+ requirements"""
    # Spatial connectivity validation
    if params.perform_spatial_analysis and "spatial_connectivities" not in adata.obsp:
        raise ValueError(
            "Spatial connectivity required for LIANA+ bivariate analysis.\n\n"
            "Run spatial neighbor computation first:\n"
            "  import squidpy as sq\n"
            "  sq.gr.spatial_neighbors(adata, coord_type='grid', n_rings=1)\n\n"
            "Platform-specific recommendations:\n"
            "  Visium: coord_type='grid', n_rings=1-2\n"
            "  MERFISH: coord_type='generic', radius=20-50\n"
            "  Slide-seq: coord_type='generic', n_neighs=10-30"
        )

    # Cell type validation
    validate_obs_column(adata, params.cell_type_key, "Cell type column")

    # Warning for resource matching
    if params.species == "mouse" and params.liana_resource == "consensus" and context:
        await context.warning(
            "Using 'consensus' for mouse data. Consider liana_resource='mouseconsensus'."
        )


async def analyze_cell_communication(
    data_id: str,
    data_store: Dict[str, Any],
    params: CellCommunicationParameters,  # No default - must be provided by caller (LLM)
    context: Optional[Context] = None,
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
    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Direct reference to avoid unnecessary copies
    adata = data_store[data_id]["adata"]

    try:
        # Apply method-specific validation
        if params.method in ["liana", "cellchat_liana"]:
            # LIANA-based methods need spatial connectivity validation
            await _validate_liana_requirements(adata, params, context)
        elif params.method == "cellphonedb":
            # Check if cell type column exists
            validate_obs_column(adata, params.cell_type_key, "Cell type column")

            # Provide data overview information
            if context:
                data_source_info = (
                    "raw layer"
                    if params.data_source == "raw" and adata.raw
                    else "current layer"
                )
                n_genes = (
                    adata.raw.n_vars
                    if params.data_source == "raw" and adata.raw
                    else adata.n_vars
                )
                await context.info(
                    f"Running CellPhoneDB with {n_genes} genes from {data_source_info} "
                    f"and {adata.n_obs} cells"
                )

                # Warnings for low counts
                if n_genes < 5000:
                    await context.warning(
                        f"Gene count ({n_genes}) is relatively low. "
                        f"This may limit the number of interactions found."
                    )

                if adata.n_obs < 100:
                    await context.warning(
                        f"Cell count ({adata.n_obs}) is relatively low. "
                        f"This may affect statistical power."
                    )

        # Handle data source selection
        if params.data_source == "raw":
            if adata.raw is None:
                raise ValueError(
                    "data_source='raw' specified but adata.raw is None. "
                    "Either use data_source='current' or load data with raw counts preserved."
                )

            # Create full dataset from raw
            adata_full = adata.raw.to_adata()

            # Transfer all metadata
            adata_full.obs = adata.obs.copy()
            adata_full.obsm = adata.obsm.copy()
            adata_full.obsp = adata.obsp.copy()  # Critical for spatial connectivity
            adata_full.uns = adata.uns.copy()

            # Check if raw data needs normalization
            if adata_full.X.max() > 100:
                import scanpy as sc

                sc.pp.normalize_total(adata_full, target_sum=1e4)
                sc.pp.log1p(adata_full)

            adata = adata_full

        # Analyze cell communication using selected method
        if params.method == "liana":
            if not LIANA_AVAILABLE:
                raise ImportError(
                    "LIANA+ is not installed. Please install it with: pip install liana"
                )
            result_data = await _analyze_communication_liana(adata, params, context)

        elif params.method == "cellphonedb":
            if not CELLPHONEDB_AVAILABLE:
                raise ImportError(
                    "CellPhoneDB is not installed. Please install it with: pip install cellphonedb"
                )
            result_data = await _analyze_communication_cellphonedb(
                adata, params, context
            )

        elif params.method == "cellchat_liana":
            if not LIANA_AVAILABLE:
                raise ImportError(
                    "LIANA+ is not installed. Please install it with: pip install liana"
                )
            result_data = await _analyze_communication_cellchat_liana(
                adata, params, context
            )

        else:
            raise ValueError(
                f"Unsupported method: {params.method}. Supported methods: 'liana', 'cellphonedb', 'cellchat_liana'"
            )

        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization
        if context:
            await context.info(
                "Cell communication analysis complete. Use visualize_data tool with plot_type='cell_communication' to visualize results"
            )

        # When data_source="raw", a new adata object is created, so copy results back
        original_adata = data_store[data_id]["adata"]

        # Copy all CellPhoneDB/LIANA results from the temporary adata to the original
        # This ensures results are preserved when data is saved
        for key in [
            "cellphonedb_deconvoluted",
            "cellphonedb_means",
            "cellphonedb_pvalues",
            "cellphonedb_significant_means",
            "cellphonedb_statistics",
            "liana_res",
            "lrdata",
            "liana_spatial",
        ]:
            if key in adata.uns:
                original_adata.uns[key] = adata.uns[key]

        # Also copy any obsm keys that were added
        for key in adata.obsm.keys():
            if key not in original_adata.obsm or key.startswith("liana_"):
                original_adata.obsm[key] = adata.obsm[key]

        # Store scientific metadata for reproducibility
        from ..utils.metadata_storage import store_analysis_metadata

        # Determine database used
        if params.method == "liana":
            database = params.liana_resource
        elif params.method == "cellphonedb":
            database = "cellphonedb"
        elif params.method == "cellchat_liana":
            database = "cellchatdb"  # Match actual LIANA resource name used in implementation (line 1801)
        else:
            database = "unknown"

        # Extract results keys
        results_keys_dict = {"obs": [], "obsm": [], "uns": []}

        if result_data.get("liana_results_key"):
            results_keys_dict["uns"].append(result_data["liana_results_key"])
        if result_data.get("liana_spatial_results_key"):
            results_keys_dict["uns"].append(result_data["liana_spatial_results_key"])
        if result_data.get("liana_spatial_scores_key"):
            results_keys_dict["obsm"].append(result_data["liana_spatial_scores_key"])
        if result_data.get("cellphonedb_results_key"):
            results_keys_dict["uns"].append(result_data["cellphonedb_results_key"])

        # Store metadata
        store_analysis_metadata(
            adata,
            analysis_name=f"cell_communication_{params.method}",
            method=params.method,
            parameters={
                "cell_type_key": params.cell_type_key,
                "n_perms": (
                    params.liana_n_perms
                    if params.method in ["liana", "cellchat_liana"]
                    else None
                ),
                "nz_prop": (
                    params.liana_nz_prop
                    if params.method in ["liana", "cellchat_liana"]
                    else None
                ),
                "min_cells": params.min_cells,
                "iterations": (
                    params.cellphonedb_iterations
                    if params.method == "cellphonedb"
                    else None
                ),
                "threshold": (
                    params.cellphonedb_threshold
                    if params.method == "cellphonedb"
                    else None
                ),
            },
            results_keys=results_keys_dict,
            statistics={
                "n_lr_pairs": result_data["n_lr_pairs"],
                "n_significant_pairs": result_data["n_significant_pairs"],
                "analysis_type": result_data.get("analysis_type"),
            },
            species=params.species,
            database=database,
        )

        # Create result
        result = CellCommunicationResult(
            data_id=data_id,
            method=params.method,
            species=params.species,
            database=database,  # Use actual database/resource determined above
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
            statistics=result_data.get("statistics", {}),
        )

        if context:
            await context.info(
                f"Successfully analyzed {result.n_significant_pairs} significant LR pairs"
            )
            if result.top_lr_pairs:
                await context.info(f"Top LR pair: {result.top_lr_pairs[0]}")

        return result

    except Exception as e:
        error_msg = f"Error in cell communication analysis: {str(e)}"
        if context:
            await context.warning(error_msg)
        raise RuntimeError(error_msg)


async def _analyze_communication_liana(
    adata: Any, params: CellCommunicationParameters, context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze cell communication using LIANA+"""
    try:
        import liana as li  # noqa: F401
    except ImportError:
        raise ImportError(
            "LIANA+ is not installed. Please install it with: pip install liana"
        )

    try:
        import time

        time.time()

        # Ensure spatial connectivity is computed
        if "spatial_connectivities" not in adata.obsp:
            # Use parameters from user or determine optimal bandwidth based on data size
            if params.liana_bandwidth is not None:
                bandwidth = params.liana_bandwidth
            elif adata.n_obs > 3000:
                bandwidth = 300  # Larger bandwidth for large datasets
            else:
                bandwidth = 200  # Standard bandwidth

            cutoff = params.liana_cutoff

            # Use Squidpy for spatial neighbor computation
            # Note: Spatial analysis requires spatial neighbors (physical coordinates), not expression neighbors
            try:
                import squidpy as sq

                # Squidpy's spatial_neighbors uses PHYSICAL coordinates
                sq.gr.spatial_neighbors(
                    adata,
                    coord_type="generic",
                    n_neighs=min(
                        30, max(6, adata.n_obs // 100)
                    ),  # Adaptive neighbor count
                    radius=bandwidth if bandwidth else None,
                    delaunay=True,  # Use Delaunay triangulation for spatial data
                    set_diag=False,  # Standard practice for spatial graphs
                )
            except ImportError:
                raise ImportError(
                    "Squidpy required for spatial analysis (computes neighbors based on physical coordinates). "
                    "Install: pip install squidpy"
                )

            if context:
                await context.info(
                    f"Spatial connectivity computed with bandwidth={bandwidth}, cutoff={cutoff}"
                )

        # Validate species parameter is specified
        if not params.species:
            raise ValueError(
                "Species parameter is required!\n\n"
                "You must explicitly specify the species of your data:\n"
                "  - species='human': For human data (genes like ACTB, GAPDH)\n"
                "  - species='mouse': For mouse data (genes like Actb, Gapdh)\n"
                "  - species='zebrafish': For zebrafish data\n\n"
                "Example usage:\n"
                "  params = {\n"
                "      'species': 'mouse',\n"
                "      'cell_type_key': 'cell_type',\n"
                "      'liana_resource': 'mouseconsensus'\n"
                "  }"
            )

        # Determine analysis type based on data characteristics
        has_clusters = params.cell_type_key in adata.obs.columns

        if has_clusters and not params.perform_spatial_analysis:
            # Single-cell style analysis with clusters
            return await _run_liana_cluster_analysis(adata, params, context)
        else:
            # Spatial bivariate analysis
            return await _run_liana_spatial_analysis(adata, params, context)

    except Exception as e:
        raise RuntimeError(f"LIANA+ analysis failed: {str(e)}")


def _get_liana_resource_name(species: str, resource_preference: str) -> str:
    """Get appropriate LIANA+ resource name based on species with enhanced resource support"""
    if species == "mouse":
        # Mouse-specific resources
        mouse_resources = ["mouseconsensus", "cellphonedb", "celltalkdb", "icellnet"]

        if resource_preference == "consensus":
            return "mouseconsensus"  # Auto-map consensus to mouseconsensus for mouse
        elif resource_preference in mouse_resources:
            return (
                resource_preference  # Use as specified if it's a valid mouse resource
            )
        else:
            # For non-mouse-specific resources, still use them but could warn
            return resource_preference
    else:
        # For human or other species, use as specified
        return resource_preference


async def _run_liana_cluster_analysis(
    adata: Any, params: CellCommunicationParameters, context: Optional[Context] = None
) -> Dict[str, Any]:
    """Run LIANA+ cluster-based analysis"""
    import liana as li

    # Use cell_type_key from params (required field, no auto-detect)
    groupby_col = params.cell_type_key

    validate_obs_column(adata, groupby_col, "Cell type column")

    # Get appropriate resource name based on species
    resource_name = _get_liana_resource_name(params.species, params.liana_resource)

    # Use parameters from user (respect user choice)
    n_perms = params.liana_n_perms

    # Run LIANA+ rank aggregate
    li.mt.rank_aggregate(
        adata,
        groupby=groupby_col,
        resource_name=resource_name,
        expr_prop=params.liana_nz_prop,
        min_cells=params.min_cells,
        n_perms=n_perms,
        verbose=False,
        use_raw=True if adata.raw is not None else False,
    )

    # Get results
    liana_res = adata.uns["liana_res"]

    # Calculate statistics using magnitude_rank (signal strength)
    # NOT specificity_rank (which has non-uniform distribution)
    n_lr_pairs = len(liana_res)
    n_significant_pairs = len(liana_res[liana_res["magnitude_rank"] <= 0.05])

    # Get top pairs
    top_lr_pairs = []
    if "magnitude_rank" in liana_res.columns:
        top_pairs_df = liana_res.nsmallest(params.plot_top_pairs, "magnitude_rank")
        top_lr_pairs = [
            f"{row['ligand_complex']}_{row['receptor_complex']}"
            for _, row in top_pairs_df.iterrows()
        ]

    # Store standardized L-R pairs for visualization
    detected_lr_pairs = []
    if "magnitude_rank" in liana_res.columns:
        for _, row in top_pairs_df.iterrows():
            ligand = row["ligand_complex"]
            receptor = row["receptor_complex"]
            detected_lr_pairs.append((ligand, receptor))

    # Store in standardized format for visualization
    adata.uns["detected_lr_pairs"] = detected_lr_pairs
    adata.uns["cell_communication_results"] = {
        "top_lr_pairs": top_lr_pairs,
        "method": "liana_cluster",
        "n_pairs": len(top_lr_pairs),
        "species": params.species,
    }

    statistics = {
        "method": "liana_cluster",
        "groupby": groupby_col,
        "n_lr_pairs_tested": n_lr_pairs,
        "n_permutations": n_perms,
        "significance_threshold": 0.05,
        "resource": params.liana_resource,
    }

    return {
        "n_lr_pairs": n_lr_pairs,
        "n_significant_pairs": n_significant_pairs,
        "top_lr_pairs": top_lr_pairs,
        # "liana_results_key": "liana_res",  # Removed to prevent potential DataFrame serialization overflow
        "analysis_type": "cluster",
        "statistics": statistics,
    }


async def _run_liana_spatial_analysis(
    adata: Any, params: CellCommunicationParameters, context: Optional[Context] = None
) -> Dict[str, Any]:
    """Run LIANA+ spatial bivariate analysis"""
    import liana as li

    # Get appropriate resource name based on species
    resource_name = _get_liana_resource_name(params.species, params.liana_resource)

    # Use parameters from user (respect user choice)
    n_perms = params.liana_n_perms
    nz_prop = params.liana_nz_prop

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
        verbose=False,
    )

    # Get results summary
    n_lr_pairs = lrdata.n_vars

    # Get top pairs based on global metric
    global_metric = params.liana_global_metric
    top_pairs_df = lrdata.var.nlargest(params.plot_top_pairs, global_metric)
    top_lr_pairs = top_pairs_df.index.tolist()

    # Count significant pairs (high spatial autocorrelation)
    # Use different thresholds based on metric
    if global_metric == "morans":
        threshold = 0.1
    else:  # lee
        threshold = 0.1

    # Both morans and lee use > threshold for significance
    n_significant_pairs = len(lrdata.var[lrdata.var[global_metric] > threshold])

    # Store results in adata
    adata.uns["liana_spatial_res"] = lrdata.var
    adata.obsm["liana_spatial_scores"] = lrdata.X.toarray()
    adata.uns["liana_spatial_interactions"] = lrdata.var.index.tolist()

    if "pvals" in lrdata.layers:
        adata.obsm["liana_spatial_pvals"] = lrdata.layers["pvals"].toarray()

    if "cats" in lrdata.layers:
        adata.obsm["liana_spatial_cats"] = lrdata.layers["cats"].toarray()

    # Store standardized L-R pairs for visualization
    detected_lr_pairs = []
    for pair_str in top_lr_pairs:
        if "^" in pair_str:
            ligand, receptor = pair_str.split("^", 1)
            detected_lr_pairs.append((ligand, receptor))
        elif "_" in pair_str:
            parts = pair_str.split("_")
            if len(parts) == 2:
                detected_lr_pairs.append((parts[0], parts[1]))

    # Store in standardized format for visualization
    adata.uns["detected_lr_pairs"] = detected_lr_pairs
    adata.uns["cell_communication_results"] = {
        "top_lr_pairs": top_lr_pairs,
        "method": "liana_spatial",
        "n_pairs": len(top_lr_pairs),
        "species": params.species,
    }

    statistics = {
        "method": "liana_spatial",
        "local_metric": params.liana_local_metric,
        "global_metric": params.liana_global_metric,
        "n_lr_pairs_tested": n_lr_pairs,
        "n_permutations": n_perms,
        "nz_proportion": nz_prop,
        "resource": params.liana_resource,
    }

    return {
        "n_lr_pairs": n_lr_pairs,
        "n_significant_pairs": n_significant_pairs,
        "top_lr_pairs": top_lr_pairs,
        "liana_spatial_results_key": "liana_spatial_res",
        "liana_spatial_scores_key": "liana_spatial_scores",
        "analysis_type": "spatial",
        "statistics": statistics,
    }


async def _ensure_cellphonedb_database(
    output_dir: str, context: Optional[Context] = None
) -> str:
    """Ensure CellPhoneDB database is available, download if not exists"""
    try:
        from cellphonedb.utils import db_utils
    except ImportError:
        raise ImportError(
            "CellPhoneDB utils not available. Please install with: pip install cellphonedb"
        )

    import os

    # Check if database file already exists
    db_path = os.path.join(output_dir, "cellphonedb.zip")

    if os.path.exists(db_path):
        if context:
            await context.info(f"Using existing CellPhoneDB database: {db_path}")
        return db_path

    try:
        if context:
            await context.info("Downloading CellPhoneDB database (v5.0.0)...")

        # Download latest database
        db_utils.download_database(output_dir, "v5.0.0")

        if context:
            await context.info(
                f"Successfully downloaded CellPhoneDB database to: {db_path}"
            )

        return db_path

    except Exception as e:
        error_msg = (
            f"Failed to download CellPhoneDB database: {str(e)}\n\n"
            "Troubleshooting:\n"
            "1. Check internet connection\n"
            "2. Verify CellPhoneDB version compatibility\n"
            "3. Try manually downloading database:\n"
            "   from cellphonedb.utils import db_utils\n"
            "   db_utils.download_database('/path/to/dir', 'v5.0.0')"
        )
        raise RuntimeError(error_msg) from e


async def _analyze_communication_cellphonedb(
    adata: Any, params: CellCommunicationParameters, context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze cell communication using CellPhoneDB"""
    try:
        import os
        import tempfile

        from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
    except ImportError:
        raise ImportError(
            "CellPhoneDB is not installed. Please install it with: pip install cellphonedb"
        )

    try:
        import time

        start_time = time.time()

        # Use cell_type_key from params (required field, no auto-detect)
        cell_type_col = params.cell_type_key

        validate_obs_column(adata, cell_type_col, "Cell type column")

        # Use original adata directly (no gene filtering needed)
        adata_for_analysis = adata

        # Import pandas for DataFrame operations
        import pandas as pd

        # Prepare counts data (genes x cells) - Use filtered data
        if hasattr(adata_for_analysis.X, "toarray"):
            counts_df = pd.DataFrame(
                adata_for_analysis.X.toarray().T,
                index=adata_for_analysis.var.index,
                columns=adata_for_analysis.obs.index,
            )
        else:
            counts_df = pd.DataFrame(
                adata_for_analysis.X.T,
                index=adata_for_analysis.var.index,
                columns=adata_for_analysis.obs.index,
            )

        # Prepare meta data
        meta_df = pd.DataFrame(
            {
                "Cell": adata_for_analysis.obs.index,
                "cell_type": adata_for_analysis.obs[cell_type_col].astype(str),
            }
        )

        if context:
            await context.info(
                f"Data prepared: {counts_df.shape[0]} genes, {counts_df.shape[1]} cells, {meta_df['cell_type'].nunique()} cell types"
            )

        # Create microenvironments file if spatial data is available and requested
        microenvs_file = None
        if (
            params.cellphonedb_use_microenvironments
            and "spatial" in adata_for_analysis.obsm
        ):
            if context:
                await context.info("Creating spatial microenvironments...")
            microenvs_file = await _create_microenvironments_file(
                adata_for_analysis, params, context
            )

        # Set random seed for reproducibility
        debug_seed = (
            params.cellphonedb_debug_seed
            if params.cellphonedb_debug_seed is not None
            else 42
        )
        np.random.seed(debug_seed)

        # Run CellPhoneDB statistical analysis
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save data to temporary files
            counts_file = os.path.join(temp_dir, "counts.txt")
            meta_file = os.path.join(temp_dir, "meta.txt")

            counts_df.to_csv(counts_file, sep="\t")
            meta_df.to_csv(meta_file, sep="\t", index=False)

            try:
                db_path = await _ensure_cellphonedb_database(temp_dir, context)
            except Exception as db_error:
                error_msg = (
                    f"Failed to setup CellPhoneDB database: {str(db_error)}\n\n"
                    "This is required for CellPhoneDB v5 API. Please:\n"
                    "1. Check internet connection for database download\n"
                    "2. Verify CellPhoneDB installation: pip install cellphonedb\n"
                    "3. Ensure write permissions to temporary directory"
                )
                raise RuntimeError(error_msg) from db_error

            # Run the analysis using CellPhoneDB v5 API with correct parameters
            try:
                # STRICT: CellPhoneDB v5 ONLY - no backward compatibility for older versions
                result = cpdb_statistical_analysis_method.call(
                    cpdb_file_path=db_path,  # Fixed: Use actual database path
                    meta_file_path=meta_file,
                    counts_file_path=counts_file,
                    counts_data="hgnc_symbol",  # Improved: Use recommended gene identifier
                    threshold=params.cellphonedb_threshold,
                    result_precision=params.cellphonedb_result_precision,
                    pvalue=params.cellphonedb_pvalue,
                    iterations=params.cellphonedb_iterations,
                    debug_seed=debug_seed,
                    output_path=temp_dir,
                    microenvs_file_path=microenvs_file,
                    score_interactions=True,  # New: Enable interaction scoring (v5 feature)
                )

                # Validate CellPhoneDB v5 format
                if not isinstance(result, dict):
                    raise RuntimeError(
                        f"CellPhoneDB returned unexpected format: {type(result).__name__}. "
                        f"Expected dict from CellPhoneDB v5. Check installation: pip install 'cellphonedb>=5.0.0'"
                    )

                # Extract results from CellPhoneDB v5 dictionary format
                deconvoluted = result.get("deconvoluted")
                means = result.get("means")
                pvalues = result.get("pvalues")
                significant_means = result.get("significant_means")

                # Validate CellPhoneDB API result completeness
                if significant_means is None and "significant_means" not in result:
                    raise RuntimeError(
                        "CellPhoneDB returned incomplete results (missing 'significant_means'). "
                        "Try method='liana' or data_source='raw' for better gene coverage."
                    )
            except Exception as api_error:
                raise RuntimeError(
                    f"CellPhoneDB analysis failed: {str(api_error)}. "
                    f"Consider using method='liana' as alternative."
                ) from api_error

            # Store results in AnnData object
            adata.uns["cellphonedb_deconvoluted"] = deconvoluted
            adata.uns["cellphonedb_means"] = means
            adata.uns["cellphonedb_pvalues"] = pvalues
            adata.uns["cellphonedb_significant_means"] = significant_means

        # Calculate statistics
        n_lr_pairs = (
            len(means) if means is not None and hasattr(means, "__len__") else 0
        )

        # Filter significant pairs based on p-values
        # CellPhoneDB v5 returns all pairs in 'significant_means', so manual filtering is needed
        if (
            pvalues is None
            or not hasattr(pvalues, "values")
            or means is None
            or not hasattr(means, "index")
        ):
            raise RuntimeError(
                "CellPhoneDB p-values unavailable - cannot identify significant interactions. "
                "Try method='liana' as alternative."
            )

        # Filter pairs where ANY cell-cell interaction has p < threshold
        # WITH multiple testing correction for cell type pairs
        threshold = params.cellphonedb_pvalue
        correction_method = params.cellphonedb_correction_method

        # Use nanmin to find minimum p-value across all cell type pairs
        # A pair is significant if its minimum p-value < threshold (after correction)
        # Convert to numeric to handle any non-numeric values
        pval_array = pvalues.select_dtypes(include=[np.number]).values
        if pval_array.shape[0] == 0:
            # No numeric columns found
            raise RuntimeError(
                "CellPhoneDB p-values are not numeric.\n"
                "This indicates a problem with CellPhoneDB analysis.\n"
                "Please check CellPhoneDB installation and data format."
            )

        # Apply multiple testing correction if requested
        # Correct p-values for each L-R pair across its cell type pairs to control FPR
        n_cell_type_pairs = pval_array.shape[1]
        n_lr_pairs_total = pval_array.shape[0]

        if correction_method == "none":
            # No correction: use minimum p-value (not recommended)
            min_pvals = np.nanmin(pval_array, axis=1)
            mask = min_pvals < threshold

            if context:
                await context.warning(
                    f"Multiple testing correction disabled. With {n_cell_type_pairs} cell type pairs, consider using 'fdr_bh' or 'bonferroni'."
                )

            # For 'none', we don't have corrected p-values per se, just use min
            min_pvals_corrected = min_pvals.copy()

        else:
            # CORRECT APPROACH: For each L-R pair, correct its cell type pair p-values
            # Then check if ANY cell type pair remains significant after correction
            from statsmodels.stats.multitest import multipletests

            mask = np.zeros(n_lr_pairs_total, dtype=bool)
            min_pvals_corrected = np.ones(
                n_lr_pairs_total
            )  # Store minimum corrected p-value

            n_uncorrected_sig = 0
            n_corrected_sig = 0

            for i in range(n_lr_pairs_total):
                # Get p-values for this L-R pair across all cell type pairs
                pvals_this_lr = pval_array[i, :]

                # Count uncorrected significance
                n_uncorrected_sig += (pvals_this_lr < threshold).any()

                # Apply correction across cell type pairs for this L-R pair
                reject_this_lr, pvals_corrected_this_lr, _, _ = multipletests(
                    pvals_this_lr,
                    alpha=threshold,
                    method=correction_method,
                    is_sorted=False,
                    returnsorted=False,
                )

                # This L-R pair is significant if ANY cell type pair is significant after correction
                if reject_this_lr.any():
                    mask[i] = True
                    n_corrected_sig += 1

                # Store minimum corrected p-value for this L-R pair
                min_pvals_corrected[i] = pvals_corrected_this_lr.min()

        n_significant_pairs = int(np.sum(mask))

        # Store minimum corrected p-values for transparency
        adata.uns["cellphonedb_pvalues_min_corrected"] = pd.Series(
            min_pvals_corrected,
            index=pvalues.index,
            name=f"min_corrected_pvalue_{correction_method}",
        )

        # Update stored significant_means to match filtered results
        if n_significant_pairs > 0:
            significant_indices = means.index[mask]
            significant_means_filtered = means.loc[significant_indices]

            # Update stored significant_means
            adata.uns["cellphonedb_significant_means"] = significant_means_filtered

            # Also update the variable for downstream use
            significant_means = significant_means_filtered

            # Log filtering results
            if context:
                await context.info(
                    f"CellPhoneDB: {n_significant_pairs}/{n_lr_pairs} significant pairs (p < {threshold}, {correction_method} correction)"
                )
        else:
            # No significant interactions found
            if context:
                await context.warning(
                    f"No significant interactions found at p < {threshold}. Consider adjusting threshold or using method='liana'."
                )

        # Get top LR pairs
        top_lr_pairs = []
        if significant_means is not None and hasattr(significant_means, "head"):
            # CellPhoneDB returns interactions in 'interacting_pair' column
            if (
                hasattr(significant_means, "columns")
                and "interacting_pair" in significant_means.columns
            ):
                top_pairs_df = significant_means.head(params.plot_top_pairs)
                top_lr_pairs = top_pairs_df["interacting_pair"].tolist()

        end_time = time.time()
        analysis_time = end_time - start_time

        if context:
            await context.info(
                f"CellPhoneDB analysis completed in {analysis_time:.2f} seconds"
            )
            await context.info(
                f"Found {n_significant_pairs} significant interactions out of {n_lr_pairs} tested"
            )

        n_cell_types = meta_df["cell_type"].nunique()
        n_cell_type_pairs = n_cell_types**2

        # Add correction statistics (useful for understanding results)
        correction_stats = {}
        if correction_method != "none" and "n_uncorrected_sig" in locals():
            correction_stats["n_uncorrected_significant"] = int(n_uncorrected_sig)
            correction_stats["n_corrected_significant"] = (
                int(n_corrected_sig) if "n_corrected_sig" in locals() else None
            )
            if (
                correction_stats["n_corrected_significant"] is not None
                and n_uncorrected_sig > 0
            ):
                correction_stats["reduction_percentage"] = round(
                    (1 - n_corrected_sig / n_uncorrected_sig) * 100, 2
                )

        statistics = {
            "method": "cellphonedb",
            "iterations": params.cellphonedb_iterations,
            "threshold": params.cellphonedb_threshold,
            "pvalue_threshold": params.cellphonedb_pvalue,
            "n_cell_types": n_cell_types,
            "n_cell_type_pairs": n_cell_type_pairs,
            "multiple_testing_correction": correction_method,
            "microenvironments_used": microenvs_file is not None,
            "analysis_time_seconds": analysis_time,
        }

        # Add correction stats if available
        if correction_stats:
            statistics["correction_statistics"] = correction_stats

        return {
            "n_lr_pairs": n_lr_pairs,
            "n_significant_pairs": n_significant_pairs,
            "top_lr_pairs": top_lr_pairs,
            "cellphonedb_results_key": "cellphonedb_means",
            "cellphonedb_pvalues_key": "cellphonedb_pvalues",
            "cellphonedb_significant_key": "cellphonedb_significant_means",
            "analysis_type": "statistical",
            "statistics": statistics,
        }

    except Exception as e:
        raise RuntimeError(f"CellPhoneDB analysis failed: {str(e)}")


async def _analyze_communication_cellchat_liana(
    adata: Any, params: CellCommunicationParameters, context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze cell communication using CellChat via LIANA"""
    try:
        import liana as li
    except ImportError:
        raise ImportError(
            "LIANA+ is not installed. Please install it with: pip install liana"
        )

    if context:
        await context.info("Running CellChat analysis via LIANA...")

    try:
        import time

        start_time = time.time()

        # Use cell_type_key from params (required field, no auto-detect)
        groupby_col = params.cell_type_key

        validate_obs_column(adata, groupby_col, "Cell type column")

        if context:
            await context.info(
                f"Running CellChat analysis grouped by '{groupby_col}'..."
            )

        # Use CellChatDB resource to match function name
        # Hardcoded to ensure consistency with method name "cellchat_liana"
        resource_name = "cellchatdb"
        if context:
            await context.info(
                f"Using LIANA resource: {resource_name} (CellChatDB database)"
            )

        # Use parameters from user (respect user choice)
        n_perms = params.liana_n_perms

        # Run LIANA with CellChat resource
        li.mt.rank_aggregate(
            adata,
            groupby=groupby_col,
            resource_name=resource_name,
            expr_prop=params.liana_nz_prop,
            min_cells=params.min_cells,
            n_perms=n_perms,
            verbose=False,
            use_raw=True if adata.raw is not None else False,
        )

        # Get results
        liana_res = adata.uns["liana_res"]

        # Calculate statistics using magnitude_rank (signal strength)
        # NOT specificity_rank (which has non-uniform distribution)
        n_lr_pairs = len(liana_res)
        n_significant_pairs = len(liana_res[liana_res["magnitude_rank"] <= 0.05])

        # Get top pairs
        top_lr_pairs = []
        if "magnitude_rank" in liana_res.columns:
            top_pairs_df = liana_res.nsmallest(params.plot_top_pairs, "magnitude_rank")
            top_lr_pairs = [
                f"{row['ligand_complex']}_{row['receptor_complex']}"
                for _, row in top_pairs_df.iterrows()
            ]

        end_time = time.time()
        analysis_time = end_time - start_time

        if context:
            await context.info(
                f"CellChat via LIANA analysis completed in {analysis_time:.2f} seconds"
            )

        statistics = {
            "method": "cellchat_liana",
            "groupby": groupby_col,
            "n_lr_pairs_tested": n_lr_pairs,
            "n_permutations": n_perms,
            "significance_threshold": 0.05,
            "resource": resource_name,  # Actual resource used (cellchatdb)
            "significance_metric": "magnitude_rank",  # Clarify which metric is used
            "analysis_time_seconds": analysis_time,
        }

        return {
            "n_lr_pairs": n_lr_pairs,
            "n_significant_pairs": n_significant_pairs,
            "top_lr_pairs": top_lr_pairs,
            # "liana_results_key": "liana_res",  # Removed to prevent DataFrame serialization overflow
            "analysis_type": "cellchat_via_liana",
            "statistics": statistics,
        }

    except Exception as e:
        raise RuntimeError(f"CellChat via LIANA analysis failed: {str(e)}")


async def _create_microenvironments_file(
    adata: Any, params: CellCommunicationParameters, context: Optional[Context] = None
) -> Optional[str]:
    """Create microenvironments file for CellPhoneDB spatial analysis"""
    try:
        import tempfile

        from sklearn.neighbors import NearestNeighbors

        if "spatial" not in adata.obsm:
            return None

        spatial_coords = adata.obsm["spatial"]

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

        # Create microenvironments using cell types
        validate_obs_column(adata, params.cell_type_key, "Cell type column")

        cell_types = adata.obs[params.cell_type_key].values

        # Create microenvironments by cell type co-occurrence
        microenv_assignments = {}
        microenv_counter = 0

        for i in range(adata.n_obs):
            neighbors = neighbor_matrix[i].indices
            if len(neighbors) > 1:  # At least one neighbor besides itself
                # Get unique cell types in this spatial neighborhood
                neighbor_cell_types = set([cell_types[j] for j in neighbors])

                # Create microenvironment signature based on co-occurring cell types
                microenv_signature = tuple(sorted(neighbor_cell_types))

                if microenv_signature not in microenv_assignments:
                    microenv_assignments[microenv_signature] = (
                        f"microenv_{microenv_counter}"
                    )
                    microenv_counter += 1

        # Generate cell_type -> microenvironment mappings
        cell_type_to_microenv = {}

        for i in range(adata.n_obs):
            neighbors = neighbor_matrix[i].indices
            if len(neighbors) > 1:
                neighbor_cell_types = set([cell_types[j] for j in neighbors])
                microenv_signature = tuple(sorted(neighbor_cell_types))
                microenv_name = microenv_assignments[microenv_signature]

                # Assign all cell types in this neighborhood to this microenvironment
                for ct in neighbor_cell_types:
                    if ct not in cell_type_to_microenv:
                        cell_type_to_microenv[ct] = set()
                    cell_type_to_microenv[ct].add(microenv_name)

        # Create final microenvironments list (cell_type, microenvironment)
        microenvs = []
        for cell_type, microenv_set in cell_type_to_microenv.items():
            for microenv in microenv_set:
                microenvs.append([cell_type, microenv])

        # Save to temporary file with CORRECT format for CellPhoneDB
        temp_file = tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix="_microenvironments.txt"
        )
        temp_file.write("cell_type\tmicroenvironment\n")  # FIXED: Correct header
        for cell_type, microenv in microenvs:
            temp_file.write(
                f"{cell_type}\t{microenv}\n"
            )  # FIXED: cell_type not cell barcode
        temp_file.close()

        if context:
            n_microenvs = len(set([m[1] for m in microenvs]))
            n_cell_types = len(set([m[0] for m in microenvs]))
            await context.info(
                f"Created microenvironments file with {n_microenvs} microenvironments covering {n_cell_types} cell types"
            )

        return temp_file.name

    except Exception as e:
        if context:
            await context.warning(f"Failed to create microenvironments file: {str(e)}")
        return None
