"""
Cell-cell communication analysis tools for spatial transcriptomics data.
"""

from typing import TYPE_CHECKING, Any, Dict, Optional

import numpy as np

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

from ..models.analysis import CellCommunicationResult
from ..models.data import CellCommunicationParameters
from ..utils import validate_obs_column
from ..utils.adata_utils import get_spatial_key, to_dense
from ..utils.dependency_manager import require, validate_r_package
from ..utils.exceptions import (
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


async def _validate_liana_requirements(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> None:
    """Validate LIANA+ requirements"""
    # Spatial connectivity validation
    if params.perform_spatial_analysis and "spatial_connectivities" not in adata.obsp:
        raise DataNotFoundError(
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
    if params.species == "mouse" and params.liana_resource == "consensus":
        await ctx.warning(
            "Using 'consensus' for mouse data. Consider liana_resource='mouseconsensus'."
        )


async def analyze_cell_communication(
    data_id: str,
    ctx: "ToolContext",
    params: CellCommunicationParameters,  # No default - must be provided by caller (LLM)
) -> CellCommunicationResult:
    """Analyze cell-cell communication in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        ctx: ToolContext for data access and logging
        params: Cell communication analysis parameters

    Returns:
        Cell communication analysis result
    """
    # Get data via ToolContext
    adata = await ctx.get_adata(data_id)

    try:
        # Apply method-specific validation
        if params.method == "liana":
            # LIANA-based methods need spatial connectivity validation
            await _validate_liana_requirements(adata, params, ctx)
        elif params.method == "cellphonedb":
            # Check if cell type column exists
            validate_obs_column(adata, params.cell_type_key, "Cell type column")

            # Check for low counts
            n_genes = adata.raw.n_vars if adata.raw is not None else adata.n_vars
            if n_genes < 5000:
                await ctx.warning(
                    f"Gene count ({n_genes}) is relatively low. "
                    f"This may limit the number of interactions found."
                )

            if adata.n_obs < 100:
                await ctx.warning(
                    f"Cell count ({adata.n_obs}) is relatively low. "
                    f"This may affect statistical power."
                )

        # Note: LIANA internally handles use_raw parameter automatically
        # No need for manual data_source switching - consistent with other tools

        # Analyze cell communication using selected method
        if params.method == "liana":
            require("liana", ctx, feature="LIANA+ cell communication analysis")
            result_data = await _analyze_communication_liana(adata, params, ctx)

        elif params.method == "cellphonedb":
            require(
                "cellphonedb", ctx, feature="CellPhoneDB cell communication analysis"
            )
            result_data = await _analyze_communication_cellphonedb(adata, params, ctx)

        elif params.method == "cellchat_r":
            validate_r_package(
                "CellChat",
                ctx,
                install_cmd="devtools::install_github('jinworks/CellChat')",
            )
            result_data = await _analyze_communication_cellchat_r(adata, params, ctx)

        else:
            raise ParameterError(
                f"Unsupported method: {params.method}. "
                f"Supported methods: 'liana', 'cellphonedb', 'cellchat_r'"
            )

        # Note: Results are already stored in adata.uns by the analysis methods
        # Since ctx.get_adata() returns a reference to the stored object,
        # modifications to adata.uns are automatically persisted

        # Store scientific metadata for reproducibility
        from ..utils.adata_utils import store_analysis_metadata

        # Determine database used
        if params.method == "liana":
            database = params.liana_resource
        elif params.method == "cellphonedb":
            database = "cellphonedb"
        elif params.method == "cellchat_liana":
            database = (
                "cellchatdb"  # Match actual LIANA resource name used in implementation
            )
        elif params.method == "cellchat_r":
            database = f"CellChatDB.{params.species}"  # Native R CellChat database
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
        if result_data.get("cellchat_r_results_key"):
            results_keys_dict["uns"].append(result_data["cellchat_r_results_key"])

        # Store metadata
        store_analysis_metadata(
            adata,
            analysis_name=f"cell_communication_{params.method}",
            method=params.method,
            parameters={
                "cell_type_key": params.cell_type_key,
                "n_perms": (params.liana_n_perms if params.method == "liana" else None),
                "nz_prop": (params.liana_nz_prop if params.method == "liana" else None),
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

        return result

    except Exception as e:
        raise ProcessingError(f"Error in cell communication analysis: {str(e)}") from e


async def _analyze_communication_liana(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> Dict[str, Any]:
    """Analyze cell communication using LIANA+"""
    # Use centralized dependency manager for consistent error handling
    require("liana")  # Raises ImportError with install instructions if missing
    import liana as li  # noqa: F401

    try:
        # Ensure spatial connectivity is computed
        if "spatial_connectivities" not in adata.obsp:
            # Use parameters from user or determine optimal bandwidth based on data size
            if params.liana_bandwidth is not None:
                bandwidth = params.liana_bandwidth
            elif adata.n_obs > 3000:
                bandwidth = 300  # Larger bandwidth for large datasets
            else:
                bandwidth = 200  # Standard bandwidth

            # Use Squidpy for spatial neighbor computation
            # Note: Spatial analysis requires spatial neighbors (physical coordinates), not expression neighbors
            # Use centralized dependency manager for consistent error handling
            require(
                "squidpy"
            )  # Raises ImportError with install instructions if missing
            import squidpy as sq

            # Squidpy's spatial_neighbors uses PHYSICAL coordinates
            sq.gr.spatial_neighbors(
                adata,
                coord_type="generic",
                n_neighs=min(30, max(6, adata.n_obs // 100)),  # Adaptive neighbor count
                radius=bandwidth if bandwidth else None,
                delaunay=True,  # Use Delaunay triangulation for spatial data
                set_diag=False,  # Standard practice for spatial graphs
            )

        # Validate species parameter is specified
        if not params.species:
            raise ParameterError(
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
            return _run_liana_cluster_analysis(adata, params, ctx)
        else:
            # Spatial bivariate analysis
            return await _run_liana_spatial_analysis(adata, params, ctx)

    except Exception as e:
        raise ProcessingError(f"LIANA+ analysis failed: {str(e)}") from e


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


def _run_liana_cluster_analysis(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
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
        use_raw=adata.raw is not None,
    )

    # Get results
    liana_res = adata.uns["liana_res"]

    # Calculate statistics using magnitude_rank (signal strength)
    # NOT specificity_rank (which has non-uniform distribution)
    n_lr_pairs = len(liana_res)
    # Use configurable significance threshold (default: 0.05)
    significance_alpha = params.liana_significance_alpha
    n_significant_pairs = len(
        liana_res[liana_res["magnitude_rank"] <= significance_alpha]
    )

    # Get top pairs using vectorized operations (faster than iterrows)
    top_lr_pairs = []
    detected_lr_pairs = []
    if "magnitude_rank" in liana_res.columns:
        top_pairs_df = liana_res.nsmallest(params.plot_top_pairs, "magnitude_rank")
        # Vectorized string concatenation
        ligands = top_pairs_df["ligand_complex"].values
        receptors = top_pairs_df["receptor_complex"].values
        top_lr_pairs = [f"{lig}_{rec}" for lig, rec in zip(ligands, receptors)]
        detected_lr_pairs = list(zip(ligands, receptors))

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
        "significance_threshold": significance_alpha,
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
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
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

    # Count significant pairs using statistical significance (p-values with FDR correction)
    #
    # P-values are ALWAYS available because:
    # 1. We always pass global_name (params.liana_global_metric, default: "morans")
    # 2. We always pass n_perms > 0 (params.liana_n_perms, default: 1000, Field(gt=0))
    # 3. LIANA computes p-values via permutation test when n_perms > 0
    #    (see liana/method/sp/_bivariate/_global_functions.py lines 104-128)
    from statsmodels.stats.multitest import multipletests

    pvals_col = f"{global_metric}_pvals"
    alpha = params.liana_significance_alpha
    pvals = lrdata.var[pvals_col]

    reject, pvals_corrected, _, _ = multipletests(
        pvals, alpha=alpha, method="fdr_bh"  # Benjamini-Hochberg FDR correction
    )

    n_significant_pairs = reject.sum()

    # Store corrected p-values and significance flags for downstream use
    lrdata.var[f"{pvals_col}_corrected"] = pvals_corrected
    lrdata.var[f"{global_metric}_significant"] = reject

    # Store results in adata
    adata.uns["liana_spatial_res"] = lrdata.var
    adata.obsm["liana_spatial_scores"] = to_dense(lrdata.X)
    adata.uns["liana_spatial_interactions"] = lrdata.var.index.tolist()

    if "pvals" in lrdata.layers:
        adata.obsm["liana_spatial_pvals"] = to_dense(lrdata.layers["pvals"])

    if "cats" in lrdata.layers:
        adata.obsm["liana_spatial_cats"] = to_dense(lrdata.layers["cats"])

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
        "significance_method": (
            "FDR-corrected p-values"
            if pvals_col in lrdata.var.columns
            else "threshold-based (deprecated)"
        ),
        "fdr_method": "Benjamini-Hochberg" if pvals_col in lrdata.var.columns else None,
        "alpha": alpha if pvals_col in lrdata.var.columns else None,
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


async def _ensure_cellphonedb_database(output_dir: str, ctx: "ToolContext") -> str:
    """Ensure CellPhoneDB database is available, download if not exists"""
    # Use centralized dependency manager for consistent error handling
    require("cellphonedb")  # Raises ImportError with install instructions if missing
    import os

    from cellphonedb.utils import db_utils

    # Check if database file already exists
    db_path = os.path.join(output_dir, "cellphonedb.zip")

    if os.path.exists(db_path):
        return db_path

    try:
        # Download latest database
        db_utils.download_database(output_dir, "v5.0.0")

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
        raise DependencyError(error_msg) from e


async def _analyze_communication_cellphonedb(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> Dict[str, Any]:
    """Analyze cell communication using CellPhoneDB"""
    # Use centralized dependency manager for consistent error handling
    require("cellphonedb")  # Raises ImportError with install instructions if missing
    import os
    import tempfile

    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

    try:
        import time

        start_time = time.time()

        # Use cell_type_key from params (required field, no auto-detect)
        cell_type_col = params.cell_type_key

        validate_obs_column(adata, cell_type_col, "Cell type column")

        # Use original adata directly (no gene filtering needed)
        adata_for_analysis = adata

        # Import pandas for DataFrame operations
        import csv

        import pandas as pd
        import scipy.sparse as sp

        # Check if data is sparse (used for efficient matrix access)
        is_sparse = sp.issparse(adata_for_analysis.X)

        # Prepare meta data
        meta_df = pd.DataFrame(
            {
                "Cell": adata_for_analysis.obs.index,
                "cell_type": adata_for_analysis.obs[cell_type_col].astype(str),
            }
        )

        # Create microenvironments file if spatial data is available and requested
        microenvs_file = None
        if (
            params.cellphonedb_use_microenvironments
            and "spatial" in adata_for_analysis.obsm
        ):
            microenvs_file = await _create_microenvironments_file(
                adata_for_analysis, params, ctx
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

            # Direct file writing: Stream sparse matrix to CSV without creating DataFrame
            # Memory-efficient approach: write gene-by-gene instead of toarray()
            with open(counts_file, "w", newline="") as f:
                writer = csv.writer(f, delimiter="\t")

                # Write header: empty first column + cell names
                header = [""] + list(adata_for_analysis.obs_names)
                writer.writerow(header)

                # Convert to CSC for efficient column access (genes)
                if is_sparse:
                    X_csc = adata_for_analysis.X.tocsc()
                else:
                    X_csc = adata_for_analysis.X

                # Write gene-by-gene (memory constant)
                for i, gene_name in enumerate(adata_for_analysis.var_names):
                    gene_expression = to_dense(X_csc[:, i]).flatten()
                    writer.writerow([gene_name] + list(gene_expression))

            meta_df.to_csv(meta_file, sep="\t", index=False)

            try:
                db_path = await _ensure_cellphonedb_database(temp_dir, ctx)
            except Exception as db_error:
                raise DependencyError(
                    f"CellPhoneDB database setup failed: {db_error}"
                ) from db_error

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
                    score_interactions=False,  # Disabled: CellPhoneDB v5 scoring has bugs
                )
            except KeyError as key_error:
                raise ProcessingError(
                    f"CellPhoneDB found no L-R interactions. "
                    f"CellPhoneDB is human-only; use method='liana' for mouse data. "
                    f"Error: {key_error}"
                ) from key_error
            except Exception as api_error:
                raise ProcessingError(
                    f"CellPhoneDB analysis failed: {str(api_error)}. "
                    f"Consider using method='liana' as alternative."
                ) from api_error

            # Validate CellPhoneDB v5 format
            if not isinstance(result, dict):
                raise ProcessingError(
                    f"CellPhoneDB returned unexpected format: {type(result).__name__}. "
                    f"Expected dict from CellPhoneDB v5. Check installation: pip install 'cellphonedb>=5.0.0'"
                )

            # Check for empty results (no interactions found)
            if not result or "significant_means" not in result:
                raise DataNotFoundError(
                    "CellPhoneDB found no L-R interactions. "
                    "CellPhoneDB is human-only; use method='liana' for mouse data."
                )

            # Extract results from CellPhoneDB v5 dictionary format
            deconvoluted = result.get("deconvoluted")
            means = result.get("means")
            pvalues = result.get("pvalues")
            significant_means = result.get("significant_means")

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
            raise DataNotFoundError(
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
            raise ProcessingError("CellPhoneDB p-values are not numeric.")

        # Apply multiple testing correction if requested
        # Correct p-values for each L-R pair across its cell type pairs to control FPR
        n_cell_type_pairs = pval_array.shape[1]
        n_lr_pairs_total = pval_array.shape[0]

        if correction_method == "none":
            # No correction: use minimum p-value (not recommended)
            min_pvals = np.nanmin(pval_array, axis=1)
            mask = min_pvals < threshold

            await ctx.warning(
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
        # Convert Series to DataFrame for H5AD compatibility (H5AD cannot store pd.Series)
        adata.uns["cellphonedb_pvalues_min_corrected"] = pd.DataFrame(
            {f"min_corrected_pvalue_{correction_method}": min_pvals_corrected},
            index=pvalues.index.astype(str),
        )

        # Update stored significant_means to match filtered results
        if n_significant_pairs > 0:
            significant_indices = means.index[mask]
            significant_means_filtered = means.loc[significant_indices]

            # Update stored significant_means
            adata.uns["cellphonedb_significant_means"] = significant_means_filtered

            # Also update the variable for downstream use
            significant_means = significant_means_filtered
        else:
            # No significant interactions found
            await ctx.warning(
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

        n_cell_types = meta_df["cell_type"].nunique()
        n_cell_type_pairs = n_cell_types**2

        # Add correction statistics (useful for understanding results)
        # When correction_method != "none", n_uncorrected_sig and n_corrected_sig
        # are always defined in the else branch above (lines 1008-1009)
        correction_stats = {}
        if correction_method != "none":
            correction_stats["n_uncorrected_significant"] = int(n_uncorrected_sig)
            correction_stats["n_corrected_significant"] = int(n_corrected_sig)
            if n_uncorrected_sig > 0:
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
        raise ProcessingError(f"CellPhoneDB analysis failed: {str(e)}") from e
    finally:
        # Cleanup: Remove temporary microenvironments file if created
        if microenvs_file is not None:
            try:
                os.remove(microenvs_file)
            except OSError:
                pass  # Cleanup failure is not critical


async def _create_microenvironments_file(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> Optional[str]:
    """Create microenvironments file for CellPhoneDB spatial analysis"""
    try:
        import tempfile

        from sklearn.neighbors import NearestNeighbors

        spatial_key = get_spatial_key(adata)
        if spatial_key is None:
            return None

        spatial_coords = adata.obsm[spatial_key]

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

        # Find spatial neighbors for each cell
        nn = NearestNeighbors(radius=radius)
        nn.fit(spatial_coords)
        neighbor_matrix = nn.radius_neighbors_graph(spatial_coords)

        # Create microenvironments using cell types
        validate_obs_column(adata, params.cell_type_key, "Cell type column")

        cell_types = adata.obs[params.cell_type_key].values

        # Create microenvironments by cell type co-occurrence
        # Optimized: Single loop to build both mappings (2x faster)
        microenv_assignments = {}
        cell_type_to_microenv = {}
        microenv_counter = 0

        for i in range(adata.n_obs):
            neighbors = neighbor_matrix[i].indices
            if len(neighbors) > 1:  # At least one neighbor besides itself
                # Get unique cell types in this spatial neighborhood (computed once)
                neighbor_cell_types = set(cell_types[j] for j in neighbors)

                # Create microenvironment signature based on co-occurring cell types
                microenv_signature = tuple(sorted(neighbor_cell_types))

                # First use: create assignment if new signature
                if microenv_signature not in microenv_assignments:
                    microenv_assignments[microenv_signature] = (
                        f"microenv_{microenv_counter}"
                    )
                    microenv_counter += 1

                # Second use: update cell_type_to_microenv mappings
                microenv_name = microenv_assignments[microenv_signature]
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

        return temp_file.name

    except Exception as e:
        await ctx.warning(f"Failed to create microenvironments file: {str(e)}")
        return None


async def _analyze_communication_cellchat_r(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> Dict[str, Any]:
    """Analyze cell communication using native R CellChat package

    This implementation uses rpy2 to call the original R CellChat package,
    which includes full features like mediator proteins and signaling pathways
    that are not available in the LIANA simplified implementation.

    Args:
        adata: AnnData object with expression data
        params: Cell communication analysis parameters
        ctx: ToolContext for logging and data access

    Returns:
        Dictionary with analysis results
    """
    import pandas as pd
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter

    try:
        import time

        start_time = time.time()

        # Validate cell type column
        validate_obs_column(adata, params.cell_type_key, "Cell type column")

        # Check for spatial data
        spatial_key = get_spatial_key(adata)
        has_spatial = spatial_key is not None

        # Prepare expression matrix (genes x cells, normalized)
        # CellChat requires normalized data with comprehensive gene coverage
        # Use adata.raw if available (contains all genes before HVG filtering)
        if adata.raw is not None:
            data_source = adata.raw
        else:
            data_source = adata

        # Run CellChat in R - start early to get gene list for pre-filtering
        with localconverter(
            ro.default_converter + pandas2ri.converter + numpy2ri.converter
        ):
            # Load CellChat
            ro.r("library(CellChat)")

            # Set species-specific database
            species_db_map = {
                "human": "CellChatDB.human",
                "mouse": "CellChatDB.mouse",
                "zebrafish": "CellChatDB.zebrafish",
            }
            db_name = species_db_map.get(params.species, "CellChatDB.human")

            # Memory optimization: Get CellChatDB gene list and pre-filter
            # This reduces memory from O(n_cells × n_all_genes) to O(n_cells × n_db_genes)
            # Typical savings: 20000 genes → 1500 genes = 13x memory reduction
            ro.r(
                f"""
                CellChatDB <- {db_name}
                # Get all genes used in CellChatDB (ligands, receptors, cofactors)
                cellchat_genes <- unique(c(
                    CellChatDB$geneInfo$Symbol,
                    unlist(strsplit(CellChatDB$interaction$ligand, "_")),
                    unlist(strsplit(CellChatDB$interaction$receptor, "_"))
                ))
                cellchat_genes <- cellchat_genes[!is.na(cellchat_genes)]
            """
            )
            cellchat_genes_r = ro.r("cellchat_genes")
            cellchat_genes = set(cellchat_genes_r)

            # Filter to genes present in both data and CellChatDB
            common_genes = data_source.var_names.intersection(cellchat_genes)

            if len(common_genes) == 0:
                raise ValueError(
                    f"No genes overlap between data and {db_name}. "
                    f"Check if species parameter matches your data."
                )

            # Create expression matrix with only CellChatDB genes (memory efficient)
            gene_indices = [data_source.var_names.get_loc(g) for g in common_genes]
            expr_matrix = pd.DataFrame(
                to_dense(data_source.X[:, gene_indices]).T,
                index=common_genes,
                columns=adata.obs_names,
            )

            # Prepare metadata
            # CellChat doesn't allow labels starting with '0', so add prefix for numeric
            cell_labels = adata.obs[params.cell_type_key].astype(str).values
            # Check if any label is '0' or starts with a digit - add 'cluster_' prefix
            if any(
                label == "0" or (label and label[0].isdigit()) for label in cell_labels
            ):
                cell_labels = [f"cluster_{label}" for label in cell_labels]
            meta_df = pd.DataFrame(
                {"labels": cell_labels},
                index=adata.obs_names,
            )

            # Prepare spatial coordinates if available
            spatial_locs = None
            if has_spatial and params.cellchat_distance_use:
                spatial_coords = adata.obsm[spatial_key]
                spatial_locs = pd.DataFrame(
                    spatial_coords[:, :2],
                    index=adata.obs_names,
                    columns=["x", "y"],
                )

            # Transfer data to R
            ro.globalenv["expr_matrix"] = expr_matrix
            ro.globalenv["meta_df"] = meta_df

            # Create CellChat object (db_name already set during gene pre-filtering)
            if (
                has_spatial
                and params.cellchat_distance_use
                and spatial_locs is not None
            ):
                # Spatial mode
                ro.globalenv["spatial_locs"] = spatial_locs

                # CellChat v2 requires spatial.factors with 'ratio' and 'tol':
                # - ratio: conversion factor from pixels to micrometers (um)
                # - tol: tolerance factor (half of spot/cell size in um)
                # Use user-configurable parameters for platform flexibility
                pixel_ratio = params.cellchat_pixel_ratio
                spatial_tol = params.cellchat_spatial_tol
                ro.globalenv["pixel_ratio"] = pixel_ratio
                ro.globalenv["spatial_tol"] = spatial_tol
                ro.r(
                    """
                    spatial.factors <- data.frame(
                        ratio = pixel_ratio,
                        tol = spatial_tol
                    )

                    cellchat <- createCellChat(
                        object = as.matrix(expr_matrix),
                        meta = meta_df,
                        group.by = "labels",
                        datatype = "spatial",
                        coordinates = as.matrix(spatial_locs),
                        spatial.factors = spatial.factors
                    )
                """
                )
            else:
                # Non-spatial mode
                ro.r(
                    """
                    cellchat <- createCellChat(
                        object = as.matrix(expr_matrix),
                        meta = meta_df,
                        group.by = "labels"
                    )
                """
                )

            # Set database
            ro.r(
                f"""
                CellChatDB <- {db_name}
            """
            )

            # Subset database by category if specified
            if params.cellchat_db_category != "All":
                ro.r(
                    f"""
                    CellChatDB.use <- subsetDB(
                        CellChatDB,
                        search = "{params.cellchat_db_category}"
                    )
                    cellchat@DB <- CellChatDB.use
                """
                )
            else:
                ro.r(
                    """
                    cellchat@DB <- CellChatDB
                """
                )

            # Preprocessing
            ro.r(
                """
                cellchat <- subsetData(cellchat)
                cellchat <- identifyOverExpressedGenes(cellchat)
                cellchat <- identifyOverExpressedInteractions(cellchat)
            """
            )

            # Project data (optional but recommended)
            ro.r(
                """
                # Project data onto PPI network (optional)
                tryCatch({
                    cellchat <- projectData(cellchat, PPI.human)
                }, error = function(e) {
                    message("Skipping data projection: ", e$message)
                })
            """
            )

            # Compute communication probability
            if has_spatial and params.cellchat_distance_use:
                # Spatial mode with distance constraints
                # CellChat v2 requires either contact.range or contact.knn.k
                if params.cellchat_contact_range is not None:
                    contact_param = f"contact.range = {params.cellchat_contact_range}"
                else:
                    contact_param = f"contact.knn.k = {params.cellchat_contact_knn_k}"

                ro.r(
                    f"""
                    cellchat <- computeCommunProb(
                        cellchat,
                        type = "{params.cellchat_type}",
                        trim = {params.cellchat_trim},
                        population.size = {str(params.cellchat_population_size).upper()},
                        distance.use = TRUE,
                        interaction.range = {params.cellchat_interaction_range},
                        scale.distance = {params.cellchat_scale_distance},
                        {contact_param}
                    )
                """
                )
            else:
                # Non-spatial mode
                ro.r(
                    f"""
                    cellchat <- computeCommunProb(
                        cellchat,
                        type = "{params.cellchat_type}",
                        trim = {params.cellchat_trim},
                        population.size = {str(params.cellchat_population_size).upper()}
                    )
                """
                )

            # Filter communication
            ro.r(
                f"""
                cellchat <- filterCommunication(cellchat, min.cells = {params.cellchat_min_cells})
            """
            )

            # Compute pathway-level communication
            ro.r(
                """
                cellchat <- computeCommunProbPathway(cellchat)
            """
            )

            # Aggregate network
            ro.r(
                """
                cellchat <- aggregateNet(cellchat)
            """
            )

            # Extract results
            ro.r(
                """
                # Get LR pairs
                lr_pairs <- cellchat@LR$LRsig

                # Get communication probabilities
                net <- cellchat@net

                # Get pathway-level probabilities
                netP <- cellchat@netP

                # Count interactions
                n_lr_pairs <- length(unique(lr_pairs$interaction_name))

                # Get significant pairs (probability > 0)
                prob_matrix <- net$prob
                n_significant <- sum(prob_matrix > 0, na.rm = TRUE)

                # Get top pathways
                pathway_names <- rownames(netP$prob)
                if (length(pathway_names) > 0) {
                    # Sum probabilities across cell type pairs for each pathway
                    pathway_sums <- rowSums(netP$prob, na.rm = TRUE)
                    top_pathway_idx <- order(pathway_sums, decreasing = TRUE)[1:min(10, length(pathway_names))]
                    top_pathways <- pathway_names[top_pathway_idx]
                } else {
                    top_pathways <- character(0)
                }

                # Get top LR pairs
                if (nrow(lr_pairs) > 0) {
                    top_lr <- head(lr_pairs$interaction_name, 10)
                } else {
                    top_lr <- character(0)
                }
            """
            )

            # Convert results back to Python
            n_lr_pairs = int(ro.r("n_lr_pairs")[0])
            n_significant_pairs = int(ro.r("n_significant")[0])
            top_pathways = list(ro.r("top_pathways"))
            top_lr_pairs = list(ro.r("top_lr"))

            # Get full results for storage
            lr_pairs_df = ro.r("as.data.frame(lr_pairs)")
            prob_matrix = ro.r("as.matrix(net$prob)")
            pval_matrix = ro.r("as.matrix(net$pval)")

            # Store in adata
            adata.uns["cellchat_r_lr_pairs"] = pd.DataFrame(lr_pairs_df)
            adata.uns["cellchat_r_prob"] = np.array(prob_matrix)
            adata.uns["cellchat_r_pval"] = np.array(pval_matrix)
            adata.uns["cellchat_r_top_pathways"] = top_pathways
            adata.uns["cellchat_r_params"] = {
                "species": params.species,
                "db_category": params.cellchat_db_category,
                "type": params.cellchat_type,
                "distance_use": params.cellchat_distance_use if has_spatial else False,
            }

            # Store detected LR pairs in standardized format for visualization
            detected_lr_pairs = []
            for pair_str in top_lr_pairs:
                if "_" in pair_str:
                    parts = pair_str.split("_", 1)
                    if len(parts) == 2:
                        detected_lr_pairs.append((parts[0], parts[1]))

            adata.uns["detected_lr_pairs"] = detected_lr_pairs
            adata.uns["cell_communication_results"] = {
                "top_lr_pairs": top_lr_pairs,
                "top_pathways": top_pathways,
                "method": "cellchat_r",
                "n_pairs": len(top_lr_pairs),
                "species": params.species,
            }

        end_time = time.time()
        analysis_time = end_time - start_time

        statistics = {
            "method": "cellchat_r",
            "species": params.species,
            "db_category": params.cellchat_db_category,
            "aggregation_type": params.cellchat_type,
            "trim": params.cellchat_trim,
            "population_size": params.cellchat_population_size,
            "min_cells": params.cellchat_min_cells,
            "spatial_mode": has_spatial and params.cellchat_distance_use,
            "n_lr_pairs_tested": n_lr_pairs,
            "analysis_time_seconds": analysis_time,
            "top_pathways": top_pathways[:5] if top_pathways else [],
        }

        return {
            "n_lr_pairs": n_lr_pairs,
            "n_significant_pairs": n_significant_pairs,
            "top_lr_pairs": top_lr_pairs,
            "cellchat_r_results_key": "cellchat_r_lr_pairs",
            "cellchat_r_prob_key": "cellchat_r_prob",
            "cellchat_r_pval_key": "cellchat_r_pval",
            "analysis_type": "cellchat_native",
            "statistics": statistics,
        }

    except Exception as e:
        raise ProcessingError(f"CellChat R analysis failed: {str(e)}") from e
