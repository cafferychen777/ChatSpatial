"""
Cell-cell communication analysis tools for spatial transcriptomics data.

Architecture:
    All CCC results are stored in a unified structure at adata.uns["ccc"].
    This ensures consistency across methods and simplifies visualization.

Storage Structure:
    adata.uns["ccc"] = {
        "method": str,           # "liana", "cellphonedb", "cellchat_r", "fastccc"
        "analysis_type": str,    # "cluster" or "spatial"
        "species": str,          # "human", "mouse", "zebrafish"
        "database": str,         # Resource/database used
        "lr_pairs": list[str],   # Standardized format: "LIGAND_RECEPTOR"
        "top_lr_pairs": list[str],
        "n_pairs": int,
        "n_significant": int,
        "results": DataFrame,    # Main results (standardized)
        "pvalues": DataFrame,    # P-values (optional)
        "autocrine": {           # Autocrine loop detection results
            "n_loops": int,
            "top_pairs": list[str],
            "results": DataFrame | None,
        },
        "statistics": dict,      # Method-specific statistics
        "method_data": dict,     # Method-specific raw data
    }

    For spatial analysis, additional data in adata.obsm:
        "ccc_spatial_scores": ndarray  # (n_spots, n_pairs)
        "ccc_spatial_pvals": ndarray   # (n_spots, n_pairs) optional
"""

from collections import defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Optional

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

from ..models.analysis import CellCommunicationResult
from ..models.data import CellCommunicationParameters
from ..utils import validate_obs_column
from ..utils.adata_utils import get_raw_data_source, get_spatial_key, to_dense
from ..utils.dependency_manager import require, validate_r_package
from ..utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)

# =============================================================================
# Unified CCC Storage Structure
# =============================================================================

# Storage keys (single source of truth)
CCC_UNS_KEY = "ccc"  # Main storage in adata.uns
CCC_SPATIAL_SCORES_KEY = "ccc_spatial_scores"  # Spatial scores in adata.obsm
CCC_SPATIAL_PVALS_KEY = "ccc_spatial_pvals"  # Spatial p-values in adata.obsm


@dataclass
class CCCAutocrine:
    """Autocrine loop detection results."""

    n_loops: int = 0
    top_pairs: list[str] = field(default_factory=list)
    results: Optional[pd.DataFrame] = None


@dataclass
class CCCStorage:
    """Unified storage structure for cell-cell communication results.

    This dataclass defines the canonical structure for all CCC results,
    regardless of the analysis method used. It ensures consistency across
    LIANA, CellPhoneDB, CellChat R, and FastCCC.

    All LR pairs are stored in standardized format: "LIGAND_RECEPTOR"
    """

    # === Core metadata ===
    method: str  # "liana", "cellphonedb", "cellchat_r", "fastccc"
    analysis_type: str  # "cluster" or "spatial"
    species: str  # "human", "mouse", "zebrafish"
    database: str  # Resource/database used

    # === Standardized LR pairs (format: "LIGAND_RECEPTOR") ===
    lr_pairs: list[str] = field(default_factory=list)
    top_lr_pairs: list[str] = field(default_factory=list)
    n_pairs: int = 0
    n_significant: int = 0

    # === Results (standardized DataFrames) ===
    results: Optional[pd.DataFrame] = None  # Main results
    pvalues: Optional[pd.DataFrame] = None  # P-values matrix

    # === Autocrine analysis ===
    autocrine: CCCAutocrine = field(default_factory=CCCAutocrine)

    # === Statistics and method-specific data ===
    statistics: dict = field(default_factory=dict)
    method_data: dict = field(default_factory=dict)  # Raw method-specific data

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for storage in adata.uns."""
        return {
            "method": self.method,
            "analysis_type": self.analysis_type,
            "species": self.species,
            "database": self.database,
            "lr_pairs": self.lr_pairs,
            "top_lr_pairs": self.top_lr_pairs,
            "n_pairs": self.n_pairs,
            "n_significant": self.n_significant,
            "results": self.results,
            "pvalues": self.pvalues,
            "autocrine": {
                "n_loops": self.autocrine.n_loops,
                "top_pairs": self.autocrine.top_pairs,
                "results": self.autocrine.results,
            },
            "statistics": self.statistics,
            "method_data": self.method_data,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "CCCStorage":
        """Create from dictionary stored in adata.uns."""
        autocrine_data = data.get("autocrine", {})
        return cls(
            method=data["method"],
            analysis_type=data["analysis_type"],
            species=data["species"],
            database=data["database"],
            lr_pairs=data.get("lr_pairs", []),
            top_lr_pairs=data.get("top_lr_pairs", []),
            n_pairs=data.get("n_pairs", 0),
            n_significant=data.get("n_significant", 0),
            results=data.get("results"),
            pvalues=data.get("pvalues"),
            autocrine=CCCAutocrine(
                n_loops=autocrine_data.get("n_loops", 0),
                top_pairs=autocrine_data.get("top_pairs", []),
                results=autocrine_data.get("results"),
            ),
            statistics=data.get("statistics", {}),
            method_data=data.get("method_data", {}),
        )


def standardize_lr_pair(pair_str: str) -> str:
    """Standardize LR pair format to 'LIGAND_RECEPTOR'.

    Handles various input formats:
    - 'LIGAND^RECEPTOR' -> 'LIGAND_RECEPTOR'
    - 'LIGAND_RECEPTOR' -> 'LIGAND_RECEPTOR' (no change)
    - 'ligand-receptor' -> 'LIGAND_RECEPTOR'
    """
    # Replace common separators with underscore
    standardized = pair_str.replace("^", "_").replace("-", "_")
    # Don't uppercase - preserve original gene names
    return standardized


def store_ccc_results(adata: Any, storage: CCCStorage) -> None:
    """Store CCC results in unified structure.

    Args:
        adata: AnnData object
        storage: CCCStorage dataclass with results
    """
    adata.uns[CCC_UNS_KEY] = storage.to_dict()


def get_ccc_results(adata: Any) -> Optional[CCCStorage]:
    """Retrieve CCC results from unified structure.

    Args:
        adata: AnnData object

    Returns:
        CCCStorage if results exist, None otherwise
    """
    if CCC_UNS_KEY not in adata.uns:
        return None
    return CCCStorage.from_dict(adata.uns[CCC_UNS_KEY])


def has_ccc_results(adata: Any) -> bool:
    """Check if CCC results exist."""
    return CCC_UNS_KEY in adata.uns


async def analyze_cell_communication(
    data_id: str,
    ctx: "ToolContext",
    params: CellCommunicationParameters,  # No default - must be provided by caller (LLM)
) -> CellCommunicationResult:
    """Analyze cell-cell communication in spatial transcriptomics data.

    All results are stored in a unified structure at adata.uns["ccc"].
    See module docstring for storage structure details.

    Args:
        data_id: Dataset ID
        ctx: ToolContext for data access and logging
        params: Cell communication analysis parameters

    Returns:
        CellCommunicationResult with analysis summary
    """
    adata = await ctx.get_adata(data_id)

    try:
        # === Method-specific validation ===
        await _validate_ccc_params(adata, params, ctx)

        # === Run analysis (returns CCCStorage) ===
        storage = await _run_ccc_analysis(adata, params, ctx)

        # === Extract autocrine loops and integrate into storage ===
        _integrate_autocrine_detection(storage, params.plot_top_pairs)

        # === Store results in unified structure ===
        store_ccc_results(adata, storage)

        # === Store spatial data in obsm if applicable ===
        if (
            storage.analysis_type == "spatial"
            and "spatial_scores" in storage.method_data
        ):
            adata.obsm[CCC_SPATIAL_SCORES_KEY] = storage.method_data["spatial_scores"]
            if "spatial_pvals" in storage.method_data:
                adata.obsm[CCC_SPATIAL_PVALS_KEY] = storage.method_data["spatial_pvals"]

        # === Store scientific metadata for reproducibility ===
        from ..utils.adata_utils import store_analysis_metadata
        from ..utils.results_export import export_analysis_result

        results_keys: dict[str, list[str]] = {
            "obs": [],
            "obsm": (
                [CCC_SPATIAL_SCORES_KEY] if storage.analysis_type == "spatial" else []
            ),
            "uns": [CCC_UNS_KEY],
        }

        store_analysis_metadata(
            adata,
            analysis_name=f"cell_communication_{params.method}",
            method=params.method,
            parameters={
                "cell_type_key": params.cell_type_key,
                "min_cells": params.min_cells,
            },
            results_keys=results_keys,
            statistics={
                "n_lr_pairs": storage.n_pairs,
                "n_significant_pairs": storage.n_significant,
                "analysis_type": storage.analysis_type,
            },
            species=params.species,
            database=storage.database,
        )

        export_analysis_result(adata, data_id, f"cell_communication_{params.method}")

        # === Create MCP response ===
        return CellCommunicationResult(
            data_id=data_id,
            method=storage.method,
            species=storage.species,
            database=storage.database,
            analysis_type=storage.analysis_type,
            n_lr_pairs=storage.n_pairs,
            n_significant_pairs=storage.n_significant,
            top_lr_pairs=storage.top_lr_pairs,
            n_autocrine_loops=storage.autocrine.n_loops,
            top_autocrine_loops=storage.autocrine.top_pairs,
            results_key=CCC_UNS_KEY,
            statistics=storage.statistics,
        )

    except Exception as e:
        raise ProcessingError(f"Error in cell communication analysis: {e}") from e


async def _validate_ccc_params(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> None:
    """Validate CCC parameters before analysis."""
    if params.method == "liana":
        if (
            params.perform_spatial_analysis
            and "spatial_connectivities" not in adata.obsp
        ):
            raise DataNotFoundError(
                "Spatial connectivity required. Run sq.gr.spatial_neighbors() first."
            )
        validate_obs_column(adata, params.cell_type_key, "Cell type")
        if params.species == "mouse" and params.liana_resource == "consensus":
            await ctx.warning(
                "Using 'consensus' for mouse data. Consider liana_resource='mouseconsensus'."
            )

    elif params.method == "cellphonedb":
        validate_obs_column(adata, params.cell_type_key, "Cell type")
        raw_result = get_raw_data_source(adata, prefer_complete_genes=True)
        if len(raw_result.var_names) < 5000:
            await ctx.warning(
                f"Gene count ({len(raw_result.var_names)}) is relatively low. "
                "This may limit the number of interactions found."
            )
        if adata.n_obs < 100:
            await ctx.warning(
                f"Cell count ({adata.n_obs}) is relatively low. "
                "This may affect statistical power."
            )

    elif params.method in ("cellchat_r", "fastccc"):
        validate_obs_column(adata, params.cell_type_key, "Cell type")


async def _run_ccc_analysis(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> CCCStorage:
    """Run CCC analysis and return unified storage structure."""
    if params.method == "liana":
        require("liana", ctx, feature="LIANA+ cell communication analysis")
        return await _analyze_communication_liana(adata, params, ctx)

    elif params.method == "cellphonedb":
        require("cellphonedb", ctx, feature="CellPhoneDB cell communication analysis")
        return await _analyze_communication_cellphonedb(adata, params, ctx)

    elif params.method == "cellchat_r":
        validate_r_package(
            "CellChat", ctx, install_cmd="devtools::install_github('jinworks/CellChat')"
        )
        return _analyze_communication_cellchat_r(adata, params, ctx)

    elif params.method == "fastccc":
        require("fastccc", ctx, feature="FastCCC cell communication analysis")
        return await _analyze_communication_fastccc(adata, params, ctx)

    else:
        raise ParameterError(
            f"Unsupported method: {params.method}. "
            "Supported methods: 'liana', 'cellphonedb', 'cellchat_r', 'fastccc'"
        )


def _integrate_autocrine_detection(storage: CCCStorage, n_top: int) -> None:
    """Extract autocrine loops and integrate into CCCStorage.

    Autocrine signaling: source == target cell type.
    This is integrated directly into the storage structure.
    """
    if storage.analysis_type == "spatial":
        # Spatial analysis doesn't have cell type pairs
        return

    if storage.results is None or len(storage.results) == 0:
        return

    results = storage.results

    # Method-specific autocrine extraction
    if storage.method == "liana" and "source" in results.columns:
        # LIANA cluster: direct source == target filter
        autocrine_df = results[results["source"] == results["target"]].copy()
        if len(autocrine_df) > 0:
            autocrine_df["lr_pair"] = (
                autocrine_df["ligand_complex"] + "_" + autocrine_df["receptor_complex"]
            )
            if "magnitude_rank" in autocrine_df.columns:
                autocrine_df = autocrine_df.nsmallest(n_top, "magnitude_rank")
            storage.autocrine = CCCAutocrine(
                n_loops=len(autocrine_df["lr_pair"].unique()),
                top_pairs=autocrine_df["lr_pair"].head(n_top).tolist(),
                results=autocrine_df,
            )

    elif storage.method in ("cellphonedb", "fastccc"):
        # Matrix format: columns are "celltype1|celltype2"
        autocrine_cols = [
            col
            for col in results.columns
            if "|" in str(col) and str(col).split("|")[0] == str(col).split("|")[1]
        ]
        if autocrine_cols:
            # Filter rows with any autocrine interaction
            numeric_cols = (
                results[autocrine_cols].select_dtypes(include=[np.number]).columns
            )
            if len(numeric_cols) > 0:
                mask = (results[numeric_cols] > 0).any(axis=1)
                autocrine_df = results[mask].copy()
                if len(autocrine_df) > 0:
                    lr_pairs = autocrine_df.index.tolist()[:n_top]
                    storage.autocrine = CCCAutocrine(
                        n_loops=len(autocrine_df),
                        top_pairs=[standardize_lr_pair(p) for p in lr_pairs],
                        results=autocrine_df,
                    )

    elif storage.method == "cellchat_r" and "prob_matrix" in storage.method_data:
        # CellChat: extract diagonal from 3D probability matrix
        prob_matrix = storage.method_data["prob_matrix"]
        if prob_matrix is not None and len(prob_matrix.shape) == 3:
            n_cell_types = prob_matrix.shape[0]
            # Sum diagonal probabilities
            autocrine_probs = np.sum(
                [prob_matrix[i, i, :] for i in range(n_cell_types)], axis=0
            )
            autocrine_mask = autocrine_probs > 0
            if autocrine_mask.any() and "interaction_name" in results.columns:
                autocrine_pairs = results[autocrine_mask]["interaction_name"].tolist()
                storage.autocrine = CCCAutocrine(
                    n_loops=len(autocrine_pairs),
                    top_pairs=[standardize_lr_pair(p) for p in autocrine_pairs[:n_top]],
                    results=None,  # Raw matrix format, not DataFrame
                )


async def _analyze_communication_liana(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> CCCStorage:
    """Analyze cell communication using LIANA+.

    Returns:
        CCCStorage with unified results structure
    """
    require("liana")
    import liana as li  # noqa: F401

    try:
        # Ensure spatial connectivity is computed for spatial analysis
        if (
            params.perform_spatial_analysis
            and "spatial_connectivities" not in adata.obsp
        ):
            bandwidth = params.liana_bandwidth or (300 if adata.n_obs > 3000 else 200)
            require("squidpy")
            import squidpy as sq

            sq.gr.spatial_neighbors(
                adata,
                coord_type="generic",
                n_neighs=min(30, max(6, adata.n_obs // 100)),
                radius=bandwidth,
                delaunay=True,
                set_diag=False,
            )

        if not params.species:
            raise ParameterError(
                "Species parameter is required! "
                "Specify species='human', 'mouse', or 'zebrafish'."
            )

        # Determine analysis type
        has_clusters = params.cell_type_key in adata.obs.columns

        if has_clusters and not params.perform_spatial_analysis:
            return _run_liana_cluster_analysis(adata, params, ctx)
        else:
            return _run_liana_spatial_analysis(adata, params, ctx)

    except Exception as e:
        raise ProcessingError(f"LIANA+ analysis failed: {e}") from e


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
) -> CCCStorage:
    """Run LIANA+ cluster-based analysis.

    Returns:
        CCCStorage with unified results structure
    """
    import liana as li

    groupby_col = params.cell_type_key
    validate_obs_column(adata, groupby_col, "Cell type")

    resource_name = _get_liana_resource_name(params.species, params.liana_resource)
    n_perms = params.liana_n_perms

    # Check for raw data
    raw_result = get_raw_data_source(adata, prefer_complete_genes=True)
    use_raw = raw_result.source == "raw"

    # Run LIANA+ rank aggregate
    li.mt.rank_aggregate(
        adata,
        groupby=groupby_col,
        resource_name=resource_name,
        expr_prop=params.liana_nz_prop,
        min_cells=params.min_cells,
        n_perms=n_perms,
        verbose=False,
        use_raw=use_raw,
    )

    # Get results from LIANA's storage location
    liana_res = adata.uns["liana_res"]
    n_lr_pairs = len(liana_res)

    # Calculate significance
    significance_alpha = params.liana_significance_alpha
    n_significant = len(liana_res[liana_res["magnitude_rank"] <= significance_alpha])

    # Extract top LR pairs (standardized format: LIGAND_RECEPTOR)
    lr_pairs: list[str] = []
    top_lr_pairs: list[str] = []

    if "magnitude_rank" in liana_res.columns and len(liana_res) > 0:
        # All LR pairs
        lr_pairs = [
            f"{row['ligand_complex']}_{row['receptor_complex']}"
            for _, row in liana_res.iterrows()
        ]
        # Top pairs by magnitude rank
        top_df = liana_res.nsmallest(params.plot_top_pairs, "magnitude_rank")
        top_lr_pairs = [
            f"{row['ligand_complex']}_{row['receptor_complex']}"
            for _, row in top_df.iterrows()
        ]

    # Build unified storage structure
    return CCCStorage(
        method="liana",
        analysis_type="cluster",
        species=params.species,
        database=resource_name,
        lr_pairs=lr_pairs,
        top_lr_pairs=top_lr_pairs,
        n_pairs=n_lr_pairs,
        n_significant=n_significant,
        results=liana_res,
        pvalues=None,  # P-values embedded in results DataFrame (magnitude_rank column)
        statistics={
            "groupby": groupby_col,
            "n_permutations": n_perms,
            "significance_threshold": significance_alpha,
            "use_raw": use_raw,
        },
        method_data={},  # All data in results field
    )


def _run_liana_spatial_analysis(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> CCCStorage:
    """Run LIANA+ spatial bivariate analysis.

    Returns:
        CCCStorage with unified results structure
    """
    import liana as li
    from statsmodels.stats.multitest import multipletests

    resource_name = _get_liana_resource_name(params.species, params.liana_resource)
    n_perms = params.liana_n_perms
    nz_prop = params.liana_nz_prop
    global_metric = params.liana_global_metric
    alpha = params.liana_significance_alpha

    # Run LIANA+ bivariate analysis
    lrdata = li.mt.bivariate(
        adata,
        resource_name=resource_name,
        local_name=params.liana_local_metric,
        global_name=global_metric,
        n_perms=n_perms,
        mask_negatives=False,
        add_categories=True,
        nz_prop=nz_prop,
        use_raw=False,
        verbose=False,
    )

    n_lr_pairs = lrdata.n_vars

    # FDR correction for significance
    pvals_col = f"{global_metric}_pvals"
    pvals = lrdata.var[pvals_col]
    reject, pvals_corrected, _, _ = multipletests(pvals, alpha=alpha, method="fdr_bh")
    n_significant = int(reject.sum())

    # Store corrected p-values in results
    lrdata.var[f"{pvals_col}_corrected"] = pvals_corrected
    lrdata.var[f"{global_metric}_significant"] = reject

    # Extract LR pairs (standardized format)
    all_lr_pairs = [standardize_lr_pair(p) for p in lrdata.var.index.tolist()]

    # Top pairs by global metric
    top_df = lrdata.var.nlargest(params.plot_top_pairs, global_metric)
    top_lr_pairs = [standardize_lr_pair(p) for p in top_df.index.tolist()]

    # Prepare spatial scores for storage
    spatial_scores = to_dense(lrdata.X)
    spatial_pvals = (
        to_dense(lrdata.layers["pvals"]) if "pvals" in lrdata.layers else None
    )

    # Build unified storage structure
    return CCCStorage(
        method="liana",
        analysis_type="spatial",
        species=params.species,
        database=resource_name,
        lr_pairs=all_lr_pairs,
        top_lr_pairs=top_lr_pairs,
        n_pairs=n_lr_pairs,
        n_significant=n_significant,
        results=lrdata.var,  # DataFrame with global metrics and p-values
        pvalues=None,  # P-values are in results DataFrame
        statistics={
            "local_metric": params.liana_local_metric,
            "global_metric": global_metric,
            "n_permutations": n_perms,
            "nz_proportion": nz_prop,
            "fdr_method": "Benjamini-Hochberg",
            "alpha": alpha,
        },
        method_data={
            "spatial_scores": spatial_scores,  # Per-spot scores for adata.obsm
            "spatial_pvals": spatial_pvals,  # Per-spot p-values for adata.obsm
        },
    )


def _ensure_cellphonedb_database(output_dir: str, ctx: "ToolContext") -> str:
    """Ensure CellPhoneDB database is available, download if not exists"""
    # Use centralized dependency manager for consistent error handling
    require("cellphonedb")  # Raises ImportError with install instructions if missing
    import os
    import ssl

    import certifi
    from cellphonedb.utils import db_utils

    # Check if database file already exists
    db_path = os.path.join(output_dir, "cellphonedb.zip")

    if os.path.exists(db_path):
        return db_path

    # Fix macOS SSL certificate issue: patch urllib to use certifi certificates
    # CellPhoneDB uses urllib.request.urlopen which fails on macOS without this fix
    original_https_context = ssl._create_default_https_context
    # Monkeypatch to use certifi certificates - type mismatch is expected
    ssl._create_default_https_context = lambda: ssl.create_default_context(  # type: ignore[misc,assignment]
        cafile=certifi.where()
    )

    try:
        # Download latest database
        db_utils.download_database(output_dir, "v5.0.0")

        return db_path

    except Exception as e:
        error_msg = (
            f"Failed to download CellPhoneDB database: {e}\n\n"
            "Troubleshooting:\n"
            "1. Check internet connection\n"
            "2. Verify CellPhoneDB version compatibility\n"
            "3. Try manually downloading database:\n"
            "   from cellphonedb.utils import db_utils\n"
            "   db_utils.download_database('/path/to/dir', 'v5.0.0')"
        )
        raise DependencyError(error_msg) from e

    finally:
        # Restore original SSL context
        ssl._create_default_https_context = original_https_context


async def _analyze_communication_cellphonedb(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> CCCStorage:
    """Analyze cell communication using CellPhoneDB.

    Returns:
        CCCStorage with unified results structure
    """
    # Use centralized dependency manager for consistent error handling
    require("cellphonedb")  # Raises ImportError with install instructions if missing
    import os
    import tempfile

    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

    try:
        import time

        start_time = time.time()

        # Initialize for finally block cleanup
        microenvs_file = None

        # Species check: CellPhoneDB is human-only
        # Reference: https://github.com/ventolab/cellphonedb
        # "CellphoneDB is a publicly available repository of HUMAN curated
        # receptors, ligands and their interactions"
        if params.species != "human":
            raise ParameterError(
                f"CellPhoneDB only supports human data. "
                f"Your data species: '{params.species}'. "
                f"For {params.species} data, please use:\n"
                f"  - method='liana' with liana_resource='mouseconsensus' (for mouse)\n"
                f"  - method='cellchat_r' (has built-in mouse/human databases)"
            )

        # Use cell_type_key from params (required field, no auto-detect)
        cell_type_col = params.cell_type_key

        validate_obs_column(adata, cell_type_col, "Cell type")

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
                db_path = _ensure_cellphonedb_database(temp_dir, ctx)
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
            # Results will be stored in unified CCCStorage.method_data

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

        # Update significant_means to match filtered results
        if n_significant_pairs > 0:
            significant_indices = means.index[mask]
            significant_means = means.loc[significant_indices]
        else:
            # No significant interactions found
            await ctx.warning(
                f"No significant interactions found at p < {threshold}. Consider adjusting threshold or using method='liana'."
            )

        # Get top LR pairs
        # CellPhoneDB returns interactions in 'interacting_pair' column
        top_lr_pairs = []
        if (
            significant_means is not None
            and hasattr(significant_means, "head")
            and hasattr(significant_means, "columns")
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
        correction_stats: dict[str, int | float] = {}
        if correction_method != "none":
            correction_stats["n_uncorrected_significant"] = int(n_uncorrected_sig)
            correction_stats["n_corrected_significant"] = int(n_corrected_sig)
            if n_uncorrected_sig > 0:
                correction_stats["reduction_percentage"] = round(
                    (1 - n_corrected_sig / n_uncorrected_sig) * 100, 2
                )

        # Extract all LR pairs (standardized format)
        all_lr_pairs: list[str] = []
        if means is not None and "interacting_pair" in means.columns:
            all_lr_pairs = [
                standardize_lr_pair(p) for p in means["interacting_pair"].tolist()
            ]
        elif means is not None:
            all_lr_pairs = [standardize_lr_pair(str(p)) for p in means.index.tolist()]

        # Standardize top LR pairs
        top_lr_standardized = [standardize_lr_pair(p) for p in top_lr_pairs]

        statistics = {
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

        # Build unified storage structure
        return CCCStorage(
            method="cellphonedb",
            analysis_type="cluster",
            species=params.species,
            database="cellphonedb",
            lr_pairs=all_lr_pairs,
            top_lr_pairs=top_lr_standardized,
            n_pairs=n_lr_pairs,
            n_significant=n_significant_pairs,
            results=means,  # Main results DataFrame
            pvalues=pvalues,  # P-values DataFrame
            statistics=statistics,
            method_data={
                "deconvoluted": deconvoluted,  # CellPhoneDB deconvoluted results
                "significant_means": significant_means,  # Filtered significant means
                "min_pvals_corrected": min_pvals_corrected,  # Corrected p-values
            },
        )

    except Exception as e:
        raise ProcessingError(f"CellPhoneDB analysis failed: {e}") from e
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
        validate_obs_column(adata, params.cell_type_key, "Cell type")

        cell_types = adata.obs[params.cell_type_key].values

        # Create microenvironments by cell type co-occurrence
        # Optimized: Vectorized neighbor extraction (6.7x faster than row-by-row)
        rows, cols = neighbor_matrix.nonzero()
        neighbor_types = cell_types[cols]

        # Aggregate cell types per row using defaultdict
        row_to_types = defaultdict(set)
        for r, ct in zip(rows, neighbor_types):
            row_to_types[r].add(ct)

        # Build microenvironment assignments
        microenv_assignments = {}
        cell_type_to_microenv = defaultdict(set)
        microenv_counter = 0

        for neighbor_cell_types in row_to_types.values():
            if len(neighbor_cell_types) > 1:  # At least 2 different cell types
                microenv_signature = tuple(sorted(neighbor_cell_types))

                if microenv_signature not in microenv_assignments:
                    microenv_assignments[microenv_signature] = (
                        f"microenv_{microenv_counter}"
                    )
                    microenv_counter += 1

                microenv_name = microenv_assignments[microenv_signature]
                for ct in neighbor_cell_types:
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
        await ctx.warning(f"Failed to create microenvironments file: {e}")
        return None


def _analyze_communication_cellchat_r(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> CCCStorage:
    """Analyze cell communication using native R CellChat package.

    This implementation uses rpy2 to call the original R CellChat package,
    which includes full features like mediator proteins and signaling pathways
    that are not available in the LIANA simplified implementation.

    Returns:
        CCCStorage with unified results structure
    """
    import pandas as pd
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter

    try:
        import time

        start_time = time.time()

        # Validate cell type column
        validate_obs_column(adata, params.cell_type_key, "Cell type")

        # Check for spatial data
        spatial_key = get_spatial_key(adata)
        has_spatial = spatial_key is not None

        # Prepare expression matrix (genes x cells, normalized)
        # CellChat requires normalized data with comprehensive gene coverage
        # Use get_raw_data_source (single source of truth) - directly use raw_result
        raw_result = get_raw_data_source(adata, prefer_complete_genes=True)

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
            common_genes = raw_result.var_names.intersection(cellchat_genes)

            if len(common_genes) == 0:
                raise DataCompatibilityError(
                    f"No genes overlap between data and {db_name}. "
                    f"Check if species parameter matches your data."
                )

            # Create expression matrix with only CellChatDB genes (memory efficient)
            gene_indices = [raw_result.var_names.get_loc(g) for g in common_genes]
            expr_matrix = pd.DataFrame(
                to_dense(raw_result.X[:, gene_indices]).T,
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

            # Convert results back to Python (force native types for h5ad compatibility)
            n_lr_pairs = int(ro.r("n_lr_pairs")[0])
            n_significant_pairs = int(ro.r("n_significant")[0])
            # Force str() to ensure Python native strings (not rpy2 objects)
            top_pathways = [str(x) for x in ro.r("top_pathways")]
            top_lr_pairs = [str(x) for x in ro.r("top_lr")]

            # Get full results for storage
            lr_pairs_df = ro.r("as.data.frame(lr_pairs)")
            # net$prob is a 3D array (n_cell_types x n_cell_types x n_lr_pairs)
            # Don't use as.matrix() which would flatten it
            prob_matrix = ro.r("net$prob")
            pval_matrix = ro.r("net$pval")
            # Get cell type names for later use
            cell_type_names = [str(x) for x in ro.r("rownames(net$prob)")]

            # Convert results to Python objects with explicit type conversion
            lr_pairs_df_py = pd.DataFrame(lr_pairs_df)
            # Ensure all string columns are Python str (not rpy2 objects)
            for col in lr_pairs_df_py.select_dtypes(include=["object"]).columns:
                lr_pairs_df_py[col] = lr_pairs_df_py[col].astype(str)
            prob_matrix_np = np.array(prob_matrix, dtype=np.float64)
            pval_matrix_np = np.array(pval_matrix, dtype=np.float64)

        end_time = time.time()
        analysis_time = end_time - start_time

        # Extract all LR pairs (standardized format)
        all_lr_pairs: list[str] = []
        if "interaction_name" in lr_pairs_df_py.columns:
            all_lr_pairs = [
                standardize_lr_pair(p)
                for p in lr_pairs_df_py["interaction_name"].tolist()
            ]

        # Standardize top LR pairs
        top_lr_standardized = [standardize_lr_pair(p) for p in top_lr_pairs]

        # Build unified storage structure
        return CCCStorage(
            method="cellchat_r",
            analysis_type="cluster",
            species=params.species,
            database=f"CellChatDB.{params.species}",
            lr_pairs=all_lr_pairs,
            top_lr_pairs=top_lr_standardized,
            n_pairs=n_lr_pairs,
            n_significant=n_significant_pairs,
            results=lr_pairs_df_py,  # LR pairs DataFrame
            pvalues=None,  # P-values are in 3D matrix format
            statistics={
                "db_category": params.cellchat_db_category,
                "aggregation_type": params.cellchat_type,
                "trim": params.cellchat_trim,
                "population_size": params.cellchat_population_size,
                "min_cells": params.cellchat_min_cells,
                "spatial_mode": has_spatial and params.cellchat_distance_use,
                "analysis_time_seconds": analysis_time,
                "top_pathways": top_pathways[:5] if top_pathways else [],
            },
            method_data={
                "prob_matrix": prob_matrix_np,  # 3D array (sources x targets x interactions)
                "pval_matrix": pval_matrix_np,  # 3D array (sources x targets x interactions)
                "cell_type_names": cell_type_names,  # Cell type names for matrix axes
                "top_pathways": top_pathways,  # Ranked pathway names
            },
        )

    except Exception as e:
        raise ProcessingError(f"CellChat R analysis failed: {e}") from e


async def _analyze_communication_fastccc(
    adata: Any, params: CellCommunicationParameters, ctx: "ToolContext"
) -> CCCStorage:
    """Analyze cell communication using FastCCC permutation-free framework.

    FastCCC uses FFT-based convolution to compute p-values analytically,
    making it extremely fast for large datasets (16M cells in minutes).

    Reference: Nature Communications 2025 (https://github.com/Svvord/FastCCC)

    Args:
        adata: AnnData object with expression data
        params: Cell communication analysis parameters
        ctx: ToolContext for logging and data access

    Returns:
        CCCStorage with unified results
    """
    import glob
    import os
    import tempfile
    import time

    import pandas as pd

    from ..utils.adata_utils import to_dense

    try:
        start_time = time.time()

        # Species check: FastCCC uses CellPhoneDB v5 which is human-only
        # Reference: https://github.com/ventolab/cellphonedb
        # "CellphoneDB is a publicly available repository of HUMAN curated
        # receptors, ligands and their interactions"
        if params.species != "human":
            raise ParameterError(
                f"FastCCC only supports human data (uses CellPhoneDB v5 database). "
                f"Your data species: '{params.species}'. "
                f"For {params.species} data, please use:\n"
                f"  - method='liana' with liana_resource='mouseconsensus' (for mouse)\n"
                f"  - method='cellchat_r' (has built-in mouse/human databases)"
            )

        # Import FastCCC
        if params.fastccc_use_cauchy:
            from fastccc import Cauchy_combination_of_statistical_analysis_methods
        else:
            from fastccc import statistical_analysis_method

        # Validate cell type column
        validate_obs_column(adata, params.cell_type_key, "Cell type")

        # Use get_raw_data_source (single source of truth) - directly use raw_result
        raw_result = get_raw_data_source(adata, prefer_complete_genes=True)

        # Create temporary directory for FastCCC I/O
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save expression data as h5ad for FastCCC
            counts_file = os.path.join(temp_dir, "counts.h5ad")

            # Create a minimal AnnData for saving (FastCCC reads h5ad directly)
            # IMPORTANT: FastCCC requires normalized log1p-transformed data
            # with max value < 14 (default threshold)
            import anndata as ad
            import scanpy as sc

            # Prepare expression matrix (cells × genes)
            # Use copy=True to ensure safe modification for normalize_total/log1p
            expr_matrix = to_dense(raw_result.X, copy=True)
            gene_names = list(raw_result.var_names)
            cell_names = list(adata.obs_names)

            # Create temporary AnnData
            temp_adata = ad.AnnData(
                X=expr_matrix,
                obs=pd.DataFrame(index=cell_names),
                var=pd.DataFrame(index=gene_names),
            )

            # Check if data needs normalization (FastCCC max threshold is 14)
            max_val = np.max(temp_adata.X)
            if max_val > 14:
                # Apply standard scanpy normalization pipeline
                sc.pp.normalize_total(temp_adata, target_sum=1e4)
                sc.pp.log1p(temp_adata)

            # Make var names unique (FastCCC requirement)
            temp_adata.var_names_make_unique()

            # Add cell type labels to obs
            temp_adata.obs[params.cell_type_key] = adata.obs[
                params.cell_type_key
            ].values

            # Save to h5ad
            temp_adata.write_h5ad(counts_file)

            # Get database directory path (FastCCC uses CellPhoneDB database format)
            # FastCCC expects a directory containing interaction_table.csv and other files
            # Check for bundled database in chatspatial package
            chatspatial_pkg_dir = os.path.dirname(os.path.dirname(__file__))
            database_dir = os.path.join(
                chatspatial_pkg_dir,
                "data",
                "cellphonedb_v5",
                "cellphonedb-data-5.0.0",
            )

            # Verify required files exist
            required_file = os.path.join(database_dir, "interaction_table.csv")
            if not os.path.exists(required_file):
                raise ProcessingError(
                    f"FastCCC requires CellPhoneDB database files. "
                    f"Expected directory: {database_dir} with interaction_table.csv. "
                    f"Please download from: https://github.com/ventolab/cellphonedb-data"
                )

            # Output directory for results
            output_dir = os.path.join(temp_dir, "results")
            os.makedirs(output_dir, exist_ok=True)

            # Run FastCCC analysis
            if params.fastccc_use_cauchy:
                # Cauchy combination method (more robust, multiple parameter combinations)
                # Note: This function saves results to files and returns None
                Cauchy_combination_of_statistical_analysis_methods(
                    database_file_path=database_dir,
                    celltype_file_path=None,  # Using meta_key instead
                    counts_file_path=counts_file,
                    convert_type="hgnc_symbol",
                    single_unit_summary_list=[
                        "Mean",
                        "Median",
                        "Q3",
                        "Quantile_0.9",
                    ],
                    complex_aggregation_list=["Minimum", "Average"],
                    LR_combination_list=["Arithmetic", "Geometric"],
                    min_percentile=params.fastccc_min_percentile,
                    save_path=output_dir,
                    meta_key=params.cell_type_key,
                    use_DEG=params.fastccc_use_deg,
                )

                # Read results from saved files (Cauchy method saves to files)
                # Find the task ID from output files
                pval_files = glob.glob(os.path.join(output_dir, "*_Cauchy_pvals.tsv"))
                if not pval_files:
                    raise ProcessingError(
                        "FastCCC Cauchy combination did not produce output files"
                    )

                # Extract task_id from filename
                pval_file = pval_files[0]
                task_id = os.path.basename(pval_file).replace("_Cauchy_pvals.tsv", "")

                # Read combined results
                pvalues = pd.read_csv(pval_file, index_col=0, sep="\t")
                strength_file = os.path.join(
                    output_dir, f"{task_id}_average_interactions_strength.tsv"
                )
                interactions_strength = pd.read_csv(
                    strength_file, index_col=0, sep="\t"
                )

                # Percentages are in individual method files, use first one
                pct_files = glob.glob(
                    os.path.join(output_dir, f"{task_id}*percents_analysis.tsv")
                )
                if pct_files:
                    percentages = pd.read_csv(pct_files[0], index_col=0, sep="\t")
                else:
                    percentages = None

            else:
                # Single method (faster)
                interactions_strength, pvalues, percentages = (
                    statistical_analysis_method(
                        database_file_path=database_dir,
                        celltype_file_path=None,  # Using meta_key instead
                        counts_file_path=counts_file,
                        convert_type="hgnc_symbol",
                        single_unit_summary=params.fastccc_single_unit_summary,
                        complex_aggregation=params.fastccc_complex_aggregation,
                        LR_combination=params.fastccc_lr_combination,
                        min_percentile=params.fastccc_min_percentile,
                        save_path=output_dir,
                        meta_key=params.cell_type_key,
                        use_DEG=params.fastccc_use_deg,
                    )
                )

        # Process results
        n_lr_pairs = len(pvalues) if pvalues is not None else 0

        # Count significant pairs
        threshold = params.fastccc_pvalue_threshold
        if pvalues is not None and hasattr(pvalues, "values"):
            # Get minimum p-value across all cell type pairs for each LR pair
            pval_array = pvalues.select_dtypes(include=[np.number]).values
            min_pvals = np.nanmin(pval_array, axis=1)
            n_significant_pairs = int(np.sum(min_pvals < threshold))
        else:
            n_significant_pairs = 0

        # Get all LR pairs and top LR pairs based on interaction strength
        all_lr_pairs = []
        top_lr_pairs_raw = []
        if interactions_strength is not None and hasattr(
            interactions_strength, "index"
        ):
            # All LR pairs (standardized format)
            all_lr_pairs = [
                standardize_lr_pair(str(p)) for p in interactions_strength.index
            ]

            # Sort by mean interaction strength across cell type pairs
            if hasattr(interactions_strength, "select_dtypes"):
                strength_array = interactions_strength.select_dtypes(
                    include=[np.number]
                ).values
                mean_strength = np.nanmean(strength_array, axis=1)
                top_indices = np.argsort(mean_strength)[::-1][: params.plot_top_pairs]
                top_lr_pairs_raw = [interactions_strength.index[i] for i in top_indices]

        # Standardize top LR pairs
        top_lr_standardized = [standardize_lr_pair(str(p)) for p in top_lr_pairs_raw]

        end_time = time.time()
        analysis_time = end_time - start_time

        # Build CCCStorage
        return CCCStorage(
            method="fastccc",
            analysis_type="cluster",  # FastCCC is cluster-level analysis
            species=params.species,
            database="CellPhoneDB_v5",
            lr_pairs=all_lr_pairs,
            top_lr_pairs=top_lr_standardized,
            n_pairs=n_lr_pairs,
            n_significant=n_significant_pairs,
            results=interactions_strength,
            pvalues=pvalues,
            statistics={
                "use_cauchy": params.fastccc_use_cauchy,
                "single_unit_summary": params.fastccc_single_unit_summary,
                "complex_aggregation": params.fastccc_complex_aggregation,
                "lr_combination": params.fastccc_lr_combination,
                "min_percentile": params.fastccc_min_percentile,
                "pvalue_threshold": threshold,
                "use_deg": params.fastccc_use_deg,
                "n_lr_pairs_tested": n_lr_pairs,
                "analysis_time_seconds": analysis_time,
                "permutation_free": True,  # Key FastCCC feature
            },
            method_data={
                "percentages": percentages,  # Expression percentages per cell type
            },
        )

    except Exception as e:
        raise ProcessingError(f"FastCCC analysis failed: {e}") from e


# Legacy autocrine functions removed - now using _integrate_autocrine_detection()
# which extracts directly from CCCStorage for consistency
