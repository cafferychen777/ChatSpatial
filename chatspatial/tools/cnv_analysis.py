"""
Copy Number Variation (CNV) analysis tools for spatial transcriptomics data.
"""

import glob
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import scanpy as sc

from ..models.analysis import CNVResult
from ..models.data import CNVParameters
from ..utils import validate_obs_column
from ..utils.adata_utils import store_analysis_metadata
from ..utils.dependency_manager import require, validate_r_environment
from ..utils.exceptions import (
    ChatSpatialError,
    DataCompatibilityError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)
from ..utils.results_export import export_analysis_result

if TYPE_CHECKING:
    import anndata as ad

    from ..spatial_mcp_adapter import ToolContext

logger = logging.getLogger(__name__)

_NUMBAT_ALLELE_COLUMNS = (
    "cell",
    "snp_id",
    "CHROM",
    "POS",
    "cM",
    "REF",
    "ALT",
    "AD",
    "DP",
    "GT",
)
_NUMBAT_CLONE_COLUMNS = ("cell", "clone_opt", "p_cnv", "compartment_opt")


def _resolve_numbat_table_path(out_dir: str, stem: str, iteration: str) -> Path | None:
    """Return a Numbat TSV output path across documented compression variants."""
    base_path = Path(out_dir) / f"{stem}_{iteration}.tsv"
    for candidate in (base_path, Path(f"{base_path}.gz")):
        if candidate.is_file():
            return candidate
    return None


# Per-spot CNV burden. The mean absolute CNV was already computed for the
# summary statistics and then discarded, while the spatial visualization —
# the reason to run CNV on spatial data at all — looks for it here and told
# the user to "Run analyze_cnv() first" when they just had.
CNV_SCORE_OBS_KEY = "cnv_score"


def _build_cnv_key(params: "CNVParameters") -> str:
    """Build a parametric analysis key for CNV results.

    Encodes method + reference categories so that runs with different
    reference baselines coexist in metadata/export/cache.
    """
    method = params.method
    if params.reference_categories:
        # Sort for deterministic key regardless of input order
        cats = "_".join(sorted(params.reference_categories))
        return f"cnv_{method}_{cats}"
    return f"cnv_{method}"


# Numbat availability is checked lazily in _infer_cnv_numbat to avoid
# import-time failures when rpy2/R is not installed


def _copy_matrix_data(data):
    """Return an independent copy of dense/sparse matrix-like data."""
    if hasattr(data, "copy"):
        return data.copy()
    return np.array(data, copy=True)


def _validate_gene_positions(adata: "ad.AnnData") -> None:
    required_position_columns = {"chromosome", "start", "end"}
    missing_position_columns = sorted(
        required_position_columns - set(adata.var.columns)
    )
    if missing_position_columns:
        raise DataCompatibilityError(
            "CNV inference requires genomic positions in adata.var. Missing "
            f"columns: {missing_position_columns} (expected chromosome, start, "
            "end). ChatSpatial does not annotate gene positions; add them from "
            "a GTF for your genome build, e.g. "
            "infercnvpy.io.genomic_position_from_gtf(gtf_path, adata), then "
            "rerun."
        )
    if not adata.var_names.is_unique:
        raise DataCompatibilityError(
            "CNV inference requires unique gene names for position-aligned results."
        )


def _expand_gene_aligned_layer(
    layer: Any,
    source_genes: pd.Index,
    target_genes: pd.Index,
    n_obs: int,
) -> Any:
    """Expand a gene-aligned layer to a target gene index without reordering values."""
    import scipy.sparse as sp_sparse

    if not source_genes.is_unique or not target_genes.is_unique:
        raise DataCompatibilityError(
            "Gene-aligned CNV results require unique source and target gene names."
        )
    if len(layer.shape) != 2 or layer.shape != (n_obs, len(source_genes)):
        raise DataCompatibilityError(
            "Gene-level CNV matrix does not match its source axes: "
            f"received {layer.shape}, expected {(n_obs, len(source_genes))}."
        )

    target_positions = target_genes.get_indexer(source_genes)
    if (target_positions < 0).any():
        missing = source_genes[target_positions < 0].tolist()
        raise DataCompatibilityError(
            f"CNV result genes are missing from the target dataset: {missing[:10]}."
        )

    if sp_sparse.issparse(layer):
        expanded = sp_sparse.lil_matrix((n_obs, len(target_genes)), dtype=layer.dtype)
        expanded[:, target_positions] = layer
        return expanded.tocsr()

    expanded = np.zeros((n_obs, len(target_genes)), dtype=layer.dtype)
    expanded[:, target_positions] = layer
    return expanded


def _validate_numbat_allele_data(data: object) -> pd.DataFrame:
    """Validate the allele-count fields consumed by Numbat before crossing into R."""
    if not isinstance(data, pd.DataFrame):
        raise ParameterError(
            "adata.uns['numbat_allele_data_raw'] must be a pandas DataFrame."
        )
    if data.empty:
        raise ParameterError("Numbat allele dataframe must contain at least one row.")

    missing_columns = [
        column for column in _NUMBAT_ALLELE_COLUMNS if column not in data.columns
    ]
    if missing_columns:
        raise ParameterError(
            f"Allele dataframe missing required columns: {missing_columns}\n"
            f"Available columns: {list(data.columns)}\n"
            "Numbat requires phased allele counts with columns: "
            f"{', '.join(_NUMBAT_ALLELE_COLUMNS)}."
        )
    return data


def _normalize_numbat_cell_table(
    table: pd.DataFrame,
    *,
    required_columns: tuple[str, ...],
    expected_cells: pd.Index,
    table_name: str,
) -> pd.DataFrame:
    """Validate a Numbat cell table and normalize its cell identifiers."""
    if not isinstance(table, pd.DataFrame):
        raise DataCompatibilityError(
            f"Numbat {table_name} output must be a pandas DataFrame."
        )
    if table.empty:
        raise DataCompatibilityError(f"Numbat {table_name} output is empty.")
    if not table.columns.is_unique:
        raise DataCompatibilityError(
            f"Numbat {table_name} output contains duplicate column labels."
        )

    missing_columns = [
        column for column in required_columns if column not in table.columns
    ]
    if missing_columns:
        raise DataCompatibilityError(
            f"Numbat {table_name} output is missing required columns: "
            f"{missing_columns}."
        )

    normalized = table.copy()
    if normalized["cell"].isna().any():
        raise DataCompatibilityError(
            f"Numbat {table_name} output contains missing cell identifiers."
        )
    normalized["cell"] = normalized["cell"].map(str)
    if normalized["cell"].duplicated().any():
        raise DataCompatibilityError(
            f"Numbat {table_name} output contains duplicate cell identifiers."
        )

    unexpected_cells = pd.Index(normalized["cell"]).difference(expected_cells)
    if len(unexpected_cells):
        raise DataCompatibilityError(
            f"Numbat {table_name} output contains cells absent from the dataset: "
            f"{unexpected_cells[:10].tolist()}."
        )
    return normalized


def _validate_and_align_numbat_outputs(
    clone_post: pd.DataFrame,
    geno: pd.DataFrame,
    obs_names: pd.Index,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """Validate Numbat output axes and align cell-level results to AnnData order."""
    expected_cells = pd.Index(obs_names.map(str), name="cell")
    if not expected_cells.is_unique:
        raise DataCompatibilityError(
            "Numbat result alignment requires unique AnnData observation names."
        )

    normalized_clone_post = _normalize_numbat_cell_table(
        clone_post,
        required_columns=_NUMBAT_CLONE_COLUMNS,
        expected_cells=expected_cells,
        table_name="clone_post",
    )
    if normalized_clone_post[["clone_opt", "compartment_opt"]].isna().any().any():
        raise DataCompatibilityError(
            "Numbat clone_post output contains missing clone assignments."
        )
    try:
        clone_probabilities = normalized_clone_post["p_cnv"].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise DataCompatibilityError(
            "Numbat clone_post p_cnv values must be numeric."
        ) from exc
    if (
        not np.isfinite(clone_probabilities).all()
        or (clone_probabilities < 0).any()
        or (clone_probabilities > 1).any()
    ):
        raise DataCompatibilityError(
            "Numbat clone_post p_cnv values must be finite probabilities in [0, 1]."
        )
    normalized_clone_post["p_cnv"] = clone_probabilities

    normalized_geno = _normalize_numbat_cell_table(
        geno,
        required_columns=("cell",),
        expected_cells=expected_cells,
        table_name="geno",
    )
    segment_columns = normalized_geno.columns.drop("cell")
    if len(segment_columns) == 0:
        raise DataCompatibilityError(
            "Numbat geno output contains no CNV segment columns."
        )
    try:
        genotype_values = normalized_geno[segment_columns].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise DataCompatibilityError(
            "Numbat geno segment values must be numeric."
        ) from exc
    if (
        not np.isfinite(genotype_values).all()
        or (genotype_values < 0).any()
        or (genotype_values > 1).any()
    ):
        raise DataCompatibilityError(
            "Numbat geno segment values must be finite probabilities in [0, 1]."
        )

    aligned_genotypes = np.full(
        (len(expected_cells), len(segment_columns)), np.nan, dtype=float
    )
    genotype_positions = expected_cells.get_indexer(normalized_geno["cell"])
    aligned_genotypes[genotype_positions] = genotype_values

    aligned_clone_post = normalized_clone_post.set_index("cell").reindex(expected_cells)
    aligned_clone_post["clone_opt"] = (
        aligned_clone_post["clone_opt"].fillna("unassigned").map(str)
    )
    aligned_clone_post["compartment_opt"] = (
        aligned_clone_post["compartment_opt"].fillna("unassigned").map(str)
    )

    return normalized_clone_post, aligned_clone_post, aligned_genotypes


def _build_infercnvpy_workspace(adata: "ad.AnnData") -> "ad.AnnData":
    """Build minimal AnnData workspace for infercnvpy without copying unrelated fields."""
    import anndata as ad

    adata_cnv = ad.AnnData(
        X=_copy_matrix_data(adata.X),
        obs=adata.obs.copy(),
        var=adata.var.copy(),
    )
    adata_cnv.obs_names = adata.obs_names.copy()
    adata_cnv.var_names = adata.var_names.copy()
    return adata_cnv


async def infer_cnv(
    data_id: str,
    ctx: "ToolContext",
    params: CNVParameters,
) -> CNVResult:
    """Infer copy number variations using selected method

    Supports two methods:
    - infercnvpy: Expression-based CNV inference (default, fast)
    - Numbat: Haplotype-aware CNV analysis (requires allele data, more accurate)

    Args:
        data_id: Dataset identifier
        ctx: Tool context for data access and logging
        params: CNV analysis parameters including method selection

    Returns:
        CNVResult containing method-specific CNV analysis results

    Raises:
        ValueError: If dataset not found or parameters are invalid
        RuntimeError: If selected method is not available
    """
    # Validate against the managed source, then give both backends an isolated
    # publication candidate. Their internal workspaces are method-specific,
    # but both write result keys incrementally before export completes.
    source_adata = await ctx.get_adata(data_id)

    # Validate common parameters
    validate_obs_column(source_adata, params.reference_key, "Reference cell type")

    available_categories = set(source_adata.obs[params.reference_key].unique())
    missing_categories = set(params.reference_categories) - available_categories
    if missing_categories:
        raise ParameterError(
            f"Reference categories {missing_categories} not found in "
            f"adata.obs['{params.reference_key}'].\n"
            f"Available categories: {sorted(available_categories)}"
        )

    if params.method not in ("infercnvpy", "numbat"):
        raise ParameterError(
            f"Unknown CNV method: {params.method}. "
            "Available methods: 'infercnvpy', 'numbat'"
        )

    adata = source_adata.copy()

    # Dispatch to appropriate method
    if params.method == "infercnvpy":
        result = await _infer_cnv_infercnvpy(data_id, adata, params, ctx)
    else:
        result = _infer_cnv_numbat(data_id, adata, params, ctx)

    await ctx.set_adata(data_id, adata)
    return result


async def _infer_cnv_infercnvpy(
    data_id: str,
    adata: "ad.AnnData",
    params: CNVParameters,
    ctx: "ToolContext",
) -> CNVResult:
    """Infer copy number variations using infercnvpy

    This function performs CNV inference on spatial transcriptomics data using
    infercnvpy, which detects chromosomal copy number alterations by comparing
    gene expression patterns across chromosomes between tumor and normal cells.

    Args:
        data_id: Dataset identifier (for result creation)
        adata: AnnData object (already retrieved via ctx.get_adata)
        params: CNV analysis parameters
        ctx: Tool context for logging

    Returns:
        CNVResult containing CNV analysis results and statistics
    """
    _validate_gene_positions(adata)

    cnv = require("infercnvpy", ctx, feature="CNV analysis")

    # Note: adata is already validated in infer_cnv() before dispatch.
    # Build a minimal workspace to avoid copying unrelated layers/obsm/uns.
    adata_cnv = _build_infercnvpy_workspace(adata)

    has_complete_position = (
        adata_cnv.var[["chromosome", "start", "end"]].notna().all(axis=1)
    )
    n_missing_positions = int((~has_complete_position).sum())
    if n_missing_positions:
        await ctx.warning(
            f"Excluding {n_missing_positions} genes with incomplete genomic positions."
        )

    genes_to_keep = has_complete_position
    if params.exclude_chromosomes is not None:
        genes_to_keep &= ~adata_cnv.var["chromosome"].isin(params.exclude_chromosomes)
    if not genes_to_keep.any():
        raise DataCompatibilityError(
            "No genes with usable genomic positions remain after chromosome filtering."
        )
    if not genes_to_keep.all():
        adata_cnv = adata_cnv[:, genes_to_keep].copy()

    cnv.tl.infercnv(
        adata_cnv,
        reference_key=params.reference_key,
        reference_cat=params.reference_categories,
        window_size=params.window_size,
        step=params.step,
        dynamic_threshold=params.dynamic_threshold,
        # Filtering is performed above so the user-facing None value means
        # "exclude nothing" instead of inheriting infercnvpy's chrX/chrY default.
        exclude_chromosomes=None,
    )

    # Optional: Cluster cells by CNV pattern
    if params.cluster_cells:
        try:
            sc.pp.neighbors(adata_cnv, use_rep="X_cnv", n_neighbors=15)
            sc.tl.leiden(adata_cnv, key_added="cnv_clusters")
        except Exception as e:
            await ctx.warning(f"Failed to cluster cells by CNV: {e}")

    # Optional: Compute dendrogram
    if params.dendrogram and params.cluster_cells:
        try:
            sc.tl.dendrogram(adata_cnv, groupby="cnv_clusters")
        except Exception as e:
            await ctx.warning(f"Failed to compute dendrogram: {e}")

    # Extract CNV statistics

    # Check what data is available
    cnv_score_key = None
    if "X_cnv" in adata_cnv.obsm:
        cnv_score_key = "X_cnv"
    elif "cnv" in adata_cnv.layers:
        cnv_score_key = "cnv"

    # Calculate statistics
    statistics = {}
    cnv_matrix = None
    if cnv_score_key == "X_cnv" and "X_cnv" in adata_cnv.obsm:
        cnv_matrix = adata_cnv.obsm["X_cnv"]
    elif cnv_score_key == "cnv" and "cnv" in adata_cnv.layers:
        cnv_matrix = adata_cnv.layers["cnv"]
    if cnv_matrix is not None:
        # ==================== OPTIMIZED: Compute statistics on sparse matrix ====================
        # Strategy: infercnvpy outputs sparse CSR matrix after noise filtering (Line 448-452)
        #           Noise filtering sets ~87% values to zero, making sparse computation efficient
        # Benefit: For 5k cells × 500 windows: save ~19 MB (50%), 1.6x faster
        # Technical: All statistics (mean, std, median, per-cell scores) can be computed
        #           directly on sparse matrices without conversion to dense

        import scipy.sparse

        if scipy.sparse.issparse(cnv_matrix):
            # Sparse matrix - compute statistics without toarray()

            # Mean: use sparse matrix's mean() method
            statistics["mean_cnv"] = float(cnv_matrix.mean())

            # Std: manual calculation using E[X^2] - E[X]^2
            mean_val = cnv_matrix.mean()
            mean_sq = cnv_matrix.multiply(cnv_matrix).mean()
            variance = float(mean_sq - mean_val**2)
            statistics["std_cnv"] = float(np.sqrt(max(0.0, variance)))

            # Median: for highly sparse matrices (>50% zeros), median is 0
            # Otherwise use approximation with non-zero values
            n_zeros = cnv_matrix.shape[0] * cnv_matrix.shape[1] - cnv_matrix.nnz
            n_total = cnv_matrix.shape[0] * cnv_matrix.shape[1]

            if n_zeros > n_total / 2:
                # Majority zeros, median is exactly 0
                statistics["median_cnv"] = 0.0
            elif n_total <= 10_000_000:
                # Small enough to compute exact median via dense conversion
                statistics["median_cnv"] = float(np.median(cnv_matrix.toarray()))
            else:
                # Large matrix: merge sorted non-zero values with
                # zero count to find exact median without densifying
                nz_sorted = np.sort(cnv_matrix.data)
                zero_pos = int(np.searchsorted(nz_sorted, 0.0))

                def _val_at(idx: int) -> float:
                    if idx < zero_pos:
                        return float(nz_sorted[idx])
                    elif idx < zero_pos + n_zeros:
                        return 0.0
                    else:
                        return float(nz_sorted[idx - n_zeros])

                if n_total % 2 == 1:
                    statistics["median_cnv"] = _val_at(n_total // 2)
                else:
                    statistics["median_cnv"] = (
                        _val_at(n_total // 2 - 1) + _val_at(n_total // 2)
                    ) / 2

            # Per-cell CNV scores: compute on sparse matrix
            # abs() preserves sparsity
            cnv_abs = cnv_matrix.copy()
            cnv_abs.data = np.abs(cnv_abs.data)
            cell_cnv_scores = np.array(cnv_abs.mean(axis=1)).flatten()
            statistics["mean_cell_cnv_score"] = float(np.mean(cell_cnv_scores))
            statistics["max_cell_cnv_score"] = float(np.max(cell_cnv_scores))
            adata_cnv.obs[CNV_SCORE_OBS_KEY] = cell_cnv_scores

        else:
            # Dense matrix - use standard numpy operations
            statistics["mean_cnv"] = float(np.mean(cnv_matrix))
            statistics["std_cnv"] = float(np.std(cnv_matrix))
            statistics["median_cnv"] = float(np.median(cnv_matrix))

            # Calculate per-cell CNV scores
            cell_cnv_scores = np.mean(np.abs(cnv_matrix), axis=1)
            statistics["mean_cell_cnv_score"] = float(np.mean(cell_cnv_scores))
            statistics["max_cell_cnv_score"] = float(np.max(cell_cnv_scores))
            adata_cnv.obs[CNV_SCORE_OBS_KEY] = np.asarray(cell_cnv_scores).ravel()

    # Count reference vs non-reference cells
    is_reference = adata_cnv.obs[params.reference_key].isin(params.reference_categories)
    statistics["n_reference_cells"] = int(is_reference.sum())
    statistics["n_non_reference_cells"] = int((~is_reference).sum())

    # Get chromosome information
    if "chromosome" in adata_cnv.var.columns:
        n_chromosomes = len(adata_cnv.var["chromosome"].unique())
    else:
        n_chromosomes = 0  # Unknown

    n_genes_analyzed = adata_cnv.n_vars

    # Transfer CNV results from the method workspace to the publication candidate.
    if cnv_score_key == "X_cnv" and "X_cnv" in adata_cnv.obsm:
        adata.obsm["X_cnv"] = adata_cnv.obsm["X_cnv"]
    elif cnv_score_key == "cnv" and "cnv" in adata_cnv.layers:
        cnv_layer = adata_cnv.layers["cnv"]
        if adata_cnv.var_names.equals(adata.var_names):
            # No gene filtering — direct copy
            adata.layers["cnv"] = cnv_layer
        else:
            # Genes were filtered (e.g., exclude_chromosomes).
            # Pad to original shape; excluded genes get CNV score of 0.
            adata.layers["cnv"] = _expand_gene_aligned_layer(
                cnv_layer,
                adata_cnv.var_names,
                adata.var_names,
                adata.n_obs,
            )

    # Store CNV metadata (required for infercnvpy plotting functions)
    if "cnv" in adata_cnv.uns:
        adata.uns["cnv"] = adata_cnv.uns["cnv"]

    if CNV_SCORE_OBS_KEY in adata_cnv.obs and adata_cnv.obs_names.equals(
        adata.obs_names
    ):
        adata.obs[CNV_SCORE_OBS_KEY] = adata_cnv.obs[CNV_SCORE_OBS_KEY]

    if params.cluster_cells and "cnv_clusters" in adata_cnv.obs:
        adata.obs["cnv_clusters"] = adata_cnv.obs["cnv_clusters"]

    if params.dendrogram and "dendrogram_cnv_clusters" in adata_cnv.uns:
        adata.uns["dendrogram_cnv_clusters"] = adata_cnv.uns["dendrogram_cnv_clusters"]

    # Store CNV analysis parameters in adata.uns for reference
    analysis_key = _build_cnv_key(params)
    cnv_summary_key = f"cnv_analysis_{analysis_key.removeprefix('cnv_')}"
    adata.uns[cnv_summary_key] = {
        "reference_key": params.reference_key,
        "reference_categories": list(params.reference_categories),  # Convert to list
        "window_size": params.window_size,
        "step": params.step,
        "cnv_score_key": cnv_score_key,
    }

    # Build results keys for metadata
    results_keys: dict = {"uns": ["cnv", cnv_summary_key]}
    if cnv_score_key == "X_cnv":
        results_keys["obsm"] = ["X_cnv"]
    elif cnv_score_key == "cnv":
        results_keys["layers"] = ["cnv"]
    if params.cluster_cells and "cnv_clusters" in adata.obs:
        results_keys.setdefault("obs", []).append("cnv_clusters")
    if params.dendrogram and "dendrogram_cnv_clusters" in adata.uns:
        results_keys["uns"].append("dendrogram_cnv_clusters")

    # Store metadata for scientific provenance tracking
    store_analysis_metadata(
        adata,
        analysis_name=analysis_key,
        method="infercnvpy",
        parameters={
            "reference_key": params.reference_key,
            "reference_categories": list(params.reference_categories),
            "window_size": params.window_size,
            "step": params.step,
        },
        results_keys=results_keys,
        statistics=statistics,
    )

    # Export results for reproducibility
    export_analysis_result(adata, data_id, analysis_key)

    return CNVResult(
        data_id=data_id,
        method="infercnvpy",
        reference_key=params.reference_key,
        reference_categories=list(params.reference_categories),  # Convert to list
        n_chromosomes=n_chromosomes,
        n_genes_analyzed=n_genes_analyzed,
        cnv_score_key=cnv_score_key,
        statistics=statistics,
        visualization_available=cnv_score_key is not None,
    )


def _infer_cnv_numbat(
    data_id: str,
    adata: "ad.AnnData",
    params: CNVParameters,
    ctx: "ToolContext",
) -> CNVResult:
    """Infer copy number variations using Numbat (haplotype-aware)

    Numbat performs haplotype-aware CNV analysis by integrating allele-specific
    counts with expression data, enabling detection of copy-neutral LOH and
    reconstruction of tumor phylogeny.

    Args:
        data_id: Dataset identifier (for result creation)
        adata: AnnData object (already retrieved via ctx.get_adata)
        params: CNV analysis parameters
        ctx: Tool context for logging

    Returns:
        CNVResult containing Numbat CNV analysis results

    Raises:
        DependencyError: If the R runtime or Numbat package is unavailable
        ChatSpatialError: If dataset or parameters are invalid
    """
    # Validate allele data exists
    # Numbat requires long-format allele dataframe (from pileup_and_phase or similar)
    # Check if we have the raw allele dataframe in adata.uns
    if "numbat_allele_data_raw" in adata.uns:
        # Use pre-prepared long-format allele data
        df_allele = _validate_numbat_allele_data(adata.uns["numbat_allele_data_raw"])

    else:
        # Fallback: try to use matrix format (less ideal for Numbat)
        raise ParameterError(
            "Numbat requires long-format allele dataframe in adata.uns['numbat_allele_data_raw'].\n"
            "This should be created during data preparation (e.g., from pileup_and_phase).\n"
            "The dataframe should contain phased SNP identifiers, genetic positions, "
            "genotypes, and allele depths.\n"
            f"Available uns keys: {list(adata.uns.keys())}"
        )

    # Get raw integer count matrix (required by Numbat)
    from ..utils.adata_utils import get_raw_data_source

    raw_result = get_raw_data_source(adata, require_integer_counts=True)
    count_mat = raw_result.X

    # Prepare metadata — gene names must match count_mat dimensions
    gene_names = list(raw_result.var_names)
    cell_barcodes = list(adata.obs_names)
    if not pd.Index(cell_barcodes).map(str).is_unique:
        raise DataCompatibilityError(
            "Numbat requires unique observation names for cell-level result alignment."
        )

    # Identify reference cells (1-indexed for R)
    ref_mask = adata.obs[params.reference_key].isin(params.reference_categories)
    ref_indices_python = [i for i, is_ref in enumerate(ref_mask) if is_ref]
    ref_indices_r = [i + 1 for i in ref_indices_python]  # R is 1-indexed

    if not ref_indices_r:
        raise ParameterError(
            f"No reference cells found with key '{params.reference_key}' and "
            f"categories {params.reference_categories}"
        )

    r_env = validate_r_environment(
        ctx,
        required_packages=["numbat"],
        require_anndata2ri=True,
        package_install_commands={
            "numbat": "install.packages('numbat', dependencies = TRUE)"
        },
    )
    ro = r_env.robjects

    # Create temporary directory for Numbat output
    import os
    import shutil
    import tempfile

    out_dir = tempfile.mkdtemp(prefix="numbat_", dir=tempfile.gettempdir())

    try:
        with r_env.local_context(anndata=True, pandas=True, numpy=True) as r_context:
            r_context["count_mat"] = count_mat.T  # R expects genes × cells
            r_context["df_allele_python"] = df_allele
            r_context["gene_names"] = gene_names
            r_context["cell_barcodes"] = cell_barcodes
            r_context["ref_indices"] = ref_indices_r
            r_context["out_dir"] = out_dir

            r_context["genome"] = params.numbat_genome
            r_context["t_param"] = params.numbat_t
            r_context["max_entropy"] = params.numbat_max_entropy
            r_context["min_cells"] = params.numbat_min_cells
            r_context["ncores"] = params.numbat_ncores
            r_context["skip_nj"] = params.numbat_skip_nj

            ro.r(
                """
                    # Keep count matrix in dgCMatrix/matrix format (do NOT convert to dataframe!)
                    # run_numbat requires dgCMatrix or matrix, not data.frame
                    # Ensure proper row/column names are set
                    rownames(count_mat) = gene_names
                    colnames(count_mat) = cell_barcodes

                    # Use allele dataframe from Python (already in correct format)
                    df_allele = df_allele_python

                    # Create cell annotation for reference cells
                    # Convert cell_barcodes to character vector (rpy2 may pass it as list)
                    cell_vec = as.character(unlist(cell_barcodes))
                    cell_annot = data.frame(
                        cell = cell_vec,
                        group = ifelse(1:length(cell_vec) %in% ref_indices, "normal", "tumor"),
                        stringsAsFactors = FALSE
                    )

                    # Aggregate reference expression profile from count matrix
                    ref_profile = aggregate_counts(count_mat, cell_annot, verbose = FALSE)

                    # Run Numbat with reference profile
                    # Note: run_numbat returns "Success" string, not results object!
                    # Results are saved to out_dir as TSV/RDS files
                    tryCatch({
                        result_status = run_numbat(
                            count_mat,         # gene x cell count matrix (dgCMatrix or matrix)
                            ref_profile,       # reference expression profile (lambdas_ref)
                            df_allele,         # allele dataframe
                            genome = genome,
                            t = t_param,
                            max_entropy = max_entropy,
                            min_cells = min_cells,
                            ncores = ncores,
                            skip_nj = skip_nj,
                            plot = FALSE,
                            out_dir = out_dir,  # Output directory for results
                            verbose = FALSE
                        )
                    }, error = function(e) {
                        stop(paste("Numbat execution failed:", e$message))
                    })
                """
            )

        # Read results from output files (Numbat saves to TSV files, not R objects)
        # The suffix (e.g., _2) is the iteration number and varies by run.
        # Discover the actual iteration by finding clone_post_*.tsv files.
        clone_post_files = glob.glob(os.path.join(out_dir, "clone_post_*.tsv"))
        if not clone_post_files:
            raise DataNotFoundError(f"No Numbat clone_post output found in: {out_dir}")

        # Extract numeric iteration suffix for correct ordering
        # (string sort would put "2" after "10")
        def _iter_num(path: str) -> int:
            stem = os.path.basename(path).replace("clone_post_", "").replace(".tsv", "")
            try:
                return int(stem)
            except ValueError:
                return -1

        clone_post_file = max(clone_post_files, key=_iter_num)
        iter_suffix = str(_iter_num(clone_post_file))

        # 1. Read clone posteriors (cell-level assignments)
        clone_post = pd.read_csv(clone_post_file, sep="\t")

        # 2. Read genotype matrix (CNV states per segment)
        geno_file = _resolve_numbat_table_path(out_dir, "geno", iter_suffix)
        if geno_file is None:
            raise DataNotFoundError(
                "Numbat genotype output not found: "
                f"geno_{iter_suffix}.tsv or geno_{iter_suffix}.tsv.gz\n"
                f"Expected output files in: {out_dir}"
            )

        geno = pd.read_csv(geno_file, sep="\t")

        # 3. Read consensus segments (optional metadata)
        segs_file = _resolve_numbat_table_path(out_dir, "segs_consensus", iter_suffix)
        segs = pd.read_csv(segs_file, sep="\t") if segs_file is not None else None

        # 4. Check for the final phylogeny tree
        tree_file = os.path.join(out_dir, f"tree_final_{iter_suffix}.rds")
        has_phylo = os.path.exists(tree_file)

        normalized_clone_post, aligned_clone_post, cnv_matrix = (
            _validate_and_align_numbat_outputs(
                clone_post, geno, pd.Index(cell_barcodes)
            )
        )

        # Track cells missing from Numbat results
        n_analyzed = int(np.count_nonzero(~np.isnan(cnv_matrix).all(axis=1)))
        n_missing = len(cell_barcodes) - n_analyzed
        if n_missing > 0:
            logger.warning(
                "%d / %d cells not in Numbat output — "
                "marked as NaN/unassigned (not analyzed, not 'no CNV')",
                n_missing,
                len(cell_barcodes),
            )

        # Store results in AnnData
        adata.obsm["X_cnv_numbat"] = cnv_matrix

        # Convert numpy types to Python native types for H5AD compatibility
        # Use "unassigned" for cells not in Numbat output (distinct from clone IDs)
        adata.obs["numbat_clone"] = aligned_clone_post["clone_opt"].to_numpy()
        adata.obs["numbat_p_cnv"] = aligned_clone_post["p_cnv"].to_numpy(dtype=float)
        adata.obs["numbat_compartment"] = aligned_clone_post[
            "compartment_opt"
        ].to_numpy()

        # Store segment information if available
        if segs is not None:
            # H5AD natively supports DataFrame storage in uns
            # However, object columns with NaN values cause serialization errors
            # Fill NaN in object columns with empty string for H5AD compatibility
            segs_clean = segs.copy()
            for col in segs_clean.columns:
                if segs_clean[col].dtype == "object":
                    segs_clean[col] = segs_clean[col].fillna("")
            adata.uns["numbat_segments"] = segs_clean

        if has_phylo:
            # The RDS lives in the temporary Numbat output directory and is
            # removed in the finally block. Record what was generated without
            # publishing a path that is guaranteed to become invalid.
            adata.uns["numbat_phylogeny"] = {
                "generated": True,
                "retained": False,
                "source_filename": os.path.basename(tree_file),
                "tree_type": "tbl_graph",
            }

        # Calculate statistics
        statistics = {
            "mean_cnv": float(np.nanmean(cnv_matrix)),
            "std_cnv": float(np.nanstd(cnv_matrix)),
            "median_cnv": float(np.nanmedian(cnv_matrix)),
            "n_clones": int(normalized_clone_post["clone_opt"].nunique()),
            "mean_p_cnv": float(normalized_clone_post["p_cnv"].mean()),
            "n_reference_cells": len(ref_indices_r),
            "n_non_reference_cells": len(cell_barcodes) - len(ref_indices_r),
            "n_cells_analyzed": n_analyzed,
            "n_cells_missing": n_missing,
            "n_segments": cnv_matrix.shape[1],
        }

        # Get clone distribution
        clone_counts = normalized_clone_post["clone_opt"].value_counts()
        statistics["clone_distribution"] = {
            str(clone): int(count) for clone, count in clone_counts.items()
        }

        # Store analysis parameters
        analysis_key = _build_cnv_key(params)
        cnv_summary_key = f"cnv_analysis_{analysis_key.removeprefix('cnv_')}"
        adata.uns[cnv_summary_key] = {
            "method": "numbat",
            "reference_key": params.reference_key,
            "reference_categories": list(params.reference_categories),
            "genome": params.numbat_genome,
            "t": params.numbat_t,
            "max_entropy": params.numbat_max_entropy,
            "min_cells": params.numbat_min_cells,
            "cnv_score_key": "X_cnv_numbat",
        }

        # Build results keys for metadata
        results_keys: dict[str, list[str]] = {
            "uns": [cnv_summary_key],
            "obsm": ["X_cnv_numbat"],
            "obs": ["numbat_clone", "numbat_p_cnv", "numbat_compartment"],
        }
        if segs is not None:
            results_keys["uns"].append("numbat_segments")
        if has_phylo:
            results_keys["uns"].append("numbat_phylogeny")

        # Store metadata for scientific provenance tracking
        store_analysis_metadata(
            adata,
            analysis_name=analysis_key,
            method="numbat",
            parameters={
                "reference_key": params.reference_key,
                "reference_categories": list(params.reference_categories),
                "genome": params.numbat_genome,
                "t": params.numbat_t,
                "max_entropy": params.numbat_max_entropy,
                "min_cells": params.numbat_min_cells,
            },
            results_keys=results_keys,
            statistics=statistics,
        )

        # Export results for reproducibility
        export_analysis_result(adata, data_id, analysis_key)

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(
            f"Numbat analysis failed: {e}\n"
            "Common issues:\n"
            "  - Allele data format incompatible\n"
            "  - Missing genomic position information\n"
            "  - Insufficient reference cells\n"
            "  - R environment configuration issues"
        ) from e
    finally:
        if os.path.exists(out_dir):
            try:
                shutil.rmtree(out_dir)
            except Exception as cleanup_error:
                debug_fn = getattr(ctx, "debug", None)
                if callable(debug_fn):
                    debug_fn(f"Numbat cleanup skipped: {cleanup_error}")

    return CNVResult(
        data_id=data_id,
        method="numbat",
        reference_key=params.reference_key,
        reference_categories=list(params.reference_categories),
        n_chromosomes=0,  # Numbat doesn't report this directly
        n_genes_analyzed=len(gene_names),
        cnv_score_key="X_cnv_numbat",
        statistics=statistics,
        visualization_available=True,
    )
