"""
Enrichment analysis tools for spatial transcriptomics data.

This module provides both standard and spatially-aware enrichment analysis methods:
- Standard methods: GSEA, ORA, ssGSEA, Enrichr (via gseapy)
- Spatial methods: EnrichMap-based spatial enrichment analysis
"""

import logging
from collections.abc import Iterable
from functools import lru_cache
from typing import TYPE_CHECKING, Optional, Union

# Core dependencies (REQUIRED - in pyproject.toml dependencies)
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ..models.analysis import EnrichmentResult
from ..utils.adata_utils import get_raw_data_source, store_analysis_metadata, to_dense
from ..utils.compute import top_n_desc_indices
from ..utils.dependency_manager import require
from ..utils.exceptions import (
    ChatSpatialError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)
from ..utils.results_export import export_analysis_result

if TYPE_CHECKING:
    import anndata as ad
    import gseapy as gp

    from ..models.data import EnrichmentParameters
    from ..spatial_mcp_adapter import ToolContext


try:
    import gseapy as gp  # optional dependency, in [full]
except ImportError:
    gp = None  # deferred until first use


def _get_gseapy():
    """Return gseapy module, raising a clear error if not installed."""
    if gp is None:
        return require(
            "gseapy",
            feature="gene set enrichment analysis",
        )
    return gp


# Subramanian et al. 2005 PNAS: "We recommend an FDR cutoff of 25%"
_GSEA_FDR_THRESHOLD = 0.25
_SPATIAL_ENRICHMAP_MAX_DATABASE_GENE_SETS = 50
_LATEST_ENRICHMENT_RESULTS_KEY = "enrichment_latest_results_key"

logger = logging.getLogger(__name__)


def _build_enrichment_key(method: str, database: str | None) -> str:
    """Build parametric analysis key for enrichment.

    Encodes method + database so multiple enrichment runs coexist.
    """
    if database:
        # Sanitize database name for use as key (replace spaces, slashes)
        db_clean = database.replace(" ", "_").replace("/", "_")
        return f"enrichment_{method}_{db_clean}"
    return f"enrichment_{method}"


def _store_enrichment_results(
    adata: "ad.AnnData",
    result_key: str,
    results: object,
    method: str,
    database: str | None,
    parameters: dict[str, object],
    statistics: dict[str, object],
    species: str | None = None,
    gene_sets: dict[str, list[str]] | None = None,
    gene_sets_shared_key: str | None = None,
    data_id: str | None = None,
) -> tuple[str, str]:
    """Store tabular enrichment results with shared, per-run, and metadata keys."""
    adata.uns[result_key] = results
    adata.uns[_LATEST_ENRICHMENT_RESULTS_KEY] = result_key

    analysis_key = _build_enrichment_key(method, database)
    per_run_key = f"{result_key}_{analysis_key.removeprefix('enrichment_')}"
    adata.uns[per_run_key] = results

    uns_keys = [per_run_key]
    if gene_sets is not None and gene_sets_shared_key is not None:
        adata.uns[gene_sets_shared_key] = gene_sets
        per_run_gs_key = (
            f"{gene_sets_shared_key}_{analysis_key.removeprefix('enrichment_')}"
        )
        adata.uns[per_run_gs_key] = gene_sets
        uns_keys.append(per_run_gs_key)

    store_analysis_metadata(
        adata,
        analysis_name=analysis_key,
        method=method,
        parameters=parameters,
        results_keys={"uns": uns_keys},
        statistics=statistics,
        species=species,
        database=database,
    )
    if data_id is not None:
        export_analysis_result(adata, data_id, analysis_key)

    return analysis_key, per_run_key


# Preferred/legacy library candidates to keep versioning deterministic while
# remaining backward-compatible with environments missing the latest release.
GENESET_LIBRARY_CANDIDATES: dict[str, tuple[str, ...]] = {
    "GO_Biological_Process": (
        "GO_Biological_Process_2025",
        "GO_Biological_Process_2023",
    ),
    "GO_Molecular_Function": (
        "GO_Molecular_Function_2025",
        "GO_Molecular_Function_2023",
    ),
    "GO_Cellular_Component": (
        "GO_Cellular_Component_2025",
        "GO_Cellular_Component_2023",
    ),
    "Reactome_Pathways": ("Reactome_Pathways_2024", "Reactome_2022"),
    "MSigDB_Hallmark": ("MSigDB_Hallmark_2020",),
    "Cell_Type_Markers": ("CellMarker_Augmented_2021",),
}
GO_ASPECT_TO_DATABASE_KEY: dict[str, str] = {
    "BP": "GO_Biological_Process",
    "MF": "GO_Molecular_Function",
    "CC": "GO_Cellular_Component",
}

# Centralized Enrichr library mapping to keep year/version consistent.
ENRICHR_LIBRARY_MAP: dict[str, str] = {
    key: candidates[0] for key, candidates in GENESET_LIBRARY_CANDIDATES.items()
}
DEFAULT_ENRICHR_DATABASES: tuple[str, ...] = (
    "GO_Biological_Process",
    "GO_Molecular_Function",
    "GO_Cellular_Component",
    "KEGG_Pathways",
    "Reactome_Pathways",
    "MSigDB_Hallmark",
)
ENRICHR_ORGANISM_DATABASE_ALIASES: dict[str, tuple[str, ...]] = {
    "Fly": ("fly", "d. melanogaster", "drosophila melanogaster"),
    "Yeast": ("yeast", "s. cerevisiae", "saccharomyces cerevisiae"),
    "Worm": ("worm", "c. elegans", "caenorhabditis elegans", "nematode"),
    "Fish": ("fish", "d. rerio", "danio rerio", "zebrafish"),
}


def _top_n_desc_indices(values: np.ndarray, n_top: int) -> np.ndarray:
    """Backward-compatible wrapper with non-finite sanitization."""
    return top_n_desc_indices(values, n_top, sanitize_nonfinite=True)


def _normalize_gene_set_library(raw_gene_sets: object) -> dict[str, list[str]]:
    """Normalize gene-set library payloads from gseapy/Enrichr boundaries."""
    if not isinstance(raw_gene_sets, dict):
        raise ProcessingError(
            f"Gene set library must be a dictionary, got {type(raw_gene_sets).__name__}"
        )

    normalized: dict[str, list[str]] = {}
    for raw_name, raw_genes in raw_gene_sets.items():
        name = raw_name.decode("utf-8") if isinstance(raw_name, bytes) else str(raw_name)
        if isinstance(raw_genes, bytes | str):
            genes_iter: Iterable[object] = [raw_genes]
        else:
            genes_iter = raw_genes

        genes = [
            gene.decode("utf-8") if isinstance(gene, bytes) else str(gene)
            for gene in genes_iter
            if gene not in (None, b"", "")
        ]
        normalized[name] = genes

    return normalized


def _enrichr_database_for_organism(organism: str) -> str:
    organism_lower = organism.lower()
    if organism_lower in {
        "human",
        "mouse",
        "hs",
        "mm",
        "homo sapiens",
        "mus musculus",
        "h. sapiens",
        "m. musculus",
    }:
        return "Enrichr"

    for key, aliases in ENRICHR_ORGANISM_DATABASE_ALIASES.items():
        if organism_lower in aliases:
            return f"{key}Enrichr"

    raise LookupError(
        "No supported database. Please input one of: "
        "Human, Mouse, Yeast, Fly, Fish, Worm"
    )


@lru_cache(maxsize=None)
def _get_enrichr_library_names(organism: str) -> tuple[str, ...]:
    return tuple(_get_gseapy().get_library_name(organism=organism))


def _is_gseapy_bytes_decode_error(error: TypeError) -> bool:
    return "bytes-like object" in str(error)


def _download_enrichr_library(library_name: str, organism: str) -> dict[str, list[str]]:
    """Download an Enrichr library without relying on gseapy's line decoding."""
    import requests

    if library_name not in _get_enrichr_library_names(organism):
        raise ValueError(
            f"Library {library_name!r} was not found for organism {organism!r}"
        )

    parser = getattr(_get_gseapy(), "parser", None)
    url_root = getattr(parser, "ENRICHR_URL", "http://maayanlab.cloud")
    url = f"{url_root}/{_enrichr_database_for_organism(organism)}/geneSetLibrary"

    with requests.get(
        url,
        params={"mode": "text", "libraryName": library_name},
        timeout=(10, 60),
        stream=True,
    ) as response:
        try:
            response.raise_for_status()
        except requests.HTTPError as e:
            raise RuntimeError(
                f"Error fetching gene set library {library_name!r} for {organism!r}."
            ) from e

        gene_sets: dict[str, list[str]] = {}
        for raw_line in response.iter_lines(chunk_size=1024, decode_unicode=True):
            if not raw_line:
                continue
            line = raw_line.decode("utf-8") if isinstance(raw_line, bytes) else raw_line
            entries = line.strip().split("\t")
            if not entries:
                continue
            gene_sets[entries[0]] = [
                entry.split(",")[0] for entry in entries[2:] if entry
            ]

    return gene_sets


def _load_enrichr_library(library_name: str, organism: str) -> dict[str, list[str]]:
    try:
        return _normalize_gene_set_library(
            _get_gseapy().get_library(library_name, organism=organism)
        )
    except TypeError as e:
        if not _is_gseapy_bytes_decode_error(e):
            raise
        return _normalize_gene_set_library(
            _download_enrichr_library(library_name, organism)
        )


def _load_library_first_available(
    database_key: str, organism: str
) -> dict[str, list[str]]:
    """Load the first available library from deterministic candidates."""
    if database_key not in GENESET_LIBRARY_CANDIDATES:
        raise ParameterError(f"Unknown database key: {database_key}")

    last_error: Exception | None = None
    for library_name in GENESET_LIBRARY_CANDIDATES[database_key]:
        try:
            return _load_enrichr_library(library_name, organism)
        except ChatSpatialError:
            raise
        except Exception as e:  # pragma: no cover - exercised via caller-level wrapping
            last_error = e
            continue

    if last_error is not None:
        raise last_error
    raise RuntimeError(
        f"No library candidates configured for database key: {database_key}"
    )


# ============================================================================
# MCP RESPONSE OPTIMIZATION
# ============================================================================


def _filter_significant_statistics(
    gene_set_statistics: dict,
    enrichment_scores: dict,
    pvalues: dict,
    adjusted_pvalues: dict,
    method: str,
    fdr_threshold: Optional[float] = None,
) -> tuple[dict, dict, dict, dict]:
    """
    Filter all enrichment result dictionaries to only include significant pathways.

    This dramatically reduces MCP response size for large gene set databases
    (e.g., KEGG 311 pathways, GO 10,000 terms) while preserving all important
    information for users.

    Args:
        gene_set_statistics: Full statistics for all gene sets
        enrichment_scores: Enrichment scores for all gene sets
        pvalues: P-values for all gene sets
        adjusted_pvalues: FDR-corrected p-values for all gene sets
        method: Enrichment method used ('gsea', 'ora', 'enrichr', 'ssgsea')
        fdr_threshold: FDR threshold for significance (default: None for method-based auto)
                       Method-based defaults (based on statistical best practices):
                       - GSEA: FDR < 0.25 (official recommendation from Subramanian et al. 2005)
                       - ORA/Enrichr: FDR < 0.05 (standard statistical threshold)
                       - ssGSEA: No filtering (no p-values produced)

    Returns:
        Tuple of (filtered_statistics, filtered_scores, filtered_pvals, filtered_adj_pvals)

    Example:
        Before: 311 pathways × 4 dicts × 100 chars = 124KB (KEGG)
        After: ~15 significant pathways × 4 dicts × 100 chars = 6KB (95% reduction)

    References:
        - GSEA: Subramanian et al. (2005) PNAS 102(43):15545-15550
          "We recommend an FDR cutoff of 25% when dealing with a single database"
        - ORA: Standard multiple testing correction threshold (Benjamini & Hochberg 1995)
    """
    if not adjusted_pvalues:
        # No p-values available (e.g., ssGSEA), return all results without filtering
        return gene_set_statistics, enrichment_scores, pvalues, adjusted_pvalues

    # Auto-determine threshold based on ANALYSIS METHOD if not specified
    # This is statistically principled: different methods have different FDR standards
    if fdr_threshold is None:
        method_lower = method.lower()
        if method_lower == "gsea":
            fdr_threshold = _GSEA_FDR_THRESHOLD
        elif method_lower in ("ora", "enrichr", "pathway_ora", "pathway_enrichr"):
            # ORA and Enrichr: standard statistical threshold
            # Based on Benjamini-Hochberg FDR control at 5%
            fdr_threshold = 0.05
        else:
            # Default fallback for unknown methods
            fdr_threshold = 0.05

    # Find significant pathways
    significant = {
        name
        for name, fdr in adjusted_pvalues.items()
        if fdr is not None and fdr < fdr_threshold
    }

    # Filter all dictionaries
    filtered_stats = {
        name: stats
        for name, stats in gene_set_statistics.items()
        if name in significant
    }

    filtered_scores = {
        name: score for name, score in enrichment_scores.items() if name in significant
    }

    filtered_pvals = {
        name: pval for name, pval in pvalues.items() if name in significant
    }

    filtered_adj_pvals = {
        name: adj_pval
        for name, adj_pval in adjusted_pvalues.items()
        if name in significant
    }

    return filtered_stats, filtered_scores, filtered_pvals, filtered_adj_pvals


# ============================================================================
# GENE SET UTILITIES
# ============================================================================


def _filter_gene_sets_by_size(
    gene_sets: dict[str, list[str]], min_size: int, max_size: int
) -> dict[str, list[str]]:
    """Filter gene sets by size constraints.

    Args:
        gene_sets: Dictionary mapping gene set names to gene lists.
        min_size: Minimum number of genes required.
        max_size: Maximum number of genes allowed.

    Returns:
        Filtered gene sets within size constraints.
    """
    return {
        name: genes
        for name, genes in gene_sets.items()
        if min_size <= len(genes) <= max_size
    }


def _match_gene_set_to_dataset(
    genes: list[str], available_genes: set[str], species: str
) -> list[str]:
    """Return genes from a gene set that are present in the dataset."""
    common_genes = [gene for gene in genes if gene in available_genes]
    if len(common_genes) < len(genes) * 0.5 and species != "unknown":
        dataset_format_genes, _ = _convert_gene_format_for_matching(
            genes, available_genes, species
        )
        if len(dataset_format_genes) > len(common_genes):
            common_genes = dataset_format_genes
    return common_genes


def _rank_gene_sets_by_dataset_overlap(
    gene_sets: dict[str, list[str]], available_genes: set[str], species: str
) -> list[tuple[str, list[str]]]:
    """Rank gene sets by overlap with the dataset for bounded spatial scoring."""
    ranked: list[tuple[int, int, str, list[str]]] = []
    for name, genes in gene_sets.items():
        common_genes = _match_gene_set_to_dataset(genes, available_genes, species)
        if len(common_genes) >= 2:
            ranked.append((len(common_genes), len(genes), name, genes))

    ranked.sort(key=lambda item: (-item[0], item[1], item[2]))
    return [(name, genes) for _, _, name, genes in ranked]


def _limit_spatial_enrichmap_gene_sets(
    gene_sets: dict[str, list[str]],
    available_genes: set[str],
    species: str,
    max_gene_sets: int = _SPATIAL_ENRICHMAP_MAX_DATABASE_GENE_SETS,
) -> dict[str, list[str]]:
    """Limit database-backed spatial EnrichMap scoring to best-overlap signatures."""
    if len(gene_sets) <= max_gene_sets:
        return gene_sets

    ranked_gene_sets = _rank_gene_sets_by_dataset_overlap(gene_sets, available_genes, species)
    return dict(ranked_gene_sets[:max_gene_sets])


# ============================================================================
# SPARSE MATRIX UTILITIES
# ============================================================================


def _compute_std_sparse_compatible(X, axis: int = 0, ddof: int = 1) -> np.ndarray:
    """
    Compute standard deviation compatible with both dense and sparse matrices.

    For sparse matrices, uses the formula: std = sqrt(E[X^2] - E[X]^2) with Bessel correction.
    For dense matrices, uses numpy's built-in std method.

    Args:
        X: Input matrix (can be sparse or dense)
        axis: Axis along which to compute std (0 for columns, 1 for rows)
        ddof: Delta Degrees of Freedom for Bessel correction (default: 1)

    Returns:
        1D numpy array of standard deviations
    """
    import scipy.sparse as sp

    if sp.issparse(X):
        # Sparse matrix: use mathematical formula
        n = X.shape[axis]
        mean = np.array(X.mean(axis=axis)).flatten()
        mean_of_squares = np.array(X.power(2).mean(axis=axis)).flatten()

        # Compute variance with Bessel correction: n/(n-ddof)
        variance = mean_of_squares - np.power(mean, 2)
        variance = np.maximum(variance, 0)  # Avoid numerical errors
        if ddof > 0:
            variance = variance * n / (n - ddof)  # Bessel correction

        return np.sqrt(variance)
    else:
        # Dense matrix: use numpy's built-in method
        return np.array(X.std(axis=axis, ddof=ddof)).flatten()


def _compute_cv_ranking(X, var_names: pd.Index) -> dict[str, float]:
    """Compute coefficient-of-variation ranking for genes."""
    mean = np.array(X.mean(axis=0)).flatten()
    std = _compute_std_sparse_compatible(X, axis=0, ddof=1)
    cv = np.zeros_like(mean)
    nonzero_mask = np.abs(mean) > 1e-10
    cv[nonzero_mask] = std[nonzero_mask] / np.abs(mean[nonzero_mask])
    return dict(zip(var_names, cv, strict=True))


def _compute_variance_ranking(X, var_names: pd.Index) -> dict[str, float]:
    """Compute variance-based ranking for genes."""
    import scipy.sparse as sp

    if sp.issparse(X):
        mean = np.array(X.mean(axis=0)).flatten()
        mean_sq = np.array(X.power(2).mean(axis=0)).flatten()
        var = mean_sq - np.power(mean, 2)
    else:
        var = np.array(X.var(axis=0)).flatten()
    return dict(zip(var_names, var, strict=True))


# ============================================================================
# GENE FORMAT CONVERSION UTILITIES
# ============================================================================


def _convert_gene_format_for_matching(
    pathway_genes: list[str], dataset_genes: set, species: str
) -> tuple[list[str], dict[str, str]]:
    """
    Rule-based gene format conversion to match dataset format.

    Handles common gene format variations between pathway databases and datasets:
    - Uppercase (GENE) vs Title case (Gene) vs lowercase (gene)
    - Species-specific formatting rules
    - Special prefixes like Gm/GM/gm for mouse genes

    Args:
        pathway_genes: Gene names from pathway database (usually uppercase from gseapy)
        dataset_genes: Available gene names in dataset
        species: Species specified by user ("mouse" or "human")

    Returns:
        (dataset_format_genes, conversion_map)
        dataset_format_genes: Gene names in dataset format that can be found
        conversion_map: Maps dataset_format -> original_pathway_format
    """
    dataset_format_genes = []
    conversion_map = {}

    for gene in pathway_genes:
        # Try direct match first
        if gene in dataset_genes:
            dataset_format_genes.append(gene)
            conversion_map[gene] = gene
            continue

        # Apply multiple format conversion rules
        format_variations = []

        if species == "mouse":
            # Mouse-specific format rules (order matters for efficiency)
            # Rule 1: Title case (most common): Cd5l, Gbp2b
            if len(gene) > 1:
                format_variations.append(gene[0].upper() + gene[1:].lower())
            # Rule 2: All lowercase: cd5l, gbp2b
            format_variations.append(gene.lower())
            # Rule 3: All uppercase: CD5L, GBP2B
            format_variations.append(gene.upper())
            # Rule 4: Capitalize first letter only
            format_variations.append(gene.capitalize())

            # Special rule for Gm-prefixed genes (common in mouse)
            if gene.upper().startswith("GM"):
                format_variations.extend(
                    [
                        "gm" + gene[2:].lower(),  # gm42418
                        "Gm" + gene[2:].lower(),  # Gm42418
                        "GM" + gene[2:].upper(),  # GM42418
                    ]
                )

        elif species == "human":
            # Human-specific format rules
            # Rule 1: All uppercase (most common): HES1, FABP4
            format_variations.append(gene.upper())
            # Rule 2: All lowercase: hes1, fabp4
            format_variations.append(gene.lower())
            # Rule 3: Capitalize first letter
            format_variations.append(gene.capitalize())

        # Remove duplicates while preserving order
        seen = set()
        unique_variations = []
        for variation in format_variations:
            if variation not in seen and variation != gene:  # Skip if same as original
                seen.add(variation)
                unique_variations.append(variation)

        # Try each format variation against dataset
        for variant in unique_variations:
            if variant in dataset_genes:
                dataset_format_genes.append(variant)  # Use dataset's actual format
                conversion_map[variant] = gene
                break  # Stop after first match

    return dataset_format_genes, conversion_map


# ============================================================================
# ENRICHR DATABASE MAPPING
# ============================================================================


def map_gene_set_database_to_enrichr_library(database_name: str, species: str) -> str:
    """Map user-friendly database names to actual Enrichr library names.

    Args:
        database_name: User-friendly database name from MCP interface
        species: Species ('human', 'mouse', or 'zebrafish')

    Returns:
        Actual Enrichr library name

    Raises:
        ValueError: If database_name is not supported
    """
    if database_name == "KEGG_Pathways":
        species_lower = species.lower()
        if species_lower == "human":
            return "KEGG_2021_Human"
        elif species_lower == "mouse":
            return "KEGG_2019_Mouse"
        else:
            raise ParameterError(
                f"KEGG_Pathways is only available for 'human' and 'mouse', "
                f"not '{species}'. Use GO or Reactome databases instead."
            )

    if database_name not in ENRICHR_LIBRARY_MAP:
        available_options = sorted([*ENRICHR_LIBRARY_MAP.keys(), "KEGG_Pathways"])
        raise ParameterError(
            f"Unknown gene set database: {database_name}. "
            f"Available options: {available_options}"
        )

    return ENRICHR_LIBRARY_MAP[database_name]


# ============================================================================
# ENRICHMENT ANALYSIS FUNCTIONS
# ============================================================================


def _get_numeric_var_columns(adata: "ad.AnnData") -> list[str]:
    return [
        col for col in adata.var.columns if pd.api.types.is_numeric_dtype(adata.var[col])
    ]


def _validate_gsea_ranking_key(adata: "ad.AnnData", ranking_key: str) -> None:
    if ranking_key in adata.var and pd.api.types.is_numeric_dtype(adata.var[ranking_key]):
        return

    available_rankings = _get_numeric_var_columns(adata)
    raise ParameterError(
        f"GSEA ranking key '{ranking_key}' was provided but is not a numeric "
        "gene-level ranking in adata.var. Provide a numeric adata.var column, or omit "
        "score_keys to let ChatSpatial compute a default expression-based ranking. "
        f"Available numeric adata.var columns: {available_rankings}"
    )


def perform_gsea(
    adata: "ad.AnnData",
    gene_sets: dict[str, list[str]],
    ranking_key: Optional[str] = None,
    method: str = "signal_to_noise",
    permutation_num: int = 1000,
    min_size: int = 10,
    max_size: int = 500,
    pvalue_cutoff: float = 0.25,
    species: Optional[str] = None,
    database: Optional[str] = None,
    ctx: Optional["ToolContext"] = None,
    data_id: Optional[str] = None,
) -> "EnrichmentResult":
    """Perform Gene Set Enrichment Analysis (GSEA).

    Args:
        adata: Annotated data matrix.
        gene_sets: Gene sets to test (name -> gene list).
        ranking_key: Key in adata.var for pre-computed ranking. Computes if None.
        method: Method for ranking genes if ranking_key is None.
        permutation_num: Number of permutations for p-value calculation.
        min_size: Minimum gene set size.
        max_size: Maximum gene set size.
        species: Species for the analysis ('mouse', 'human').
        database: Gene set database used ('KEGG_Pathways', 'GO_Biological_Process').
        ctx: MCP tool context for logging.
        data_id: Dataset identifier for result tracking.

    Returns:
        EnrichmentResult with enrichment scores and statistics.
    """
    # gseapy imported at module level (required dependency)

    ranking_method = method.strip().lower()

    # Prepare ranking
    if ranking_key:
        _validate_gsea_ranking_key(adata, ranking_key)
        ranking = adata.var[ranking_key].to_dict()
    else:
        # Compute ranking from expression data
        # Use get_raw_data_source (single source of truth) for complete gene coverage
        raw_result = get_raw_data_source(adata, prefer_complete_genes=True)
        X = raw_result.X
        var_names = raw_result.var_names

        # Compute gene ranking metric
        # IMPORTANT: GSEA requires biologically meaningful ranking, not just variance
        # Reference: Subramanian et al. (2005) PNAS, GSEA-MSIGDB documentation

        if ranking_method not in {
            "signal_to_noise",
            "coefficient_of_variation",
            "variance",
            "highly_variable_rank",
            "dispersions_norm",
        }:
            raise ParameterError(
                "Unsupported GSEA ranking method: "
                f"{method}. Supported: signal_to_noise, coefficient_of_variation, "
                "variance, highly_variable_rank, dispersions_norm"
            )

        if ranking_method == "signal_to_noise":
            group_key = "condition" if "condition" in adata.obs else "group"
            groups = adata.obs[group_key].unique() if group_key in adata.obs else []

            if len(groups) == 2:
                # Binary comparison: Use Signal-to-Noise Ratio (GSEA default)
                # S2N = (μ1 - μ2) / (σ1 + σ2)
                # This captures both differential expression AND expression stability
                group1_mask = adata.obs[group_key] == groups[0]
                group2_mask = adata.obs[group_key] == groups[1]

                # Compute means
                mean1 = np.array(X[group1_mask, :].mean(axis=0)).flatten()
                mean2 = np.array(X[group2_mask, :].mean(axis=0)).flatten()

                # Compute standard deviations (sparse-compatible)
                std1 = _compute_std_sparse_compatible(X[group1_mask, :], axis=0, ddof=1)
                std2 = _compute_std_sparse_compatible(X[group2_mask, :], axis=0, ddof=1)

                # Apply minimum std threshold (GSEA standard: 0.2 * |mean|)
                # This prevents division by zero and reduces noise from low-variance genes
                min_std_factor = 0.2
                std1 = np.maximum(std1, min_std_factor * np.abs(mean1))
                std2 = np.maximum(std2, min_std_factor * np.abs(mean2))

                # Compute Signal-to-Noise Ratio
                s2n = (mean1 - mean2) / (std1 + std2)
                ranking = dict(zip(var_names, s2n, strict=True))
            else:
                # S2N requires exactly 2 groups; fall back to CV when unavailable.
                ranking = _compute_cv_ranking(X, var_names)
        elif ranking_method == "coefficient_of_variation":
            ranking = _compute_cv_ranking(X, var_names)
        elif ranking_method == "variance":
            ranking = _compute_variance_ranking(X, var_names)
        elif ranking_method == "highly_variable_rank":
            if "highly_variable_rank" not in adata.var:
                raise DataNotFoundError(
                    "Ranking method 'highly_variable_rank' requested but "
                    "adata.var['highly_variable_rank'] is missing."
                )
            ranking = adata.var["highly_variable_rank"].to_dict()
        elif ranking_method == "dispersions_norm":
            if "dispersions_norm" not in adata.var:
                raise DataNotFoundError(
                    "Ranking method 'dispersions_norm' requested but "
                    "adata.var['dispersions_norm'] is missing."
                )
            ranking = adata.var["dispersions_norm"].to_dict()
        else:  # pragma: no cover - defensive fallback, guarded by method validation
            # Defensive fallback (should be unreachable due validation above)
            ranking = _compute_cv_ranking(X, var_names)

    # Run GSEA preranked
    try:
        # Convert ranking dict to DataFrame for gseapy
        ranking_df = pd.DataFrame.from_dict(ranking, orient="index", columns=["score"])
        ranking_df.index.name = "gene"
        ranking_df = ranking_df.sort_values("score", ascending=False)

        res = _get_gseapy().prerank(
            rnk=ranking_df,  # Pass DataFrame instead of dict
            gene_sets=gene_sets,
            processes=1,
            permutation_num=permutation_num,
            min_size=min_size,
            max_size=max_size,
            seed=42,
            verbose=False,
            no_plot=True,
            outdir=None,
        )

        # Extract results
        results_df = res.res2d

        # Prepare output - OPTIMIZED: vectorized dict + array iteration (16x faster)
        enrichment_scores = dict(
            zip(results_df["Term"], results_df["ES"], strict=True)
        )
        pvalues = dict(
            zip(results_df["Term"], results_df["NOM p-val"], strict=True)
        )
        adjusted_pvalues = dict(
            zip(results_df["Term"], results_df["FDR q-val"], strict=True)
        )

        # Pre-extract arrays for fast iteration
        terms = results_df["Term"].values
        es_vals = results_df["ES"].values
        nes_vals = results_df["NES"].values
        pval_vals = results_df["NOM p-val"].values
        fdr_vals = results_df["FDR q-val"].values
        has_matched_size = "Matched_size" in results_df.columns
        has_lead_genes = "Lead_genes" in results_df.columns
        size_vals = (
            results_df["Matched_size"].values
            if has_matched_size
            else np.zeros(len(terms))
        )
        lead_genes_vals = (
            results_df["Lead_genes"].values if has_lead_genes else [""] * len(terms)
        )

        gene_set_statistics = {}
        for i in range(len(terms)):
            lead_genes_str = lead_genes_vals[i] if has_lead_genes else ""
            gene_set_statistics[terms[i]] = {
                "es": float(es_vals[i]),
                "nes": float(nes_vals[i]),
                "pval": float(pval_vals[i]),
                "fdr": float(fdr_vals[i]),
                "size": int(size_vals[i]) if has_matched_size else 0,
                "lead_genes": lead_genes_str.split(";")[:10] if lead_genes_str else [],
            }

        # Get top enriched and depleted
        results_df_sorted = results_df.sort_values("NES", ascending=False)
        top_enriched = (
            results_df_sorted[results_df_sorted["NES"] > 0].head(10)["Term"].tolist()
        )
        top_depleted = (
            results_df_sorted[results_df_sorted["NES"] < 0].head(10)["Term"].tolist()
        )

        _store_enrichment_results(
            adata=adata,
            result_key="gsea_results",
            results=results_df,
            method="gsea",
            database=database,
            parameters={
                "permutation_num": permutation_num,
                "ranking_method": ranking_method,
                "min_size": min_size,
                "max_size": max_size,
                "ranking_key": ranking_key,
            },
            statistics={
                "n_gene_sets": len(gene_sets),
                "n_significant": int((results_df["FDR q-val"] < pvalue_cutoff).sum()),
            },
            species=species,
            gene_sets=gene_sets,
            gene_sets_shared_key="enrichment_gsea_gene_sets",
            data_id=data_id,
        )

        # Filter all result dictionaries to only significant pathways (reduces MCP response size)
        # Uses method-based FDR threshold: GSEA = 0.25 (Subramanian et al. 2005)
        (
            filtered_statistics,
            filtered_scores,
            filtered_pvals,
            filtered_adj_pvals,
        ) = _filter_significant_statistics(
            gene_set_statistics,
            enrichment_scores,
            pvalues,
            adjusted_pvalues,
            method="gsea",
            fdr_threshold=pvalue_cutoff,
        )

        return EnrichmentResult(
            method="gsea",
            n_gene_sets=len(gene_sets),
            n_significant=len(results_df[results_df["FDR q-val"] < pvalue_cutoff]),
            enrichment_scores=filtered_scores,
            pvalues=filtered_pvals,
            adjusted_pvalues=filtered_adj_pvals,
            gene_set_statistics=filtered_statistics,
            top_gene_sets=top_enriched,
            top_depleted_sets=top_depleted,
        )

    except Exception as e:
        logger.error(f"GSEA failed: {e}")
        raise


def perform_ora(
    adata: "ad.AnnData",
    gene_sets: dict[str, list[str]],
    gene_list: Optional[list[str]] = None,
    pvalue_threshold: float = 0.05,
    significance_threshold: Optional[float] = None,
    min_size: int = 10,
    max_size: int = 500,
    adjust_method: str = "fdr",
    species: Optional[str] = None,
    database: Optional[str] = None,
    ctx: Optional["ToolContext"] = None,
    data_id: Optional[str] = None,
) -> "EnrichmentResult":
    """Perform Over-Representation Analysis (ORA).

    Args:
        adata: Annotated data matrix.
        gene_sets: Gene sets to test (name -> gene list).
        gene_list: Genes to test. Uses DEGs from rank_genes_groups if None.
        pvalue_threshold: P-value threshold for selecting DEGs from
            rank_genes_groups (only used when gene_list is None).
        significance_threshold: Adjusted p-value threshold for counting and
            filtering significant pathways in ORA results. If None, defaults
            to 0.05 (standard Benjamini-Hochberg FDR control).
        min_size: Minimum gene set size.
        max_size: Maximum gene set size.
        adjust_method: Multiple testing correction ('fdr', 'bonferroni', 'none').
        species: Species for the analysis ('mouse', 'human').
        database: Gene set database used ('KEGG_Pathways', 'GO_Biological_Process').
        ctx: MCP tool context for logging.
        data_id: Dataset identifier for result tracking.

    Returns:
        EnrichmentResult with enrichment scores and statistics.

    Note:
        LogFC filtering removed. ORA should use genes pre-filtered by find_markers.
        Gene filtering is the responsibility of differential expression analysis.
    """
    # Default significance threshold for ORA: standard FDR 5%
    if significance_threshold is None:
        significance_threshold = 0.05
    # Get gene list if not provided
    if gene_list is None:
        # Try to get DEGs from adata
        if "rank_genes_groups" in adata.uns:
            # Get DEGs
            result = adata.uns["rank_genes_groups"]
            names = result["names"]

            # Check if pvals exist (not all rank_genes_groups have pvals)
            pvals = None
            if "pvals_adj" in result:
                pvals = result["pvals_adj"]
            elif "pvals" in result:
                pvals = result["pvals"]

            # Get DEGs from all groups and merge
            # IMPORTANT: names is a numpy recarray with shape (n_genes,)
            # and dtype.names contains group names as fields
            # Access genes by group name: names[group_name][i]
            degs_seen: set[str] = set()
            degs_ordered: list[str] = []

            # Iterate over all groups
            for group_name in names.dtype.names:
                for i in range(len(names)):
                    # Skip genes that don't pass filter criteria
                    if pvals is not None and pvals[group_name][i] >= pvalue_threshold:
                        continue
                    if pvals is None and i >= 100:  # Top 100 genes when no pvals
                        continue

                    gene_name = str(names[group_name][i])
                    if gene_name not in degs_seen:
                        degs_seen.add(gene_name)
                        degs_ordered.append(gene_name)

            gene_list = degs_ordered
        else:
            # Use highly variable genes
            if "highly_variable" in adata.var:
                gene_list = adata.var_names[adata.var["highly_variable"]].tolist()
            else:
                # Use top variable genes (based on Coefficient of Variation)
                # CV = σ/μ is more appropriate than raw variance
                mean = np.array(adata.X.mean(axis=0)).flatten()
                std = _compute_std_sparse_compatible(adata.X, axis=0, ddof=1)

                # Compute CV (avoid division by zero)
                cv = np.zeros_like(mean)
                nonzero_mask = np.abs(mean) > 1e-10
                cv[nonzero_mask] = std[nonzero_mask] / np.abs(mean[nonzero_mask])

                top_indices = top_n_desc_indices(cv, 500, sanitize_nonfinite=True)
                gene_list = adata.var_names[top_indices].tolist()

    # Background genes
    # Use get_raw_data_source (single source of truth) to get complete gene set
    # This handles gene name casing differences between raw and filtered data
    bg_result = get_raw_data_source(adata, prefer_complete_genes=True)
    background_genes = set(bg_result.var_names)

    # Case-insensitive matching as fallback for gene name format differences
    # (e.g., MT.CO1 vs MT-CO1, uppercase vs lowercase)
    query_genes = set(gene_list) & background_genes

    # If no direct matches, try case-insensitive matching
    if len(query_genes) == 0 and len(gene_list) > 0:
        # Create case-insensitive lookup
        gene_name_map = {g.upper(): g for g in background_genes}
        query_genes = set()
        for gene in gene_list:
            if gene.upper() in gene_name_map:
                query_genes.add(gene_name_map[gene.upper()])

    # Perform hypergeometric test for each gene set
    enrichment_scores = {}
    pvalues = {}
    gene_set_statistics = {}

    for gs_name, gs_genes in gene_sets.items():
        gs_genes_set = set(gs_genes) & background_genes

        if len(gs_genes_set) < min_size or len(gs_genes_set) > max_size:
            continue

        # Hypergeometric test
        # a: genes in both query and gene set
        # b: genes in query but not in gene set
        # c: genes in gene set but not in query
        # d: genes in neither

        a = len(query_genes & gs_genes_set)
        b = len(query_genes - gs_genes_set)
        c = len(gs_genes_set - query_genes)
        d = len(background_genes - query_genes - gs_genes_set)

        # Fisher's exact test
        odds_ratio, p_value = stats.fisher_exact(
            [[a, b], [c, d]], alternative="greater"
        )

        enrichment_scores[gs_name] = odds_ratio
        pvalues[gs_name] = p_value

        gene_set_statistics[gs_name] = {
            "odds_ratio": odds_ratio,
            "pval": p_value,
            "overlap": a,
            "query_size": len(query_genes),
            "gs_size": len(gs_genes_set),
            "overlapping_genes": list(query_genes & gs_genes_set)[:20],  # Top 20
        }

    # Multiple testing correction
    _method_map = {"fdr": "fdr_bh", "bonferroni": "bonferroni"}
    if pvalues and adjust_method != "none":
        stats_method = _method_map.get(adjust_method, "fdr_bh")
        pval_array = np.array(list(pvalues.values()))
        _, adjusted_pvals, _, _ = multipletests(pval_array, method=stats_method)
        adjusted_pvalues = dict(zip(pvalues.keys(), adjusted_pvals, strict=True))
    else:
        # No correction: adjusted = raw
        adjusted_pvalues = dict(pvalues) if pvalues else {}

    # Get top results
    sorted_by_pval = sorted(pvalues.items(), key=lambda x: x[1])
    top_gene_sets = [x[0] for x in sorted_by_pval[:10]]

    # Save results to adata.uns for visualization
    # Create DataFrame for visualization compatibility
    ora_df = pd.DataFrame(
        {
            "pathway": list(enrichment_scores),
            "odds_ratio": list(enrichment_scores.values()),
            "pvalue": [pvalues.get(k, 1.0) for k in enrichment_scores],
            "adjusted_pvalue": [
                adjusted_pvalues.get(k, 1.0) for k in enrichment_scores
            ],
        }
    )
    ora_df["NES"] = ora_df["odds_ratio"]  # Use odds_ratio as score for visualization
    ora_df = ora_df.sort_values("pvalue")

    # Count significant pathways using the pathway significance threshold
    # (distinct from the DEG selection threshold used for input gene filtering)
    n_significant = sum(
        1
        for p in adjusted_pvalues.values()
        if p is not None and p < significance_threshold
    )

    _store_enrichment_results(
        adata=adata,
        result_key="ora_results",
        results=ora_df,
        method="ora",
        database=database,
        parameters={
            "deg_pvalue_threshold": pvalue_threshold,
            "significance_threshold": significance_threshold,
            "adjust_method": adjust_method,
            "min_size": min_size,
            "max_size": max_size,
            "n_query_genes": len(query_genes),
        },
        statistics={
            "n_gene_sets": len(gene_sets),
            "n_significant": n_significant,
            "n_query_genes": len(query_genes),
        },
        species=species,
        gene_sets=gene_sets,
        gene_sets_shared_key="enrichment_ora_gene_sets",
        data_id=data_id,
    )

    # Filter result dicts to significant pathways (reduces MCP response size)
    (
        filtered_statistics,
        filtered_scores,
        filtered_pvals,
        filtered_adj_pvals,
    ) = _filter_significant_statistics(
        gene_set_statistics,
        enrichment_scores,
        pvalues,
        adjusted_pvalues,
        method="ora",
        fdr_threshold=significance_threshold,
    )

    return EnrichmentResult(
        method="ora",
        n_gene_sets=len(gene_sets),
        n_significant=n_significant,
        enrichment_scores=filtered_scores,
        pvalues=filtered_pvals,
        adjusted_pvalues=filtered_adj_pvals,
        gene_set_statistics=filtered_statistics,
        top_gene_sets=top_gene_sets,
        top_depleted_sets=[],  # ORA does not produce depleted gene sets
    )


def perform_ssgsea(
    adata: "ad.AnnData",
    gene_sets: dict[str, list[str]],
    min_size: int = 10,
    max_size: int = 500,
    species: Optional[str] = None,
    database: Optional[str] = None,
    ctx: Optional["ToolContext"] = None,
    data_id: Optional[str] = None,
) -> "EnrichmentResult":
    """Perform single-sample Gene Set Enrichment Analysis (ssGSEA).

    Calculates enrichment scores for each sample independently.

    Args:
        adata: Annotated data matrix.
        gene_sets: Gene sets to test (name -> gene list).
        min_size: Minimum gene set size.
        max_size: Maximum gene set size.
        species: Species for the analysis ('mouse', 'human').
        database: Gene set database used ('KEGG_Pathways', 'GO_Biological_Process').
        ctx: MCP tool context for logging.
        data_id: Dataset identifier for result tracking.

    Returns:
        EnrichmentResult with per-sample enrichment scores.
    """
    # gseapy imported at module level (required dependency)

    # Memory-efficient batch processing for large datasets
    # Threshold: process in batches if > 1000 samples to avoid OOM
    BATCH_SIZE = 500
    n_samples = adata.n_obs

    # Run ssGSEA (with batch processing for large datasets)
    try:
        if n_samples <= BATCH_SIZE:
            # Small dataset: process all at once (original behavior)
            expr_df = pd.DataFrame(
                to_dense(adata.X).T, index=adata.var_names, columns=adata.obs_names
            )
            res = _get_gseapy().ssgsea(
                data=expr_df,
                gene_sets=gene_sets,
                min_size=min_size,
                max_size=max_size,
                permutation_num=0,
                no_plot=True,
                threads=1,
                seed=42,
            )
        else:
            # Large dataset: batch processing to reduce peak memory
            # Memory reduction: O(n_genes × n_samples) -> O(n_genes × batch_size)
            all_batch_results = []

            for batch_start in range(0, n_samples, BATCH_SIZE):
                batch_end = min(batch_start + BATCH_SIZE, n_samples)
                batch_indices = list(range(batch_start, batch_end))

                # Extract batch - only convert this batch to dense
                batch_X = to_dense(adata.X[batch_indices, :])
                batch_df = pd.DataFrame(
                    batch_X.T,
                    index=adata.var_names,
                    columns=adata.obs_names[batch_indices],
                )

                batch_res = _get_gseapy().ssgsea(
                    data=batch_df,
                    gene_sets=gene_sets,
                    min_size=min_size,
                    max_size=max_size,
                    permutation_num=0,
                    no_plot=True,
                    threads=1,
                    seed=42,
                )

                if hasattr(batch_res, "results"):
                    all_batch_results.append(batch_res.results)

                # Free batch memory
                del batch_X, batch_df

            # Merge batch results into unified format
            # Create a mock result object with combined results
            class CombinedResult:
                def __init__(self, results_list):
                    self.results = {}
                    for batch_results in results_list:
                        if isinstance(batch_results, dict):
                            self.results.update(batch_results)

            res = CombinedResult(all_batch_results)

        if hasattr(res, "res2d") and isinstance(res.res2d, pd.DataFrame):
            res2d = res.res2d
            if {"Name", "Term", "ES"}.issubset(res2d.columns):
                scores_df = res2d.pivot(index="Term", columns="Name", values="ES")
            else:
                scores_df = pd.DataFrame()
        elif hasattr(res, "results") and isinstance(res.results, dict):
            all_samples = list(res.results.keys())
            all_gene_sets = set()

            for sample_result in res.results.values():
                if isinstance(sample_result, pd.DataFrame) and "Term" in sample_result.columns:
                    all_gene_sets.update(sample_result["Term"].values)
                elif isinstance(sample_result, dict):
                    all_gene_sets.update(sample_result)

            all_gene_sets_list = list(all_gene_sets)
            scores_matrix = pd.DataFrame(
                index=all_gene_sets_list, columns=all_samples, dtype=float
            )

            for sample, sample_result in res.results.items():
                if (
                    isinstance(sample_result, pd.DataFrame)
                    and "Term" in sample_result.columns
                    and "ES" in sample_result.columns
                ):
                    sample_scores = sample_result.set_index("Term")["ES"]
                    scores_matrix[sample] = sample_scores.reindex(all_gene_sets_list)
                elif isinstance(sample_result, dict):
                    sample_scores = {
                        term: values.get("es")
                        for term, values in sample_result.items()
                        if isinstance(values, dict) and "es" in values
                    }
                    scores_matrix[sample] = pd.Series(sample_scores).reindex(
                        all_gene_sets_list
                    )

            scores_df = scores_matrix
        else:
            error_msg = "ssGSEA results format not recognized."
            logger.error(error_msg)
            raise ProcessingError(error_msg)

        # Calculate statistics - vectorized, NaN-aware (50x faster than row-by-row)
        enrichment_scores = {}
        gene_set_statistics = {}

        if not scores_df.empty:
            values = scores_df.values
            means = np.nanmean(values, axis=1)
            stds = np.nanstd(values, axis=1)
            mins = np.nanmin(values, axis=1)
            maxs = np.nanmax(values, axis=1)

            enrichment_scores = dict(
                zip(scores_df.index, means.astype(float), strict=True)
            )

            for i, gs_name in enumerate(scores_df.index):
                gene_set_statistics[gs_name] = {
                    "mean_score": float(means[i]),
                    "std_score": float(stds[i]),
                    "min_score": float(mins[i]),
                    "max_score": float(maxs[i]),
                    "size": len(gene_sets.get(gs_name, [])),
                }

            # Store per-cell scores in adata.obs, preserving NaN to
            # distinguish "not estimable" from "zero enrichment".
            scores_T = scores_df.T
            for gs_name in scores_df.index:
                adata.obs[f"ssgsea_{gs_name}"] = scores_T[gs_name].values

            # Per-run parametrized key so multiple database runs coexist
            analysis_key = _build_enrichment_key("ssgsea", database)

            # Store gene set membership (shared for compat + per-run for provenance)
            adata.uns["enrichment_ssgsea_gene_sets"] = gene_sets
            per_run_gs_key = f"enrichment_ssgsea_gene_sets_{analysis_key.removeprefix('enrichment_')}"
            adata.uns[per_run_gs_key] = gene_sets

            adata.uns["enrichment_latest_score_key"] = analysis_key

            # Store metadata for scientific provenance tracking
            obs_keys = [f"ssgsea_{gs_name}" for gs_name in scores_df.index]
            store_analysis_metadata(
                adata,
                analysis_name=analysis_key,
                method="ssgsea",
                parameters={
                    "min_size": min_size,
                    "max_size": max_size,
                },
                results_keys={
                    "obs": obs_keys,
                    "uns": [per_run_gs_key],
                },
                statistics={
                    "n_gene_sets": len(gene_sets),
                    "n_samples": adata.n_obs,
                },
                species=species,
                database=database,
            )

            # Export results for reproducibility
            if data_id is not None:
                export_analysis_result(adata, data_id, analysis_key)

        # Get top gene sets by mean enrichment
        sorted_by_mean = sorted(
            enrichment_scores.items(), key=lambda x: x[1], reverse=True
        )
        top_gene_sets = [x[0] for x in sorted_by_mean[:10]]

        # ssGSEA doesn't provide p-values, so return empty gene_set_statistics
        # to reduce MCP response size (no significance filtering possible)
        pvalues = None
        adjusted_pvalues = None

        return EnrichmentResult(
            method="ssgsea",
            n_gene_sets=len(gene_sets),
            # IMPORTANT: ssGSEA does NOT perform significance testing
            # Setting n_significant=0 is honest: no pathways are "statistically significant"
            # All gene sets receive enrichment scores, but these are sample-level metrics
            # without associated p-values. Use GSEA or ORA for significance testing.
            n_significant=0,  # ssGSEA doesn't test significance - no p-values produced
            enrichment_scores=enrichment_scores,  # Mean scores per gene set
            pvalues=pvalues,
            adjusted_pvalues=adjusted_pvalues,
            gene_set_statistics={},  # Empty to reduce response size (no p-values available)
            top_gene_sets=top_gene_sets,
            top_depleted_sets=[],  # ssGSEA doesn't produce depleted sets
        )

    except Exception as e:
        logger.error(f"ssGSEA failed: {e}")
        raise


def perform_enrichr(
    gene_list: list[str],
    gene_sets: Optional[str] = None,
    organism: str = "human",
    pvalue_cutoff: float = 0.05,
    ctx: Optional["ToolContext"] = None,
    adata: Optional["ad.AnnData"] = None,
    species: Optional[str] = None,
    database: Optional[str] = None,
    data_id: Optional[str] = None,
) -> "EnrichmentResult":
    """Perform enrichment analysis using Enrichr web service.

    Args:
        gene_list: List of genes to analyze.
        gene_sets: Enrichr library name. If None, use default libraries.
        organism: Organism ('human' or 'mouse').
        pvalue_cutoff: Adjusted p-value threshold for significance.
        ctx: MCP tool context for logging.

    Returns:
        EnrichmentResult containing enrichment results.
    """
    # gseapy imported at module level (required dependency)

    # Default gene set libraries - use separate variable for list
    gene_sets_list: list[str]
    if gene_sets is None:
        gene_sets_list = [
            map_gene_set_database_to_enrichr_library(db_name, organism)
            for db_name in DEFAULT_ENRICHR_DATABASES
        ]
    else:
        # Map user-friendly database name to actual Enrichr library name
        enrichr_library = map_gene_set_database_to_enrichr_library(gene_sets, organism)
        gene_sets_list = [enrichr_library]

    # Run Enrichr
    try:
        enr = _get_gseapy().enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets_list,
            organism=organism.capitalize(),
            outdir=None,
            cutoff=pvalue_cutoff,
        )

        # Get results - enr.results is already a DataFrame
        all_results = enr.results

        # Prepare output - OPTIMIZED: vectorized dict + array iteration (12x faster)
        enrichment_scores = dict(
            zip(
                all_results["Term"],
                all_results["Combined Score"],
                strict=True,
            )
        )
        pvalues = dict(
            zip(all_results["Term"], all_results["P-value"], strict=True)
        )
        adjusted_pvalues = dict(
            zip(
                all_results["Term"],
                all_results["Adjusted P-value"],
                strict=True,
            )
        )

        # Pre-extract arrays for fast iteration
        terms = all_results["Term"].values
        combined_scores = all_results["Combined Score"].values
        p_values = all_results["P-value"].values
        adj_p_values = all_results["Adjusted P-value"].values
        z_scores = (
            all_results["Z-score"].values
            if "Z-score" in all_results.columns
            else np.full(len(terms), np.nan)
        )
        overlaps = all_results["Overlap"].values
        genes_strs = all_results["Genes"].values
        odds_ratios = (
            all_results["Odds Ratio"].values
            if "Odds Ratio" in all_results.columns
            else np.ones(len(terms))
        )

        # Pre-split all genes once
        genes_split = [
            genes_str.split(";") if isinstance(genes_str, str) else []
            for genes_str in genes_strs
        ]

        # Build gene_set_statistics with array indexing
        gene_set_statistics = {}
        for i in range(len(terms)):
            gene_set_statistics[terms[i]] = {
                "combined_score": float(combined_scores[i]),
                "pval": float(p_values[i]),
                "adjusted_pval": float(adj_p_values[i]),
                "z_score": float(z_scores[i]),
                "overlap": overlaps[i],
                "genes": genes_split[i],
                "odds_ratio": float(odds_ratios[i]),
            }

        # Get top results
        all_results_sorted = all_results.sort_values("Combined Score", ascending=False)
        top_gene_sets = all_results_sorted.head(10)["Term"].tolist()

        n_significant = int((all_results["Adjusted P-value"] < pvalue_cutoff).sum())
        if adata is not None:
            _store_enrichment_results(
                adata=adata,
                result_key="enrichr_results",
                results=all_results,
                method="enrichr",
                database=database,
                parameters={
                    "pvalue_cutoff": pvalue_cutoff,
                    "n_query_genes": len(gene_list),
                },
                statistics={
                    "n_gene_sets": len(all_results),
                    "n_significant": n_significant,
                },
                species=species or organism,
                data_id=data_id,
            )

        # Filter all result dictionaries to only significant pathways (reduces MCP response size)
        # Uses method-based FDR threshold: Enrichr = 0.05 (same as ORA, hypergeometric-based)
        (
            filtered_statistics,
            filtered_scores,
            filtered_pvals,
            filtered_adj_pvals,
        ) = _filter_significant_statistics(
            gene_set_statistics,
            enrichment_scores,
            pvalues,
            adjusted_pvalues,
            method="enrichr",
            fdr_threshold=pvalue_cutoff,
        )

        return EnrichmentResult(
            method="enrichr",
            n_gene_sets=len(all_results),
            n_significant=n_significant,
            enrichment_scores=filtered_scores,
            pvalues=filtered_pvals,
            adjusted_pvalues=filtered_adj_pvals,
            gene_set_statistics=filtered_statistics,
            top_gene_sets=top_gene_sets,
            top_depleted_sets=[],  # Enrichr doesn't produce depleted sets
        )

    except Exception as e:
        logger.error(f"Enrichr failed: {e}")
        raise


# ============================================================================
# Spatial Enrichment Analysis Functions (EnrichMap-based)
# ============================================================================


async def perform_spatial_enrichment(
    data_id: str,
    ctx: "ToolContext",
    gene_sets: Union[list[str], dict[str, list[str]]],
    score_keys: Optional[Union[str, list[str]]] = None,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    smoothing: bool = True,
    correct_spatial_covariates: bool = True,
    batch_key: Optional[str] = None,
    species: str = "unknown",
    database: Optional[str] = None,
) -> "EnrichmentResult":
    """Perform spatially-aware gene set enrichment analysis using EnrichMap.

    Args:
        data_id: Identifier for the spatial data in the data store.
        ctx: MCP tool context for data access and logging.
        gene_sets: Either a single gene list or a dictionary of gene sets where
            keys are signature names and values are lists of genes.
        score_keys: Names for the gene signatures if gene_sets is a list.
            Ignored if gene_sets is already a dictionary.
        spatial_key: Key in adata.obsm containing spatial coordinates.
        n_neighbors: Number of nearest spatial neighbors for smoothing.
        smoothing: Whether to perform spatial smoothing.
        correct_spatial_covariates: Whether to correct for spatial covariates using GAM.
        batch_key: Column in adata.obs for batch-wise normalization.
        species: Species for the analysis (e.g., 'mouse', 'human').
        database: Gene set database used (e.g., 'KEGG_Pathways', 'GO_Biological_Process').

    Returns:
        EnrichmentResult containing enrichment scores and statistics.
    """
    em = require("enrichmap", ctx, feature="spatial enrichment analysis")

    # Get data using standard ctx pattern
    adata = await ctx.get_adata(data_id)

    # Validate spatial coordinates
    if spatial_key not in adata.obsm:
        raise ProcessingError(
            f"Spatial coordinates '{spatial_key}' not found in adata.obsm"
        )

    # Convert single gene list to dictionary format
    gene_sets_dict: dict[str, list[str]]
    if isinstance(gene_sets, list):
        # For a single gene list, score_keys should be a string name
        if score_keys is None:
            sig_name = "enrichmap_signature"
        elif isinstance(score_keys, str):
            sig_name = score_keys
        else:
            # If score_keys is a list, use the first element
            sig_name = score_keys[0] if score_keys else "enrichmap_signature"
        gene_sets_dict = {sig_name: gene_sets}
    else:
        gene_sets_dict = gene_sets

    # Validate gene sets with format conversion
    available_genes = set(adata.var_names)
    validated_gene_sets = {}

    for sig_name, genes in gene_sets_dict.items():
        common_genes = _match_gene_set_to_dataset(genes, available_genes, species)

        if len(common_genes) < 2:
            await ctx.warning(
                f"Signature '{sig_name}' has {len(common_genes)} genes in the dataset. Skipping."
            )
            continue
        validated_gene_sets[sig_name] = common_genes
        await ctx.info(
            f"Signature '{sig_name}': {len(common_genes)}/{len(genes)} genes found"
        )

    if not validated_gene_sets:
        raise ProcessingError(
            f"No valid gene signatures found (≥2 genes). "
            f"Dataset: {len(available_genes)} genes, requested: {len(gene_sets_dict)} signatures. "
            f"Check species (human/mouse) and gene name format."
        )

    # Run EnrichMap scoring - process each gene set individually
    failed_signatures = []
    successful_signatures = []

    for sig_name, genes in validated_gene_sets.items():
        try:
            em.tl.score(
                adata=adata,
                gene_set=genes,  # Fixed: use gene_set (correct API parameter name)
                score_key=sig_name,  # Fixed: provide explicit score_key
                spatial_key=spatial_key,
                n_neighbors=n_neighbors,
                smoothing=smoothing,
                correct_spatial_covariates=correct_spatial_covariates,
                batch_key=batch_key,
            )
            successful_signatures.append(sig_name)

        except Exception as e:
            await ctx.warning(f"EnrichMap failed for '{sig_name}': {e}")
            failed_signatures.append((sig_name, str(e)))

    # Check if any signatures were processed successfully
    if not successful_signatures:
        error_details = "; ".join(
            [f"{name}: {error}" for name, error in failed_signatures]
        )
        raise ProcessingError(
            f"All EnrichMap scoring failed. This may indicate:\n"
            f"1. EnrichMap package installation issues\n"
            f"2. Incompatible gene names or data format\n"
            f"3. Insufficient spatial information\n"
            f"Details: {error_details}"
        )

    # Update validated_gene_sets to only include successful ones
    validated_gene_sets = {
        sig: validated_gene_sets[sig] for sig in successful_signatures
    }

    if ctx and failed_signatures:
        await ctx.warning(
            f"Failed to process {len(failed_signatures)} gene sets: {[name for name, _ in failed_signatures]}"
        )

    # Collect results
    score_columns = [f"{sig}_score" for sig in validated_gene_sets]

    # Calculate summary statistics
    summary_stats = {}
    for sig_name in validated_gene_sets:
        score_col = f"{sig_name}_score"
        scores = adata.obs[score_col]

        summary_stats[sig_name] = {
            "mean": float(scores.mean()),
            "std": float(scores.std()),
            "min": float(scores.min()),
            "max": float(scores.max()),
            "median": float(scores.median()),
            "q25": float(scores.quantile(0.25)),
            "q75": float(scores.quantile(0.75)),
            "n_genes": len(validated_gene_sets[sig_name]),
        }

    # Per-run parametrized key so multiple database runs coexist
    analysis_key = _build_enrichment_key("spatial", database)

    # Store gene set membership (shared for compat + per-run for provenance)
    adata.uns["enrichment_spatial_gene_sets"] = validated_gene_sets
    per_run_gs_key = (
        f"enrichment_spatial_gene_sets_{analysis_key.removeprefix('enrichment_')}"
    )
    adata.uns[per_run_gs_key] = validated_gene_sets

    # Store metadata for scientific provenance tracking
    store_analysis_metadata(
        adata,
        analysis_name=analysis_key,
        method="spatial_enrichmap",
        parameters={
            "spatial_key": spatial_key,
            "n_neighbors": n_neighbors,
            "smoothing": smoothing,
            "correct_spatial_covariates": correct_spatial_covariates,
            "batch_key": batch_key,
        },
        results_keys={
            "obs": score_columns,
            "uns": [per_run_gs_key],
        },
        statistics={
            "n_gene_sets": len(validated_gene_sets),
            "n_successful_signatures": len(successful_signatures),
            "n_failed_signatures": len(failed_signatures),
        },
        species=species,
        database=database,
    )

    # Export results for reproducibility
    export_analysis_result(adata, data_id, analysis_key)

    # Create enrichment scores (use max score per gene set)
    enrichment_scores = {
        sig_name: float(stats["max"]) for sig_name, stats in summary_stats.items()
    }

    # Sort by enrichment score to get top gene sets
    sorted_sigs = sorted(enrichment_scores.items(), key=lambda x: x[1], reverse=True)
    top_gene_sets = [sig_name for sig_name, _ in sorted_sigs[:10]]

    # Spatial enrichment doesn't provide p-values, so return empty gene_set_statistics
    # to reduce MCP response size (no significance filtering possible)
    pvalues = None
    adjusted_pvalues = None

    return EnrichmentResult(
        method="spatial_enrichmap",
        n_gene_sets=len(validated_gene_sets),
        # No significance testing in spatial enrichment — n_significant=0
        # (successful computation != statistical significance)
        n_significant=0,
        n_successful_signatures=len(successful_signatures),
        enrichment_scores=enrichment_scores,
        pvalues=pvalues,
        adjusted_pvalues=adjusted_pvalues,
        gene_set_statistics={},  # Empty to reduce response size (no p-values available)
        spatial_scores_key=None,  # Scores are in obs columns, not obsm
        top_gene_sets=top_gene_sets,
        top_depleted_sets=[],  # Spatial enrichment doesn't produce depleted sets
    )


# ============================================================================
# Gene Set Loading Functions
# ============================================================================
# Simplified from GeneSetLoader class - no need for class overhead when
# functions are only called once from load_gene_sets()


def _get_organism_name(species: str) -> str:
    """Get organism name for gseapy from species code."""
    return "Homo sapiens" if species.lower() == "human" else "Mus musculus"


def load_msigdb_gene_sets(
    species: str,
    collection: str = "H",
    subcollection: Optional[str] = None,
    min_size: int = 10,
    max_size: int = 500,
) -> dict[str, list[str]]:
    """Load gene sets from MSigDB using gseapy.

    Args:
        species: Species for gene sets ('human' or 'mouse').
        collection: MSigDB collection name. Options include:
            H (hallmark gene sets), C1 (positional gene sets),
            C2 (curated gene sets, e.g., CGP, CP:KEGG, CP:REACTOME),
            C3 (motif gene sets), C4 (computational gene sets),
            C5 (GO gene sets: CC, BP, MF), C6 (oncogenic signatures),
            C7 (immunologic signatures), C8 (cell type signatures).
        subcollection: Subcollection for specific databases (e.g., 'CP:KEGG', 'GO:BP').
        min_size: Minimum gene set size.
        max_size: Maximum gene set size.

    Returns:
        Dictionary of gene sets.
    """
    # gseapy imported at module level (required dependency)
    try:
        organism = _get_organism_name(species)
        gene_sets_dict = {}

        if collection == "H":
            # Hallmark gene sets
            gene_sets_dict = _load_library_first_available("MSigDB_Hallmark", organism)

        elif collection == "C2" and subcollection == "CP:KEGG":
            # KEGG pathways
            if species.lower() == "human":
                gene_sets_dict = _load_enrichr_library("KEGG_2021_Human", organism)
            else:
                gene_sets_dict = _load_enrichr_library("KEGG_2019_Mouse", organism)

        elif collection == "C2" and subcollection == "CP:REACTOME":
            # Reactome pathways
            gene_sets_dict = _load_library_first_available(
                "Reactome_Pathways", organism
            )

        elif collection == "C5":
            # GO gene sets
            if subcollection == "GO:BP" or subcollection is None:
                gene_sets_dict.update(
                    _load_library_first_available("GO_Biological_Process", organism)
                )
            if subcollection == "GO:MF" or subcollection is None:
                gene_sets_dict.update(
                    _load_library_first_available("GO_Molecular_Function", organism)
                )
            if subcollection == "GO:CC" or subcollection is None:
                gene_sets_dict.update(
                    _load_library_first_available("GO_Cellular_Component", organism)
                )

        elif collection == "C8":
            # Cell type signatures
            gene_sets_dict = _load_library_first_available(
                "Cell_Type_Markers", organism
            )

        # Filter by size
        filtered_sets = _filter_gene_sets_by_size(gene_sets_dict, min_size, max_size)
        return filtered_sets

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to load MSigDB gene sets: {e}") from e


def load_go_gene_sets(
    species: str,
    aspect: str = "BP",
    min_size: int = 10,
    max_size: int = 500,
) -> dict[str, list[str]]:
    """Load GO terms using gseapy.

    Args:
        species: Species for gene sets ('human' or 'mouse').
        aspect: GO aspect. Options: 'BP' (biological process),
            'MF' (molecular function), 'CC' (cellular component).
        min_size: Minimum gene set size.
        max_size: Maximum gene set size.

    Returns:
        Dictionary of GO gene sets.
    """
    if aspect not in GO_ASPECT_TO_DATABASE_KEY:
        raise ParameterError(f"Invalid GO aspect: {aspect}")

    # gseapy imported at module level (required dependency)
    try:
        organism = _get_organism_name(species)
        database_key = GO_ASPECT_TO_DATABASE_KEY[aspect]
        gene_sets = _load_library_first_available(database_key, organism)

        # Filter by size
        filtered_sets = _filter_gene_sets_by_size(gene_sets, min_size, max_size)
        return filtered_sets

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to load GO gene sets: {e}") from e


def load_kegg_gene_sets(
    species: str, min_size: int = 10, max_size: int = 500
) -> dict[str, list[str]]:
    """Load KEGG pathways using gseapy.

    Args:
        species: Species for gene sets ('human' or 'mouse').
        min_size: Minimum gene set size.
        max_size: Maximum gene set size.

    Returns:
        Dictionary of KEGG pathway gene sets.
    """
    # gseapy imported at module level (required dependency)
    try:
        organism = _get_organism_name(species)

        species_lower = species.lower()
        if species_lower == "human":
            gene_sets = _load_enrichr_library("KEGG_2021_Human", organism)
        elif species_lower == "mouse":
            gene_sets = _load_enrichr_library("KEGG_2019_Mouse", organism)
        else:
            raise ParameterError(
                f"KEGG_Pathways is only available for 'human' and 'mouse', "
                f"not '{species}'. Use GO or Reactome databases instead."
            )

        # Filter by size
        filtered_sets = _filter_gene_sets_by_size(gene_sets, min_size, max_size)
        return filtered_sets

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to load KEGG pathways: {e}") from e


def load_reactome_gene_sets(
    species: str, min_size: int = 10, max_size: int = 500
) -> dict[str, list[str]]:
    """Load Reactome pathways using gseapy.

    Args:
        species: Species for gene sets ('human' or 'mouse').
        min_size: Minimum gene set size.
        max_size: Maximum gene set size.

    Returns:
        Dictionary of Reactome pathway gene sets.
    """
    # gseapy imported at module level (required dependency)
    try:
        organism = _get_organism_name(species)
        gene_sets = _load_library_first_available("Reactome_Pathways", organism)

        # Filter by size (use shared utility for consistency)
        filtered_sets = _filter_gene_sets_by_size(gene_sets, min_size, max_size)
        return filtered_sets

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to load Reactome pathways: {e}") from e


def load_cell_marker_gene_sets(
    species: str, min_size: int = 5, max_size: int = 200
) -> dict[str, list[str]]:
    """Load cell type marker gene sets using gseapy.

    Args:
        species: Species for gene sets ('human' or 'mouse').
        min_size: Minimum gene set size.
        max_size: Maximum gene set size.

    Returns:
        Dictionary of cell type marker gene sets.
    """
    # gseapy imported at module level (required dependency)
    try:
        organism = _get_organism_name(species)
        gene_sets = _load_enrichr_library("CellMarker_Augmented_2021", organism)

        # Filter by size
        filtered_sets = _filter_gene_sets_by_size(gene_sets, min_size, max_size)
        return filtered_sets

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to load cell markers: {e}") from e


def load_gene_sets(
    database: str,
    species: str = "human",
    min_genes: int = 10,
    max_genes: int = 500,
    ctx: Optional["ToolContext"] = None,
) -> dict[str, list[str]]:
    """Load gene sets from specified database.

    Args:
        database: Database name. Options: GO_Biological_Process, GO_Molecular_Function,
            GO_Cellular_Component, KEGG_Pathways, Reactome_Pathways, MSigDB_Hallmark,
            Cell_Type_Markers.
        species: Species ('human' or 'mouse').
        min_genes: Minimum gene set size.
        max_genes: Maximum gene set size.
        ctx: MCP tool context for logging.

    Returns:
        Dictionary of gene sets.
    """
    # Direct function calls - no class overhead
    database_map = {
        "GO_Biological_Process": lambda: load_go_gene_sets(
            species, "BP", min_genes, max_genes
        ),
        "GO_Molecular_Function": lambda: load_go_gene_sets(
            species, "MF", min_genes, max_genes
        ),
        "GO_Cellular_Component": lambda: load_go_gene_sets(
            species, "CC", min_genes, max_genes
        ),
        "KEGG_Pathways": lambda: load_kegg_gene_sets(species, min_genes, max_genes),
        "Reactome_Pathways": lambda: load_reactome_gene_sets(
            species, min_genes, max_genes
        ),
        "MSigDB_Hallmark": lambda: load_msigdb_gene_sets(
            species, "H", None, min_genes, max_genes
        ),
        "Cell_Type_Markers": lambda: load_cell_marker_gene_sets(
            species, min_genes, max_genes
        ),
    }

    if database not in database_map:
        raise ParameterError(
            f"Unknown database: {database}. Available: {list(database_map)}"
        )

    gene_sets = database_map[database]()
    return gene_sets


# ============================================================================
# UNIFIED ENRICHMENT ANALYSIS ENTRY POINT
# ============================================================================


async def analyze_enrichment(
    data_id: str,
    ctx: "ToolContext",
    params: "EnrichmentParameters",
) -> EnrichmentResult:
    """
    Unified entry point for gene set enrichment analysis.

    This function handles all enrichment methods with a consistent interface:
    - Gene set loading from databases
    - Method dispatch (GSEA, ORA, ssGSEA, Enrichr, spatial)
    - Error handling with clear messages

    Args:
        data_id: Dataset ID
        ctx: ToolContext for data access and logging
        params: EnrichmentParameters with method, species, database, etc.

    Returns:
        EnrichmentResult with enrichment scores and statistics

    Raises:
        ParameterError: If params is None or invalid
        ProcessingError: If gene set loading or analysis fails
    """
    # Import here to avoid circular imports
    from ..utils.adata_utils import get_highly_variable_genes

    # Validate params
    if params is None:
        raise ParameterError(
            "params parameter is required for enrichment analysis.\n"
            "You must provide EnrichmentParameters with at least 'species' specified.\n"
            "Example: params={'species': 'mouse', 'method': 'pathway_ora'}"
        )

    # Get adata
    adata = await ctx.get_adata(data_id)

    ranking_key: str | None = None
    if params.score_keys is not None:
        ranking_key = (
            params.score_keys[0]
            if isinstance(params.score_keys, list)
            else params.score_keys
        )
        if params.method == "pathway_gsea":
            _validate_gsea_ranking_key(adata, ranking_key)

    # Load gene sets
    loaded_from_database = False
    gene_sets = params.gene_sets
    if gene_sets is None and params.gene_set_database:
        try:
            gene_sets = load_gene_sets(
                database=params.gene_set_database,
                species=params.species,
                min_genes=params.min_genes,
                max_genes=params.max_genes,
                ctx=ctx,
            )
            loaded_from_database = True
        except ChatSpatialError:
            raise
        except Exception as e:
            await ctx.error(f"Gene set database loading failed: {e}")
            raise ProcessingError(
                f"Failed to load gene sets from {params.gene_set_database}: {e}\n\n"
                f"SOLUTIONS:\n"
                f"1. Check your internet connection\n"
                f"2. Verify species parameter: '{params.species}'\n"
                f"3. Try a different database (KEGG_Pathways, GO_Biological_Process)\n"
                f"4. Provide custom gene sets via 'gene_sets' parameter"
            ) from e

    # Validate gene sets
    if gene_sets is None or len(gene_sets) == 0:
        raise ProcessingError(
            "No valid gene sets available. "
            "Please provide gene sets via 'gene_sets' parameter or "
            "specify a valid 'gene_set_database'."
        )

    # Normalize gene_sets to dict format (convert list to single gene set dict)
    gene_sets_dict: dict[str, list[str]]
    if isinstance(gene_sets, list):
        gene_sets_dict = {"user_genes": gene_sets}
    else:
        gene_sets_dict = gene_sets

    if params.method == "spatial_enrichmap" and loaded_from_database:
        n_loaded_gene_sets = len(gene_sets_dict)
        gene_sets_dict = _limit_spatial_enrichmap_gene_sets(
            gene_sets_dict,
            available_genes=set(adata.var_names),
            species=params.species,
        )
        if len(gene_sets_dict) < n_loaded_gene_sets:
            await ctx.warning(
                f"Spatial EnrichMap scored the {len(gene_sets_dict)} database gene sets "
                f"with highest dataset overlap out of {n_loaded_gene_sets} loaded sets. "
                "Provide custom gene_sets to score specific signatures."
            )
        gene_sets = gene_sets_dict

    # Dispatch to appropriate method
    if params.method == "spatial_enrichmap":
        result = await perform_spatial_enrichment(
            data_id=data_id,
            ctx=ctx,
            gene_sets=gene_sets,
            score_keys=params.score_keys,
            spatial_key=params.spatial_key,
            n_neighbors=params.n_neighbors,
            smoothing=params.smoothing,
            correct_spatial_covariates=params.correct_spatial_covariates,
            batch_key=params.batch_key,
            species=params.species,
            database=params.gene_set_database,
        )

    elif params.method == "pathway_gsea":
        result = perform_gsea(
            adata=adata,
            gene_sets=gene_sets_dict,
            ranking_key=ranking_key,
            permutation_num=params.n_permutations,
            min_size=params.min_genes,
            max_size=params.max_genes,
            pvalue_cutoff=params.pvalue_cutoff,
            species=params.species,
            database=params.gene_set_database,
            ctx=ctx,
            data_id=data_id,
        )

    elif params.method == "pathway_ora":
        result = perform_ora(
            adata=adata,
            gene_sets=gene_sets_dict,
            pvalue_threshold=params.pvalue_cutoff,
            min_size=params.min_genes,
            max_size=params.max_genes,
            adjust_method=params.adjust_method,
            species=params.species,
            database=params.gene_set_database,
            ctx=ctx,
            data_id=data_id,
        )

    elif params.method == "pathway_ssgsea":
        result = perform_ssgsea(
            adata=adata,
            gene_sets=gene_sets_dict,
            min_size=params.min_genes,
            max_size=params.max_genes,
            species=params.species,
            database=params.gene_set_database,
            ctx=ctx,
            data_id=data_id,
        )

    elif params.method == "pathway_enrichr":
        gene_list = get_highly_variable_genes(adata, max_genes=500)
        result = perform_enrichr(
            gene_list=gene_list,
            gene_sets=params.gene_set_database,
            organism=params.species,
            pvalue_cutoff=params.pvalue_cutoff,
            ctx=ctx,
            adata=adata,
            species=params.species,
            database=params.gene_set_database,
            data_id=data_id,
        )

    else:
        raise ParameterError(f"Unknown enrichment method: {params.method}")

    await ctx.set_adata(data_id, adata)
    return result
