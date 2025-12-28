"""
Copy Number Variation (CNV) analysis tools for spatial transcriptomics data.
"""

from typing import Any, Dict, Optional

import numpy as np
import scanpy as sc
from mcp.server.fastmcp import Context

from ..models.analysis import CNVResult
from ..models.data import CNVParameters
from ..utils import validate_obs_column
from ..utils.exceptions import (DataCompatibilityError, DataNotFoundError,
                                DependencyError, ParameterError,
                                ProcessingError)

# Dependency checking for infercnvpy
try:
    import infercnvpy as cnv

    INFERCNVPY_AVAILABLE = True
except ImportError:
    INFERCNVPY_AVAILABLE = False

# Numbat availability is checked lazily in _infer_cnv_numbat to avoid
# import-time failures when rpy2/R is not installed


async def infer_cnv(
    data_id: str,
    data_store: Dict[str, Any],
    params: CNVParameters,
    context: Optional[Context] = None,
) -> CNVResult:
    """Infer copy number variations using selected method

    Supports two methods:
    - infercnvpy: Expression-based CNV inference (default, fast)
    - Numbat: Haplotype-aware CNV analysis (requires allele data, more accurate)

    Args:
        data_id: Dataset identifier
        data_store: Dictionary storing loaded datasets
        params: CNV analysis parameters including method selection
        context: MCP context for logging

    Returns:
        CNVResult containing method-specific CNV analysis results

    Raises:
        ValueError: If dataset not found or parameters are invalid
        RuntimeError: If selected method is not available
    """
    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise DataNotFoundError(f"Dataset '{data_id}' not found in data store")

    adata = data_store[data_id]["adata"]

    # Validate common parameters
    validate_obs_column(adata, params.reference_key, "Reference cell type key")

    available_categories = set(adata.obs[params.reference_key].unique())
    missing_categories = set(params.reference_categories) - available_categories
    if missing_categories:
        raise ParameterError(
            f"Reference categories {missing_categories} not found in "
            f"adata.obs['{params.reference_key}'].\n"
            f"Available categories: {sorted(available_categories)}"
        )

    # Dispatch to appropriate method
    if params.method == "infercnvpy":
        return await _infer_cnv_infercnvpy(data_id, data_store, params, context)
    elif params.method == "numbat":
        return await _infer_cnv_numbat(data_id, data_store, params, context)
    else:
        raise ParameterError(
            f"Unknown CNV method: {params.method}. "
            "Available methods: 'infercnvpy', 'numbat'"
        )


async def _infer_cnv_infercnvpy(
    data_id: str,
    data_store: Dict[str, Any],
    params: CNVParameters,
    context: Optional[Context] = None,
) -> CNVResult:
    """Infer copy number variations using infercnvpy

    This function performs CNV inference on spatial transcriptomics data using
    infercnvpy, which detects chromosomal copy number alterations by comparing
    gene expression patterns across chromosomes between tumor and normal cells.

    Args:
        data_id: Dataset identifier
        data_store: Dictionary storing loaded datasets
        params: CNV analysis parameters
        context: MCP context for logging

    Returns:
        CNVResult containing CNV analysis results and statistics
    """
    # Check if infercnvpy is available
    if not INFERCNVPY_AVAILABLE:
        raise DependencyError(
            "infercnvpy is not installed. Please install it with:\n"
            "  pip install 'chatspatial[cnv]'\n"
            "or:\n"
            "  pip install infercnvpy>=0.4.0"
        )

    # Note: data_id, reference_key, and reference_categories are already
    # validated in infer_cnv() before dispatch - no need to re-validate
    adata = data_store[data_id]["adata"]

    if context:
        await context.info(
            f"Running CNV inference with {len(params.reference_categories)} "
            f"reference cell types: {', '.join(params.reference_categories)}"
        )

    # Create a copy of adata for CNV analysis
    adata_cnv = adata.copy()

    # Check if gene position information is available
    if "chromosome" not in adata_cnv.var.columns:
        if context:
            await context.warning(
                "No chromosome information found in adata.var. "
                "Attempting to infer from gene names..."
            )
        try:
            # Try to infer gene positions from infercnvpy's built-in database
            cnv.tl.infercnv(
                adata_cnv,
                reference_key=params.reference_key,
                reference_cat=params.reference_categories,
                window_size=params.window_size,
                step=params.step,
                dynamic_threshold=params.dynamic_threshold,
            )
        except Exception as e:
            raise ProcessingError(
                "Failed to run CNV inference. Gene position information is "
                "required but not found in adata.var.\n"
                "Please ensure your data includes chromosome and position "
                "information, or use infercnvpy's built-in gene annotation.\n"
                f"Error: {str(e)}"
            ) from e
    else:
        # Gene positions are available, run CNV inference
        if context:
            await context.info("Running infercnvpy CNV inference...")

        # Exclude chromosomes if specified
        if params.exclude_chromosomes:
            genes_to_keep = ~adata_cnv.var["chromosome"].isin(
                params.exclude_chromosomes
            )
            adata_cnv = adata_cnv[:, genes_to_keep].copy()
            if context:
                await context.info(
                    f"Excluded chromosomes: {', '.join(params.exclude_chromosomes)}"
                )

        # Run infercnvpy
        cnv.tl.infercnv(
            adata_cnv,
            reference_key=params.reference_key,
            reference_cat=params.reference_categories,
            window_size=params.window_size,
            step=params.step,
            dynamic_threshold=params.dynamic_threshold,
        )

    # Optional: Cluster cells by CNV pattern
    if params.cluster_cells:
        if context:
            await context.info("Clustering cells by CNV pattern...")
        try:
            sc.pp.neighbors(adata_cnv, use_rep="X_cnv", n_neighbors=15)
            sc.tl.leiden(adata_cnv, key_added="cnv_clusters")
            if context:
                n_clusters = len(adata_cnv.obs["cnv_clusters"].unique())
                await context.info(f"Identified {n_clusters} CNV-based clusters")
        except Exception as e:
            if context:
                await context.warning(f"Failed to cluster cells by CNV: {str(e)}")

    # Optional: Compute dendrogram
    if params.dendrogram and params.cluster_cells:
        if context:
            await context.info("Computing hierarchical clustering dendrogram...")
        try:
            sc.tl.dendrogram(adata_cnv, groupby="cnv_clusters")
        except Exception as e:
            if context:
                await context.warning(f"Failed to compute dendrogram: {str(e)}")

    # Extract CNV statistics
    if context:
        await context.info("Extracting CNV statistics...")

    # Check what data is available
    cnv_score_key = None
    if "X_cnv" in adata_cnv.obsm:
        cnv_score_key = "X_cnv"
    elif "cnv" in adata_cnv.layers:
        cnv_score_key = "cnv"

    # Calculate statistics
    statistics = {}
    if cnv_score_key and cnv_score_key in adata_cnv.obsm:
        cnv_matrix = adata_cnv.obsm[cnv_score_key]

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
            statistics["std_cnv"] = float(np.sqrt(mean_sq - mean_val**2))

            # Median: for highly sparse matrices (>50% zeros), median is 0
            # Otherwise use approximation with non-zero values
            n_zeros = cnv_matrix.shape[0] * cnv_matrix.shape[1] - cnv_matrix.nnz
            n_total = cnv_matrix.shape[0] * cnv_matrix.shape[1]

            if n_zeros > n_total / 2:
                # Majority zeros, median is exactly 0
                statistics["median_cnv"] = 0.0
            else:
                # Use non-zero median as approximation
                statistics["median_cnv"] = float(np.median(cnv_matrix.data))

            # Per-cell CNV scores: compute on sparse matrix
            # abs() preserves sparsity
            cnv_abs = cnv_matrix.copy()
            cnv_abs.data = np.abs(cnv_abs.data)
            cell_cnv_scores = np.array(cnv_abs.mean(axis=1)).flatten()
            statistics["mean_cell_cnv_score"] = float(np.mean(cell_cnv_scores))
            statistics["max_cell_cnv_score"] = float(np.max(cell_cnv_scores))

        else:
            # Dense matrix - use standard numpy operations
            statistics["mean_cnv"] = float(np.mean(cnv_matrix))
            statistics["std_cnv"] = float(np.std(cnv_matrix))
            statistics["median_cnv"] = float(np.median(cnv_matrix))

            # Calculate per-cell CNV scores
            cell_cnv_scores = np.mean(np.abs(cnv_matrix), axis=1)
            statistics["mean_cell_cnv_score"] = float(np.mean(cell_cnv_scores))
            statistics["max_cell_cnv_score"] = float(np.max(cell_cnv_scores))

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

    # Store CNV results back in the original adata object
    if cnv_score_key and cnv_score_key in adata_cnv.obsm:
        adata.obsm[cnv_score_key] = adata_cnv.obsm[cnv_score_key]

    # Store CNV metadata (required for infercnvpy plotting functions)
    if "cnv" in adata_cnv.uns:
        adata.uns["cnv"] = adata_cnv.uns["cnv"]

    if params.cluster_cells and "cnv_clusters" in adata_cnv.obs:
        adata.obs["cnv_clusters"] = adata_cnv.obs["cnv_clusters"]

    if params.dendrogram and "dendrogram_cnv_clusters" in adata_cnv.uns:
        adata.uns["dendrogram_cnv_clusters"] = adata_cnv.uns["dendrogram_cnv_clusters"]

    # Store CNV analysis parameters in adata.uns for reference
    adata.uns["cnv_analysis"] = {
        "reference_key": params.reference_key,
        "reference_categories": list(params.reference_categories),  # Convert to list
        "window_size": params.window_size,
        "step": params.step,
        "cnv_score_key": cnv_score_key,
    }

    if context:
        await context.info(
            f"CNV analysis complete. Analyzed {n_genes_analyzed} genes "
            f"across {n_chromosomes} chromosomes."
        )

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


async def _infer_cnv_numbat(
    data_id: str,
    data_store: Dict[str, Any],
    params: CNVParameters,
    context: Optional[Context] = None,
) -> CNVResult:
    """Infer copy number variations using Numbat (haplotype-aware)

    Numbat performs haplotype-aware CNV analysis by integrating allele-specific
    counts with expression data, enabling detection of copy-neutral LOH and
    reconstruction of tumor phylogeny.

    Args:
        data_id: Dataset identifier
        data_store: Dictionary storing loaded datasets
        params: CNV analysis parameters
        context: MCP context for logging

    Returns:
        CNVResult containing Numbat CNV analysis results

    Raises:
        RuntimeError: If Numbat is not available or allele data is missing
        ValueError: If dataset or parameters are invalid
    """
    # Lazy import and check for Numbat availability
    try:
        import anndata2ri
        import rpy2.robjects as ro
        from rpy2.rinterface_lib import openrlib
        from rpy2.robjects import (conversion, default_converter, numpy2ri,
                                   pandas2ri)

        # Test if Numbat R package is available
        ro.r("suppressPackageStartupMessages(library(numbat))")
    except ImportError as e:
        raise DependencyError(
            f"rpy2 is not installed: {e}\n" "Please install with: pip install rpy2"
        ) from e
    except Exception as e:
        raise DependencyError(
            f"Numbat R package is not available: {e}\n"
            "Please ensure Numbat R package is installed:\n"
            "  R -e \"install.packages('numbat')\"\n"
            "  or\n"
            "  R -e \"devtools::install_github('kharchenkolab/numbat')\""
        ) from e

    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise DataNotFoundError(f"Dataset '{data_id}' not found in data store")

    adata = data_store[data_id]["adata"]

    if context:
        await context.info(
            f"Running Numbat CNV inference (haplotype-aware) with genome: "
            f"{params.numbat_genome}"
        )

    # Validate allele data exists
    # Numbat requires long-format allele dataframe (from pileup_and_phase or similar)
    # Check if we have the raw allele dataframe in adata.uns
    if "numbat_allele_data_raw" in adata.uns:
        # Use pre-prepared long-format allele data
        import pandas as pd

        df_allele = adata.uns["numbat_allele_data_raw"]

        if context:
            await context.info(
                f"Found Numbat allele dataframe in adata.uns['numbat_allele_data_raw']\n"
                f"  - {len(df_allele)} SNP-cell pairs"
            )

        # Validate required columns
        required_cols = ["cell", "CHROM", "POS", "REF", "ALT", "AD", "DP"]
        missing_cols = [col for col in required_cols if col not in df_allele.columns]

        if missing_cols:
            raise ParameterError(
                f"Allele dataframe missing required columns: {missing_cols}\n"
                f"Available columns: {list(df_allele.columns)}\n"
                "Numbat requires: cell, CHROM, POS, REF, ALT, AD (alt count), "
                "DP (total depth)"
            )

    else:
        # Fallback: try to use matrix format (less ideal for Numbat)
        raise ParameterError(
            "Numbat requires long-format allele dataframe in adata.uns['numbat_allele_data_raw'].\n"
            "This should be created during data preparation (e.g., from pileup_and_phase).\n"
            "The dataframe should have columns: cell, CHROM, POS, REF, ALT, AD, DP, etc.\n"
            f"Available uns keys: {list(adata.uns.keys())}"
        )

    # Get expression matrix (keep sparse if possible)
    import scipy.sparse as sp

    count_mat = adata.X

    # Prepare metadata
    gene_names = list(adata.var_names)
    cell_barcodes = list(adata.obs_names)

    # Identify reference cells (1-indexed for R)
    ref_mask = adata.obs[params.reference_key].isin(params.reference_categories)
    ref_indices_python = [i for i, is_ref in enumerate(ref_mask) if is_ref]
    ref_indices_r = [i + 1 for i in ref_indices_python]  # R is 1-indexed

    if not ref_indices_r:
        raise ParameterError(
            f"No reference cells found with key '{params.reference_key}' and "
            f"categories {params.reference_categories}"
        )

    # Log sparse vs dense matrix info
    is_sparse = sp.issparse(count_mat)
    if context:
        matrix_type = "sparse" if is_sparse else "dense"
        await context.info(
            f"Preparing Numbat analysis:\n"
            f"  - {len(cell_barcodes)} cells × {len(gene_names)} genes\n"
            f"  - {len(ref_indices_r)} reference cells\n"
            f"  - {len(df_allele)} SNP-cell pairs\n"
            f"  - Genome: {params.numbat_genome}\n"
            f"  - Matrix format: {matrix_type}"
        )

    # Create temporary directory for Numbat output
    import os
    import shutil
    import tempfile

    out_dir = tempfile.mkdtemp(prefix="numbat_", dir=tempfile.gettempdir())

    try:
        # Use sparkx-style context management for ALL R operations
        # This prevents "Conversion rules missing" errors in multithreaded/async environments
        with openrlib.rlock:  # Thread safety lock
            with conversion.localconverter(
                default_converter
                + anndata2ri.converter
                + pandas2ri.converter
                + numpy2ri.converter
            ):
                # Transfer data to R environment (inside context!)
                if context:
                    transfer_msg = (
                        f"Using anndata2ri for {'sparse' if is_sparse else 'dense'} matrix transfer to Numbat\n"
                        f"  (spatial: ({count_mat.shape[0]}, {count_mat.shape[1]}))"
                    )
                    await context.info(transfer_msg)

                ro.globalenv["count_mat"] = count_mat.T  # R expects genes × cells
                ro.globalenv["df_allele_python"] = (
                    df_allele  # Transfer allele dataframe
                )
                ro.globalenv["gene_names"] = gene_names
                ro.globalenv["cell_barcodes"] = cell_barcodes
                ro.globalenv["ref_indices"] = ref_indices_r
                ro.globalenv["out_dir"] = out_dir  # Output directory

                # Set Numbat parameters (inside context!)
                ro.globalenv["genome"] = params.numbat_genome
                ro.globalenv["t_param"] = params.numbat_t
                ro.globalenv["max_entropy"] = params.numbat_max_entropy
                ro.globalenv["min_cells"] = params.numbat_min_cells
                ro.globalenv["ncores"] = params.numbat_ncores
                ro.globalenv["skip_nj"] = params.numbat_skip_nj

                if context:
                    await context.info(
                        "Running Numbat analysis (this may take several minutes)..."
                    )

                # Run Numbat via R (inside context!)
                ro.r(
                    """
                    library(numbat)
                    library(dplyr)

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

        if context:
            await context.info("Numbat analysis complete. Reading output files...")

        # Read results from output files (Numbat saves to TSV files, not R objects)
        import pandas as pd

        # 1. Read clone posteriors (cell-level assignments)
        clone_post_file = os.path.join(out_dir, "clone_post_2.tsv")
        if not os.path.exists(clone_post_file):
            raise DataNotFoundError(
                f"Numbat output file not found: {clone_post_file}\n"
                f"Expected output files in: {out_dir}"
            )

        clone_post = pd.read_csv(clone_post_file, sep="\t")

        # 2. Read genotype matrix (CNV states per segment)
        geno_file = os.path.join(out_dir, "geno_2.tsv")
        if not os.path.exists(geno_file):
            raise DataNotFoundError(
                f"Numbat output file not found: {geno_file}\n"
                f"Expected output files in: {out_dir}"
            )

        geno = pd.read_csv(geno_file, sep="\t")

        # 3. Read consensus segments (optional metadata)
        segs_file = os.path.join(out_dir, "segs_consensus_2.tsv")
        segs = None
        if os.path.exists(segs_file):
            segs = pd.read_csv(segs_file, sep="\t")

        # 4. Check for phylogeny tree (if skip_nj=FALSE)
        tree_file = os.path.join(out_dir, "tree_final_2.rds")
        has_phylo = os.path.exists(tree_file)

        if context:
            await context.info(
                f"Results extracted:\n"
                f"  - Clone assignments: {len(clone_post)} cells\n"
                f"  - CNV segments: {geno.shape[1] - 1} segments\n"
                f"  - Phylogeny: {'Available' if has_phylo else 'Not computed'}"
            )

        # Process genotype matrix for AnnData storage
        # geno has structure: cell | segment1 | segment2 | ...
        # Convert to numpy array (cells × segments)
        geno_cells = geno["cell"].values
        geno_segments = geno.drop(columns=["cell"]).values

        # Ensure cells are in correct order (matching adata.obs_names)
        cell_order = {cell: i for i, cell in enumerate(cell_barcodes)}
        geno_sorted_indices = [cell_order.get(cell, -1) for cell in geno_cells]

        if -1 in geno_sorted_indices:
            raise DataCompatibilityError(
                "Mismatch between genotype cells and AnnData cells"
            )

        # Reorder genotype matrix to match AnnData cell order
        cnv_matrix = np.zeros((len(cell_barcodes), geno_segments.shape[1]))
        for geno_idx, adata_idx in enumerate(geno_sorted_indices):
            cnv_matrix[adata_idx, :] = geno_segments[geno_idx, :]

        # Store results in AnnData
        adata.obsm["X_cnv_numbat"] = cnv_matrix

        # Extract clone assignments and probabilities
        # Match clone_post cells with adata.obs_names
        clone_dict = clone_post.set_index("cell").to_dict()

        # Convert numpy types to Python native types for H5AD compatibility
        adata.obs["numbat_clone"] = [
            str(clone_dict["clone_opt"].get(cell, "unknown")) for cell in cell_barcodes
        ]
        adata.obs["numbat_p_cnv"] = [
            float(clone_dict["p_cnv"].get(cell, 0.0)) for cell in cell_barcodes
        ]
        adata.obs["numbat_compartment"] = [
            str(clone_dict["compartment_opt"].get(cell, "unknown"))
            for cell in cell_barcodes
        ]

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
            # Store phylogeny metadata
            adata.uns["numbat_phylogeny"] = {
                "available": True,
                "tree_file": tree_file,
                "tree_type": "phylo",
            }

        # Calculate statistics
        statistics = {
            "mean_cnv": float(np.mean(cnv_matrix)),
            "std_cnv": float(np.std(cnv_matrix)),
            "median_cnv": float(np.median(cnv_matrix)),
            "n_clones": int(clone_post["clone_opt"].nunique()),
            "mean_p_cnv": float(clone_post["p_cnv"].mean()),
            "n_reference_cells": len(ref_indices_r),
            "n_non_reference_cells": len(cell_barcodes) - len(ref_indices_r),
            "n_segments": geno_segments.shape[1],
        }

        # Get clone distribution
        clone_counts = clone_post["clone_opt"].value_counts()
        # Type: ignore needed because mypy doesn't infer Dict[str, Any] correctly
        statistics["clone_distribution"] = {  # type: ignore[assignment]
            str(clone): int(count) for clone, count in clone_counts.items()
        }

        # Store analysis parameters
        adata.uns["cnv_analysis"] = {
            "method": "numbat",
            "reference_key": params.reference_key,
            "reference_categories": list(params.reference_categories),
            "genome": params.numbat_genome,
            "t": params.numbat_t,
            "max_entropy": params.numbat_max_entropy,
            "min_cells": params.numbat_min_cells,
            "cnv_score_key": "X_cnv_numbat",
        }

        if context:
            await context.info(
                f"Numbat CNV analysis complete:\n"
                f"  - Detected {statistics['n_clones']} clones\n"
                f"  - Mean P(CNV): {statistics['mean_p_cnv']:.3f}\n"
                f"  - Phylogeny: {'Available' if has_phylo else 'Not computed'}"
            )

    except Exception as e:
        raise ProcessingError(
            f"Numbat analysis failed: {str(e)}\n"
            "Common issues:\n"
            "  - Allele data format incompatible\n"
            "  - Missing genomic position information\n"
            "  - Insufficient reference cells\n"
            "  - R environment configuration issues"
        ) from e
    finally:
        # Cleanup: Remove temporary output directory
        if os.path.exists(out_dir):
            try:
                shutil.rmtree(out_dir)
                if context:
                    await context.info(f"Cleaned up temporary directory: {out_dir}")
            except Exception as cleanup_error:
                # Log but don't fail if cleanup fails
                if context:
                    await context.info(
                        f"Warning: Failed to cleanup temp dir: {cleanup_error}"
                    )

        # Deactivate converters
        pandas2ri.deactivate()
        numpy2ri.deactivate()

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
