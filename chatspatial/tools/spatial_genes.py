"""
Spatial Variable Genes (SVG) identification for ChatSpatial MCP.

This module provides implementations for multiple SVG detection methods including GASTON,
SpatialDE, and SPARK-X, enabling comprehensive spatial transcriptomics analysis. Each method
offers distinct advantages for identifying genes with spatial expression patterns.

Methods Overview:
    - SPARK-X (default): Non-parametric statistical method, fastest execution, requires R
    - GASTON: Deep learning topographic mapping with isodepth analysis, GPU-accelerated
    - SpatialDE: Gaussian process-based kernel method, statistically rigorous

The module integrates these tools into the ChatSpatial MCP framework, handling data preparation,
execution, result formatting, and error management across different computational backends.
"""

import logging
import os
import shutil
import tempfile
from typing import TYPE_CHECKING, Any, Dict, List, Tuple

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)

# GASTON import will be done at runtime
GASTON_AVAILABLE = None
GASTON_IMPORT_ERROR = None

from ..models.analysis import SpatialVariableGenesResult
from ..models.data import SpatialVariableGenesParameters
from ..utils.error_handling import suppress_output


async def identify_spatial_genes(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None,
) -> SpatialVariableGenesResult:
    """
    Identify spatial variable genes using various computational methods.

    This is the main entry point for spatial gene detection, routing to the appropriate
    method based on params.method. Each method has different strengths:

    Method Selection Guide:
        - SPARK-X: Best for quick screening, handles large datasets efficiently
        - GASTON: Best for discovering complex spatial patterns (gradients, domains)
        - SpatialDE: Best for statistical rigor in publication-ready analyses

    Data Requirements:
        - SPARK-X: Works with raw counts or normalized data
        - GASTON: Requires normalized log-transformed data (will validate)
        - SpatialDE: Flexible, can handle both raw and normalized data

    Args:
        data_id: Dataset identifier in data store
        data_store: Dictionary containing loaded datasets
        params: Method-specific parameters (see SpatialVariableGenesParameters)
        context: MCP context for logging and status updates

    Returns:
        SpatialVariableGenesResult containing:
            - List of significant spatial genes
            - Statistical metrics (p-values, q-values)
            - Method-specific results (e.g., GASTON domains, isodepth)

    Raises:
        ValueError: If dataset not found or spatial coordinates missing
        ImportError: If required method dependencies not installed

    Performance Notes:
        - SPARK-X: ~2-5 min for 3000 spots × 20000 genes
        - GASTON: ~10-20 min (GPU recommended for larger datasets)
        - SpatialDE: ~15-30 min (scales with spot count squared)
    """
    if context:
        await context.info(
            f"Starting spatial variable genes identification using {params.method}"
        )

    # Get data
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Work with the original adata object, not a copy
    adata = data_store[data_id]["adata"]

    # Validate spatial coordinates
    if params.spatial_key not in adata.obsm:
        raise ValueError(
            f"Spatial coordinates not found in adata.obsm['{params.spatial_key}']"
        )

    # Extract spatial coordinates
    spatial_coords = adata.obsm[params.spatial_key]
    if spatial_coords.shape[1] != 2:
        raise ValueError("Spatial coordinates must be 2D (x, y)")

    # Log data information
    if context:
        await context.info(
            f"Processing data: {adata.n_obs} spots, {adata.n_vars} genes"
        )

    # Route to appropriate method
    if params.method == "gaston":
        return await _identify_spatial_genes_gaston(
            data_id, data_store, params, context
        )
    elif params.method == "spatialde":
        return await _identify_spatial_genes_spatialde(
            data_id, data_store, params, context
        )
    elif params.method == "sparkx":
        return await _identify_spatial_genes_sparkx(
            data_id, data_store, params, context
        )
    else:
        raise ValueError(
            f"Unsupported method: {params.method}. Available methods: gaston, spatialde, sparkx"
        )


async def _identify_spatial_genes_gaston(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None,
) -> SpatialVariableGenesResult:
    """
    Identify spatial variable genes using the GASTON deep learning method.

    GASTON (Generative Adversarial Spatial Transcriptomics Optimization Network) uses
    neural networks to learn topographic mappings from spatial coordinates to gene
    expression patterns. The method discovers both continuous spatial gradients and
    discrete spatial domains through isodepth analysis.

    Workflow:
        1. Feature engineering: Reduce expression to low-dimensional space (GLM-PCA/Pearson)
        2. Neural network training: Learn mapping from coordinates to expression features
        3. Isodepth extraction: Derive continuous depth values and discrete domains
        4. Gene classification: Identify genes with continuous vs discontinuous patterns
        5. Statistical analysis: Bin expression along isodepth and perform piecewise fitting

    Key Parameters:
        - preprocessing_method: 'glmpca' for count data, 'pearson_residuals' for normalized
        - n_domains: Number of spatial domains to identify (default: 5)
        - epochs: Training iterations (default: 10000)
        - continuous_quantile: Threshold for continuous gene classification (default: 0.9)

    Returns:
        Results including:
            - List of continuous gradient genes
            - List of discontinuous (domain-specific) genes
            - Isodepth values for visualization
            - Spatial domain assignments
            - Model performance metrics

    Requirements:
        - gaston-spatial package installed
        - Normalized, log-transformed data (validates automatically)
        - PyTorch for neural network operations
    """
    # Import dependencies at runtime
    import numpy as np

    # Check GASTON availability at runtime
    global GASTON_AVAILABLE, GASTON_IMPORT_ERROR
    if GASTON_AVAILABLE is None:
        try:
            import gaston  # noqa: F401
            from gaston import (
                binning_and_plotting,
                dp_related,  # noqa: F401
                neural_net,
                process_NN_output,
                segmented_fit,
                spatial_gene_classification,
            )

            GASTON_AVAILABLE = True
            GASTON_IMPORT_ERROR = None
        except ImportError as e:
            GASTON_AVAILABLE = False
            GASTON_IMPORT_ERROR = str(e)

    if not GASTON_AVAILABLE:
        error_msg = (
            f"GASTON is not available: {GASTON_IMPORT_ERROR}\n\n"
            "To use GASTON, please install it using one of these methods:\n"
            "1. pip install gaston-spatial\n"
            "2. Clone and install from source:\n"
            "   git clone https://github.com/raphael-group/GASTON.git\n"
            "   cd GASTON\n"
            "   pip install -e .\n\n"
            "Documentation: https://gaston.readthedocs.io/"
        )
        raise ImportError(error_msg)

    adata = data_store[data_id]["adata"]
    spatial_coords = adata.obsm[params.spatial_key]

    # Validate data preprocessing state
    data_max = adata.X.max() if hasattr(adata.X, "max") else np.max(adata.X)
    if data_max > 100:
        raise ValueError(
            "GASTON requires preprocessed data but raw counts detected. "
            "Please run basic preprocessing in preprocessing.py: "
            "1) sc.pp.normalize_total(adata, target_sum=1e4) "
            "2) sc.pp.log1p(adata)"
        )

    # GASTON-specific feature engineering (algorithm requirement)
    if context:
        await context.info(
            f"Applying GASTON-specific feature engineering using {params.preprocessing_method}"
        )

    if params.preprocessing_method == "glmpca":
        expression_features = await _gaston_feature_engineering_glmpca(
            adata, params.n_components, context
        )
    else:
        expression_features = await _gaston_feature_engineering_pearson(
            adata, params.n_components, context
        )

    # Create temporary directory for GASTON outputs
    temp_dir = tempfile.mkdtemp(prefix="gaston_")

    try:
        # Prepare data for GASTON
        coords_file = os.path.join(temp_dir, "spatial_coords.npy")
        expression_file = os.path.join(temp_dir, "expression_features.npy")

        np.save(coords_file, spatial_coords)
        np.save(expression_file, expression_features)

        if context:
            await context.info("Training GASTON neural network")

        # Train GASTON model
        model, loss_list, final_loss = await _train_gaston_model(
            coords_file, expression_file, params, temp_dir, context
        )

        if context:
            await context.info("Analyzing spatial patterns and gene classifications")

        # Analyze spatial patterns
        spatial_analysis = await _analyze_spatial_patterns(
            model, spatial_coords, expression_features, adata, params, context
        )

        # Store results in adata
        results_key = f"gaston_results_{params.random_seed}"
        await _store_results_in_adata(adata, spatial_analysis, results_key, context)

        # Create standardized result
        continuous_genes = spatial_analysis["continuous_genes"]  # Already a list now
        discontinuous_genes = spatial_analysis[
            "discontinuous_genes"
        ]  # Already a list now
        all_spatial_genes = list(set(continuous_genes + discontinuous_genes))

        # Create gene statistics (using isodepth correlation as primary statistic)
        gene_statistics = {}
        p_values = {}
        q_values = {}

        for gene in all_spatial_genes:
            # GASTON provides spatial classifications but not traditional statistical metrics
            # Be honest about unavailable statistics rather than using misleading placeholders
            gene_statistics[gene] = (
                None  # GASTON doesn't provide traditional statistics
            )
            p_values[gene] = None  # GASTON doesn't perform statistical testing
            q_values[gene] = None  # GASTON doesn't provide adjusted p-values

        # Create GASTON-specific results
        gaston_results = {
            "preprocessing_method": params.preprocessing_method,
            "n_components": params.n_components,
            "n_epochs_trained": params.epochs,
            "final_loss": final_loss,
            "spatial_hidden_layers": params.spatial_hidden_layers,
            "expression_hidden_layers": params.expression_hidden_layers,
            "n_spatial_domains": spatial_analysis["n_domains"],
            "continuous_gradient_genes": spatial_analysis["continuous_genes"],
            "discontinuous_genes": spatial_analysis["discontinuous_genes"],
            "model_performance": spatial_analysis["performance_metrics"],
        }

        result = SpatialVariableGenesResult(
            data_id=data_id,
            method="gaston",
            n_genes_analyzed=adata.n_vars,
            n_significant_genes=len(all_spatial_genes),
            spatial_genes=all_spatial_genes,
            gene_statistics=gene_statistics,
            p_values=p_values,
            q_values=q_values,
            results_key=results_key,
            gaston_results=gaston_results,
            isodepth_visualization={"plot_type": "gaston_isodepth"},
            spatial_domains_visualization={"plot_type": "gaston_domains"},
            top_genes_visualization={"plot_type": "gaston_genes"},
        )

        if context:
            await context.info("GASTON analysis completed successfully")
            await context.info(
                f"Found {len(continuous_genes)} genes with continuous gradients"
            )
            await context.info(
                f"Found {len(discontinuous_genes)} genes with discontinuities"
            )

        return result

    finally:
        # Clean up temporary directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


async def _identify_spatial_genes_spatialde(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None,
) -> SpatialVariableGenesResult:
    """
    Identify spatial variable genes using the SpatialDE statistical framework.

    SpatialDE employs Gaussian process regression with spatial kernels to decompose
    gene expression variance into spatial and non-spatial components. It provides
    rigorous statistical testing for spatial expression patterns with multiple
    testing correction.

    Method Details:
        - Models spatial correlation using squared exponential kernel
        - Tests significance via likelihood ratio test
        - Applies FDR correction for multiple testing
        - Returns both raw and adjusted p-values

    Key Parameters:
        - spatialde_normalized: Whether input is pre-normalized (default: True)
        - spatialde_kernel: Kernel type for spatial modeling (default: 'SE')
        - n_top_genes: Return only top N significant genes (optional)

    Data Handling:
        - If normalized=True: Expects log-transformed normalized data
        - If normalized=False: Applies SpatialDE-specific normalization to raw counts

    Returns:
        Results including:
            - List of significant spatial genes (q-value < 0.05)
            - Log-likelihood ratios as test statistics
            - Raw p-values and FDR-corrected q-values
            - Spatial correlation length scale per gene

    Requirements:
        - SpatialDE package installed
        - 2D spatial coordinates
        - Expression data (raw counts or normalized)
    """
    # Import dependencies at runtime
    import numpy as np
    import pandas as pd

    try:
        import SpatialDE
        from SpatialDE.util import qvalue
    except ImportError:
        raise ImportError(
            "SpatialDE not installed. Install with: pip install spatialde"
        )

    adata = data_store[data_id]["adata"]

    if context:
        await context.info("Running SpatialDE analysis")

    # Prepare data
    coords = pd.DataFrame(
        adata.obsm[params.spatial_key][:, :2],  # Ensure 2D coordinates
        columns=["x", "y"],
        index=adata.obs_names,
    )

    # Validate and get expression data for SpatialDE
    if params.spatialde_normalized:
        # Use pre-normalized data
        if "log1p" not in adata.uns_keys():
            raise ValueError(
                "SpatialDE normalized mode requires log-transformed data but log1p not found. "
                "Please run preprocessing in preprocessing.py: "
                "1) sc.pp.normalize_total(adata, target_sum=1e4) "
                "2) sc.pp.log1p(adata)"
            )

        counts = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
            columns=adata.var_names,
            index=adata.obs_names,
        )

        if context:
            await context.info("Using pre-normalized data for SpatialDE")
    else:
        # SpatialDE-specific normalization (algorithm requirement)
        if context:
            await context.info(
                "Applying SpatialDE-specific normalization (algorithm requirement)"
            )

        if adata.raw is not None:
            raw_counts = pd.DataFrame(
                (
                    adata.raw.X.toarray()
                    if hasattr(adata.raw.X, "toarray")
                    else adata.raw.X
                ),
                columns=adata.raw.var_names,
                index=adata.obs_names,
            )
        else:
            # Check if current data appears to be raw counts
            data_max = adata.X.max() if hasattr(adata.X, "max") else np.max(adata.X)
            if data_max <= 10:  # Likely already normalized
                raise ValueError(
                    "SpatialDE raw mode requires raw count data but normalized data detected. "
                    "Either switch to normalized mode or provide raw count data."
                )

            raw_counts = pd.DataFrame(
                adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
                columns=adata.var_names,
                index=adata.obs_names,
            )

        # SpatialDE-specific normalization and log transformation (algorithm requirement)
        total_counts = raw_counts.sum(axis=1)
        norm_counts = raw_counts.div(total_counts, axis=0) * np.median(total_counts)
        counts = np.log1p(norm_counts)

    # Run SpatialDE
    results = SpatialDE.run(coords.values, counts)

    # Multiple testing correction
    results["qval"] = qvalue(results["pval"].values, pi0=0.1)

    # Sort by q-value
    results = results.sort_values("qval")

    # Filter significant genes
    significant_genes = results[results["qval"] < 0.05]["g"].tolist()

    # Get top genes if requested
    if params.n_top_genes is not None:
        results = results.head(params.n_top_genes)
        significant_genes = results["g"].tolist()

    # Store results in adata
    results_key = f"spatialde_results_{data_id}"
    adata.var["spatialde_pval"] = results.set_index("g")["pval"]
    adata.var["spatialde_qval"] = results.set_index("g")["qval"]
    adata.var["spatialde_l"] = results.set_index("g")["l"]

    # Store scientific metadata for reproducibility
    from ..utils.metadata_storage import store_analysis_metadata

    store_analysis_metadata(
        adata,
        analysis_name="spatial_genes_spatialde",
        method="spatialde",
        parameters={
            "kernel": params.spatialde_kernel,
            "normalized": params.spatialde_normalized,
        },
        results_keys={
            "var": ["spatialde_pval", "spatialde_qval", "spatialde_l"],
            "obs": [],
            "obsm": [],
            "uns": [],
        },
        statistics={
            "n_genes_analyzed": len(results),
            "n_significant_genes": len(
                results[results["qval"] < params.pvalue_threshold]
            ),
        },
    )

    # Create gene statistics dictionaries
    gene_statistics = dict(zip(results["g"], results["LLR"]))  # Log-likelihood ratio
    p_values = dict(zip(results["g"], results["pval"]))
    q_values = dict(zip(results["g"], results["qval"]))

    # Create SpatialDE-specific results
    # Only return summary statistics (top 10 genes) to avoid exceeding MCP token limit
    top_results = results.head(10)
    spatialde_results = {
        "top_genes_summary": {
            "genes": top_results["g"].tolist(),
            "pvalues": top_results["pval"].tolist(),
            "qvalues": top_results["qval"].tolist(),
            "log_likelihood_ratios": top_results["LLR"].tolist(),
        },
        "kernel": params.spatialde_kernel,
        "normalized": params.spatialde_normalized,
        "n_significant_genes": len(significant_genes),
        "note": "Full results stored in adata.var['spatialde_pval', 'spatialde_qval', 'spatialde_l']",
    }

    result = SpatialVariableGenesResult(
        data_id=data_id,
        method="spatialde",
        n_genes_analyzed=len(results),
        n_significant_genes=len(significant_genes),
        spatial_genes=significant_genes,
        gene_statistics=gene_statistics,
        p_values=p_values,
        q_values=q_values,
        results_key=results_key,
        spatialde_results=spatialde_results,
    )

    if context:
        await context.info("SpatialDE analysis completed")
        await context.info(f"Found {len(significant_genes)} significant spatial genes")

    return result


async def _gaston_feature_engineering_glmpca(adata, n_components: int, context):
    """
    Perform GASTON-specific feature engineering using GLM-PCA.

    GLM-PCA (Generalized Linear Model PCA) is designed for count data and provides
    a better low-dimensional representation than standard PCA for UMI-based data.
    This preprocessing step is crucial for GASTON's neural network to learn
    meaningful spatial-expression relationships.

    Why GLM-PCA:
        - Designed for count data (Poisson/negative binomial)
        - Avoids log-transformation artifacts
        - Preserves biological signal better than standard PCA

    Args:
        adata: AnnData object with count matrix
        n_components: Number of latent dimensions (typically 10-20)
        context: MCP context for logging

    Returns:
        numpy.ndarray: Factor loadings of shape (n_obs, n_components)

    Note:
        Requires glmpca package: pip install glmpca
    """
    # Import dependencies at runtime
    import numpy as np

    try:
        from glmpca.glmpca import glmpca
    except ImportError:
        try:
            from glmpca import glmpca
        except ImportError:
            raise ImportError(
                "glmpca package is required for GASTON GLM-PCA feature engineering"
            )

    if context:
        await context.info(
            "Running GASTON GLM-PCA feature engineering (algorithm requirement)"
        )

    # Get count matrix
    if hasattr(adata.X, "toarray"):
        counts = adata.X.toarray()
    else:
        counts = adata.X.copy()

    # Ensure counts are integers and non-negative
    counts = np.round(np.maximum(counts, 0)).astype(int)

    if context:
        await context.info(
            f"Running GLM-PCA with {n_components} components on {counts.shape[0]} spots and {counts.shape[1]} genes"
        )

    # Run GLM-PCA
    glmpca_result = glmpca(counts.T, L=n_components, fam="poi")

    # Return the factors (PCs)
    return glmpca_result["factors"]


async def _gaston_feature_engineering_pearson(adata, n_components: int, context):
    """
    Perform GASTON-specific feature engineering using Pearson residuals.

    This alternative preprocessing computes Pearson residuals to normalize for
    technical variation, then applies PCA for dimensionality reduction. This
    approach is recommended for data that has already undergone quality control.

    Why Pearson Residuals:
        - Variance stabilization across expression levels
        - Removes mean-variance relationship
        - Better handles overdispersion in count data

    Args:
        adata: AnnData object (will be modified in-place)
        n_components: Number of principal components
        context: MCP context for logging

    Returns:
        numpy.ndarray: PCA coordinates of shape (n_obs, n_components)

    Note:
        Modifies adata in-place by adding normalized values
    """
    # Import dependencies at runtime
    import scanpy as sc

    if context:
        await context.info(
            "Computing GASTON Pearson residuals feature engineering (algorithm requirement)"
        )

    # GASTON-specific Pearson residuals computation (algorithm requirement)
    sc.experimental.pp.normalize_pearson_residuals(adata)

    # GASTON-specific PCA on residuals (algorithm requirement)
    sc.tl.pca(adata, n_comps=n_components)

    return adata.obsm["X_pca"]


async def _train_gaston_model(
    coords_file: str,
    expression_file: str,
    params: SpatialVariableGenesParameters,
    output_dir: str,
    context,
) -> Tuple[Any, List[float], float]:
    """
    Train the GASTON neural network model.

    The model architecture consists of two branches:
        - Spatial branch (f_S): Maps 2D coordinates to embedding space
        - Expression branch (f_A): Maps embedding to expression features

    Training optimizes the combined network to predict expression features
    from spatial coordinates, learning the topographic organization of the tissue.

    Args:
        coords_file: Path to numpy file with spatial coordinates (n_spots × 2)
        expression_file: Path to numpy file with expression features (n_spots × n_components)
        params: Training configuration including architecture and hyperparameters
        output_dir: Directory for saving model checkpoints
        context: MCP context for logging

    Returns:
        Tuple of:
            - Trained PyTorch model
            - List of loss values per epoch
            - Final training loss

    Architecture Parameters:
        - spatial_hidden_layers: Hidden units for spatial branch
        - expression_hidden_layers: Hidden units for expression branch
        - embedding_size: Dimension of latent space connecting branches
    """
    # Import dependencies at runtime
    import numpy as np
    from gaston import neural_net

    # Load data
    S = np.load(coords_file)
    A = np.load(expression_file)

    # Convert to torch tensors and rescale
    S_torch, A_torch = neural_net.load_rescale_input_data(S, A)

    # Train model using official parameters
    checkpoint_interval = params.checkpoint_interval  # Official default: 500

    model, loss_list = neural_net.train(
        S_torch,
        A_torch,
        S_hidden_list=params.spatial_hidden_layers,
        A_hidden_list=params.expression_hidden_layers,
        epochs=params.epochs,
        checkpoint=checkpoint_interval,
        save_dir=output_dir,
        optim=params.optimizer,
        lr=params.learning_rate,
        seed=params.random_seed,
        save_final=True,
        embed_size=params.embedding_size,
        sigma=params.sigma,
        batch_size=params.batch_size,
    )

    final_loss = float(loss_list[-1]) if len(loss_list) > 0 else 0.0

    return model, loss_list, final_loss


async def _analyze_spatial_patterns(
    model,
    spatial_coords,
    expression_features,
    adata,
    params: SpatialVariableGenesParameters,
    context,
) -> Dict[str, Any]:
    """
    Analyze spatial patterns from trained GASTON model.

    This executes the complete GASTON post-training analysis pipeline:

    1. Isodepth Extraction:
        - Derives continuous depth values representing position in expression manifold
        - Identifies discrete spatial domains from depth discontinuities

    2. Gene Binning:
        - Bins gene expression along isodepth axis within each domain
        - Applies UMI threshold filtering for reliable genes

    3. Piecewise Linear Fitting:
        - Fits segmented regression models to binned expression
        - Tests for significant slopes and breakpoints

    4. Gene Classification:
        - Continuous: Genes with smooth gradients along isodepth
        - Discontinuous: Genes with sharp transitions between domains

    Args:
        model: Trained GASTON PyTorch model
        spatial_coords: Original spatial coordinates
        expression_features: Low-dimensional expression features
        adata: AnnData object for accessing gene expression
        params: Analysis parameters including thresholds
        context: MCP context for logging

    Returns:
        Dictionary with:
            - isodepth: Continuous depth values per spot
            - spatial_domains: Discrete domain labels
            - continuous_genes: List of gradient genes
            - discontinuous_genes: List of domain-specific genes
            - performance_metrics: Model evaluation statistics
    """
    # Import dependencies at runtime
    import numpy as np
    import pandas as pd
    import torch
    from gaston import (
        binning_and_plotting,
        dp_related,
        segmented_fit,
        spatial_gene_classification,
    )

    if context:
        await context.info("Processing neural network output following GASTON tutorial")

    # Step 1: Use dp_related.get_isodepth_labels to get isodepth and labels
    # This is the official GASTON way according to the tutorial
    S_torch = torch.tensor(spatial_coords, dtype=torch.float32)
    A_torch = torch.tensor(expression_features, dtype=torch.float32)

    # Use GASTON's official method to get isodepth and labels
    gaston_isodepth, gaston_labels = dp_related.get_isodepth_labels(
        model, A_torch, S_torch, params.n_domains
    )

    if context:
        actual_domains = len(np.unique(gaston_labels))
        await context.info(f"GASTON identified {actual_domains} spatial domains")
        await context.info(
            f"Isodepth range: [{gaston_isodepth.min():.3f}, {gaston_isodepth.max():.3f}]"
        )

    # Step 2: Get model predictions for performance metrics
    with torch.no_grad():
        predictions = model(S_torch).numpy()

    # Step 3: Prepare data for GASTON analysis
    # Get original count matrix (need raw counts for Poisson regression)
    if hasattr(adata.X, "toarray"):
        counts_mat = adata.X.toarray()  # GASTON expects N x G matrix
    else:
        counts_mat = adata.X

    # Ensure counts are non-negative integers
    counts_mat = np.maximum(counts_mat, 0).astype(int)
    gene_labels = adata.var_names.values

    # Create dummy cell type dataframe (for all cell types analysis)
    cell_type_df = pd.DataFrame(
        {"All": np.ones(len(spatial_coords))}, index=adata.obs_names
    )

    # Step 4: Perform binning using GASTON's binning function
    if context:
        await context.info(f"Running GASTON binning with {len(gene_labels)} genes")
        await context.info(f"UMI threshold: {params.umi_threshold}")

    binning_output = binning_and_plotting.bin_data(
        counts_mat=counts_mat,
        gaston_labels=gaston_labels,
        gaston_isodepth=gaston_isodepth,
        cell_type_df=cell_type_df,
        gene_labels=gene_labels,
        num_bins=params.num_bins,
        umi_threshold=params.umi_threshold,
        pc=0,  # No pseudocount
        pc_exposure=True,
    )

    if context:
        await context.info(
            f"Binning completed. Analyzing {len(binning_output['gene_labels_idx'])} genes"
        )

    # Step 5: Perform piecewise linear fitting using GASTON's segmented fit
    if context:
        await context.info("Running piecewise linear fitting")

    pw_fit_dict = segmented_fit.pw_linear_fit(
        counts_mat=counts_mat,
        gaston_labels=gaston_labels,
        gaston_isodepth=gaston_isodepth,
        cell_type_df=cell_type_df,
        ct_list=[],  # Empty list for all cell types analysis only
        umi_threshold=params.umi_threshold,
        t=params.pvalue_threshold,
        isodepth_mult_factor=params.isodepth_mult_factor,
        reg=params.regularization,
        zero_fit_threshold=params.zero_fit_threshold,
    )

    if context:
        await context.info("Piecewise linear fitting completed")
        await context.info(
            "Classifying genes into continuous and discontinuous patterns"
        )

    # Step 6: Classify genes using GASTON's classification functions
    continuous_genes = spatial_gene_classification.get_cont_genes(
        pw_fit_dict, binning_output, q=params.continuous_quantile
    )

    discontinuous_genes = spatial_gene_classification.get_discont_genes(
        pw_fit_dict, binning_output, q=params.discontinuous_quantile
    )

    if context:
        await context.info(
            f"Found {len(continuous_genes)} genes with continuous gradients"
        )
        await context.info(
            f"Found {len(discontinuous_genes)} genes with discontinuities"
        )

    # Performance metrics
    mse = np.mean((predictions - expression_features) ** 2)
    r2 = 1 - (
        np.sum((expression_features - predictions) ** 2)
        / np.sum((expression_features - np.mean(expression_features, axis=0)) ** 2)
    )

    performance_metrics = {
        "mse": float(mse),
        "r2": float(r2),
        "isodepth_range": [float(gaston_isodepth.min()), float(gaston_isodepth.max())],
        "isodepth_std": float(gaston_isodepth.std()),
    }

    # Convert complex structures to simple lists for MCP output
    continuous_genes_list = []
    discontinuous_genes_list = []

    if continuous_genes:
        continuous_genes_list = (
            list(continuous_genes.keys()) if hasattr(continuous_genes, "keys") else []
        )

    if discontinuous_genes:
        discontinuous_genes_list = (
            list(discontinuous_genes.keys())
            if hasattr(discontinuous_genes, "keys")
            else []
        )

    return {
        "isodepth": gaston_isodepth,
        "spatial_domains": gaston_labels,
        "n_domains": params.n_domains,
        "continuous_genes": continuous_genes_list,
        "discontinuous_genes": discontinuous_genes_list,
        "performance_metrics": performance_metrics,
    }


# Note: Visualization functions have been moved to visualization.py
# Use visualize_data tool with plot_type="gaston_isodepth", "gaston_domains", or "gaston_genes"


async def _store_results_in_adata(
    adata, spatial_analysis: Dict[str, Any], results_key: str, context
):
    """
    Store GASTON analysis results in the AnnData object.

    Persists key results for downstream analysis and visualization:

    Stored Annotations:
        - .obs[f'{results_key}_isodepth']: Continuous depth values
        - .obs[f'{results_key}_spatial_domains']: Discrete domain labels
        - .obsm[f'{results_key}_embedding']: Isodepth as 1D embedding
        - .uns[f'{results_key}_metadata']: Analysis metadata and metrics

    Args:
        adata: AnnData object to annotate
        spatial_analysis: Results from _analyze_spatial_patterns
        results_key: Prefix for result keys (includes random seed)
        context: MCP context for logging

    Note:
        Results are stored with unique keys to allow multiple GASTON runs
        with different parameters on the same dataset.
    """

    # Store isodepth values
    adata.obs[f"{results_key}_isodepth"] = spatial_analysis["isodepth"]

    # Store spatial domains
    adata.obs[f"{results_key}_spatial_domains"] = spatial_analysis[
        "spatial_domains"
    ].astype(str)

    # Store predictions (if available - removed from spatial results to reduce MCP output size)
    # adata.obsm[f"{results_key}_predictions"] = spatial_analysis.get('predictions', None)

    # Store spatial embedding (isodepth as 1D embedding)
    adata.obsm[f"{results_key}_embedding"] = spatial_analysis["isodepth"].reshape(-1, 1)

    # Store metadata
    adata.uns[f"{results_key}_metadata"] = {
        "method": "GASTON",
        "n_domains": spatial_analysis["n_domains"],
        "performance_metrics": spatial_analysis["performance_metrics"],
    }

    if context:
        await context.info(f"Results stored in adata with key prefix: {results_key}")


async def _identify_spatial_genes_sparkx(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None,
) -> SpatialVariableGenesResult:
    """
    Identify spatial variable genes using the SPARK-X non-parametric method.

    SPARK-X is an efficient non-parametric method for detecting spatially variable
    genes without assuming specific distribution models. It uses spatial covariance
    testing and is particularly effective for large-scale datasets. The method is
    implemented in R and accessed via rpy2.

    Method Advantages:
        - Non-parametric: No distributional assumptions required
        - Computationally efficient: Scales well with gene count
        - Robust: Handles various spatial patterns effectively
        - Flexible: Works with both single and mixture spatial kernels

    Key Parameters:
        - sparkx_option: 'single' or 'mixture' kernel (default: 'mixture')
        - sparkx_percentage: Min percentage of cells expressing gene (default: 0.1)
        - sparkx_min_total_counts: Min total counts per gene (default: 10)
        - sparkx_num_core: Number of CPU cores for parallel processing

    Data Processing:
        - Automatically filters low-expression genes based on parameters
        - Uses raw counts when available (adata.raw), otherwise current matrix
        - Handles duplicate gene names by adding suffixes

    Returns:
        Results including:
            - List of significant spatial genes (adjusted p-value < 0.05)
            - Raw p-values from spatial covariance test
            - Bonferroni-adjusted p-values
            - Results dataframe with all tested genes

    Requirements:
        - R installation with SPARK package
        - rpy2 Python package for R integration
        - Raw count data preferred (will use adata.raw if available)

    Performance:
        - Fastest among the three methods
        - ~2-5 minutes for typical datasets (3000 spots × 20000 genes)
        - Memory efficient through gene filtering
    """
    # Import dependencies at runtime
    import numpy as np
    import pandas as pd

    try:
        from rpy2 import robjects as ro
        from rpy2.robjects import conversion, default_converter
        from rpy2.robjects.packages import importr
        from rpy2.rinterface_lib import openrlib  # For thread safety
    except ImportError:
        raise ImportError("rpy2 not installed. Install with: pip install rpy2")

    if context:
        await context.info("Running SPARK-X non-parametric analysis")

    adata = data_store[data_id]["adata"]

    # Prepare spatial coordinates - SPARK needs data.frame format
    coords_array = adata.obsm[params.spatial_key][:, :2].astype(float)
    n_spots, n_genes = adata.shape

    if context:
        await context.info(f"Preparing data: {n_spots} spots × {n_genes} genes")

    # Get count matrix - use raw counts if available, otherwise current matrix
    if adata.raw is not None:
        # Use raw counts for SPARK
        if hasattr(adata.raw.X, "toarray"):
            counts_matrix = adata.raw.X.toarray()
        else:
            counts_matrix = adata.raw.X.copy()
        # Use raw gene names
        gene_names = [str(name) for name in adata.raw.var_names]
        n_genes = len(gene_names)
    else:
        # Fallback to current matrix
        if hasattr(adata.X, "toarray"):
            counts_matrix = adata.X.toarray()
        else:
            counts_matrix = adata.X.copy()
        gene_names = [str(name) for name in adata.var_names]
        n_genes = len(gene_names)

    # Ensure gene names are unique (required for SPARK-X R rownames)
    if len(gene_names) != len(set(gene_names)):
        from collections import Counter

        gene_counts = Counter(gene_names)
        unique_names = []
        seen_counts = {}
        for gene in gene_names:
            if gene_counts[gene] > 1:
                # Add suffix for duplicates
                if gene not in seen_counts:
                    seen_counts[gene] = 0
                    unique_names.append(gene)
                else:
                    seen_counts[gene] += 1
                    unique_names.append(f"{gene}_{seen_counts[gene]}")
            else:
                unique_names.append(gene)
        gene_names = unique_names
        if context:
            await context.info(
                f"Made duplicate gene names unique (found {sum(1 for c in gene_counts.values() if c > 1)} duplicates)"
            )

    # Ensure counts are non-negative integers
    counts_matrix = np.maximum(counts_matrix, 0).astype(int)

    if context:
        await context.info(
            f"Using {'raw' if adata.raw is not None else 'current'} count matrix"
        )

    # Apply gene filtering based on SPARK-X parameters (like CreateSPARKObject in R)
    percentage = params.sparkx_percentage
    min_total_counts = params.sparkx_min_total_counts

    # Calculate total counts per gene
    gene_totals = counts_matrix.sum(axis=0)
    n_expressed = (counts_matrix > 0).sum(axis=0)

    # Filter genes: must be expressed in at least percentage of cells AND have min total counts
    min_cells = int(np.ceil(n_spots * percentage))
    keep_genes = (n_expressed >= min_cells) & (gene_totals >= min_total_counts)

    if keep_genes.sum() < len(gene_names):
        # Apply filtering
        counts_matrix = counts_matrix[:, keep_genes]
        gene_names = [gene for gene, keep in zip(gene_names, keep_genes) if keep]
        n_filtered = keep_genes.sum()

        if context:
            await context.info(
                f"Filtered to {n_filtered}/{len(keep_genes)} genes (>{percentage*100:.0f}% cells, >{min_total_counts} counts)"
            )

    # Update gene count after filtering
    n_genes = len(gene_names)

    # Transpose for SPARK format (genes × spots)
    counts_transposed = counts_matrix.T

    if context:
        await context.info(
            f"Count matrix shape: {counts_transposed.shape} (genes × spots)"
        )
        await context.info(
            f"Passing {n_genes} genes to SPARK-X for analysis"
        )

    # Create spot names
    spot_names = [str(name) for name in adata.obs_names]

    # Wrap ALL R operations in thread lock and localconverter for proper contextvars handling
    # This prevents "Conversion rules missing" errors in multithreaded/async environments
    with openrlib.rlock:  # Thread safety lock
        with conversion.localconverter(default_converter):  # Conversion context
            # Import SPARK package inside context (FIX for contextvars issue)
            try:
                spark = importr("SPARK")
            except Exception as e:
                raise ImportError(
                    f"SPARK not installed in R. Install with: install.packages('SPARK'). Error: {e}"
                )

            # Convert to R format (already in context)
            # Count matrix: genes × spots
            r_counts = ro.r.matrix(
                ro.IntVector(counts_transposed.flatten()),
                nrow=n_genes,
                ncol=n_spots,
                byrow=True,
            )
            r_counts.rownames = ro.StrVector(gene_names)
            r_counts.colnames = ro.StrVector(spot_names)

            # Coordinates as data.frame (SPARK requirement)
            coords_df = pd.DataFrame(coords_array, columns=["x", "y"], index=spot_names)
            r_coords = ro.r["data.frame"](
                x=ro.FloatVector(coords_df["x"]),
                y=ro.FloatVector(coords_df["y"]),
                row_names=ro.StrVector(coords_df.index),
            )

            if context:
                await context.info(
                    "Running SPARK-X analysis using sparkx function (output suppressed for MCP compatibility)"
                )

            try:
                # Execute SPARK-X analysis inside context (FIX for contextvars issue)
                # Keep suppress_output for MCP communication compatibility
                with suppress_output():
                    results = spark.sparkx(
                        count_in=r_counts,
                        locus_in=r_coords,
                        X_in=ro.NULL,  # No additional covariates (could be extended in future)
                        numCores=params.sparkx_num_core,
                        option=params.sparkx_option,
                        verbose=False,  # Ensure verbose is off for cleaner MCP communication
                    )

                if context:
                    await context.info("SPARK-X analysis completed successfully")

                # Extract p-values from results (inside context for proper conversion)
                try:
                    pvals = results.rx2("res_mtest")
                    if pvals:
                        # SPARK-X returns res_mtest as a data.frame with columns:
                        # - combinedPval: combined p-values across kernels
                        # - adjustedPval: adjusted p-values
                        # We need to extract the combinedPval column

                        # Check if it's a data.frame (which it should be for SPARK-X)
                        is_dataframe = ro.r["is.data.frame"](pvals)[0]

                        if is_dataframe:
                            # Extract combinedPval column
                            combined_pvals = ro.r["$"](pvals, "combinedPval")
                            pval_list = [float(p) for p in combined_pvals]

                            # Also extract adjustedPval if available
                            adjusted_pvals = ro.r["$"](pvals, "adjustedPval")
                            adjusted_pval_list = [float(p) for p in adjusted_pvals]

                            # Create results dataframe
                            results_df = pd.DataFrame(
                                {
                                    "gene": gene_names[: len(pval_list)],
                                    "pvalue": pval_list,
                                    "adjusted_pvalue": adjusted_pval_list,
                                }
                            )
                        else:
                            # Fallback for older format (numeric vector)
                            pval_list = []
                            try:
                                pvals_numeric = ro.r["as.numeric"](pvals)
                            except:
                                pvals_numeric = pvals

                            for i in range(len(pvals_numeric)):
                                val = pvals_numeric[i]
                                if hasattr(val, "__len__") and hasattr(
                                    val, "__getitem__"
                                ):
                                    try:
                                        pval_list.append(float(val[0]))
                                    except:
                                        pval_list.append(float(val))
                                else:
                                    pval_list.append(float(val))

                            # Create results dataframe
                            results_df = pd.DataFrame(
                                {
                                    "gene": gene_names[: len(pval_list)],
                                    "pvalue": pval_list,
                                }
                            )

                            # Add adjusted p-values (Bonferroni correction)
                            n_tests = len(pval_list)
                            results_df["adjusted_pvalue"] = results_df[
                                "pvalue"
                            ].apply(lambda p: min(p * n_tests, 1.0))

                        if context:
                            await context.info(
                                f"Extracted results for {len(results_df)} genes"
                            )
                            # Warn if returned genes much fewer than input genes
                            if len(results_df) < n_genes * 0.5:
                                await context.warning(
                                    f"SPARK-X returned results for only {len(results_df)}/{n_genes} genes. "
                                    f"This may indicate a problem with the R environment, SPARK package, or input data. "
                                    f"Consider checking R logs or trying alternative methods (GASTON, SpatialDE)."
                                )
                    else:
                        # SPARK-X results format not recognized - fail honestly instead of fake results
                        error_msg = (
                            "SPARK-X results format not recognized. Expected 'res_mtest' component with p-values. "
                            "This may indicate an issue with the SPARK-X R package, rpy2 integration, or input data. "
                            "Please check the R environment and SPARK-X installation."
                        )
                        if context:
                            await context.error(error_msg)
                        raise RuntimeError(error_msg)

                except Exception as e:
                    # P-value extraction failed - fail honestly instead of creating fake results
                    error_msg = (
                        f"SPARK-X p-value extraction failed: {e}. "
                        f"This indicates an issue with R-Python communication or SPARK-X result format. "
                        f"Please check rpy2 installation and R package versions."
                    )
                    if context:
                        await context.error(error_msg)
                    raise RuntimeError(error_msg)

            except Exception as e:
                if context:
                    await context.info(f"SPARK-X analysis failed: {e}")
                raise RuntimeError(f"SPARK-X analysis failed: {e}")

    # Sort by adjusted p-value
    results_df = results_df.sort_values("adjusted_pvalue")

    # Filter significant genes
    significant_genes = results_df[results_df["adjusted_pvalue"] < 0.05][
        "gene"
    ].tolist()

    # Get top genes if requested
    if params.n_top_genes is not None:
        results_df = results_df.head(params.n_top_genes)
        significant_genes = results_df["gene"].tolist()

    # Store results in adata
    results_key = f"sparkx_results_{data_id}"
    adata.var["sparkx_pval"] = pd.Series(
        dict(zip(results_df["gene"], results_df["pvalue"])), name="sparkx_pval"
    ).reindex(adata.var_names, fill_value=1.0)

    adata.var["sparkx_qval"] = pd.Series(
        dict(zip(results_df["gene"], results_df["adjusted_pvalue"])), name="sparkx_qval"
    ).reindex(adata.var_names, fill_value=1.0)

    # Store scientific metadata for reproducibility
    from ..utils.metadata_storage import store_analysis_metadata

    store_analysis_metadata(
        adata,
        analysis_name="spatial_genes_sparkx",
        method="sparkx",
        parameters={
            "num_core": params.sparkx_num_core,
            "percentage": params.sparkx_percentage,
            "min_total_counts": params.sparkx_min_total_counts,
            "option": params.sparkx_option,
        },
        results_keys={
            "var": ["sparkx_pval", "sparkx_qval"],
            "obs": [],
            "obsm": [],
            "uns": [],
        },
        statistics={
            "n_genes_analyzed": len(results_df),
            "n_significant_genes": len(significant_genes),
        },
    )

    # Create gene statistics dictionaries
    gene_statistics = dict(zip(results_df["gene"], results_df["pvalue"]))
    p_values = dict(zip(results_df["gene"], results_df["pvalue"]))
    q_values = dict(zip(results_df["gene"], results_df["adjusted_pvalue"]))

    # Create SPARK-X specific results
    # Only return summary statistics (top 10 genes) to avoid exceeding MCP token limit
    top_results = results_df.head(10)
    sparkx_results = {
        "top_genes_summary": {
            "genes": top_results["gene"].tolist(),
            "pvalues": top_results["pvalue"].tolist(),
            "adjusted_pvalues": top_results["adjusted_pvalue"].tolist(),
            "combined_pvalues": (
                top_results["combined_pvalue"].tolist()
                if "combined_pvalue" in top_results.columns
                else None
            ),
        },
        "method": "sparkx",
        "num_core": params.sparkx_num_core,
        "option": params.sparkx_option,
        "n_significant_genes": len(significant_genes),
        "data_format": "genes_x_spots",
        "note": "Full results stored in adata.var['sparkx_pval', 'sparkx_qval']",
    }

    result = SpatialVariableGenesResult(
        data_id=data_id,
        method="sparkx",
        n_genes_analyzed=len(results_df),
        n_significant_genes=len(significant_genes),
        spatial_genes=significant_genes,
        gene_statistics=gene_statistics,
        p_values=p_values,
        q_values=q_values,
        results_key=results_key,
        sparkx_results=sparkx_results,
    )

    if context:
        await context.info("SPARK-X analysis completed successfully")
        await context.info(f"Analyzed {len(results_df)} genes")
        await context.info(
            f"Found {len(significant_genes)} significant spatial genes (q < 0.05)"
        )

    return result
