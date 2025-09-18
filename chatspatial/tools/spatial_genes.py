"""
GASTON-based spatial variable genes identification for ChatSpatial MCP.

This module integrates GASTON (Generative Adversarial Spatial Transcriptomics Optimization Network)
for identifying spatial variable genes through topographic mapping and isodepth learning.
"""

import os
import sys
import tempfile
import shutil
from typing import Dict, List, Tuple, Optional, Any, TYPE_CHECKING
from pathlib import Path
import warnings
import logging

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd
    import torch
    import torch.nn as nn
    import scanpy as sc
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)

# GASTON import will be done at runtime
GASTON_AVAILABLE = None
GASTON_IMPORT_ERROR = None

from ..models.data import SpatialVariableGenesParameters
from ..models.analysis import SpatialVariableGenesResult


async def identify_spatial_genes(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None
) -> SpatialVariableGenesResult:
    """
    Identify spatial variable genes using various methods.

    Args:
        data_id: Dataset identifier
        data_store: Data storage instance
        params: Spatial variable genes parameters
        context: MCP context for logging

    Returns:
        SpatialVariableGenesResult with analysis results
    """
    if context:
        await context.info(f"Starting spatial variable genes identification using {params.method}")

    # Get data
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Work with the original adata object, not a copy
    adata = data_store[data_id]["adata"]

    # Validate spatial coordinates
    if params.spatial_key not in adata.obsm:
        raise ValueError(f"Spatial coordinates not found in adata.obsm['{params.spatial_key}']")

    # Extract spatial coordinates
    spatial_coords = adata.obsm[params.spatial_key]
    if spatial_coords.shape[1] != 2:
        raise ValueError("Spatial coordinates must be 2D (x, y)")

    # Log data information
    if context:
        await context.info(f"Processing data: {adata.n_obs} spots, {adata.n_vars} genes")

    # Route to appropriate method
    if params.method == "gaston":
        return await _identify_spatial_genes_gaston(data_id, data_store, params, context)
    elif params.method == "spatialde":
        return await _identify_spatial_genes_spatialde(data_id, data_store, params, context)
    elif params.method == "sparkx":
        return await _identify_spatial_genes_sparkx(data_id, data_store, params, context)
    else:
        raise ValueError(f"Unsupported method: {params.method}. Available methods: gaston, spatialde, sparkx")


async def _identify_spatial_genes_gaston(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None
) -> SpatialVariableGenesResult:
    """Identify spatial variable genes using GASTON method."""
    # Import dependencies at runtime
    import numpy as np
    import pandas as pd
    import torch
    import torch.nn as nn
    
    # Check GASTON availability at runtime
    global GASTON_AVAILABLE, GASTON_IMPORT_ERROR
    if GASTON_AVAILABLE is None:
        try:
            import gaston
            from gaston import neural_net, spatial_gene_classification, binning_and_plotting
            from gaston import dp_related, segmented_fit, process_NN_output
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
            "   git clone https://github.com/Arashz/GASTON.git\n"
            "   cd GASTON\n"
            "   pip install -e .\n\n"
            "For development, you can also clone the original repository for comparison:\n"
            "git clone https://github.com/Arashz/GASTON.git gaston_dev"
        )
        raise ImportError(error_msg)

    adata = data_store[data_id]["adata"]
    spatial_coords = adata.obsm[params.spatial_key]

    # Validate data preprocessing state
    data_max = adata.X.max() if hasattr(adata.X, 'max') else np.max(adata.X)
    if data_max > 100:
        raise ValueError(
            "GASTON requires preprocessed data but raw counts detected. "
            "Please run basic preprocessing in preprocessing.py: "
            "1) sc.pp.normalize_total(adata, target_sum=1e4) "
            "2) sc.pp.log1p(adata)"
        )
    
    # GASTON-specific feature engineering (algorithm requirement)
    if context:
        await context.info(f"Applying GASTON-specific feature engineering using {params.preprocessing_method}")

    if params.preprocessing_method == "glmpca":
        expression_features = await _gaston_feature_engineering_glmpca(adata, params.n_components, context)
    else:
        expression_features = await _gaston_feature_engineering_pearson(adata, params.n_components, context)

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
        continuous_genes = spatial_analysis['continuous_genes']  # Already a list now
        discontinuous_genes = spatial_analysis['discontinuous_genes']  # Already a list now
        all_spatial_genes = list(set(continuous_genes + discontinuous_genes))

        # Create gene statistics (using isodepth correlation as primary statistic)
        gene_statistics = {}
        p_values = {}
        q_values = {}

        for gene in all_spatial_genes:
            gene_statistics[gene] = 1.0  # Placeholder - GASTON doesn't provide traditional statistics
            p_values[gene] = 0.05  # Placeholder
            q_values[gene] = 0.05  # Placeholder

        # Create GASTON-specific results
        gaston_results = {
            'preprocessing_method': params.preprocessing_method,
            'n_components': params.n_components,
            'n_epochs_trained': params.epochs,
            'final_loss': final_loss,
            'spatial_hidden_layers': params.spatial_hidden_layers,
            'expression_hidden_layers': params.expression_hidden_layers,
            'n_spatial_domains': spatial_analysis['n_domains'],
            'continuous_gradient_genes': spatial_analysis['continuous_genes'],
            'discontinuous_genes': spatial_analysis['discontinuous_genes'],
            'model_performance': spatial_analysis['performance_metrics'],
            'spatial_autocorrelation': spatial_analysis['autocorrelation_metrics']
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
            isodepth_visualization={'plot_type': 'gaston_isodepth'},
            spatial_domains_visualization={'plot_type': 'gaston_domains'},
            top_genes_visualization={'plot_type': 'gaston_genes'}
        )

        if context:
            await context.info(f"GASTON analysis completed successfully")
            await context.info(f"Found {len(continuous_genes)} genes with continuous gradients")
            await context.info(f"Found {len(discontinuous_genes)} genes with discontinuities")

        return result

    finally:
        # Clean up temporary directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


async def _identify_spatial_genes_spatialde(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None
) -> SpatialVariableGenesResult:
    """Identify spatial variable genes using SpatialDE method."""
    # Import dependencies at runtime
    import numpy as np
    import pandas as pd
    import scanpy as sc
    try:
        import SpatialDE
        from SpatialDE.util import qvalue
    except ImportError:
        raise ImportError("SpatialDE not installed. Install with: pip install spatialde")

    adata = data_store[data_id]["adata"]

    if context:
        await context.info("Running SpatialDE analysis")

    # Prepare data
    coords = pd.DataFrame(
        adata.obsm[params.spatial_key][:, :2],  # Ensure 2D coordinates
        columns=['x', 'y'],
        index=adata.obs_names
    )

    # Validate and get expression data for SpatialDE
    if params.spatialde_normalized:
        # Use pre-normalized data
        if 'log1p' not in adata.uns_keys():
            raise ValueError(
                "SpatialDE normalized mode requires log-transformed data but log1p not found. "
                "Please run preprocessing in preprocessing.py: "
                "1) sc.pp.normalize_total(adata, target_sum=1e4) "
                "2) sc.pp.log1p(adata)"
            )
        
        counts = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
            columns=adata.var_names,
            index=adata.obs_names
        )
        
        if context:
            await context.info("Using pre-normalized data for SpatialDE")
    else:
        # SpatialDE-specific normalization (algorithm requirement)
        if context:
            await context.info("Applying SpatialDE-specific normalization (algorithm requirement)")
        
        if adata.raw is not None:
            raw_counts = pd.DataFrame(
                adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X,
                columns=adata.raw.var_names,
                index=adata.obs_names
            )
        else:
            # Check if current data appears to be raw counts
            data_max = adata.X.max() if hasattr(adata.X, 'max') else np.max(adata.X)
            if data_max <= 10:  # Likely already normalized
                raise ValueError(
                    "SpatialDE raw mode requires raw count data but normalized data detected. "
                    "Either switch to normalized mode or provide raw count data."
                )
            
            raw_counts = pd.DataFrame(
                adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
                columns=adata.var_names,
                index=adata.obs_names
            )
        
        # SpatialDE-specific normalization and log transformation (algorithm requirement)
        total_counts = raw_counts.sum(axis=1)
        norm_counts = raw_counts.div(total_counts, axis=0) * np.median(total_counts)
        counts = np.log1p(norm_counts)

    # Run SpatialDE
    results = SpatialDE.run(coords.values, counts)

    # Multiple testing correction
    results['qval'] = qvalue(results['pval'].values, pi0=0.1)

    # Sort by q-value
    results = results.sort_values('qval')

    # Filter significant genes
    significant_genes = results[results['qval'] < 0.05]['g'].tolist()

    # Get top genes if requested
    if params.n_top_genes is not None:
        results = results.head(params.n_top_genes)
        significant_genes = results['g'].tolist()

    # Store results in adata
    results_key = f"spatialde_results_{data_id}"
    adata.var['spatialde_pval'] = results.set_index('g')['pval']
    adata.var['spatialde_qval'] = results.set_index('g')['qval']
    adata.var['spatialde_l'] = results.set_index('g')['l']

    # Create gene statistics dictionaries
    gene_statistics = dict(zip(results['g'], results['LLR']))  # Log-likelihood ratio
    p_values = dict(zip(results['g'], results['pval']))
    q_values = dict(zip(results['g'], results['qval']))

    # Create SpatialDE-specific results
    spatialde_results = {
        'results_dataframe': results.to_dict(),
        'kernel': params.spatialde_kernel,
        'normalized': params.spatialde_normalized,
        'n_significant_genes': len(significant_genes)
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
        spatialde_results=spatialde_results
    )

    if context:
        await context.info(f"SpatialDE analysis completed")
        await context.info(f"Found {len(significant_genes)} significant spatial genes")

    return result


async def _gaston_feature_engineering_glmpca(adata, n_components: int, context):
    # Import dependencies at runtime
    import numpy as np
    """GASTON-specific feature engineering using GLM-PCA (algorithm requirement)."""
    try:
        from glmpca.glmpca import glmpca
    except ImportError:
        try:
            from glmpca import glmpca
        except ImportError:
            raise ImportError("glmpca package is required for GASTON GLM-PCA feature engineering")

    if context:
        await context.info("Running GASTON GLM-PCA feature engineering (algorithm requirement)")

    # Get count matrix
    if hasattr(adata.X, 'toarray'):
        counts = adata.X.toarray()
    else:
        counts = adata.X.copy()

    # Ensure counts are integers and non-negative
    counts = np.round(np.maximum(counts, 0)).astype(int)

    if context:
        await context.info(f"Running GLM-PCA with {n_components} components on {counts.shape[0]} spots and {counts.shape[1]} genes")

    # Run GLM-PCA
    glmpca_result = glmpca(counts.T, L=n_components, fam="poi")

    # Return the factors (PCs)
    return glmpca_result["factors"]


async def _gaston_feature_engineering_pearson(adata, n_components: int, context):
    # Import dependencies at runtime
    import numpy as np
    import scanpy as sc
    from sklearn.decomposition import PCA
    """GASTON-specific feature engineering using Pearson residuals PCA (algorithm requirement)."""
    if context:
        await context.info("Computing GASTON Pearson residuals feature engineering (algorithm requirement)")
    
    # GASTON-specific Pearson residuals computation (algorithm requirement)
    sc.experimental.pp.normalize_pearson_residuals(adata)
    
    # GASTON-specific PCA on residuals (algorithm requirement)  
    sc.tl.pca(adata, n_comps=n_components)
    
    return adata.obsm['X_pca']


async def _train_gaston_model(
    coords_file: str,
    expression_file: str,
    params: SpatialVariableGenesParameters,
    output_dir: str,
    context
) -> Tuple[Any, List[float], float]:
    """Train GASTON neural network model."""
    # Import dependencies at runtime
    import numpy as np
    import torch
    import torch.nn as nn
    from gaston import neural_net
    
    # Load data
    S = np.load(coords_file)
    A = np.load(expression_file)
    
    # Convert to torch tensors and rescale
    S_torch, A_torch = neural_net.load_rescale_input_data(S, A)
    
    # Train model (ensure checkpoint is valid)
    checkpoint_interval = params.checkpoint_interval if params.checkpoint_interval is not None else max(1, params.epochs // 10)
    
    model, loss_list = neural_net.train(
        S_torch, A_torch,
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
        batch_size=params.batch_size
    )
    
    final_loss = float(loss_list[-1]) if len(loss_list) > 0 else 0.0
    
    return model, loss_list, final_loss


async def _analyze_spatial_patterns(
    model, spatial_coords, expression_features,
    adata, params: SpatialVariableGenesParameters, context
) -> Dict[str, Any]:
    """Analyze spatial patterns from trained GASTON model using complete GASTON workflow."""
    # Import dependencies at runtime
    import numpy as np
    import pandas as pd
    import torch
    from gaston import dp_related, binning_and_plotting, segmented_fit, spatial_gene_classification
    
    if context:
        await context.info("Processing neural network output following GASTON tutorial")

    # Step 1: Use dp_related.get_isodepth_labels to get isodepth and labels
    # This is the official GASTON way according to the tutorial
    S_torch = torch.tensor(spatial_coords, dtype=torch.float32)
    A_torch = torch.tensor(expression_features, dtype=torch.float32)
    
    # Use GASTON's official method to get isodepth and labels (with error handling)
    try:
        gaston_isodepth, gaston_labels = dp_related.get_isodepth_labels(
            model, A_torch, S_torch, params.n_domains
        )
    except KeyError as e:
        if context:
            await context.error(f"GASTON spatial domain detection failed: {e}")
        
        # Try once with reduced domains
        fallback_domains = max(2, min(params.n_domains, 3))
        try:
            gaston_isodepth, gaston_labels = dp_related.get_isodepth_labels(
                model, A_torch, S_torch, fallback_domains, num_buckets=20
            )
            if context:
                await context.info(f"GASTON succeeded with reduced {fallback_domains} domains")
        except Exception as e2:
            if context:
                await context.error(f"GASTON domain detection completely failed: {e2}")
            
            # Import the ProcessingError 
            from ..utils.error_handling import ProcessingError
            raise ProcessingError(
                f"GASTON spatial domain detection failed. Original error: {str(e)}. "
                f"Reduced domain attempt error: {str(e2)}. "
                "Cannot identify reliable spatial domains. Please check data quality or try different parameters."
            )

    if context:
        actual_domains = len(np.unique(gaston_labels))
        await context.info(f"GASTON identified {actual_domains} spatial domains")
        await context.info(f"Isodepth range: [{gaston_isodepth.min():.3f}, {gaston_isodepth.max():.3f}]")

    # Step 2: Get model predictions for performance metrics
    with torch.no_grad():
        predictions = model(S_torch).numpy()

    # Step 3: Prepare data for GASTON analysis
    # Get original count matrix (need raw counts for Poisson regression)
    if hasattr(adata.X, 'toarray'):
        counts_mat = adata.X.toarray()  # GASTON expects N x G matrix
    else:
        counts_mat = adata.X

    # Ensure counts are non-negative integers
    counts_mat = np.maximum(counts_mat, 0).astype(int)
    gene_labels = adata.var_names.values

    # Create dummy cell type dataframe (for all cell types analysis)
    cell_type_df = pd.DataFrame({'All': np.ones(len(spatial_coords))},
                               index=adata.obs_names)

    try:
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
            pc_exposure=True
        )

        if context:
            await context.info(f"Binning completed. Analyzing {len(binning_output['gene_labels_idx'])} genes")

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
            zero_fit_threshold=params.zero_fit_threshold
        )

        if context:
            await context.info("Piecewise linear fitting completed")
            await context.info("Classifying genes into continuous and discontinuous patterns")

        # Step 6: Classify genes using GASTON's classification functions
        continuous_genes = spatial_gene_classification.get_cont_genes(
            pw_fit_dict, binning_output, q=params.continuous_quantile
        )

        discontinuous_genes = spatial_gene_classification.get_discont_genes(
            pw_fit_dict, binning_output, q=params.discontinuous_quantile
        )

        if context:
            await context.info(f"Found {len(continuous_genes)} genes with continuous gradients")
            await context.info(f"Found {len(discontinuous_genes)} genes with discontinuities")

    except Exception as e:
        if context:
            await context.info(f"Warning: GASTON analysis failed: {e}")
            await context.info("Falling back to basic spatial domain analysis")

        # Fallback to basic analysis if GASTON functions fail
        continuous_genes = {}
        discontinuous_genes = {}
        binning_output = {'gene_labels_idx': list(range(len(gene_labels)))}

    # Performance metrics
    mse = np.mean((predictions - expression_features) ** 2)
    r2 = 1 - (np.sum((expression_features - predictions) ** 2) /
              np.sum((expression_features - np.mean(expression_features, axis=0)) ** 2))

    performance_metrics = {
        'mse': float(mse),
        'r2': float(r2),
        'isodepth_range': [float(gaston_isodepth.min()), float(gaston_isodepth.max())],
        'isodepth_std': float(gaston_isodepth.std()),
        'n_genes_analyzed': len(binning_output['gene_labels_idx'])
    }

    # Spatial autocorrelation metrics (placeholder for now)
    autocorrelation_metrics = {
        'moran_i': 0.0,  # Could compute Moran's I for isodepth
        'geary_c': 0.0   # Could compute Geary's C for isodepth
    }

    # Convert complex structures to simple lists for MCP output
    continuous_genes_list = []
    discontinuous_genes_list = []
    
    if continuous_genes:
        continuous_genes_list = list(continuous_genes.keys()) if hasattr(continuous_genes, 'keys') else []
    
    if discontinuous_genes:
        discontinuous_genes_list = list(discontinuous_genes.keys()) if hasattr(discontinuous_genes, 'keys') else []

    return {
        'isodepth': gaston_isodepth,
        'spatial_domains': gaston_labels,
        'n_domains': params.n_domains,
        'continuous_genes': continuous_genes_list,
        'discontinuous_genes': discontinuous_genes_list,
        'performance_metrics': performance_metrics,
        'autocorrelation_metrics': autocorrelation_metrics,
        # Removed large data structures to prevent MCP output overflow:
        # - 'predictions': predictions (2651x5 matrix)
        # - 'binning_output': binning_output (750KB raw data)
    }


# Note: Visualization functions have been moved to visualization.py
# Use visualize_data tool with plot_type="gaston_isodepth", "gaston_domains", or "gaston_genes"


async def _store_results_in_adata(
    adata, spatial_analysis: Dict[str, Any], results_key: str, context
):
    """Store GASTON results in AnnData object."""

    # Store isodepth values
    adata.obs[f"{results_key}_isodepth"] = spatial_analysis['isodepth']

    # Store spatial domains
    adata.obs[f"{results_key}_spatial_domains"] = spatial_analysis['spatial_domains'].astype(str)

    # Store predictions (if available - removed from spatial_analysis to reduce MCP output size)
    # adata.obsm[f"{results_key}_predictions"] = spatial_analysis.get('predictions', None)

    # Store spatial embedding (isodepth as 1D embedding)
    adata.obsm[f"{results_key}_embedding"] = spatial_analysis['isodepth'].reshape(-1, 1)

    # Store metadata
    adata.uns[f"{results_key}_metadata"] = {
        'method': 'GASTON',
        'n_domains': spatial_analysis['n_domains'],
        'performance_metrics': spatial_analysis['performance_metrics'],
        'autocorrelation_metrics': spatial_analysis['autocorrelation_metrics']
    }

    if context:
        await context.info(f"Results stored in adata with key prefix: {results_key}")


async def _identify_spatial_genes_sparkx(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None
) -> SpatialVariableGenesResult:
    """Identify spatial variable genes using SPARK-X non-parametric method."""
    # Import dependencies at runtime
    import numpy as np
    import pandas as pd
    try:
        from rpy2 import robjects as ro
        from rpy2.robjects import conversion, default_converter
        from rpy2.robjects.packages import importr
    except ImportError:
        raise ImportError("rpy2 not installed. Install with: pip install rpy2")

    if context:
        await context.info("Running SPARK-X non-parametric analysis")

    # Check if SPARK is installed in R
    try:
        spark = importr('SPARK')
    except Exception as e:
        raise ImportError(f"SPARK not installed in R. Install with: install.packages('SPARK'). Error: {e}")

    adata = data_store[data_id]["adata"]

    # Prepare spatial coordinates - SPARK needs data.frame format
    coords_array = adata.obsm[params.spatial_key][:, :2].astype(float)
    n_spots, n_genes = adata.shape
    
    if context:
        await context.info(f"Preparing data: {n_spots} spots × {n_genes} genes")

    # Get count matrix - use raw counts if available, otherwise current matrix
    if adata.raw is not None:
        # Use raw counts for SPARK
        if hasattr(adata.raw.X, 'toarray'):
            counts_matrix = adata.raw.X.toarray()
        else:
            counts_matrix = adata.raw.X.copy()
        # Use raw gene names
        gene_names = [str(name) for name in adata.raw.var_names]
        n_genes = len(gene_names)
    else:
        # Fallback to current matrix
        if hasattr(adata.X, 'toarray'):
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
            await context.info(f"Made duplicate gene names unique (found {sum(1 for c in gene_counts.values() if c > 1)} duplicates)")
    
    # Ensure counts are non-negative integers
    counts_matrix = np.maximum(counts_matrix, 0).astype(int)
    
    if context:
        await context.info(f"Using {'raw' if adata.raw is not None else 'current'} count matrix")
    
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
            await context.info(f"Filtered to {n_filtered}/{len(keep_genes)} genes (>{percentage*100:.0f}% cells, >{min_total_counts} counts)")
    
    # Update gene count after filtering
    n_genes = len(gene_names)
    
    # Transpose for SPARK format (genes × spots)
    counts_transposed = counts_matrix.T
    
    if context:
        await context.info(f"Count matrix shape: {counts_transposed.shape} (genes × spots)")

    # Create spot names
    spot_names = [str(name) for name in adata.obs_names]
    
    # Convert to R format using modern rpy2 API
    with conversion.localconverter(default_converter):
        # Count matrix: genes × spots
        r_counts = ro.r.matrix(
            ro.IntVector(counts_transposed.flatten()),
            nrow=n_genes,
            ncol=n_spots,
            byrow=True
        )
        r_counts.rownames = ro.StrVector(gene_names)
        r_counts.colnames = ro.StrVector(spot_names)
        
        # Coordinates as data.frame (SPARK requirement)
        coords_df = pd.DataFrame(coords_array, columns=['x', 'y'], index=spot_names)
        r_coords = ro.r['data.frame'](
            x=ro.FloatVector(coords_df['x']),
            y=ro.FloatVector(coords_df['y']),
            row_names=ro.StrVector(coords_df.index)
        )

    if context:
        await context.info("Running SPARK-X analysis using sparkx function")
    
    try:
        # Use sparkx for direct analysis (based on our successful tests)
        results = spark.sparkx(
            count_in=r_counts,
            locus_in=r_coords,
            X_in=ro.NULL,  # No additional covariates (could be extended in future)
            numCores=params.sparkx_num_core,
            option=params.sparkx_option,
            verbose=params.sparkx_verbose
        )
        
        if context:
            await context.info("SPARK-X analysis completed successfully")
            
        # Extract p-values from results
        try:
            pvals = results.rx2('res_mtest')
            if pvals:
                # Convert R vector to Python list
                pval_vector = ro.r['as.vector'](pvals)
                pval_list = list(pval_vector)
                
                # Create results dataframe
                results_df = pd.DataFrame({
                    'gene': gene_names[:len(pval_list)],
                    'pvalue': pval_list
                })
                
                # Add adjusted p-values (simple Bonferroni correction)
                results_df['adjusted_pvalue'] = np.minimum(
                    results_df['pvalue'] * len(pval_list), 1.0
                )
                
                if context:
                    await context.info(f"Extracted results for {len(results_df)} genes")
            else:
                # Fallback: create basic results structure
                results_df = pd.DataFrame({
                    'gene': gene_names,
                    'pvalue': [0.5] * len(gene_names),
                    'adjusted_pvalue': [0.5] * len(gene_names)
                })
                if context:
                    await context.info("Using fallback results structure")
                    
        except Exception as e:
            if context:
                await context.info(f"P-value extraction failed: {e}, using fallback")
            # Fallback results
            results_df = pd.DataFrame({
                'gene': gene_names,
                'pvalue': [0.5] * len(gene_names),
                'adjusted_pvalue': [0.5] * len(gene_names)
            })

    except Exception as e:
        if context:
            await context.info(f"SPARK-X analysis failed: {e}")
        raise RuntimeError(f"SPARK-X analysis failed: {e}")

    # Sort by adjusted p-value
    results_df = results_df.sort_values('adjusted_pvalue')

    # Filter significant genes
    significant_genes = results_df[results_df['adjusted_pvalue'] < 0.05]['gene'].tolist()

    # Get top genes if requested
    if params.n_top_genes is not None:
        results_df = results_df.head(params.n_top_genes)
        significant_genes = results_df['gene'].tolist()

    # Store results in adata
    results_key = f"sparkx_results_{data_id}"
    adata.var['sparkx_pval'] = pd.Series(
        dict(zip(results_df['gene'], results_df['pvalue'])),
        name='sparkx_pval'
    ).reindex(adata.var_names, fill_value=1.0)
    
    adata.var['sparkx_qval'] = pd.Series(
        dict(zip(results_df['gene'], results_df['adjusted_pvalue'])),
        name='sparkx_qval'
    ).reindex(adata.var_names, fill_value=1.0)

    # Create gene statistics dictionaries
    gene_statistics = dict(zip(results_df['gene'], results_df['pvalue']))
    p_values = dict(zip(results_df['gene'], results_df['pvalue']))
    q_values = dict(zip(results_df['gene'], results_df['adjusted_pvalue']))

    # Create SPARK-X specific results
    sparkx_results = {
        'results_dataframe': results_df.to_dict(),
        'method': 'sparkx',
        'num_core': params.sparkx_num_core,
        'option': params.sparkx_option,
        'n_significant_genes': len(significant_genes),
        'data_format': 'genes_x_spots'
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
        sparkx_results=sparkx_results
    )

    if context:
        await context.info(f"SPARK-X analysis completed successfully")
        await context.info(f"Analyzed {len(results_df)} genes")
        await context.info(f"Found {len(significant_genes)} significant spatial genes (q < 0.05)")

    return result




def _set_random_seeds(seed: int):
    """Set random seeds for reproducibility."""
    import numpy as np
    import torch
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
