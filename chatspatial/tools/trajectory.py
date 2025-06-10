"""
Tools for RNA velocity and trajectory analysis in spatial transcriptomics data.
"""

from typing import Dict, Any, Optional, List, Tuple
import numpy as np
import pandas as pd

from mcp.server.fastmcp import Context

from ..models.data import RNAVelocityParameters, TrajectoryParameters
from ..models.analysis import RNAVelocityResult, TrajectoryResult
from ..utils.output_utils import suppress_output, ProcessingError


def preprocess_for_velocity(adata):
    """Prepare data for RNA velocity analysis"""
    import scvelo as scv

    # Check if data contains unspliced and spliced counts
    if 'spliced' not in adata.layers or 'unspliced' not in adata.layers:
        raise ValueError("Data lacks unspliced and spliced counts")

    # Standard preprocessing
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)

    return adata


def compute_rna_velocity(adata, mode="stochastic"):
    """Compute RNA velocity"""
    import scvelo as scv

    # Compute RNA velocity based on selected mode
    if mode == "deterministic":
        scv.tl.velocity(adata)
    elif mode == "stochastic":
        scv.tl.velocity(adata, mode="stochastic")
    elif mode == "dynamical":
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode="dynamical")
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    # Compute velocity graph
    scv.tl.velocity_graph(adata)

    return adata


def prepare_for_sirv(spatial_adata, scrna_adata):
    """Prepare spatial data and scRNA-seq data for SIRV analysis"""
    import scvelo as scv

    # Ensure scRNA-seq data contains unspliced and spliced counts
    if 'spliced' not in scrna_adata.layers or 'unspliced' not in scrna_adata.layers:
        raise ValueError("scRNA-seq data lacks unspliced and spliced counts")

    # Ensure spatial data contains spatial coordinates
    if 'spatial' not in spatial_adata.obsm:
        raise ValueError("Spatial data lacks spatial coordinates")

    # Preprocess scRNA-seq data
    scv.pp.filter_and_normalize(scrna_adata)
    scv.pp.moments(scrna_adata)

    # Compute RNA velocity for scRNA-seq data
    scv.tl.velocity(scrna_adata)
    scv.tl.velocity_graph(scrna_adata)

    return spatial_adata, scrna_adata


def run_sirv(spatial_adata, scrna_adata, n_pcs=30, labels=None):
    """Run SIRV analysis"""
    try:
        from SIRV.main import SIRV

        # Run SIRV
        spatial_adata_with_velocity = SIRV(
            spatial_data=spatial_adata,
            scrna_data=scrna_adata,
            n_pcs=n_pcs,
            labels=labels
        )

        return spatial_adata_with_velocity
    except ImportError:
        raise ProcessingError("SIRV package not found. Install it with: pip install git+https://github.com/tabdelaal/SIRV.git")


def infer_spatial_trajectory_cellrank(adata, spatial_weight=0.5, kernel_weights=(0.8, 0.2), n_states=5):
    """Integrate RNA velocity and spatial information for trajectory inference using CellRank.
    
    Raises exceptions directly if CellRank fails - no fallback logic.
    """
    import cellrank as cr
    from scipy.spatial.distance import pdist, squareform
    from scipy.sparse import csr_matrix

    # Get spatial coordinates
    if 'spatial' not in adata.obsm:
        raise ValueError("Data lacks spatial coordinate information")

    spatial_coords = adata.obsm['spatial']

    # Create RNA velocity kernel
    vk = cr.kernels.VelocityKernel(adata)
    vk.compute_transition_matrix()

    # Create expression similarity-based kernel
    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()

    # Create custom spatial kernel
    spatial_dist = squareform(pdist(spatial_coords))
    spatial_sim = np.exp(-spatial_dist / spatial_dist.mean())
    spatial_kernel = csr_matrix(spatial_sim)

    # Create spatial kernel
    sk = cr.kernels.PrecomputedKernel(spatial_kernel, adata)
    sk.compute_transition_matrix()

    # Combine kernels using configurable weights
    vk_weight, ck_weight = kernel_weights
    combined_kernel = (1 - spatial_weight) * (vk_weight * vk + ck_weight * ck) + spatial_weight * sk

    # Use GPCCA for analysis
    g = cr.estimators.GPCCA(combined_kernel)

    # Compute eigendecomposition
    g.compute_eigendecomposition()

    # Compute macrostates
    g.compute_macrostates(n_states=n_states)

    # Predict terminal states
    g.predict_terminal_states(method='stability')

    # Compute absorption probabilities
    g.compute_absorption_probabilities()

    # Get absorption probabilities
    absorption_probs = g.absorption_probabilities

    # Compute pseudotime
    terminal_states = list(g.terminal_states.cat.categories)
    if not terminal_states:
        raise RuntimeError("CellRank could not identify any terminal states")
        
    root_state = terminal_states[0]
    pseudotime = 1 - absorption_probs[root_state].X.flatten()

    # Store results
    adata.obs['pseudotime'] = pseudotime
    adata.obsm['absorption_probabilities'] = absorption_probs

    return adata


def spatial_aware_embedding(adata, spatial_weight=0.3):
    """Generate spatially-aware low-dimensional embedding"""
    from sklearn.metrics.pairwise import euclidean_distances
    from umap import UMAP

    # Get spatial coordinates
    if 'spatial' not in adata.obsm:
        raise ValueError("Data lacks spatial coordinates")

    spatial_coords = adata.obsm['spatial']

    # Ensure PCA has been computed
    if 'X_pca' not in adata.obsm:
        raise ValueError("PCA has not been computed. Run preprocessing first.")

    # Calculate expression-based distance matrix
    expr_dist = euclidean_distances(adata.obsm['X_pca'])

    # Calculate spatial distance matrix
    spatial_dist = euclidean_distances(spatial_coords)

    # Combine distance matrices
    combined_dist = (1 - spatial_weight) * expr_dist + spatial_weight * spatial_dist

    # Use UMAP for dimensionality reduction
    umap_op = UMAP(metric='precomputed')
    embedding = umap_op.fit_transform(combined_dist)

    # Store results
    adata.obsm['X_spatial_umap'] = embedding

    return adata


def infer_pseudotime_palantir(adata, root_cells=None):
    """Infer pseudotime based on spatial pattern using Palantir.
    
    Raises exceptions directly if Palantir fails - no fallback logic.
    """
    import palantir
    import pandas as pd
    import numpy as np

    # Make sure we have PCA computed
    if 'X_pca' not in adata.obsm:
        import scanpy as sc
        sc.pp.pca(adata)

    # Convert PCA projections to DataFrame for Palantir
    pca_df = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)

    # Run diffusion maps
    dm_res = palantir.utils.run_diffusion_maps(pca_df, n_components=10)

    # Determine multiscale space
    ms_data = pd.DataFrame(
        dm_res['EigenVectors'],
        index=pca_df.index
    )

    # Determine early cell (start cell)
    if root_cells is not None and len(root_cells) > 0:
        if root_cells[0] not in ms_data.index:
            raise ValueError(f"Root cell '{root_cells[0]}' not found in data")
        start_cell = root_cells[0]
    else:
        # Use the cell at the extreme of the first diffusion component
        start_cell = ms_data.iloc[:, 0].idxmax()

    # Run Palantir
    pr_res = palantir.core.run_palantir(ms_data, start_cell)

    # Store pseudotime and branch probabilities
    adata.obs['palantir_pseudotime'] = pr_res.pseudotime
    adata.obsm['palantir_branch_probs'] = pr_res.branch_probs

    return adata


def compute_dpt_fallback(adata, root_cells=None):
    """Compute Diffusion Pseudotime as fallback method."""
    import scanpy as sc
    import numpy as np

    # Make sure we have PCA computed
    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata)

    # Compute neighbors if not already computed
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, use_rep='X_pca')

    # Compute diffusion map
    sc.tl.diffmap(adata)

    # Set root cell if provided
    if root_cells is not None and len(root_cells) > 0:
        if root_cells[0] in adata.obs_names:
            adata.uns['iroot'] = np.where(adata.obs_names == root_cells[0])[0][0]

    # Compute diffusion pseudotime
    sc.tl.dpt(adata)
    
    # Handle any NaN values
    adata.obs['dpt_pseudotime'] = adata.obs['dpt_pseudotime'].fillna(0)

    return adata


async def analyze_rna_velocity(
    data_id: str,
    data_store: Dict[str, Any],
    params: RNAVelocityParameters = RNAVelocityParameters(),
    context: Optional[Context] = None
) -> RNAVelocityResult:
    """Analyze RNA velocity in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing datasets
        params: RNA velocity analysis parameters
        context: MCP context

    Returns:
        RNA velocity analysis result (computation metadata only, no visualization)
    """
    try:
        import scvelo as scv
    except ImportError:
        raise ProcessingError("scvelo package not found. Install it with: pip install scvelo>=0.2.5")

    if context:
        await context.info(f"Analyzing RNA velocity using {params.mode} mode")

    # Get AnnData object
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Copy AnnData to avoid modifying original data
    adata = data_store[data_id]["adata"].copy()

    # Check if data has velocity information
    has_velocity_data = 'spliced' in adata.layers and 'unspliced' in adata.layers

    velocity_computed = False

    with suppress_output():
        try:
            if has_velocity_data:
                if context:
                    await context.info("Found spliced/unspliced layers. Computing velocity directly.")
                
                # Preprocess and compute velocity
                adata = preprocess_for_velocity(adata)
                adata = compute_rna_velocity(adata, mode=params.mode)
                velocity_computed = True
                
            elif params.reference_data_id:
                if context:
                    await context.info(f"No layers found. Using reference data '{params.reference_data_id}' for SIRV.")
                
                if params.reference_data_id not in data_store:
                    raise ValueError(f"Reference dataset {params.reference_data_id} not found in data store")
                
                ref_adata = data_store[params.reference_data_id]["adata"]
                
                # Prepare data for SIRV
                spatial_adata, scrna_adata = prepare_for_sirv(adata, ref_adata)
                
                # Run SIRV
                adata = run_sirv(
                    spatial_adata=spatial_adata,
                    scrna_adata=scrna_adata,
                    n_pcs=params.n_pcs,
                    labels=params.labels
                )
                velocity_computed = True
                
            else:
                if context:
                    await context.warning("Data lacks spliced/unspliced counts and no reference data was provided. Cannot compute RNA velocity.")

        except Exception as e:
            raise ProcessingError(f"RNA velocity analysis failed: {str(e)}") from e

    # Update data store
    data_store[data_id]["adata"] = adata

    # Return result with metadata only (no visualization)
    return RNAVelocityResult(
        data_id=data_id,
        velocity_computed=velocity_computed,
        velocity_graph_key='velocity_graph' if velocity_computed else None,
        mode=params.mode
    )


async def analyze_trajectory(
    data_id: str,
    data_store: Dict[str, Any],
    params: TrajectoryParameters = TrajectoryParameters(),
    context: Optional[Context] = None
) -> TrajectoryResult:
    """Analyze trajectory and cell state transitions in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing datasets
        params: Trajectory analysis parameters
        context: MCP context

    Returns:
        Trajectory analysis result (computation metadata only, no visualization)
    """
    if context:
        await context.info(f"Analyzing trajectory using {params.method} method")

    # Get AnnData object
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Copy AnnData to avoid modifying original data
    adata = data_store[data_id]["adata"].copy()

    # Check if RNA velocity has been computed
    has_velocity = 'velocity_graph' in adata.uns

    pseudotime_key = None
    method_used = params.method

    # Strategy 1: Try CellRank if requested and velocity is available
    if params.method == "cellrank" and has_velocity:
        if context:
            await context.info("Attempting trajectory inference with CellRank...")
        
        try:
            import cellrank as cr
        except ImportError:
            raise ProcessingError("cellrank package not found. Install it with: pip install cellrank>=2.0.0")

        try:
            with suppress_output():
                adata = infer_spatial_trajectory_cellrank(
                    adata,
                    spatial_weight=params.spatial_weight,
                    kernel_weights=params.cellrank_kernel_weights,
                    n_states=params.cellrank_n_states
                )
            pseudotime_key = 'pseudotime'
            if context:
                await context.info("CellRank analysis completed successfully.")
        except Exception as cellrank_error:
            if context:
                await context.warning(f"CellRank analysis failed: {cellrank_error}")
                await context.info("Falling back to Palantir method...")
            # Force switch to Palantir
            method_used = "palantir"

    # Strategy 2: Try Palantir (either as primary method or fallback)
    if method_used == "palantir" or not has_velocity:
        if context:
            await context.info("Attempting trajectory inference with Palantir...")
        
        try:
            with suppress_output():
                # Run spatially-aware embedding
                adata = spatial_aware_embedding(adata, spatial_weight=params.spatial_weight)
                
                # Run Palantir
                adata = infer_pseudotime_palantir(adata, root_cells=params.root_cells)
                
            pseudotime_key = 'palantir_pseudotime'
            method_used = "palantir"
            if context:
                await context.info("Palantir analysis completed successfully.")
                
        except Exception as palantir_error:
            if not params.allow_fallback_to_dpt:
                raise ProcessingError(f"Palantir analysis failed: {palantir_error}") from palantir_error
            
            if context:
                await context.warning(f"Palantir analysis failed: {palantir_error}")
                await context.info("Falling back to Diffusion Pseudotime (DPT)...")
            
            # Final fallback: DPT
            try:
                with suppress_output():
                    adata = compute_dpt_fallback(adata, root_cells=params.root_cells)
                    
                pseudotime_key = 'dpt_pseudotime'
                method_used = "dpt"
                if context:
                    await context.info("DPT fallback analysis completed successfully.")
                    
            except Exception as dpt_error:
                raise ProcessingError(f"All trajectory inference methods failed. Last error (DPT): {dpt_error}") from dpt_error

    # Ensure pseudotime key exists
    if pseudotime_key is None or pseudotime_key not in adata.obs.columns:
        raise ProcessingError("Failed to compute pseudotime with any available method")

    # Update data store
    data_store[data_id]["adata"] = adata

    # Return result with metadata only (no visualization)
    return TrajectoryResult(
        data_id=data_id,
        pseudotime_computed=True,
        velocity_computed=has_velocity,
        pseudotime_key=pseudotime_key,
        method=method_used,
        spatial_weight=params.spatial_weight
    )
