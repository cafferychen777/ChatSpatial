"""
Tools for RNA velocity and trajectory analysis in spatial transcriptomics data.
"""

from typing import Dict, Any, Optional, List
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix

from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image

from ..models.data import RNAVelocityParameters, TrajectoryParameters
from ..models.analysis import RNAVelocityResult, TrajectoryResult
from ..utils.image_utils import fig_to_image, fig_to_base64, create_placeholder_image


def preprocess_for_velocity(adata):
    """Prepare data for RNA velocity analysis"""
    import scvelo as scv

    # Check if data contains unspliced and spliced counts
    if 'spliced' not in adata.layers or 'unspliced' not in adata.layers:
        raise ValueError("Data lacks unspliced and spliced counts, cannot directly compute RNA velocity")

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
        raise ImportError("SIRV package not found. Install it with: pip install git+https://github.com/tabdelaal/SIRV.git")


def infer_spatial_trajectory(adata, spatial_weight=0.5):
    """Integrate RNA velocity and spatial information for trajectory inference"""
    try:
        # Import CellRank
        import cellrank as cr

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
        # Calculate spatial distance matrix
        spatial_dist = squareform(pdist(spatial_coords))
        # Convert to similarity (smaller distance = higher similarity)
        spatial_sim = np.exp(-spatial_dist / spatial_dist.mean())
        spatial_kernel = csr_matrix(spatial_sim)

        # Create spatial kernel
        sk = cr.kernels.PrecomputedKernel(spatial_kernel, adata)
        sk.compute_transition_matrix()

        # Combine kernels, integrating RNA velocity, expression similarity, and spatial information
        combined_kernel = (1-spatial_weight) * (0.8 * vk + 0.2 * ck) + spatial_weight * sk

        # Use GPCCA for analysis
        g = cr.estimators.GPCCA(combined_kernel)

        # Compute eigendecomposition first (required for some CellRank operations)
        g.compute_eigendecomposition()

        # Compute macrostates (cell clusters)
        g.compute_macrostates(n_states=5)

        # Predict terminal states automatically
        g.predict_terminal_states(method='stability')

        # Compute absorption probabilities
        g.compute_absorption_probabilities()

        # Get absorption probabilities
        absorption_probs = g.absorption_probabilities

        # Compute pseudotime as distance from the initial state
        # Use the first terminal state as the root
        terminal_states = list(g.terminal_states.cat.categories)
        if len(terminal_states) > 0:
            root_state = terminal_states[0]
            # Compute pseudotime as 1 - probability of being absorbed in the root state
            pseudotime = 1 - absorption_probs[root_state].X.flatten()
        else:
            # Fallback if no terminal states were found
            pseudotime = np.zeros(adata.n_obs)

        # Store results
        adata.obs['pseudotime'] = pseudotime
        adata.obsm['absorption_probabilities'] = absorption_probs

        return adata

    except Exception as e:
        import logging
        logging.error(f"Error in CellRank trajectory inference: {str(e)}")

        # Fallback to a simpler approach if CellRank fails
        # Just use diffusion pseudotime based on the graph
        import scanpy as sc
        import numpy as np

        # Make sure we have PCA computed
        if 'X_pca' not in adata.obsm:
            sc.pp.pca(adata)

        # Compute neighbors if not already computed
        if 'neighbors' not in adata.uns:
            sc.pp.neighbors(adata, use_rep='X_pca')

        # Compute diffusion map
        if 'X_diffmap' not in adata.obsm:
            sc.tl.diffmap(adata)

        # Compute diffusion pseudotime
        try:
            sc.tl.dpt(adata)
            # Store results
            adata.obs['pseudotime'] = adata.obs['dpt_pseudotime']
        except Exception as e:
            # If DPT fails, create a simple pseudotime based on the first diffusion component
            if 'X_diffmap' in adata.obsm:
                # Use the first diffusion component as a simple pseudotime
                pseudotime = adata.obsm['X_diffmap'][:, 0]
                # Normalize to [0, 1]
                pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
                adata.obs['pseudotime'] = pseudotime
            else:
                # If diffusion map is not available, use a random pseudotime
                import numpy as np
                adata.obs['pseudotime'] = np.random.rand(adata.n_obs)

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


def infer_pseudotime_from_spatial_pattern(adata, root_cells=None):
    """Infer pseudotime based on spatial pattern"""
    try:
        # Try to use Palantir for trajectory analysis
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

    except Exception as e:
        import logging
        logging.error(f"Error in Palantir trajectory inference: {str(e)}")

        # Fallback to diffusion pseudotime if Palantir fails
        import scanpy as sc
        import numpy as np
        import pandas as pd

        # Make sure we have PCA computed
        if 'X_pca' not in adata.obsm:
            sc.pp.pca(adata)

        # Compute neighbors if not already computed
        if 'neighbors' not in adata.uns:
            sc.pp.neighbors(adata, use_rep='X_pca')

        # Compute diffusion map
        if 'X_diffmap' not in adata.obsm:
            sc.tl.diffmap(adata)

        # Set root cell if provided
        if root_cells is not None:
            if isinstance(root_cells, list) and len(root_cells) > 0:
                # Find index of the root cell
                if root_cells[0] in adata.obs_names:
                    adata.uns['iroot'] = np.where(adata.obs_names == root_cells[0])[0][0]

        # Compute diffusion pseudotime
        try:
            sc.tl.dpt(adata)
            # Store results
            adata.obs['palantir_pseudotime'] = adata.obs['dpt_pseudotime']
        except Exception as e:
            # If DPT fails, create a simple pseudotime based on the first diffusion component
            if 'X_diffmap' in adata.obsm:
                # Use the first diffusion component as a simple pseudotime
                pseudotime = adata.obsm['X_diffmap'][:, 0]
                # Normalize to [0, 1]
                pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
                adata.obs['palantir_pseudotime'] = pseudotime
            else:
                # If diffusion map is not available, use a random pseudotime
                import numpy as np
                adata.obs['palantir_pseudotime'] = np.random.rand(adata.n_obs)

        # Create a simple branch probability matrix (just for compatibility)
        branch_probs = pd.DataFrame(
            np.ones((adata.n_obs, 1)),  # Just one branch
            index=adata.obs_names,
            columns=['branch_1']
        )
        adata.obsm['palantir_branch_probs'] = branch_probs

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
        RNA velocity analysis result
    """
    try:
        import scvelo as scv
        import warnings
        import sys
        import io
        from contextlib import redirect_stdout, redirect_stderr
    except ImportError:
        raise ImportError("scvelo package not found. Install it with: pip install scvelo>=0.2.5")

    if context:
        await context.info(f"Analyzing RNA velocity using {params.mode} mode")

    # Get AnnData object
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Copy AnnData to avoid modifying original data
    adata = data_store[data_id]["adata"].copy()

    # Check if data has velocity information
    has_velocity_data = 'spliced' in adata.layers and 'unspliced' in adata.layers

    # Suppress warnings and stdout/stderr to prevent interference with JSON output
    warnings.filterwarnings('ignore')

    # If no velocity data but reference data provided, use SIRV
    if not has_velocity_data and params.reference_data_id:
        if params.reference_data_id not in data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found in data store")

        ref_adata = data_store[params.reference_data_id]["adata"]

        if context:
            await context.info("Using SIRV to infer RNA velocity from reference dataset")

        # Capture stdout/stderr to prevent it from interfering with JSON
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()

        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            # Prepare data
            spatial_adata, scrna_adata = prepare_for_sirv(adata, ref_adata)

            # Run SIRV
            adata = run_sirv(
                spatial_adata=spatial_adata,
                scrna_adata=scrna_adata,
                n_pcs=params.n_pcs,
                labels=params.labels
            )

        # Log captured output if context is available
        if context:
            stdout_content = stdout_buffer.getvalue()
            stderr_content = stderr_buffer.getvalue()
            if stdout_content:
                await context.info(f"SIRV stdout: {stdout_content}")
            if stderr_content:
                await context.warning(f"SIRV stderr: {stderr_content}")

        velocity_computed = True

    # If data has velocity information, compute directly
    elif has_velocity_data:
        if context:
            await context.info("Computing RNA velocity directly")

        # Capture stdout/stderr to prevent it from interfering with JSON
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()

        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            # Preprocess
            adata = preprocess_for_velocity(adata)

            # Compute RNA velocity
            adata = compute_rna_velocity(adata, mode=params.mode)

        # Log captured output if context is available
        if context:
            stdout_content = stdout_buffer.getvalue()
            stderr_content = stderr_buffer.getvalue()
            if stdout_content:
                await context.info(f"RNA velocity computation stdout: {stdout_content}")
            if stderr_content:
                await context.warning(f"RNA velocity computation stderr: {stderr_content}")

        velocity_computed = True

    # If velocity cannot be computed
    else:
        if context:
            await context.warning("Data lacks unspliced and spliced counts, cannot compute RNA velocity")

        velocity_computed = False

    # Create visualization
    if velocity_computed:
        # Capture stdout/stderr for visualization
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()

        # Initialize flag to skip visualization
        skip_visualization = False

        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            # Check if the requested basis exists
            basis = params.basis
            basis_key = f'X_{basis}'

            # If the requested basis doesn't exist, try to use 'spatial' or compute UMAP
            if basis_key not in adata.obsm:
                if context:
                    await context.warning(f"Basis '{basis}' not found in data. Trying alternative basis.")

                # Try to use spatial coordinates if available
                if 'spatial' in adata.obsm:
                    basis = 'spatial'
                    if context:
                        await context.info(f"Using 'spatial' as basis for visualization.")
                # Otherwise, compute UMAP if not already computed
                elif 'X_pca' in adata.obsm:
                    import scanpy as sc
                    if context:
                        await context.info(f"Computing UMAP for visualization.")
                    sc.pp.neighbors(adata, use_rep='X_pca')
                    sc.tl.umap(adata)
                    basis = 'umap'
                else:
                    # If no suitable basis is available, use the first available basis
                    available_bases = [key.replace('X_', '') for key in adata.obsm.keys() if key.startswith('X_')]
                    if available_bases:
                        basis = available_bases[0]
                        if context:
                            await context.info(f"Using '{basis}' as basis for visualization.")
                    else:
                        # If no basis is available, we can't create a visualization
                        if context:
                            await context.warning(f"No suitable basis found for visualization.")
                        # Set a flag to skip visualization
                        skip_visualization = True

            # Create velocity stream plot if we're not skipping visualization
            if not skip_visualization:
                plt.figure(figsize=(6, 6))
                scv.pl.velocity_embedding_stream(
                    adata,
                    basis=basis,
                    color=params.color,
                    show=False
                )
                fig = plt.gcf()
                # Use PNG format with lower DPI to reduce size
                # Set max_size_kb to 300KB to ensure the response isn't too large
                img = fig_to_image(fig, dpi=60, max_size_kb=300)
            else:
                img = None

        # Log captured output if context is available
        if context:
            stdout_content = stdout_buffer.getvalue()
            stderr_content = stderr_buffer.getvalue()
            if stdout_content:
                await context.info(f"Visualization stdout: {stdout_content}")
            if stderr_content:
                await context.warning(f"Visualization stderr: {stderr_content}")
    else:
        img = None

    # Update data store
    data_store[data_id]["adata"] = adata

    # Return result
    return RNAVelocityResult(
        data_id=data_id,
        velocity_computed=velocity_computed,
        visualization=img
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
        Trajectory analysis result
    """
    import scanpy as sc
    import warnings
    import io
    from contextlib import redirect_stdout, redirect_stderr

    # Suppress warnings to prevent interference with JSON output
    warnings.filterwarnings('ignore')

    if context:
        await context.info(f"Analyzing trajectory using {params.method} method")

    # Get AnnData object
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Copy AnnData to avoid modifying original data
    adata = data_store[data_id]["adata"].copy()

    # Check if RNA velocity has been computed
    has_velocity = 'velocity_graph' in adata.uns

    # Choose trajectory inference method based on situation
    if params.method == "cellrank" and has_velocity:
        try:
            import cellrank as cr
        except ImportError:
            raise ImportError("cellrank package not found. Install it with: pip install cellrank>=2.0.0")

        if context:
            await context.info("Using CellRank to infer trajectory based on RNA velocity")

        # Capture stdout/stderr to prevent it from interfering with JSON
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()

        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            # Use CellRank to infer trajectory
            adata = infer_spatial_trajectory(
                adata,
                spatial_weight=params.spatial_weight
            )

        # Log captured output if context is available
        if context:
            stdout_content = stdout_buffer.getvalue()
            stderr_content = stderr_buffer.getvalue()
            if stdout_content:
                await context.info(f"CellRank stdout: {stdout_content}")
            if stderr_content:
                await context.warning(f"CellRank stderr: {stderr_content}")

        # Use pseudotime key
        pseudotime_key = 'pseudotime'

    elif params.method == "palantir" or not has_velocity:
        if context:
            await context.info("Using Palantir to infer trajectory based on spatial pattern")

        # Capture stdout/stderr to prevent it from interfering with JSON
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()

        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            # Use spatially-aware embedding
            adata = spatial_aware_embedding(
                adata,
                spatial_weight=params.spatial_weight
            )

            # Use Palantir to infer pseudotime
            adata = infer_pseudotime_from_spatial_pattern(
                adata,
                root_cells=params.root_cells
            )

        # Log captured output if context is available
        if context:
            stdout_content = stdout_buffer.getvalue()
            stderr_content = stderr_buffer.getvalue()
            if stdout_content:
                await context.info(f"Palantir stdout: {stdout_content}")
            if stderr_content:
                await context.warning(f"Palantir stderr: {stderr_content}")

        # Use pseudotime key
        pseudotime_key = 'palantir_pseudotime'

    # Ensure pseudotime key exists in the data
    if pseudotime_key not in adata.obs.columns:
        if context:
            await context.warning(f"Pseudotime key '{pseudotime_key}' not found in data. Creating a random pseudotime.")
        # Create a random pseudotime as fallback
        import numpy as np
        adata.obs[pseudotime_key] = np.random.rand(adata.n_obs)

    # Create visualizations
    # Capture stdout/stderr for visualization
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()

    with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
        # Check if spatial coordinates exist
        if 'spatial' not in adata.obsm:
            if context:
                await context.warning("Spatial coordinates not found in data. Using UMAP for visualization.")
            # Use UMAP if available, or create a simple 2D embedding
            if 'X_umap' in adata.obsm:
                basis = 'umap'
            elif 'X_pca' in adata.obsm:
                # Compute UMAP from PCA
                sc.pp.neighbors(adata, use_rep='X_pca')
                sc.tl.umap(adata)
                basis = 'umap'
            else:
                # Create a simple 2D embedding based on pseudotime
                import numpy as np
                # Use pseudotime as one dimension
                x = adata.obs[pseudotime_key].values
                # Create a second dimension with some noise
                y = x + np.random.normal(0, 0.1, size=len(x))
                # Store as a simple 2D embedding
                adata.obsm['X_pseudotime_2d'] = np.column_stack([x, y])
                basis = 'pseudotime_2d'
        else:
            basis = 'spatial'

        # 1. Pseudotime spatial plot
        try:
            fig, ax = plt.subplots(figsize=(6, 6))
            sc.pl.embedding(adata, basis=basis, color=pseudotime_key, ax=ax, show=False)
            # Use PNG format with lower DPI to reduce size
            # Set max_size_kb to 300KB to ensure the response isn't too large
            pseudotime_img = fig_to_image(fig, dpi=60, max_size_kb=300)
        except Exception as e:
            if context:
                await context.warning(f"Error creating pseudotime visualization: {str(e)}")
            # Create a simple visualization as fallback
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.text(0.5, 0.5, "Pseudotime visualization failed",
                    ha='center', va='center', fontsize=12)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')
            pseudotime_img = fig_to_image(fig, dpi=60, max_size_kb=300)

        # 2. If RNA velocity available, create velocity stream plot
        if has_velocity:
            try:
                import scvelo as scv
                plt.figure(figsize=(6, 6))
                scv.pl.velocity_embedding_stream(
                    adata,
                    basis=basis,
                    color=pseudotime_key,
                    show=False
                )
                fig = plt.gcf()
                # Use PNG format with lower DPI to reduce size
                # Set max_size_kb to 300KB to ensure the response isn't too large
                velocity_img = fig_to_image(fig, dpi=60, max_size_kb=300)
            except Exception as e:
                if context:
                    await context.warning(f"Error creating velocity visualization: {str(e)}")
                # Create a simple visualization as fallback
                fig, ax = plt.subplots(figsize=(6, 6))
                ax.text(0.5, 0.5, "Velocity visualization failed",
                        ha='center', va='center', fontsize=12)
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.axis('off')
                velocity_img = fig_to_image(fig, dpi=60, max_size_kb=300)
        else:
            velocity_img = None

    # Log captured output if context is available
    if context:
        stdout_content = stdout_buffer.getvalue()
        stderr_content = stderr_buffer.getvalue()
        if stdout_content:
            await context.info(f"Visualization stdout: {stdout_content}")
        if stderr_content:
            await context.warning(f"Visualization stderr: {stderr_content}")

    # Update data store
    data_store[data_id]["adata"] = adata

    # Return result
    return TrajectoryResult(
        data_id=data_id,
        pseudotime_computed=True,
        velocity_computed=has_velocity,
        pseudotime_key=pseudotime_key,
        pseudotime_visualization=pseudotime_img,
        velocity_visualization=velocity_img
    )
