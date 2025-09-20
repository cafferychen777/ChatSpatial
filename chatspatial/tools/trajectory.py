"""
A module for RNA velocity and trajectory inference in spatial transcriptomics.

This module provides a suite of tools for analyzing cellular dynamics by
integrating single-cell RNA velocity with spatial information. It includes
functions for preprocessing data with spliced and unspliced counts, computing
RNA velocity using various models (e.g., stochastic, dynamical), and inferring
cellular trajectories.

Key functionalities are organized into two main analysis pipelines:
1. `analyze_rna_velocity`: Computes RNA velocity directly from spatial data
   that contains both spliced and unspliced count layers.
2. `analyze_trajectory`: Infers developmental trajectories and pseudotime by
   combining velocity information with spatial context. It provides access to
   advanced methods like CellRank and Palantir, with a fallback to Diffusion
   Pseudotime (DPT) for robustness.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from mcp.server.fastmcp import Context

from ..models.analysis import RNAVelocityResult, TrajectoryResult
from ..models.data import RNAVelocityParameters, TrajectoryParameters
from ..utils.error_handling import (DataNotFoundError, ProcessingError,
                                    suppress_output, validate_adata)

logger = logging.getLogger(__name__)

# Import scvi-tools for advanced trajectory analysis
try:
    import scvi
    from scvi.external import VELOVI
except ImportError:
    scvi = None
    VELOVI = None


def validate_velocity_data(adata) -> Tuple[bool, List[str]]:
    """
    Validate if AnnData object has necessary data for RNA velocity analysis.

    This function now uses the unified validation system for consistency,
    but maintains the same interface for backward compatibility.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix to validate

    Returns
    -------
    Tuple[bool, List[str]]
        Tuple of (is_valid, list_of_issues)
    """
    try:
        validate_adata(adata, {}, check_velocity=True)
        return True, []
    except DataNotFoundError as e:
        # Extract individual issues from the error message
        error_msg = str(e)
        if "Validation failed: " in error_msg:
            issues = error_msg.replace("Validation failed: ", "").split(", ")
        else:
            issues = [error_msg]
        return False, issues


def validate_spatial_data(
    adata, spatial_key: str = "spatial"
) -> Tuple[bool, List[str]]:
    """
    Validate if AnnData object has necessary spatial information.

    This function now uses the unified validation system for consistency,
    but maintains the same interface for backward compatibility.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix to validate
    spatial_key : str, default 'spatial'
        Key for spatial coordinates in adata.obsm

    Returns
    -------
    Tuple[bool, List[str]]
        Tuple of (is_valid, list_of_issues)
    """
    try:
        validate_adata(adata, {}, check_spatial=True, spatial_key=spatial_key)
        return True, []
    except DataNotFoundError as e:
        # Extract individual issues from the error message
        error_msg = str(e)
        if "Validation failed: " in error_msg:
            issues = error_msg.replace("Validation failed: ", "").split(", ")
        else:
            issues = [error_msg]
        return False, issues


def preprocess_for_velocity(
    adata, min_shared_counts=30, n_top_genes=2000, n_pcs=30, n_neighbors=30, params=None
):
    """
    Prepares an AnnData object for RNA velocity analysis using the scVelo pipeline.

    This function performs the standard scVelo preprocessing workflow, which is a
    prerequisite for computing RNA velocity. The steps include:
    1. Filtering genes based on minimum shared counts between spliced and
       unspliced layers.
    2. Normalizing the data.
    3. Selecting a subset of highly variable genes.
    4. Computing first and second-order moments (mean and uncentered variance)
       across nearest neighbors in PCA space.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix, which must contain 'spliced' and 'unspliced' layers.
    min_shared_counts : int, default 30
        Minimum number of counts shared between spliced and unspliced layers.
    n_top_genes : int, default 2000
        Number of highly variable genes to use for analysis.
    n_pcs : int, default 30
        Number of principal components to compute.
    n_neighbors : int, default 30
        Number of nearest neighbors for moment computation.
    params : RNAVelocityParameters, optional
        If provided, this object's attributes will override the individual parameters.
    """
    import scvelo as scv

    # If params object is provided, use its values
    if params is not None:
        from ..models.data import RNAVelocityParameters

        if isinstance(params, RNAVelocityParameters):
            min_shared_counts = params.min_shared_counts
            n_top_genes = params.n_top_genes
            n_pcs = params.n_pcs
            n_neighbors = params.n_neighbors

    # Validate velocity data
    is_valid, issues = validate_velocity_data(adata)
    if not is_valid:
        raise ValueError(f"Invalid velocity data: {'; '.join(issues)}")

    # Standard preprocessing with configurable parameters
    scv.pp.filter_and_normalize(
        adata, min_shared_counts=min_shared_counts, n_top_genes=n_top_genes
    )
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    return adata


def compute_rna_velocity(adata, mode="stochastic", params=None):
    """
    Computes RNA velocity to infer the direction of cellular differentiation.

    This function executes the core RNA velocity workflow. It first ensures that
    the necessary preprocessing steps (e.g., moment computation) have been performed.
    It then estimates RNA velocity, which represents the rate of change of gene
    expression, by modeling the balance of spliced and unspliced mRNA counts.
    Finally, it constructs a velocity graph to represent the inferred cell-to-cell
    transitions.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix with 'spliced' and 'unspliced' layers.
    mode : str, default 'stochastic'
        The model for velocity estimation.
        - 'stochastic': A likelihood-based model that accounts for noise.
        - 'deterministic': A simpler model based on steady-state assumptions.
        - 'dynamical': A model that solves the full transcriptional dynamics.
    params : RNAVelocityParameters, optional
        An object containing parameters for both preprocessing and velocity computation.

    Returns
    -------
    AnnData
        The AnnData object updated with the computed velocity vectors and graph.
    """
    import scvelo as scv

    # Use params for mode if provided
    if params is not None:
        from ..models.data import RNAVelocityParameters

        if isinstance(params, RNAVelocityParameters):
            mode = params.mode

    # Check if preprocessing is needed
    if "Ms" not in adata.layers or "Mu" not in adata.layers:
        # Run preprocessing
        adata = preprocess_for_velocity(adata, params=params)

    # Compute velocity based on mode
    if mode == "dynamical":
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode="dynamical")
    else:
        scv.tl.velocity(adata, mode=mode)

    # Compute velocity graph
    scv.tl.velocity_graph(adata)

    return adata


def validate_rna_velocity_computation(adata, mode="stochastic"):
    """Validate RNA velocity computation requirements (no computation performed)"""

    # Check if velocity has been computed
    velocity_key = "velocity"
    if velocity_key not in adata.layers:
        raise ValueError(
            f"RNA velocity ({mode} mode) not found but required for analysis. "
            "Please run velocity computation in preprocessing.py: "
            f"scv.tl.velocity(adata, mode='{mode}')"
        )

    # Check if velocity graph has been computed
    if "velocity_graph" not in adata.uns:
        raise ValueError(
            "Velocity graph not found but required for trajectory analysis. "
            "Please run in preprocessing.py: scv.tl.velocity_graph(adata)"
        )

    return adata


def infer_spatial_trajectory_cellrank(
    adata, spatial_weight=0.5, kernel_weights=(0.8, 0.2), n_states=5
):
    """
    Infers cellular trajectories by combining RNA velocity and spatial data with CellRank.

    This function uses CellRank to model cell-state transitions by constructing
    a transition matrix from multiple kernels. It combines:
    1. A velocity kernel, derived from RNA velocity.
    2. A connectivity kernel, based on transcriptomic similarity.
    3. A custom spatial kernel, based on physical proximity.

    By analyzing the eigenvectors of the combined transition matrix, CellRank
    identifies macrostates (representative cell states), terminal states, and computes
    fate probabilities to map the paths of cellular development.

    Raises exceptions directly if CellRank fails - no fallback logic.
    """
    import cellrank as cr
    from scipy.sparse import csr_matrix
    from scipy.spatial.distance import pdist, squareform

    # Validate spatial data
    is_valid, issues = validate_spatial_data(adata)
    if not is_valid:
        raise ValueError(f"Invalid spatial data: {'; '.join(issues)}")

    spatial_coords = adata.obsm["spatial"]

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
    combined_kernel = (1 - spatial_weight) * (
        vk_weight * vk + ck_weight * ck
    ) + spatial_weight * sk

    # Use GPCCA for analysis
    g = cr.estimators.GPCCA(combined_kernel)

    # Compute eigendecomposition
    g.compute_eigendecomposition()

    # Compute macrostates
    try:
        g.compute_macrostates(n_states=n_states)
    except Exception as e:
        # If automatic n_states fails, try with fewer states
        for alt_n_states in [n_states - 1, n_states - 2, 3, 2]:
            if alt_n_states < 2:
                break
            try:
                g.compute_macrostates(n_states=alt_n_states)
                break
            except:
                continue
        else:
            raise RuntimeError(
                f"Failed to compute macrostates with any number of states: {e}"
            )

    # Predict terminal states
    try:
        g.predict_terminal_states(method="stability")
    except ValueError as e:
        if "No macrostates have been selected" in str(e):
            # Skip terminal state prediction if no macrostates
            pass
        else:
            raise

    # Check if we have terminal states
    has_terminal_states = (
        hasattr(g, "terminal_states") and g.terminal_states is not None
    )

    if has_terminal_states and len(g.terminal_states.cat.categories) > 0:
        # Compute fate probabilities (renamed from absorption_probabilities in newer versions)
        g.compute_fate_probabilities()

        # Get fate probabilities (renamed from absorption_probabilities)
        absorption_probs = g.fate_probabilities

        # Compute pseudotime
        terminal_states = list(g.terminal_states.cat.categories)
        root_state = terminal_states[0]
        pseudotime = 1 - absorption_probs[root_state].X.flatten()

        # Store results
        adata.obs["pseudotime"] = pseudotime
        adata.obsm["fate_probabilities"] = absorption_probs
        adata.obs["terminal_states"] = g.terminal_states
    else:
        # No terminal states, use macrostates for pseudotime
        if hasattr(g, "macrostates") and g.macrostates is not None:
            # Use the macrostate memberships as a proxy for pseudotime
            macrostate_probs = g.macrostates_memberships
            # Use the first macrostate as the "early" state
            pseudotime = 1 - macrostate_probs[:, 0].X.flatten()
            adata.obs["pseudotime"] = pseudotime
        else:
            raise RuntimeError(
                "CellRank could not compute either terminal states or macrostates"
            )

    # Always store macrostates if available
    if hasattr(g, "macrostates") and g.macrostates is not None:
        adata.obs["macrostates"] = g.macrostates

    return adata


def spatial_aware_embedding(adata, spatial_weight=0.3):
    """Generate spatially-aware low-dimensional embedding"""
    from sklearn.metrics.pairwise import euclidean_distances
    from umap import UMAP

    # Validate spatial data
    is_valid, issues = validate_spatial_data(adata)
    if not is_valid:
        raise ValueError(f"Invalid spatial data: {'; '.join(issues)}")

    spatial_coords = adata.obsm["spatial"]

    # Ensure PCA has been computed
    if "X_pca" not in adata.obsm:
        raise ValueError("PCA has not been computed. Run preprocessing first.")

    # Calculate expression-based distance matrix
    expr_dist = euclidean_distances(adata.obsm["X_pca"])

    # Calculate spatial distance matrix
    spatial_dist = euclidean_distances(spatial_coords)

    # Combine distance matrices
    combined_dist = (1 - spatial_weight) * expr_dist + spatial_weight * spatial_dist

    # Use UMAP for dimensionality reduction
    umap_op = UMAP(metric="precomputed")
    embedding = umap_op.fit_transform(combined_dist)

    # Store results
    adata.obsm["X_spatial_umap"] = embedding

    return adata


def infer_pseudotime_palantir(
    adata, root_cells=None, n_diffusion_components=10, num_waypoints=500
):
    """
    Infers cellular trajectories and pseudotime using the Palantir algorithm.

    Palantir models cellular differentiation as a stochastic process on a graph.
    It uses diffusion maps to capture the geometry of the data and then determines
    cell fate probabilities and differentiation potential by simulating random walks
    from a specified root cell. This produces a continuous pseudotime ordering of
    cells along developmental trajectories.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix, which must contain PCA results in `adata.obsm['X_pca']`.
    root_cells : list of str, optional
        A list of cell identifiers to be used as the starting point(s) of the trajectory.
        If not provided, Palantir will select a root cell automatically.
    n_diffusion_components : int, default 10
        The number of diffusion components to compute.
    num_waypoints : int, default 500
        The number of waypoints to use for modeling the trajectory, affecting granularity.

    Raises exceptions directly if Palantir fails - no fallback logic.
    """
    import palantir

    # Validate PCA requirement
    if "X_pca" not in adata.obsm:
        raise ValueError(
            "Palantir trajectory analysis requires PCA but X_pca not found. "
            "Please run PCA in preprocessing.py: sc.tl.pca(adata)"
        )

    # Convert PCA projections to DataFrame for Palantir
    pca_df = pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names)

    # Run diffusion maps with configurable components
    dm_res = palantir.utils.run_diffusion_maps(
        pca_df, n_components=n_diffusion_components
    )

    # Determine multiscale space
    ms_data = pd.DataFrame(dm_res["EigenVectors"], index=pca_df.index)

    # Determine early cell (start cell)
    if root_cells is not None and len(root_cells) > 0:
        if root_cells[0] not in ms_data.index:
            raise ValueError(f"Root cell '{root_cells[0]}' not found in data")
        start_cell = root_cells[0]
    else:
        # Use the cell at the extreme of the first diffusion component
        start_cell = ms_data.iloc[:, 0].idxmax()

    # Run Palantir with configurable waypoints
    pr_res = palantir.core.run_palantir(
        ms_data, start_cell, num_waypoints=num_waypoints
    )

    # Store pseudotime and branch probabilities
    adata.obs["palantir_pseudotime"] = pr_res.pseudotime
    adata.obsm["palantir_branch_probs"] = pr_res.branch_probs

    return adata


def compute_dpt_fallback(adata, root_cells=None):
    """Compute Diffusion Pseudotime as fallback method."""
    import numpy as np
    import scanpy as sc

    # Validate trajectory analysis prerequisites
    if "X_pca" not in adata.obsm:
        raise ValueError(
            "Diffusion pseudotime requires PCA but X_pca not found. "
            "Please run PCA in preprocessing.py: sc.tl.pca(adata)"
        )

    if "neighbors" not in adata.uns:
        raise ValueError(
            "Diffusion pseudotime requires neighborhood graph but neighbors not found. "
            "Please run in preprocessing.py: sc.pp.neighbors(adata, use_rep='X_pca')"
        )

    # Check if diffusion map has been computed
    if "X_diffmap" not in adata.obsm:
        raise ValueError(
            "Diffusion pseudotime requires diffusion map but X_diffmap not found. "
            "Please run in preprocessing.py: sc.tl.diffmap(adata)"
        )

    # Set root cell if provided
    if root_cells is not None and len(root_cells) > 0:
        if root_cells[0] in adata.obs_names:
            adata.uns["iroot"] = np.where(adata.obs_names == root_cells[0])[0][0]
        else:
            # If provided root cell not found, use first cell and warn
            # Root cell not found - using first cell as fallback
            adata.uns["iroot"] = 0
    else:
        # If no root cell specified, set the first cell as root
        adata.uns["iroot"] = 0

    # Validate or compute diffusion pseudotime
    if "dpt_pseudotime" not in adata.obs:
        # DPT needs to be computed - this is the core algorithm, not preprocessing
        try:
            import scanpy as sc

            sc.tl.dpt(adata)
        except Exception as e:
            # DPT computation failed - do not create fake pseudotime
            raise RuntimeError(
                f"Standard DPT computation failed: {e}. "
                "This indicates a problem with the data preprocessing or diffusion map computation. "
                "Please check that PCA, neighbors, and diffusion map were computed correctly."
            )

    # Check if dpt_pseudotime was created
    if "dpt_pseudotime" not in adata.obs.columns:
        raise RuntimeError("DPT computation did not create 'dpt_pseudotime' column")

    # Handle any NaN values
    adata.obs["dpt_pseudotime"] = adata.obs["dpt_pseudotime"].fillna(0)

    return adata


async def analyze_rna_velocity(
    data_id: str,
    data_store: Dict[str, Any],
    params: RNAVelocityParameters = RNAVelocityParameters(),
    context: Optional[Context] = None,
) -> RNAVelocityResult:
    """
    Computes RNA velocity for spatial transcriptomics data.

    This function requires the input dataset to contain 'spliced' and 'unspliced'
    count layers. It uses the scVelo library to preprocess the data and compute
    velocity vectors, which indicate the predicted future state of individual cells.

    Args:
        data_id: The identifier for the dataset.
        data_store: A dictionary that stores the loaded datasets.
        params: An object containing parameters for the RNA velocity analysis.
        context: The MCP context for logging and communication.

    Returns:
        An RNAVelocityResult object containing metadata about the computation.

    Raises:
        ProcessingError: If scVelo is not installed or if the velocity computation fails.
        DataNotFoundError: If the input data lacks the required 'spliced'/'unspliced' layers.
    """
    try:
        import scvelo as scv
    except ImportError:
        raise ProcessingError(
            "scvelo package not found. Install it with: pip install scvelo>=0.2.5"
        )

    if context:
        await context.info(f"Analyzing RNA velocity using {params.mode} mode")

    # Get AnnData object
    if data_id not in data_store:
        raise DataNotFoundError(f"Dataset {data_id} not found in data store")

    # Copy AnnData to avoid modifying original data
    adata = data_store[data_id]["adata"].copy()

    # Validate data for velocity analysis
    is_valid, issues = validate_velocity_data(adata)
    if not is_valid:
        error_message = (
            "The dataset is missing required data for RNA velocity analysis. "
            f"Specific issues: {'; '.join(issues)}. Please ensure the data "
            "contains both 'spliced' and 'unspliced' count layers."
        )
        if context:
            await context.error(error_message)
        raise DataNotFoundError(error_message)

    velocity_computed = False
    with suppress_output():
        try:
            if context:
                await context.info(
                    "Found spliced/unspliced layers. Computing velocity directly."
                )

            # Use unified velocity computation
            adata = compute_rna_velocity(adata, params=params)
            velocity_computed = True

        except Exception as e:
            error_msg = f"RNA velocity analysis failed: {str(e)}"
            if context:
                await context.error(error_msg)
            raise ProcessingError(error_msg) from e

    # Update data store
    data_store[data_id]["adata"] = adata

    # Return result with metadata only
    return RNAVelocityResult(
        data_id=data_id,
        velocity_computed=velocity_computed,
        velocity_graph_key="velocity_graph" if velocity_computed else None,
        mode=params.mode,
    )


async def analyze_trajectory(
    data_id: str,
    data_store: Dict[str, Any],
    params: TrajectoryParameters = TrajectoryParameters(),
    context: Optional[Context] = None,
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
    has_velocity = "velocity_graph" in adata.uns

    pseudotime_key = None
    method_used = params.method

    # Strategy 1: Try CellRank if requested and velocity is available
    if params.method == "cellrank" and has_velocity:
        if context:
            await context.info("Attempting trajectory inference with CellRank...")

        try:
            import cellrank as cr
        except ImportError:
            raise ProcessingError(
                "cellrank package not found. Install it with: pip install cellrank>=2.0.0"
            )

        try:
            with suppress_output():
                adata = infer_spatial_trajectory_cellrank(
                    adata,
                    spatial_weight=params.spatial_weight,
                    kernel_weights=params.cellrank_kernel_weights,
                    n_states=params.cellrank_n_states,
                )
            pseudotime_key = "pseudotime"
            if context:
                await context.info("CellRank analysis completed successfully.")
        except Exception as cellrank_error:
            if context:
                await context.error(f"CellRank analysis failed: {cellrank_error}")
            # Don't silently fallback - let the LLM decide what to do
            raise ProcessingError(
                f"CellRank trajectory inference failed: {cellrank_error}"
            ) from cellrank_error

    # Strategy 2: Try DPT if explicitly requested
    if params.method == "dpt":
        if context:
            await context.info("Attempting trajectory inference with DPT...")
        try:
            with suppress_output():
                adata = compute_dpt_fallback(adata, root_cells=params.root_cells)
            pseudotime_key = "dpt_pseudotime"
            method_used = "dpt"
            if context:
                await context.info("DPT analysis completed successfully.")
        except Exception as dpt_error:
            raise ProcessingError(f"DPT analysis failed: {dpt_error}") from dpt_error

    # Strategy 3: Try Palantir (either as primary method or fallback)
    elif method_used == "palantir" or not has_velocity:
        if context:
            await context.info("Attempting trajectory inference with Palantir...")

        try:
            with suppress_output():
                # Run spatially-aware embedding
                adata = spatial_aware_embedding(
                    adata, spatial_weight=params.spatial_weight
                )

                # Run Palantir with configurable parameters
                adata = infer_pseudotime_palantir(
                    adata,
                    root_cells=params.root_cells,
                    n_diffusion_components=params.palantir_n_diffusion_components,
                    num_waypoints=params.palantir_num_waypoints,
                )

            pseudotime_key = "palantir_pseudotime"
            method_used = "palantir"
            if context:
                await context.info("Palantir analysis completed successfully.")

        except Exception as palantir_error:
            if context:
                await context.error(f"Palantir analysis failed: {palantir_error}")
            # Don't silently fallback - let the LLM decide what to do
            # The LLM can explicitly request DPT method if it wants to try that
            raise ProcessingError(
                f"Palantir trajectory inference failed: {palantir_error}"
            ) from palantir_error

    # Ensure pseudotime key exists
    if pseudotime_key is None or pseudotime_key not in adata.obs.columns:
        raise ProcessingError("Failed to compute pseudotime with any available method")

    # Update data store
    data_store[data_id]["adata"] = adata

    # Determine if fallback was used
    requested_method = params.method if hasattr(params, "method") else "cellrank"
    fallback_used = (method_used == "dpt" and requested_method != "dpt") or (
        method_used == "palantir" and requested_method == "cellrank"
    )

    # Return result with metadata only (no visualization)
    return TrajectoryResult(
        data_id=data_id,
        pseudotime_computed=True,
        velocity_computed=has_velocity,
        pseudotime_key=pseudotime_key,
        method=method_used,
        spatial_weight=params.spatial_weight,
        method_fallback_used=fallback_used,
    )


async def analyze_velocity_with_velovi(
    adata,
    n_epochs: int = 1000,
    n_hidden: int = 128,
    n_latent: int = 10,
    use_gpu: bool = False,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """
    Analyzes RNA velocity using the deep learning model VELOVI.

    VELOVI (Velocity Variational Inference) is a probabilistic deep generative model
    that estimates transcriptional dynamics from spliced and unspliced mRNA counts.
    It provides a more detailed view of cellular dynamics by not only estimating
    velocity vectors but also quantifying the uncertainty associated with these
    estimates. This method is part of the scvi-tools ecosystem.

    Args:
        adata: The AnnData object, which must contain 'spliced' and 'unspliced' layers.
        n_epochs: The number of training epochs for the model.
        n_hidden: The number of hidden units in the neural network layers.
        n_latent: The dimensionality of the latent space.
        use_gpu: If True, training will be performed on a GPU if available.
        context: The MCP context for logging.

    Returns:
        A dictionary containing the results and metadata from the VELOVI analysis.

    Raises:
        ImportError: If the `scvi-tools` package is not installed.
        ValueError: If the input data is missing the required 'spliced' or 'unspliced' layers.
        RuntimeError: If an error occurs during model training or result extraction.
    """
    try:
        if scvi is None or VELOVI is None:
            raise ImportError(
                "scvi-tools package is required for VELOVI analysis. Install with 'pip install scvi-tools'"
            )

        if context:
            await context.info("Starting VELOVI velocity analysis...")

        # Validate required layers
        is_valid, issues = validate_velocity_data(adata)
        if not is_valid:
            raise ValueError(f"Invalid data for velocity analysis: {'; '.join(issues)}")

        if context:
            await context.info(
                f"Analyzing velocity for {adata.n_obs} cells and {adata.n_vars} genes with VELOVI"
            )

        # Setup VELOVI
        VELOVI.setup_anndata(
            adata, spliced_layer="spliced", unspliced_layer="unspliced"
        )

        # Create VELOVI model
        velovi_model = VELOVI(adata, n_hidden=n_hidden, n_latent=n_latent)

        if context:
            await context.info("Training VELOVI model...")

        # Train model
        if use_gpu:
            velovi_model.train(max_epochs=n_epochs, accelerator="gpu")
        else:
            velovi_model.train(max_epochs=n_epochs)

        if context:
            await context.info("VELOVI training completed")

        # Get results
        if context:
            await context.info("Extracting VELOVI velocity results...")

        # Get velocity estimates
        velocity = velovi_model.get_velocity()

        # Handle different return types from get_velocity
        if isinstance(velocity, tuple):
            # Sometimes get_velocity returns (velocity, something_else)
            velocity = velocity[0]

        # Ensure velocity is numpy array
        if hasattr(velocity, "detach"):
            # It's a tensor, convert to numpy
            velocity = velocity.detach().cpu().numpy()
        elif hasattr(velocity, "toarray"):
            # It's a sparse matrix
            velocity = velocity.toarray()

        # Ensure velocity is 2D array
        if velocity.ndim == 1:
            velocity = velocity.reshape(-1, 1)

        # Get latent representation
        latent = velovi_model.get_latent_representation()

        # Get velocity uncertainty (if available)
        try:
            velocity_uncertainty = velovi_model.get_velocity_uncertainty()
            uncertainty_computed = True
        except AttributeError:
            try:
                uncertainty_result = velovi_model.get_directional_uncertainty()
                # Handle case where method returns a tuple
                if isinstance(uncertainty_result, tuple):
                    velocity_uncertainty = uncertainty_result[
                        0
                    ]  # Usually the first element is the uncertainty
                else:
                    velocity_uncertainty = uncertainty_result
                uncertainty_computed = True
            except (AttributeError, IndexError):
                velocity_uncertainty = None
                uncertainty_computed = False

        # Store results in adata
        adata.layers["velocity_velovi"] = velocity
        adata.obsm["X_velovi_latent"] = latent
        if velocity_uncertainty is not None:
            try:
                adata.layers["velocity_velovi_uncertainty"] = velocity_uncertainty
            except Exception:
                # If uncertainty storage fails, just skip it
                uncertainty_computed = False

        # Calculate velocity statistics
        velocity_norm = np.linalg.norm(velocity, axis=1)

        # Store velocity norm in obs
        adata.obs["velocity_velovi_norm"] = velocity_norm

        # Calculate summary statistics
        results = {
            "method": "VELOVI",
            "n_latent_dims": n_latent,
            "n_epochs": n_epochs,
            "velocity_computed": True,
            "velocity_mean_norm": velocity_norm.mean(),
            "velocity_std_norm": velocity_norm.std(),
            "velocity_shape": velocity.shape,
            "latent_shape": latent.shape,
            "uncertainty_computed": uncertainty_computed,
            "training_completed": True,
            "device": "GPU" if use_gpu else "CPU",
        }

        if context:
            await context.info("VELOVI velocity analysis completed successfully")
            await context.info("Stored velocity in adata.layers['velocity_velovi']")
            await context.info(
                "Stored latent representation in adata.obsm['X_velovi_latent']"
            )
            await context.info(
                "Stored uncertainty in adata.layers['velocity_velovi_uncertainty']"
            )

        return results

    except Exception as e:
        error_msg = f"VELOVI velocity analysis failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg) from e
