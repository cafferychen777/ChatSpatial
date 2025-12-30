"""
A module for RNA velocity and trajectory inference in spatial transcriptomics.

This module provides a suite of tools for analyzing cellular dynamics by
integrating single-cell RNA velocity with spatial information. It includes
functions for preprocessing data with spliced and unspliced counts, computing
RNA velocity using various models (e.g., stochastic, dynamical), and inferring
cellular trajectories.

Key functionalities are organized into two main analysis pipelines:
1. `analyze_rna_velocity`: Computes RNA velocity directly from spatial data
   that contains both spliced and unspliced count layers. Supports multiple
   velocity methods including scVelo (standard), VELOVI (deep learning),
   and SIRV (reference-based).
2. `analyze_trajectory`: Infers developmental trajectories and pseudotime by
   combining velocity information with spatial context. It provides access to
   advanced methods including CellRank (velocity-based), Palantir (expression-based),
   and DPT (diffusion-based). Each method has specific data requirements and
   should be explicitly selected.
"""

from typing import TYPE_CHECKING, Any, Dict, Optional

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

from ..models.analysis import RNAVelocityResult, TrajectoryResult
from ..models.data import RNAVelocityParameters, TrajectoryParameters
from ..utils.adata_utils import (
    get_spatial_key,
    require_spatial_coords,
    validate_adata,
    validate_obs_column,
)
from ..utils.compute import ensure_diffmap, ensure_neighbors, ensure_pca
from ..utils.dependency_manager import require
from ..utils.exceptions import (
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)
from ..utils.mcp_utils import suppress_output


def prepare_gam_model_for_visualization(
    adata,
    genes: list,
    time_key: str = "latent_time",
    fate_key: str = "lineages_fwd",
):
    """
    Prepare a GAM model for CellRank gene trends visualization.

    This function handles the computation logic needed for CellRank 2.0 gene trends
    and fate heatmap visualizations. Requires data to be analyzed in memory via
    analyze_rna_velocity (dynamical mode) and analyze_trajectory (cellrank method).

    Note: CellRank's cr.pl.gene_trends() handles model.prepare() internally,
    so we just need to validate data and create the GAM model.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix with CellRank results (fate probabilities).
    genes : list
        List of gene names to prepare the model for.
    time_key : str, default 'latent_time'
        Key in adata.obs for pseudotime/latent time values.
    fate_key : str, default 'lineages_fwd'
        Key in adata.obsm for fate probabilities (CellRank Lineage object).

    Returns
    -------
    tuple
        (model, lineage_names) - The GAM model and list of lineage names.

    Raises
    ------
    ValueError
        If required data (time_key, fate_key, genes) is missing from adata.
    """
    # Use centralized dependency manager for consistent error handling
    require("cellrank")  # Raises ImportError with install instructions if missing
    from cellrank.models import GAM

    # Validate required data
    validate_obs_column(adata, time_key, "Time key")

    if fate_key not in adata.obsm:
        raise DataNotFoundError(
            f"Fate probabilities '{fate_key}' not found. Run analyze_trajectory first."
        )

    # Validate Lineage object has names (must be from in-memory analysis)
    fate_probs = adata.obsm[fate_key]
    if not hasattr(fate_probs, "names") or fate_probs.names is None:
        raise DataError(
            "Fate probabilities must be a CellRank Lineage object with names. "
            "This requires running the full analysis pipeline in memory:\n"
            "1. analyze_rna_velocity(data_id, params={'scvelo_mode': 'dynamical'})\n"
            "2. analyze_trajectory(data_id, params={'method': 'cellrank'})\n"
            "3. Then visualize with plot_type='trajectory', subtype='gene_trends'"
        )
    lineage_names = list(fate_probs.names)

    # Validate genes exist
    missing_genes = [g for g in genes if g not in adata.var_names]
    if missing_genes:
        raise DataNotFoundError(
            f"Genes not found in data: {missing_genes}. "
            f"Available genes: {list(adata.var_names[:10])}..."
        )

    # Create GAM model - don't call prepare(), let cr.pl functions handle it
    model = GAM(adata)

    return model, lineage_names


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
    try:
        validate_adata(adata, {}, check_velocity=True)
    except DataNotFoundError as e:
        raise DataError(f"Invalid velocity data: {e}") from e

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
            mode = params.scvelo_mode

    # Check if preprocessing is needed
    if "Ms" not in adata.layers or "Mu" not in adata.layers:
        # Run preprocessing
        adata = preprocess_for_velocity(adata, params=params)

    # Compute velocity based on mode
    if mode == "dynamical":
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode="dynamical")
        # Compute latent time (required for gene_trends visualization)
        scv.tl.latent_time(adata)
    else:
        scv.tl.velocity(adata, mode=mode)

    # Compute velocity graph
    scv.tl.velocity_graph(adata)

    return adata


def infer_spatial_trajectory_cellrank(
    adata, spatial_weight=0.5, kernel_weights=(0.8, 0.2), n_states=5
):
    """
    Infers cellular trajectories by combining RNA velocity with CellRank.

    This function uses CellRank to model cell-state transitions by constructing
    a transition matrix from multiple kernels. It combines:
    1. A velocity kernel, derived from RNA velocity.
    2. A connectivity kernel, based on transcriptomic similarity.
    3. (Optional) A custom spatial kernel, based on physical proximity.

    By analyzing the eigenvectors of the combined transition matrix, CellRank
    identifies macrostates (representative cell states), terminal states, and computes
    fate probabilities to map the paths of cellular development.

    Note: Spatial kernel is optional. CellRank can work with velocity + connectivity
    kernels only. If spatial coordinates are available and spatial_weight > 0, a
    spatial kernel will be added to enhance the analysis.

    Raises ProcessingError if CellRank computation fails.
    """
    import cellrank as cr
    from scipy.sparse import csr_matrix
    from scipy.spatial.distance import pdist, squareform

    # Check if spatial data is available
    spatial_key = get_spatial_key(adata)
    has_spatial = spatial_key is not None

    # Adjust spatial_weight if no spatial data
    if not has_spatial and spatial_weight > 0:
        # No spatial coordinates found, fall back to velocity + connectivity kernels only
        spatial_weight = 0

    # Create RNA velocity kernel
    # Handle different velocity methods
    if "velocity_method" in adata.uns and adata.uns["velocity_method"] == "velovi":
        # For VELOVI, we need to use the preprocessed data that contains Ms/Mu layers
        if "velovi_adata" in adata.uns:
            # Use the VELOVI preprocessed data which has the proper layers
            adata_for_cellrank = adata.uns["velovi_adata"]
            # Transfer spatial coordinates if available
            if has_spatial:
                adata_for_cellrank.obsm["spatial"] = adata.obsm[spatial_key]

            # VELOVI stores velocity as 'velocity_velovi', but CellRank expects 'velocity'
            # Create a reference to the standard velocity layer name
            if "velocity_velovi" in adata_for_cellrank.layers:
                adata_for_cellrank.layers["velocity"] = adata_for_cellrank.layers[
                    "velocity_velovi"
                ]

            # Use VELOVI data for velocity kernel
            vk = cr.kernels.VelocityKernel(adata_for_cellrank)
            vk.compute_transition_matrix()
        else:
            raise ProcessingError("VELOVI velocity data not found")
    else:
        # Standard velocity (scVelo) - use original adata
        adata_for_cellrank = adata
        vk = cr.kernels.VelocityKernel(adata_for_cellrank)
        vk.compute_transition_matrix()

    # Create expression similarity-based kernel
    ck = cr.kernels.ConnectivityKernel(adata_for_cellrank)
    ck.compute_transition_matrix()

    # Combine kernels based on spatial data availability
    vk_weight, ck_weight = kernel_weights

    if has_spatial and spatial_weight > 0:
        # Create custom spatial kernel
        spatial_coords = adata.obsm[spatial_key]
        spatial_dist = squareform(pdist(spatial_coords))
        spatial_sim = np.exp(-spatial_dist / spatial_dist.mean())
        spatial_kernel = csr_matrix(spatial_sim)

        # Create spatial kernel
        sk = cr.kernels.PrecomputedKernel(spatial_kernel, adata_for_cellrank)
        sk.compute_transition_matrix()

        # Combine all three kernels
        combined_kernel = (1 - spatial_weight) * (
            vk_weight * vk + ck_weight * ck
        ) + spatial_weight * sk
    else:
        # No spatial data or spatial_weight=0: use only velocity + connectivity
        combined_kernel = vk_weight * vk + ck_weight * ck

    # Use GPCCA for analysis
    g = cr.estimators.GPCCA(combined_kernel)

    # Compute eigendecomposition
    g.compute_eigendecomposition()

    # Compute macrostates
    try:
        g.compute_macrostates(n_states=n_states)
    except Exception as e:
        raise ProcessingError(
            f"CellRank failed with n_states={n_states}: {e}. "
            f"Try reducing n_states or use method='palantir'/'dpt'."
        ) from e

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

        # Store results in adata_for_cellrank (which could be velovi_adata)
        adata_for_cellrank.obs["pseudotime"] = pseudotime
        adata_for_cellrank.obsm["fate_probabilities"] = absorption_probs
        adata_for_cellrank.obs["terminal_states"] = g.terminal_states
    else:
        # No terminal states, use macrostates for pseudotime
        if hasattr(g, "macrostates") and g.macrostates is not None:
            # Use the macrostate memberships as a proxy for pseudotime
            macrostate_probs = g.macrostates_memberships
            # Use the first macrostate as the "early" state
            pseudotime = 1 - macrostate_probs[:, 0].X.flatten()
            # Store in adata_for_cellrank first for consistency
            adata_for_cellrank.obs["pseudotime"] = pseudotime
        else:
            raise ProcessingError(
                "CellRank could not compute either terminal states or macrostates"
            )

    # Always store macrostates if available
    if hasattr(g, "macrostates") and g.macrostates is not None:
        adata_for_cellrank.obs["macrostates"] = g.macrostates

    # Transfer all CellRank results back to original adata
    # This ensures consistency regardless of whether VELOVI or standard velocity was used
    if "pseudotime" in adata_for_cellrank.obs:
        adata.obs["pseudotime"] = adata_for_cellrank.obs["pseudotime"]
    if "terminal_states" in adata_for_cellrank.obs:
        adata.obs["terminal_states"] = adata_for_cellrank.obs["terminal_states"]
    if "macrostates" in adata_for_cellrank.obs:
        adata.obs["macrostates"] = adata_for_cellrank.obs["macrostates"]
    if "fate_probabilities" in adata_for_cellrank.obsm:
        adata.obsm["fate_probabilities"] = adata_for_cellrank.obsm["fate_probabilities"]

    # Update velovi_adata if it was used (to preserve CellRank results for future use)
    if "velocity_method" in adata.uns and adata.uns["velocity_method"] == "velovi":
        if "velovi_adata" in adata.uns:
            adata.uns["velovi_adata"] = adata_for_cellrank

    return adata


def spatial_aware_embedding(adata, spatial_weight=0.3):
    """Generate spatially-aware low-dimensional embedding"""
    from sklearn.metrics.pairwise import euclidean_distances
    from umap import UMAP

    # Validate and get spatial coordinates
    spatial_coords = require_spatial_coords(adata)

    # Ensure PCA has been computed (lazy computation)
    ensure_pca(adata)

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

    Raises ProcessingError if Palantir computation fails.
    """
    import palantir

    # Ensure PCA has been computed (lazy computation)
    ensure_pca(adata)

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
            raise ParameterError(f"Root cell '{root_cells[0]}' not found in data")
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


def compute_dpt_trajectory(adata, root_cells=None, ctx: Optional["ToolContext"] = None):
    """Compute Diffusion Pseudotime trajectory analysis."""
    import numpy as np
    import scanpy as sc

    # Ensure trajectory analysis prerequisites (lazy computation)
    ensure_pca(adata)
    ensure_neighbors(adata)
    ensure_diffmap(adata)

    # Set root cell if provided
    if root_cells is not None and len(root_cells) > 0:
        if root_cells[0] in adata.obs_names:
            adata.uns["iroot"] = np.where(adata.obs_names == root_cells[0])[0][0]
        else:
            raise ParameterError(
                f"Root cell '{root_cells[0]}' not found. "
                f"Use valid cell ID from adata.obs_names or omit to auto-select."
            )
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
            raise ProcessingError(f"DPT computation failed: {e}") from e

    # Check if dpt_pseudotime was created
    if "dpt_pseudotime" not in adata.obs.columns:
        raise ProcessingError("DPT computation did not create 'dpt_pseudotime' column")

    # Handle any NaN values
    adata.obs["dpt_pseudotime"] = adata.obs["dpt_pseudotime"].fillna(0)

    return adata


async def analyze_rna_velocity(
    data_id: str,
    ctx: "ToolContext",
    params: RNAVelocityParameters = RNAVelocityParameters(),
) -> RNAVelocityResult:
    """
    Computes RNA velocity for spatial transcriptomics data.

    This function requires the input dataset to contain 'spliced' and 'unspliced'
    count layers. It uses the scVelo library to preprocess the data and compute
    velocity vectors, which indicate the predicted future state of individual cells.

    Args:
        data_id: The identifier for the dataset.
        ctx: ToolContext for data access and logging.
        params: An object containing parameters for the RNA velocity analysis.

    Returns:
        An RNAVelocityResult object containing metadata about the computation.

    Raises:
        ProcessingError: If scVelo is not installed or if the velocity computation fails.
        DataNotFoundError: If the input data lacks the required 'spliced'/'unspliced' layers.
    """
    # Use centralized dependency manager for consistent error handling
    require("scvelo")  # Raises ImportError with install instructions if missing
    import scvelo as scv  # noqa: F401

    # Get AnnData object via ToolContext
    adata = await ctx.get_adata(data_id)

    # Validate data for velocity analysis
    try:
        validate_adata(adata, {}, check_velocity=True)
    except DataNotFoundError as e:
        raise DataNotFoundError(
            f"Missing velocity data: {e}. Requires 'spliced' and 'unspliced' layers."
        ) from e

    velocity_computed = False
    velocity_method_used = params.method

    # Branch based on velocity computation method
    if params.method == "scvelo":
        with suppress_output():
            try:
                adata = compute_rna_velocity(
                    adata, mode=params.scvelo_mode, params=params
                )
                velocity_computed = True

            except Exception as e:
                raise ProcessingError(
                    f"scVelo RNA velocity analysis failed: {str(e)}"
                ) from e

    elif params.method == "velovi":
        # VELOVI deep learning velocity computation
        require("scvi", feature="VELOVI velocity analysis")

        try:
            velovi_results = await analyze_velocity_with_velovi(
                adata,
                n_epochs=params.velovi_n_epochs,
                n_hidden=params.velovi_n_hidden,
                n_latent=params.velovi_n_latent,
                use_gpu=params.velovi_use_gpu,
                ctx=ctx,
            )

            # Check if velocity was successfully computed
            if velovi_results.get("velocity_computed", False):
                velocity_computed = True
                # VELOVI stores results in adata.uns["velovi_adata"] and adata.obs
                if "velovi_adata" in adata.uns:
                    adata.uns["velocity_graph"] = True  # Marker that velocity exists
                    adata.uns["velocity_method"] = "velovi"
            else:
                raise ProcessingError("VELOVI failed to compute velocity")

        except Exception as e:
            raise ProcessingError(f"VELOVI velocity analysis failed: {str(e)}") from e

    else:
        raise ParameterError(f"Unknown velocity method: {params.method}")

    # COW FIX: No need to update data_store - changes already reflected via direct reference
    # All modifications to adata.obs/uns/obsm are in-place and preserved

    # Return result with metadata only
    return RNAVelocityResult(
        data_id=data_id,
        velocity_computed=velocity_computed,
        velocity_graph_key="velocity_graph" if velocity_computed else None,
        mode=velocity_method_used if params.method == "scvelo" else params.method,
    )


async def analyze_trajectory(
    data_id: str,
    ctx: "ToolContext",
    params: TrajectoryParameters = TrajectoryParameters(),
) -> TrajectoryResult:
    """Analyze trajectory and cell state transitions in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        ctx: ToolContext for data access and logging
        params: Trajectory analysis parameters

    Returns:
        Trajectory analysis result (computation metadata only, no visualization)
    """
    # Get AnnData object via ToolContext
    adata = await ctx.get_adata(data_id)

    # Check if RNA velocity has been computed (by any method)
    has_velocity = (
        "velocity_graph" in adata.uns  # scVelo velocity
        or "velovi_adata" in adata.uns  # VELOVI velocity (stored separately)
        or "velocity_method" in adata.uns  # Generic velocity marker
    )

    pseudotime_key = None
    method_used = params.method

    # Execute requested method with explicit requirements checking
    if params.method == "cellrank":
        if not has_velocity:
            raise ProcessingError(
                "CellRank requires velocity data. Run velocity analysis first or use palantir/dpt."
            )

        require("cellrank")  # Raises ImportError with install instructions if missing
        import cellrank as cr  # noqa: F401

        try:
            with suppress_output():
                adata = infer_spatial_trajectory_cellrank(
                    adata,
                    spatial_weight=params.spatial_weight,
                    kernel_weights=params.cellrank_kernel_weights,
                    n_states=params.cellrank_n_states,
                )
            pseudotime_key = "pseudotime"
            method_used = "cellrank"
        except Exception as cellrank_error:
            raise ProcessingError(
                f"CellRank trajectory inference failed: {cellrank_error}"
            ) from cellrank_error

    elif params.method == "palantir":
        try:
            with suppress_output():
                # Optional: Run spatially-aware embedding if spatial data available
                has_spatial = get_spatial_key(adata) is not None
                if has_spatial and params.spatial_weight > 0:
                    adata = spatial_aware_embedding(
                        adata, spatial_weight=params.spatial_weight
                    )
                elif not has_spatial and params.spatial_weight > 0:
                    await ctx.warning(
                        f"Spatial weight {params.spatial_weight} specified but no spatial "
                        "coordinates found. Using expression-only Palantir."
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

        except Exception as palantir_error:
            raise ProcessingError(
                f"Palantir trajectory inference failed: {palantir_error}"
            ) from palantir_error

    elif params.method == "dpt":
        try:
            with suppress_output():
                adata = compute_dpt_trajectory(
                    adata, root_cells=params.root_cells, ctx=ctx
                )
            pseudotime_key = "dpt_pseudotime"
            method_used = "dpt"
        except Exception as dpt_error:
            raise ProcessingError(f"DPT analysis failed: {dpt_error}") from dpt_error

    else:
        raise ParameterError(f"Unknown trajectory method: {params.method}")

    # Ensure pseudotime key exists
    if pseudotime_key is None or pseudotime_key not in adata.obs.columns:
        raise ProcessingError("Failed to compute pseudotime with any available method")

    # COW FIX: No need to update data_store - changes already reflected via direct reference
    # All modifications to adata.obs/uns/obsm are in-place and preserved

    # Store scientific metadata for reproducibility
    from ..utils.adata_utils import store_analysis_metadata

    # Determine results keys based on method
    results_keys_dict = {"obs": [pseudotime_key], "obsm": [], "uns": []}

    if method_used == "cellrank":
        results_keys_dict["obs"].extend(["terminal_states", "macrostates"])
        results_keys_dict["obsm"].append("fate_probabilities")
        results_keys_dict["uns"].append("velocity_method")
    elif method_used == "palantir":
        results_keys_dict["obsm"].append("palantir_branch_probs")
    elif method_used == "dpt":
        results_keys_dict["uns"].append("iroot")

    # Prepare parameters dict
    parameters_dict: Dict[str, Any] = {"spatial_weight": params.spatial_weight}
    if method_used == "cellrank":
        parameters_dict.update(
            {
                "kernel_weights": params.cellrank_kernel_weights,
                "n_states": params.cellrank_n_states,
            }
        )
    elif method_used == "palantir":
        parameters_dict.update(
            {
                "n_diffusion_components": params.palantir_n_diffusion_components,
                "num_waypoints": params.palantir_num_waypoints,
            }
        )

    if params.root_cells:
        parameters_dict["root_cells"] = params.root_cells

    # Prepare statistics dict
    statistics_dict = {
        "velocity_computed": has_velocity,
        "pseudotime_key": pseudotime_key,
    }

    # Store metadata
    store_analysis_metadata(
        adata,
        analysis_name=f"trajectory_{method_used}",
        method=method_used,
        parameters=parameters_dict,
        results_keys=results_keys_dict,
        statistics=statistics_dict,
    )

    # Return result with metadata only (no visualization)
    return TrajectoryResult(
        data_id=data_id,
        pseudotime_computed=True,
        velocity_computed=has_velocity,
        pseudotime_key=pseudotime_key,
        method=method_used,
        spatial_weight=params.spatial_weight,
    )


async def _prepare_velovi_data(adata, ctx: "ToolContext"):
    """Prepare data for VELOVI according to official standards"""
    import scvelo as scv

    adata_velovi = adata.copy()

    # 1. Convert layer names to VELOVI standards
    if "spliced" in adata_velovi.layers and "unspliced" in adata_velovi.layers:
        adata_velovi.layers["Ms"] = adata_velovi.layers["spliced"]
        adata_velovi.layers["Mu"] = adata_velovi.layers["unspliced"]
    else:
        raise DataNotFoundError("Missing required 'spliced' and 'unspliced' layers")

    # 2. scvelo preprocessing
    try:
        scv.pp.filter_and_normalize(
            adata_velovi, min_shared_counts=30, n_top_genes=2000, enforce=False
        )
    except Exception as e:
        await ctx.warning(f"scvelo preprocessing warning: {e}")

    # 3. Compute moments
    try:
        scv.pp.moments(adata_velovi, n_pcs=30, n_neighbors=30)
    except Exception as e:
        await ctx.warning(f"moments computation warning: {e}")

    return adata_velovi


def _validate_velovi_data(adata):
    """VELOVI-specific data validation"""
    if "Ms" not in adata.layers or "Mu" not in adata.layers:
        raise DataNotFoundError("Missing required layers 'Ms' and 'Mu' for VELOVI")

    ms_data = adata.layers["Ms"]
    mu_data = adata.layers["Mu"]

    if ms_data.shape != mu_data.shape:
        raise DataError(f"Shape mismatch: Ms {ms_data.shape} vs Mu {mu_data.shape}")

    if ms_data.ndim != 2 or mu_data.ndim != 2:
        raise DataError(
            f"Expected 2D arrays, got Ms:{ms_data.ndim}D, Mu:{mu_data.ndim}D"
        )

    return True


async def analyze_velocity_with_velovi(
    adata,
    n_epochs: int = 1000,
    n_hidden: int = 128,
    n_latent: int = 10,
    use_gpu: bool = False,
    ctx: Optional["ToolContext"] = None,
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
        ctx: ToolContext for data access and logging.

    Returns:
        A dictionary containing the results and metadata from the VELOVI analysis.

    Raises:
        ImportError: If the `scvi-tools` package is not installed.
        ValueError: If the input data is missing the required 'spliced' or 'unspliced' layers.
        RuntimeError: If an error occurs during model training or result extraction.
    """
    try:
        require("scvi", feature="VELOVI velocity analysis")
        from scvi.external import VELOVI

        # Data preprocessing
        adata_prepared = await _prepare_velovi_data(adata, ctx)

        # Data validation
        _validate_velovi_data(adata_prepared)

        # VELOVI setup
        VELOVI.setup_anndata(
            adata_prepared,
            spliced_layer="Ms",
            unspliced_layer="Mu",
        )

        # Model creation
        velovi_model = VELOVI(adata_prepared, n_hidden=n_hidden, n_latent=n_latent)

        # Model training
        if use_gpu:
            velovi_model.train(max_epochs=n_epochs, accelerator="gpu")
        else:
            velovi_model.train(max_epochs=n_epochs)

        # Result extraction
        latent_time = velovi_model.get_latent_time(n_samples=25)
        velocities = velovi_model.get_velocity(n_samples=25, velo_statistic="mean")
        latent_repr = velovi_model.get_latent_representation()

        # Handle pandas/numpy compatibility issues
        if hasattr(latent_time, "values"):
            latent_time = latent_time.values
        if hasattr(velocities, "values"):
            velocities = velocities.values

        # Ensure numpy array format
        latent_time = np.asarray(latent_time)
        velocities = np.asarray(velocities)
        latent_repr = np.asarray(latent_repr)

        # Safe scaling calculation
        t = latent_time
        if t.ndim > 1:
            # Multi-dimensional array, calculate max value for each feature
            t_max = np.max(t, axis=0)
            # Safe condition check: use all() to handle arrays
            if np.all(t_max > 0):
                scaling = 20 / t_max
            else:
                # If some features have max value of 0, use 1 as default scaling
                scaling = np.where(t_max > 0, 20 / t_max, 1.0)
        else:
            # 1-dimensional array
            t_max = np.max(t)
            scaling = 20 / t_max if t_max > 0 else 1.0

        # Ensure scaling is numpy array
        if hasattr(scaling, "to_numpy"):
            scaling = scaling.to_numpy()
        scaling = np.asarray(scaling)

        # Calculate final scaled velocities
        if scaling.ndim == 0:
            scaled_velocities = velocities / scaling
        elif scaling.ndim == 1 and velocities.ndim == 2:
            scaled_velocities = velocities / scaling[np.newaxis, :]
        else:
            scaled_velocities = velocities / scaling

        # Key fix: store results in preprocessed data object (dimension matching)
        adata_prepared.layers["velocity_velovi"] = scaled_velocities
        adata_prepared.layers["latent_time_velovi"] = latent_time
        adata_prepared.obsm["X_velovi_latent"] = latent_repr

        # Calculate velocity statistics
        velocity_norm = np.linalg.norm(scaled_velocities, axis=1)
        adata_prepared.obs["velocity_velovi_norm"] = velocity_norm

        # Solution: transfer key information back to original data object, avoiding dimension conflicts
        # Only transfer cell-related information (obs info), not gene-related layer information
        adata.obs["velocity_velovi_norm"] = velocity_norm
        adata.obsm["X_velovi_latent"] = latent_repr

        # Store preprocessed data object in uns for future use
        adata.uns["velovi_adata"] = adata_prepared
        adata.uns["velovi_gene_names"] = adata_prepared.var_names.tolist()

        # Build results dictionary
        results = {
            "method": "VELOVI",
            "velocity_computed": True,
            "n_epochs": n_epochs,
            "n_hidden": n_hidden,
            "n_latent": n_latent,
            "velocity_shape": scaled_velocities.shape,
            "latent_time_shape": latent_time.shape,
            "latent_repr_shape": latent_repr.shape,
            "velocity_mean_norm": float(velocity_norm.mean()),
            "velocity_std_norm": float(velocity_norm.std()),
            "n_genes_analyzed": adata_prepared.n_vars,
            "original_n_genes": adata.n_vars,
            "scaling_shape": scaling.shape,
            "training_completed": True,
            "device": "GPU" if use_gpu else "CPU",
            "storage_note": "Full velocity results stored in adata.uns['velovi_adata']",
        }

        return results

    except Exception as e:
        raise ProcessingError(f"VELOVI velocity analysis failed: {str(e)}") from e
