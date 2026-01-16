"""
RNA velocity analysis for spatial transcriptomics.

This module computes RNA velocity to infer the direction of cellular state changes
by analyzing the balance of spliced and unspliced mRNA counts.

Key functionality:
- `analyze_rna_velocity`: Main MCP entry point for velocity analysis
- Supports scVelo (standard) and VELOVI (deep learning) methods
"""

from typing import TYPE_CHECKING, Any, Optional

import numpy as np

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

from ..models.analysis import RNAVelocityResult
from ..models.data import RNAVelocityParameters
from ..utils.adata_utils import store_velovi_essential_data, validate_adata
from ..utils.dependency_manager import require
from ..utils.exceptions import (
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)
from ..utils.mcp_utils import suppress_output


def preprocess_for_velocity(
    adata, min_shared_counts=30, n_top_genes=2000, n_pcs=30, n_neighbors=30, params=None
):
    """
    Prepares an AnnData object for RNA velocity analysis using the scVelo pipeline.

    This function performs the standard scVelo preprocessing workflow:
    1. Filtering genes based on minimum shared counts between spliced and
       unspliced layers.
    2. Normalizing the data.
    3. Selecting a subset of highly variable genes.
    4. Computing first and second-order moments across nearest neighbors.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix with 'spliced' and 'unspliced' layers.
    min_shared_counts : int, default 30
        Minimum number of counts shared between spliced and unspliced layers.
    n_top_genes : int, default 2000
        Number of highly variable genes to use.
    n_pcs : int, default 30
        Number of principal components to compute.
    n_neighbors : int, default 30
        Number of nearest neighbors for moment computation.
    params : RNAVelocityParameters, optional
        If provided, overrides the individual parameters.
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
    # enforce=True ensures scvelo recomputes everything even if data was pre-normalized
    # This is important when running after MCP's general preprocessing step
    scv.pp.filter_and_normalize(
        adata, min_shared_counts=min_shared_counts, n_top_genes=n_top_genes, enforce=True
    )
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    return adata


def compute_rna_velocity(adata, mode="stochastic", params=None):
    """
    Computes RNA velocity to infer the direction of cellular differentiation.

    This function executes the core RNA velocity workflow:
    1. Ensures preprocessing (moment computation) is complete.
    2. Estimates RNA velocity using the specified model.
    3. Constructs a velocity graph for cell-to-cell transitions.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix with 'spliced' and 'unspliced' layers.
    mode : str, default 'stochastic'
        The model for velocity estimation:
        - 'stochastic': Likelihood-based model accounting for noise.
        - 'deterministic': Simpler steady-state model.
        - 'dynamical': Full transcriptional dynamics with ODE fitting.
    params : RNAVelocityParameters, optional
        Parameter object (mode will be extracted from params.scvelo_mode).

    Returns
    -------
    AnnData
        Updated with velocity vectors and graph.
    """
    import scvelo as scv

    # Use params for mode if provided
    if params is not None:
        from ..models.data import RNAVelocityParameters

        if isinstance(params, RNAVelocityParameters):
            mode = params.scvelo_mode

    # Check if preprocessing is needed
    if "Ms" not in adata.layers or "Mu" not in adata.layers:
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


async def _prepare_velovi_data(adata, ctx: "ToolContext"):
    """Prepare data for VELOVI according to official standards."""
    import scvelo as scv

    adata_velovi = adata.copy()

    # Convert layer names to VELOVI standards
    if "spliced" in adata_velovi.layers and "unspliced" in adata_velovi.layers:
        adata_velovi.layers["Ms"] = adata_velovi.layers["spliced"]
        adata_velovi.layers["Mu"] = adata_velovi.layers["unspliced"]
    else:
        raise DataNotFoundError("Missing required 'spliced' and 'unspliced' layers")

    # scvelo preprocessing
    # enforce=True ensures scvelo recomputes everything even if data was pre-normalized
    try:
        scv.pp.filter_and_normalize(
            adata_velovi, min_shared_counts=30, n_top_genes=2000, enforce=True
        )
    except Exception as e:
        await ctx.warning(f"scvelo preprocessing warning: {e}")

    # Compute moments
    try:
        scv.pp.moments(adata_velovi, n_pcs=30, n_neighbors=30)
    except Exception as e:
        await ctx.warning(f"moments computation warning: {e}")

    return adata_velovi


def _validate_velovi_data(adata):
    """VELOVI-specific data validation."""
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
) -> dict[str, Any]:
    """
    Analyzes RNA velocity using the deep learning model VELOVI.

    VELOVI (Velocity Variational Inference) is a probabilistic deep generative model
    that estimates transcriptional dynamics from spliced and unspliced mRNA counts.
    It provides velocity vectors with uncertainty quantification.

    Args:
        adata: AnnData with 'spliced' and 'unspliced' layers.
        n_epochs: Number of training epochs.
        n_hidden: Number of hidden units in neural network layers.
        n_latent: Dimensionality of the latent space.
        use_gpu: If True, use GPU for training.
        ctx: ToolContext for logging.

    Returns:
        Dictionary with VELOVI results and metadata.
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

        # Handle pandas/numpy compatibility
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
            t_max = np.max(t, axis=0)
            if np.all(t_max > 0):
                scaling = 20 / t_max
            else:
                scaling = np.where(t_max > 0, 20 / t_max, 1.0)
        else:
            t_max = np.max(t)
            scaling = 20 / t_max if t_max > 0 else 1.0

        if hasattr(scaling, "to_numpy"):
            scaling = scaling.to_numpy()
        scaling = np.asarray(scaling)

        # Calculate scaled velocities
        if scaling.ndim == 0:
            scaled_velocities = velocities / scaling
        elif scaling.ndim == 1 and velocities.ndim == 2:
            scaled_velocities = velocities / scaling[np.newaxis, :]
        else:
            scaled_velocities = velocities / scaling

        # Store results in preprocessed data object
        adata_prepared.layers["velocity_velovi"] = scaled_velocities
        adata_prepared.layers["latent_time_velovi"] = latent_time
        adata_prepared.obsm["X_velovi_latent"] = latent_repr

        # Calculate velocity statistics
        velocity_norm = np.linalg.norm(scaled_velocities, axis=1)
        adata_prepared.obs["velocity_velovi_norm"] = velocity_norm

        # Transfer key information back to original adata
        adata.obs["velocity_velovi_norm"] = velocity_norm
        adata.obsm["X_velovi_latent"] = latent_repr

        # Store essential data for CellRank (optimized: ~78% memory savings)
        # Instead of storing full adata (~160 MB), stores only velocity/neighbors (~35 MB)
        store_velovi_essential_data(adata, adata_prepared)

        return {
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
            "training_completed": True,
            "device": "GPU" if use_gpu else "CPU",
        }

    except Exception as e:
        raise ProcessingError(f"VELOVI velocity analysis failed: {e}") from e


async def analyze_rna_velocity(
    data_id: str,
    ctx: "ToolContext",
    params: RNAVelocityParameters = RNAVelocityParameters(),
) -> RNAVelocityResult:
    """
    Computes RNA velocity for spatial transcriptomics data.

    This is the main MCP entry point for velocity analysis. It requires
    'spliced' and 'unspliced' count layers in the input dataset.

    Args:
        data_id: Dataset identifier.
        ctx: ToolContext for data access and logging.
        params: RNA velocity parameters.

    Returns:
        RNAVelocityResult with computation metadata.

    Raises:
        DataNotFoundError: If data lacks required layers.
        ProcessingError: If velocity computation fails.
    """
    require("scvelo")
    import scvelo as scv  # noqa: F401

    # Get AnnData object
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

    # Dispatch based on method
    if params.method == "scvelo":
        with suppress_output():
            try:
                adata = compute_rna_velocity(
                    adata, mode=params.scvelo_mode, params=params
                )
                velocity_computed = True
            except Exception as e:
                raise ProcessingError(
                    f"scVelo RNA velocity analysis failed: {e}"
                ) from e

    elif params.method == "velovi":
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

            if velovi_results.get("velocity_computed", False):
                velocity_computed = True
                adata.uns["velocity_graph"] = True
                adata.uns["velocity_method"] = "velovi"
            else:
                raise ProcessingError("VELOVI failed to compute velocity")

        except Exception as e:
            raise ProcessingError(f"VELOVI velocity analysis failed: {e}") from e

    else:
        raise ParameterError(f"Unknown velocity method: {params.method}")

    return RNAVelocityResult(
        data_id=data_id,
        velocity_computed=velocity_computed,
        velocity_graph_key="velocity_graph" if velocity_computed else None,
        mode=velocity_method_used if params.method == "scvelo" else params.method,
    )
