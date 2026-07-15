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
import pandas as pd
from scipy import sparse

if TYPE_CHECKING:
    import anndata as ad

    from ..spatial_mcp_adapter import ToolContext

from ..models.analysis import RNAVelocityResult
from ..models.data import RNAVelocityParameters
from ..utils.adata_utils import (
    store_analysis_metadata,
    store_velovi_essential_data,
    validate_adata,
)
from ..utils.dependency_manager import require
from ..utils.device_utils import get_accelerator
from ..utils.exceptions import (
    ChatSpatialError,
    DataCompatibilityError,
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)
from ..utils.mcp_utils import suppress_output
from ..utils.results_export import export_analysis_result


def _build_velocity_key(params: "RNAVelocityParameters") -> str:
    """Build a parametric analysis key for RNA velocity results.

    Encodes method + mode so that different scvelo modes (stochastic,
    dynamical) or scvelo vs velovi coexist in metadata/export/cache.
    """
    if params.method == "scvelo":
        return f"velocity_scvelo_{params.scvelo_mode}"
    return f"velocity_{params.method}"


def _velocity_parameters_for_metadata(
    params: "RNAVelocityParameters",
) -> dict[str, Any]:
    """Return exactly the parameters that affect the selected velocity method."""
    metadata: dict[str, Any] = {
        "n_top_genes": params.n_top_genes,
        "n_pcs": params.n_pcs,
        "n_neighbors": params.n_neighbors,
        "min_shared_counts": params.min_shared_counts,
    }
    if params.method == "scvelo":
        metadata["mode"] = params.scvelo_mode
    else:
        metadata.update(
            {
                "n_epochs": params.velovi_n_epochs,
                "n_hidden": params.velovi_n_hidden,
                "n_latent": params.velovi_n_latent,
                "n_layers": params.velovi_n_layers,
                "dropout_rate": params.velovi_dropout_rate,
                "learning_rate": params.velovi_learning_rate,
                "use_gpu": params.velovi_use_gpu,
            }
        )
    return metadata


def _copy_matrix_data(data: Any) -> Any:
    """Return an independent copy of dense/sparse matrix-like data."""
    if hasattr(data, "copy"):
        return data.copy()
    return np.array(data, copy=True)


def _compute_velocity_moments(
    adata: "ad.AnnData", *, n_pcs: int, n_neighbors: int
) -> None:
    """Compute an explicit Scanpy neighbor graph and scVelo moments."""
    import scanpy as sc
    import scvelo as scv

    sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    # Passing None makes scVelo consume the graph above instead of invoking its
    # deprecated implicit neighbor calculation.
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
    if "Ms" not in adata.layers or "Mu" not in adata.layers:
        raise ProcessingError("scVelo did not produce the required Ms and Mu moments.")


def _coerce_velovi_numeric_array(result: Any, output_name: str) -> np.ndarray:
    """Convert a VELOVI output to a finite real-valued float array."""
    raw_values = np.asarray(result)
    if np.iscomplexobj(raw_values):
        raise DataCompatibilityError(
            f"VELOVI {output_name} contains complex values."
        )
    try:
        values = np.asarray(raw_values, dtype=np.float32)
    except (TypeError, ValueError) as exc:
        raise DataCompatibilityError(
            f"VELOVI {output_name} must be numeric."
        ) from exc
    if not np.isfinite(values).all():
        raise DataCompatibilityError(
            f"VELOVI {output_name} contains non-finite values."
        )
    return values


def _normalize_velovi_cell_gene_output(
    result: Any,
    adata: "ad.AnnData",
    output_name: str,
    *,
    require_nonnegative: bool = False,
) -> np.ndarray:
    """Validate a VELOVI cell-by-gene output and restore canonical axis order."""
    expected_cells = pd.Index(adata.obs_names.map(str), name="cell")
    expected_genes = pd.Index(adata.var_names.map(str), name="gene")
    if not expected_cells.is_unique:
        raise DataCompatibilityError(
            "VELOVI workspace contains duplicate observation IDs."
        )
    if not expected_genes.is_unique:
        raise DataCompatibilityError("VELOVI workspace contains duplicate gene IDs.")

    expected_shape = (adata.n_obs, adata.n_vars)
    if isinstance(result, pd.DataFrame):
        normalized = result.copy()
        normalized.index = normalized.index.map(str)
        normalized.columns = normalized.columns.map(str)
        if not normalized.index.is_unique:
            raise DataCompatibilityError(
                f"VELOVI {output_name} contains duplicate observation IDs."
            )
        if not normalized.columns.is_unique:
            raise DataCompatibilityError(
                f"VELOVI {output_name} contains duplicate gene IDs."
            )

        missing_cells = expected_cells.difference(normalized.index)
        unexpected_cells = normalized.index.difference(expected_cells)
        missing_genes = expected_genes.difference(normalized.columns)
        unexpected_genes = normalized.columns.difference(expected_genes)
        if any(
            len(labels)
            for labels in (
                missing_cells,
                unexpected_cells,
                missing_genes,
                unexpected_genes,
            )
        ):
            raise DataCompatibilityError(
                f"VELOVI {output_name} axes do not match the model workspace: "
                f"{len(missing_cells)} missing and {len(unexpected_cells)} "
                "unexpected observations; "
                f"{len(missing_genes)} missing and {len(unexpected_genes)} "
                "unexpected genes."
            )
        raw_values = normalized.reindex(
            index=expected_cells, columns=expected_genes
        ).to_numpy()
    else:
        raw_values = np.asarray(result)
        if raw_values.ndim != 2 or raw_values.shape != expected_shape:
            raise DataCompatibilityError(
                f"VELOVI {output_name} has shape {raw_values.shape}; "
                f"expected {expected_shape}."
            )

    values = _coerce_velovi_numeric_array(raw_values, output_name)
    if values.shape != expected_shape:
        raise DataCompatibilityError(
            f"VELOVI {output_name} has shape {values.shape}; "
            f"expected {expected_shape}."
        )
    if require_nonnegative and np.any(values < 0):
        raise DataCompatibilityError(
            f"VELOVI {output_name} contains negative values."
        )
    return values


def _normalize_velovi_latent_representation(
    result: Any, adata: "ad.AnnData"
) -> np.ndarray:
    """Validate the positional cell-by-latent output from VELOVI."""
    values = _coerce_velovi_numeric_array(result, "latent representation")
    if values.ndim != 2 or values.shape[0] != adata.n_obs or values.shape[1] == 0:
        raise DataCompatibilityError(
            "VELOVI latent representation has shape "
            f"{values.shape}; expected ({adata.n_obs}, n_latent) with n_latent > 0."
        )
    return values


def _validate_velovi_workspace_axis(
    original: "ad.AnnData", prepared: "ad.AnnData"
) -> None:
    """Ensure preprocessing preserved the cell axis used by stored graph matrices."""
    original_cells = pd.Index(original.obs_names.map(str))
    prepared_cells = pd.Index(prepared.obs_names.map(str))
    if not original_cells.is_unique or not prepared_cells.is_unique:
        raise DataCompatibilityError(
            "VELOVI requires unique observation IDs before preprocessing."
        )
    if not prepared_cells.equals(original_cells):
        raise DataCompatibilityError(
            "VELOVI preprocessing changed the observation identities or order."
        )
    prepared_genes = pd.Index(prepared.var_names.map(str))
    if len(prepared_genes) == 0:
        raise DataCompatibilityError("VELOVI preprocessing retained no genes.")
    if not prepared_genes.is_unique:
        raise DataCompatibilityError(
            "VELOVI preprocessing produced duplicate gene IDs."
        )


def preprocess_for_velocity(
    adata: "ad.AnnData",
    min_shared_counts: int = 30,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 30,
    params: Optional[RNAVelocityParameters] = None,
) -> "ad.AnnData":
    """Prepare AnnData for RNA velocity analysis using scVelo pipeline.

    Performs the standard scVelo preprocessing workflow:
    1. Filtering genes based on minimum shared counts between spliced/unspliced layers
    2. Normalizing the data
    3. Selecting highly variable genes
    4. Computing first and second-order moments across nearest neighbors

    Args:
        adata: Annotated data matrix with 'spliced' and 'unspliced' layers.
        min_shared_counts: Minimum counts shared between spliced and unspliced layers.
        n_top_genes: Number of highly variable genes to use.
        n_pcs: Number of principal components to compute.
        n_neighbors: Number of nearest neighbors for moment computation.
        params: If provided, overrides the individual parameters.

    Returns:
        Preprocessed AnnData with velocity-ready layers and moments.
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
        adata,
        min_shared_counts=min_shared_counts,
        n_top_genes=n_top_genes,
        enforce=True,
    )
    _compute_velocity_moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    return adata


def compute_rna_velocity(
    adata: "ad.AnnData",
    mode: str = "stochastic",
    params: Optional[RNAVelocityParameters] = None,
) -> "ad.AnnData":
    """Compute RNA velocity to infer cellular differentiation direction.

    Executes the core RNA velocity workflow:
    1. Ensures preprocessing (moment computation) is complete
    2. Estimates RNA velocity using the specified model
    3. Constructs a velocity graph for cell-to-cell transitions

    Args:
        adata: Annotated data matrix with 'spliced' and 'unspliced' layers.
        mode: Model for velocity estimation:
            - 'stochastic': Likelihood-based model accounting for noise
            - 'deterministic': Simpler steady-state model
            - 'dynamical': Full transcriptional dynamics with ODE fitting
        params: Parameter object (mode extracted from params.scvelo_mode).

    Returns:
        AnnData updated with velocity vectors and graph.
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


async def _prepare_velovi_data(
    adata: "ad.AnnData",
    ctx: Optional["ToolContext"],
    params: Optional["RNAVelocityParameters"] = None,
) -> "ad.AnnData":
    """Prepare data for VELOVI according to official standards.

    Args:
        adata: Input AnnData with spliced/unspliced layers.
        ctx: ToolContext for logging warnings.
        params: Velocity parameters for preprocessing settings.

    Returns:
        Preprocessed AnnData copy ready for VELOVI.
    """
    import anndata as ad
    import scvelo as scv

    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise DataNotFoundError("Missing required 'spliced' and 'unspliced' layers")

    # Build a minimal working AnnData to avoid copying unrelated layers/obsm/uns.
    spliced = _copy_matrix_data(adata.layers["spliced"])
    unspliced = _copy_matrix_data(adata.layers["unspliced"])
    adata_velovi = ad.AnnData(
        X=spliced,
        obs=adata.obs.copy(),
        var=adata.var.copy(),
    )
    adata_velovi.obs_names = adata.obs_names.copy()
    adata_velovi.var_names = adata.var_names.copy()
    adata_velovi.layers["spliced"] = spliced
    adata_velovi.layers["unspliced"] = unspliced

    # Convert layer names to VELOVI standards.
    adata_velovi.layers["Ms"] = adata_velovi.layers["spliced"]
    adata_velovi.layers["Mu"] = adata_velovi.layers["unspliced"]

    # Resolve preprocessing parameters from params (same defaults as scvelo path)
    min_shared_counts = 30
    n_top_genes = 2000
    n_pcs = 30
    n_neighbors = 30
    if params is not None:
        min_shared_counts = params.min_shared_counts
        n_top_genes = params.n_top_genes
        n_pcs = params.n_pcs
        n_neighbors = params.n_neighbors

    # scvelo preprocessing
    # enforce=True ensures scvelo recomputes everything even if data was pre-normalized
    try:
        scv.pp.filter_and_normalize(
            adata_velovi,
            min_shared_counts=min_shared_counts,
            n_top_genes=n_top_genes,
            enforce=True,
        )
    except Exception as e:
        raise ProcessingError(
            f"VELOVI preprocessing failed: {e}. "
            "Check that spliced/unspliced layers exist "
            "and contain valid data."
        ) from e

    # Compute moments
    try:
        _compute_velocity_moments(
            adata_velovi, n_pcs=n_pcs, n_neighbors=n_neighbors
        )
    except Exception as e:
        raise ProcessingError(
            f"Moments computation failed: {e}. "
            "This may indicate insufficient cells or genes "
            "after filtering."
        ) from e

    return adata_velovi


def _validate_velovi_data(adata: "ad.AnnData") -> None:
    """Validate data for VELOVI requirements.

    Args:
        adata: AnnData to validate for VELOVI compatibility.

    Raises:
        DataNotFoundError: If required layers are missing.
        DataError: If layer shapes are incompatible.
    """
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

    for layer_name, layer_data in (("Ms", ms_data), ("Mu", mu_data)):
        raw_values = layer_data.data if sparse.issparse(layer_data) else layer_data
        values = _coerce_velovi_numeric_array(raw_values, f"{layer_name} layer")
        if np.any(values < 0):
            raise DataError(f"VELOVI {layer_name} layer contains negative values.")
        if not np.any(values > 0):
            raise DataError(f"VELOVI {layer_name} layer contains no positive values.")


async def analyze_velocity_with_velovi(
    adata: "ad.AnnData",
    params: Optional["RNAVelocityParameters"] = None,
    *,
    n_epochs: int = 1000,
    n_hidden: int = 128,
    n_latent: int = 10,
    n_layers: int = 1,
    dropout_rate: float = 0.1,
    learning_rate: float = 1e-3,
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
        params: Full parameter object (preprocessing params forwarded to data prep).
        n_epochs: Number of training epochs.
        n_hidden: Number of hidden units in neural network layers.
        n_latent: Dimensionality of the latent space.
        n_layers: Number of hidden layers in encoder/decoder.
        dropout_rate: Dropout rate for regularization.
        learning_rate: Learning rate for Adam optimizer.
        use_gpu: If True, use GPU for training.
        ctx: ToolContext for logging.

    Returns:
        Dictionary with VELOVI results and metadata.
    """
    try:
        require("scvelo", feature="VELOVI preprocessing")
        require("scvi", feature="VELOVI velocity analysis")
        from scvi.external import VELOVI

        # Data preprocessing (forward params so user's preprocessing settings apply)
        adata_prepared = await _prepare_velovi_data(adata, ctx, params=params)
        _validate_velovi_workspace_axis(adata, adata_prepared)

        # Data validation
        _validate_velovi_data(adata_prepared)

        # VELOVI setup
        VELOVI.setup_anndata(
            adata_prepared,
            spliced_layer="Ms",
            unspliced_layer="Mu",
        )

        # Model creation (all architecture params forwarded)
        velovi_model = VELOVI(
            adata_prepared,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
        )

        # Model training
        accelerator = get_accelerator(prefer_gpu=use_gpu)
        if use_gpu and accelerator == "cpu" and ctx is not None:
            await ctx.warning("GPU requested but unavailable; VELOVI is using CPU.")
        velovi_model.train(
            max_epochs=n_epochs,
            accelerator=accelerator,
            lr=learning_rate,
        )

        # Result extraction
        latent_time_result = velovi_model.get_latent_time(
            n_samples=25, return_numpy=False
        )
        velocity_result = velovi_model.get_velocity(
            n_samples=25,
            velo_statistic="mean",
            return_numpy=False,
        )
        latent_repr_result = velovi_model.get_latent_representation()

        latent_time = _normalize_velovi_cell_gene_output(
            latent_time_result,
            adata_prepared,
            "latent time",
            require_nonnegative=True,
        )
        velocities = _normalize_velovi_cell_gene_output(
            velocity_result,
            adata_prepared,
            "velocity",
        )
        latent_repr = _normalize_velovi_latent_representation(
            latent_repr_result, adata_prepared
        )

        # Match the official per-gene scaling while retaining a defined value
        # for genes whose inferred latent time is identically zero.
        t_max = latent_time.max(axis=0).astype(np.float64, copy=False)
        scaling = np.ones_like(t_max)
        positive_time = t_max > 0
        scaling[positive_time] = 20.0 / t_max[positive_time]
        scaled_velocities = np.asarray(
            velocities / scaling[np.newaxis, :], dtype=np.float32
        )
        if not np.isfinite(scaled_velocities).all():
            raise DataCompatibilityError(
                "VELOVI scaled velocity contains non-finite values."
            )

        # Store results in preprocessed data object.
        adata_prepared.layers["velocity_velovi"] = scaled_velocities
        adata_prepared.layers["latent_time_velovi"] = latent_time
        adata_prepared.obsm["X_velovi_latent"] = latent_repr

        # Calculate velocity statistics
        velocity_norm = np.linalg.norm(scaled_velocities, axis=1)
        adata_prepared.obs["velocity_velovi_norm"] = velocity_norm

        # Validate and persist the complete CellRank payload before publishing
        # user-facing result keys on the original object.
        store_velovi_essential_data(adata, adata_prepared)

        # Transfer key information back to original adata
        adata.obs["velocity_velovi_norm"] = velocity_norm
        adata.obsm["X_velovi_latent"] = latent_repr
        adata.obsm["velovi_latent_time"] = pd.DataFrame(
            latent_time,
            index=adata.obs_names.copy(),
            columns=adata_prepared.var_names.copy(),
        )

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
            "device": accelerator.upper(),
        }

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"VELOVI velocity analysis failed: {e}") from e


async def analyze_rna_velocity(
    data_id: str,
    ctx: "ToolContext",
    params: RNAVelocityParameters | None = None,
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
    if params is None:
        params = RNAVelocityParameters()

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
    # Dispatch based on method
    if params.method == "scvelo":
        require("scvelo", feature="scVelo RNA velocity analysis")
        with suppress_output():
            try:
                adata = compute_rna_velocity(
                    adata, mode=params.scvelo_mode, params=params
                )
                velocity_computed = True
                adata.uns["velocity_method"] = f"scvelo_{params.scvelo_mode}"
            except Exception as e:
                raise ProcessingError(
                    f"scVelo RNA velocity analysis failed: {e}"
                ) from e

    elif params.method == "velovi":
        try:
            velovi_results = await analyze_velocity_with_velovi(
                adata,
                params=params,
                n_epochs=params.velovi_n_epochs,
                n_hidden=params.velovi_n_hidden,
                n_latent=params.velovi_n_latent,
                n_layers=params.velovi_n_layers,
                dropout_rate=params.velovi_dropout_rate,
                learning_rate=params.velovi_learning_rate,
                use_gpu=params.velovi_use_gpu,
                ctx=ctx,
            )

            if velovi_results.get("velocity_computed", False):
                velocity_computed = True
                # Note: do NOT set adata.uns["velocity_graph"] here.
                # VELOVI does not produce a scVelo-compatible transition graph;
                # setting a boolean placeholder would mislead downstream tools
                # (scv.pl.velocity_embedding_stream, scv.tl.velocity_pseudotime)
                # that expect a sparse matrix.
                adata.uns["velocity_method"] = "velovi"
            else:
                raise ProcessingError("VELOVI failed to compute velocity")

        except ChatSpatialError:
            raise
        except Exception as e:
            raise ProcessingError(f"VELOVI velocity analysis failed: {e}") from e

    else:
        raise ParameterError(f"Unknown velocity method: {params.method}")

    # Build results keys based on what was computed
    # Note: velocity layers NOT exported (too large for CSV)
    results_keys: dict[str, list[str]] = {
        "uns": [],
        "obs": [],
        "obsm": [],
    }

    # Claim only results produced by the selected method in this run.
    if params.method == "velovi":
        if "velocity_velovi_norm" in adata.obs:
            results_keys["obs"].append("velocity_velovi_norm")
        if "X_velovi_latent" in adata.obsm:
            results_keys["obsm"].append("X_velovi_latent")
        if "velovi_latent_time" in adata.obsm:
            results_keys["obsm"].append("velovi_latent_time")

    # Only claim latent_time if THIS run computed it (scvelo dynamical only)
    if (
        params.method == "scvelo"
        and params.scvelo_mode == "dynamical"
        and "latent_time" in adata.obs
    ):
        results_keys["obs"].append("latent_time")

    # Store metadata for scientific provenance tracking
    analysis_key = _build_velocity_key(params)
    statistics: dict[str, Any] = {"velocity_computed": velocity_computed}
    if params.method == "velovi":
        for key in (
            "n_genes_analyzed",
            "original_n_genes",
            "velocity_mean_norm",
            "velocity_std_norm",
            "device",
        ):
            if key in velovi_results:
                statistics[key] = velovi_results[key]
    store_analysis_metadata(
        adata,
        analysis_name=analysis_key,
        method=params.method,
        parameters=_velocity_parameters_for_metadata(params),
        results_keys=results_keys,
        statistics=statistics,
    )

    # Export results for reproducibility
    export_analysis_result(adata, data_id, analysis_key)

    return RNAVelocityResult(
        data_id=data_id,
        velocity_computed=velocity_computed,
        velocity_graph_key=(
            "velocity_graph"
            if params.method == "scvelo" and "velocity_graph" in adata.uns
            else None
        ),
        mode=params.scvelo_mode if params.method == "scvelo" else params.method,
    )
