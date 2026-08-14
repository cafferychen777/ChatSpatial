"""
Trajectory inference for spatial transcriptomics.

This module infers cellular trajectories and pseudotime by combining
expression patterns with optional velocity and spatial information.

Key functionality:
- `analyze_trajectory`: Main MCP entry point for trajectory analysis
- Supports CellRank (velocity-based), Palantir (expression-based), and DPT (diffusion-based)
"""

import logging
from typing import TYPE_CHECKING, Any, Optional

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

    from ..spatial_mcp_adapter import ToolContext

from ..models.analysis import TrajectoryResult
from ..models.data import TrajectoryParameters
from ..utils.adata_utils import (
    get_spatial_key,
    has_velovi_essential_data,
    reconstruct_velovi_adata,
    validate_obs_column,
)
from ..utils.compute import ensure_diffmap, ensure_neighbors, ensure_pca
from ..utils.dependency_manager import require, require_module
from ..utils.exceptions import (
    ChatSpatialError,
    DataCompatibilityError,
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)
from ..utils.mcp_utils import suppress_output

# Every method resolves one starting cell, and that choice sets the direction
# of the whole pseudotime, so it is recorded in one place regardless of method.
TRAJECTORY_ROOT_KEY = "trajectory_root_cell"


def _build_trajectory_key(params: "TrajectoryParameters") -> str:
    """Build a parametric analysis key for trajectory results.

    Encodes method + the parameter most likely to change between runs:
    - cellrank: spatial_weight (0=expression only, 1=spatial only)
    - palantir/dpt: first root cell (different roots → different trajectories)
    """
    method = params.method
    if method == "cellrank":
        sw = f"{params.spatial_weight:.2f}".replace(".", "_")
        return f"trajectory_cellrank_sw{sw}"
    if params.root_cells:
        root_tag = params.root_cells[0]
        return f"trajectory_{method}_{root_tag}"
    return f"trajectory_{method}"


def prepare_gam_model_for_visualization(
    adata: "ad.AnnData",
    genes: list[str],
    time_key: str = "latent_time",
    fate_key: str = "lineages_fwd",
) -> tuple[Any, list[str]]:
    """Prepare a GAM model for CellRank gene trends visualization.

    Handles the computation logic needed for CellRank 2.0 gene trends
    and fate heatmap visualizations. Requires data analyzed via
    analyze_rna_velocity (dynamical mode) and analyze_trajectory (cellrank).

    Args:
        adata: Annotated data matrix with CellRank results.
        genes: List of gene names to prepare the model for.
        time_key: Key in adata.obs for pseudotime/latent time values.
        fate_key: Key in adata.obsm for fate probabilities.

    Returns:
        Tuple of (GAM model, list of lineage names).

    Raises:
        DataNotFoundError: If fate probabilities or genes not found.
        DataError: If fate probabilities lack proper Lineage names.
    """
    # Validate required data
    validate_obs_column(adata, time_key, "Time")

    if fate_key not in adata.obsm:
        raise DataNotFoundError(
            f"Fate probabilities '{fate_key}' not found. Run analyze_trajectory first."
        )

    # Validate Lineage object has names
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

    cellrank_models = require_module(
        "cellrank",
        "cellrank.models",
        feature="CellRank gene trend visualization",
    )
    model = cellrank_models.GAM(adata)
    return model, lineage_names


def infer_spatial_trajectory_cellrank(
    adata: "ad.AnnData",
    spatial_weight: float = 0.5,
    kernel_weights: tuple[float, float] = (0.8, 0.2),
    n_states: int = 5,
    stability_threshold: float = 0.96,
) -> "ad.AnnData":
    """Infer cellular trajectories by combining RNA velocity with CellRank.

    Uses CellRank to model cell-state transitions by constructing
    a transition matrix from multiple kernels:
    1. A velocity kernel from RNA velocity
    2. A connectivity kernel based on transcriptomic similarity
    3. (Optional) A spatial kernel based on physical proximity

    Args:
        adata: AnnData with velocity data computed.
        spatial_weight: Weight for spatial kernel (0=no spatial, 1=spatial only).
        kernel_weights: Tuple of (velocity_weight, connectivity_weight).
        n_states: Number of macrostates to compute.

    Returns:
        AnnData with pseudotime, terminal states, and fate probabilities.

    Raises:
        ProcessingError: If CellRank computation fails.
    """
    cr = require("cellrank", feature="CellRank trajectory inference")
    from scipy.sparse import csr_matrix
    from scipy.spatial.distance import pdist, squareform

    # Check if spatial data is available
    spatial_key = get_spatial_key(adata)
    has_spatial = spatial_key is not None

    if not has_spatial and spatial_weight > 0:
        spatial_weight = 0

    # Handle different velocity methods
    if "velocity_method" in adata.uns and adata.uns["velocity_method"] == "velovi":
        # Reconstruct velovi adata from essential data stored in uns
        if not has_velovi_essential_data(adata):
            raise ProcessingError(
                "VELOVI velocity data not found. Run analyze_velocity_data first."
            )
        adata_for_cellrank = reconstruct_velovi_adata(adata)
        vk = cr.kernels.VelocityKernel(adata_for_cellrank)
        vk.compute_transition_matrix()
    else:
        adata_for_cellrank = adata
        vk = cr.kernels.VelocityKernel(adata_for_cellrank)
        vk.compute_transition_matrix()

    # Create connectivity kernel
    ck = cr.kernels.ConnectivityKernel(adata_for_cellrank)
    ck.compute_transition_matrix()

    # Combine kernels
    vk_weight, ck_weight = kernel_weights

    if has_spatial and spatial_weight > 0:
        spatial_coords = adata.obsm[spatial_key]
        spatial_dist = squareform(pdist(spatial_coords))
        dist_mean = spatial_dist.mean()
        if dist_mean < 1e-10:
            logging.getLogger(__name__).warning(
                "Spatial coordinates are nearly identical; disabling spatial kernel."
            )
            # Fall back to non-spatial kernel combination
            combined_kernel = vk_weight * vk + ck_weight * ck
        else:
            spatial_sim = np.exp(-spatial_dist / dist_mean)
            spatial_kernel = csr_matrix(spatial_sim)

            sk = cr.kernels.PrecomputedKernel(spatial_kernel, adata_for_cellrank)
            sk.compute_transition_matrix()

            combined_kernel = (1 - spatial_weight) * (
                vk_weight * vk + ck_weight * ck
            ) + spatial_weight * sk
    else:
        combined_kernel = vk_weight * vk + ck_weight * ck

    # GPCCA analysis
    g = cr.estimators.GPCCA(combined_kernel)
    g.compute_eigendecomposition()

    # Compute macrostates
    try:
        g.compute_macrostates(n_states=n_states)
    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(
            f"CellRank macrostate computation failed: {e}. "
            f"Try reducing n_states or use method='palantir'/'dpt'."
        ) from e

    # Terminal states are optional: when none qualify, the macrostate
    # fallback below still yields a pseudotime. CellRank reports this as a
    # ValueError whose wording depends on which criterion failed, so the
    # authoritative test is whether terminal states exist afterwards —
    # matching on the message turns a recoverable run into a hard failure
    # as soon as upstream rephrases it.
    try:
        g.predict_terminal_states(
            method="stability", stability_threshold=stability_threshold
        )
    except ValueError as e:
        logging.getLogger(__name__).warning(
            "CellRank selected no terminal states (%s). Falling back to "
            "macrostate-based pseudotime; lower "
            "cellrank_stability_threshold to select terminal states.",
            e,
        )

    # Check terminal states and compute fate probabilities
    has_terminal_states = (
        hasattr(g, "terminal_states")
        and g.terminal_states is not None
        and len(g.terminal_states.cat.categories) > 0
    )

    if has_terminal_states:
        try:
            g.compute_fate_probabilities()
            absorption_probs = g.fate_probabilities

            # Derive pseudotime from fate probability entropy:
            # high entropy = multipotent/undifferentiated = early
            # low entropy = committed to one fate = late
            fate_matrix = np.asarray(absorption_probs)

            # Check for NaN in fate probabilities
            nan_mask = np.isnan(fate_matrix).any(axis=1)
            if nan_mask.all():
                raise ProcessingError(
                    "CellRank fate probabilities are all NaN. "
                    "This indicates numerical instability in "
                    "GPCCA. Try reducing n_states or adjusting "
                    "kernel weights."
                )

            fate_matrix = np.clip(fate_matrix, 1e-10, None)  # avoid log(0)
            entropy = -np.sum(fate_matrix * np.log(fate_matrix), axis=1)

            # NaN cells get pseudotime = NaN (not 0)
            max_entropy = np.nanmax(entropy)
            if max_entropy > 0:
                pseudotime = 1 - entropy / max_entropy
            else:
                pseudotime = np.zeros_like(entropy)

            if nan_mask.any():
                # Cells without fate probabilities get no pseudotime; the
                # caller reports the count from the pseudotime column.
                pseudotime[nan_mask] = np.nan

            adata_for_cellrank.obs["pseudotime"] = pseudotime
            adata_for_cellrank.obsm["fate_probabilities"] = absorption_probs
            adata_for_cellrank.obs["terminal_states"] = g.terminal_states
        except Exception as e:
            raise ProcessingError(
                f"CellRank fate probability computation failed: {e}. "
                f"This often indicates numerical instability. "
                f"Try method='palantir' or 'dpt' instead."
            ) from e
    else:
        # Fall back to macrostates-based pseudotime. Guard on the
        # memberships actually read, not on the macrostates attribute.
        macrostate_probs = getattr(g, "macrostates_memberships", None)
        if macrostate_probs is not None:
            pseudotime = 1 - macrostate_probs[:, 0].X.flatten()
            adata_for_cellrank.obs["pseudotime"] = pseudotime
        else:
            raise ProcessingError(
                "CellRank could not compute terminal states or macrostates. "
                "Try method='palantir' or 'dpt' instead."
            )

    if hasattr(g, "macrostates") and g.macrostates is not None:
        adata_for_cellrank.obs["macrostates"] = g.macrostates

    # Transfer results back to original adata
    if "pseudotime" in adata_for_cellrank.obs:
        adata.obs["pseudotime"] = adata_for_cellrank.obs["pseudotime"]
    if "terminal_states" in adata_for_cellrank.obs:
        adata.obs["terminal_states"] = adata_for_cellrank.obs["terminal_states"]
    if "macrostates" in adata_for_cellrank.obs:
        adata.obs["macrostates"] = adata_for_cellrank.obs["macrostates"]
    if "fate_probabilities" in adata_for_cellrank.obsm:
        adata.obsm["fate_probabilities"] = adata_for_cellrank.obsm["fate_probabilities"]
        # Also write CellRank-standard alias so viz doesn't need to
        adata.obsm["to_terminal_states"] = adata.obsm["fate_probabilities"]
    if (
        "to_terminal_states" in adata_for_cellrank.obsm
        and "to_terminal_states" not in adata.obsm
    ):
        adata.obsm["to_terminal_states"] = adata_for_cellrank.obsm["to_terminal_states"]

    # Note: With optimized storage, velovi data is stored as individual arrays
    # in uns (velovi_velocity, velovi_Ms, etc.) rather than a full adata copy.
    # Results are already transferred to original adata above.

    return adata


def _normalize_palantir_matrix(
    matrix: object,
    obs_names: pd.Index,
    *,
    name: str,
    allow_empty: bool = False,
) -> pd.DataFrame:
    """Validate and align a Palantir cell-by-feature matrix.

    ``allow_empty`` distinguishes a matrix with no columns that is a corrupt
    result (an empty diffusion space, which nothing downstream can use) from
    one that is a legitimate finding (no terminal states, which still leaves
    a usable pseudotime).
    """
    expected_index = pd.Index(obs_names)
    if not expected_index.is_unique:
        raise DataCompatibilityError(
            "Palantir requires unique observation names for result alignment."
        )

    if isinstance(matrix, pd.DataFrame):
        if not matrix.index.is_unique:
            raise DataCompatibilityError(
                f"Palantir {name} contains duplicate cell identifiers."
            )
        missing = expected_index.difference(matrix.index)
        unexpected = matrix.index.difference(expected_index)
        if len(missing) or len(unexpected):
            raise DataCompatibilityError(
                f"Palantir {name} cells do not match the dataset: "
                f"{len(missing)} missing and {len(unexpected)} unexpected."
            )
        normalized = matrix.reindex(expected_index)
    else:
        try:
            values = np.asarray(matrix, dtype=float)
        except (TypeError, ValueError) as exc:
            raise DataCompatibilityError(f"Palantir {name} must be numeric.") from exc
        if values.ndim != 2 or values.shape[0] != len(expected_index):
            raise DataCompatibilityError(
                f"Palantir {name} has shape {values.shape}; expected "
                f"({len(expected_index)}, n_features)."
            )
        normalized = pd.DataFrame(values, index=expected_index)

    if normalized.shape[1] == 0 and not allow_empty:
        raise DataCompatibilityError(f"Palantir {name} contains no components.")
    try:
        values = normalized.to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise DataCompatibilityError(f"Palantir {name} must be numeric.") from exc
    if not np.isfinite(values).all():
        raise DataCompatibilityError(f"Palantir {name} contains non-finite values.")
    return pd.DataFrame(values, index=expected_index, columns=normalized.columns)


def _normalize_palantir_pseudotime(
    pseudotime: object, obs_names: pd.Index
) -> pd.Series:
    """Validate and align Palantir pseudotime to the AnnData observation axis."""
    expected_index = pd.Index(obs_names)
    if isinstance(pseudotime, pd.Series):
        if not pseudotime.index.is_unique:
            raise DataCompatibilityError(
                "Palantir pseudotime contains duplicate cell identifiers."
            )
        missing = expected_index.difference(pseudotime.index)
        unexpected = pseudotime.index.difference(expected_index)
        if len(missing) or len(unexpected):
            raise DataCompatibilityError(
                "Palantir pseudotime cells do not match the dataset."
            )
        aligned = pseudotime.reindex(expected_index)
    else:
        try:
            values = np.asarray(pseudotime, dtype=float)
        except (TypeError, ValueError) as exc:
            raise DataCompatibilityError(
                "Palantir pseudotime must be numeric."
            ) from exc
        if values.ndim != 1 or len(values) != len(expected_index):
            raise DataCompatibilityError(
                "Palantir pseudotime must contain exactly one value per cell."
            )
        aligned = pd.Series(values, index=expected_index)

    try:
        values = aligned.to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise DataCompatibilityError("Palantir pseudotime must be numeric.") from exc
    if (
        not np.isfinite(values).all()
        or (values < -1e-10).any()
        or (values > 1 + 1e-10).any()
    ):
        raise DataCompatibilityError(
            "Palantir pseudotime must contain finite values in [0, 1]."
        )
    return pd.Series(np.clip(values, 0, 1), index=expected_index)


def infer_pseudotime_palantir(
    adata: "ad.AnnData",
    root_cells: Optional[list[str]] = None,
    n_diffusion_components: int = 10,
    num_waypoints: int = 500,
    ctx: Optional["ToolContext"] = None,
    palantir_module: Any | None = None,
) -> "ad.AnnData":
    """Infer cellular trajectories and pseudotime using Palantir.

    Palantir models differentiation as a stochastic process on a graph,
    using diffusion maps to capture data geometry and computing fate
    probabilities via random walks from a root cell.

    Args:
        adata: Annotated data matrix with PCA results.
        root_cells: Cell identifiers as starting points. Auto-selected if None.
        n_diffusion_components: Number of diffusion components.
        num_waypoints: Number of waypoints for trajectory granularity.

    Returns:
        AnnData with palantir_pseudotime and palantir_branch_probs.

    Raises:
        ParameterError: If specified root cell not found in data.
    """
    if root_cells and root_cells[0] not in adata.obs_names:
        raise ParameterError(f"Root cell '{root_cells[0]}' not found in data")

    ensure_pca(adata)
    palantir = (
        palantir_module
        if palantir_module is not None
        else require("palantir", ctx, feature="Palantir trajectory inference")
    )

    pca_df = pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names)
    dm_res = palantir.utils.run_diffusion_maps(
        pca_df, n_components=n_diffusion_components
    )
    ms_data = _normalize_palantir_matrix(
        palantir.utils.determine_multiscale_space(dm_res),
        adata.obs_names,
        name="multiscale diffusion space",
    )

    if root_cells is not None and len(root_cells) > 0:
        start_cell = root_cells[0]
    else:
        # Sign-invariant: pick cell with largest absolute value in first DC
        start_cell = ms_data.iloc[:, 0].abs().idxmax()
    adata.uns[TRAJECTORY_ROOT_KEY] = str(start_cell)

    pr_res = palantir.core.run_palantir(
        ms_data, start_cell, num_waypoints=num_waypoints
    )

    pseudotime = _normalize_palantir_pseudotime(pr_res.pseudotime, adata.obs_names)
    branch_probs = _normalize_palantir_matrix(
        pr_res.branch_probs,
        adata.obs_names,
        name="branch probabilities",
        allow_empty=True,
    )
    branch_values = branch_probs.to_numpy()
    if (branch_values < -1e-10).any() or (branch_values > 1 + 1e-10).any():
        raise DataCompatibilityError(
            "Palantir branch probabilities must lie in [0, 1]."
        )

    adata.obs["palantir_pseudotime"] = pseudotime
    # No terminal states is a statement about the data, not a failed run:
    # pseudotime is the primary output and stays valid without fate
    # probabilities, so the empty matrix is left unwritten and reported.
    if branch_probs.shape[1] > 0:
        adata.obsm["palantir_branch_probs"] = branch_probs.clip(0, 1)

    return adata


def compute_dpt_trajectory(
    adata: "ad.AnnData",
    root_cells: Optional[list[str]] = None,
) -> "ad.AnnData":
    """Compute Diffusion Pseudotime (DPT) trajectory analysis.

    Args:
        adata: Annotated data matrix (will compute PCA/neighbors/diffmap if needed).
        root_cells: Cell identifiers as starting points. Uses first cell if None.

    Returns:
        AnnData with dpt_pseudotime in obs.

    Raises:
        ParameterError: If specified root cell not found.
        ProcessingError: If DPT computation fails.
    """
    import scanpy as sc

    ensure_pca(adata)
    ensure_neighbors(adata)
    ensure_diffmap(adata)

    if root_cells is not None and len(root_cells) > 0:
        if root_cells[0] in adata.obs_names:
            adata.uns["iroot"] = np.where(adata.obs_names == root_cells[0])[0][0]
        else:
            raise ParameterError(
                f"Root cell '{root_cells[0]}' not found. "
                f"Use valid cell ID from adata.obs_names or omit to auto-select."
            )
    else:
        # Use first diffusion component for principled root selection
        # (consistent with Palantir heuristic, sign-invariant)
        if "X_diffmap" in adata.obsm:
            dc1 = adata.obsm["X_diffmap"][:, 0]
            adata.uns["iroot"] = int(np.argmax(np.abs(dc1)))
        else:
            adata.uns["iroot"] = 0
    adata.uns[TRAJECTORY_ROOT_KEY] = str(adata.obs_names[adata.uns["iroot"]])

    try:
        sc.tl.dpt(adata)
    except Exception as e:
        raise ProcessingError(f"DPT computation failed: {e}") from e

    if "dpt_pseudotime" not in adata.obs.columns:
        raise ProcessingError("DPT computation did not create 'dpt_pseudotime' column")

    # NaN pseudotime is preserved rather than filled: those cells are
    # unreachable from the root, which the caller reports.
    return adata


def has_velocity_data(adata: "ad.AnnData") -> bool:
    """Check if RNA velocity has been computed (by any method).

    Args:
        adata: Annotated data matrix to check for velocity data.

    Returns:
        True if velocity data exists (scVelo or VELOVI), False otherwise.
    """
    return (
        "velocity_graph" in adata.uns
        or "velocity_method" in adata.uns
        or has_velovi_essential_data(adata)
    )


async def analyze_trajectory(
    data_id: str,
    ctx: "ToolContext",
    params: TrajectoryParameters | None = None,
) -> TrajectoryResult:
    """
    Analyze trajectory and cell state transitions in spatial transcriptomics data.

    This is the main MCP entry point for trajectory inference. It supports:
    - CellRank: Requires pre-computed velocity data
    - Palantir: Expression-based, no velocity required
    - DPT: Diffusion-based, no velocity required

    Args:
        data_id: Dataset identifier.
        ctx: ToolContext for data access and logging.
        params: Trajectory analysis parameters.

    Returns:
        TrajectoryResult with pseudotime and method metadata.
    """
    if params is None:
        params = TrajectoryParameters()

    source_adata = await ctx.get_adata(data_id)

    velocity_available = has_velocity_data(source_adata)
    pseudotime_key = None
    method_used = params.method

    # All supported backends mutate AnnData incrementally. Run the complete
    # workflow on a candidate so failures cannot leave stale pseudotime,
    # transition states, or metadata in the managed source dataset.
    adata = source_adata.copy()

    # Execute requested method
    if params.method == "cellrank":
        if not velocity_available:
            raise ProcessingError(
                "CellRank requires velocity data. Run velocity analysis first or use palantir/dpt."
            )

        try:
            with suppress_output():
                adata = infer_spatial_trajectory_cellrank(
                    adata,
                    spatial_weight=params.spatial_weight,
                    kernel_weights=(
                        params.cellrank_kernel_weights[0],
                        params.cellrank_kernel_weights[1],
                    ),
                    n_states=params.cellrank_n_states,
                    stability_threshold=params.cellrank_stability_threshold,
                )
            pseudotime_key = "pseudotime"
            method_used = "cellrank"

            # Per-run suffix for CellRank result coexistence
            sw_str = f"{params.spatial_weight:.2f}".replace(".", "_")
            cellrank_suffix = f"cellrank_sw{sw_str}"

            # Save per-run copies for provenance (shared keys kept for viz)
            if "pseudotime" in adata.obs.columns:
                adata.obs[f"pseudotime_{cellrank_suffix}"] = adata.obs["pseudotime"]
            for _ck in ("terminal_states", "macrostates"):
                if _ck in adata.obs.columns:
                    adata.obs[f"{_ck}_{cellrank_suffix}"] = adata.obs[_ck]
            if "fate_probabilities" in adata.obsm:
                adata.obsm[f"fate_probabilities_{cellrank_suffix}"] = adata.obsm[
                    "fate_probabilities"
                ]
        except ChatSpatialError:
            raise
        except Exception as e:
            raise ProcessingError(f"CellRank trajectory inference failed: {e}") from e

    elif params.method == "palantir":
        try:
            with suppress_output():
                adata = infer_pseudotime_palantir(
                    adata,
                    root_cells=params.root_cells,
                    n_diffusion_components=params.palantir_n_diffusion_components,
                    num_waypoints=params.palantir_n_waypoints,
                    ctx=ctx,
                )

            pseudotime_key = "palantir_pseudotime"
            method_used = "palantir"

        except ChatSpatialError:
            raise
        except Exception as e:
            raise ProcessingError(f"Palantir trajectory inference failed: {e}") from e

    elif params.method == "dpt":
        try:
            with suppress_output():
                adata = compute_dpt_trajectory(adata, root_cells=params.root_cells)
            pseudotime_key = "dpt_pseudotime"
            method_used = "dpt"
        except ChatSpatialError:
            raise
        except Exception as e:
            raise ProcessingError(f"DPT analysis failed: {e}") from e

    else:
        raise ParameterError(f"Unknown trajectory method: {params.method}")

    if pseudotime_key is None or pseudotime_key not in adata.obs.columns:
        raise ProcessingError("Failed to compute pseudotime with any available method")

    # The root cell orders the whole trajectory, so an automatic choice is a
    # result the caller has to see, not a log line.
    root_cell = adata.uns.get(TRAJECTORY_ROOT_KEY)
    if root_cell is not None and not params.root_cells:
        await ctx.warning(
            f"No root_cells given, so {method_used} started from '{root_cell}', "
            "chosen as the extreme of the first diffusion component. Pseudotime "
            "is measured from that cell; pass root_cells to set the origin."
        )

    n_unreachable = int(adata.obs[pseudotime_key].isna().sum())
    if n_unreachable:
        await ctx.warning(
            f"{n_unreachable} cells have no pseudotime because they are "
            "unreachable from the root. They are kept as NaN rather than filled."
        )

    if method_used == "palantir" and "palantir_branch_probs" not in adata.obsm:
        await ctx.warning(
            "Palantir identified no terminal states, so pseudotime is reported "
            "without fate probabilities and fate visualizations are unavailable. "
            "Tissue without distinct differentiation endpoints gives this result; "
            "pass root_cells if the automatic root sits mid-trajectory."
        )

    # Store scientific metadata
    from ..utils.adata_utils import store_analysis_metadata
    from ..utils.results_export import export_analysis_result

    results_keys_dict: dict[str, Any] = {"obs": [pseudotime_key], "obsm": [], "uns": []}

    # For CellRank, use per-run keys in results_keys for provenance accuracy
    if method_used == "cellrank":
        sw_str = f"{params.spatial_weight:.2f}".replace(".", "_")
        cellrank_suffix = f"cellrank_sw{sw_str}"
        results_keys_dict["obs"] = [f"pseudotime_{cellrank_suffix}"]

    if method_used == "cellrank":
        for obs_key in ("terminal_states", "macrostates"):
            suffixed = f"{obs_key}_{cellrank_suffix}"
            if suffixed in adata.obs.columns:
                results_keys_dict["obs"].append(suffixed)
        fate_suffixed = f"fate_probabilities_{cellrank_suffix}"
        if fate_suffixed in adata.obsm:
            results_keys_dict["obsm"].append(fate_suffixed)
    elif method_used == "palantir":
        if "palantir_branch_probs" in adata.obsm:
            results_keys_dict["obsm"].append("palantir_branch_probs")
    elif method_used == "dpt":
        results_keys_dict["uns"].append("iroot")

    if TRAJECTORY_ROOT_KEY in adata.uns:
        results_keys_dict["uns"].append(TRAJECTORY_ROOT_KEY)

    parameters_dict: dict[str, Any] = {"spatial_weight": params.spatial_weight}
    if method_used == "cellrank":
        parameters_dict.update(
            {
                "kernel_weights": list(params.cellrank_kernel_weights),
                "n_states": params.cellrank_n_states,
            }
        )
    elif method_used == "palantir":
        parameters_dict.update(
            {
                "n_diffusion_components": params.palantir_n_diffusion_components,
                "num_waypoints": params.palantir_n_waypoints,
            }
        )

    if params.root_cells:
        parameters_dict["root_cells"] = params.root_cells

    statistics_dict = {
        "velocity_computed": velocity_available,
        "pseudotime_key": pseudotime_key,
    }

    analysis_key = _build_trajectory_key(params)
    store_analysis_metadata(
        adata,
        analysis_name=analysis_key,
        method=method_used,
        parameters=parameters_dict,
        results_keys=results_keys_dict,
        statistics=statistics_dict,
    )

    # Export results for reproducibility
    export_analysis_result(adata, data_id, analysis_key)

    result = TrajectoryResult(
        data_id=data_id,
        pseudotime_computed=True,
        velocity_computed=velocity_available,
        pseudotime_key=pseudotime_key,
        method=method_used,
        spatial_weight=params.spatial_weight,
    )
    await ctx.set_adata(data_id, adata)
    return result
