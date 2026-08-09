"""Integration contracts for tools.velocity and tools.trajectory entrypoints."""

from __future__ import annotations

import numpy as np
import pytest

from chatspatial.models.analysis import RNAVelocityResult, TrajectoryResult
from chatspatial.models.data import RNAVelocityParameters, TrajectoryParameters
from chatspatial.tools import trajectory as traj_module
from chatspatial.tools import velocity as vel_module
from chatspatial.tools.trajectory import analyze_trajectory
from chatspatial.tools.velocity import analyze_rna_velocity
from chatspatial.utils.exceptions import ProcessingError


class DummyCtx:
    def __init__(self, adata):
        self._adata = adata

    async def get_adata(self, data_id: str):
        return self._adata

    async def set_adata(self, data_id: str, adata):
        self._adata = adata

    async def info(self, msg: str):
        return None

    async def warning(self, msg: str):
        return None


@pytest.mark.integration
@pytest.mark.asyncio
async def test_analyze_rna_velocity_scvelo_branch_returns_expected_contract(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = adata.X.copy()
    adata.layers["unspliced"] = adata.X.copy()
    ctx = DummyCtx(adata)

    fake_scv = object()
    monkeypatch.setattr(
        vel_module,
        "require",
        lambda *args, **kwargs: fake_scv,
    )
    monkeypatch.setattr(vel_module, "validate_adata", lambda *args, **kwargs: None)

    def _fake_compute(adata, mode, params, *, scv):
        assert scv is fake_scv
        # Simulate scv.tl.velocity_graph() output (sparse transition matrix)
        import scipy.sparse

        adata.uns["velocity_graph"] = scipy.sparse.eye(adata.n_obs, format="csr")
        return adata

    monkeypatch.setattr(vel_module, "compute_rna_velocity", _fake_compute)
    monkeypatch.setattr(
        vel_module, "store_analysis_metadata", lambda *args, **kwargs: None
    )
    monkeypatch.setattr(
        vel_module, "export_analysis_result", lambda *args, **kwargs: None
    )

    params = RNAVelocityParameters(method="scvelo", scvelo_mode="stochastic")
    result = await analyze_rna_velocity("d1", ctx, params)

    assert isinstance(result, RNAVelocityResult)
    assert result.velocity_computed is True
    assert result.mode == "stochastic"
    assert result.velocity_graph_key == "velocity_graph"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_analyze_rna_velocity_velovi_branch_sets_uns_method(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = adata.X.copy()
    adata.layers["unspliced"] = adata.X.copy()
    ctx = DummyCtx(adata)

    monkeypatch.setattr(vel_module, "require", lambda *args, **kwargs: None)
    monkeypatch.setattr(vel_module, "validate_adata", lambda *args, **kwargs: None)

    async def fake_velovi(*args, **kwargs):
        return {"velocity_computed": True}

    monkeypatch.setattr(vel_module, "analyze_velocity_with_velovi", fake_velovi)
    monkeypatch.setattr(
        vel_module, "store_analysis_metadata", lambda *args, **kwargs: None
    )
    monkeypatch.setattr(
        vel_module, "export_analysis_result", lambda *args, **kwargs: None
    )

    params = RNAVelocityParameters(method="velovi")
    result = await analyze_rna_velocity("d2", ctx, params)

    assert isinstance(result, RNAVelocityResult)
    assert result.mode == "velovi"
    assert ctx._adata is not adata
    assert ctx._adata.uns["velocity_method"] == "velovi"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_analyze_trajectory_cellrank_requires_velocity_data(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    params = TrajectoryParameters(method="cellrank")
    with pytest.raises(ProcessingError, match="requires velocity data"):
        await analyze_trajectory("d3", ctx, params)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_analyze_trajectory_palantir_and_dpt_dispatch_contract(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    def fake_palantir(
        adata, root_cells, n_diffusion_components, num_waypoints, **_kwargs
    ):
        adata.obs["palantir_pseudotime"] = 0.1
        adata.obsm["palantir_branch_probs"] = np.ones((adata.n_obs, 1), dtype=float)
        return adata

    def fake_dpt(adata, root_cells):
        adata.obs["dpt_pseudotime"] = 0.2
        adata.uns["iroot"] = 0
        return adata

    monkeypatch.setattr(traj_module, "infer_pseudotime_palantir", fake_palantir)
    monkeypatch.setattr(traj_module, "compute_dpt_trajectory", fake_dpt)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *args, **kwargs: None,
    )

    palantir_result = await analyze_trajectory(
        "d4", ctx, TrajectoryParameters(method="palantir")
    )
    dpt_result = await analyze_trajectory("d4", ctx, TrajectoryParameters(method="dpt"))

    assert isinstance(palantir_result, TrajectoryResult)
    assert palantir_result.method == "palantir"
    assert palantir_result.pseudotime_key == "palantir_pseudotime"

    assert isinstance(dpt_result, TrajectoryResult)
    assert dpt_result.method == "dpt"
    assert dpt_result.pseudotime_key == "dpt_pseudotime"
