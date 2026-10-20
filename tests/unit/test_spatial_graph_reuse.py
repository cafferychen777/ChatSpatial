"""Regression tests: an explicit n_neighbors must not silently reuse a different spatial graph."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest
import squidpy as sq

from chatspatial.models.data import SpatialStatisticsParameters
from chatspatial.tools import spatial_statistics as spstats
from chatspatial.utils import compute

pytestmark = pytest.mark.unit


class DummyCtx:
    def __init__(self, adata):
        self.adata = adata
        self.warnings: list[str] = []

    async def get_adata(self, _data_id):
        return self.adata

    async def set_adata(self, _data_id, adata):
        self.adata = adata

    async def warning(self, msg):
        self.warnings.append(msg)


@pytest.fixture(autouse=True)
def no_export(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(spstats, "export_analysis_result", lambda *_a, **_k: [])


def _grid_adata() -> ad.AnnData:
    """A 12 x 12 lattice with one smooth and several noisy genes."""
    rng = np.random.default_rng(7)
    xs, ys = np.meshgrid(np.arange(12), np.arange(12))
    coords = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    smooth = np.sin(coords[:, 0] / 2.0) + rng.normal(0, 0.3, len(coords))
    X = np.column_stack([smooth] + [rng.normal(size=len(coords)) for _ in range(5)])
    adata = ad.AnnData(X.astype(np.float32))
    adata.var_names = ["SMOOTH"] + [f"N{i}" for i in range(5)]
    adata.obsm["spatial"] = coords
    adata.var["highly_variable"] = True
    return adata


def _with_grid_graph() -> ad.AnnData:
    adata = _grid_adata()
    # What compute_embeddings(spatial_coord_type="grid") leaves behind.
    sq.gr.spatial_neighbors(adata, coord_type="grid", n_rings=1)
    return adata


async def _moran(adata, **params):
    ctx = DummyCtx(adata.copy())
    result = await spstats.analyze_spatial_statistics(
        "d",
        ctx,
        SpatialStatisticsParameters(analysis_type="moran", genes=["SMOOTH"], **params),
    )
    return result, ctx


@pytest.mark.asyncio
async def test_explicit_n_neighbors_replaces_a_different_stored_graph():
    adata = _with_grid_graph()

    r6, _ = await _moran(adata, n_neighbors=6)
    r10, ctx10 = await _moran(adata, n_neighbors=10)

    i6 = r6.statistics["mean_morans_i"]
    i10 = r10.statistics["mean_morans_i"]
    assert i6 != pytest.approx(i10)
    assert r10.spatial_graph["coord_type"] == "generic"
    assert r10.spatial_graph["n_neighbors"] == 10
    assert r10.spatial_graph["source"] == "rebuilt"
    assert any("Replaced the stored spatial graph" in w for w in ctx10.warnings)
    # Each spot of a kNN graph has exactly k neighbors.
    graph = ctx10.adata.obsp["spatial_connectivities"]
    assert set(np.asarray((graph > 0).sum(axis=1)).ravel()) == {10}


@pytest.mark.asyncio
async def test_default_n_neighbors_reuses_the_stored_graph():
    adata = _with_grid_graph()

    result, ctx = await _moran(adata)

    assert result.spatial_graph["coord_type"] == "grid"
    assert result.spatial_graph["source"] == "existing"
    assert not ctx.warnings
    stored = adata.obsp["spatial_connectivities"]
    used = ctx.adata.obsp["spatial_connectivities"]
    assert (stored != used).nnz == 0


@pytest.mark.asyncio
async def test_explicit_n_neighbors_matching_the_stored_graph_reuses_it():
    adata = _grid_adata()
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=10)

    result, ctx = await _moran(adata, n_neighbors=10)

    assert result.spatial_graph["source"] == "existing"
    assert not ctx.warnings


@pytest.mark.asyncio
async def test_graph_is_computed_when_missing():
    result, _ = await _moran(_grid_adata(), n_neighbors=7)

    assert result.spatial_graph["source"] == "computed"
    assert result.spatial_graph["n_neighbors"] == 7


def test_ensure_spatial_neighbors_rebuild_flag():
    adata = _with_grid_graph()

    assert not compute.ensure_spatial_neighbors(adata, n_neighs=10)
    assert compute.spatial_graph_params(adata)["coord_type"] == "grid"

    assert compute.ensure_spatial_neighbors(
        adata, n_neighs=10, rebuild_if_different=True
    )
    assert compute.spatial_graph_params(adata)["n_neighbors"] == 10
    assert not compute.ensure_spatial_neighbors(
        adata, n_neighs=10, rebuild_if_different=True
    )

    # A graph without recorded parameters cannot be shown to match.
    del adata.uns["spatial_neighbors"]
    assert compute.spatial_graph_params(adata) is None
    assert compute.ensure_spatial_neighbors(
        adata, n_neighs=10, rebuild_if_different=True
    )
