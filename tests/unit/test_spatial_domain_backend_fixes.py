"""Regression tests for spatial-domain backend behaviour that the schema promises.

Covers: a STAGATE radius that follows coordinate units, GraphST seeding inside
its training worker, GraphST Louvain running Louvain, cluster_n_neighbors being
honoured over a stored graph, per-method timeouts, and schema bounds.
"""

from __future__ import annotations

import sys
import types

import anndata as ad
import numpy as np
import pytest
import scanpy as sc
from pydantic import ValidationError

from chatspatial.models.data import SpatialDomainParameters
from chatspatial.tools import spatial_domains as sd
from chatspatial.utils.compute import ensure_neighbors
from chatspatial.utils.exceptions import DependencyError

pytestmark = pytest.mark.unit


def _visium_array_grid(n_rows: int = 30, n_cols: int = 60) -> np.ndarray:
    """Visium array indices: spots at (row, col) with row + col even."""
    return np.array(
        [(c, r) for r in range(n_rows) for c in range(n_cols) if (r + c) % 2 == 0],
        dtype=float,
    )


def _visium_pixels(spacing: float = 138.0) -> np.ndarray:
    """The same lattice in physical (pixel) units: a regular hexagonal grid."""
    grid = _visium_array_grid()
    return np.column_stack(
        [grid[:, 0] * spacing / 2, grid[:, 1] * spacing * 3**0.5 / 2]
    )


class _Ctx:
    def __init__(self) -> None:
        self.warnings: list[str] = []

    async def warning(self, msg: str) -> None:
        self.warnings.append(msg)

    def debug(self, _msg: str) -> None:
        return None


async def _cpu(prefer_gpu: bool, ctx) -> str:
    del prefer_gpu, ctx
    return "cpu"


async def _inline(function, **_kwargs):
    return function()


def _adata_with_coords(coords: np.ndarray) -> ad.AnnData:
    rng = np.random.default_rng(0)
    adata = ad.AnnData(rng.poisson(3, size=(coords.shape[0], 8)).astype(np.float32))
    adata.obsm["spatial"] = coords
    return adata


# ---------------------------------------------------------------- STAGATE radius


@pytest.mark.parametrize("scale", [1.0, 138.0 / 2, 1000.0])
def test_stagate_auto_radius_links_first_ring_in_any_unit(scale: float) -> None:
    pixels = _visium_pixels() / 138.0 * scale
    radius = sd._stagate_auto_radius(pixels)
    mean = sd._mean_radius_neighbors(pixels, radius)
    # Interior spots have exactly the six first-ring neighbours; edges fewer.
    assert 5.0 <= mean <= 6.0
    assert np.isclose(radius, scale * 3**0.25, rtol=1e-6)


def test_stagate_auto_radius_on_array_indices_stays_local() -> None:
    grid = _visium_array_grid()
    radius = sd._stagate_auto_radius(grid)
    mean = sd._mean_radius_neighbors(grid, radius)
    # The old fixed radius of 50 linked thousands of spots on this grid.
    assert sd._mean_radius_neighbors(grid, 50.0) > 500
    assert 4.0 <= mean <= 8.0


def _fake_stagate(captured: dict):
    class _FakeSTAGATE:
        @staticmethod
        def Cal_Spatial_Net(a, rad_cutoff=None):
            captured["rad_cutoff"] = rad_cutoff

        @staticmethod
        def Stats_Spatial_Net(_a):
            return None

        @staticmethod
        def train_STAGATE(a, device=None, random_seed=0):
            del device
            captured["random_seed"] = random_seed
            a.obsm["STAGATE"] = np.random.default_rng(0).normal(size=(a.n_obs, 4))
            return a

    return _FakeSTAGATE


def _patch_stagate(monkeypatch: pytest.MonkeyPatch, captured: dict) -> None:
    torch_mod = types.ModuleType("torch")
    torch_mod.device = lambda x: x
    deps = {"torch": torch_mod, "STAGATE_pyG": _fake_stagate(captured)}
    monkeypatch.setattr(sd, "require", lambda name, *_a, **_k: deps[name])
    monkeypatch.setattr(sd, "resolve_device_async", _cpu)
    monkeypatch.setattr(sd, "run_sync_with_timeout", _inline)


@pytest.mark.asyncio
@pytest.mark.parametrize("coords", [_visium_array_grid(), _visium_pixels()])
async def test_stagate_default_radius_follows_coordinate_units(
    monkeypatch: pytest.MonkeyPatch, coords: np.ndarray
) -> None:
    captured: dict = {}
    _patch_stagate(monkeypatch, captured)
    ctx = _Ctx()
    _, _, stats = await sd._identify_domains_stagate(
        _adata_with_coords(coords),
        SpatialDomainParameters(method="stagate", n_domains=3, stagate_random_seed=7),
        ctx,
    )
    assert stats["rad_cutoff_source"] == "auto (spot spacing)"
    assert captured["rad_cutoff"] == pytest.approx(stats["rad_cutoff"])
    assert 4.0 <= stats["mean_neighbors"] <= 8.0
    assert captured["random_seed"] == 7
    assert ctx.warnings == []


@pytest.mark.asyncio
async def test_stagate_explicit_radius_is_kept_and_degenerate_graph_warns(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict = {}
    _patch_stagate(monkeypatch, captured)
    ctx = _Ctx()
    _, _, stats = await sd._identify_domains_stagate(
        _adata_with_coords(_visium_pixels()),
        SpatialDomainParameters(method="stagate", n_domains=3, stagate_rad_cutoff=50),
        ctx,
    )
    assert captured["rad_cutoff"] == 50
    assert stats["rad_cutoff_source"] == "user"
    assert stats["mean_neighbors"] == 0
    assert any("neighbours per spot" in w for w in ctx.warnings)


# ---------------------------------------------------------------- GraphST


def _install_fake_graphst(monkeypatch: pytest.MonkeyPatch, seeds: list) -> None:
    class _FakeGraphST:
        def __init__(self, adata_graphst, device=None, random_seed=0):
            del device, random_seed
            self._adata = adata_graphst

        def train(self):
            # Mimics GraphST: weights and per-epoch permutations use the global
            # NumPy generator, which is unseeded in a spawned worker.
            self._adata.obsm["emb"] = np.random.normal(size=(self._adata.n_obs, 24))
            return self._adata

    sub = types.ModuleType("GraphST.GraphST")
    sub.GraphST = _FakeGraphST
    pkg = types.ModuleType("GraphST")
    pkg.GraphST = sub
    pre = types.ModuleType("GraphST.preprocess")

    def _fix_seed(seed):
        seeds.append(seed)
        np.random.seed(seed)

    pre.fix_seed = _fix_seed
    monkeypatch.setitem(sys.modules, "GraphST", pkg)
    monkeypatch.setitem(sys.modules, "GraphST.GraphST", sub)
    monkeypatch.setitem(sys.modules, "GraphST.preprocess", pre)
    torch_mod = types.ModuleType("torch")
    torch_mod.device = lambda x: x
    monkeypatch.setitem(sys.modules, "torch", torch_mod)
    monkeypatch.setattr(
        sd,
        "require_module",
        lambda _n, module_name, *_a, **_k: sys.modules[module_name],
    )
    monkeypatch.setattr(sd, "resolve_device_async", _cpu)


@pytest.mark.asyncio
async def test_graphst_training_worker_is_seeded_and_repeatable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seeds: list = []
    _install_fake_graphst(monkeypatch, seeds)
    monkeypatch.setattr(sd, "require", lambda name, *_a, **_k: sys.modules[name])

    async def _unseeded_worker(function, **_kwargs):
        np.random.seed(None)  # a fresh process starts from OS entropy
        return function()

    monkeypatch.setattr(sd, "run_sync_with_timeout", _unseeded_worker)
    params = SpatialDomainParameters(
        method="graphst",
        n_domains=3,
        graphst_clustering_method="mclust",
        graphst_refinement=False,
        graphst_random_seed=11,
    )
    coords = _visium_pixels()[:200]
    runs = []
    for _ in range(2):
        labels, _, stats = await sd._identify_domains_graphst(
            _adata_with_coords(coords), params, _Ctx()
        )
        runs.append(labels.tolist())
    assert seeds == [11, 11]
    assert runs[0] == runs[1]
    assert stats["random_seed"] == 11


@pytest.mark.asyncio
async def test_graphst_louvain_runs_louvain_not_leiden(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_graphst(monkeypatch, [])
    monkeypatch.setattr(sd, "require", lambda name, *_a, **_k: sys.modules.get(name))
    monkeypatch.setattr(sd, "run_sync_with_timeout", _inline)
    monkeypatch.setattr(sd.sc.pp, "neighbors", lambda *_a, **_k: None)
    called: list[str] = []

    def _louvain(a, resolution=1.0, random_state=0):
        del resolution, random_state
        called.append("louvain")
        a.obs["louvain"] = [str(i % 3) for i in range(a.n_obs)]

    def _leiden(*_a, **_k):
        raise AssertionError(
            "Leiden must not run for graphst_clustering_method='louvain'"
        )

    monkeypatch.setattr(sd.sc.tl, "louvain", _louvain)
    monkeypatch.setattr(sd.sc.tl, "leiden", _leiden)
    labels, _, stats = await sd._identify_domains_graphst(
        _adata_with_coords(_visium_pixels()[:90]),
        SpatialDomainParameters(
            method="graphst",
            n_domains=3,
            graphst_clustering_method="louvain",
            graphst_refinement=False,
        ),
        _Ctx(),
    )
    assert called and set(called) == {"louvain"}
    assert stats["clustering_method"] == "louvain"
    assert labels.nunique() == 3


@pytest.mark.asyncio
async def test_graphst_louvain_without_package_fails_before_training(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_graphst(monkeypatch, [])

    def _require(name, *_a, **_k):
        if name == "louvain":
            raise DependencyError("louvain is not installed")
        return sys.modules[name]

    async def _no_training(*_a, **_k):
        raise AssertionError("training must not start")

    monkeypatch.setattr(sd, "require", _require)
    monkeypatch.setattr(sd, "run_sync_with_timeout", _no_training)
    with pytest.raises(DependencyError, match="louvain"):
        await sd._identify_domains_graphst(
            _adata_with_coords(_visium_pixels()[:90]),
            SpatialDomainParameters(
                method="graphst", graphst_clustering_method="louvain"
            ),
            _Ctx(),
        )


# ---------------------------------------------------------------- neighbours


def _pca_adata() -> ad.AnnData:
    rng = np.random.default_rng(1)
    adata = ad.AnnData(rng.normal(size=(80, 10)).astype(np.float32))
    adata.obsm["X_pca"] = rng.normal(size=(80, 5))
    return adata


def test_ensure_neighbors_rebuilds_graph_when_k_is_required_and_differs() -> None:
    adata = _pca_adata()
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
    assert ensure_neighbors(adata, n_neighbors=30) is False  # default: reuse
    assert ensure_neighbors(adata, n_neighbors=30, require_n_neighbors=True) is True
    assert adata.uns["neighbors"]["params"]["n_neighbors"] == 30
    assert ensure_neighbors(adata, n_neighbors=30, require_n_neighbors=True) is False


@pytest.mark.asyncio
async def test_cluster_n_neighbors_is_honoured_over_stored_graph(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    adata = _pca_adata()
    adata.obsm["spatial"] = _visium_pixels()[:80]
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")

    class _Sq:
        class gr:
            @staticmethod
            def spatial_neighbors(a, **_k):
                a.obsp["spatial_connectivities"] = a.obsp["connectivities"].copy()

    monkeypatch.setattr(sd, "require", lambda *_a, **_k: _Sq)
    monkeypatch.setattr(
        sd.sc.tl,
        "leiden",
        lambda a, resolution, key_added: a.obs.__setitem__(key_added, ["0"] * a.n_obs),
    )
    _, _, stats = await sd._identify_domains_clustering(
        adata, SpatialDomainParameters(method="leiden", cluster_n_neighbors=8), _Ctx()
    )
    assert stats["n_neighbors"] == 8
    assert adata.uns["neighbors"]["params"]["n_neighbors"] == 8

    reused = _pca_adata()
    reused.obsm["spatial"] = _visium_pixels()[:80]
    sc.pp.neighbors(reused, n_neighbors=12, use_rep="X_pca")
    _, _, stats = await sd._identify_domains_clustering(
        reused, SpatialDomainParameters(method="leiden"), _Ctx()
    )
    assert stats["n_neighbors"] == 12  # reports the graph that was used


# ---------------------------------------------------------------- timeouts/schema


def test_training_backends_get_a_longer_default_timeout() -> None:
    assert sd._resolve_timeout(SpatialDomainParameters(method="stagate")) == 3600
    assert sd._resolve_timeout(SpatialDomainParameters(method="graphst")) == 3600
    assert sd._resolve_timeout(SpatialDomainParameters(method="spagcn")) == 600
    assert (
        sd._resolve_timeout(SpatialDomainParameters(method="graphst", timeout=90)) == 90
    )


@pytest.mark.parametrize(
    "field, value",
    [
        ("resolution", 0.0),
        ("resolution", 50.0),
        ("stagate_rad_cutoff", 0.0),
        ("graphst_radius", 0),
        ("cluster_n_neighbors", 1),
        ("graphst_random_seed", -1),
        ("graphst_n_clusters", 0),
    ],
)
def test_schema_rejects_out_of_range_values(field: str, value) -> None:
    with pytest.raises(ValidationError):
        SpatialDomainParameters(**{field: value})


def test_schema_states_effective_defaults_and_applicability() -> None:
    fields = SpatialDomainParameters.model_fields
    for name in (
        "cluster_n_neighbors",
        "cluster_spatial_weight",
        "stagate_rad_cutoff",
        "stagate_random_seed",
        "graphst_n_clusters",
        "timeout",
    ):
        assert "None" in (fields[name].description or ""), name
    assert "graphst_refinement=True" in fields["graphst_radius"].description
    assert "Ignored by leiden" in fields["n_domains"].description
    assert "Not used by graphst" in fields["resolution"].description
