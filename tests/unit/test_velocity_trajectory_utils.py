"""Unit tests for velocity/trajectory utility contracts."""

from __future__ import annotations

from types import ModuleType, SimpleNamespace
import warnings

import numpy as np
import pandas as pd
import pytest

from chatspatial.tools import trajectory as traj
from chatspatial.tools import velocity as vel
from chatspatial.utils.exceptions import (
    DataCompatibilityError,
    DataError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


def _patch_trajectory_dependency(
    monkeypatch: pytest.MonkeyPatch,
    dependency_name: str,
    dependency: object,
) -> None:
    def _require(requested_name: str, *_args, **_kwargs):
        if requested_name != dependency_name:
            raise AssertionError(f"Unexpected dependency request: {requested_name}")
        return dependency

    monkeypatch.setattr(traj, "require", _require)


def test_validate_velovi_data_contracts(minimal_spatial_adata):
    adata = SimpleNamespace(layers={})
    with pytest.raises(DataNotFoundError, match="Missing required layers"):
        vel._validate_velovi_data(adata)

    adata.layers["Ms"] = np.ones((10, 4))
    adata.layers["Mu"] = np.ones((10, 5))
    with pytest.raises(DataError, match="Shape mismatch"):
        vel._validate_velovi_data(adata)

    adata.layers["Mu"] = np.ones((10, 4))
    assert vel._validate_velovi_data(adata) is None


def test_compute_velocity_moments_uses_explicit_scanpy_neighbors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    calls: list[tuple[str, dict[str, int | None]]] = []

    def _neighbors(_adata, **kwargs):
        calls.append(("neighbors", kwargs))

    def _moments(adata_obj, **kwargs):
        calls.append(("moments", kwargs))
        adata_obj.layers["Ms"] = np.ones(adata_obj.shape)
        adata_obj.layers["Mu"] = np.ones(adata_obj.shape)

    fake_scanpy = ModuleType("scanpy")
    fake_scanpy.pp = SimpleNamespace(neighbors=_neighbors)
    fake_scvelo = ModuleType("scvelo")
    fake_scvelo.pp = SimpleNamespace(moments=_moments)
    monkeypatch.setitem(__import__("sys").modules, "scanpy", fake_scanpy)
    vel._compute_velocity_moments(
        adata,
        n_pcs=11,
        n_neighbors=7,
        scv=fake_scvelo,
    )

    assert calls == [
        ("neighbors", {"n_pcs": 11, "n_neighbors": 7}),
        ("moments", {"n_pcs": None, "n_neighbors": None}),
    ]


def test_preprocess_for_velocity_maps_validation_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    def _raise_validate(*_args, **_kwargs):
        raise DataNotFoundError("velocity layers missing")

    monkeypatch.setattr(vel, "validate_adata", _raise_validate)
    with pytest.raises(DataError, match="Invalid velocity data"):
        vel.preprocess_for_velocity(adata)


def test_compute_rna_velocity_dynamical_calls_expected_steps(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["Ms"] = np.ones((adata.n_obs, adata.n_vars))
    adata.layers["Mu"] = np.ones((adata.n_obs, adata.n_vars))
    called: list[str] = []

    fake_scv = ModuleType("scvelo")
    fake_scv.tl = SimpleNamespace(
        recover_dynamics=lambda *_args, **_kwargs: called.append("recover"),
        velocity=lambda *_args, **kwargs: called.append(f"velocity:{kwargs['mode']}"),
        latent_time=lambda *_args, **_kwargs: called.append("latent_time"),
        velocity_graph=lambda *_args, **_kwargs: called.append("graph"),
    )
    out = vel.compute_rna_velocity(adata, mode="dynamical", scv=fake_scv)
    assert out is adata
    assert called == ["recover", "velocity:dynamical", "latent_time", "graph"]


def test_compute_rna_velocity_preprocess_fallback(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    called: dict[str, bool] = {}

    fake_scv = ModuleType("scvelo")
    fake_scv.tl = SimpleNamespace(
        velocity=lambda *_args, **_kwargs: None,
        velocity_graph=lambda *_args, **_kwargs: None,
    )

    def _fake_preprocess(_adata, params=None, *, scv=None):
        del params
        called["preprocess"] = True
        called["shared_dependency"] = scv is fake_scv
        return _adata

    monkeypatch.setattr(vel, "preprocess_for_velocity", _fake_preprocess)

    vel.compute_rna_velocity(adata, mode="stochastic", scv=fake_scv)
    assert called["preprocess"] is True
    assert called["shared_dependency"] is True


def test_compute_rna_velocity_validates_mode_before_dependency_resolution(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(
        vel,
        "require",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("dependency resolution must follow mode validation")
        ),
    )

    with pytest.raises(ParameterError, match="Unsupported scVelo mode"):
        vel.compute_rna_velocity(minimal_spatial_adata.copy(), mode="invalid")


def test_velocity_metadata_contains_only_selected_method_parameters():
    scvelo_params = vel.RNAVelocityParameters(
        method="scvelo", scvelo_mode="dynamical"
    )
    scvelo_metadata = vel._velocity_parameters_for_metadata(scvelo_params)
    assert scvelo_metadata["mode"] == "dynamical"
    assert "n_epochs" not in scvelo_metadata

    velovi_params = vel.RNAVelocityParameters(
        method="velovi",
        velovi_n_epochs=12,
        velovi_n_hidden=64,
        velovi_n_latent=7,
        velovi_n_layers=2,
        velovi_dropout_rate=0.2,
        velovi_learning_rate=0.005,
        velovi_use_gpu=True,
    )
    velovi_metadata = vel._velocity_parameters_for_metadata(velovi_params)
    assert velovi_metadata == {
        "n_top_genes": 2000,
        "n_pcs": 30,
        "n_neighbors": 30,
        "min_shared_counts": 30,
        "n_epochs": 12,
        "n_hidden": 64,
        "n_latent": 7,
        "n_layers": 2,
        "dropout_rate": 0.2,
        "learning_rate": 0.005,
        "use_gpu": True,
    }


def test_prepare_gam_model_for_visualization_errors_and_success(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = SimpleNamespace(
        obs=pd.DataFrame({"latent_time": np.linspace(0, 1, 6)}),
        obsm={},
        var_names=["gene_0", "gene_1"],
    )

    monkeypatch.setattr(traj, "validate_obs_column", lambda *_args, **_kwargs: None)

    fake_cellrank_models = ModuleType("cellrank.models")
    fake_cellrank_models.GAM = lambda ad: ("GAM", len(ad.obs))
    monkeypatch.setattr(
        traj,
        "require_module",
        lambda *_args, **_kwargs: fake_cellrank_models,
    )

    with pytest.raises(DataNotFoundError, match="Fate probabilities"):
        traj.prepare_gam_model_for_visualization(adata, genes=["gene_0"])

    class _LineageNoNames:
        names = None

    adata.obsm["lineages_fwd"] = _LineageNoNames()
    with pytest.raises(DataError, match="must be a CellRank Lineage object"):
        traj.prepare_gam_model_for_visualization(adata, genes=["gene_0"])

    class _Lineage:
        names = ["state_a", "state_b"]

    adata.obsm["lineages_fwd"] = _Lineage()
    with pytest.raises(DataNotFoundError, match="Genes not found"):
        traj.prepare_gam_model_for_visualization(adata, genes=["MISSING_GENE"])

    model, lineages = traj.prepare_gam_model_for_visualization(adata, genes=["gene_0"])
    assert model[0] == "GAM"
    assert lineages == ["state_a", "state_b"]


def test_infer_pseudotime_palantir_root_validation(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(
        traj,
        "require",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("Root validation must precede dependency resolution")
        ),
    )

    with pytest.raises(ParameterError, match="Root cell 'missing' not found"):
        traj.infer_pseudotime_palantir(adata, root_cells=["missing"])


def test_compute_dpt_trajectory_root_validation_and_error_wrap(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(traj, "ensure_pca", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(traj, "ensure_neighbors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(traj, "ensure_diffmap", lambda *_args, **_kwargs: None)

    with pytest.raises(ParameterError, match="Root cell 'missing' not found"):
        traj.compute_dpt_trajectory(adata, root_cells=["missing"])

    fake_scanpy = ModuleType("scanpy")
    fake_scanpy.tl = SimpleNamespace(
        dpt=lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("fail"))
    )
    monkeypatch.setitem(__import__("sys").modules, "scanpy", fake_scanpy)

    with pytest.raises(ProcessingError, match="DPT computation failed"):
        traj.compute_dpt_trajectory(adata, root_cells=None)


class _VelCtx:
    def __init__(self, adata):
        self._adata = adata
        self.warnings: list[str] = []
        self.set_calls = 0

    async def get_adata(self, _data_id: str):
        return self._adata

    async def set_adata(self, _data_id: str, adata) -> None:
        self._adata = adata
        self.set_calls += 1

    async def warning(self, message: str) -> None:
        self.warnings.append(message)


@pytest.mark.asyncio
async def test_analyze_rna_velocity_missing_layers_raises_data_not_found(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(vel, "require", lambda *_a, **_k: None)
    monkeypatch.setattr(
        vel,
        "validate_adata",
        lambda *_a, **_k: (_ for _ in ()).throw(DataNotFoundError("spliced missing")),
    )
    monkeypatch.setitem(__import__("sys").modules, "scvelo", ModuleType("scvelo"))

    with pytest.raises(DataNotFoundError, match="Missing velocity data"):
        await vel.analyze_rna_velocity("d1", _VelCtx(adata))


@pytest.mark.asyncio
async def test_analyze_rna_velocity_scvelo_wraps_compute_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars))
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars))

    monkeypatch.setattr(vel, "require", lambda *_a, **_k: None)
    monkeypatch.setitem(__import__("sys").modules, "scvelo", ModuleType("scvelo"))
    monkeypatch.setattr(vel, "validate_adata", lambda *_a, **_k: None)

    def _partially_mutate_then_fail(candidate, *_args, **_kwargs):
        candidate.obs["partial_velocity"] = 1.0
        raise RuntimeError("cv boom")

    monkeypatch.setattr(vel, "compute_rna_velocity", _partially_mutate_then_fail)

    ctx = _VelCtx(adata)

    with pytest.raises(ProcessingError, match="scVelo RNA velocity analysis failed"):
        await vel.analyze_rna_velocity("d2", ctx)

    assert "partial_velocity" not in adata.obs
    assert ctx._adata is adata
    assert ctx.set_calls == 0


@pytest.mark.asyncio
async def test_analyze_rna_velocity_scvelo_preserves_domain_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars))
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars))
    fake_scv = ModuleType("scvelo")

    monkeypatch.setattr(vel, "require", lambda *_args, **_kwargs: fake_scv)
    monkeypatch.setattr(vel, "validate_adata", lambda *_args, **_kwargs: None)
    error = ParameterError("managed velocity error")
    monkeypatch.setattr(
        vel,
        "compute_rna_velocity",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(error),
    )

    with pytest.raises(ParameterError) as exc_info:
        await vel.analyze_rna_velocity("d2", _VelCtx(adata))
    assert exc_info.value is error


@pytest.mark.asyncio
async def test_analyze_rna_velocity_velovi_success_stores_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars))
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars))
    adata.uns["velocity_graph"] = "stale-scvelo-graph"

    captured = {}

    monkeypatch.setattr(vel, "require", lambda *_a, **_k: None)
    monkeypatch.setitem(__import__("sys").modules, "scvelo", ModuleType("scvelo"))
    monkeypatch.setattr(vel, "validate_adata", lambda *_a, **_k: None)

    async def _fake_velovi(*_a, **_k):
        return {
            "velocity_computed": True,
            "velocity_shape": (adata.n_obs, adata.n_vars),
            "n_genes_analyzed": 12,
            "original_n_genes": adata.n_vars,
            "velocity_mean_norm": 0.4,
            "velocity_std_norm": 0.1,
            "device": "CPU",
        }

    monkeypatch.setattr(vel, "analyze_velocity_with_velovi", _fake_velovi)
    monkeypatch.setattr(
        vel,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(vel, "export_analysis_result", lambda *_a, **_k: [])

    out = await vel.analyze_rna_velocity(
        "d3",
        _VelCtx(adata),
        vel.RNAVelocityParameters(method="velovi"),
    )

    assert out.velocity_computed is True
    assert out.mode == "velovi"
    assert out.velocity_graph_key is None
    assert captured["analysis_name"] == "velocity_velovi"
    assert captured["statistics"]["velocity_computed"] is True
    assert captured["statistics"]["n_genes_analyzed"] == 12
    assert captured["statistics"]["device"] == "CPU"


@pytest.mark.asyncio
async def test_analyze_trajectory_cellrank_requires_velocity_data(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()

    with pytest.raises(ProcessingError, match="CellRank requires velocity data"):
        await traj.analyze_trajectory(
            "t1",
            _VelCtx(adata),
            traj.TrajectoryParameters(method="cellrank"),
        )


@pytest.mark.asyncio
async def test_analyze_trajectory_palantir_success_records_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    captured = {}

    def _fake_palantir(_adata, **_kwargs):
        _adata.obs["palantir_pseudotime"] = np.linspace(0, 1, _adata.n_obs)
        _adata.obsm["palantir_branch_probs"] = np.ones((_adata.n_obs, 2), dtype=float)
        return _adata

    monkeypatch.setattr(traj, "infer_pseudotime_palantir", _fake_palantir)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_a, **_k: [],
    )

    out = await traj.analyze_trajectory(
        "t2",
        _VelCtx(adata),
        traj.TrajectoryParameters(method="palantir"),
    )

    assert out.pseudotime_computed is True
    assert out.method == "palantir"
    assert out.pseudotime_key == "palantir_pseudotime"
    assert captured["analysis_name"] == "trajectory_palantir"
    assert captured["results_keys"]["obsm"] == ["palantir_branch_probs"]


@pytest.mark.asyncio
@pytest.mark.parametrize("method", ["palantir", "dpt"])
async def test_analyze_trajectory_preserves_root_cell_parameter_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, method: str
):
    adata = minimal_spatial_adata.copy()
    helper_name = {
        "palantir": "infer_pseudotime_palantir",
        "dpt": "compute_dpt_trajectory",
    }[method]

    monkeypatch.setattr(
        traj,
        helper_name,
        lambda *_a, **_k: (_ for _ in ()).throw(
            ParameterError("Root cell 'missing' not found")
        ),
    )

    with pytest.raises(ParameterError, match="Root cell 'missing' not found"):
        await traj.analyze_trajectory(
            "t3",
            _VelCtx(adata),
            traj.TrajectoryParameters(method=method, root_cells=["missing"]),
        )


@pytest.mark.asyncio
async def test_analyze_trajectory_dpt_wraps_internal_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    def _partially_mutate_then_fail(candidate, *_args, **_kwargs):
        candidate.obs["partial_pseudotime"] = 0.5
        raise RuntimeError("dpt boom")

    monkeypatch.setattr(traj, "compute_dpt_trajectory", _partially_mutate_then_fail)

    ctx = _VelCtx(adata)

    with pytest.raises(ProcessingError, match="DPT analysis failed"):
        await traj.analyze_trajectory(
            "t3",
            ctx,
            traj.TrajectoryParameters(method="dpt"),
        )

    assert "partial_pseudotime" not in adata.obs
    assert ctx._adata is adata
    assert ctx.set_calls == 0


@pytest.mark.asyncio
async def test_analyze_trajectory_unknown_method_raises_parameter_error(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    params = traj.TrajectoryParameters(method="dpt").model_copy(
        update={"method": "unknown"}
    )

    with pytest.raises(ParameterError, match="Unknown trajectory method"):
        await traj.analyze_trajectory("t4", _VelCtx(adata), params)


def _install_fake_cellrank(
    monkeypatch: pytest.MonkeyPatch,
    *,
    n_obs: int,
    terminal_categories: list[str] | None = None,
    fail_macrostates: bool = False,
):
    class _Kernel:
        def __init__(self, *_a, **_k):
            self.computed = False

        def compute_transition_matrix(self):
            self.computed = True
            return self

        def __mul__(self, _other):
            return self

        def __rmul__(self, _other):
            return self

        def __add__(self, _other):
            return self

    class _Memberships:
        def __init__(self, n: int):
            self._x = np.linspace(0.1, 0.9, n, dtype=float).reshape(n, 1)

        def __getitem__(self, idx):
            return SimpleNamespace(X=self._x[idx])

    class _FateProbabilities(np.ndarray):
        def __new__(cls, n: int):
            arr = np.linspace(0.2, 0.8, n, dtype=float).reshape(n, 1)
            return arr.view(cls)

        def __getitem__(self, key):
            if isinstance(key, str):
                return SimpleNamespace(X=np.asarray(self))
            return super().__getitem__(key)

    class _GPCCA:
        def __init__(self, _kernel):
            self.terminal_states = None
            self.fate_probabilities = None
            self.macrostates = pd.Series(
                pd.Categorical(["M0"] * n_obs, categories=["M0"])
            )
            self.macrostates_memberships = _Memberships(n_obs)

        def compute_eigendecomposition(self):
            return None

        def compute_macrostates(self, n_states: int):
            if fail_macrostates:
                raise RuntimeError("macro failed")
            self.n_states = n_states

        def predict_terminal_states(self, method: str = "stability"):
            _ = method
            cats = terminal_categories or []
            if cats:
                values = [cats[0]] * n_obs
                self.terminal_states = pd.Series(
                    pd.Categorical(values, categories=cats)
                )
            else:
                self.terminal_states = pd.Series(pd.Categorical([], categories=[]))

        def compute_fate_probabilities(self):
            self.fate_probabilities = _FateProbabilities(n_obs)

    fake_cr = ModuleType("cellrank")
    fake_cr.kernels = SimpleNamespace(
        VelocityKernel=_Kernel,
        ConnectivityKernel=_Kernel,
        PrecomputedKernel=_Kernel,
    )
    fake_cr.estimators = SimpleNamespace(GPCCA=_GPCCA)
    _patch_trajectory_dependency(monkeypatch, "cellrank", fake_cr)


def test_infer_spatial_trajectory_cellrank_velovi_path_transfers_results(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.uns["velocity_method"] = "velovi"
    reconstructed = SimpleNamespace(
        obs=pd.DataFrame(index=adata.obs_names.copy()),
        obsm={},
        uns={},
    )
    cleanup_called = {"v": False}

    monkeypatch.setattr(
        traj,
        "ensure_cellrank_compat",
        lambda: lambda: cleanup_called.__setitem__("v", True),
    )
    monkeypatch.setattr(traj, "get_spatial_key", lambda _a: "spatial")
    monkeypatch.setattr(traj, "has_velovi_essential_data", lambda _a: True)
    monkeypatch.setattr(traj, "reconstruct_velovi_adata", lambda _a: reconstructed)

    _install_fake_cellrank(
        monkeypatch,
        n_obs=adata.n_obs,
        terminal_categories=["T0"],
        fail_macrostates=False,
    )

    out = traj.infer_spatial_trajectory_cellrank(adata, spatial_weight=0.5, n_states=3)
    assert out is adata
    assert "pseudotime" in adata.obs
    assert "terminal_states" in adata.obs
    assert "macrostates" in adata.obs
    assert "fate_probabilities" in adata.obsm
    # Analysis must also write CellRank-standard alias
    assert "to_terminal_states" in adata.obsm
    np.testing.assert_array_equal(
        adata.obsm["to_terminal_states"],
        adata.obsm["fate_probabilities"],
    )
    assert cleanup_called["v"] is True


def test_infer_spatial_trajectory_cellrank_wraps_macrostate_errors_and_cleans_up(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    cleanup_called = {"v": False}

    monkeypatch.setattr(
        traj,
        "ensure_cellrank_compat",
        lambda: lambda: cleanup_called.__setitem__("v", True),
    )
    monkeypatch.setattr(traj, "get_spatial_key", lambda _a: None)
    _install_fake_cellrank(
        monkeypatch,
        n_obs=adata.n_obs,
        terminal_categories=["T0"],
        fail_macrostates=True,
    )

    with pytest.raises(ProcessingError, match="CellRank macrostate computation failed"):
        traj.infer_spatial_trajectory_cellrank(adata, spatial_weight=0.0, n_states=9)

    assert cleanup_called["v"] is True


def test_infer_spatial_trajectory_cellrank_falls_back_to_macrostate_pseudotime(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    monkeypatch.setattr(traj, "ensure_cellrank_compat", lambda: (lambda: None))
    monkeypatch.setattr(traj, "get_spatial_key", lambda _a: None)
    _install_fake_cellrank(
        monkeypatch,
        n_obs=adata.n_obs,
        terminal_categories=[],
        fail_macrostates=False,
    )

    out = traj.infer_spatial_trajectory_cellrank(adata, spatial_weight=0.0, n_states=4)
    assert out is adata
    assert "pseudotime" in adata.obs
    assert "macrostates" in adata.obs


@pytest.mark.asyncio
async def test_analyze_velocity_with_velovi_success_handles_zero_latent_time(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata_prepared = adata.copy()
    adata_prepared.layers["Ms"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata_prepared.layers["Mu"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    called = {"setup": False, "train": False, "stored": False}
    canonical_time = np.arange(
        adata.n_obs * adata.n_vars, dtype=float
    ).reshape(adata.n_obs, adata.n_vars)
    canonical_time[:, 0] = 0.0
    canonical_velocity = np.arange(
        1, adata.n_obs * adata.n_vars + 1, dtype=float
    ).reshape(adata.n_obs, adata.n_vars)

    class _FakeVELOVI:
        @staticmethod
        def setup_anndata(_adata, spliced_layer: str, unspliced_layer: str):
            assert spliced_layer == "Ms"
            assert unspliced_layer == "Mu"
            called["setup"] = True

        def __init__(
            self,
            _adata,
            n_hidden: int,
            n_latent: int,
            n_layers: int = 1,
            dropout_rate: float = 0.1,
        ):
            assert n_hidden == 32
            assert n_latent == 5

        def train(self, max_epochs: int, accelerator: str, lr: float = 1e-3):
            assert max_epochs == 2
            assert accelerator == "cpu"
            called["train"] = True

        def get_latent_time(self, n_samples: int, return_numpy: bool):
            assert n_samples == 25
            assert return_numpy is False
            return pd.DataFrame(
                canonical_time,
                index=adata_prepared.obs_names,
                columns=adata_prepared.var_names,
            ).iloc[::-1, ::-1]

        def get_velocity(
            self, n_samples: int, velo_statistic: str, return_numpy: bool
        ):
            assert n_samples == 25
            assert velo_statistic == "mean"
            assert return_numpy is False
            return pd.DataFrame(
                canonical_velocity,
                index=adata_prepared.obs_names,
                columns=adata_prepared.var_names,
            ).iloc[::-1, ::-1]

        def get_latent_representation(self):
            return np.ones((adata.n_obs, 3), dtype=float)

    fake_scvi_external = ModuleType("scvi.external")
    fake_scvi_external.VELOVI = _FakeVELOVI
    monkeypatch.setitem(__import__("sys").modules, "scvi.external", fake_scvi_external)

    monkeypatch.setattr(vel, "require", lambda *_a, **_k: None)
    monkeypatch.setattr(
        vel,
        "require_module",
        lambda *_args, **_kwargs: fake_scvi_external,
    )

    async def _fake_prepare(*_a, **_k):
        return adata_prepared

    monkeypatch.setattr(vel, "_prepare_velovi_data", _fake_prepare)
    monkeypatch.setattr(vel, "get_accelerator", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(
        vel,
        "store_velovi_essential_data",
        lambda *_a, **_k: called.__setitem__("stored", True),
    )

    ctx = _VelCtx(adata)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = await vel.analyze_velocity_with_velovi(
            adata, n_epochs=2, n_hidden=32, n_latent=5, use_gpu=True, ctx=ctx
        )
    assert result["velocity_computed"] is True
    assert result["velocity_shape"] == (adata.n_obs, adata.n_vars)
    assert called == {"setup": True, "train": True, "stored": True}
    assert "velocity_velovi_norm" in adata.obs
    assert "X_velovi_latent" in adata.obsm
    assert "latent_time_velovi" not in adata.obs
    pd.testing.assert_frame_equal(
        adata.obsm["velovi_latent_time"],
        pd.DataFrame(
            canonical_time.astype(np.float32),
            index=adata.obs_names,
            columns=adata.var_names,
        ),
    )
    scaling = np.ones(adata.n_vars, dtype=float)
    positive = canonical_time.max(axis=0) > 0
    scaling[positive] = 20.0 / canonical_time.max(axis=0)[positive]
    np.testing.assert_allclose(
        adata_prepared.layers["velocity_velovi"],
        canonical_velocity / scaling[np.newaxis, :],
    )
    assert result["device"] == "CPU"
    assert ctx.warnings == ["GPU requested but unavailable; VELOVI is using CPU."]
    assert not any("divide by zero" in str(w.message) for w in caught)


@pytest.mark.asyncio
async def test_analyze_velocity_with_velovi_wraps_model_failures(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata_prepared = adata.copy()
    adata_prepared.layers["Ms"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata_prepared.layers["Mu"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)

    class _FailVELOVI:
        @staticmethod
        def setup_anndata(*_a, **_k):
            return None

        def __init__(self, *_a, **_k):
            return None

        def train(self, *_a, **_k):
            raise RuntimeError("train failed")

    fake_scvi_external = ModuleType("scvi.external")
    fake_scvi_external.VELOVI = _FailVELOVI
    monkeypatch.setitem(__import__("sys").modules, "scvi.external", fake_scvi_external)

    monkeypatch.setattr(vel, "require", lambda *_a, **_k: None)
    monkeypatch.setattr(
        vel,
        "require_module",
        lambda *_args, **_kwargs: fake_scvi_external,
    )

    async def _fake_prepare(*_a, **_k):
        return adata_prepared

    monkeypatch.setattr(vel, "_prepare_velovi_data", _fake_prepare)
    monkeypatch.setattr(vel, "get_accelerator", lambda prefer_gpu=False: "cpu")

    with pytest.raises(ProcessingError, match="VELOVI velocity analysis failed"):
        await vel.analyze_velocity_with_velovi(adata, n_epochs=1, ctx=None)


def test_preprocess_for_velocity_uses_params_object_values(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)

    calls: dict[str, dict[str, int | bool]] = {}
    fake_scv = ModuleType("scvelo")
    fake_scv.pp = SimpleNamespace(
        filter_and_normalize=lambda _adata, **kwargs: calls.__setitem__(
            "filter", kwargs
        ),
    )
    monkeypatch.setattr(vel, "validate_adata", lambda *_a, **_k: None)
    monkeypatch.setattr(
        vel,
        "_compute_velocity_moments",
        lambda _adata, **kwargs: calls.__setitem__("moments", kwargs),
    )

    params = vel.RNAVelocityParameters(
        min_shared_counts=11,
        n_top_genes=123,
        n_pcs=17,
        n_neighbors=9,
    )
    out = vel.preprocess_for_velocity(adata, params=params, scv=fake_scv)
    assert out is adata
    assert calls["filter"] == {
        "min_shared_counts": 11,
        "n_top_genes": 123,
        "enforce": True,
    }
    assert calls["moments"] == {
        "n_pcs": 17,
        "n_neighbors": 9,
        "scv": fake_scv,
    }


def test_compute_rna_velocity_uses_params_mode_override(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["Ms"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["Mu"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    seen: dict[str, str] = {}

    fake_scv = ModuleType("scvelo")
    fake_scv.tl = SimpleNamespace(
        velocity=lambda *_a, **kwargs: seen.__setitem__("mode", kwargs["mode"]),
        velocity_graph=lambda *_a, **_k: None,
    )
    params = vel.RNAVelocityParameters(scvelo_mode="deterministic")
    vel.compute_rna_velocity(
        adata,
        mode="stochastic",
        params=params,
        scv=fake_scv,
    )
    assert seen["mode"] == "deterministic"


@pytest.mark.asyncio
async def test_prepare_velovi_data_raises_on_scv_preprocessing_failure(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)

    def _raise_pp(*_a, **_k):
        raise RuntimeError("pp fail")

    fake_scv = ModuleType("scvelo")
    fake_scv.pp = SimpleNamespace(
        filter_and_normalize=_raise_pp,
    )
    with pytest.raises(
        ProcessingError, match="VELOVI preprocessing failed"
    ):
        await vel._prepare_velovi_data(adata, None, scv=fake_scv)


@pytest.mark.asyncio
async def test_prepare_velovi_data_raises_on_moments_failure(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)

    def _raise_moments(*_a, **_k):
        raise RuntimeError("moments fail")

    fake_scv = ModuleType("scvelo")
    fake_scv.pp = SimpleNamespace(filter_and_normalize=lambda *_a, **_k: None)
    monkeypatch.setattr(vel, "_compute_velocity_moments", _raise_moments)

    with pytest.raises(
        ProcessingError, match="Moments computation failed"
    ):
        await vel._prepare_velovi_data(adata, None, scv=fake_scv)


@pytest.mark.asyncio
async def test_prepare_velovi_data_requires_spliced_and_unspliced_layers(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers.clear()
    monkeypatch.setattr(
        vel,
        "require",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("layer validation must precede dependency resolution")
        ),
    )

    with pytest.raises(
        DataNotFoundError, match="Missing required 'spliced' and 'unspliced' layers"
    ):
        await vel._prepare_velovi_data(adata, ctx=None)


@pytest.mark.asyncio
async def test_prepare_velovi_data_keeps_input_unmodified_and_builds_independent_workspace(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=np.float32)
    adata.layers["unspliced"] = (
        np.ones((adata.n_obs, adata.n_vars), dtype=np.float32) * 2
    )
    spliced_before = adata.layers["spliced"].copy()
    unspliced_before = adata.layers["unspliced"].copy()

    def _fake_filter_and_normalize(adata_obj, **_kwargs):
        adata_obj.X[0, 0] = 999.0

    def _fake_moments(adata_obj, **_kwargs):
        adata_obj.obsp["connectivities"] = np.eye(adata_obj.n_obs, dtype=np.float32)

    fake_scv = ModuleType("scvelo")
    fake_scv.pp = SimpleNamespace(
        filter_and_normalize=_fake_filter_and_normalize,
    )
    monkeypatch.setattr(vel, "_compute_velocity_moments", _fake_moments)

    out = await vel._prepare_velovi_data(adata, ctx=None, scv=fake_scv)

    # Working object is independent.
    assert out is not adata
    assert out.layers["spliced"] is not adata.layers["spliced"]
    assert out.layers["unspliced"] is not adata.layers["unspliced"]
    assert "Ms" in out.layers and "Mu" in out.layers

    # Input adata remains unchanged by working-copy preprocessing.
    np.testing.assert_allclose(adata.layers["spliced"], spliced_before)
    np.testing.assert_allclose(adata.layers["unspliced"], unspliced_before)
    assert "Ms" not in adata.layers
    assert "Mu" not in adata.layers


def test_validate_velovi_data_rejects_non_2d_layers():
    adata = SimpleNamespace(layers={"Ms": np.ones(5), "Mu": np.ones(5)})
    with pytest.raises(DataError, match="Expected 2D arrays"):
        vel._validate_velovi_data(adata)


@pytest.mark.parametrize(
    ("layer_name", "values", "message"),
    [
        ("Ms", np.array([[np.nan, 1.0]]), "non-finite"),
        ("Mu", np.array([[-1.0, 1.0]]), "negative"),
        ("Ms", np.zeros((1, 2)), "no positive values"),
    ],
)
def test_validate_velovi_data_rejects_invalid_layer_values(
    layer_name: str, values: np.ndarray, message: str
):
    layers = {
        "Ms": np.ones((1, 2), dtype=float),
        "Mu": np.ones((1, 2), dtype=float),
    }
    layers[layer_name] = values
    with pytest.raises(DataError, match=message):
        vel._validate_velovi_data(SimpleNamespace(layers=layers))


def test_normalize_velovi_cell_gene_output_aligns_labeled_axes(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata[:, :3].copy()
    expected = np.arange(adata.n_obs * adata.n_vars, dtype=float).reshape(
        adata.n_obs, adata.n_vars
    )
    result = pd.DataFrame(
        expected,
        index=adata.obs_names,
        columns=adata.var_names,
    ).iloc[::-1, ::-1]

    normalized = vel._normalize_velovi_cell_gene_output(
        result, adata, "velocity"
    )

    np.testing.assert_array_equal(normalized, expected.astype(np.float32))


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (
            lambda frame: frame.rename(index={frame.index[0]: "unexpected"}),
            "axes do not match",
        ),
        (
            lambda frame: frame.set_axis(
                [frame.index[0], *frame.index[1:-1], frame.index[0]], axis=0
            ),
            "duplicate observation IDs",
        ),
        (
            lambda frame: frame.set_axis(
                [frame.columns[0], *frame.columns[1:-1], frame.columns[0]], axis=1
            ),
            "duplicate gene IDs",
        ),
    ],
)
def test_normalize_velovi_cell_gene_output_rejects_invalid_labels(
    minimal_spatial_adata, mutate, message: str
):
    adata = minimal_spatial_adata[:, :3].copy()
    result = pd.DataFrame(
        np.ones((adata.n_obs, adata.n_vars)),
        index=adata.obs_names,
        columns=adata.var_names,
    )

    with pytest.raises(DataCompatibilityError, match=message):
        vel._normalize_velovi_cell_gene_output(
            mutate(result), adata, "velocity"
        )


def test_normalize_velovi_cell_gene_output_rejects_wrong_array_shape(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata[:, :1].copy()
    result = np.ones((adata.n_obs - 1, adata.n_vars))
    with pytest.raises(DataCompatibilityError, match="has shape"):
        vel._normalize_velovi_cell_gene_output(
            result,
            adata,
            "latent time",
        )


@pytest.mark.parametrize(
    ("value", "message"),
    [(np.nan, "non-finite"), (-1.0, "negative")],
)
def test_normalize_velovi_latent_time_rejects_invalid_values(
    minimal_spatial_adata, value: float, message: str
):
    adata = minimal_spatial_adata[:, :1].copy()
    result = np.full((adata.n_obs, adata.n_vars), value)
    with pytest.raises(DataCompatibilityError, match=message):
        vel._normalize_velovi_cell_gene_output(
            result,
            adata,
            "latent time",
            require_nonnegative=True,
        )


@pytest.mark.parametrize(
    "case",
    ["one_dimensional", "empty_latent_axis", "nonfinite"],
)
def test_normalize_velovi_latent_representation_rejects_invalid_output(
    minimal_spatial_adata, case: str
):
    if case == "one_dimensional":
        result = np.ones(minimal_spatial_adata.n_obs)
    elif case == "empty_latent_axis":
        result = np.ones((minimal_spatial_adata.n_obs, 0))
    else:
        result = np.full((minimal_spatial_adata.n_obs, 2), np.nan)
    with pytest.raises(DataCompatibilityError, match="latent representation"):
        vel._normalize_velovi_latent_representation(
            result, minimal_spatial_adata
        )


def test_validate_velovi_workspace_axis_rejects_reordered_cells(
    minimal_spatial_adata,
):
    prepared = minimal_spatial_adata[::-1].copy()
    with pytest.raises(DataCompatibilityError, match="identities or order"):
        vel._validate_velovi_workspace_axis(minimal_spatial_adata, prepared)


@pytest.mark.asyncio
async def test_analyze_rna_velocity_scvelo_success_includes_latent_time_result_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars))
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars))
    adata.obs["velocity_velovi_norm"] = 99.0
    adata.obsm["X_velovi_latent"] = np.full((adata.n_obs, 2), 99.0)
    adata.obsm["velovi_latent_time"] = pd.DataFrame(
        np.full((adata.n_obs, 2), 99.0),
        index=adata.obs_names,
        columns=["stale_g1", "stale_g2"],
    )
    captured: dict[str, object] = {}

    fake_scv = ModuleType("scvelo")
    monkeypatch.setattr(vel, "require", lambda *_a, **_k: fake_scv)
    monkeypatch.setattr(vel, "validate_adata", lambda *_a, **_k: None)

    def _fake_compute(_adata, mode: str, params, *, scv):
        del mode, params
        assert scv is fake_scv
        _adata.obs["latent_time"] = np.linspace(0, 1, _adata.n_obs)
        return _adata

    monkeypatch.setattr(vel, "compute_rna_velocity", _fake_compute)
    monkeypatch.setattr(
        vel, "store_analysis_metadata", lambda _adata, **kwargs: captured.update(kwargs)
    )
    monkeypatch.setattr(vel, "export_analysis_result", lambda *_a, **_k: [])

    out = await vel.analyze_rna_velocity(
        "vel-s1",
        _VelCtx(adata),
        vel.RNAVelocityParameters(method="scvelo", scvelo_mode="dynamical"),
    )
    assert out.velocity_computed is True
    assert out.mode == "dynamical"
    assert "latent_time" in captured["results_keys"]["obs"]
    assert "velocity_velovi_norm" not in captured["results_keys"]["obs"]
    assert captured["results_keys"]["obsm"] == []
    assert captured["parameters"]["mode"] == "dynamical"


@pytest.mark.asyncio
async def test_analyze_rna_velocity_velovi_records_velovi_result_keys(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars))
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars))
    captured: dict[str, object] = {}

    monkeypatch.setattr(vel, "require", lambda *_a, **_k: None)
    monkeypatch.setitem(__import__("sys").modules, "scvelo", ModuleType("scvelo"))
    monkeypatch.setattr(vel, "validate_adata", lambda *_a, **_k: None)

    async def _fake_velovi(_adata, **_kwargs):
        _adata.obs["velocity_velovi_norm"] = np.ones(_adata.n_obs, dtype=float)
        _adata.obsm["X_velovi_latent"] = np.ones((_adata.n_obs, 2), dtype=float)
        _adata.obsm["velovi_latent_time"] = pd.DataFrame(
            np.ones((_adata.n_obs, 2), dtype=float),
            index=_adata.obs_names,
            columns=["g1", "g2"],
        )
        return {"velocity_computed": True}

    monkeypatch.setattr(vel, "analyze_velocity_with_velovi", _fake_velovi)
    monkeypatch.setattr(
        vel, "store_analysis_metadata", lambda _adata, **kwargs: captured.update(kwargs)
    )
    monkeypatch.setattr(vel, "export_analysis_result", lambda *_a, **_k: [])

    out = await vel.analyze_rna_velocity(
        "vel-s2",
        _VelCtx(adata),
        vel.RNAVelocityParameters(method="velovi"),
    )
    assert out.velocity_computed is True
    assert "velocity_velovi_norm" in captured["results_keys"]["obs"]
    assert "X_velovi_latent" in captured["results_keys"]["obsm"]
    assert "velovi_latent_time" in captured["results_keys"]["obsm"]


@pytest.mark.asyncio
async def test_analyze_rna_velocity_velovi_false_result_raises_processing(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars))
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars))

    monkeypatch.setattr(vel, "require", lambda *_a, **_k: None)
    monkeypatch.setitem(__import__("sys").modules, "scvelo", ModuleType("scvelo"))
    monkeypatch.setattr(vel, "validate_adata", lambda *_a, **_k: None)

    async def _fake_velovi_false(*_a, **_k):
        return {"velocity_computed": False}

    monkeypatch.setattr(vel, "analyze_velocity_with_velovi", _fake_velovi_false)

    with pytest.raises(ProcessingError, match="VELOVI failed to compute velocity"):
        await vel.analyze_rna_velocity(
            "vel-s3",
            _VelCtx(adata),
            vel.RNAVelocityParameters(method="velovi"),
        )


@pytest.mark.asyncio
async def test_analyze_rna_velocity_unknown_method_raises_parameter_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars))
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars))

    monkeypatch.setattr(
        vel,
        "require",
        lambda *_a, **_k: (_ for _ in ()).throw(
            AssertionError("dependency checks must follow method validation")
        ),
    )
    monkeypatch.setitem(__import__("sys").modules, "scvelo", ModuleType("scvelo"))
    monkeypatch.setattr(vel, "validate_adata", lambda *_a, **_k: None)

    params = vel.RNAVelocityParameters(method="scvelo").model_copy(
        update={"method": "unknown"}
    )
    with pytest.raises(ParameterError, match="Unknown velocity method"):
        await vel.analyze_rna_velocity("vel-s4", _VelCtx(adata), params)


def test_infer_spatial_trajectory_cellrank_disables_spatial_when_coords_missing(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    cleanup_called = {"v": False}

    monkeypatch.setattr(
        traj,
        "ensure_cellrank_compat",
        lambda: lambda: cleanup_called.__setitem__("v", True),
    )
    monkeypatch.setattr(traj, "get_spatial_key", lambda _a: None)
    _install_fake_cellrank(
        monkeypatch,
        n_obs=adata.n_obs,
        terminal_categories=["T0"],
        fail_macrostates=False,
    )

    out = traj.infer_spatial_trajectory_cellrank(adata, spatial_weight=0.9, n_states=3)
    assert out is adata
    assert "pseudotime" in adata.obs
    assert cleanup_called["v"] is True


def test_infer_spatial_trajectory_cellrank_velovi_missing_essential_data_raises(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.uns["velocity_method"] = "velovi"

    monkeypatch.setattr(traj, "ensure_cellrank_compat", lambda: (lambda: None))
    monkeypatch.setattr(traj, "get_spatial_key", lambda _a: "spatial")
    monkeypatch.setattr(traj, "has_velovi_essential_data", lambda _a: False)
    _install_fake_cellrank(
        monkeypatch,
        n_obs=adata.n_obs,
        terminal_categories=["T0"],
        fail_macrostates=False,
    )

    with pytest.raises(ProcessingError, match="VELOVI velocity data not found"):
        traj.infer_spatial_trajectory_cellrank(adata)


def test_infer_spatial_trajectory_cellrank_predict_terminal_value_error_bubbles(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    cleanup_called = {"v": False}

    monkeypatch.setattr(
        traj,
        "ensure_cellrank_compat",
        lambda: lambda: cleanup_called.__setitem__("v", True),
    )
    monkeypatch.setattr(traj, "get_spatial_key", lambda _a: None)

    class _Kernel:
        def __init__(self, *_a, **_k):
            pass

        def compute_transition_matrix(self):
            return self

        def __mul__(self, _other):
            return self

        def __rmul__(self, _other):
            return self

        def __add__(self, _other):
            return self

    class _GPCCA:
        def __init__(self, _kernel):
            self.macrostates = pd.Series(pd.Categorical(["M0"] * adata.n_obs))

        def compute_eigendecomposition(self):
            return None

        def compute_macrostates(self, n_states: int):
            self.n_states = n_states

        def predict_terminal_states(self, method: str = "stability"):
            del method
            raise ValueError("unexpected state error")

    fake_cr = ModuleType("cellrank")
    fake_cr.kernels = SimpleNamespace(
        VelocityKernel=_Kernel,
        ConnectivityKernel=_Kernel,
        PrecomputedKernel=_Kernel,
    )
    fake_cr.estimators = SimpleNamespace(GPCCA=_GPCCA)
    _patch_trajectory_dependency(monkeypatch, "cellrank", fake_cr)

    with pytest.raises(ValueError, match="unexpected state error"):
        traj.infer_spatial_trajectory_cellrank(adata, spatial_weight=0.0)
    assert cleanup_called["v"] is True


def test_infer_spatial_trajectory_cellrank_wraps_fate_probability_failures(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(traj, "ensure_cellrank_compat", lambda: (lambda: None))
    monkeypatch.setattr(traj, "get_spatial_key", lambda _a: None)

    class _Kernel:
        def __init__(self, *_a, **_k):
            pass

        def compute_transition_matrix(self):
            return self

        def __mul__(self, _other):
            return self

        def __rmul__(self, _other):
            return self

        def __add__(self, _other):
            return self

    class _TerminalSeries:
        cat = SimpleNamespace(categories=["T0"])

    class _GPCCA:
        def __init__(self, _kernel):
            self.terminal_states = _TerminalSeries()
            self.macrostates = pd.Series(pd.Categorical(["M0"] * adata.n_obs))

        def compute_eigendecomposition(self):
            return None

        def compute_macrostates(self, n_states: int):
            self.n_states = n_states

        def predict_terminal_states(self, method: str = "stability"):
            del method
            return None

        def compute_fate_probabilities(self):
            raise RuntimeError("fate failed")

    fake_cr = ModuleType("cellrank")
    fake_cr.kernels = SimpleNamespace(
        VelocityKernel=_Kernel,
        ConnectivityKernel=_Kernel,
        PrecomputedKernel=_Kernel,
    )
    fake_cr.estimators = SimpleNamespace(GPCCA=_GPCCA)
    _patch_trajectory_dependency(monkeypatch, "cellrank", fake_cr)

    with pytest.raises(ProcessingError, match="fate probability computation failed"):
        traj.infer_spatial_trajectory_cellrank(adata, spatial_weight=0.0)


def test_infer_spatial_trajectory_cellrank_raises_when_no_terminal_or_macrostates(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(traj, "ensure_cellrank_compat", lambda: (lambda: None))
    monkeypatch.setattr(traj, "get_spatial_key", lambda _a: None)

    class _Kernel:
        def __init__(self, *_a, **_k):
            pass

        def compute_transition_matrix(self):
            return self

        def __mul__(self, _other):
            return self

        def __rmul__(self, _other):
            return self

        def __add__(self, _other):
            return self

    class _EmptyTerminalSeries:
        cat = SimpleNamespace(categories=[])

    class _GPCCA:
        def __init__(self, _kernel):
            self.terminal_states = _EmptyTerminalSeries()
            self.macrostates = None

        def compute_eigendecomposition(self):
            return None

        def compute_macrostates(self, n_states: int):
            self.n_states = n_states

        def predict_terminal_states(self, method: str = "stability"):
            del method
            return None

    fake_cr = ModuleType("cellrank")
    fake_cr.kernels = SimpleNamespace(
        VelocityKernel=_Kernel,
        ConnectivityKernel=_Kernel,
        PrecomputedKernel=_Kernel,
    )
    fake_cr.estimators = SimpleNamespace(GPCCA=_GPCCA)
    _patch_trajectory_dependency(monkeypatch, "cellrank", fake_cr)

    with pytest.raises(ProcessingError, match="could not compute terminal states"):
        traj.infer_spatial_trajectory_cellrank(adata, spatial_weight=0.0)


def test_infer_pseudotime_palantir_valid_root_populates_outputs(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_pca"] = np.ones((adata.n_obs, 5))
    monkeypatch.setattr(traj, "ensure_pca", lambda *_a, **_k: None)
    captured: dict[str, object] = {}

    class _PR:
        pseudotime = pd.Series(np.linspace(0, 1, adata.n_obs), index=adata.obs_names)
        branch_probs = np.ones((adata.n_obs, 2))

    fake_palantir = ModuleType("palantir")
    fake_palantir.utils = SimpleNamespace(
        run_diffusion_maps=lambda *_a, **_k: {"EigenVectors": adata.obsm["X_pca"]},
        determine_multiscale_space=lambda result: pd.DataFrame(
            result["EigenVectors"], index=adata.obs_names
        ),
    )

    def _fake_run_palantir(_ms_data, start_cell, num_waypoints):
        captured["start_cell"] = start_cell
        captured["num_waypoints"] = num_waypoints
        return _PR()

    fake_palantir.core = SimpleNamespace(run_palantir=_fake_run_palantir)

    out = traj.infer_pseudotime_palantir(
        adata,
        root_cells=[adata.obs_names[1]],
        num_waypoints=123,
        palantir_module=fake_palantir,
    )
    assert out is adata
    assert captured["start_cell"] == adata.obs_names[1]
    assert captured["num_waypoints"] == 123
    assert "palantir_pseudotime" in adata.obs
    assert "palantir_branch_probs" in adata.obsm


def test_infer_pseudotime_palantir_auto_selects_root_from_first_component(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_pca"] = np.ones((adata.n_obs, 5))
    monkeypatch.setattr(traj, "ensure_pca", lambda *_a, **_k: None)
    captured: dict[str, object] = {}

    eigen = np.zeros((adata.n_obs, 3), dtype=float)
    eigen[1, 0] = 10.0
    multiscale = np.zeros((adata.n_obs, 2), dtype=float)
    multiscale[5, 0] = -4.0

    class _PR:
        pseudotime = pd.Series(np.linspace(0, 1, adata.n_obs), index=adata.obs_names)
        branch_probs = np.ones((adata.n_obs, 2))

    fake_palantir = ModuleType("palantir")

    def _determine_multiscale_space(_result):
        captured["used_multiscale"] = True
        return pd.DataFrame(multiscale, index=adata.obs_names)

    fake_palantir.utils = SimpleNamespace(
        run_diffusion_maps=lambda *_a, **_k: {"EigenVectors": eigen},
        determine_multiscale_space=_determine_multiscale_space,
    )

    def _fake_run_palantir(_ms_data, start_cell, num_waypoints):
        captured["start_cell"] = start_cell
        captured["num_waypoints"] = num_waypoints
        return _PR()

    fake_palantir.core = SimpleNamespace(run_palantir=_fake_run_palantir)

    out = traj.infer_pseudotime_palantir(
        adata,
        root_cells=None,
        num_waypoints=77,
        palantir_module=fake_palantir,
    )
    assert out is adata
    assert captured["used_multiscale"] is True
    assert captured["start_cell"] == adata.obs_names[5]
    assert captured["num_waypoints"] == 77


def test_normalize_palantir_outputs_aligns_labels_and_rejects_invalid_values():
    obs_names = pd.Index(["c1", "c2", "c3"])
    matrix = pd.DataFrame(
        {"fate": [0.3, 0.1, 0.2]}, index=["c3", "c1", "c2"]
    )

    aligned = traj._normalize_palantir_matrix(
        matrix, obs_names, name="branch probabilities"
    )
    np.testing.assert_allclose(aligned["fate"], [0.1, 0.2, 0.3])

    with pytest.raises(DataCompatibilityError, match=r"values in \[0, 1\]"):
        traj._normalize_palantir_pseudotime(
            pd.Series([0.0, 0.5, 1.2], index=obs_names), obs_names
        )

    with pytest.raises(DataCompatibilityError, match="cells do not match"):
        traj._normalize_palantir_matrix(
            matrix.iloc[:-1], obs_names, name="branch probabilities"
        )


def test_compute_dpt_trajectory_valid_root_and_fillna(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(traj, "ensure_pca", lambda *_a, **_k: None)
    monkeypatch.setattr(traj, "ensure_neighbors", lambda *_a, **_k: None)
    monkeypatch.setattr(traj, "ensure_diffmap", lambda *_a, **_k: None)

    fake_scanpy = ModuleType("scanpy")
    fake_scanpy.tl = SimpleNamespace(
        dpt=lambda _a: _a.obs.__setitem__(
            "dpt_pseudotime",
            pd.Series([0.2, np.nan] + [0.1] * (_a.n_obs - 2), index=_a.obs_names),
        )
    )
    monkeypatch.setitem(__import__("sys").modules, "scanpy", fake_scanpy)

    root = adata.obs_names[3]
    out = traj.compute_dpt_trajectory(adata, root_cells=[root])
    assert out is adata
    assert int(adata.uns["iroot"]) == 3
    # NaN pseudotime is preserved (unreachable cells), not filled with 0
    assert adata.obs["dpt_pseudotime"].isna().sum() == 1


def test_compute_dpt_trajectory_raises_when_dpt_column_missing(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(traj, "ensure_pca", lambda *_a, **_k: None)
    monkeypatch.setattr(traj, "ensure_neighbors", lambda *_a, **_k: None)
    monkeypatch.setattr(traj, "ensure_diffmap", lambda *_a, **_k: None)

    fake_scanpy = ModuleType("scanpy")
    fake_scanpy.tl = SimpleNamespace(dpt=lambda *_a, **_k: None)
    monkeypatch.setitem(__import__("sys").modules, "scanpy", fake_scanpy)

    with pytest.raises(ProcessingError, match="did not create 'dpt_pseudotime'"):
        traj.compute_dpt_trajectory(adata)


@pytest.mark.asyncio
async def test_analyze_trajectory_cellrank_success_records_cellrank_specific_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.uns["velocity_graph"] = True
    captured: dict[str, object] = {}

    def _fake_infer(_adata, **_kwargs):
        _adata.obs["pseudotime"] = np.linspace(0, 1, _adata.n_obs)
        _adata.obs["terminal_states"] = pd.Categorical(["T0"] * _adata.n_obs)
        _adata.obs["macrostates"] = pd.Categorical(["M0"] * _adata.n_obs)
        _adata.obsm["fate_probabilities"] = np.ones((_adata.n_obs, 1), dtype=float)
        _adata.uns["velocity_method"] = "scvelo"
        return _adata

    monkeypatch.setattr(traj, "infer_spatial_trajectory_cellrank", _fake_infer)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_a, **_k: [],
    )

    out = await traj.analyze_trajectory(
        "t-cellrank",
        _VelCtx(adata),
        traj.TrajectoryParameters(
            method="cellrank",
            root_cells=[adata.obs_names[0]],
            cellrank_n_states=4,
            cellrank_kernel_weights=(0.7, 0.3),
            spatial_weight=0.25,
        ),
    )

    assert out.method == "cellrank"
    assert out.pseudotime_key == "pseudotime"
    assert "terminal_states_cellrank_sw0_25" in captured["results_keys"]["obs"]
    assert "macrostates_cellrank_sw0_25" in captured["results_keys"]["obs"]
    assert "fate_probabilities_cellrank_sw0_25" in captured["results_keys"]["obsm"]
    assert captured["parameters"]["kernel_weights"] == [0.7, 0.3]
    assert captured["parameters"]["n_states"] == 4
    assert captured["parameters"]["root_cells"] == [adata.obs_names[0]]


@pytest.mark.asyncio
async def test_cellrank_metadata_omits_absent_keys(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """CellRank metadata only registers keys that actually exist in adata."""
    adata = minimal_spatial_adata.copy()
    adata.uns["velocity_graph"] = True
    captured: dict[str, object] = {}

    def _fake_infer_minimal(_adata, **_kwargs):
        # Only pseudotime, no terminal_states / macrostates / fate_probabilities
        _adata.obs["pseudotime"] = np.linspace(0, 1, _adata.n_obs)
        return _adata

    monkeypatch.setattr(traj, "infer_spatial_trajectory_cellrank", _fake_infer_minimal)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_a, **_k: [],
    )

    await traj.analyze_trajectory(
        "t-cellrank-min",
        _VelCtx(adata),
        traj.TrajectoryParameters(method="cellrank"),
    )

    rk = captured["results_keys"]
    assert "terminal_states" not in rk["obs"]
    assert "macrostates" not in rk["obs"]
    assert "fate_probabilities" not in rk["obsm"]
    assert "velocity_method" not in rk["uns"]


@pytest.mark.asyncio
async def test_analyze_trajectory_cellrank_wraps_inference_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.uns["velocity_graph"] = True
    monkeypatch.setattr(
        traj,
        "infer_spatial_trajectory_cellrank",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("cellrank boom")),
    )

    with pytest.raises(ProcessingError, match="CellRank trajectory inference failed"):
        await traj.analyze_trajectory(
            "t-cellrank-fail",
            _VelCtx(adata),
            traj.TrajectoryParameters(method="cellrank"),
        )


@pytest.mark.asyncio
async def test_analyze_trajectory_palantir_wraps_inference_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(
        traj,
        "infer_pseudotime_palantir",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("palantir boom")),
    )

    with pytest.raises(ProcessingError, match="Palantir trajectory inference failed"):
        await traj.analyze_trajectory(
            "t-palantir-fail",
            _VelCtx(adata),
            traj.TrajectoryParameters(method="palantir"),
        )


@pytest.mark.asyncio
async def test_analyze_trajectory_palantir_preserves_dependency_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(
        traj,
        "infer_pseudotime_palantir",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("palantir is unavailable")
        ),
    )

    with pytest.raises(DependencyError, match="palantir is unavailable"):
        await traj.analyze_trajectory(
            "t-palantir-missing",
            _VelCtx(adata),
            traj.TrajectoryParameters(method="palantir"),
        )


@pytest.mark.asyncio
async def test_analyze_trajectory_raises_when_method_returns_without_pseudotime(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(
        traj, "compute_dpt_trajectory", lambda _adata, root_cells=None: _adata
    )

    with pytest.raises(ProcessingError, match="Failed to compute pseudotime"):
        await traj.analyze_trajectory(
            "t-no-pt",
            _VelCtx(adata),
            traj.TrajectoryParameters(method="dpt"),
        )
