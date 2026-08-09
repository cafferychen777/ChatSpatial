"""Additional branch-contract tests for annotation methods."""

from __future__ import annotations

import sys
import warnings
from types import ModuleType, SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from chatspatial.models.data import AnnotationParameters
from chatspatial.tools import annotation as ann
from chatspatial.utils.exceptions import (
    DataError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


def _required_module(name: str, *_args, **_kwargs):
    """Return the fake backend installed by each test."""
    return sys.modules.get(name, object())


class DummyWarnCtx:
    def __init__(self, datasets: dict[str, object] | None = None):
        self.datasets = datasets or {}
        self.warnings: list[str] = []
        self.calls: list[str] = []

    async def get_adata(self, data_id: str):
        self.calls.append(data_id)
        return self.datasets[data_id]

    async def set_adata(self, data_id: str, adata):
        self.datasets[data_id] = adata

    async def warning(self, msg: str):
        self.warnings.append(msg)

    async def info(self, _msg: str):
        return None


@pytest.mark.asyncio
async def test_singler_celldex_reference_without_labels_raises(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    class _Ref:
        def get_column_data(self):
            return self

        def column(self, _name: str):
            raise KeyError("missing")

    fake_celldex = ModuleType("celldex")
    fake_celldex.fetch_reference = lambda *_a, **_k: _Ref()

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(
        ann,
        "get",
        lambda name, *_args, **_kwargs: fake_celldex if name == "celldex" else None,
    )
    monkeypatch.setitem(__import__("sys").modules, "singler", ModuleType("singler"))
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=False: SimpleNamespace(
            X=_adata.X, is_integer_counts=True
        ),
    )

    with pytest.raises(DataNotFoundError, match="Could not find labels"):
        await ann._annotate_with_singler(
            adata,
            AnnotationParameters(method="singler", singler_reference="hpca"),
            DummyWarnCtx(),
            "cell_type_singler",
            "confidence_singler",
            reference_adata=None,
        )


@pytest.mark.asyncio
async def test_singler_uses_normalized_layers_and_warns_for_low_delta_confidence(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    adata.uns["log1p"] = {"base": None}
    ref.uns["log1p"] = {"base": None}
    adata.layers["X_normalized"] = np.asarray(adata.X, dtype=float) + 1.0
    ref.layers["X_normalized"] = np.asarray(ref.X, dtype=float) + 2.0
    ref.obs["ctype"] = pd.Categorical(
        ["B"] * (ref.n_obs // 2) + ["T"] * (ref.n_obs - ref.n_obs // 2)
    )

    captured: dict[str, object] = {}

    class _SingleResults:
        def __init__(self, n_obs: int):
            self._best = ["B" if i % 2 == 0 else "T" for i in range(n_obs)]
            self._scores = pd.DataFrame(
                {
                    "B": [0.7 if i % 2 == 0 else 0.2 for i in range(n_obs)],
                    "T": [0.3 if i % 2 == 0 else 0.8 for i in range(n_obs)],
                }
            )
            self._delta = [0.01 if i < int(n_obs * 0.5) else 0.2 for i in range(n_obs)]

        def column(self, name: str):
            if name == "best":
                return self._best
            if name == "scores":
                return self._scores
            if name == "delta":
                return self._delta
            raise KeyError(name)

    fake_singler = ModuleType("singler")

    def _annotate_single(**kwargs):
        captured["kwargs"] = kwargs
        return _SingleResults(adata.n_obs)

    fake_singler.annotate_single = _annotate_single
    monkeypatch.setitem(__import__("sys").modules, "singler", fake_singler)

    async def _fake_ensure_unique_var_names_async(_ad, _ctx, label: str):
        return 1 if label == "query data" else 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "get", lambda *_a, **_k: None)
    monkeypatch.setattr(
        ann, "ensure_unique_var_names_async", _fake_ensure_unique_var_names_async
    )
    monkeypatch.setattr(
        ann,
        "find_common_genes",
        lambda test_features, ref_features: [
            g for g in test_features if g in set(ref_features)
        ],
    )

    ctx = DummyWarnCtx()
    out = await ann._annotate_with_singler(
        adata,
        AnnotationParameters(
            method="singler",
            reference_data_id="r",
            cell_type_key="ctype",
            n_threads=3,
            singler_fine_tune=False,
        ),
        ctx,
        "cell_type_singler",
        "confidence_singler",
        reference_adata=ref,
    )

    assert np.array_equal(
        captured["kwargs"]["test_data"], adata.layers["X_normalized"].T
    )
    assert np.array_equal(captured["kwargs"]["ref_data"], ref.layers["X_normalized"].T)
    assert captured["kwargs"]["num_threads"] == 3
    assert captured["kwargs"]["classify_args"] == {"use_fine_tune": False}
    assert "n_threads" not in captured["kwargs"]
    assert set(out.cell_types) == {"B", "T"}
    assert any("low confidence scores" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_singler_delta_extraction_failure_falls_back_to_score_confidence(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(
        ["B"] * (ref.n_obs // 2) + ["T"] * (ref.n_obs - ref.n_obs // 2)
    )

    class _SingleResults:
        def __init__(self, n_obs: int):
            self._best = ["B" if i % 2 == 0 else "T" for i in range(n_obs)]
            self._scores = pd.DataFrame(
                {
                    "B": [0.9 if i % 2 == 0 else -0.1 for i in range(n_obs)],
                    "T": [0.1 if i % 2 == 0 else 0.8 for i in range(n_obs)],
                }
            )

        def column(self, name: str):
            if name == "best":
                return self._best
            if name == "scores":
                return self._scores
            if name == "delta":
                return 1  # non-subscriptable -> delta processing exception path
            raise KeyError(name)

    fake_singler = ModuleType("singler")
    fake_singler.annotate_single = lambda **_kwargs: _SingleResults(adata.n_obs)
    monkeypatch.setitem(__import__("sys").modules, "singler", fake_singler)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "get", lambda *_a, **_k: None)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(
        ann,
        "find_common_genes",
        lambda test_features, ref_features: [
            g for g in test_features if g in set(ref_features)
        ],
    )
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=False: SimpleNamespace(
            X=_adata.X, is_integer_counts=True
        ),
    )

    out = await ann._annotate_with_singler(
        adata,
        AnnotationParameters(
            method="singler", reference_data_id="r", cell_type_key="ctype"
        ),
        DummyWarnCtx(),
        "cell_type_singler",
        "confidence_singler",
        reference_adata=ref,
    )

    assert "B" in out.confidence and "T" in out.confidence


@pytest.mark.asyncio
async def test_singler_delta_confidence_computation_failure_falls_back_to_scores(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(
        ["B"] * (ref.n_obs // 2) + ["T"] * (ref.n_obs - ref.n_obs // 2)
    )

    class _DeltaToken:
        def __init__(self, value: float):
            self.value = value

        def __bool__(self):
            return True

        def __lt__(self, other):
            return self.value < other

    class _SingleResults:
        def __init__(self, n_obs: int):
            self._best = ["B" if i % 2 == 0 else "T" for i in range(n_obs)]
            self._scores = pd.DataFrame(
                {
                    "B": [0.8 if i % 2 == 0 else 0.1 for i in range(n_obs)],
                    "T": [0.2 if i % 2 == 0 else 0.9 for i in range(n_obs)],
                }
            )
            self._delta = [_DeltaToken(0.1)] * n_obs

        def column(self, name: str):
            if name == "best":
                return self._best
            if name == "scores":
                return self._scores
            if name == "delta":
                return self._delta
            raise KeyError(name)

    fake_singler = ModuleType("singler")
    fake_singler.annotate_single = lambda **_kwargs: _SingleResults(adata.n_obs)
    monkeypatch.setitem(__import__("sys").modules, "singler", fake_singler)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "get", lambda *_a, **_k: None)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(
        ann,
        "find_common_genes",
        lambda test_features, ref_features: [
            g for g in test_features if g in set(ref_features)
        ],
    )
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=False: SimpleNamespace(
            X=_adata.X, is_integer_counts=True
        ),
    )

    out = await ann._annotate_with_singler(
        adata,
        AnnotationParameters(
            method="singler", reference_data_id="r", cell_type_key="ctype"
        ),
        DummyWarnCtx(),
        "cell_type_singler",
        "confidence_singler",
        reference_adata=ref,
    )

    assert out.confidence
    assert "B" in out.confidence and "T" in out.confidence


@pytest.mark.asyncio
async def test_tangram_clusters_mode_without_detectable_cluster_label_raises(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_a, **_k: None
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_cell_type_key", lambda *_a, **_k: None)
    monkeypatch.setattr(ann, "get_cluster_key", lambda *_a, **_k: None)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    with pytest.raises(ParameterError, match="No cluster label found"):
        await ann._annotate_with_tangram(
            adata,
            AnnotationParameters(
                method="tangram",
                tangram_mode="clusters",
                cluster_label=None,
                training_genes=["gene_0", "gene_1"],
            ),
            DummyWarnCtx(),
            "cell_type_tangram",
            "confidence_tangram",
            reference_adata=ref,
        )


@pytest.mark.asyncio
async def test_tangram_optional_mapping_args_validation_warnings_and_low_score_warning(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.raw = adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(
        ["B"] * (ref.n_obs // 2) + ["T"] * (ref.n_obs - ref.n_obs // 2)
    )

    captured: dict[str, object] = {}

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_a, **_k: None

    def _map_cells_to_space(_adata_sc, _adata_sp, **kwargs):
        captured["kwargs"] = kwargs
        return SimpleNamespace(uns={"training_history": {"main_loss": ["bad-score"]}})

    def _project_cell_annotations(_ad_map, adata_sp, annotation):
        captured["annotation"] = annotation
        pattern = np.array([[2.0, -1.0], [0.0, 0.0], [-1.0, 2.0]], dtype=float)
        vals = np.vstack([pattern[i % 3] for i in range(adata_sp.n_obs)])
        adata_sp.obsm["tangram_ct_pred"] = pd.DataFrame(
            vals, columns=["B", "T"], index=adata_sp.obs_names
        )

    fake_tg.map_cells_to_space = _map_cells_to_space
    fake_tg.compare_spatial_geneexp = lambda *_a, **_k: (_ for _ in ()).throw(
        RuntimeError("validation failed")
    )
    fake_tg.project_genes = lambda *_a, **_k: (_ for _ in ()).throw(
        RuntimeError("project genes failed")
    )
    fake_tg.project_cell_annotations = _project_cell_annotations
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    ctx = DummyWarnCtx()
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="The behavior of DataFrame.idxmax with all-NA values",
            category=FutureWarning,
        )
        out = await ann._annotate_with_tangram(
            adata,
            AnnotationParameters(
                method="tangram",
                tangram_mode="clusters",
                cluster_label="ctype",
                tangram_use_gpu=True,
                tangram_lambda_r=0.4,
                tangram_lambda_neighborhood=0.2,
                tangram_compute_validation=True,
                tangram_project_genes=True,
                training_genes=["gene_0", "gene_1"],
            ),
            ctx,
            "cell_type_tangram",
            "confidence_tangram",
            reference_adata=ref,
        )

    assert captured["kwargs"]["lambda_r"] == 0.4
    assert captured["kwargs"]["lambda_neighborhood"] == 0.2
    assert captured["kwargs"]["cluster_label"] == "ctype"
    assert captured["annotation"] == "ctype"
    assert out.tangram_mapping_score == 0.0
    assert any("GPU requested but not available" in msg for msg in ctx.warnings)
    assert any("Could not compute validation metrics" in msg for msg in ctx.warnings)
    assert any("Could not project genes" in msg for msg in ctx.warnings)
    assert any("normalized probabilities are negative" in msg for msg in ctx.warnings)
    assert any("normalized probabilities exceed 1.0" in msg for msg in ctx.warnings)
    assert any("score is suspiciously low" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_tangram_score_extraction_unknown_history_format_raises(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_a, **_k: None
    fake_tg.map_cells_to_space = lambda *_a, **_k: SimpleNamespace(
        uns={"training_history": [1, 2, 3]}
    )
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    with pytest.raises(ProcessingError, match="history format not recognized"):
        await ann._annotate_with_tangram(
            adata,
            AnnotationParameters(
                method="tangram",
                cell_type_key="ctype",
                training_genes=["gene_0", "gene_1"],
            ),
            DummyWarnCtx(),
            "cell_type_tangram",
            "confidence_tangram",
            reference_adata=ref,
        )


@pytest.mark.asyncio
async def test_tangram_copies_gene_predictions_back_to_original_adata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.raw = adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_a, **_k: None
    fake_tg.map_cells_to_space = lambda *_a, **_k: SimpleNamespace(
        uns={"training_history": {"main_loss": [1.0]}}
    )
    fake_tg.project_genes = lambda *_a, **_k: SimpleNamespace(
        X=np.ones((adata.n_obs, adata.n_vars), dtype=float)
    )
    fake_tg.project_cell_annotations = (
        lambda _ad_map, adata_sp, annotation: adata_sp.obsm.__setitem__(
            "tangram_ct_pred",
            pd.DataFrame(
                {"B": np.ones(adata_sp.n_obs, dtype=float)}, index=adata_sp.obs_names
            ),
        )
    )
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    await ann._annotate_with_tangram(
        adata,
        AnnotationParameters(
            method="tangram",
            cell_type_key="ctype",
            tangram_project_genes=True,
            training_genes=["gene_0", "gene_1"],
        ),
        DummyWarnCtx(),
        "cell_type_tangram",
        "confidence_tangram",
        reference_adata=ref,
    )

    # Default tangram_ct_pred_key is "tangram_ct_pred", so gene_pred_key
    # is "tangram_gene_predictions" (no suffix when called directly).
    assert "tangram_gene_predictions" in adata.obsm


@pytest.mark.asyncio
async def test_tangram_gene_predictions_uses_parametrized_key_with_suffix(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """When tangram_ct_pred_key has a suffix, gene_predictions key follows."""
    adata = minimal_spatial_adata.copy()
    adata.raw = adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_a, **_k: None
    fake_tg.map_cells_to_space = lambda *_a, **_k: SimpleNamespace(
        uns={"training_history": {"main_loss": [1.0]}}
    )
    fake_tg.project_genes = lambda *_a, **_k: SimpleNamespace(
        X=np.ones((adata.n_obs, adata.n_vars), dtype=float)
    )
    fake_tg.project_cell_annotations = (
        lambda _ad_map, adata_sp, annotation: adata_sp.obsm.__setitem__(
            "tangram_ct_pred",
            pd.DataFrame(
                {"B": np.ones(adata_sp.n_obs, dtype=float)},
                index=adata_sp.obs_names,
            ),
        )
    )
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    await ann._annotate_with_tangram(
        adata,
        AnnotationParameters(
            method="tangram",
            cell_type_key="ctype",
            tangram_project_genes=True,
            training_genes=["gene_0", "gene_1"],
        ),
        DummyWarnCtx(),
        "cell_type_tangram_myref",
        "confidence_tangram_myref",
        reference_adata=ref,
        tangram_ct_pred_key="tangram_ct_pred_tangram_myref",
    )

    # The gene_predictions key should use the same suffix
    assert "tangram_gene_predictions_tangram_myref" in adata.obsm
    # The shared key should NOT be present (parametrized replaces it)
    assert "tangram_gene_predictions" not in adata.obsm


@pytest.mark.asyncio
async def test_scanvi_pretrain_subsetting_and_missing_batch_key_are_handled(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(
        ["B"] * (ref.n_obs // 2) + ["T"] * (ref.n_obs - ref.n_obs // 2)
    )

    seen: dict[str, object] = {}

    class _FakeSCVI:
        @staticmethod
        def setup_anndata(*_a, **_k):
            seen["scvi_setup"] = _k

        def __init__(self, adata_obj, **_kwargs):
            self.adata_obj = adata_obj

        def train(self, **_kwargs):
            seen["scvi_train"] = _kwargs

    class _FakeSCANVI:
        @staticmethod
        def setup_anndata(*_a, **_k):
            seen.setdefault("scanvi_setups", []).append(_k)

        def __init__(self, adata_obj, **_kwargs):
            self._adata = adata_obj

        @classmethod
        def from_scvi_model(cls, scvi_model, _unlabeled):
            return cls(scvi_model.adata_obj)

        @staticmethod
        def load_query_data(adata_subset, _model):
            seen["query_batch_values"] = list(adata_subset.obs["batch"].unique())
            return _FakeSCANVI(adata_subset)

        def train(self, **_kwargs):
            seen.setdefault("scanvi_trains", []).append(_kwargs)

        def predict(self, soft: bool = False):
            if soft:
                return pd.DataFrame(
                    {
                        "B": [
                            0.9 if i % 2 == 0 else 0.1 for i in range(self._adata.n_obs)
                        ],
                        "T": [
                            0.1 if i % 2 == 0 else 0.9 for i in range(self._adata.n_obs)
                        ],
                    },
                    index=self._adata.obs_names,
                )
            return pd.Categorical(
                ["B" if i % 2 == 0 else "T" for i in range(self._adata.n_obs)]
            )

    fake_scvi = SimpleNamespace(
        model=SimpleNamespace(SCVI=_FakeSCVI, SCANVI=_FakeSCANVI)
    )

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "validate_scvi_tools", lambda *_a, **_k: fake_scvi)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(
        ann,
        "find_common_genes",
        lambda ref_genes, _qry_genes: list(ref_genes)[:15],
    )
    monkeypatch.setattr(ann, "ensure_counts_layer", lambda *_a, **_k: None)
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x.copy())

    ctx = DummyWarnCtx()
    out = await ann._annotate_with_scanvi(
        adata,
        AnnotationParameters(
            method="scanvi",
            reference_data_id="r",
            cell_type_key="ctype",
            scanvi_use_scvi_pretrain=True,
            batch_key="batch",
            scanvi_query_epochs=2,
            scanvi_scvi_epochs=2,
            scanvi_scanvi_epochs=2,
        ),
        ctx,
        "cell_type_scanvi",
        "confidence_scanvi",
        reference_adata=ref,
    )

    assert set(out.cell_types) == {"B", "T"}
    assert seen["query_batch_values"] == ["query_batch"]
    assert any("Subsetting to 15 common genes" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_cellassign_non_raw_source_data_cleaning_and_index_predictions(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.X = sp.csr_matrix(np.asarray(adata.X, dtype=float))

    class _FakeCellAssign:
        @staticmethod
        def setup_anndata(*_a, **_k):
            return None

        def __init__(self, adata_subset, marker_gene_matrix, **_kwargs):
            self._adata = adata_subset
            self._n_types = len(marker_gene_matrix.columns)

        def train(self, **_kwargs):
            return None

        def predict(self):
            return np.array([i % self._n_types for i in range(self._adata.n_obs)])

    fake_scvi_external = ModuleType("scvi.external")
    fake_scvi_external.CellAssign = _FakeCellAssign
    monkeypatch.setitem(__import__("sys").modules, "scvi.external", fake_scvi_external)

    def _fake_to_dense(X):
        dense = np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=float)
        if dense.shape[1] >= 3:
            dense[:, 0] = 1.0  # zero-variance gene
            dense[0, 1] = np.nan
            dense[1, 1] = np.inf
            dense[2, 2] = -5.0
        return dense

    monkeypatch.setattr(
        ann,
        "validate_scvi_tools",
        lambda *_a, **_k: SimpleNamespace(external=fake_scvi_external),
    )
    monkeypatch.setattr(ann, "to_dense", _fake_to_dense)
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True, require_integer_counts=False: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="X",
        ),
    )

    params = AnnotationParameters(
        method="cellassign",
        marker_genes={
            "B": ["gene_0", "MISSING_1", "MISSING_2", "MISSING_3"],
            "T": ["gene_1", "gene_2"],
        },
    )

    ctx = DummyWarnCtx()
    out = await ann._annotate_with_cellassign(
        adata,
        params,
        ctx,
        "cell_type_cellassign",
        "confidence_cellassign",
    )

    assert set(out.cell_types) == {"B", "T"}
    assert out.confidence == {}
    assert any("Using X data for marker gene validation" in msg for msg in ctx.warnings)
    assert any("Missing most markers for B" in msg for msg in ctx.warnings)
    assert any("zero variance" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_cellassign_retries_without_n_hidden_for_newer_scvi_api(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    init_kwargs = []

    class _FakeCellAssign:
        @staticmethod
        def setup_anndata(*_a, **_k):
            return None

        def __init__(self, adata_subset, marker_gene_matrix, **kwargs):
            init_kwargs.append(kwargs)
            if "n_hidden" in kwargs:
                raise TypeError(
                    "CellAssignModule.__init__() got an unexpected keyword argument 'n_hidden'"
                )
            self._adata = adata_subset
            self._labels = list(reversed(marker_gene_matrix.columns.tolist()))

        def train(self, **_kwargs):
            return None

        def predict(self):
            values = np.zeros((self._adata.n_obs, len(self._labels)))
            for i in range(self._adata.n_obs):
                values[i, i % len(self._labels)] = 1.0
            return pd.DataFrame(
                values, index=self._adata.obs_names, columns=self._labels
            )

    fake_scvi_external = ModuleType("scvi.external")
    fake_scvi_external.CellAssign = _FakeCellAssign
    monkeypatch.setitem(__import__("sys").modules, "scvi.external", fake_scvi_external)
    monkeypatch.setattr(
        ann,
        "validate_scvi_tools",
        lambda *_a, **_k: SimpleNamespace(external=fake_scvi_external),
    )
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True, require_integer_counts=False: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="X",
        ),
    )

    params = AnnotationParameters(
        method="cellassign",
        marker_genes={"B": ["gene_0"], "T": ["gene_1"]},
    )

    ctx = DummyWarnCtx()
    out = await ann._annotate_with_cellassign(
        adata,
        params,
        ctx,
        "cell_type_cellassign",
        "confidence_cellassign",
    )

    assert out.cell_types == ["B", "T"]
    assert set(adata.obs["cell_type_cellassign"]) == {"B", "T"}
    assert init_kwargs == [{"n_hidden": params.cellassign_n_hidden}, {}]
    assert any("does not accept n_hidden" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_cellassign_rejects_unknown_prediction_labels(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    class _FakeCellAssign:
        @staticmethod
        def setup_anndata(*_a, **_k):
            return None

        def __init__(self, adata_subset, _marker_gene_matrix, **_kwargs):
            self._adata = adata_subset

        def train(self, **_kwargs):
            return None

        def predict(self):
            return pd.DataFrame(
                {"Unknown": np.ones(self._adata.n_obs)}, index=self._adata.obs_names
            )

    fake_scvi_external = ModuleType("scvi.external")
    fake_scvi_external.CellAssign = _FakeCellAssign
    monkeypatch.setitem(__import__("sys").modules, "scvi.external", fake_scvi_external)
    monkeypatch.setattr(
        ann,
        "validate_scvi_tools",
        lambda *_a, **_k: SimpleNamespace(external=fake_scvi_external),
    )
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True, require_integer_counts=False: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="X",
        ),
    )

    with pytest.raises(ann.DataError, match="not present in the marker matrix"):
        await ann._annotate_with_cellassign(
            adata,
            AnnotationParameters(
                method="cellassign", marker_genes={"B": ["gene_0"], "T": ["gene_1"]}
            ),
            DummyWarnCtx(),
            "cell_type_cellassign",
            "confidence_cellassign",
        )


@pytest.mark.asyncio
async def test_mllmcelltype_consensus_path_handles_unmapped_clusters(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(
        ["0"] * 20 + ["1"] * 20 + ["2"] * (adata.n_obs - 40)
    )

    def _fake_rank_genes_groups(adata_obj, *_a, **_k):
        adata_obj.uns["rank_genes_groups"] = {
            "names": {
                "0": np.array(["gene_0", "gene_1"], dtype=object),
                "1": np.array(["gene_2", "gene_3"], dtype=object),
                "2": np.array(["gene_4", "gene_5"], dtype=object),
            }
        }

    fake_mllm = ModuleType("mllmcelltype")
    fake_mllm.interactive_consensus_annotation = lambda **_kwargs: {
        "consensus": {"Cluster_0": "T", "Cluster_1": "B"}
    }
    monkeypatch.setitem(__import__("sys").modules, "mllmcelltype", fake_mllm)
    monkeypatch.setattr(ann.sc.tl, "rank_genes_groups", _fake_rank_genes_groups)
    monkeypatch.setattr(ann, "require", _required_module)

    ctx = DummyWarnCtx()
    out = await ann._annotate_with_mllmcelltype(
        adata,
        AnnotationParameters(
            method="mllmcelltype",
            cluster_label="leiden",
            mllm_use_consensus=True,
            mllm_models=["gpt-5"],
        ),
        ctx,
        "cell_type_mllmcelltype",
        "confidence_mllmcelltype",
    )

    assert "Unknown" in out.counts
    assert any("unmapped clusters" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("backend_error", "expected_error", "message"),
    [
        (
            RuntimeError("api boom"),
            ProcessingError,
            "mLLMCellType annotation failed: api boom",
        ),
        (
            DependencyError("credentials missing"),
            DependencyError,
            "credentials missing",
        ),
    ],
)
async def test_mllmcelltype_single_model_error_contracts(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
    backend_error,
    expected_error,
    message,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(
        ["0"] * (adata.n_obs // 2) + ["1"] * (adata.n_obs - adata.n_obs // 2)
    )

    def _fake_rank_genes_groups(adata_obj, *_a, **_k):
        adata_obj.uns["rank_genes_groups"] = {
            "names": {
                "0": np.array(["gene_0", "gene_1"], dtype=object),
                "1": np.array(["gene_2", "gene_3"], dtype=object),
            }
        }

    fake_mllm = ModuleType("mllmcelltype")
    fake_mllm.annotate_clusters = lambda **_kwargs: (_ for _ in ()).throw(backend_error)
    monkeypatch.setitem(__import__("sys").modules, "mllmcelltype", fake_mllm)
    monkeypatch.setattr(ann.sc.tl, "rank_genes_groups", _fake_rank_genes_groups)
    monkeypatch.setattr(ann, "require", _required_module)

    with pytest.raises(expected_error, match=message):
        await ann._annotate_with_mllmcelltype(
            adata,
            AnnotationParameters(method="mllmcelltype", cluster_label="leiden"),
            DummyWarnCtx(),
            "cell_type_mllmcelltype",
            "confidence_mllmcelltype",
        )


@pytest.mark.asyncio
async def test_annotate_cell_types_dispatch_metadata_for_tangram_scanvi_and_mllm(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ctx = DummyWarnCtx({"q": adata, "r": ref})
    captured_calls: list[dict[str, object]] = []

    def _capture_metadata(_adata, **kwargs):
        captured_calls.append(kwargs)

    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata", _capture_metadata
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result", lambda *_a, **_k: []
    )

    async def _fake_tangram(
        _adata,
        _params,
        _ctx,
        output_key,
        confidence_key,
        reference_adata,
        tangram_ct_pred_key="tangram_ct_pred",
    ):
        del reference_adata
        _adata.obs[output_key] = ["B"] * _adata.n_obs
        _adata.obs[confidence_key] = [0.9] * _adata.n_obs
        _adata.obsm[tangram_ct_pred_key] = pd.DataFrame(
            {"B": np.ones(_adata.n_obs)}, index=_adata.obs_names
        )
        return ann.AnnotationMethodOutput(
            cell_types=["B"],
            counts={"B": _adata.n_obs},
            confidence={"B": 0.9},
            tangram_mapping_score=0.77,
        )

    async def _fake_scanvi(
        _adata, _params, _ctx, output_key, confidence_key, reference_adata
    ):
        del reference_adata
        _adata.obs[output_key] = ["T"] * _adata.n_obs
        _adata.obs[confidence_key] = [0.8] * _adata.n_obs
        return ann.AnnotationMethodOutput(
            cell_types=["T"],
            counts={"T": _adata.n_obs},
            confidence={"T": 0.8},
            tangram_mapping_score=None,
        )

    async def _fake_mllm(_adata, _params, _ctx, output_key, confidence_key):
        del confidence_key
        _adata.obs[output_key] = ["M"] * _adata.n_obs
        return ann.AnnotationMethodOutput(
            cell_types=["M"],
            counts={"M": _adata.n_obs},
            confidence={},
            tangram_mapping_score=None,
        )

    monkeypatch.setattr(ann, "_annotate_with_tangram", _fake_tangram)
    monkeypatch.setattr(ann, "_annotate_with_scanvi", _fake_scanvi)
    monkeypatch.setattr(ann, "_annotate_with_mllmcelltype", _fake_mllm)

    out_t = await ann.annotate_cell_types(
        "q",
        ctx,
        AnnotationParameters(method="tangram", reference_data_id="r"),
    )
    out_s = await ann.annotate_cell_types(
        "q",
        ctx,
        AnnotationParameters(method="scanvi", reference_data_id="r"),
    )
    out_m = await ann.annotate_cell_types(
        "q",
        ctx,
        AnnotationParameters(method="mllmcelltype"),
    )

    assert out_t.tangram_mapping_score == 0.77
    # scanvi with reference_data_id="r" → suffix "scanvi_r"
    assert out_s.confidence_key == "confidence_scanvi_r"
    assert out_m.confidence_key is None

    tangram_meta = next(c for c in captured_calls if c["method"] == "tangram")
    scanvi_meta = next(c for c in captured_calls if c["method"] == "scanvi")
    mllm_meta = next(c for c in captured_calls if c["method"] == "mllmcelltype")

    assert tangram_meta["results_keys"]["obsm"] == ["tangram_ct_pred_tangram_r"]
    assert "learning_rate" in tangram_meta["parameters"]
    assert "mapping_score" in tangram_meta["statistics"]
    assert "n_latent" in scanvi_meta["parameters"]
    assert "n_marker_genes" in mllm_meta["parameters"]


@pytest.mark.asyncio
async def test_singler_default_celldex_reference_path_runs_successfully(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    class _Ref:
        def get_column_data(self):
            return self

        def column(self, name: str):
            if name == "label.main":
                return ["B"] * adata.n_obs
            raise KeyError(name)

    class _SingleResults:
        def column(self, name: str):
            if name == "best":
                return ["B"] * adata.n_obs
            if name == "scores":
                return pd.DataFrame({"B": [0.8] * adata.n_obs})
            if name == "delta":
                return [0.2] * adata.n_obs
            raise KeyError(name)

    fake_celldex = ModuleType("celldex")
    fake_celldex.fetch_reference = lambda *_a, **_k: _Ref()
    fake_singler = ModuleType("singler")
    fake_singler.annotate_single = lambda **_kwargs: _SingleResults()

    monkeypatch.setitem(__import__("sys").modules, "singler", fake_singler)
    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(
        ann,
        "get",
        lambda name, *_args, **_kwargs: fake_celldex if name == "celldex" else None,
    )
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=False: SimpleNamespace(
            X=_adata.X, is_integer_counts=True
        ),
    )

    out = await ann._annotate_with_singler(
        adata,
        AnnotationParameters(method="singler"),
        DummyWarnCtx(),
        "cell_type_singler",
        "confidence_singler",
        reference_adata=None,
    )

    assert out.cell_types == ["B"]


@pytest.mark.asyncio
async def test_singler_celldex_label_fallback_and_integrated_branch(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    class _Ref(list):
        def get_column_data(self):
            return self

        def column(self, name: str):
            if name == "label.main":
                raise KeyError("missing main label")
            if name == "label.fine":
                return ["T"] * adata.n_obs
            raise KeyError(name)

    class _Integrated:
        def column(self, name: str):
            if name == "best_label":
                return ["T"] * adata.n_obs
            if name == "scores":
                return pd.DataFrame({"T": [0.9] * adata.n_obs})
            raise KeyError(name)

    fake_celldex = ModuleType("celldex")
    fake_celldex.fetch_reference = lambda *_a, **_k: _Ref([1, 2, 3])
    captured: dict[str, object] = {}
    fake_singler = ModuleType("singler")

    def _annotate_integrated(*_args, **kwargs):
        captured.update(kwargs)
        return None, _Integrated()

    fake_singler.annotate_integrated = _annotate_integrated
    monkeypatch.setitem(__import__("sys").modules, "singler", fake_singler)
    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(
        ann,
        "get",
        lambda name, *_args, **_kwargs: fake_celldex if name == "celldex" else None,
    )
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=False: SimpleNamespace(
            X=_adata.X, is_integer_counts=True
        ),
    )

    out = await ann._annotate_with_singler(
        adata,
        AnnotationParameters(
            method="singler",
            singler_integrated=True,
            singler_reference="hpca",
            n_threads=2,
            singler_fine_tune=False,
        ),
        DummyWarnCtx(),
        "cell_type_singler",
        "confidence_singler",
        reference_adata=None,
    )

    assert out.cell_types == ["T"]
    assert out.confidence == {"T": 0.9}
    assert captured["num_threads"] == 2
    assert captured["classify_single_args"] == {"use_fine_tune": False}
    assert captured["classify_integrated_args"] == {"use_fine_tune": False}
    assert "n_threads" not in captured


@pytest.mark.asyncio
async def test_singler_custom_reference_requires_cell_type_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "get", lambda *_a, **_k: None)
    monkeypatch.setitem(__import__("sys").modules, "singler", ModuleType("singler"))
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(
        ann,
        "find_common_genes",
        lambda test_features, ref_features: [
            g for g in test_features if g in set(ref_features)
        ],
    )
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=False: SimpleNamespace(
            X=_adata.X, is_integer_counts=True
        ),
    )

    params = AnnotationParameters(method="singler", reference_data_id="r").model_copy(
        update={"cell_type_key": None}
    )
    with pytest.raises(ParameterError, match="cell_type_key is required"):
        await ann._annotate_with_singler(
            adata,
            params,
            DummyWarnCtx(),
            "cell_type_singler",
            "confidence_singler",
            reference_adata=ref,
        )


@pytest.mark.asyncio
async def test_singler_score_fallback_handles_non_dataframe_score_objects(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(
        [0] * (ref.n_obs // 2) + [1] * (ref.n_obs - ref.n_obs // 2)
    )

    class _ScoresArray:
        def __init__(self, arr):
            self._arr = arr

        def to_numpy(self):
            return self._arr

    class _SingleResults:
        def __init__(self, n_obs: int):
            self._best = [0 if i % 2 == 0 else 1 for i in range(n_obs)]
            self._scores = _ScoresArray(
                np.array(
                    [[0.9, 0.1] if i % 2 == 0 else [0.2, 0.8] for i in range(n_obs)],
                    dtype=float,
                )
            )

        def column(self, name: str):
            if name == "best":
                return self._best
            if name == "scores":
                return self._scores
            if name == "delta":
                raise KeyError("no delta")
            raise KeyError(name)

    fake_singler = ModuleType("singler")
    fake_singler.annotate_single = lambda **_kwargs: _SingleResults(adata.n_obs)
    monkeypatch.setitem(__import__("sys").modules, "singler", fake_singler)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "get", lambda *_a, **_k: None)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(
        ann,
        "find_common_genes",
        lambda test_features, ref_features: [
            g for g in test_features if g in set(ref_features)
        ],
    )
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=False: SimpleNamespace(
            X=_adata.X, is_integer_counts=True
        ),
    )

    out = await ann._annotate_with_singler(
        adata,
        AnnotationParameters(
            method="singler", reference_data_id="r", cell_type_key="ctype"
        ),
        DummyWarnCtx(),
        "cell_type_singler",
        "confidence_singler",
        reference_adata=ref,
    )

    assert out.confidence


@pytest.mark.asyncio
async def test_tangram_requires_hvg_when_training_genes_and_markers_not_provided(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    if "highly_variable" in ref.var.columns:
        del ref.var["highly_variable"]

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_a, **_k: None
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    with pytest.raises(DataNotFoundError, match="HVGs not found"):
        await ann._annotate_with_tangram(
            adata,
            AnnotationParameters(
                method="tangram", training_genes=None, marker_genes=None
            ),
            DummyWarnCtx(),
            "cell_type_tangram",
            "confidence_tangram",
            reference_adata=ref,
        )


@pytest.mark.asyncio
async def test_tangram_uses_marker_genes_when_training_genes_not_given(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)
    seen: dict[str, object] = {}

    fake_tg = ModuleType("tangram")

    def _pp_adatas(_adata_sc, _adata_sp, genes):
        seen["genes"] = genes

    fake_tg.pp_adatas = _pp_adatas
    fake_tg.map_cells_to_space = lambda *_a, **_k: SimpleNamespace(
        uns={"training_history": {"main_loss": [1.0]}}
    )
    fake_tg.project_cell_annotations = (
        lambda _ad_map, adata_sp, annotation: adata_sp.obsm.__setitem__(
            "tangram_ct_pred",
            pd.DataFrame(
                {"B": np.ones(adata_sp.n_obs, dtype=float)}, index=adata_sp.obs_names
            ),
        )
    )
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    await ann._annotate_with_tangram(
        adata,
        AnnotationParameters(
            method="tangram",
            cell_type_key="ctype",
            training_genes=None,
            marker_genes={"B": ["gene_0", "gene_1"], "T": ["gene_1", "gene_2"]},
        ),
        DummyWarnCtx(),
        "cell_type_tangram",
        "confidence_tangram",
        reference_adata=ref,
    )

    assert seen["genes"] == ["gene_0", "gene_1", "gene_2"]


@pytest.mark.asyncio
async def test_tangram_uses_reference_hvgs_when_training_genes_not_provided(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)
    ref.var["highly_variable"] = False
    ref.var.loc[ref.var_names[:2], "highly_variable"] = True
    seen: dict[str, object] = {}

    fake_tg = ModuleType("tangram")

    def _pp_adatas(_adata_sc, _adata_sp, genes):
        seen["genes"] = genes

    fake_tg.pp_adatas = _pp_adatas
    fake_tg.map_cells_to_space = lambda *_a, **_k: SimpleNamespace(
        uns={"training_history": {"main_loss": [1.0]}}
    )
    fake_tg.project_cell_annotations = (
        lambda _ad_map, adata_sp, annotation: adata_sp.obsm.__setitem__(
            "tangram_ct_pred",
            pd.DataFrame(
                {"B": np.ones(adata_sp.n_obs, dtype=float)}, index=adata_sp.obs_names
            ),
        )
    )
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    await ann._annotate_with_tangram(
        adata,
        AnnotationParameters(
            method="tangram",
            cell_type_key="ctype",
            training_genes=None,
            marker_genes=None,
        ),
        DummyWarnCtx(),
        "cell_type_tangram",
        "confidence_tangram",
        reference_adata=ref,
    )

    assert set(seen["genes"]) == set(ref.var_names[:2])


@pytest.mark.asyncio
async def test_tangram_stores_validation_scores_on_success(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_a, **_k: None
    fake_tg.map_cells_to_space = lambda *_a, **_k: SimpleNamespace(
        uns={"training_history": {"main_loss": [1.0]}}
    )
    fake_tg.compare_spatial_geneexp = lambda *_a, **_k: {"pearson": 0.9}
    fake_tg.project_cell_annotations = (
        lambda _ad_map, adata_sp, annotation: adata_sp.obsm.__setitem__(
            "tangram_ct_pred",
            pd.DataFrame(
                {"B": np.ones(adata_sp.n_obs, dtype=float)}, index=adata_sp.obs_names
            ),
        )
    )
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    await ann._annotate_with_tangram(
        adata,
        AnnotationParameters(
            method="tangram",
            cell_type_key="ctype",
            tangram_compute_validation=True,
            training_genes=["gene_0"],
        ),
        DummyWarnCtx(),
        "cell_type_tangram",
        "confidence_tangram",
        reference_adata=ref,
    )

    assert "tangram_validation_scores" in adata.uns


@pytest.mark.asyncio
async def test_tangram_validates_cell_type_key_before_mapping(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["available_col"] = pd.Categorical(["B"] * ref.n_obs)

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_a, **_k: None
    mapping_calls = 0

    def _map_cells_to_space(*_args, **_kwargs):
        nonlocal mapping_calls
        mapping_calls += 1
        return SimpleNamespace(uns={"training_history": {"main_loss": [1.0]}})

    fake_tg.map_cells_to_space = _map_cells_to_space
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", _required_module)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    with pytest.raises(
        ParameterError, match="Cell type column 'missing_key' not found"
    ):
        await ann._annotate_with_tangram(
            adata,
            AnnotationParameters(
                method="tangram",
                cell_type_key="missing_key",
                training_genes=["gene_0"],
            ),
            DummyWarnCtx(),
            "cell_type_tangram",
            "confidence_tangram",
            reference_adata=ref,
        )

    assert mapping_calls == 0


@pytest.mark.asyncio
async def test_scanvi_raises_on_insufficient_gene_overlap(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    fake_scvi = SimpleNamespace(model=SimpleNamespace(SCANVI=object))
    monkeypatch.setattr(ann, "validate_scvi_tools", lambda *_a, **_k: fake_scvi)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "find_common_genes", lambda *_a, **_k: [])

    with pytest.raises(DataError, match="Insufficient gene overlap"):
        await ann._annotate_with_scanvi(
            adata,
            AnnotationParameters(
                method="scanvi", reference_data_id="r", cell_type_key="ctype"
            ),
            DummyWarnCtx(),
            "cell_type_scanvi",
            "confidence_scanvi",
            reference_adata=ref,
        )


@pytest.mark.asyncio
async def test_annotate_cell_types_rejects_unsupported_method_runtime_guard(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    params = AnnotationParameters(method="sctype", sctype_tissue="Brain").model_copy(
        update={"method": "unknown"}
    )

    with pytest.raises(ParameterError, match="Unsupported method"):
        await ann.annotate_cell_types("q", DummyWarnCtx({"q": adata}), params)
