"""Unit tests for annotate_cell_types routing and error contracts."""

from __future__ import annotations

from types import ModuleType, SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import AnnotationParameters
from chatspatial.tools import annotation as ann
from chatspatial.utils.exceptions import (
    DataError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


class DummyCtx:
    def __init__(self, datasets: dict[str, object]):
        self.datasets = datasets
        self.calls: list[str] = []

    async def get_adata(self, data_id: str):
        self.calls.append(data_id)
        if data_id not in self.datasets:
            raise DataNotFoundError(f"missing: {data_id}")
        return self.datasets[data_id]

    async def set_adata(self, data_id: str, adata):
        self.datasets[data_id] = adata

    async def warning(self, _msg: str):
        return None

    async def info(self, _msg: str):
        return None


def test_extract_singler_reference_labels_only_falls_back_for_missing_columns():
    class _Columns:
        def column(self, name: str):
            if name == "label.fine":
                return ["B", "T"]
            raise AttributeError(name)

    reference = SimpleNamespace(get_column_data=lambda: _Columns())

    assert ann._extract_singler_reference_labels(reference, "hpca") == ["B", "T"]

    broken_reference = SimpleNamespace(
        get_column_data=lambda: (_ for _ in ()).throw(RuntimeError("corrupt frame"))
    )
    with pytest.raises(RuntimeError, match="corrupt frame"):
        ann._extract_singler_reference_labels(broken_reference, "hpca")


def test_singler_confidence_preserves_valid_zeros_and_fills_missing_delta():
    scores = pd.DataFrame({"B": [0.8, 0.1], "T": [0.2, 0.9]})

    zero_confidence = ann._singler_confidence_values(
        ["B", "T"],
        [0.0, 0.0],
        scores,
    )
    np.testing.assert_array_equal(zero_confidence, [0.0, 0.0])
    assert np.isfinite(zero_confidence).all()

    partial_delta = ann._singler_confidence_values(
        ["B", "T"],
        [None, 0.2],
        scores,
    )
    assert partial_delta[0] == 0.8
    assert partial_delta[1] == pytest.approx(round(-np.expm1(-0.2), 3))

    adata = SimpleNamespace(obs=pd.DataFrame(index=range(2)))
    assert ann._store_confidence_if_available(
        adata,
        "confidence",
        zero_confidence,
    )
    np.testing.assert_array_equal(adata.obs["confidence"], [0.0, 0.0])
    assert ann._summarize_confidence(
        ["B", "T"],
        ["B", "T"],
        zero_confidence,
        digits=2,
    ) == {"B": 0.0, "T": 0.0}


@pytest.mark.asyncio
async def test_annotate_cell_types_sctype_happy_path_records_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx({"d1": adata})
    captured: dict[str, object] = {}

    async def _fake_sctype(_adata, _params, _ctx, output_key, confidence_key):
        _adata.obs[output_key] = ["T"] * _adata.n_obs
        return ann.AnnotationMethodOutput(
            cell_types=["T"],
            counts={"T": _adata.n_obs},
            confidence={},
            tangram_mapping_score=None,
        )

    monkeypatch.setattr(ann, "_annotate_with_sctype", _fake_sctype)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await ann.annotate_cell_types(
        "d1",
        ctx,
        AnnotationParameters(method="sctype", sctype_tissue="Brain"),
    )
    assert out.method == "sctype"
    assert out.output_key == "cell_type_sctype"
    assert out.confidence_key is None
    assert out.cell_types == ["T"]
    assert captured["analysis_name"] == "annotation_sctype"
    assert captured["results_keys"] == {
        "obs": ["cell_type_sctype"],
        "obsm": [],
        "uns": ["cell_type_sctype_counts"],
    }
    assert captured["statistics"] == {"n_cell_types": 1}
    assert ctx.datasets["d1"] is not adata
    assert "cell_type_sctype" not in adata.obs
    assert "cell_type_sctype" in ctx.datasets["d1"].obs


@pytest.mark.asyncio
async def test_annotate_empty_counts_does_not_register_uns_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """When annotation returns empty counts, no uns key is registered."""
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx({"d1": adata})
    captured: dict[str, object] = {}

    async def _fake_sctype(_adata, _params, _ctx, output_key, confidence_key):
        _adata.obs[output_key] = ["T"] * _adata.n_obs
        return ann.AnnotationMethodOutput(
            cell_types=["T"],
            counts={},  # empty counts
            confidence={},
            tangram_mapping_score=None,
        )

    monkeypatch.setattr(ann, "_annotate_with_sctype", _fake_sctype)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    await ann.annotate_cell_types(
        "d1",
        ctx,
        AnnotationParameters(method="sctype", sctype_tissue="Brain"),
    )

    assert captured["results_keys"]["uns"] == []


@pytest.mark.asyncio
async def test_annotate_cell_types_loads_reference_for_singler(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ctx = DummyCtx({"query": adata, "ref": ref})
    seen: dict[str, object] = {}

    async def _fake_singler(
        _adata, _params, _ctx, output_key, confidence_key, reference_adata
    ):
        seen["reference_adata"] = reference_adata
        _adata.obs[output_key] = ["B"] * _adata.n_obs
        return ann.AnnotationMethodOutput(
            cell_types=["B"],
            counts={"B": _adata.n_obs},
            confidence={"B": 0.88},
            tangram_mapping_score=None,
        )

    monkeypatch.setattr(ann, "_annotate_with_singler", _fake_singler)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await ann.annotate_cell_types(
        "query",
        ctx,
        AnnotationParameters(method="singler", reference_data_id="ref"),
    )
    assert seen["reference_adata"] is not ref
    assert seen["reference_adata"].obs_names.equals(ref.obs_names)
    # singler + reference_data_id="ref" → suffix "singler_ref"
    assert out.confidence_key == "confidence_singler_ref"
    assert ctx.calls == ["query", "ref"]


@pytest.mark.asyncio
async def test_annotate_cell_types_passes_through_parameter_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx({"d1": adata})

    async def _raise_param(*_args, **_kwargs):
        raise ParameterError("invalid marker_genes")

    monkeypatch.setattr(ann, "_annotate_with_cellassign", _raise_param)

    with pytest.raises(ParameterError, match="invalid marker_genes"):
        await ann.annotate_cell_types(
            "d1", ctx, AnnotationParameters(method="cellassign")
        )


@pytest.mark.asyncio
async def test_annotate_cell_types_converts_raw_import_error_to_dependency_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx({"d1": adata})

    async def _raise_import(*_args, **_kwargs):
        raise ImportError("mllmcelltype is required")

    monkeypatch.setattr(ann, "_annotate_with_mllmcelltype", _raise_import)

    with pytest.raises(
        DependencyError, match="Annotation dependency import failed"
    ) as exc:
        await ann.annotate_cell_types(
            "d1", ctx, AnnotationParameters(method="mllmcelltype")
        )
    assert isinstance(exc.value.__cause__, ImportError)


@pytest.mark.asyncio
async def test_annotate_cell_types_wraps_unexpected_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx({"d1": adata})

    async def _raise_unexpected(*_args, **_kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(ann, "_annotate_with_mllmcelltype", _raise_unexpected)

    with pytest.raises(ProcessingError, match="Annotation failed: boom"):
        await ann.annotate_cell_types(
            "d1", ctx, AnnotationParameters(method="mllmcelltype")
        )


@pytest.mark.asyncio
async def test_annotate_cell_types_late_failure_does_not_publish_partial_results(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx({"d1": adata})

    async def _fake_sctype(_adata, _params, _ctx, output_key, _confidence_key):
        _adata.obs[output_key] = ["T"] * _adata.n_obs
        return ann.AnnotationMethodOutput(
            cell_types=["T"],
            counts={"T": _adata.n_obs},
            confidence={},
        )

    monkeypatch.setattr(ann, "_annotate_with_sctype", _fake_sctype)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **_kwargs: _adata.uns.update(partial_metadata=True),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("export failed")),
    )

    with pytest.raises(RuntimeError, match="export failed"):
        await ann.annotate_cell_types(
            "d1",
            ctx,
            AnnotationParameters(method="sctype", sctype_tissue="Brain"),
        )

    assert ctx.datasets["d1"] is adata
    assert "cell_type_sctype" not in adata.obs
    assert "partial_metadata" not in adata.uns


@pytest.mark.asyncio
async def test_annotate_with_tangram_requires_reference_data_id(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(ann, "require", lambda *_args, **_kwargs: object())
    monkeypatch.setitem(__import__("sys").modules, "tangram", ModuleType("tangram"))
    ctx = DummyCtx({"d1": minimal_spatial_adata.copy()})

    with pytest.raises(ParameterError, match="Tangram requires reference_data_id"):
        await ann._annotate_with_tangram(
            minimal_spatial_adata.copy(),
            AnnotationParameters(method="tangram"),
            ctx,
            "cell_type_tangram",
            "confidence_tangram",
            reference_adata=None,
        )


@pytest.mark.asyncio
async def test_annotate_with_scanvi_requires_reference_data_id(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(
        ann, "validate_scvi_tools", lambda *_args, **_kwargs: SimpleNamespace()
    )
    ctx = DummyCtx({"d1": minimal_spatial_adata.copy()})

    with pytest.raises(ParameterError, match="scANVI requires reference_data_id"):
        await ann._annotate_with_scanvi(
            minimal_spatial_adata.copy(),
            AnnotationParameters(method="scanvi"),
            ctx,
            "cell_type_scanvi",
            "confidence_scanvi",
            reference_adata=None,
        )


@pytest.mark.asyncio
async def test_annotate_with_cellassign_requires_marker_genes(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ctx = DummyCtx({"d1": minimal_spatial_adata.copy()})
    with pytest.raises(ParameterError, match="CellAssign requires marker genes"):
        await ann._annotate_with_cellassign(
            minimal_spatial_adata.copy(),
            AnnotationParameters(method="cellassign", marker_genes=None),
            ctx,
            "cell_type_cellassign",
            "confidence_cellassign",
        )


@pytest.mark.asyncio
async def test_annotate_with_mllmcelltype_requires_cluster_label(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(ann, "require", lambda *_args, **_kwargs: object())
    monkeypatch.setitem(
        __import__("sys").modules,
        "mllmcelltype",
        ModuleType("mllmcelltype"),
    )
    ctx = DummyCtx({"d1": minimal_spatial_adata.copy()})

    with pytest.raises(ParameterError, match="cluster_label parameter is required"):
        await ann._annotate_with_mllmcelltype(
            minimal_spatial_adata.copy(),
            AnnotationParameters(method="mllmcelltype", cluster_label=None),
            ctx,
            "cell_type_mllmcelltype",
            "confidence_mllmcelltype",
        )


@pytest.mark.asyncio
async def test_annotate_with_mllmcelltype_consensus_without_models_raises_parameter_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(
        ["0"] * (adata.n_obs // 2) + ["1"] * (adata.n_obs - adata.n_obs // 2)
    )

    monkeypatch.setattr(
        ann,
        "require",
        lambda name, *_args, **_kwargs: __import__("sys").modules.get(name, object()),
    )
    monkeypatch.setitem(
        __import__("sys").modules,
        "mllmcelltype",
        ModuleType("mllmcelltype"),
    )

    def _fake_rank_genes_groups(adata_obj, *_args, **_kwargs):
        adata_obj.uns["rank_genes_groups"] = {
            "names": {
                "0": np.array(["gene_0", "gene_1"], dtype=object),
                "1": np.array(["gene_2", "gene_3"], dtype=object),
            }
        }

    monkeypatch.setattr(ann.sc.tl, "rank_genes_groups", _fake_rank_genes_groups)

    ctx = DummyCtx({"d1": adata})
    params = AnnotationParameters(
        method="mllmcelltype",
        cluster_label="leiden",
        mllm_use_consensus=True,
        mllm_models=None,
    )
    with pytest.raises(ParameterError, match="mllm_models parameter is required"):
        await ann._annotate_with_mllmcelltype(
            adata,
            params,
            ctx,
            "cell_type_mllmcelltype",
            "confidence_mllmcelltype",
        )


class DummyWarnCtx(DummyCtx):
    def __init__(self, datasets: dict[str, object]):
        super().__init__(datasets)
        self.warnings: list[str] = []

    async def warning(self, msg: str):
        self.warnings.append(msg)


@pytest.mark.asyncio
async def test_annotate_with_singler_custom_reference_success_deterministic_order(
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
                    "B": [0.8 if i % 2 == 0 else 0.2 for i in range(n_obs)],
                    "T": [0.2 if i % 2 == 0 else 0.8 for i in range(n_obs)],
                }
            )
            self._delta = [0.2] * n_obs

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

    monkeypatch.setattr(
        ann,
        "require",
        lambda name, *_args, **_kwargs: fake_singler if name == "singler" else object(),
    )
    monkeypatch.setattr(ann, "get", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=False: SimpleNamespace(
            X=_adata.X, is_integer_counts=True
        ),
    )

    ctx = DummyWarnCtx({"q": adata, "r": ref})
    out = await ann._annotate_with_singler(
        adata,
        AnnotationParameters(
            method="singler", reference_data_id="r", cell_type_key="ctype"
        ),
        ctx,
        "cell_type_singler",
        "confidence_singler",
        reference_adata=ref,
    )

    assert out.cell_types == ["B", "T"]
    assert set(out.counts.keys()) == {"B", "T"}
    assert 0.0 < out.confidence["B"] <= 1.0
    assert 0.0 < out.confidence["T"] <= 1.0
    assert adata.obs["cell_type_singler"].dtype.name == "category"
    assert "confidence_singler" in adata.obs.columns


@pytest.mark.asyncio
async def test_annotate_with_tangram_success_extracts_score_and_copies_back(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.raw = adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(
        ["B"] * (ref.n_obs // 2) + ["T"] * (ref.n_obs - ref.n_obs // 2)
    )

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_args, **_kwargs: None
    fake_tg.map_cells_to_space = lambda *_args, **_kwargs: SimpleNamespace(
        uns={"training_history": {"main_loss": ["tensor(0.9050, grad_fn=<x>)"]}}
    )

    def _project_cell_annotations(_ad_map, adata_sp, annotation):
        assert annotation == "ctype"
        adata_sp.obsm["tangram_ct_pred"] = pd.DataFrame(
            {
                "B": [0.8 if i % 2 == 0 else 0.2 for i in range(adata_sp.n_obs)],
                "T": [0.2 if i % 2 == 0 else 0.8 for i in range(adata_sp.n_obs)],
            },
            index=adata_sp.obs_names,
        )

    fake_tg.project_cell_annotations = _project_cell_annotations
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", lambda *_args, **_kwargs: fake_tg)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    ctx = DummyWarnCtx({"q": adata, "r": ref})
    out = await ann._annotate_with_tangram(
        adata,
        AnnotationParameters(
            method="tangram",
            reference_data_id="r",
            cell_type_key="ctype",
            training_genes=["gene_0", "gene_1"],
        ),
        ctx,
        "cell_type_tangram",
        "confidence_tangram",
        reference_adata=ref,
    )

    assert pytest.approx(out.tangram_mapping_score, rel=1e-6) == 0.905
    assert set(out.cell_types) == {"B", "T"}
    assert "cell_type_tangram" in adata.obs.columns
    # Default tangram_ct_pred_key when no suffix override is given
    assert "tangram_ct_pred" in adata.obsm


@pytest.mark.asyncio
async def test_annotate_with_tangram_raises_when_no_predictions_generated(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_args, **_kwargs: None
    fake_tg.map_cells_to_space = lambda *_args, **_kwargs: SimpleNamespace(
        uns={"training_history": {"main_loss": [1.0]}}
    )
    fake_tg.project_cell_annotations = lambda *_args, **_kwargs: None
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", lambda *_args, **_kwargs: fake_tg)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    with pytest.raises(ProcessingError, match="no cell type predictions"):
        await ann._annotate_with_tangram(
            adata,
            AnnotationParameters(
                method="tangram",
                reference_data_id="r",
                cell_type_key="ctype",
                training_genes=["gene_0", "gene_1"],
            ),
            DummyWarnCtx({"q": adata, "r": ref}),
            "cell_type_tangram",
            "confidence_tangram",
            reference_adata=ref,
        )


@pytest.mark.asyncio
async def test_annotate_with_tangram_treats_projection_as_core_step(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_args, **_kwargs: None
    fake_tg.map_cells_to_space = lambda *_args, **_kwargs: SimpleNamespace(
        uns={"training_history": {"main_loss": [1.0]}}
    )
    fake_tg.project_cell_annotations = lambda *_args, **_kwargs: (_ for _ in ()).throw(
        RuntimeError("projection boom")
    )

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", lambda *_args, **_kwargs: fake_tg)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    with pytest.raises(
        ProcessingError,
        match="Tangram cell annotation projection failed: projection boom",
    ):
        await ann._annotate_with_tangram(
            adata,
            AnnotationParameters(
                method="tangram",
                reference_data_id="r",
                cell_type_key="ctype",
                training_genes=["gene_0", "gene_1"],
            ),
            DummyWarnCtx({"q": adata, "r": ref}),
            "cell_type_tangram",
            "confidence_tangram",
            reference_adata=ref,
        )


@pytest.mark.asyncio
async def test_annotate_with_scanvi_direct_training_success(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(
        ["B"] * (ref.n_obs // 2) + ["T"] * (ref.n_obs - ref.n_obs // 2)
    )

    class _FakeSCANVIModel:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, adata_obj, **_kwargs):
            self._adata = adata_obj

        def train(self, **_kwargs):
            return None

        def predict(self, soft: bool = False):
            if soft:
                return pd.DataFrame(
                    {
                        "B": [
                            0.8 if i % 2 == 0 else 0.2 for i in range(self._adata.n_obs)
                        ],
                        "T": [
                            0.2 if i % 2 == 0 else 0.8 for i in range(self._adata.n_obs)
                        ],
                    },
                    index=self._adata.obs_names,
                )
            return pd.Categorical(
                ["B" if i % 2 == 0 else "T" for i in range(self._adata.n_obs)]
            )

        @staticmethod
        def load_query_data(adata_subset, _model):
            return _FakeSCANVIModel(adata_subset)

    fake_scvi = SimpleNamespace(model=SimpleNamespace(SCANVI=_FakeSCANVIModel))

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "validate_scvi_tools", lambda *_args, **_kwargs: fake_scvi)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(
        ann,
        "find_common_genes",
        lambda ref_genes, qry_genes: list(ref_genes)[
            : min(len(ref_genes), len(qry_genes))
        ],
    )
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x.copy())
    monkeypatch.setattr(ann, "ensure_counts_layer", lambda *_args, **_kwargs: None)

    ctx = DummyWarnCtx({"q": adata, "r": ref})
    out = await ann._annotate_with_scanvi(
        adata,
        AnnotationParameters(
            method="scanvi",
            reference_data_id="r",
            cell_type_key="ctype",
            scanvi_use_scvi_pretrain=False,
            scanvi_query_epochs=2,
            scanvi_scanvi_epochs=2,
        ),
        ctx,
        "cell_type_scanvi",
        "confidence_scanvi",
        reference_adata=ref,
    )

    assert set(out.cell_types) == {"B", "T"}
    assert set(out.counts.keys()) == {"B", "T"}
    assert "cell_type_scanvi" in adata.obs.columns
    assert "confidence_scanvi" in adata.obs.columns


@pytest.mark.asyncio
async def test_annotate_with_cellassign_probability_output_success(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    class _FakeCellAssign:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, adata_subset, marker_gene_matrix, **_kwargs):
            self._adata = adata_subset
            self._marker_cols = list(marker_gene_matrix.columns)

        def train(self, **_kwargs):
            return None

        def predict(self):
            n = self._adata.n_obs
            # Return probability dataframe to exercise confidence path.
            return pd.DataFrame(
                {
                    self._marker_cols[0]: [
                        0.8 if i % 2 == 0 else 0.2 for i in range(n)
                    ],
                    self._marker_cols[1]: [
                        0.2 if i % 2 == 0 else 0.8 for i in range(n)
                    ],
                },
                index=self._adata.obs_names,
            )

    fake_scvi_external = ModuleType("scvi.external")
    fake_scvi_external.CellAssign = _FakeCellAssign
    monkeypatch.setitem(__import__("sys").modules, "scvi.external", fake_scvi_external)

    monkeypatch.setattr(
        ann,
        "validate_scvi_tools",
        lambda *_args, **_kwargs: SimpleNamespace(external=fake_scvi_external),
    )
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata,
        prefer_complete_genes=True,
        require_integer_counts=True: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="raw",
        ),
    )

    params = AnnotationParameters(
        method="cellassign",
        marker_genes={"B": ["gene_0", "gene_1"], "T": ["gene_2", "gene_3"]},
    )

    out = await ann._annotate_with_cellassign(
        adata,
        params,
        DummyWarnCtx({"d1": adata}),
        "cell_type_cellassign",
        "confidence_cellassign",
    )

    assert out.cell_types == ["B", "T"]
    assert set(out.counts.keys()) == {"B", "T"}
    assert "cell_type_cellassign" in adata.obs.columns
    assert "confidence_cellassign" in adata.obs.columns


@pytest.mark.asyncio
async def test_annotate_with_cellassign_raises_when_all_markers_missing(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    monkeypatch.setattr(
        ann,
        "validate_scvi_tools",
        lambda *_args, **_kwargs: SimpleNamespace(external=fake_scvi_external),
    )
    fake_scvi_external = ModuleType("scvi.external")
    fake_scvi_external.CellAssign = object
    fake_scvi_pkg = ModuleType("scvi")
    fake_scvi_pkg.external = fake_scvi_external
    monkeypatch.setitem(__import__("sys").modules, "scvi", fake_scvi_pkg)
    monkeypatch.setitem(__import__("sys").modules, "scvi.external", fake_scvi_external)
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata,
        prefer_complete_genes=True,
        require_integer_counts=True: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="raw",
        ),
    )

    with pytest.raises(DataError, match="No valid marker genes found"):
        await ann._annotate_with_cellassign(
            adata,
            AnnotationParameters(
                method="cellassign",
                marker_genes={"B": ["NOT_A_GENE"], "T": ["ANOTHER_MISSING"]},
            ),
            DummyWarnCtx({"d1": adata}),
            "cell_type_cellassign",
            "confidence_cellassign",
        )


@pytest.mark.asyncio
async def test_annotate_with_singler_raises_when_no_reference_available(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    monkeypatch.setattr(ann, "require", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(ann, "get", lambda *_args, **_kwargs: None)
    monkeypatch.setitem(__import__("sys").modules, "singler", ModuleType("singler"))
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=False: SimpleNamespace(
            X=_adata.X, is_integer_counts=True
        ),
    )

    with pytest.raises(DataNotFoundError, match="No reference data"):
        await ann._annotate_with_singler(
            adata,
            AnnotationParameters(method="singler", reference_data_id=None),
            DummyWarnCtx({"q": adata}),
            "cell_type_singler",
            "confidence_singler",
            reference_adata=None,
        )


@pytest.mark.asyncio
async def test_annotate_with_singler_raises_on_insufficient_gene_overlap(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(["B"] * ref.n_obs)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(ann, "get", lambda *_args, **_kwargs: None)
    monkeypatch.setitem(__import__("sys").modules, "singler", ModuleType("singler"))
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "find_common_genes", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(
        ann,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=False: SimpleNamespace(
            X=_adata.X, is_integer_counts=True
        ),
    )

    with pytest.raises(DataError, match="Insufficient gene overlap"):
        await ann._annotate_with_singler(
            adata,
            AnnotationParameters(
                method="singler", reference_data_id="r", cell_type_key="ctype"
            ),
            DummyWarnCtx({"q": adata, "r": ref}),
            "cell_type_singler",
            "confidence_singler",
            reference_adata=ref,
        )


@pytest.mark.asyncio
async def test_annotate_with_scanvi_uses_ndarray_probability_fallback(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(
        ["B"] * (ref.n_obs // 2) + ["T"] * (ref.n_obs - ref.n_obs // 2)
    )

    class _FakeSCANVIModel:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, adata_obj, **_kwargs):
            self._adata = adata_obj

        def train(self, **_kwargs):
            return None

        def predict(self, soft: bool = False):
            if soft:
                return np.array(
                    [
                        [0.9, 0.1] if i % 2 == 0 else [0.1, 0.9]
                        for i in range(self._adata.n_obs)
                    ],
                    dtype=float,
                )
            return pd.Categorical(
                ["B" if i % 2 == 0 else "T" for i in range(self._adata.n_obs)]
            )

        @staticmethod
        def load_query_data(adata_subset, _model):
            return _FakeSCANVIModel(adata_subset)

    fake_scvi = SimpleNamespace(model=SimpleNamespace(SCANVI=_FakeSCANVIModel))

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "validate_scvi_tools", lambda *_args, **_kwargs: fake_scvi)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(
        ann,
        "find_common_genes",
        lambda ref_genes, qry_genes: list(ref_genes)[
            : min(len(ref_genes), len(qry_genes))
        ],
    )
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x.copy())
    monkeypatch.setattr(ann, "ensure_counts_layer", lambda *_args, **_kwargs: None)

    out = await ann._annotate_with_scanvi(
        adata,
        AnnotationParameters(
            method="scanvi",
            reference_data_id="r",
            cell_type_key="ctype",
            scanvi_use_scvi_pretrain=False,
            scanvi_query_epochs=2,
            scanvi_scanvi_epochs=2,
        ),
        DummyWarnCtx({"q": adata, "r": ref}),
        "cell_type_scanvi",
        "confidence_scanvi",
        reference_adata=ref,
    )

    assert out.confidence
    assert "confidence_scanvi" in adata.obs.columns


@pytest.mark.asyncio
async def test_annotate_with_scanvi_warns_when_probability_extraction_fails(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(
        ["B"] * (ref.n_obs // 2) + ["T"] * (ref.n_obs - ref.n_obs // 2)
    )

    class _FakeSCANVIModel:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, adata_obj, **_kwargs):
            self._adata = adata_obj

        def train(self, **_kwargs):
            return None

        def predict(self, soft: bool = False):
            if soft:
                raise RuntimeError("probabilities unavailable")
            return pd.Categorical(
                ["B" if i % 2 == 0 else "T" for i in range(self._adata.n_obs)]
            )

        @staticmethod
        def load_query_data(adata_subset, _model):
            return _FakeSCANVIModel(adata_subset)

    fake_scvi = SimpleNamespace(model=SimpleNamespace(SCANVI=_FakeSCANVIModel))

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "validate_scvi_tools", lambda *_args, **_kwargs: fake_scvi)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(
        ann,
        "find_common_genes",
        lambda ref_genes, qry_genes: list(ref_genes)[
            : min(len(ref_genes), len(qry_genes))
        ],
    )
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x.copy())
    monkeypatch.setattr(ann, "ensure_counts_layer", lambda *_args, **_kwargs: None)

    ctx = DummyWarnCtx({"q": adata, "r": ref})
    out = await ann._annotate_with_scanvi(
        adata,
        AnnotationParameters(
            method="scanvi",
            reference_data_id="r",
            cell_type_key="ctype",
            scanvi_use_scvi_pretrain=False,
            scanvi_query_epochs=2,
            scanvi_scanvi_epochs=2,
        ),
        ctx,
        "cell_type_scanvi",
        "confidence_scanvi",
        reference_adata=ref,
    )

    assert out.confidence == {}
    assert "confidence_scanvi" not in adata.obs.columns
    assert any("Could not get confidence scores" in msg for msg in ctx.warnings)


# =============================================================================
# Issue 2 regression: Tangram zero-sum rows → "unassigned" in cell_types
# =============================================================================


@pytest.mark.asyncio
async def test_annotate_with_tangram_zero_sum_rows_include_unassigned_in_cell_types(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """When some spots have zero-sum tangram_ct_pred, the returned cell_types
    and counts must include 'unassigned' to match the actual label space."""
    adata = minimal_spatial_adata.copy()
    adata.raw = adata.copy()
    ref = minimal_spatial_adata.copy()
    ref.obs["ctype"] = pd.Categorical(
        ["B"] * (ref.n_obs // 2) + ["T"] * (ref.n_obs - ref.n_obs // 2)
    )

    fake_tg = ModuleType("tangram")
    fake_tg.pp_adatas = lambda *_args, **_kwargs: None
    fake_tg.map_cells_to_space = lambda *_args, **_kwargs: SimpleNamespace(
        uns={"training_history": {"main_loss": ["tensor(0.9, grad_fn=<x>)"]}}
    )

    def _project_with_zero_rows(_ad_map, adata_sp, annotation):
        n = adata_sp.n_obs
        b_vals = [0.8 if i % 2 == 0 else 0.2 for i in range(n)]
        t_vals = [0.2 if i % 2 == 0 else 0.8 for i in range(n)]
        # Make first row all-zero → "unassigned"
        b_vals[0] = 0.0
        t_vals[0] = 0.0
        adata_sp.obsm["tangram_ct_pred"] = pd.DataFrame(
            {"B": b_vals, "T": t_vals}, index=adata_sp.obs_names
        )

    fake_tg.project_cell_annotations = _project_with_zero_rows
    monkeypatch.setitem(__import__("sys").modules, "tangram", fake_tg)

    async def _no_dupes(*_args, **_kwargs):
        return 0

    monkeypatch.setattr(ann, "require", lambda *_args, **_kwargs: fake_tg)
    monkeypatch.setattr(ann, "ensure_unique_var_names_async", _no_dupes)
    monkeypatch.setattr(ann, "get_device", lambda prefer_gpu=False: "cpu")
    monkeypatch.setattr(ann, "shallow_copy_adata", lambda x: x)

    ctx = DummyWarnCtx({"q": adata, "r": ref})
    out = await ann._annotate_with_tangram(
        adata,
        AnnotationParameters(
            method="tangram",
            reference_data_id="r",
            cell_type_key="ctype",
            training_genes=["gene_0", "gene_1"],
        ),
        ctx,
        "cell_type_tangram",
        "confidence_tangram",
        reference_adata=ref,
    )

    # "unassigned" must appear in the returned metadata
    assert "unassigned" in out.cell_types
    assert "unassigned" in out.counts
    assert out.counts["unassigned"] >= 1

    # Actual obs labels must match
    assert "unassigned" in adata.obs["cell_type_tangram"].values
    # Confidence for unassigned spots must be 0
    unassigned_mask = adata.obs["cell_type_tangram"] == "unassigned"
    assert (adata.obs.loc[unassigned_mask, "confidence_tangram"] == 0.0).all()
