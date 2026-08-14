"""Unit tests for lightweight utilities in cell_communication module."""

from __future__ import annotations

import io
import ssl
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse as sp

from chatspatial.models.data import CellCommunicationParameters
from chatspatial.tools import cell_communication as ccc
from chatspatial.tools.cell_communication import (
    CCCAutocrine,
    CCCStorage,
    _integrate_autocrine_detection,
    get_ccc_results,
    has_ccc_results,
    standardize_lr_pair,
    store_ccc_results,
)
from chatspatial.utils.dependency_manager import REnvironment
from chatspatial.utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


class DummyCtx:
    def __init__(self):
        self.warnings: list[str] = []
        self.committed = None

    async def warning(self, msg: str):
        self.warnings.append(msg)

    async def get_adata(self, _data_id: str):
        return None

    async def set_adata(self, _data_id: str, adata):
        self.committed = adata


def _required_dependency(name: str, *_args, **_kwargs):
    return __import__("sys").modules.get(name, object())


def _fail_dependency_load(*_args, **_kwargs):
    raise AssertionError("Dependency should not be loaded on this path")


def _patch_required_module(
    monkeypatch: pytest.MonkeyPatch,
    module_name: str,
    module: object,
) -> None:
    def _require_module(
        _dependency: str,
        requested_module: str,
        *_args,
        **_kwargs,
    ):
        if requested_module != module_name:
            raise AssertionError(f"Unexpected module request: {requested_module}")
        return module

    monkeypatch.setattr(ccc, "require_module", _require_module)


def _patch_required_dependency(
    monkeypatch: pytest.MonkeyPatch,
    dependency_name: str,
    dependency: object,
    *,
    configure_fastccc_database: bool = True,
) -> None:
    def _require(
        requested_dependency: str,
        *_args,
        **_kwargs,
    ):
        if requested_dependency != dependency_name:
            raise AssertionError(
                f"Unexpected dependency request: {requested_dependency}"
            )
        return dependency

    monkeypatch.setattr(ccc, "require", _require)
    if dependency_name == "fastccc" and configure_fastccc_database:
        monkeypatch.setattr(
            ccc,
            "_ensure_fastccc_database_dir",
            lambda _output_dir: "/fake/cpdb",
        )
        monkeypatch.setattr(
            ccc,
            "_build_fastccc_lri_pair_map",
            lambda _database_dir: {
                "L1^R1": "L1^R1",
                "L2-R2": "L2^R2",
            },
        )


def test_standardize_lr_pair_normalizes_separators():
    # ^ is the canonical separator — preserved as-is
    assert standardize_lr_pair("LIG^REC") == "LIG^REC"
    # Single underscore → treated as ligand/receptor separator → ^
    assert standardize_lr_pair("LIG_REC") == "LIG^REC"
    # Multi-underscore without ^ is ambiguous → returned as-is
    assert standardize_lr_pair("ITGAL_ICAM1_ITGB2") == "ITGAL_ICAM1_ITGB2"
    # Complex names with ^ are unambiguous
    assert standardize_lr_pair("ITGAL^ICAM1_ITGB2") == "ITGAL^ICAM1_ITGB2"
    # External matrix outputs may expose non-string index labels.
    assert standardize_lr_pair(123) == "123"


def test_store_and_get_ccc_results_roundtrip(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    storage = CCCStorage(
        method="liana",
        analysis_type="cluster",
        species="human",
        database="consensus",
        lr_pairs=["L1^R1"],
        top_lr_pairs=["L1^R1"],
        n_pairs=1,
        n_significant=1,
        autocrine=CCCAutocrine(n_loops=1, top_pairs=["L1^R1"]),
    )

    store_ccc_results(adata, storage)
    assert has_ccc_results(adata) is True

    restored = get_ccc_results(adata)
    assert restored is not None
    assert restored.method == "liana"
    assert restored.n_pairs == 1
    assert restored.autocrine.n_loops == 1


def test_integrate_autocrine_detection_for_liana_cluster_results():
    results = pd.DataFrame(
        {
            "source": ["T", "B", "T", "T", "B"],
            "target": ["T", "B", "T", "B", "B"],
            "ligand_complex": ["L1", "L3", "L1", "L2", "L3"],
            "receptor_complex": ["R1", "R3", "R1", "R2", "R3"],
            "magnitude_rank": [0.03, 0.02, 0.01, 0.2, 0.04],
        }
    )
    storage = CCCStorage(
        method="liana",
        analysis_type="cluster",
        species="human",
        database="consensus",
        results=results,
    )

    _integrate_autocrine_detection(storage, n_top=5)

    assert storage.autocrine.n_loops == 2
    assert storage.autocrine.top_pairs == ["L1^R1", "L3^R3"]
    assert len(storage.autocrine.results) == 4


def test_integrate_autocrine_detection_for_matrix_based_methods():
    # Simulate cellphonedb/fastccc matrix format: columns are source|target
    results = pd.DataFrame(
        {
            "interacting_pair": ["L1^R1", "L1_R1", "L2^R2"],
            "T|T": [0.4, 0.2, 0.0],
            "T|B": [0.1, 0.2, 0.2],
            "B|B": [0.0, 0.0, 0.9],
        },
        index=[10, 11, 12],
    )
    storage = CCCStorage(
        method="cellphonedb",
        analysis_type="cluster",
        species="human",
        database="cellphonedb",
        results=results,
    )

    _integrate_autocrine_detection(storage, n_top=5)

    assert storage.autocrine.n_loops == 2
    # Ranked by autocrine strength: L2^R2 scores 0.9 in B|B, L1^R1 0.4 in T|T.
    assert storage.autocrine.top_pairs == ["L2^R2", "L1^R1"]
    assert len(storage.autocrine.results) == 3


def test_integrate_autocrine_detection_counts_only_significant_loops():
    # Interaction strength is non-zero for nearly every catalogued pair, so a
    # "strength > 0" filter would call all three rows autocrine loops.
    results = pd.DataFrame(
        {
            "interacting_pair": ["L1^R1", "L2^R2", "L3^R3"],
            "T|T": [0.4, 0.6, 0.2],
            "T|B": [0.1, 0.2, 0.2],
            "B|B": [0.3, 0.1, 0.5],
        },
        index=[10, 11, 12],
    )
    pvalues = pd.DataFrame(
        {
            "interacting_pair": ["L1^R1", "L2^R2", "L3^R3"],
            "T|T": [0.9, 0.01, 0.8],
            "T|B": [0.001, 0.9, 0.9],
            "B|B": [0.02, 0.9, 0.7],
        },
        index=[10, 11, 12],
    )
    storage = CCCStorage(
        method="fastccc",
        analysis_type="cluster",
        species="human",
        database="CellPhoneDB_v5",
        results=results,
        pvalues=pvalues,
        statistics={"pvalue_threshold": 0.05},
    )

    _integrate_autocrine_detection(storage, n_top=5)

    # L1^R1 is significant in B|B, L2^R2 in T|T; L3^R3 is significant nowhere,
    # and L1^R1's T|B hit is paracrine, not a loop.
    assert storage.autocrine.n_loops == 2
    assert storage.autocrine.top_pairs == ["L2^R2", "L1^R1"]


def test_integrate_autocrine_detection_respects_backend_significance_call():
    """Loops must not outnumber the significant pairs they are drawn from.

    The per-column threshold is uncorrected, while the backend applies
    multiple-testing correction across cell-type pairs. Reusing the backend's
    own call keeps the two numbers consistent.
    """
    results = pd.DataFrame(
        {
            "interacting_pair": ["L1^R1", "L2^R2"],
            "T|T": [0.4, 0.6],
            "B|B": [0.3, 0.1],
        },
        index=[10, 11],
    )
    pvalues = pd.DataFrame(
        {
            "interacting_pair": ["L1^R1", "L2^R2"],
            "T|T": [0.04, 0.01],
            "B|B": [0.9, 0.9],
        },
        index=[10, 11],
    )
    storage = CCCStorage(
        method="fastccc",
        analysis_type="cluster",
        species="human",
        database="CellPhoneDB_v5",
        results=results,
        pvalues=pvalues,
        n_significant=1,
        statistics={"pvalue_threshold": 0.05},
        # L1^R1 clears the raw threshold but not the correction.
        method_data={"significant_lr_mask": np.array([False, True])},
    )

    _integrate_autocrine_detection(storage, n_top=5)

    assert storage.autocrine.n_loops == 1
    assert storage.autocrine.top_pairs == ["L2^R2"]


def test_integrate_autocrine_detection_filters_liana_by_significance():
    results = pd.DataFrame(
        {
            "source": ["T", "T", "B"],
            "target": ["T", "B", "B"],
            "ligand_complex": ["L1", "L2", "L3"],
            "receptor_complex": ["R1", "R2", "R3"],
            "magnitude_rank": [0.01, 0.02, 0.4],
        }
    )
    storage = CCCStorage(
        method="liana",
        analysis_type="cluster",
        species="human",
        database="consensus",
        results=results,
        statistics={"significance_threshold": 0.05},
    )

    _integrate_autocrine_detection(storage, n_top=5)

    # B->B exists but ranks far below the threshold, so it is not a loop.
    assert storage.autocrine.n_loops == 1
    assert storage.autocrine.top_pairs == ["L1^R1"]


@pytest.mark.asyncio
async def test_validate_ccc_params_requires_spatial_connectivity_for_liana(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    params = CellCommunicationParameters(
        method="liana",
        species="human",
        cell_type_key="cell_type",
        perform_spatial_analysis=True,
    )
    with pytest.raises(DataNotFoundError, match="Spatial connectivity required"):
        await ccc._validate_ccc_params(adata, params, DummyCtx())


@pytest.mark.asyncio
async def test_validate_ccc_params_warns_for_mouse_consensus_resource(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    ctx = DummyCtx()
    params = CellCommunicationParameters(
        method="liana",
        species="mouse",
        cell_type_key="cell_type",
        perform_spatial_analysis=False,
        liana_resource="consensus",
    )
    await ccc._validate_ccc_params(adata, params, ctx)
    assert any("mouseconsensus" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_run_ccc_analysis_dispatches_to_fastccc(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    params = CellCommunicationParameters(
        method="fastccc",
        species="human",
        cell_type_key="cell_type",
    )
    expected = CCCStorage(
        method="fastccc",
        analysis_type="cluster",
        species="human",
        database="fastccc",
        n_pairs=2,
        n_significant=1,
    )

    async def _fake_fastccc(*_args, **_kwargs):
        return expected

    monkeypatch.setattr(ccc, "require", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(ccc, "_analyze_communication_fastccc", _fake_fastccc)

    out = await ccc._run_ccc_analysis(adata, params, DummyCtx())
    assert out is expected


@pytest.mark.asyncio
async def test_analyze_cell_communication_happy_path_cluster(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    params = CellCommunicationParameters(
        method="fastccc",
        species="human",
        cell_type_key="cell_type",
    )
    ctx = DummyCtx()
    ctx.get_adata = lambda _data_id: adata  # type: ignore[method-assign]
    captured: dict[str, object] = {}

    async def _ok_get_adata(_data_id: str):
        return adata

    ctx.get_adata = _ok_get_adata  # type: ignore[method-assign]

    async def _fake_validate(*_args, **_kwargs):
        return None

    async def _fake_run(*_args, **_kwargs):
        return CCCStorage(
            method="fastccc",
            analysis_type="cluster",
            species="human",
            database="fastccc",
            lr_pairs=["L1^R1"],
            top_lr_pairs=["L1^R1"],
            n_pairs=1,
            n_significant=1,
            statistics={"ok": 1.0},
        )

    monkeypatch.setattr(ccc, "_validate_ccc_params", _fake_validate)
    monkeypatch.setattr(ccc, "_run_ccc_analysis", _fake_run)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await ccc.analyze_cell_communication("d1", ctx, params)
    assert out.method == "fastccc"
    assert out.n_lr_pairs == 1
    assert out.results_key == "ccc_fastccc"
    assert captured["analysis_name"] == "cell_communication_fastccc"
    assert captured["results_keys"] == {"obs": [], "obsm": [], "uns": ["ccc_fastccc"]}
    assert ctx.committed is not adata
    assert "ccc_fastccc" not in adata.uns
    assert "ccc_fastccc" in ctx.committed.uns


@pytest.mark.asyncio
async def test_analyze_cell_communication_wraps_unexpected_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    params = CellCommunicationParameters(
        method="fastccc",
        species="human",
        cell_type_key="cell_type",
    )

    async def _get_adata(_data_id: str):
        return adata

    ctx = DummyCtx()
    ctx.get_adata = _get_adata  # type: ignore[method-assign]

    async def _raise_validate(*_args, **_kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(ccc, "_validate_ccc_params", _raise_validate)

    with pytest.raises(
        ProcessingError, match="Error in cell communication analysis: boom"
    ):
        await ccc.analyze_cell_communication("d1", ctx, params)


@pytest.mark.asyncio
async def test_analyze_cell_communication_passes_through_data_not_found(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    params = CellCommunicationParameters(
        method="fastccc",
        species="human",
        cell_type_key="cell_type",
    )

    async def _get_adata(_data_id: str):
        return adata

    async def _raise_not_found(*_args, **_kwargs):
        raise DataNotFoundError("missing neighbor graph")

    ctx = DummyCtx()
    ctx.get_adata = _get_adata  # type: ignore[method-assign]
    monkeypatch.setattr(ccc, "_validate_ccc_params", _raise_not_found)

    with pytest.raises(DataNotFoundError, match="missing neighbor graph"):
        await ccc.analyze_cell_communication("d1", ctx, params)


def test_get_liana_resource_name_maps_mouse_consensus_and_passthrough():
    assert ccc._get_liana_resource_name("mouse", "consensus") == "mouseconsensus"
    assert ccc._get_liana_resource_name("mouse", "cellphonedb") == "cellphonedb"
    assert ccc._get_liana_resource_name("human", "consensus") == "consensus"


@pytest.mark.asyncio
async def test_validate_ccc_params_cellphonedb_warns_low_gene_and_cell_counts(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    ctx = DummyCtx()

    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True: type(
            "Raw", (), {"var_names": list(range(300))}
        )(),
    )

    await ccc._validate_ccc_params(
        adata,
        CellCommunicationParameters(
            method="cellphonedb",
            species="human",
            cell_type_key="cell_type",
        ),
        ctx,
    )

    assert any("Gene count" in w for w in ctx.warnings)
    assert any("Cell count" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_run_ccc_analysis_unknown_method_raises_parameter_error(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    params = CellCommunicationParameters(
        method="fastccc", species="human", cell_type_key="cell_type"
    ).model_copy(update={"method": "unknown"})

    with pytest.raises(ParameterError, match="Unsupported method"):
        await ccc._run_ccc_analysis(adata, params, DummyCtx())


@pytest.mark.asyncio
async def test_analyze_cell_communication_spatial_writes_obsm_scores(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    params = CellCommunicationParameters(
        method="liana", species="human", cell_type_key="cell_type"
    )

    async def _get_adata(_data_id: str):
        return adata

    ctx = DummyCtx()
    ctx.get_adata = _get_adata  # type: ignore[method-assign]
    captured: dict[str, object] = {}

    async def _fake_validate(*_args, **_kwargs):
        return None

    async def _fake_run(*_args, **_kwargs):
        return CCCStorage(
            method="liana",
            analysis_type="spatial",
            species="human",
            database="consensus",
            lr_pairs=["L1^R1"],
            top_lr_pairs=["L1^R1"],
            n_pairs=1,
            n_significant=1,
            method_data={
                "spatial_scores": pd.DataFrame(
                    {"x": [0.1] * adata.n_obs}, index=adata.obs_names
                ),
                "spatial_pvals": pd.DataFrame(
                    {"x": [0.05] * adata.n_obs}, index=adata.obs_names
                ),
            },
        )

    monkeypatch.setattr(ccc, "_validate_ccc_params", _fake_validate)
    monkeypatch.setattr(ccc, "_run_ccc_analysis", _fake_run)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await ccc.analyze_cell_communication("d1", ctx, params)

    assert out.analysis_type == "spatial"
    committed = ctx.committed
    assert committed is not adata
    # Shared keys still written for viz backward compatibility
    assert ccc.CCC_SPATIAL_SCORES_KEY in committed.obsm
    assert ccc.CCC_SPATIAL_PVALS_KEY in committed.obsm
    # Per-method keys also written
    assert "ccc_spatial_scores_liana" in committed.obsm
    assert "ccc_spatial_pvals_liana" in committed.obsm
    # Metadata points to per-method keys (not shared)
    assert "ccc_spatial_scores_liana" in captured["results_keys"]["obsm"]
    assert "ccc_spatial_pvals_liana" in captured["results_keys"]["obsm"]
    assert ccc.CCC_SPATIAL_SCORES_KEY not in captured["results_keys"]["obsm"]
    assert ccc.CCC_SPATIAL_PVALS_KEY not in captured["results_keys"]["obsm"]


@pytest.mark.asyncio
async def test_analyze_cell_communication_late_failure_keeps_source_unchanged(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    ctx = DummyCtx()

    async def _get_adata(_data_id: str):
        return adata

    async def _validate(*_args, **_kwargs):
        return None

    async def _run(candidate, *_args, **_kwargs):
        candidate.uns["backend_partial"] = True
        return CCCStorage(
            method="fastccc",
            analysis_type="cluster",
            species="human",
            database="fastccc",
            n_pairs=1,
            n_significant=1,
        )

    ctx.get_adata = _get_adata  # type: ignore[method-assign]
    monkeypatch.setattr(ccc, "_validate_ccc_params", _validate)
    monkeypatch.setattr(ccc, "_run_ccc_analysis", _run)
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("export failed")),
    )

    params = CellCommunicationParameters(
        method="fastccc", species="human", cell_type_key="cell_type"
    )
    with pytest.raises(ProcessingError, match="export failed"):
        await ccc.analyze_cell_communication("d1", ctx, params)

    assert ctx.committed is None
    assert "backend_partial" not in adata.uns
    assert "ccc_fastccc" not in adata.uns


@pytest.mark.asyncio
async def test_analyze_cell_communication_passes_through_data_compatibility_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    params = CellCommunicationParameters(
        method="fastccc", species="human", cell_type_key="cell_type"
    )

    async def _get_adata(_data_id: str):
        return adata

    async def _raise_compat(*_args, **_kwargs):
        raise DataCompatibilityError("incompatible matrix")

    ctx = DummyCtx()
    ctx.get_adata = _get_adata  # type: ignore[method-assign]
    monkeypatch.setattr(ccc, "_validate_ccc_params", _raise_compat)

    with pytest.raises(DataCompatibilityError, match="incompatible matrix"):
        await ccc.analyze_cell_communication("d1", ctx, params)


def test_run_liana_cluster_analysis_builds_expected_storage(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(
        ["T"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )

    def _rank_aggregate(_adata, **_kwargs):
        _adata.uns["liana_res"] = pd.DataFrame(
            {
                "source": ["T", "B", "T", "B"],
                "target": ["B", "T", "T", "B"],
                "ligand_complex": ["L1", "L1", "L2", "L3"],
                "receptor_complex": ["R1", "R1", "R2", "R3"],
                "magnitude_rank": [0.01, 0.02, 0.03, 0.2],
            }
        )

    fake_liana = type(
        "L",
        (),
        {"mt": type("MT", (), {"rank_aggregate": staticmethod(_rank_aggregate)})},
    )()
    monkeypatch.setitem(__import__("sys").modules, "liana", fake_liana)
    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True: type("Raw", (), {"source": "raw"})(),
    )

    params = CellCommunicationParameters(
        method="liana",
        species="human",
        cell_type_key="cell_type",
        perform_spatial_analysis=False,
        plot_top_pairs=3,
        liana_significance_alpha=0.05,
    )

    out = ccc._run_liana_cluster_analysis(adata, params, DummyCtx(), fake_liana)

    assert out.analysis_type == "cluster"
    assert out.n_pairs == 4
    assert out.n_significant == 3
    assert out.lr_pairs == ["L1^R1", "L1^R1", "L2^R2", "L3^R3"]
    # L1^R1 is scored twice (T->B and B->T); deduplicating before truncation is
    # what lets plot_top_pairs=3 actually return three distinct pairs.
    assert out.top_lr_pairs == ["L1^R1", "L2^R2", "L3^R3"]
    assert len(out.results) == 4
    assert out.statistics["use_raw"] is True


def test_run_liana_spatial_analysis_builds_expected_storage(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    var = pd.DataFrame(
        {
            "morans": [0.9, 0.8, 0.3],
            "morans_pvals": [0.001, 0.002, 0.3],
        },
        index=["L1^R1", "L1^R1", "L2-R2"],
    )
    lrdata = type(
        "LR",
        (),
        {
            "n_vars": 3,
            "var": var,
            "X": np.array([[0.5, 0.4, 0.2]] * adata.n_obs),
            "layers": {"pvals": np.array([[0.01, 0.02, 0.2]] * adata.n_obs)},
        },
    )()

    fake_liana = type(
        "L",
        (),
        {
            "mt": type(
                "MT", (), {"bivariate": staticmethod(lambda *_args, **_kwargs: lrdata)}
            )
        },
    )()
    monkeypatch.setitem(__import__("sys").modules, "liana", fake_liana)

    params = CellCommunicationParameters(
        method="liana",
        species="human",
        liana_global_metric="morans",
        cell_type_key="cell_type",
        plot_top_pairs=3,
    )

    out = ccc._run_liana_spatial_analysis(adata, params, DummyCtx(), fake_liana)

    assert out.analysis_type == "spatial"
    assert out.n_pairs == 3
    assert out.lr_pairs == ["L1^R1", "L1^R1", "L2^R2"]
    assert out.top_lr_pairs == ["L1^R1", "L2^R2"]
    assert len(out.results) == 3
    assert out.method_data["spatial_scores"].shape == (adata.n_obs, 3)
    assert out.method_data["spatial_pvals"].shape == (adata.n_obs, 3)
    assert "morans_pvals_corrected" in out.results.columns
    assert "morans_significant" in out.results.columns


@pytest.mark.asyncio
async def test_analyze_communication_liana_auto_builds_spatial_neighbors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    calls: dict[str, object] = {}

    def _spatial_neighbors(_adata, **kwargs):
        calls["neighbors_kwargs"] = kwargs
        _adata.obsp["spatial_connectivities"] = np.eye(_adata.n_obs)

    fake_sq = type(
        "SQ",
        (),
        {"gr": type("GR", (), {"spatial_neighbors": staticmethod(_spatial_neighbors)})},
    )()
    monkeypatch.setitem(__import__("sys").modules, "squidpy", fake_sq)
    monkeypatch.setitem(__import__("sys").modules, "liana", type("L", (), {})())
    monkeypatch.setattr(
        ccc,
        "require",
        lambda name, *_args, **_kwargs: fake_sq if name == "squidpy" else object(),
    )

    expected = CCCStorage(
        method="liana",
        analysis_type="spatial",
        species="human",
        database="consensus",
        n_pairs=1,
        n_significant=1,
    )

    monkeypatch.setattr(
        ccc, "_run_liana_spatial_analysis", lambda *_args, **_kwargs: expected
    )

    params = CellCommunicationParameters(
        method="liana",
        species="human",
        cell_type_key="cell_type",
        perform_spatial_analysis=True,
    )

    out = await ccc._analyze_communication_liana(adata, params, DummyCtx())

    assert out is expected
    assert "neighbors_kwargs" in calls
    assert calls["neighbors_kwargs"]["coord_type"] == "generic"


@pytest.mark.asyncio
async def test_analyze_communication_liana_without_cluster_column_dispatches_spatial(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    if "cell_type" in adata.obs:
        del adata.obs["cell_type"]

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setitem(__import__("sys").modules, "liana", type("L", (), {})())

    expected = CCCStorage(
        method="liana",
        analysis_type="spatial",
        species="human",
        database="consensus",
        n_pairs=0,
        n_significant=0,
    )

    monkeypatch.setattr(
        ccc, "_run_liana_spatial_analysis", lambda *_args, **_kwargs: expected
    )

    params = CellCommunicationParameters(
        method="liana",
        species="human",
        cell_type_key="cell_type",
        perform_spatial_analysis=False,
    )

    out = await ccc._analyze_communication_liana(adata, params, DummyCtx())
    assert out is expected


@pytest.mark.asyncio
async def test_analyze_communication_liana_cluster_dispatch_when_cluster_present(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setitem(__import__("sys").modules, "liana", type("L", (), {})())

    expected = CCCStorage(
        method="liana",
        analysis_type="cluster",
        species="human",
        database="consensus",
        n_pairs=1,
        n_significant=1,
    )

    monkeypatch.setattr(
        ccc, "_run_liana_cluster_analysis", lambda *_args, **_kwargs: expected
    )

    params = CellCommunicationParameters(
        method="liana",
        species="human",
        cell_type_key="cell_type",
        perform_spatial_analysis=False,
    )

    out = await ccc._analyze_communication_liana(adata, params, DummyCtx())
    assert out is expected


@pytest.mark.asyncio
async def test_analyze_communication_cellphonedb_rejects_non_human_early(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    monkeypatch.setattr(ccc, "require_module", _fail_dependency_load)

    with pytest.raises(ParameterError, match="CellPhoneDB only supports human"):
        await ccc._analyze_communication_cellphonedb(
            adata,
            CellCommunicationParameters(
                method="cellphonedb",
                species="mouse",
                cell_type_key="cell_type",
            ),
            DummyCtx(),
        )


@pytest.mark.asyncio
async def test_analyze_communication_fastccc_rejects_non_human_early(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    with pytest.raises(ParameterError, match="FastCCC only supports human"):
        await ccc._analyze_communication_fastccc(
            adata,
            CellCommunicationParameters(
                method="fastccc",
                species="mouse",
                cell_type_key="cell_type",
            ),
            DummyCtx(),
        )


@pytest.mark.asyncio
async def test_run_ccc_analysis_dispatches_to_liana(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    expected = CCCStorage(
        method="liana",
        analysis_type="cluster",
        species="human",
        database="consensus",
        n_pairs=1,
        n_significant=1,
    )

    async def _fake_liana(*_args, **_kwargs):
        return expected

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(ccc, "_analyze_communication_liana", _fake_liana)

    out = await ccc._run_ccc_analysis(
        adata,
        CellCommunicationParameters(
            method="liana", species="human", cell_type_key="cell_type"
        ),
        DummyCtx(),
    )
    assert out is expected


@pytest.mark.asyncio
async def test_run_ccc_analysis_dispatches_to_cellphonedb_and_cellchat_r(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    cpdb_expected = CCCStorage(
        method="cellphonedb",
        analysis_type="cluster",
        species="human",
        database="cellphonedb",
        n_pairs=2,
        n_significant=1,
    )
    chat_expected = CCCStorage(
        method="cellchat_r",
        analysis_type="cluster",
        species="human",
        database="cellchatdb",
        n_pairs=2,
        n_significant=1,
    )

    async def _fake_cpdb(*_args, **_kwargs):
        return cpdb_expected

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(ccc, "_analyze_communication_cellphonedb", _fake_cpdb)
    monkeypatch.setattr(
        ccc,
        "_analyze_communication_cellchat_r",
        lambda *_args, **_kwargs: chat_expected,
    )

    cpdb_out = await ccc._run_ccc_analysis(
        adata,
        CellCommunicationParameters(
            method="cellphonedb", species="human", cell_type_key="cell_type"
        ),
        DummyCtx(),
    )
    chat_out = await ccc._run_ccc_analysis(
        adata,
        CellCommunicationParameters(
            method="cellchat_r", species="human", cell_type_key="cell_type"
        ),
        DummyCtx(),
    )

    assert cpdb_out is cpdb_expected
    assert chat_out is chat_expected


@pytest.mark.asyncio
async def test_create_microenvironments_file_returns_none_without_spatial(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    if "spatial" in adata.obsm:
        del adata.obsm["spatial"]
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    out = await ccc._create_microenvironments_file(
        adata,
        CellCommunicationParameters(
            method="cellphonedb",
            species="human",
            cell_type_key="cell_type",
            cellphonedb_spatial_radius=10.0,
        ),
        DummyCtx(),
    )
    assert out is None


@pytest.mark.asyncio
async def test_create_microenvironments_file_writes_expected_format(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(
        ["T"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )
    # Force all cells into one connected neighborhood signature.
    adata.obsm["spatial"] = np.column_stack(
        [
            np.linspace(0.0, 1.0, adata.n_obs),
            np.linspace(0.0, 1.0, adata.n_obs),
        ]
    )

    path = await ccc._create_microenvironments_file(
        adata,
        CellCommunicationParameters(
            method="cellphonedb",
            species="human",
            cell_type_key="cell_type",
            cellphonedb_spatial_radius=10.0,
        ),
        DummyCtx(),
    )

    assert path is not None
    with open(path, encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    assert lines[0] == "cell_type\tmicroenvironment"
    data_lines = set(lines[1:])
    # Both cell types should map into at least one shared microenvironment.
    assert any(ln.startswith("T\tmicroenv_") for ln in data_lines)
    assert any(ln.startswith("B\tmicroenv_") for ln in data_lines)

    import os

    os.remove(path)


@pytest.mark.asyncio
async def test_create_microenvironments_file_skips_empty_assignments(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    ctx = DummyCtx()

    path = await ccc._create_microenvironments_file(
        adata,
        CellCommunicationParameters(
            method="cellphonedb",
            species="human",
            cell_type_key="cell_type",
            cellphonedb_spatial_radius=10.0,
        ),
        ctx,
    )

    assert path is None
    assert any(
        "No mixed cell-type neighborhoods" in warning for warning in ctx.warnings
    )


def test_integrate_autocrine_detection_cellchat_prob_matrix_branch():
    results = pd.DataFrame({"interaction_name": ["L1^R1", "L1_R1", "L2-R2", "L3^R3"]})
    prob = np.zeros((2, 2, 4), dtype=float)
    prob[0, 0, :] = [1.0, 0.5, 0.0, 2.0]
    prob[1, 1, :] = [0.0, 0.0, 0.0, 0.0]

    storage = CCCStorage(
        method="cellchat_r",
        analysis_type="cluster",
        species="human",
        database="cellchatdb",
        results=results,
        lr_pairs=["L1^R1", "L1^R1", "L2^R2", "L3^R3"],
        method_data={"prob_matrix": prob},
    )

    ccc._integrate_autocrine_detection(storage, n_top=2)

    assert storage.autocrine.n_loops == 2
    assert storage.autocrine.top_pairs == ["L1^R1", "L3^R3"]


def _install_fake_cellphonedb_download(
    monkeypatch: pytest.MonkeyPatch, download_impl
) -> None:
    monkeypatch.setattr(ccc, "_download_cellphonedb_release", download_impl)


def _write_valid_cellphonedb_archive(path: Path) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        for name in ccc.CELLPHONEDB_REQUIRED_TABLES:
            archive.writestr(name, "id\n")


def _valid_cellphonedb_archive_bytes() -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("gene_table.csv", "id\n")
    return buffer.getvalue()


def test_download_cellphonedb_release_extracts_database_with_scoped_tls(
    tmp_path, monkeypatch: pytest.MonkeyPatch
):
    release_buffer = io.BytesIO()
    with zipfile.ZipFile(release_buffer, "w") as release:
        release.writestr(
            "cellphonedb-data-v5.0.0/cellphonedb.zip",
            _valid_cellphonedb_archive_bytes(),
        )

    ssl_context = ssl.create_default_context()
    request: dict[str, object] = {}

    def _urlopen(url: str, *, context: ssl.SSLContext):
        request.update(url=url, context=context)
        return io.BytesIO(release_buffer.getvalue())

    monkeypatch.setattr(ccc.urllib.request, "urlopen", _urlopen)

    ccc._download_cellphonedb_release(tmp_path, "v5.0.0", ssl_context)

    assert request == {
        "url": (
            "https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v5.0.0.zip"
        ),
        "context": ssl_context,
    }
    assert zipfile.is_zipfile(tmp_path / "cellphonedb.zip")


@pytest.mark.parametrize(
    "members",
    [
        {},
        {
            "first/cellphonedb.zip": b"first",
            "second/cellphonedb.zip": b"second",
        },
    ],
    ids=["missing", "ambiguous"],
)
def test_download_cellphonedb_release_rejects_invalid_database_members(
    tmp_path, monkeypatch: pytest.MonkeyPatch, members: dict[str, bytes]
):
    release_buffer = io.BytesIO()
    with zipfile.ZipFile(release_buffer, "w") as release:
        for name, content in members.items():
            release.writestr(name, content)

    monkeypatch.setattr(
        ccc.urllib.request,
        "urlopen",
        lambda *_args, **_kwargs: io.BytesIO(release_buffer.getvalue()),
    )

    with pytest.raises(
        DependencyError,
        match="must contain exactly one cellphonedb.zip archive",
    ):
        ccc._download_cellphonedb_release(
            tmp_path,
            "v5.0.0",
            ssl.create_default_context(),
        )


def test_ensure_cellphonedb_database_returns_existing_file(
    tmp_path, monkeypatch: pytest.MonkeyPatch
):
    import os

    monkeypatch.setattr(ccc, "require_module", _fail_dependency_load)

    db_file = tmp_path / "cellphonedb.zip"
    _write_valid_cellphonedb_archive(db_file)

    out = ccc._ensure_cellphonedb_database(str(tmp_path))

    assert out == os.path.join(str(tmp_path), "cellphonedb.zip")


def test_ensure_cellphonedb_database_downloads_atomically_with_scoped_tls(
    tmp_path, monkeypatch: pytest.MonkeyPatch
):
    import os

    calls: dict[str, object] = {}
    scoped_context = object()

    def _download(output_dir, version, ssl_context):
        calls["download_args"] = (output_dir, version, ssl_context)
        _write_valid_cellphonedb_archive(output_dir / "cellphonedb.zip")

    _install_fake_cellphonedb_download(monkeypatch, _download)
    monkeypatch.setattr("certifi.where", lambda: str(tmp_path / "ca.pem"))
    monkeypatch.setattr(
        ccc.ssl,
        "create_default_context",
        lambda *, cafile: calls.update(cafile=cafile) or scoped_context,
    )

    original_https_context = ssl._create_default_https_context

    out = ccc._ensure_cellphonedb_database(str(tmp_path))

    assert out == os.path.join(str(tmp_path), "cellphonedb.zip")
    download_dir, version, passed_context = calls["download_args"]
    assert Path(download_dir).parent == tmp_path
    assert Path(download_dir) != tmp_path
    assert version == ccc.CELLPHONEDB_DATABASE_VERSION
    assert passed_context is scoped_context
    assert calls["cafile"] == str(tmp_path / "ca.pem")
    assert zipfile.is_zipfile(out)
    assert ssl._create_default_https_context is original_https_context


def test_ensure_cellphonedb_database_replaces_corrupt_cache(
    tmp_path, monkeypatch: pytest.MonkeyPatch
):
    db_file = tmp_path / "cellphonedb.zip"
    db_file.write_bytes(b"partial download")
    calls = {"count": 0}

    def _download(output_dir, _version, _ssl_context):
        calls["count"] += 1
        _write_valid_cellphonedb_archive(output_dir / "cellphonedb.zip")

    _install_fake_cellphonedb_download(monkeypatch, _download)

    out = ccc._ensure_cellphonedb_database(str(tmp_path))

    assert out == str(db_file)
    assert calls["count"] == 1
    assert zipfile.is_zipfile(db_file)


def test_ensure_fastccc_database_extracts_shared_archive(
    tmp_path, monkeypatch: pytest.MonkeyPatch
):
    archive_path = tmp_path / "cellphonedb.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        for name in ccc.CELLPHONEDB_REQUIRED_TABLES:
            archive.writestr(name, f"content-{name}")

    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda _output_dir: str(archive_path),
    )

    database_dir = Path(ccc._ensure_fastccc_database_dir(str(tmp_path)))

    assert database_dir == tmp_path / "fastccc"
    for name in ccc.FASTCCC_REQUIRED_TABLES:
        assert (database_dir / name).read_text() == f"content-{name}"


def test_ensure_fastccc_database_reuses_complete_extraction(
    tmp_path, monkeypatch: pytest.MonkeyPatch
):
    database_dir = tmp_path / "fastccc"
    database_dir.mkdir()
    for name in ccc.FASTCCC_REQUIRED_TABLES:
        (database_dir / name).write_text("cached")

    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda _output_dir: (_ for _ in ()).throw(
            AssertionError("complete extraction must not re-open the archive")
        ),
    )

    assert ccc._ensure_fastccc_database_dir(str(tmp_path)) == str(database_dir)


@pytest.mark.asyncio
async def test_analyze_communication_fastccc_success_single_method(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(
        ["T"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )

    interactions_strength = pd.DataFrame(
        {
            "T|T": [0.8, 0.1],
            "T|B": [0.2, 0.9],
        },
        index=["L1^R1", "L2-R2"],
    )
    pvalues = pd.DataFrame(
        {
            "T|T": [0.01, 0.2],
            "T|B": [0.03, 0.8],
        },
        index=["L1^R1", "L2-R2"],
    )
    percentages = pd.DataFrame({"pct": [0.5, 0.6]}, index=["L1^R1", "L2-R2"])

    def _fake_statistical_analysis_method(**kwargs):
        saved = __import__("anndata").read_h5ad(kwargs["counts_file_path"])
        assert sp.issparse(saved.X)
        return interactions_strength, pvalues, percentages

    from types import ModuleType

    fastccc_mod = ModuleType("fastccc")
    fastccc_mod.statistical_analysis_method = _fake_statistical_analysis_method
    _patch_required_dependency(monkeypatch, "fastccc", fastccc_mod)

    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True: type(
            "Raw",
            (),
            {
                "X": sp.csr_matrix(np.asarray(_adata.X)),
                "var_names": _adata.var_names,
            },
        )(),
    )

    params = CellCommunicationParameters(
        method="fastccc",
        species="human",
        cell_type_key="cell_type",
        fastccc_use_cauchy=False,
        fastccc_pvalue_threshold=0.05,
        plot_top_pairs=1,
    )

    out = await ccc._analyze_communication_fastccc(adata, params, DummyCtx())

    assert out.method == "fastccc"
    assert out.n_pairs == 2
    assert out.n_significant == 1
    assert out.lr_pairs == ["L1^R1", "L2^R2"]
    assert out.top_lr_pairs == ["L2^R2"]
    assert out.method_data["percentages"].index.tolist() == ["L1^R1", "L2^R2"]


def test_ensure_cellphonedb_database_raises_dependency_error_on_download_failure(
    tmp_path, monkeypatch: pytest.MonkeyPatch
):
    def _download_fail(*_args, **_kwargs):
        raise RuntimeError("network down")

    _install_fake_cellphonedb_download(monkeypatch, _download_fail)
    monkeypatch.setattr("certifi.where", lambda: str(tmp_path / "ca.pem"))
    monkeypatch.setattr(
        ccc.ssl,
        "create_default_context",
        lambda *, cafile: object(),
    )

    with pytest.raises(
        DependencyError, match="Failed to download CellPhoneDB database"
    ):
        ccc._ensure_cellphonedb_database(str(tmp_path))


@pytest.mark.asyncio
async def test_create_microenvironments_file_warns_and_returns_none_on_exception(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    class _BrokenNN:
        def __init__(self, *args, **kwargs):
            pass

        def fit(self, *_args, **_kwargs):
            raise RuntimeError("knn fail")

    monkeypatch.setattr("sklearn.neighbors.NearestNeighbors", _BrokenNN)

    ctx = DummyCtx()
    out = await ccc._create_microenvironments_file(
        adata,
        CellCommunicationParameters(
            method="cellphonedb",
            species="human",
            cell_type_key="cell_type",
        ),
        ctx,
    )

    assert out is None
    assert any("Failed to create microenvironments file" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_create_microenvironments_file_removes_partial_file_on_write_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import tempfile
    from pathlib import Path

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(
        ["T"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )
    real_named_temporary_file = tempfile.NamedTemporaryFile
    created_paths: list[str] = []

    class _BrokenTemporaryFile:
        def __init__(self, **kwargs):
            self._file = real_named_temporary_file(**kwargs)
            self.name = self._file.name
            created_paths.append(self.name)

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            self._file.close()

        def write(self, _value):
            raise OSError("disk full")

    monkeypatch.setattr(tempfile, "NamedTemporaryFile", _BrokenTemporaryFile)

    ctx = DummyCtx()
    out = await ccc._create_microenvironments_file(
        adata,
        CellCommunicationParameters(
            method="cellphonedb",
            species="human",
            cell_type_key="cell_type",
            cellphonedb_spatial_radius=10.0,
        ),
        ctx,
    )

    assert out is None
    assert len(created_paths) == 1
    assert not Path(created_paths[0]).exists()
    assert any("disk full" in warning for warning in ctx.warnings)


@pytest.mark.asyncio
async def test_analyze_communication_liana_wraps_internal_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setitem(__import__("sys").modules, "liana", type("L", (), {})())

    def _boom(*_args, **_kwargs):
        raise RuntimeError("cluster fail")

    monkeypatch.setattr(ccc, "_run_liana_cluster_analysis", _boom)

    with pytest.raises(ProcessingError, match=r"LIANA\+ analysis failed: cluster fail"):
        await ccc._analyze_communication_liana(
            adata,
            CellCommunicationParameters(
                method="liana",
                species="human",
                cell_type_key="cell_type",
                perform_spatial_analysis=False,
            ),
            DummyCtx(),
        )


def test_run_liana_spatial_analysis_without_pvals_layer_sets_none(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    var = pd.DataFrame(
        {
            "morans": [0.9],
            "morans_pvals": [0.001],
        },
        index=["L1^R1"],
    )
    lrdata = type(
        "LR",
        (),
        {
            "n_vars": 1,
            "var": var,
            "X": np.array([[0.5]] * adata.n_obs),
            "layers": {},
        },
    )()

    fake_liana = type(
        "L",
        (),
        {
            "mt": type(
                "MT", (), {"bivariate": staticmethod(lambda *_args, **_kwargs: lrdata)}
            )
        },
    )()
    monkeypatch.setitem(__import__("sys").modules, "liana", fake_liana)

    out = ccc._run_liana_spatial_analysis(
        adata,
        CellCommunicationParameters(
            method="liana",
            species="human",
            cell_type_key="cell_type",
            liana_global_metric="morans",
            plot_top_pairs=1,
        ),
        DummyCtx(),
        fake_liana,
    )

    assert out.method_data["spatial_pvals"] is None
    assert out.lr_pairs == ["L1^R1"]


@pytest.mark.asyncio
@pytest.mark.parametrize("debug_seed", [None, 17])
async def test_analyze_communication_cellphonedb_success_with_correction_stats(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, debug_seed
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(
        ["T"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )

    means = pd.DataFrame(
        {
            "interacting_pair": ["L1^R1", "L2-R2"],
            "T|B": [0.8, 0.2],
            "B|T": [0.6, 0.1],
        },
        index=["pair1", "pair2"],
    )
    pvals = pd.DataFrame(
        {
            "T|B": [0.001, 0.8],
            "B|T": [0.01, 0.7],
        },
        index=["pair1", "pair2"],
    )

    captured_kwargs: dict[str, object] = {}

    def _call(**kwargs):
        captured_kwargs.update(kwargs)
        if "debug_seed" in kwargs:
            # Simulate CellPhoneDB's process-global RNG mutation.
            np.random.seed(kwargs["debug_seed"])
        return {
            "deconvoluted": pd.DataFrame(),
            "means": means,
            "pvalues": pvals,
            "significant_means": means.copy(),
        }

    fake_cpdb_method = types.SimpleNamespace(call=_call)
    _patch_required_module(
        monkeypatch,
        "cellphonedb.src.core.methods.cpdb_statistical_analysis_method",
        fake_cpdb_method,
    )

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda *_args, **_kwargs: "/tmp/fake_cpdb.zip",
    )

    params = CellCommunicationParameters(
        method="cellphonedb",
        species="human",
        cell_type_key="cell_type",
        cellphonedb_use_microenvironments=False,
        cellphonedb_iterations=5,
        cellphonedb_pvalue=0.05,
        cellphonedb_correction_method="fdr_bh",
        cellphonedb_debug_seed=debug_seed,
        plot_top_pairs=2,
    )

    rng_state_before = np.random.get_state()
    out = await ccc._analyze_communication_cellphonedb(adata, params, DummyCtx())
    rng_state_after = np.random.get_state()

    assert out.method == "cellphonedb"
    assert out.n_pairs == 2
    assert out.n_significant == 1
    assert out.lr_pairs == ["L1^R1", "L2^R2"]
    assert out.top_lr_pairs == ["L1^R1"]
    assert "correction_statistics" in out.statistics
    assert rng_state_after[0] == rng_state_before[0]
    np.testing.assert_array_equal(rng_state_after[1], rng_state_before[1])
    assert rng_state_after[2:] == rng_state_before[2:]
    if debug_seed is None:
        assert "debug_seed" not in captured_kwargs
        assert "threads" not in captured_kwargs
    else:
        assert captured_kwargs["debug_seed"] == debug_seed
        assert captured_kwargs["threads"] == 1


@pytest.mark.asyncio
async def test_analyze_communication_cellphonedb_ranks_top_pairs_by_strength(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """CellPhoneDB returns interactions in database order, not by strength.

    Taking the first rows would make "top pairs" mean "first in the database".
    """
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(
        ["T"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )

    means = pd.DataFrame(
        {
            "interacting_pair": ["L1^R1", "L2^R2", "L3^R3"],
            "T|B": [0.1, 0.5, 0.9],
            "B|T": [0.1, 0.5, 0.9],
        },
        index=["pair1", "pair2", "pair3"],
    )
    pvals = pd.DataFrame(
        {
            "T|B": [0.001, 0.001, 0.001],
            "B|T": [0.001, 0.001, 0.001],
        },
        index=["pair1", "pair2", "pair3"],
    )

    fake_cpdb_method = types.SimpleNamespace(
        call=lambda **_kwargs: {
            "deconvoluted": pd.DataFrame(),
            "means": means,
            "pvalues": pvals,
            "significant_means": means.copy(),
        }
    )
    _patch_required_module(
        monkeypatch,
        "cellphonedb.src.core.methods.cpdb_statistical_analysis_method",
        fake_cpdb_method,
    )
    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda *_args, **_kwargs: "/tmp/fake_cpdb.zip",
    )

    params = CellCommunicationParameters(
        method="cellphonedb",
        species="human",
        cell_type_key="cell_type",
        cellphonedb_use_microenvironments=False,
        cellphonedb_iterations=5,
        cellphonedb_correction_method="none",
        plot_top_pairs=2,
    )

    out = await ccc._analyze_communication_cellphonedb(adata, params, DummyCtx())

    assert out.n_significant == 3
    assert out.top_lr_pairs == ["L3^R3", "L2^R2"]


@pytest.mark.asyncio
async def test_analyze_communication_fastccc_success_standard_path(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(
        ["T"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )

    interactions_strength = pd.DataFrame(
        {"T|B": [0.7, 0.1], "B|T": [0.5, 0.2]}, index=["L1^R1", "L2-R2"]
    )
    pvalues = pd.DataFrame(
        {"T|B": [0.01, 0.8], "B|T": [0.03, 0.6]}, index=["L1^R1", "L2-R2"]
    )
    percentages = pd.DataFrame(
        {"T": [0.8, 0.2], "B": [0.6, 0.4]}, index=["L1^R1", "L2-R2"]
    )

    fake_fastccc = types.ModuleType("fastccc")
    fake_fastccc.statistical_analysis_method = lambda **_kwargs: (
        interactions_strength,
        pvalues,
        percentages,
    )
    _patch_required_dependency(monkeypatch, "fastccc", fake_fastccc)
    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda *_a, **_k: type(
            "Raw",
            (),
            {
                "X": adata.X,
                "var_names": adata.var_names,
                "is_integer_counts": True,
                "has_negatives": False,
                "source": "X",
            },
        )(),
    )

    params = CellCommunicationParameters(
        method="fastccc",
        species="human",
        cell_type_key="cell_type",
        fastccc_use_cauchy=False,
        fastccc_pvalue_threshold=0.05,
        plot_top_pairs=1,
    )

    out = await ccc._analyze_communication_fastccc(adata, params, DummyCtx())

    assert out.method == "fastccc"
    assert out.n_pairs == 2
    assert out.n_significant == 1
    assert out.lr_pairs == ["L1^R1", "L2^R2"]
    assert out.top_lr_pairs == ["L1^R1"]
    assert out.method_data["percentages"].shape == (2, 2)


@pytest.mark.asyncio
async def test_analyze_communication_fastccc_cauchy_without_outputs_raises(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    fake_fastccc = types.ModuleType("fastccc")
    fake_fastccc.Cauchy_combination_of_statistical_analysis_methods = (
        lambda **_kwargs: None
    )
    _patch_required_dependency(monkeypatch, "fastccc", fake_fastccc)
    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda *_a, **_k: type(
            "Raw",
            (),
            {
                "X": adata.X,
                "var_names": adata.var_names,
                "is_integer_counts": True,
                "has_negatives": False,
                "source": "X",
            },
        )(),
    )

    monkeypatch.setattr("glob.glob", lambda _pattern: [])

    params = CellCommunicationParameters(
        method="fastccc",
        species="human",
        cell_type_key="cell_type",
        fastccc_use_cauchy=True,
    )

    with pytest.raises(
        ProcessingError, match="Cauchy combination did not produce output files"
    ):
        await ccc._analyze_communication_fastccc(adata, params, DummyCtx())


@pytest.mark.asyncio
async def test_analyze_communication_cellphonedb_rejects_unexpected_return_format(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    fake_cpdb_method = types.SimpleNamespace(call=lambda **_kwargs: "bad-result")
    _patch_required_module(
        monkeypatch,
        "cellphonedb.src.core.methods.cpdb_statistical_analysis_method",
        fake_cpdb_method,
    )

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda *_args, **_kwargs: "/tmp/fake_cpdb.zip",
    )

    params = CellCommunicationParameters(
        method="cellphonedb",
        species="human",
        cell_type_key="cell_type",
        cellphonedb_use_microenvironments=False,
    )

    with pytest.raises(ProcessingError, match="returned unexpected format"):
        await ccc._analyze_communication_cellphonedb(adata, params, DummyCtx())


@pytest.mark.asyncio
async def test_analyze_communication_cellphonedb_rejects_missing_significant_means(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    means = pd.DataFrame({"interacting_pair": ["L1^R1"]})
    pvals = pd.DataFrame({"T|T": [0.01]})

    fake_cpdb_method = types.SimpleNamespace(
        call=lambda **_kwargs: {
            "means": means,
            "pvalues": pvals,
            # intentionally missing 'significant_means'
        }
    )
    _patch_required_module(
        monkeypatch,
        "cellphonedb.src.core.methods.cpdb_statistical_analysis_method",
        fake_cpdb_method,
    )

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda *_args, **_kwargs: "/tmp/fake_cpdb.zip",
    )

    params = CellCommunicationParameters(
        method="cellphonedb",
        species="human",
        cell_type_key="cell_type",
        cellphonedb_use_microenvironments=False,
    )

    with pytest.raises(DataNotFoundError, match="found no L-R interactions"):
        await ccc._analyze_communication_cellphonedb(adata, params, DummyCtx())


@pytest.mark.asyncio
async def test_analyze_communication_cellphonedb_rejects_non_numeric_pvalues(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    means = pd.DataFrame({"interacting_pair": ["L1^R1"]})
    pvals = pd.DataFrame({"T|T": ["not_numeric"]})

    fake_cpdb_method = types.SimpleNamespace(
        call=lambda **_kwargs: {
            "deconvoluted": pd.DataFrame(),
            "means": means,
            "pvalues": pvals,
            "significant_means": means.copy(),
        }
    )
    _patch_required_module(
        monkeypatch,
        "cellphonedb.src.core.methods.cpdb_statistical_analysis_method",
        fake_cpdb_method,
    )

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda *_args, **_kwargs: "/tmp/fake_cpdb.zip",
    )

    params = CellCommunicationParameters(
        method="cellphonedb",
        species="human",
        cell_type_key="cell_type",
        cellphonedb_use_microenvironments=False,
    )

    with pytest.raises(ProcessingError, match="p-values are not numeric"):
        await ccc._analyze_communication_cellphonedb(adata, params, DummyCtx())


@pytest.mark.asyncio
@pytest.mark.parametrize("method", ["cellchat_r", "fastccc"])
async def test_validate_ccc_params_cellchat_and_fastccc_require_cell_type_key(
    minimal_spatial_adata, method: str
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    await ccc._validate_ccc_params(
        adata,
        CellCommunicationParameters(
            method=method,  # type: ignore[arg-type]
            species="human",
            cell_type_key="cell_type",
        ),
        DummyCtx(),
    )


@pytest.mark.asyncio
async def test_analyze_communication_liana_requires_species_parameter(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)
    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setitem(__import__("sys").modules, "liana", type("L", (), {})())

    params = CellCommunicationParameters(
        method="liana",
        species="human",
        cell_type_key="cell_type",
        perform_spatial_analysis=False,
    ).model_copy(update={"species": None})

    with pytest.raises(ParameterError, match="Species parameter is required"):
        await ccc._analyze_communication_liana(adata, params, DummyCtx())


def test_get_liana_resource_name_mouse_unknown_resource_passthrough():
    assert ccc._get_liana_resource_name("mouse", "custom_resource") == "custom_resource"


@pytest.mark.asyncio
async def test_analyze_communication_cellphonedb_sparse_microenv_and_cleanup_oseror(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.X = sp.csr_matrix(np.asarray(adata.X))
    adata.obsm["X_spatial"] = adata.obsm.pop("spatial")
    adata.obs["cell_type"] = pd.Categorical(
        ["T"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )

    means = pd.DataFrame(
        {"interacting_pair": ["L1^R1"], "T|B": [0.5], "B|T": [0.4]},
        index=["pair1"],
    )
    pvals = pd.DataFrame({"T|B": [0.01], "B|T": [0.02]}, index=["pair1"])

    fake_cpdb_method = types.SimpleNamespace(
        call=lambda **_kwargs: {
            "deconvoluted": pd.DataFrame(),
            "means": means,
            "pvalues": pvals,
            "significant_means": means.copy(),
        }
    )
    _patch_required_module(
        monkeypatch,
        "cellphonedb.src.core.methods.cpdb_statistical_analysis_method",
        fake_cpdb_method,
    )

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda *_args, **_kwargs: "/tmp/fake_cpdb.zip",
    )
    monkeypatch.setattr(
        ccc,
        "_create_microenvironments_file",
        lambda *_args, **_kwargs: __import__("asyncio").sleep(
            0, result="/tmp/fake_microenv.txt"
        ),
    )
    monkeypatch.setattr(
        "os.remove", lambda _p: (_ for _ in ()).throw(OSError("cleanup fail"))
    )

    params = CellCommunicationParameters(
        method="cellphonedb",
        species="human",
        cell_type_key="cell_type",
        cellphonedb_use_microenvironments=True,
        cellphonedb_pvalue=0.05,
        plot_top_pairs=1,
    )

    out = await ccc._analyze_communication_cellphonedb(adata, params, DummyCtx())
    assert out.method == "cellphonedb"
    assert out.statistics["microenvironments_used"] is True


@pytest.mark.asyncio
async def test_analyze_communication_cellphonedb_database_setup_error_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    fake_cpdb_method = types.SimpleNamespace(call=lambda **_kwargs: {})
    _patch_required_module(
        monkeypatch,
        "cellphonedb.src.core.methods.cpdb_statistical_analysis_method",
        fake_cpdb_method,
    )

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("db boom")),
    )

    with pytest.raises(DependencyError, match="database setup failed"):
        await ccc._analyze_communication_cellphonedb(
            adata,
            CellCommunicationParameters(
                method="cellphonedb",
                species="human",
                cell_type_key="cell_type",
                cellphonedb_use_microenvironments=False,
            ),
            DummyCtx(),
        )


@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("exc", "msg"),
    [
        (KeyError("no interactions"), "found no L-R interactions"),
        (RuntimeError("api boom"), "analysis failed"),
    ],
)
async def test_analyze_communication_cellphonedb_api_error_branches(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, exc: Exception, msg: str
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    def _raise_call(**_kwargs):
        raise exc

    fake_cpdb_method = types.SimpleNamespace(call=_raise_call)
    _patch_required_module(
        monkeypatch,
        "cellphonedb.src.core.methods.cpdb_statistical_analysis_method",
        fake_cpdb_method,
    )

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda *_args, **_kwargs: "/tmp/fake_cpdb.zip",
    )

    with pytest.raises(ProcessingError, match=msg):
        await ccc._analyze_communication_cellphonedb(
            adata,
            CellCommunicationParameters(
                method="cellphonedb",
                species="human",
                cell_type_key="cell_type",
                cellphonedb_use_microenvironments=False,
            ),
            DummyCtx(),
        )


@pytest.mark.asyncio
async def test_analyze_communication_cellphonedb_handles_missing_pvalues_and_index_lr_pairs(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    means_missing_pvals = pd.DataFrame({"T|T": [0.1]}, index=["L1^R1"])
    pvals_high = pd.DataFrame({"T|T": [0.9], "B|B": [0.8]}, index=["L1^R1"])

    calls = {"i": 0}

    def _call(**_kwargs):
        calls["i"] += 1
        if calls["i"] == 1:
            return {
                "deconvoluted": pd.DataFrame(),
                "means": means_missing_pvals,
                "pvalues": None,
                "significant_means": means_missing_pvals.copy(),
            }
        return {
            "deconvoluted": pd.DataFrame(),
            "means": means_missing_pvals,
            "pvalues": pvals_high,
            "significant_means": means_missing_pvals.copy(),
        }

    fake_cpdb_method = types.SimpleNamespace(call=_call)
    _patch_required_module(
        monkeypatch,
        "cellphonedb.src.core.methods.cpdb_statistical_analysis_method",
        fake_cpdb_method,
    )

    monkeypatch.setattr(ccc, "require", _required_dependency)
    monkeypatch.setattr(
        ccc,
        "_ensure_cellphonedb_database",
        lambda *_args, **_kwargs: "/tmp/fake_cpdb.zip",
    )

    with pytest.raises(DataNotFoundError, match="p-values unavailable"):
        await ccc._analyze_communication_cellphonedb(
            adata,
            CellCommunicationParameters(
                method="cellphonedb",
                species="human",
                cell_type_key="cell_type",
                cellphonedb_use_microenvironments=False,
            ),
            DummyCtx(),
        )

    ctx = DummyCtx()
    out = await ccc._analyze_communication_cellphonedb(
        adata,
        CellCommunicationParameters(
            method="cellphonedb",
            species="human",
            cell_type_key="cell_type",
            cellphonedb_use_microenvironments=False,
            cellphonedb_correction_method="none",
            cellphonedb_pvalue=0.05,
        ),
        ctx,
    )
    assert out.lr_pairs == ["L1^R1"]
    assert out.n_significant == 0
    assert any("Multiple testing correction disabled" in w for w in ctx.warnings)
    assert any("No significant interactions found" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_create_microenvironments_file_auto_radius_branch(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(
        ["T"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )

    path = await ccc._create_microenvironments_file(
        adata,
        CellCommunicationParameters(
            method="cellphonedb",
            species="human",
            cell_type_key="cell_type",
            cellphonedb_spatial_radius=None,
        ),
        DummyCtx(),
    )
    assert path is not None
    __import__("os").remove(path)


@pytest.mark.asyncio
async def test_create_microenvironments_file_auto_radius_handles_small_dataset(
    minimal_spatial_adata,
):
    from pathlib import Path

    adata = minimal_spatial_adata[:4].copy()
    adata.obs["cell_type"] = pd.Categorical(["T", "T", "B", "B"])
    adata.obsm["spatial"] = np.column_stack(
        [np.arange(4, dtype=float), np.zeros(4, dtype=float)]
    )

    path = await ccc._create_microenvironments_file(
        adata,
        CellCommunicationParameters(
            method="cellphonedb",
            species="human",
            cell_type_key="cell_type",
            cellphonedb_spatial_radius=None,
        ),
        DummyCtx(),
    )

    assert path is not None
    try:
        contents = Path(path).read_text(encoding="utf-8")
        assert "T\tmicroenv_" in contents
        assert "B\tmicroenv_" in contents
    finally:
        Path(path).unlink(missing_ok=True)


@pytest.mark.asyncio
async def test_analyze_communication_fastccc_normalizes_when_values_too_large(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import sys
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    interactions_strength = pd.DataFrame({"T|T": [0.9]}, index=["L1^R1"])
    pvalues = pd.DataFrame({"T|T": [0.01]}, index=["L1^R1"])

    fake_fastccc = types.ModuleType("fastccc")
    fake_fastccc.statistical_analysis_method = lambda **_kwargs: (
        interactions_strength,
        pvalues,
        None,
    )
    _patch_required_dependency(monkeypatch, "fastccc", fake_fastccc)

    calls = {"norm": 0, "log1p": 0}
    fake_scanpy = types.ModuleType("scanpy")
    fake_scanpy.pp = types.SimpleNamespace(
        normalize_total=lambda *_a, **_k: calls.__setitem__("norm", calls["norm"] + 1),
        log1p=lambda *_a, **_k: calls.__setitem__("log1p", calls["log1p"] + 1),
    )
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda *_a, **_k: type(
            "Raw",
            (),
            {
                "X": np.full((adata.n_obs, adata.n_vars), 100.0),
                "var_names": adata.var_names,
            },
        )(),
    )
    out = await ccc._analyze_communication_fastccc(
        adata,
        CellCommunicationParameters(
            method="fastccc",
            species="human",
            cell_type_key="cell_type",
            fastccc_use_cauchy=False,
        ),
        DummyCtx(),
    )

    assert out.method == "fastccc"
    assert calls["norm"] == 1
    assert calls["log1p"] == 1


@pytest.mark.asyncio
async def test_analyze_communication_fastccc_missing_database_files_raises(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    fake_fastccc = types.ModuleType("fastccc")
    fake_fastccc.statistical_analysis_method = lambda **_kwargs: (
        pd.DataFrame(),
        pd.DataFrame(),
        None,
    )
    _patch_required_dependency(
        monkeypatch,
        "fastccc",
        fake_fastccc,
        configure_fastccc_database=False,
    )

    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda *_a, **_k: type(
            "Raw",
            (),
            {
                "X": np.asarray(adata.X),
                "var_names": adata.var_names,
                "is_integer_counts": True,
                "has_negatives": False,
                "source": "X",
            },
        )(),
    )
    monkeypatch.setattr(
        ccc,
        "_ensure_fastccc_database_dir",
        lambda _output_dir: (_ for _ in ()).throw(
            ProcessingError("FastCCC requires a complete CellPhoneDB database")
        ),
    )

    with pytest.raises(
        ProcessingError, match="requires a complete CellPhoneDB database"
    ):
        await ccc._analyze_communication_fastccc(
            adata,
            CellCommunicationParameters(
                method="fastccc",
                species="human",
                cell_type_key="cell_type",
                fastccc_use_cauchy=False,
            ),
            DummyCtx(),
        )


@pytest.mark.asyncio
async def test_analyze_communication_fastccc_cauchy_reads_saved_outputs_and_none_percentages(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import os
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    def _fake_cauchy(**kwargs):
        save_path = kwargs["save_path"]
        task_id = "task_x"
        pd.DataFrame({"T|T": [0.01]}, index=["L1^R1"]).to_csv(
            os.path.join(save_path, f"{task_id}_Cauchy_pvals.tsv"), sep="\t"
        )
        pd.DataFrame({"T|T": [0.9]}, index=["L1^R1"]).to_csv(
            os.path.join(save_path, f"{task_id}_average_interactions_strength.tsv"),
            sep="\t",
        )

    fake_fastccc = types.ModuleType("fastccc")
    fake_fastccc.Cauchy_combination_of_statistical_analysis_methods = _fake_cauchy
    _patch_required_dependency(monkeypatch, "fastccc", fake_fastccc)

    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda *_a, **_k: type(
            "Raw",
            (),
            {
                "X": np.asarray(adata.X),
                "var_names": adata.var_names,
                "is_integer_counts": True,
                "has_negatives": False,
                "source": "X",
            },
        )(),
    )
    import glob as glob_mod

    _orig_glob = glob_mod.glob
    monkeypatch.setattr(
        "glob.glob",
        lambda pattern: [] if "percents_analysis" in pattern else _orig_glob(pattern),
    )

    out = await ccc._analyze_communication_fastccc(
        adata,
        CellCommunicationParameters(
            method="fastccc",
            species="human",
            cell_type_key="cell_type",
            fastccc_use_cauchy=True,
            fastccc_pvalue_threshold=0.05,
            plot_top_pairs=1,
        ),
        DummyCtx(),
    )

    assert out.method == "fastccc"
    assert out.n_significant == 1
    assert out.method_data["percentages"] is None


@pytest.mark.asyncio
async def test_analyze_communication_fastccc_none_pvalues_sets_zero_significant(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    fake_fastccc = types.ModuleType("fastccc")
    fake_fastccc.statistical_analysis_method = lambda **_kwargs: (
        pd.DataFrame({"T|T": [0.5]}, index=["L1^R1"]),
        None,
        None,
    )
    _patch_required_dependency(monkeypatch, "fastccc", fake_fastccc)

    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda *_a, **_k: type(
            "Raw",
            (),
            {
                "X": np.asarray(adata.X),
                "var_names": adata.var_names,
                "is_integer_counts": True,
                "has_negatives": False,
                "source": "X",
            },
        )(),
    )
    out = await ccc._analyze_communication_fastccc(
        adata,
        CellCommunicationParameters(
            method="fastccc",
            species="human",
            cell_type_key="cell_type",
            fastccc_use_cauchy=False,
        ),
        DummyCtx(),
    )

    assert out.n_significant == 0


def _install_fake_rpy2(
    monkeypatch: pytest.MonkeyPatch,
    *,
    cellchat_genes: list[str],
    n_lr_pairs: int = 2,
    n_significant: int = 3,
    top_pathways: list[str] | None = None,
    top_lr: list[str] | None = None,
):
    import sys
    import types
    from contextlib import contextmanager

    class _Conv:
        def __add__(self, _other):
            return self

    class _LCtx:
        def __enter__(self):
            return None

        def __exit__(self, exc_type, exc, tb):
            return False

    class _RExec:
        def __init__(self):
            self.scripts: list[str] = []
            self.lookup = {
                "cellchat_genes": list(cellchat_genes),
                "n_lr_pairs": np.array([n_lr_pairs]),
                "n_significant": np.array([n_significant]),
                "top_pathways": list(top_pathways or ["PathA", "PathB"]),
                "top_lr": list(top_lr or ["L1^R1", "L2-R2"]),
                "as.data.frame(lr_pairs)": pd.DataFrame(
                    {"interaction_name": ["L1^R1", "L2-R2"]}
                ),
                "net$prob": np.ones((2, 2, 2), dtype=float),
                "net$pval": np.full((2, 2, 2), 0.05, dtype=float),
                "dimnames(net$prob)[[3]]": ["L1^R1", "L2-R2"],
                "rownames(net$prob)": ["T", "B"],
            }

        def __call__(self, code: str):
            key = code.strip()
            if key in self.lookup:
                return self.lookup[key]
            self.scripts.append(code)
            return None

    r_exec = _RExec()
    fake_ro = types.ModuleType("rpy2.robjects")
    fake_ro.r = r_exec
    fake_ro.globalenv = {}
    fake_ro.default_converter = _Conv()
    fake_ro.numpy2ri = types.SimpleNamespace(converter=_Conv())
    fake_ro.pandas2ri = types.SimpleNamespace(converter=_Conv())

    @contextmanager
    def _local_context(*, use_rlock: bool = True):
        del use_rlock
        global_environment = fake_ro.globalenv
        local_environment: dict[str, object] = {}
        fake_ro.globalenv = local_environment
        fake_ro.local_context_calls += 1
        try:
            yield local_environment
        finally:
            fake_ro.last_localenv = local_environment
            fake_ro.globalenv = global_environment
            fake_ro.local_context_exits += 1

    fake_ro.local_context = _local_context
    fake_ro.local_context_calls = 0
    fake_ro.local_context_exits = 0
    fake_ro.last_localenv = None

    fake_conversion = types.ModuleType("rpy2.robjects.conversion")
    fake_conversion.localconverter = lambda _converter: _LCtx()

    class _RLock:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    environment = REnvironment(
        robjects=fake_ro,
        pandas2ri=fake_ro.pandas2ri,
        numpy2ri=fake_ro.numpy2ri,
        packages=types.SimpleNamespace(importr=lambda _package: object()),
        conversion=fake_conversion,
        openrlib=types.SimpleNamespace(rlock=_RLock()),
        anndata2ri=None,
    )
    monkeypatch.setattr(
        ccc,
        "validate_r_environment",
        lambda *_args, **_kwargs: environment,
    )

    fake_rpy2 = types.ModuleType("rpy2")
    fake_rpy2.robjects = fake_ro

    monkeypatch.setitem(sys.modules, "rpy2", fake_rpy2)
    monkeypatch.setitem(sys.modules, "rpy2.robjects", fake_ro)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.conversion", fake_conversion)

    return r_exec, fake_ro


def test_analyze_communication_cellchat_r_non_spatial_success(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(
        ["T"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )

    r_exec, fake_ro = _install_fake_rpy2(
        monkeypatch,
        cellchat_genes=["gene_0", "gene_1", "gene_2", "gene_3"],
    )

    monkeypatch.setattr(ccc, "validate_obs_column", lambda *_a, **_k: None)
    monkeypatch.setattr(ccc, "get_spatial_key", lambda *_a, **_k: None)
    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True: type(
            "Raw",
            (),
            {
                "X": np.asarray(_adata.X),
                "var_names": pd.Index(_adata.var_names),
                "is_integer_counts": True,
                "has_negatives": False,
                "source": "X",
            },
        )(),
    )

    params = CellCommunicationParameters(
        method="cellchat_r",
        species="human",
        cell_type_key="cell_type",
    )

    out = ccc._analyze_communication_cellchat_r(adata, params, DummyCtx())

    assert out.method == "cellchat_r"
    assert out.n_pairs == 2
    assert out.n_significant == 3
    assert out.lr_pairs == ["L1^R1", "L2^R2"]
    assert out.method_data["prob_matrix"].shape == (2, 2, 2)
    assert out.method_data["lr_names"] == ["L1^R1", "L2-R2"]
    assert "createCellChat" in "\n".join(r_exec.scripts)
    assert fake_ro.globalenv == {}
    assert "spatial_locs" not in fake_ro.last_localenv


def test_analyze_communication_cellchat_r_spatial_contact_range_and_db_category(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    # Numeric-like labels trigger CellChat-safe prefixing branch.
    adata.obs["cell_type"] = pd.Categorical(
        ["0"] * (adata.n_obs // 2) + ["1"] * (adata.n_obs - adata.n_obs // 2)
    )

    r_exec, fake_ro = _install_fake_rpy2(
        monkeypatch,
        cellchat_genes=["gene_0", "gene_1", "gene_2", "gene_3"],
    )

    monkeypatch.setattr(ccc, "validate_obs_column", lambda *_a, **_k: None)
    monkeypatch.setattr(ccc, "get_spatial_key", lambda *_a, **_k: "spatial")
    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True: type(
            "Raw",
            (),
            {
                "X": np.asarray(_adata.X),
                "var_names": pd.Index(_adata.var_names),
                "is_integer_counts": True,
                "has_negatives": False,
                "source": "X",
            },
        )(),
    )

    params = CellCommunicationParameters(
        method="cellchat_r",
        species="human",
        cell_type_key="cell_type",
        cellchat_db_category="Secreted Signaling",
        cellchat_contact_range=100.0,
        cellchat_distance_use=True,
    )

    out = ccc._analyze_communication_cellchat_r(adata, params, DummyCtx())

    scripts = "\n".join(r_exec.scripts)
    assert out.statistics["spatial_mode"] is True
    assert "contact.range = 100.0" in scripts
    assert "subsetDB" in scripts
    assert fake_ro.globalenv == {}
    assert "spatial_locs" in fake_ro.last_localenv
    # Labels should be prefixed when numeric-like.
    assert set(fake_ro.last_localenv["meta_df"]["labels"].unique()) == {
        "cluster_0",
        "cluster_1",
    }


def test_analyze_communication_cellchat_r_spatial_uses_contact_knn_when_range_is_none(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    r_exec, _fake_ro = _install_fake_rpy2(
        monkeypatch,
        cellchat_genes=["gene_0", "gene_1", "gene_2", "gene_3"],
    )

    monkeypatch.setattr(ccc, "validate_obs_column", lambda *_a, **_k: None)
    monkeypatch.setattr(ccc, "get_spatial_key", lambda *_a, **_k: "spatial")
    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True: type(
            "Raw",
            (),
            {
                "X": np.asarray(_adata.X),
                "var_names": pd.Index(_adata.var_names),
                "is_integer_counts": True,
                "has_negatives": False,
                "source": "X",
            },
        )(),
    )

    out = ccc._analyze_communication_cellchat_r(
        adata,
        CellCommunicationParameters(
            method="cellchat_r",
            species="human",
            cell_type_key="cell_type",
            cellchat_contact_range=None,
            cellchat_contact_knn_k=9,
        ),
        DummyCtx(),
    )

    assert out.statistics["spatial_mode"] is True
    assert "contact.knn.k = 9" in "\n".join(r_exec.scripts)


def test_analyze_communication_cellchat_r_preserves_no_gene_overlap_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    _r_exec, fake_ro = _install_fake_rpy2(monkeypatch, cellchat_genes=["NOT_IN_DATA"])

    monkeypatch.setattr(ccc, "validate_obs_column", lambda *_a, **_k: None)
    monkeypatch.setattr(ccc, "get_spatial_key", lambda *_a, **_k: None)
    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True: type(
            "Raw",
            (),
            {
                "X": np.asarray(_adata.X),
                "var_names": pd.Index(_adata.var_names),
                "is_integer_counts": True,
                "has_negatives": False,
                "source": "X",
            },
        )(),
    )

    with pytest.raises(DataCompatibilityError, match="No genes overlap"):
        ccc._analyze_communication_cellchat_r(
            adata,
            CellCommunicationParameters(
                method="cellchat_r",
                species="human",
                cell_type_key="cell_type",
            ),
            DummyCtx(),
        )

    assert fake_ro.globalenv == {}
    assert fake_ro.local_context_calls == 1
    assert fake_ro.local_context_exits == 1


@pytest.mark.asyncio
async def test_analyze_communication_fastccc_cauchy_reads_percentages_when_present(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import os
    import types

    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = pd.Categorical(["T"] * adata.n_obs)

    def _fake_cauchy(**kwargs):
        save_path = kwargs["save_path"]
        task_id = "task_y"
        pd.DataFrame({"T|T": [0.01]}, index=["L1^R1"]).to_csv(
            os.path.join(save_path, f"{task_id}_Cauchy_pvals.tsv"), sep="\t"
        )
        pd.DataFrame({"T|T": [0.9]}, index=["L1^R1"]).to_csv(
            os.path.join(save_path, f"{task_id}_average_interactions_strength.tsv"),
            sep="\t",
        )
        pd.DataFrame({"T|T": [0.7]}, index=["L1^R1"]).to_csv(
            os.path.join(save_path, f"{task_id}_percents_analysis.tsv"), sep="\t"
        )

    fake_fastccc = types.ModuleType("fastccc")
    fake_fastccc.Cauchy_combination_of_statistical_analysis_methods = _fake_cauchy
    _patch_required_dependency(monkeypatch, "fastccc", fake_fastccc)

    monkeypatch.setattr(
        ccc,
        "get_raw_data_source",
        lambda *_a, **_k: type(
            "Raw",
            (),
            {"X": np.asarray(adata.X), "var_names": adata.var_names},
        )(),
    )
    out = await ccc._analyze_communication_fastccc(
        adata,
        CellCommunicationParameters(
            method="fastccc",
            species="human",
            cell_type_key="cell_type",
            fastccc_use_cauchy=True,
            fastccc_pvalue_threshold=0.05,
        ),
        DummyCtx(),
    )

    assert out.method_data["percentages"] is not None


@pytest.mark.unit
def test_cellchat_reads_the_matrix_its_own_contract_asks_for():
    """createCellChat documents "a normalized (NOT count) data matrix".

    Nothing downstream normalizes — identifyOverExpressedGenes runs
    presto::wilcoxauc straight off object@data.signaling — and computeCommunProb
    divides by max(data) before a Hill function with a fixed Kh = 0.5. Counts
    put the divisor on one extreme gene and leave every ligand-receptor product
    orders of magnitude below Kh, so the saturation the model exists to express
    never engages: on a lymph node the median interaction probability moved from
    1.3e-07 to 4.9e-05 and TGFB1 went from absent to the top interaction.
    """
    import inspect

    from chatspatial.tools import cell_communication as cc

    source = inspect.getsource(cc._analyze_communication_cellchat_r)

    assert "build_log_normalized_view(" in source
    # The old path prepared counts and warned in the inverted direction: it
    # flagged the input as suspect exactly when it was already normalized.
    assert "get_raw_data_source" not in source
    assert "double-normalization" not in source
