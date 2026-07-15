"""Unit tests for cell communication visualization utilities and routing."""

from __future__ import annotations

import sys
from types import ModuleType

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.cell_communication import CCCStorage, store_ccc_results
from chatspatial.tools.visualization import cell_comm as viz_cc
from chatspatial.utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)


class DummyCtx:
    def __init__(self):
        self.infos: list[str] = []

    async def info(self, msg: str):
        self.infos.append(msg)


@pytest.fixture(autouse=True)
def _close_figures_after_test():
    yield
    plt.close("all")


def _store_cluster_ccc(adata, *, method: str = "cellphonedb") -> None:
    results = pd.DataFrame(
        {
            "T|B": [0.8, 0.2],
            "B|T": [0.1, 0.6],
        },
        index=["L1_R1", "L2-R2"],
    )
    pvals = pd.DataFrame(
        {
            "T|B": [0.05, 0.2],
            "B|T": [0.01, 0.3],
        },
        index=results.index,
    )
    storage = CCCStorage(
        method=method,
        analysis_type="cluster",
        species="human",
        database="cellphonedb",
        lr_pairs=["L1^R1", "L2^R2"],
        top_lr_pairs=["L1^R1"],
        n_pairs=2,
        n_significant=1,
        results=results,
        pvalues=pvals,
    )
    store_ccc_results(adata, storage)


def _mock_liana_results(include_scores: bool = True) -> pd.DataFrame:
    data: dict[str, list] = {
        "source": ["T", "B", "T"],
        "target": ["B", "T", "T"],
        "ligand_complex": ["L1", "L2", "L1"],
        "receptor_complex": ["R1", "R2", "R1"],
    }
    if include_scores:
        data["magnitude_rank"] = [0.1, 0.4, 0.2]
        data["lr_means"] = [0.9, 0.4, 0.7]
    return pd.DataFrame(data)


def _mock_cc_data(
    results: pd.DataFrame,
    *,
    analysis_type: str = "cluster",
    lr_pairs: list[str] | None = None,
    spatial_scores: np.ndarray | None = None,
):
    return viz_cc.CellCommunicationData(
        results=results,
        method="cellphonedb",
        analysis_type=analysis_type,
        lr_pairs=lr_pairs or ["L1^R1", "L2^R2"],
        pvalues=None,
        spatial_scores=spatial_scores,
        spatial_pvals=None,
        source_labels=["T", "B"],
        target_labels=["T", "B"],
        method_data=None,
    )


class _FakePlotnine:
    def __init__(self, should_fail: bool = False):
        self.should_fail = should_fail

    def draw(self):
        if self.should_fail:
            raise RuntimeError("draw failed")
        fig, _ = plt.subplots()
        return fig


def test_parse_lr_pair_supports_common_formats():
    # Canonical "^" separator — always safe
    assert viz_cc._parse_lr_pair("LIG^REC") == ("LIG", "REC")
    # Single underscore — unambiguous
    assert viz_cc._parse_lr_pair("LIG_REC") == ("LIG", "REC")
    # Prefix removal + single dash
    assert viz_cc._parse_lr_pair("complex:LIG-REC") == ("LIG", "REC")
    # Multi-underscore (complex names) — cannot safely split
    assert viz_cc._parse_lr_pair("ITGAL_ICAM1_ITGB2") == (
        "ITGAL_ICAM1_ITGB2",
        "ITGAL_ICAM1_ITGB2",
    )
    # "^" with complex names — safe because ^ is the canonical separator
    assert viz_cc._parse_lr_pair("ITGAL^ICAM1_ITGB2") == ("ITGAL", "ICAM1_ITGB2")


def test_parse_lr_pair_falls_back_when_no_separator():
    assert viz_cc._parse_lr_pair("NOSPLIT") == ("NOSPLIT", "NOSPLIT")


def test_convert_to_liana_format_matrix_contract():
    results = pd.DataFrame(
        {
            "interacting_pair": ["CXCL12_CXCR4", "CCL5^CCR5"],
            "A|B": [1.2, 0.0],
            "B|A": [0.5, np.nan],
        }
    )
    pvals = pd.DataFrame({"A|B": [0.01, 0.2], "B|A": [0.03, 0.4]})

    out = viz_cc._convert_to_liana_format(results, pvals, method="cellphonedb")

    assert not out.empty
    assert {"source", "target", "ligand_complex", "receptor_complex"}.issubset(
        out.columns
    )
    assert set(out["source"]) == {"A", "B"}
    assert set(out["target"]) == {"A", "B"}
    assert out["lr_means"].min() >= 0


def test_convert_to_liana_format_handles_empty_and_passthrough_cases():
    empty = viz_cc._convert_to_liana_format(
        pd.DataFrame(), pvalues=None, method="cellphonedb"
    )
    assert empty.empty

    liana_ready = _mock_liana_results()
    out_liana = viz_cc._convert_to_liana_format(
        liana_ready, pvalues=None, method="liana"
    )
    out_passthrough = viz_cc._convert_to_liana_format(
        liana_ready, pvalues=None, method="cellphonedb"
    )
    assert out_liana.equals(liana_ready)
    assert out_passthrough.equals(liana_ready)
    assert out_liana is not liana_ready
    assert out_passthrough is not liana_ready


def test_convert_to_liana_format_cellchat_3d_contract():
    results = pd.DataFrame(
        {
            "interaction_name": ["L1_R1", "L2_R2"],
            "ligand": ["L1", "L2"],
            "receptor": ["R1", "R2"],
        }
    )
    prob = np.zeros((2, 2, 2), dtype=float)
    prob[0, 1, 0] = 0.6
    prob[1, 0, 1] = 0.2
    pval = np.ones((2, 2, 2), dtype=float)
    pval[0, 1, 0] = 0.01
    pval[1, 0, 1] = 0.2

    out = viz_cc._convert_to_liana_format(
        results=results,
        pvalues=None,
        method="cellchat_r",
        method_data={
            "prob_matrix": prob,
            "pval_matrix": pval,
            "cell_type_names": ["A", "B"],
        },
    )

    assert len(out) == 2
    assert set(out["ligand_complex"]) == {"L1", "L2"}
    assert out["magnitude_rank"].between(0, 1).all()


def test_cellchat_3d_uses_matrix_position_for_lr_metadata():
    results = pd.DataFrame(
        {
            "interaction_name": ["CXCL12_CXCR4"],
            "ligand": ["CXCL12"],
            "receptor": ["CXCR4"],
        },
        index=["interaction_42"],
    )

    out = viz_cc._cellchat_3d_to_liana_format(
        results,
        {
            "prob_matrix": np.ones((1, 1, 1), dtype=float),
            "cell_type_names": "T cell",
        },
    )

    assert out[["ligand_complex", "receptor_complex"]].iloc[0].tolist() == [
        "CXCL12",
        "CXCR4",
    ]


def test_cellchat_3d_prefers_named_interaction_axis_over_result_row_order():
    results = pd.DataFrame(
        {
            "interaction_name": ["L2_R2", "L1_R1"],
            "ligand": ["L2", "L1"],
            "receptor": ["R2", "R1"],
        }
    )
    probability = np.zeros((1, 1, 2), dtype=float)
    probability[0, 0, 0] = 0.8

    out = viz_cc._cellchat_3d_to_liana_format(
        results,
        {
            "prob_matrix": probability,
            "cell_type_names": ["T"],
            "lr_names": ["L1_R1", "L2_R2"],
        },
    )

    assert out[["ligand_complex", "receptor_complex"]].iloc[0].tolist() == [
        "L1",
        "R1",
    ]


def test_cellchat_3d_validates_matrix_invariants():
    results = pd.DataFrame({"interaction_name": ["L1_R1"]})

    with pytest.raises(ProcessingError, match="must have shape"):
        viz_cc._cellchat_3d_to_liana_format(
            results,
            {"prob_matrix": np.ones((1, 1)), "cell_type_names": ["A"]},
        )

    with pytest.raises(ProcessingError, match="labels do not match"):
        viz_cc._cellchat_3d_to_liana_format(
            results,
            {"prob_matrix": np.ones((2, 2, 1)), "cell_type_names": ["A"]},
        )

    with pytest.raises(ProcessingError, match="p-value matrix must match"):
        viz_cc._cellchat_3d_to_liana_format(
            results,
            {
                "prob_matrix": np.ones((1, 1, 1)),
                "pval_matrix": np.ones((1, 1, 2)),
                "cell_type_names": ["A"],
            },
        )

    with pytest.raises(ProcessingError, match="metadata does not match"):
        viz_cc._cellchat_3d_to_liana_format(
            results,
            {"prob_matrix": np.ones((1, 1, 2)), "cell_type_names": ["A"]},
        )


def test_cell_communication_data_uses_immutable_empty_top_pairs():
    data = viz_cc.CellCommunicationData(
        results=pd.DataFrame(),
        method="liana",
        analysis_type="cluster",
        lr_pairs=(),
    )

    assert data.top_lr_pairs == ()


def test_cellchat_3d_to_liana_format_handles_nan_pval_and_empty_rows():
    results = pd.DataFrame(
        {
            "interaction_name": ["L1_R1"],
            "ligand": ["L1"],
            "receptor": ["R1"],
        }
    )
    prob = np.full((1, 1, 1), 0.7, dtype=float)
    pval = np.full((1, 1, 1), np.nan, dtype=float)

    out = viz_cc._cellchat_3d_to_liana_format(
        results=results,
        method_data={
            "prob_matrix": prob,
            "pval_matrix": pval,
            "cell_type_names": ["A"],
        },
    )
    assert len(out) == 1
    assert out.iloc[0]["ligand_complex"] == "L1"
    assert out.iloc[0]["receptor_complex"] == "R1"
    assert out.iloc[0]["magnitude_rank"] == 1.0

    empty = viz_cc._cellchat_3d_to_liana_format(
        results=results,
        method_data={
            "prob_matrix": np.zeros((1, 1, 1), dtype=float),
            "pval_matrix": np.ones((1, 1, 1), dtype=float),
            "cell_type_names": ["A"],
        },
    )
    assert empty.empty


def test_matrix_to_liana_format_fallback_to_all_numeric_and_handles_nan_and_keyerror():
    # Use "|" separator for cell type pairs (safe separator).
    # Columns without "|" are skipped to avoid wrong splits on cell type names.
    nan_rank = viz_cc._matrix_to_liana_format(
        results=pd.DataFrame({"A|B": [1.0]}, index=["PAIR_1"]),
        pvalues=pd.DataFrame({"A|B": [np.nan]}, index=["PAIR_1"]),
        method="cellphonedb",
    )
    assert len(nan_rank) == 1
    assert nan_rank.iloc[0]["source"] == "A"
    assert nan_rank.iloc[0]["target"] == "B"
    assert nan_rank.iloc[0]["magnitude_rank"] == 1.0

    keyerr_rank = viz_cc._matrix_to_liana_format(
        results=pd.DataFrame({"A|B": [0.8]}, index=["PAIR_1"]),
        pvalues=pd.DataFrame({"A|B": [0.2]}, index=["OTHER_PAIR"]),
        method="cellphonedb",
    )
    assert len(keyerr_rank) == 1
    assert keyerr_rank.iloc[0]["magnitude_rank"] == 1.0


def test_matrix_to_liana_format_skips_underscore_cell_pair_columns():
    """Columns with '_' but no '|' are skipped — '_' in cell type names is ambiguous."""
    underscore_cols = viz_cc._matrix_to_liana_format(
        results=pd.DataFrame({"T_cell_B_cell": [0.5]}, index=["L1_R1"]),
        pvalues=None,
        method="cellphonedb",
    )
    assert underscore_cols.empty


def test_matrix_to_liana_format_returns_empty_for_no_numeric_and_unparsable_pairs():
    no_numeric = viz_cc._matrix_to_liana_format(
        results=pd.DataFrame({"interacting_pair": ["L1_R1"]}),
        pvalues=None,
        method="cellphonedb",
    )
    assert no_numeric.empty

    unparsable_pairs = viz_cc._matrix_to_liana_format(
        results=pd.DataFrame({"AB": [0.5]}, index=["L1_R1"]),
        pvalues=None,
        method="cellphonedb",
    )
    assert unparsable_pairs.empty


def test_matrix_to_liana_format_skips_invalid_pipe_split():
    invalid_pipe = viz_cc._matrix_to_liana_format(
        results=pd.DataFrame({"A|B|C": [0.4]}, index=["L1_R1"]),
        pvalues=None,
        method="cellphonedb",
    )
    assert invalid_pipe.empty


@pytest.mark.asyncio
async def test_get_cell_communication_data_requires_stored_results(minimal_spatial_adata):
    with pytest.raises(DataNotFoundError, match="No cell communication results found"):
        await viz_cc.get_cell_communication_data(minimal_spatial_adata)


@pytest.mark.asyncio
async def test_get_cell_communication_data_converts_and_extracts_labels(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    _store_cluster_ccc(adata)
    ctx = DummyCtx()

    data = await viz_cc.get_cell_communication_data(adata, context=ctx)

    assert data.method == "cellphonedb"
    assert not data.results.empty
    assert set(data.source_labels) == {"T", "B"}
    assert set(data.target_labels) == {"T", "B"}
    assert data.spatial_scores is None
    assert any("converted to LIANA format" in m for m in ctx.infos)


@pytest.mark.asyncio
async def test_get_cell_communication_data_rejects_missing_result_table(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    _store_cluster_ccc(adata)
    adata.uns["ccc"]["results"] = None

    with pytest.raises(DataCompatibilityError, match="DataFrame result table"):
        await viz_cc.get_cell_communication_data(adata)


@pytest.mark.asyncio
async def test_get_cell_communication_data_ignores_stale_spatial_data_for_cluster_result(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    _store_cluster_ccc(adata)
    adata.obsm["ccc_spatial_scores"] = np.ones((adata.n_obs, 2))

    data = await viz_cc.get_cell_communication_data(adata)

    assert data.analysis_type == "cluster"
    assert data.spatial_scores is None


@pytest.mark.asyncio
async def test_get_cell_communication_data_validates_spatial_matrix_axes(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    _store_cluster_ccc(adata)
    adata.uns["ccc"]["analysis_type"] = "spatial"
    adata.obsm["ccc_spatial_scores"] = np.ones((adata.n_obs, 1))

    expected_shape = rf"expected \({adata.n_obs}, 2\)"
    with pytest.raises(DataCompatibilityError, match=expected_shape):
        await viz_cc.get_cell_communication_data(adata)


@pytest.mark.asyncio
async def test_create_cell_communication_visualization_default_routes_to_dotplot(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    _store_cluster_ccc(adata)

    sentinel = object()
    called: dict[str, bool] = {}

    def fake_dotplot(*_args, **_kwargs):
        called["dotplot"] = True
        return sentinel

    monkeypatch.setattr(viz_cc, "_create_unified_dotplot", fake_dotplot)

    out = await viz_cc.create_cell_communication_visualization(
        adata, VisualizationParameters(plot_type="communication")
    )

    assert out is sentinel
    assert called["dotplot"] is True


@pytest.mark.asyncio
async def test_create_cell_communication_visualization_defaults_to_dotplot_when_subtype_none_for_cluster(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    data = _mock_cc_data(_mock_liana_results(), analysis_type="cluster")
    sentinel = object()

    async def _fake_get(*_args, **_kwargs):
        return data

    monkeypatch.setattr(viz_cc, "get_cell_communication_data", _fake_get)
    monkeypatch.setattr(viz_cc, "_create_unified_dotplot", lambda *_a, **_k: sentinel)

    params = VisualizationParameters(plot_type="communication", subtype="tileplot")
    params.subtype = None
    out = await viz_cc.create_cell_communication_visualization(adata, params)
    assert out is sentinel


@pytest.mark.asyncio
async def test_create_cell_communication_visualization_spatial_requires_spatial_analysis(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    _store_cluster_ccc(adata, method="liana")

    with pytest.raises(ParameterError, match="requires spatial analysis"):
        await viz_cc.create_cell_communication_visualization(
            adata,
            VisualizationParameters(plot_type="communication", subtype="spatial"),
        )


@pytest.mark.asyncio
async def test_cluster_only_communication_visualizations_reject_spatial_results(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    data = _mock_cc_data(
        pd.DataFrame({"morans": [0.8]}, index=["L1^R1"]),
        analysis_type="spatial",
        spatial_scores=np.ones((adata.n_obs, 1), dtype=float),
    )

    async def _fake_get(*_args, **_kwargs):
        return data

    monkeypatch.setattr(viz_cc, "get_cell_communication_data", _fake_get)

    with pytest.raises(ParameterError, match="requires cluster-level communication results"):
        await viz_cc.create_cell_communication_visualization(
            adata,
            VisualizationParameters(plot_type="communication", subtype="circle_plot"),
        )


@pytest.mark.asyncio
async def test_create_cell_communication_visualization_unknown_subtype_error(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    _store_cluster_ccc(adata)

    with pytest.raises(ParameterError, match="Unknown visualization type"):
        await viz_cc.create_cell_communication_visualization(
            adata,
            VisualizationParameters(plot_type="communication", subtype="unknown"),
        )


@pytest.mark.asyncio
async def test_create_cell_communication_visualization_logs_context_and_spatial_unknown_lists_available(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    data = _mock_cc_data(
        _mock_liana_results(),
        analysis_type="spatial",
        spatial_scores=np.ones((adata.n_obs, 2), dtype=float),
    )
    ctx = DummyCtx()

    async def _fake_get(*_args, **_kwargs):
        return data

    monkeypatch.setattr(viz_cc, "get_cell_communication_data", _fake_get)

    with pytest.raises(ParameterError, match="Available: spatial, dotplot, tileplot, circle_plot"):
        await viz_cc.create_cell_communication_visualization(
            adata,
            VisualizationParameters(plot_type="communication", subtype="bad"),
            context=ctx,
        )

    assert any("Creating cell communication visualization" in m for m in ctx.infos)
    assert any("Using cellphonedb results" in m for m in ctx.infos)


def test_cellchat_3d_to_liana_format_handles_empty_and_normalization():
    empty = viz_cc._cellchat_3d_to_liana_format(
        results=pd.DataFrame({"interaction_name": ["L1_R1"]}),
        method_data={"pval_matrix": np.ones((1, 1, 1))},
    )
    assert empty.empty

    prob = np.zeros((2, 2, 1), dtype=float)
    prob[0, 1, 0] = 0.8
    prob[1, 0, 0] = 0.3
    pval = np.zeros((2, 2, 1), dtype=float)
    pval[0, 1, 0] = 2.0
    pval[1, 0, 0] = 4.0

    out = viz_cc._cellchat_3d_to_liana_format(
        results=pd.DataFrame(
            {
                "interaction_name": ["L1_R1"],
                "ligand": ["L1"],
                "receptor": ["R1"],
            }
        ),
        method_data={
            "prob_matrix": prob,
            "pval_matrix": pval,
            "cell_type_names": np.array(["A", "B"]),
        },
    )
    assert len(out) == 2
    assert out["magnitude_rank"].max() == pytest.approx(1.0)


def test_matrix_to_liana_format_parses_cell_pairs_and_normalizes():
    results = pd.DataFrame(
        {
            "interacting_pair": ["L1_R1", "L2-R2"],
            "A|B": [0.7, np.nan],
            "B|A": [0.5, 0.4],
        }
    )
    pvals = pd.DataFrame({"A|B": [3.0, 1.0], "B|A": [9.0, 6.0]})

    out = viz_cc._matrix_to_liana_format(results, pvals, method="cellphonedb")
    assert not out.empty
    assert set(out["source"]) == {"A", "B"}
    assert set(out["target"]) == {"A", "B"}
    assert out["magnitude_rank"].max() == pytest.approx(1.0)


@pytest.mark.asyncio
async def test_create_cell_communication_visualization_defaults_to_spatial_for_spatial_analysis(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    sentinel = object()
    data = _mock_cc_data(
        _mock_liana_results(),
        analysis_type="spatial",
        spatial_scores=np.ones((adata.n_obs, 2), dtype=float),
    )

    async def _fake_get(*_args, **_kwargs):
        return data

    monkeypatch.setattr(viz_cc, "get_cell_communication_data", _fake_get)
    monkeypatch.setattr(viz_cc, "_create_spatial_lr_visualization", lambda *_a, **_k: sentinel)

    params = VisualizationParameters(plot_type="communication", subtype="dotplot")
    params.subtype = None
    out = await viz_cc.create_cell_communication_visualization(
        adata,
        params,
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_create_cell_communication_visualization_routes_tile_and_circle(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    data = _mock_cc_data(_mock_liana_results(), analysis_type="cluster")

    async def _fake_get(*_args, **_kwargs):
        return data

    sentinel_tile = object()
    sentinel_circle = object()

    monkeypatch.setattr(viz_cc, "get_cell_communication_data", _fake_get)
    monkeypatch.setattr(viz_cc, "_create_unified_tileplot", lambda *_a, **_k: sentinel_tile)
    monkeypatch.setattr(viz_cc, "_create_unified_circle_plot", lambda *_a, **_k: sentinel_circle)

    tile = await viz_cc.create_cell_communication_visualization(
        adata, VisualizationParameters(plot_type="communication", subtype="tileplot")
    )
    circle = await viz_cc.create_cell_communication_visualization(
        adata, VisualizationParameters(plot_type="communication", subtype="circle_plot")
    )
    assert tile is sentinel_tile
    assert circle is sentinel_circle


def test_create_spatial_lr_visualization_requires_scores(minimal_spatial_adata):
    data = _mock_cc_data(_mock_liana_results(), analysis_type="spatial", spatial_scores=None)
    with pytest.raises(DataNotFoundError, match="No spatial communication scores found"):
        viz_cc._create_spatial_lr_visualization(
            minimal_spatial_adata.copy(),
            data,
            VisualizationParameters(plot_type="communication", subtype="spatial"),
        )


def test_create_spatial_lr_visualization_uses_metric_and_handles_missing_pair_column(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    results = pd.DataFrame(
        {"morans": [0.3, 0.9]},
        index=["L1_R1", "L2_R2"],
    )
    data = _mock_cc_data(
        results,
        analysis_type="spatial",
        lr_pairs=["L1^R1", "L2^R2"],
        spatial_scores=np.ones((adata.n_obs, 1), dtype=float),
    )
    fig = viz_cc._create_spatial_lr_visualization(
        adata,
        data,
        VisualizationParameters(
            plot_type="communication",
            subtype="spatial",
            plot_top_pairs=2,
        ),
    )
    assert fig._suptitle is not None
    assert "Spatial Cell Communication" in fig._suptitle.get_text()
    fig.clf()


def test_create_spatial_lr_visualization_handles_no_metric_with_single_panel(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    data = _mock_cc_data(
        pd.DataFrame({"note": [1]}, index=["L1_R1"]),
        analysis_type="spatial",
        lr_pairs=["L1^R1"],
        spatial_scores=np.ones((adata.n_obs, 1), dtype=float),
    )
    fig = viz_cc._create_spatial_lr_visualization(
        adata,
        data,
        VisualizationParameters(plot_type="communication", subtype="spatial"),
    )
    assert len(fig.axes) >= 1
    fig.clf()


def test_create_spatial_lr_visualization_falls_back_when_top_pairs_not_in_lr_list_and_hides_extra_axes(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    data = _mock_cc_data(
        pd.DataFrame({"morans": [0.9, 0.8, 0.7, 0.6]}, index=["X1", "X2", "X3", "X4"]),
        analysis_type="spatial",
        lr_pairs=["L1^R1", "L2^R2", "L3^R3", "L4^R4"],
        spatial_scores=np.ones((adata.n_obs, 4), dtype=float),
    )
    fig = viz_cc._create_spatial_lr_visualization(
        adata,
        data,
        VisualizationParameters(
            plot_type="communication",
            subtype="spatial",
            plot_top_pairs=4,
        ),
    )
    assert sum(not ax.get_visible() for ax in fig.axes) >= 2
    fig.clf()


def test_create_spatial_lr_visualization_errors_for_empty_selected_pairs_and_missing_scores(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    data = viz_cc.CellCommunicationData(
        results=pd.DataFrame(),
        method="cellphonedb",
        analysis_type="spatial",
        lr_pairs=["L1^R1"],
        pvalues=None,
        spatial_scores=np.ones((adata.n_obs, 1), dtype=float),
        spatial_pvals=None,
        source_labels=None,
        target_labels=None,
        method_data=None,
    )
    params = VisualizationParameters(
        plot_type="communication",
        subtype="spatial",
        plot_top_pairs=1,
    )
    params.plot_top_pairs = -1
    with pytest.raises(DataNotFoundError, match="No LR pairs found in spatial results."):
        viz_cc._create_spatial_lr_visualization(
            adata,
            data,
            params,
        )

    no_scores_data = _mock_cc_data(
        pd.DataFrame(),
        analysis_type="spatial",
        lr_pairs=[],
        spatial_scores=None,
    )
    with pytest.raises(DataNotFoundError, match="No spatial communication scores found"):
        viz_cc._create_spatial_lr_visualization(
            adata,
            no_scores_data,
            VisualizationParameters(plot_type="communication", subtype="spatial"),
        )


def test_create_unified_dotplot_validates_required_columns(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(viz_cc, "require", lambda *_a, **_k: None)
    fake_liana = ModuleType("liana")
    fake_liana.pl = ModuleType("liana.pl")
    fake_liana.pl.dotplot = lambda **_kwargs: _FakePlotnine()
    monkeypatch.setitem(sys.modules, "liana", fake_liana)
    bad = pd.DataFrame({"source": ["A"], "target": ["B"]})
    data = _mock_cc_data(bad)
    with pytest.raises(DataNotFoundError, match="Missing required columns"):
        viz_cc._create_unified_dotplot(
            data,
            VisualizationParameters(plot_type="communication", subtype="dotplot"),
        )


def test_create_unified_dotplot_rejects_when_no_score_columns(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(viz_cc, "require", lambda *_a, **_k: None)
    df = _mock_liana_results(include_scores=False)
    data = _mock_cc_data(df)

    fake_liana = ModuleType("liana")
    fake_liana.pl = ModuleType("liana.pl")
    fake_liana.pl.dotplot = lambda **_kwargs: _FakePlotnine()
    monkeypatch.setitem(sys.modules, "liana", fake_liana)

    with pytest.raises(DataNotFoundError, match="No suitable columns for visualization"):
        viz_cc._create_unified_dotplot(
            data,
            VisualizationParameters(plot_type="communication", subtype="dotplot"),
        )


def test_create_unified_dotplot_falls_back_on_liana_error(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(viz_cc, "require", lambda *_a, **_k: None)
    data = _mock_cc_data(_mock_liana_results())
    sentinel = object()

    fake_liana = ModuleType("liana")
    fake_liana.pl = ModuleType("liana.pl")

    def _boom(**_kwargs):
        raise RuntimeError("liana plot failure")

    fake_liana.pl.dotplot = _boom
    monkeypatch.setitem(sys.modules, "liana", fake_liana)
    monkeypatch.setattr(viz_cc, "_create_fallback_dotplot", lambda *_a, **_k: sentinel)

    out = viz_cc._create_unified_dotplot(
        data,
        VisualizationParameters(plot_type="communication", subtype="dotplot"),
    )
    assert out is sentinel


def test_create_unified_dotplot_success_and_empty_contract(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(viz_cc, "require", lambda *_a, **_k: None)

    fake_liana = ModuleType("liana")
    fake_liana.pl = ModuleType("liana.pl")
    fake_liana.pl.dotplot = lambda **_kwargs: _FakePlotnine()
    monkeypatch.setitem(sys.modules, "liana", fake_liana)

    good = _mock_cc_data(_mock_liana_results())
    fig = viz_cc._create_unified_dotplot(
        good,
        VisualizationParameters(plot_type="communication", subtype="dotplot", dpi=111),
    )
    assert fig.get_dpi() == 111
    fig.clf()

    empty = _mock_cc_data(pd.DataFrame())
    with pytest.raises(DataNotFoundError, match="No cellphonedb results found"):
        viz_cc._create_unified_dotplot(
            empty,
            VisualizationParameters(plot_type="communication", subtype="dotplot"),
        )


def test_create_fallback_dotplot_zero_rank_branch_sets_constant_size():
    df = _mock_liana_results()
    df["magnitude_rank"] = 0.0
    data = _mock_cc_data(df)

    fig = viz_cc._create_fallback_dotplot(
        data,
        VisualizationParameters(plot_type="communication", subtype="dotplot"),
    )
    sizes = fig.axes[0].collections[0].get_sizes()
    assert np.allclose(sizes, 100.0)
    fig.clf()


def test_create_fallback_dotplot_requires_rank_and_scales_positive_rank():
    with pytest.raises(DataNotFoundError, match="No ranking column found"):
        viz_cc._create_fallback_dotplot(
            _mock_cc_data(_mock_liana_results(include_scores=False)),
            VisualizationParameters(plot_type="communication", subtype="dotplot"),
        )

    fig = viz_cc._create_fallback_dotplot(
        _mock_cc_data(_mock_liana_results()),
        VisualizationParameters(plot_type="communication", subtype="dotplot"),
    )
    sizes = fig.axes[0].collections[0].get_sizes()
    assert np.max(sizes) > np.min(sizes)
    fig.clf()


def test_create_unified_tileplot_requires_value_columns(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(viz_cc, "require", lambda *_a, **_k: None)
    df = _mock_liana_results(include_scores=False)
    data = _mock_cc_data(df)

    fake_liana = ModuleType("liana")
    fake_liana.pl = ModuleType("liana.pl")
    fake_liana.pl.tileplot = lambda **_kwargs: _FakePlotnine()
    monkeypatch.setitem(sys.modules, "liana", fake_liana)

    with pytest.raises(DataNotFoundError, match="No suitable columns for tileplot"):
        viz_cc._create_unified_tileplot(
            data,
            VisualizationParameters(plot_type="communication", subtype="tileplot"),
        )


def test_create_unified_tileplot_falls_back_on_liana_error(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(viz_cc, "require", lambda *_a, **_k: None)
    data = _mock_cc_data(_mock_liana_results())
    sentinel = object()

    fake_liana = ModuleType("liana")
    fake_liana.pl = ModuleType("liana.pl")

    def _boom(**_kwargs):
        raise RuntimeError("tileplot failure")

    fake_liana.pl.tileplot = _boom
    monkeypatch.setitem(sys.modules, "liana", fake_liana)
    monkeypatch.setattr(viz_cc, "_create_fallback_tileplot", lambda *_a, **_k: sentinel)

    out = viz_cc._create_unified_tileplot(
        data,
        VisualizationParameters(plot_type="communication", subtype="tileplot"),
    )
    assert out is sentinel


def test_create_unified_tileplot_success_and_empty_contract(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(viz_cc, "require", lambda *_a, **_k: None)

    fake_liana = ModuleType("liana")
    fake_liana.pl = ModuleType("liana.pl")
    fake_liana.pl.tileplot = lambda **_kwargs: _FakePlotnine()
    monkeypatch.setitem(sys.modules, "liana", fake_liana)

    good = _mock_cc_data(_mock_liana_results())
    fig = viz_cc._create_unified_tileplot(
        good,
        VisualizationParameters(plot_type="communication", subtype="tileplot", dpi=109),
    )
    assert fig.get_dpi() == 109
    fig.clf()

    empty = _mock_cc_data(pd.DataFrame())
    with pytest.raises(DataNotFoundError, match="No cellphonedb results found"):
        viz_cc._create_unified_tileplot(
            empty,
            VisualizationParameters(plot_type="communication", subtype="tileplot"),
        )


def test_create_fallback_tileplot_requires_value_column():
    data = _mock_cc_data(_mock_liana_results(include_scores=False))
    with pytest.raises(DataNotFoundError, match="No suitable value column found"):
        viz_cc._create_fallback_tileplot(
            data,
            VisualizationParameters(plot_type="communication", subtype="tileplot"),
        )


def test_create_fallback_tileplot_success_heatmap_path():
    data = _mock_cc_data(_mock_liana_results())
    fig = viz_cc._create_fallback_tileplot(
        data,
        VisualizationParameters(
            plot_type="communication",
            subtype="tileplot",
            colormap="coolwarm",
            plot_top_pairs=2,
        ),
    )
    ax = fig.axes[0]
    assert "Cell Type Pairs" in ax.get_xlabel()
    assert "Ligand-Receptor Pairs" in ax.get_ylabel()
    fig.clf()


def test_create_unified_circle_plot_supports_weighted_and_count_modes():
    weighted = _mock_cc_data(_mock_liana_results())
    fig = viz_cc._create_unified_circle_plot(
        weighted,
        VisualizationParameters(plot_type="communication", subtype="circle_plot"),
    )
    assert fig.axes[0].name == "polar"
    fig.clf()

    counts_df = _mock_liana_results(include_scores=False)
    counts = _mock_cc_data(counts_df)
    fig2 = viz_cc._create_unified_circle_plot(
        counts,
        VisualizationParameters(plot_type="communication", subtype="circle_plot"),
    )
    assert fig2.axes[0].name == "polar"
    labels = [text.get_text() for text in fig2.axes[0].texts]
    assert labels[:2] == ["T", "B"]
    fig2.clf()

    with pytest.raises(DataNotFoundError, match="No cellphonedb results found"):
        viz_cc._create_unified_circle_plot(
            _mock_cc_data(pd.DataFrame()),
            VisualizationParameters(plot_type="communication", subtype="circle_plot"),
        )

    with pytest.raises(DataNotFoundError, match="Missing source/target columns"):
        viz_cc._create_unified_circle_plot(
            _mock_cc_data(pd.DataFrame({"ligand_complex": ["L"], "receptor_complex": ["R"]})),
            VisualizationParameters(plot_type="communication", subtype="circle_plot"),
        )


def test_plotnine_to_matplotlib_success_and_error_wrap():
    fig = viz_cc._plotnine_to_matplotlib(
        _FakePlotnine(should_fail=False),
        VisualizationParameters(plot_type="communication", dpi=123),
    )
    assert fig.get_dpi() == 123
    fig.clf()

    with pytest.raises(ProcessingError, match="Failed to convert plotnine figure"):
        viz_cc._plotnine_to_matplotlib(
            _FakePlotnine(should_fail=True),
            VisualizationParameters(plot_type="communication"),
        )


@pytest.mark.asyncio
async def test_get_cell_communication_data_uses_per_method_key(minimal_spatial_adata):
    """When method is specified, per-method uns snapshot should be used."""
    adata = minimal_spatial_adata.copy()

    # Store two methods: liana (shared/latest) and cellphonedb (older, per-method)
    cpdb_storage = CCCStorage(
        method="cellphonedb",
        analysis_type="cluster",
        species="human",
        database="cellphonedb",
        lr_pairs=["L1^R1"],
        top_lr_pairs=["L1^R1"],
        n_pairs=1,
        n_significant=1,
        results=pd.DataFrame({"T|B": [0.8]}, index=["L1_R1"]),
        pvalues=pd.DataFrame({"T|B": [0.05]}, index=["L1_R1"]),
    )
    liana_storage = CCCStorage(
        method="liana",
        analysis_type="cluster",
        species="human",
        database="consensus",
        lr_pairs=["L2^R2"],
        top_lr_pairs=["L2^R2"],
        n_pairs=1,
        n_significant=1,
        results=pd.DataFrame({"T|B": [0.5]}, index=["L2_R2"]),
        pvalues=pd.DataFrame({"T|B": [0.1]}, index=["L2_R2"]),
    )

    # Simulate dual-write: shared slot has liana (latest)
    store_ccc_results(adata, liana_storage)
    # Per-method snapshots
    adata.uns["ccc_cellphonedb"] = cpdb_storage.to_dict()
    adata.uns["ccc_liana"] = liana_storage.to_dict()

    # Without method → shared (liana)
    data_latest = await viz_cc.get_cell_communication_data(adata)
    assert data_latest.method == "liana"

    # With method → per-method (cellphonedb)
    data_cpdb = await viz_cc.get_cell_communication_data(adata, method="cellphonedb")
    assert data_cpdb.method == "cellphonedb"

    # With method → per-method (liana)
    data_liana = await viz_cc.get_cell_communication_data(adata, method="liana")
    assert data_liana.method == "liana"


@pytest.mark.asyncio
async def test_create_cell_communication_visualization_passes_analysis_method(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """create_cell_communication_visualization must forward analysis_method to get_cell_communication_data."""
    adata = minimal_spatial_adata.copy()
    captured: dict[str, object] = {}

    async def _fake_get(_adata, method=None, context=None):
        captured["method"] = method
        return _mock_cc_data(_mock_liana_results())

    monkeypatch.setattr(viz_cc, "get_cell_communication_data", _fake_get)
    monkeypatch.setattr(
        viz_cc, "_create_unified_dotplot", lambda *_a, **_k: plt.figure()
    )

    # With analysis_method set
    params = VisualizationParameters(
        plot_type="communication", analysis_method="cellphonedb"
    )
    await viz_cc.create_cell_communication_visualization(adata, params)
    assert captured["method"] == "cellphonedb"

    # Without analysis_method (None)
    params_none = VisualizationParameters(plot_type="communication")
    await viz_cc.create_cell_communication_visualization(adata, params_none)
    assert captured["method"] is None


@pytest.mark.asyncio
async def test_get_cell_communication_data_per_method_spatial_obsm(minimal_spatial_adata):
    """Per-method obsm keys should be used when method is specified."""
    adata = minimal_spatial_adata.copy()
    _store_cluster_ccc(adata)
    adata.uns["ccc"]["analysis_type"] = "spatial"
    adata.uns["ccc_cellphonedb"] = adata.uns["ccc"].copy()

    # Shared obsm: ones
    adata.obsm["ccc_spatial_scores"] = np.ones((adata.n_obs, 2))
    # Per-method obsm: twos
    adata.obsm["ccc_spatial_scores_cellphonedb"] = np.full((adata.n_obs, 2), 2.0)

    data_shared = await viz_cc.get_cell_communication_data(adata)
    assert data_shared.spatial_scores is not None
    assert data_shared.spatial_scores[0, 0] == pytest.approx(1.0)

    data_method = await viz_cc.get_cell_communication_data(
        adata, method="cellphonedb"
    )
    assert data_method.spatial_scores is not None
    assert data_method.spatial_scores[0, 0] == pytest.approx(2.0)
