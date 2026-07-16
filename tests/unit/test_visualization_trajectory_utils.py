"""Unit tests for trajectory visualization routing and guard contracts."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools import trajectory as traj_tool
from chatspatial.tools.visualization import trajectory as viz_traj
from chatspatial.utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    ParameterError,
)


def _required_dependency(name: str, *_args, **_kwargs):
    return sys.modules.get(name, object())


class DummyCtx:
    def __init__(self):
        self.infos: list[str] = []

    async def info(self, msg: str):
        self.infos.append(msg)


def _install_fake_cellrank(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_cr = type("FakeCellRank", (), {})()
    fake_cr.pl = type("FakeCRPlt", (), {})()
    monkeypatch.setitem(__import__("sys").modules, "cellrank", fake_cr)


@pytest.mark.asyncio
async def test_create_trajectory_visualization_rejects_unsupported_subtype(
    minimal_spatial_adata,
):
    with pytest.raises(ParameterError, match="Unsupported subtype for trajectory"):
        await viz_traj.create_trajectory_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="trajectory", subtype="bad"),
        )


@pytest.mark.asyncio
async def test_create_trajectory_visualization_routes_pseudotime(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    sentinel = object()

    async def _fake(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_traj, "_create_trajectory_pseudotime_plot", _fake)
    out = await viz_traj.create_trajectory_visualization(
        minimal_spatial_adata,
        VisualizationParameters(plot_type="trajectory", subtype="pseudotime"),
        context=DummyCtx(),
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_pseudotime_plot_requires_pseudotime_column(minimal_spatial_adata):
    with pytest.raises(DataNotFoundError, match="No pseudotime found"):
        await viz_traj._create_trajectory_pseudotime_plot(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="trajectory", subtype="pseudotime"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_pseudotime_plot_requires_valid_basis(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["dpt_pseudotime"] = np.linspace(0, 1, adata.n_obs)
    monkeypatch.setattr(viz_traj, "infer_basis", lambda *_args, **_kwargs: None)
    with pytest.raises(DataCompatibilityError, match="No valid embedding basis found"):
        await viz_traj._create_trajectory_pseudotime_plot(
            adata,
            VisualizationParameters(plot_type="trajectory", subtype="pseudotime"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_cellrank_circular_requires_fate_probabilities(minimal_spatial_adata, monkeypatch):
    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    _install_fake_cellrank(monkeypatch)
    with pytest.raises(DataNotFoundError, match="CellRank fate probabilities not found"):
        await viz_traj._create_cellrank_circular_projection(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="trajectory", subtype="circular"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_cellrank_fate_map_requires_cluster_key_when_no_categorical(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2))
    # remove categorical fallback
    adata.obs = pd.DataFrame(index=adata.obs.index)
    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    _install_fake_cellrank(monkeypatch)

    with pytest.raises(ParameterError, match="cluster_key is required for fate map"):
        await viz_traj._create_cellrank_fate_map(
            adata,
            VisualizationParameters(plot_type="trajectory", subtype="fate_map"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_cellrank_gene_trends_requires_time_key(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2))
    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    _install_fake_cellrank(monkeypatch)

    with pytest.raises(DataNotFoundError, match="No pseudotime found"):
        await viz_traj._create_cellrank_gene_trends(
            adata,
            VisualizationParameters(plot_type="trajectory", subtype="gene_trends"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_palantir_results_requires_pseudotime(minimal_spatial_adata):
    with pytest.raises(DataNotFoundError, match="Palantir results not found"):
        await viz_traj._create_palantir_results(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="trajectory", subtype="palantir"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_palantir_results_nonempty_fate_probs_branch(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["palantir_pseudotime"] = np.linspace(0, 1, adata.n_obs)
    adata.obs["palantir_entropy"] = np.linspace(1, 0, adata.n_obs)
    adata.obsm["palantir_fate_probs"] = pd.DataFrame(
        {"f1": np.linspace(0.6, 0.3, adata.n_obs), "f2": np.linspace(0.4, 0.7, adata.n_obs)},
        index=adata.obs_names,
    )

    monkeypatch.setattr(viz_traj, "infer_basis", lambda *_args, **_kwargs: "spatial")

    calls: list[str] = []

    def _embedding(*_args, **kwargs):
        calls.append(str(kwargs.get("color")))
        return None

    monkeypatch.setattr(viz_traj.sc.pl, "embedding", _embedding)

    fig = await viz_traj._create_palantir_results(
        adata,
        VisualizationParameters(plot_type="trajectory", subtype="palantir"),
        context=DummyCtx(),
    )
    assert fig is not None
    assert "palantir_pseudotime" in calls
    assert "palantir_entropy" in calls
    assert "_dominant_fate" in calls
    assert "_dominant_fate" not in adata.obs.columns
    plt.close("all")


@pytest.mark.asyncio
async def test_create_trajectory_visualization_routes_all_supported_subtypes(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    sentinel = object()

    async def _fake(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_traj, "_create_cellrank_circular_projection", _fake)
    monkeypatch.setattr(viz_traj, "_create_cellrank_fate_map", _fake)
    monkeypatch.setattr(viz_traj, "_create_cellrank_gene_trends", _fake)
    monkeypatch.setattr(viz_traj, "_create_cellrank_fate_heatmap", _fake)
    monkeypatch.setattr(viz_traj, "_create_palantir_results", _fake)

    for subtype in ["circular", "fate_map", "gene_trends", "fate_heatmap", "palantir"]:
        out = await viz_traj.create_trajectory_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="trajectory", subtype=subtype),
            context=DummyCtx(),
        )
        assert out is sentinel


@pytest.mark.asyncio
async def test_pseudotime_success_with_velocity_stream_panel(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["dpt_pseudotime"] = np.linspace(0, 1, adata.n_obs)
    adata.uns["velocity_graph"] = np.eye(adata.n_obs)

    monkeypatch.setattr(viz_traj, "infer_basis", lambda *_a, **_k: "spatial")

    def _embedding(*_args, **kwargs):
        kwargs["ax"].scatter([0], [0], c=[1.0])

    monkeypatch.setattr(viz_traj.sc.pl, "embedding", _embedding)

    fake_scv = ModuleType("scvelo")

    def _stream(*_args, **kwargs):
        kwargs["ax"].scatter([0], [0], c=[1.0])

    fake_scv.pl = SimpleNamespace(velocity_embedding_stream=_stream)
    monkeypatch.setattr(
        viz_traj,
        "require",
        lambda dependency, *_args, **_kwargs: (
            fake_scv
            if dependency == "scvelo"
            else pytest.fail(f"Unexpected dependency: {dependency}")
        ),
    )

    ctx = DummyCtx()
    fig = await viz_traj._create_trajectory_pseudotime_plot(
        adata,
        VisualizationParameters(plot_type="trajectory", subtype="pseudotime"),
        context=ctx,
    )
    assert len(fig.axes) >= 2
    assert any("Found pseudotime column: dpt_pseudotime" in msg for msg in ctx.infos)
    for ax in fig.axes[:2]:
        assert ax.yaxis_inverted()
    fig.clf()


@pytest.mark.asyncio
async def test_cellrank_circular_projection_success(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)

    monkeypatch.setattr(viz_traj, "require", _required_dependency)

    fake_cr = ModuleType("cellrank")

    def _circular(*_args, **_kwargs):
        plt.figure()

    fake_cr.pl = SimpleNamespace(circular_projection=_circular)
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    fig = await viz_traj._create_cellrank_circular_projection(
        adata,
        VisualizationParameters(
            plot_type="trajectory",
            subtype="circular",
            title="Circular Title",
        ),
        context=DummyCtx(),
    )
    assert fig._suptitle is not None
    fig.clf()


@pytest.mark.asyncio
async def test_cellrank_fate_map_success_with_auto_cluster_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)

    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    captured: dict[str, object] = {}

    fake_cr = ModuleType("cellrank")

    def _aggregate(*_args, **kwargs):
        captured.update(kwargs)
        plt.figure()

    fake_cr.pl = SimpleNamespace(aggregate_fate_probabilities=_aggregate)
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    ctx = DummyCtx()
    fig = await viz_traj._create_cellrank_fate_map(
        adata,
        VisualizationParameters(plot_type="trajectory", subtype="fate_map"),
        context=ctx,
    )
    assert captured["cluster_key"] == "group"
    assert captured["mode"] == "heatmap"
    assert any("Using cluster_key: 'group'" in msg for msg in ctx.infos)
    assert any("Creating CellRank fate map for 'group'" in msg for msg in ctx.infos)
    fig.clf()


@pytest.mark.asyncio
async def test_cellrank_gene_trends_success_with_filtered_features(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)
    adata.obs["latent_time"] = np.linspace(0, 1, adata.n_obs)

    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    monkeypatch.setattr(
        traj_tool,
        "prepare_gam_model_for_visualization",
        lambda *_a, **_k: ("fake_model", ["L1", "L2"]),
    )

    captured: dict[str, object] = {}
    fake_cr = ModuleType("cellrank")

    def _gene_trends(*_args, **kwargs):
        captured["genes"] = kwargs.get("genes")
        captured["time_key"] = kwargs.get("time_key")
        plt.figure()

    fake_cr.pl = SimpleNamespace(gene_trends=_gene_trends)
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    ctx = DummyCtx()
    fig = await viz_traj._create_cellrank_gene_trends(
        adata,
        VisualizationParameters(
            plot_type="trajectory",
            subtype="gene_trends",
            feature=["gene_0", "missing_gene", "gene_1"],
        ),
        context=ctx,
    )
    assert captured["genes"] == ["gene_0", "gene_1"]
    assert captured["time_key"] == "latent_time"
    assert any("Creating gene trends for: ['gene_0', 'gene_1']" in msg for msg in ctx.infos)
    assert any("Lineages: ['L1', 'L2']" in msg for msg in ctx.infos)
    fig.clf()


@pytest.mark.asyncio
async def test_cellrank_fate_heatmap_success_and_missing_genes_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)
    adata.obs["latent_time"] = np.linspace(0, 1, adata.n_obs)

    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    monkeypatch.setattr(
        traj_tool,
        "prepare_gam_model_for_visualization",
        lambda *_a, **_k: ("fake_model", ["L1"]),
    )

    captured: dict[str, object] = {}
    fake_cr = ModuleType("cellrank")

    def _heatmap(*_args, **kwargs):
        captured["genes"] = kwargs.get("genes")
        captured["time_key"] = kwargs.get("time_key")
        plt.figure()

    fake_cr.pl = SimpleNamespace(heatmap=_heatmap)
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    ctx = DummyCtx()
    fig = await viz_traj._create_cellrank_fate_heatmap(
        adata,
        VisualizationParameters(
            plot_type="trajectory",
            subtype="fate_heatmap",
            feature=["gene_0", "gene_1"],
        ),
        context=ctx,
    )
    assert captured["genes"] == ["gene_0", "gene_1"]
    assert captured["time_key"] == "latent_time"
    assert any("Creating fate heatmap with 2 genes" in msg for msg in ctx.infos)
    assert any("Lineages: ['L1']" in msg for msg in ctx.infos)
    fig.clf()

    with pytest.raises(DataNotFoundError, match="None of the genes found"):
        await viz_traj._create_cellrank_fate_heatmap(
            adata,
            VisualizationParameters(
                plot_type="trajectory",
                subtype="fate_heatmap",
                feature=["missing_a", "missing_b"],
            ),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_pseudotime_plot_feature_list_single_panel_and_title(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["pt"] = np.linspace(0, 1, adata.n_obs)
    monkeypatch.setattr(viz_traj, "infer_basis", lambda *_a, **_k: "umap")

    def _embedding(*_args, **kwargs):
        kwargs["ax"].scatter([0], [0], c=[1.0])

    monkeypatch.setattr(viz_traj.sc.pl, "embedding", _embedding)

    fig = await viz_traj._create_trajectory_pseudotime_plot(
        adata,
        VisualizationParameters(
            plot_type="trajectory",
            subtype="pseudotime",
            feature=["pt"],
            show_colorbar=False,
            title="PT Title",
        ),
        context=DummyCtx(),
    )
    assert fig._suptitle is not None
    assert fig._suptitle.get_text() == "PT Title"
    fig.clf()


@pytest.mark.asyncio
async def test_cellrank_fate_map_missing_fate_and_title_branch(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    _install_fake_cellrank(monkeypatch)

    with pytest.raises(DataNotFoundError, match="CellRank fate probabilities not found"):
        await viz_traj._create_cellrank_fate_map(
            adata,
            VisualizationParameters(plot_type="trajectory", subtype="fate_map"),
            context=DummyCtx(),
        )

    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)
    fake_cr = ModuleType("cellrank")
    fake_cr.pl = SimpleNamespace(
        aggregate_fate_probabilities=lambda *_a, **_k: plt.figure()
    )
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    fig = await viz_traj._create_cellrank_fate_map(
        adata,
        VisualizationParameters(
            plot_type="trajectory",
            subtype="fate_map",
            cluster_key="group",
            title="Fate Title",
        ),
        context=DummyCtx(),
    )
    assert fig._suptitle is not None
    assert fig._suptitle.get_text() == "Fate Title"
    fig.clf()


@pytest.mark.asyncio
async def test_cellrank_gene_trends_edge_branches(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["latent_time"] = np.linspace(0, 1, adata.n_obs)
    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    _install_fake_cellrank(monkeypatch)

    with pytest.raises(DataNotFoundError, match="fate probabilities not found"):
        await viz_traj._create_cellrank_gene_trends(
            adata,
            VisualizationParameters(plot_type="trajectory", subtype="gene_trends"),
            context=DummyCtx(),
        )

    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)
    with pytest.raises(DataNotFoundError, match="None of the specified genes found"):
        await viz_traj._create_cellrank_gene_trends(
            adata,
            VisualizationParameters(
                plot_type="trajectory",
                subtype="gene_trends",
                feature="missing_gene",
            ),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_cellrank_gene_trends_defaults_and_title(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)
    adata.obs["latent_time"] = np.linspace(0, 1, adata.n_obs)
    adata.var["highly_variable"] = [i < 4 for i in range(adata.n_vars)]

    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    monkeypatch.setattr(
        traj_tool,
        "prepare_gam_model_for_visualization",
        lambda *_a, **_k: ("fake_model", ["L1"]),
    )

    fake_cr = ModuleType("cellrank")
    captured: dict[str, object] = {}

    def _gene_trends(*_args, **kwargs):
        captured["genes"] = kwargs.get("genes")
        plt.figure()

    fake_cr.pl = SimpleNamespace(gene_trends=_gene_trends)
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    fig = await viz_traj._create_cellrank_gene_trends(
        adata,
        VisualizationParameters(
            plot_type="trajectory",
            subtype="gene_trends",
            title="GT Title",
        ),
        context=DummyCtx(),
    )
    assert captured["genes"] == ["gene_0", "gene_1", "gene_2", "gene_3"]
    assert fig._suptitle is not None
    assert fig._suptitle.get_text() == "GT Title"
    fig.clf()


@pytest.mark.asyncio
async def test_cellrank_fate_heatmap_edge_branches(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    _install_fake_cellrank(monkeypatch)

    with pytest.raises(DataNotFoundError, match="fate probabilities not found"):
        await viz_traj._create_cellrank_fate_heatmap(
            adata,
            VisualizationParameters(plot_type="trajectory", subtype="fate_heatmap"),
            context=DummyCtx(),
        )

    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)
    with pytest.raises(DataNotFoundError, match="No pseudotime found for fate heatmap"):
        await viz_traj._create_cellrank_fate_heatmap(
            adata,
            VisualizationParameters(plot_type="trajectory", subtype="fate_heatmap"),
            context=DummyCtx(),
        )

    adata.obs["latent_time"] = np.linspace(0, 1, adata.n_obs)
    monkeypatch.setattr(
        traj_tool,
        "prepare_gam_model_for_visualization",
        lambda *_a, **_k: ("fake_model", ["L1"]),
    )
    fake_cr = ModuleType("cellrank")
    captured: dict[str, object] = {}

    def _heatmap(*_args, **kwargs):
        captured["genes"] = kwargs.get("genes")
        plt.figure()

    fake_cr.pl = SimpleNamespace(heatmap=_heatmap)
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    fig = await viz_traj._create_cellrank_fate_heatmap(
        adata,
        VisualizationParameters(
            plot_type="trajectory",
            subtype="fate_heatmap",
            feature="gene_0",
            title="FH Title",
        ),
        context=DummyCtx(),
    )
    assert captured["genes"] == ["gene_0"]
    assert fig._suptitle is not None
    assert fig._suptitle.get_text() == "FH Title"
    fig.clf()


@pytest.mark.asyncio
async def test_cellrank_fate_heatmap_defaults_without_hvg_uses_var_names(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)
    adata.obs["latent_time"] = np.linspace(0, 1, adata.n_obs)

    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    monkeypatch.setattr(
        traj_tool,
        "prepare_gam_model_for_visualization",
        lambda *_a, **_k: ("fake_model", ["L1"]),
    )

    fake_cr = ModuleType("cellrank")
    captured: dict[str, object] = {}

    def _heatmap(*_args, **kwargs):
        captured["genes"] = kwargs.get("genes")
        plt.figure()

    fake_cr.pl = SimpleNamespace(heatmap=_heatmap)
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    fig = await viz_traj._create_cellrank_fate_heatmap(
        adata,
        VisualizationParameters(plot_type="trajectory", subtype="fate_heatmap"),
        context=DummyCtx(),
    )
    assert captured["genes"] == list(adata.var_names[:50])
    fig.clf()


@pytest.mark.asyncio
async def test_palantir_single_panel_and_ndarray_fates(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["palantir_pseudotime"] = np.linspace(0, 1, adata.n_obs)
    monkeypatch.setattr(viz_traj, "infer_basis", lambda *_a, **_k: "umap")
    calls: list[str] = []

    def _embedding(*_args, **kwargs):
        calls.append(str(kwargs.get("color")))
        kwargs["ax"].scatter([0], [0], c=[1.0])

    monkeypatch.setattr(viz_traj.sc.pl, "embedding", _embedding)

    fig_single = await viz_traj._create_palantir_results(
        adata,
        VisualizationParameters(plot_type="trajectory", subtype="palantir"),
        context=DummyCtx(),
    )
    assert len(fig_single.axes) == 1
    fig_single.clf()

    adata.obsm["palantir_branch_probs"] = np.column_stack(
        [np.linspace(0.6, 0.2, adata.n_obs), np.linspace(0.4, 0.8, adata.n_obs)]
    )
    fig_nd = await viz_traj._create_palantir_results(
        adata,
        VisualizationParameters(plot_type="trajectory", subtype="palantir"),
        context=DummyCtx(),
    )
    assert "_dominant_fate" in calls
    assert "_dominant_fate" not in adata.obs.columns
    fig_nd.clf()


@pytest.mark.asyncio
async def test_cellrank_gene_trends_default_without_hvg_uses_var_names(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)
    adata.obs["latent_time"] = np.linspace(0, 1, adata.n_obs)
    if "highly_variable" in adata.var.columns:
        del adata.var["highly_variable"]

    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    monkeypatch.setattr(
        traj_tool,
        "prepare_gam_model_for_visualization",
        lambda *_a, **_k: ("fake_model", ["L1"]),
    )

    fake_cr = ModuleType("cellrank")
    captured: dict[str, object] = {}

    def _gene_trends(*_args, **kwargs):
        captured["genes"] = kwargs.get("genes")
        plt.figure()

    fake_cr.pl = SimpleNamespace(gene_trends=_gene_trends)
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    fig = await viz_traj._create_cellrank_gene_trends(
        adata,
        VisualizationParameters(plot_type="trajectory", subtype="gene_trends"),
        context=DummyCtx(),
    )
    assert captured["genes"] == list(adata.var_names[:6])
    fig.clf()


@pytest.mark.asyncio
async def test_cellrank_fate_heatmap_defaults_with_hvg_uses_hvg_genes(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)
    adata.obs["latent_time"] = np.linspace(0, 1, adata.n_obs)
    adata.var["highly_variable"] = [i < 3 for i in range(adata.n_vars)]

    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    monkeypatch.setattr(
        traj_tool,
        "prepare_gam_model_for_visualization",
        lambda *_a, **_k: ("fake_model", ["L1"]),
    )

    fake_cr = ModuleType("cellrank")
    captured: dict[str, object] = {}

    def _heatmap(*_args, **kwargs):
        captured["genes"] = kwargs.get("genes")
        plt.figure()

    fake_cr.pl = SimpleNamespace(heatmap=_heatmap)
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    fig = await viz_traj._create_cellrank_fate_heatmap(
        adata,
        VisualizationParameters(plot_type="trajectory", subtype="fate_heatmap"),
        context=DummyCtx(),
    )
    assert captured["genes"] == ["gene_0", "gene_1", "gene_2"]
    fig.clf()


# ---------------------------------------------------------------------------
# Regression: fate_probabilities alias → CellRank-recognised key
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_circular_projection_finds_fate_probabilities_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """When only 'fate_probabilities' is present, _resolve_cellrank_fate_key
    must find it without modifying adata.obsm (read-only)."""
    adata = minimal_spatial_adata.copy()
    adata.obsm["fate_probabilities"] = np.ones((adata.n_obs, 2), dtype=float)
    obsm_keys_before = set(adata.obsm.keys())

    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    fake_cr = ModuleType("cellrank")
    fake_cr.pl = SimpleNamespace(circular_projection=lambda *_a, **_k: plt.figure())
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    fig = await viz_traj._create_cellrank_circular_projection(
        adata,
        VisualizationParameters(plot_type="trajectory", subtype="circular"),
        context=DummyCtx(),
    )
    # _resolve_cellrank_fate_key should NOT write new keys
    assert set(adata.obsm.keys()) == obsm_keys_before
    fig.clf()


@pytest.mark.asyncio
async def test_fate_map_finds_fate_probabilities_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """fate_map with only 'fate_probabilities' must find it read-only."""
    adata = minimal_spatial_adata.copy()
    adata.obsm["fate_probabilities"] = np.ones((adata.n_obs, 2), dtype=float)
    obsm_keys_before = set(adata.obsm.keys())

    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    fake_cr = ModuleType("cellrank")
    fake_cr.pl = SimpleNamespace(
        aggregate_fate_probabilities=lambda *_a, **_k: plt.figure()
    )
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    fig = await viz_traj._create_cellrank_fate_map(
        adata,
        VisualizationParameters(
            plot_type="trajectory", subtype="fate_map", cluster_key="group"
        ),
        context=DummyCtx(),
    )
    # _resolve_cellrank_fate_key should NOT write new keys
    assert set(adata.obsm.keys()) == obsm_keys_before
    fig.clf()


# ---------------------------------------------------------------------------
# _resolve_cellrank_fate_key is read-only (no side-effects)
# ---------------------------------------------------------------------------


def test_resolve_cellrank_fate_key_is_read_only(minimal_spatial_adata):
    """_resolve_cellrank_fate_key must not modify adata.obsm."""
    adata = minimal_spatial_adata.copy()
    adata.obsm["fate_probabilities"] = np.ones((adata.n_obs, 2), dtype=float)
    obsm_keys_before = set(adata.obsm.keys())

    key = viz_traj._resolve_cellrank_fate_key(adata)
    assert key == "fate_probabilities"
    assert set(adata.obsm.keys()) == obsm_keys_before


def test_resolve_cellrank_fate_key_prefers_lineages_fwd(minimal_spatial_adata):
    """Priority order: lineages_fwd > to_terminal_states > fate_probabilities."""
    adata = minimal_spatial_adata.copy()
    adata.obsm["fate_probabilities"] = np.ones((adata.n_obs, 2), dtype=float)
    adata.obsm["to_terminal_states"] = np.ones((adata.n_obs, 2), dtype=float)
    adata.obsm["lineages_fwd"] = np.ones((adata.n_obs, 2), dtype=float)

    key = viz_traj._resolve_cellrank_fate_key(adata)
    assert key == "lineages_fwd"


# ---------------------------------------------------------------------------
# Regression: bbknn empty results_keys still records in _index.json
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_export_empty_results_keys_still_writes_index(
    minimal_spatial_adata, tmp_path, monkeypatch
):
    """bbknn-like analysis with empty results_keys must appear in _index.json."""
    import json

    from chatspatial.utils import results_export

    adata = minimal_spatial_adata.copy()
    # Simulate metadata written by integration (bbknn)
    adata.uns["integration_bbknn_metadata"] = {
        "method": "bbknn",
        "parameters": {"batch_key": "batch"},
        "results_keys": {},
        "statistics": {"n_batches": 2},
    }

    monkeypatch.setattr(results_export, "_is_export_enabled", lambda: True)
    monkeypatch.setattr(
        results_export,
        "get_results_dir",
        lambda data_id: tmp_path,
    )

    files = results_export.export_analysis_result(adata, "d1", "integration_bbknn")
    assert files == []  # no CSVs exported

    index_path = tmp_path / "_index.json"
    assert index_path.exists()
    index = json.loads(index_path.read_text())
    assert "integration_bbknn" in index["analyses"]
    assert index["analyses"]["integration_bbknn"]["method"] == "bbknn"


# ---------------------------------------------------------------------------
# Regression: fate_map must use a PAGA-independent mode
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_fate_map_uses_paga_independent_heatmap_mode(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["to_terminal_states"] = np.ones((adata.n_obs, 2), dtype=float)

    monkeypatch.setattr(viz_traj, "require", _required_dependency)
    captured: dict[str, object] = {}

    fake_cr = ModuleType("cellrank")

    def _aggregate(*_args, **kwargs):
        captured.update(kwargs)
        plt.figure()

    fake_cr.pl = SimpleNamespace(aggregate_fate_probabilities=_aggregate)
    monkeypatch.setitem(sys.modules, "cellrank", fake_cr)

    await viz_traj._create_cellrank_fate_map(
        adata,
        VisualizationParameters(
            plot_type="trajectory", subtype="fate_map", cluster_key="group"
        ),
        context=DummyCtx(),
    )
    assert captured["mode"] == "heatmap"
    plt.close("all")
