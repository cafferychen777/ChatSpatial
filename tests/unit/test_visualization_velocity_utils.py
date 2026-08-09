"""Unit tests for RNA velocity visualization contracts."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import velocity as viz_vel
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


@pytest.mark.asyncio
async def test_create_rna_velocity_visualization_rejects_unsupported_subtype(
    minimal_spatial_adata,
):
    with pytest.raises(ParameterError, match="Unsupported subtype for rna_velocity"):
        await viz_vel.create_rna_velocity_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="velocity", subtype="unknown"),
        )


@pytest.mark.asyncio
async def test_create_rna_velocity_visualization_routes_stream_by_default(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    sentinel = object()

    async def _fake_stream(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_vel, "_create_velocity_stream_plot", _fake_stream)
    out = await viz_vel.create_rna_velocity_visualization(
        minimal_spatial_adata,
        VisualizationParameters(plot_type="velocity", subtype="stream"),
        context=DummyCtx(),
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_stream_requires_velocity_graph(minimal_spatial_adata, monkeypatch):
    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    with pytest.raises(DataNotFoundError, match="RNA velocity not computed"):
        await viz_vel._create_velocity_stream_plot(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="velocity", subtype="stream"),
        )


@pytest.mark.asyncio
async def test_stream_requires_valid_basis(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.uns["velocity_graph"] = np.eye(adata.n_obs)
    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    monkeypatch.setattr(viz_vel, "infer_basis", lambda *_args, **_kwargs: None)

    with pytest.raises(DataCompatibilityError, match="No valid embedding basis found"):
        await viz_vel._create_velocity_stream_plot(
            adata,
            VisualizationParameters(plot_type="velocity", subtype="stream"),
        )


@pytest.mark.asyncio
async def test_phase_requires_required_layers(minimal_spatial_adata, monkeypatch):
    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    with pytest.raises(DataNotFoundError, match="Missing layers for phase plot"):
        await viz_vel._create_velocity_phase_plot(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="velocity", subtype="phase"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_proportions_requires_velocity_layers(minimal_spatial_adata, monkeypatch):
    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    with pytest.raises(
        DataNotFoundError, match="Spliced and unspliced layers are required"
    ):
        await viz_vel._create_velocity_proportions_plot(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="velocity", subtype="proportions"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_heatmap_requires_time_or_velocity_graph(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    with pytest.raises(DataNotFoundError, match="No time ordering available"):
        await viz_vel._create_velocity_heatmap(
            adata,
            VisualizationParameters(plot_type="velocity", subtype="heatmap"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_heatmap_computes_velocity_pseudotime_when_graph_exists(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.uns["velocity_graph"] = np.eye(adata.n_obs)
    called: dict[str, bool] = {}

    monkeypatch.setattr(viz_vel, "require", _required_dependency)

    fake_scv = ModuleType("scvelo")

    def _vp(adata_obj):
        adata_obj.obs["velocity_pseudotime"] = np.linspace(0, 1, adata_obj.n_obs)
        called["vp"] = True

    def _heatmap(*_args, **_kwargs):
        plt.figure()
        return None

    fake_scv.tl = SimpleNamespace(velocity_pseudotime=_vp)
    fake_scv.pl = SimpleNamespace(heatmap=_heatmap)
    monkeypatch.setitem(__import__("sys").modules, "scvelo", fake_scv)

    fig = await viz_vel._create_velocity_heatmap(
        adata,
        VisualizationParameters(plot_type="velocity", subtype="heatmap"),
        context=DummyCtx(),
    )
    assert fig is not None
    assert called["vp"] is True


@pytest.mark.asyncio
async def test_paga_requires_cluster_key(minimal_spatial_adata):
    with pytest.raises(ParameterError, match="cluster_key is required for PAGA plot"):
        await viz_vel._create_velocity_paga_plot(
            minimal_spatial_adata,
            VisualizationParameters(
                plot_type="velocity", subtype="paga", cluster_key="missing_key"
            ),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_create_rna_velocity_visualization_routes_remaining_subtypes(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    sentinel = object()

    async def _fake(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_vel, "_create_velocity_phase_plot", _fake)
    monkeypatch.setattr(viz_vel, "_create_velocity_proportions_plot", _fake)
    monkeypatch.setattr(viz_vel, "_create_velocity_heatmap", _fake)
    monkeypatch.setattr(viz_vel, "_create_velocity_paga_plot", _fake)

    for subtype in ["phase", "proportions", "heatmap", "paga"]:
        out = await viz_vel.create_rna_velocity_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="velocity", subtype=subtype),
            context=DummyCtx(),
        )
        assert out is sentinel


@pytest.mark.asyncio
async def test_stream_success_uses_inferred_basis_and_auto_feature(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.uns["velocity_graph"] = np.eye(adata.n_obs)

    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    monkeypatch.setattr(viz_vel, "infer_basis", lambda *_a, **_k: "spatial")

    fig, ax = plt.subplots()
    monkeypatch.setattr(
        viz_vel, "create_figure_from_params", lambda *_a, **_k: (fig, [ax])
    )

    captured: dict[str, object] = {}
    fake_scv = ModuleType("scvelo")

    def _stream(*_args, **kwargs):
        captured["basis"] = kwargs.get("basis")
        captured["color"] = kwargs.get("color")
        captured["legend_loc"] = kwargs.get("legend_loc")
        kwargs["ax"].scatter([0], [0], c=[1.0])

    fake_scv.pl = SimpleNamespace(velocity_embedding_stream=_stream)
    monkeypatch.setitem(sys.modules, "scvelo", fake_scv)

    ctx = DummyCtx()
    await viz_vel._create_velocity_stream_plot(
        adata,
        VisualizationParameters(plot_type="velocity", subtype="stream", basis="umap"),
        context=ctx,
    )

    assert captured["basis"] == "spatial"
    assert captured["color"] == "group"
    assert captured["legend_loc"] == "right margin"
    assert any("Using 'spatial' as basis" in msg for msg in ctx.infos)
    assert any("Using 'group' for coloring" in msg for msg in ctx.infos)
    assert ax.yaxis_inverted()
    fig.clf()


@pytest.mark.asyncio
async def test_phase_success_uses_velocity_genes_and_context(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["velocity"] = np.zeros((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["Ms"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["Mu"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.var["velocity_genes"] = [True] * 5 + [False] * (adata.n_vars - 5)

    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    monkeypatch.setattr(viz_vel, "infer_basis", lambda *_a, **_k: "umap")

    captured: dict[str, object] = {}
    fake_scv = ModuleType("scvelo")

    def _velocity(*_args, **kwargs):
        captured["var_names"] = kwargs.get("var_names")
        captured["basis"] = kwargs.get("basis")
        captured["color"] = kwargs.get("color")
        plt.figure()

    fake_scv.pl = SimpleNamespace(velocity=_velocity)
    monkeypatch.setitem(sys.modules, "scvelo", fake_scv)

    ctx = DummyCtx()
    fig = await viz_vel._create_velocity_phase_plot(
        adata,
        VisualizationParameters(
            plot_type="velocity",
            subtype="phase",
            cluster_key="group",
            title="Phase Title",
        ),
        context=ctx,
    )
    assert captured["basis"] == "umap"
    assert captured["color"] == "group"
    assert len(captured["var_names"]) <= 4
    assert all(g in adata.var_names for g in captured["var_names"])
    assert any("Creating phase plot for genes" in msg for msg in ctx.infos)
    assert fig._suptitle is not None
    fig.clf()


@pytest.mark.asyncio
async def test_phase_supports_string_feature_and_default_genes_without_velocity_flag(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["velocity"] = np.zeros((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["Ms"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["Mu"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)

    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    monkeypatch.setattr(viz_vel, "infer_basis", lambda *_a, **_k: "umap")
    captured: dict[str, object] = {}
    fake_scv = ModuleType("scvelo")

    def _velocity(*_args, **kwargs):
        captured["var_names"] = kwargs.get("var_names")
        plt.figure()

    fake_scv.pl = SimpleNamespace(velocity=_velocity)
    monkeypatch.setitem(sys.modules, "scvelo", fake_scv)

    fig1 = await viz_vel._create_velocity_phase_plot(
        adata,
        VisualizationParameters(
            plot_type="velocity",
            subtype="phase",
            feature="gene_0",
        ),
        context=DummyCtx(),
    )
    assert captured["var_names"] == ["gene_0"]
    fig1.clf()

    fig2 = await viz_vel._create_velocity_phase_plot(
        adata,
        VisualizationParameters(plot_type="velocity", subtype="phase"),
        context=DummyCtx(),
    )
    assert captured["var_names"] == list(adata.var_names[:4])
    fig2.clf()


@pytest.mark.asyncio
async def test_phase_raises_when_requested_genes_are_missing(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["velocity"] = np.zeros((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["Ms"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["Mu"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)

    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    monkeypatch.setitem(sys.modules, "scvelo", ModuleType("scvelo"))

    with pytest.raises(DataNotFoundError, match="None of the specified genes found"):
        await viz_vel._create_velocity_phase_plot(
            adata,
            VisualizationParameters(
                plot_type="velocity",
                subtype="phase",
                feature=["missing_a", "missing_b"],
            ),
        )


@pytest.mark.asyncio
async def test_proportions_success_and_cluster_auto_selection(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)

    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    captured: dict[str, object] = {}

    fake_scv = ModuleType("scvelo")

    def _props(*_args, **kwargs):
        captured["groupby"] = kwargs.get("groupby")
        plt.figure()

    fake_scv.pl = SimpleNamespace(proportions=_props)
    monkeypatch.setitem(sys.modules, "scvelo", fake_scv)

    ctx = DummyCtx()
    fig = await viz_vel._create_velocity_proportions_plot(
        adata,
        VisualizationParameters(plot_type="velocity", subtype="proportions"),
        context=ctx,
    )
    assert captured["groupby"] == "group"
    assert any("Using cluster_key: 'group'" in msg for msg in ctx.infos)
    assert any(
        "Creating proportions plot grouped by 'group'" in msg for msg in ctx.infos
    )
    fig.clf()

    fig_titled = await viz_vel._create_velocity_proportions_plot(
        adata,
        VisualizationParameters(
            plot_type="velocity",
            subtype="proportions",
            title="Props Title",
        ),
        context=DummyCtx(),
    )
    assert fig_titled._suptitle is not None
    assert fig_titled._suptitle.get_text() == "Props Title"
    fig_titled.clf()


@pytest.mark.asyncio
async def test_proportions_requires_cluster_key_when_no_categorical(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.layers["unspliced"] = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    adata.obs = pd.DataFrame(index=adata.obs.index)

    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    monkeypatch.setitem(sys.modules, "scvelo", ModuleType("scvelo"))

    with pytest.raises(
        ParameterError, match="cluster_key is required for proportions plot"
    ):
        await viz_vel._create_velocity_proportions_plot(
            adata,
            VisualizationParameters(plot_type="velocity", subtype="proportions"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_heatmap_success_uses_latent_time_and_hvg_fallback(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["latent_time"] = np.linspace(0.0, 1.0, adata.n_obs)
    adata.var["highly_variable"] = [True] * 8 + [False] * (adata.n_vars - 8)

    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    captured: dict[str, object] = {}
    fake_scv = ModuleType("scvelo")

    def _heatmap(*_args, **kwargs):
        captured["sortby"] = kwargs.get("sortby")
        captured["var_names"] = kwargs.get("var_names")
        captured["col_color"] = kwargs.get("col_color")
        plt.figure()

    fake_scv.pl = SimpleNamespace(heatmap=_heatmap)
    fake_scv.tl = SimpleNamespace(velocity_pseudotime=lambda *_args, **_kwargs: None)
    monkeypatch.setitem(sys.modules, "scvelo", fake_scv)

    ctx = DummyCtx()
    fig = await viz_vel._create_velocity_heatmap(
        adata,
        VisualizationParameters(
            plot_type="velocity",
            subtype="heatmap",
            cluster_key="group",
            title="Heatmap Title",
        ),
        context=ctx,
    )
    assert captured["sortby"] == "latent_time"
    assert len(captured["var_names"]) == 8
    assert captured["col_color"] == "group"
    assert any("Using 'latent_time' for heatmap ordering" in msg for msg in ctx.infos)
    assert any("Creating velocity heatmap with 8 genes" in msg for msg in ctx.infos)
    assert fig._suptitle is not None
    fig.clf()


@pytest.mark.asyncio
async def test_heatmap_feature_string_and_velocity_genes_branches(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["latent_time"] = np.linspace(0.0, 1.0, adata.n_obs)
    adata.var["velocity_genes"] = [True] * 6 + [False] * (adata.n_vars - 6)

    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    captured: dict[str, object] = {}
    fake_scv = ModuleType("scvelo")

    def _heatmap(*_args, **kwargs):
        captured["var_names"] = kwargs.get("var_names")
        plt.figure()

    fake_scv.pl = SimpleNamespace(heatmap=_heatmap)
    fake_scv.tl = SimpleNamespace(velocity_pseudotime=lambda *_args, **_kwargs: None)
    monkeypatch.setitem(sys.modules, "scvelo", fake_scv)

    fig_feature = await viz_vel._create_velocity_heatmap(
        adata,
        VisualizationParameters(
            plot_type="velocity",
            subtype="heatmap",
            feature="gene_0",
        ),
        context=DummyCtx(),
    )
    assert captured["var_names"] == ["gene_0"]
    fig_feature.clf()

    fig_default = await viz_vel._create_velocity_heatmap(
        adata,
        VisualizationParameters(plot_type="velocity", subtype="heatmap"),
        context=DummyCtx(),
    )
    assert captured["var_names"] == list(adata.var_names[:6])
    fig_default.clf()


@pytest.mark.asyncio
async def test_heatmap_rejects_missing_requested_genes(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["latent_time"] = np.linspace(0.0, 1.0, adata.n_obs)

    monkeypatch.setattr(viz_vel, "require", _required_dependency)
    fake_scv = ModuleType("scvelo")
    fake_scv.pl = SimpleNamespace(heatmap=lambda *_args, **_kwargs: None)
    fake_scv.tl = SimpleNamespace(velocity_pseudotime=lambda *_args, **_kwargs: None)
    monkeypatch.setitem(sys.modules, "scvelo", fake_scv)

    with pytest.raises(DataNotFoundError, match="None of the specified genes found"):
        await viz_vel._create_velocity_heatmap(
            adata,
            VisualizationParameters(
                plot_type="velocity",
                subtype="heatmap",
                feature=["missing_gene"],
            ),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_paga_success_recompute_and_uns_group_shortcut(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    calls: dict[str, int] = {"paga": 0, "plot": 0}

    import scanpy as sc

    def _paga(*_args, **_kwargs):
        calls["paga"] += 1

    def _plot_paga(*_args, **kwargs):
        calls["plot"] += 1
        kwargs["ax"].scatter([0], [0], c=[1.0])

    monkeypatch.setattr(sc.tl, "paga", _paga)
    monkeypatch.setattr(sc.pl, "paga", _plot_paga)

    fig, ax = plt.subplots()
    monkeypatch.setattr(
        viz_vel, "create_figure_from_params", lambda *_a, **_k: (fig, [ax])
    )

    ctx = DummyCtx()
    out = await viz_vel._create_velocity_paga_plot(
        adata,
        VisualizationParameters(
            plot_type="velocity", subtype="paga", cluster_key="group"
        ),
        context=ctx,
    )
    assert out is fig
    assert calls["paga"] == 1
    assert calls["plot"] == 1
    assert any("Computing PAGA for cluster_key='group'" in msg for msg in ctx.infos)
    assert any("Creating PAGA plot for 'group'" in msg for msg in ctx.infos)
    fig.clf()

    adata.uns["paga"] = {"groups": "group"}
    fig2, ax2 = plt.subplots()
    monkeypatch.setattr(
        viz_vel, "create_figure_from_params", lambda *_a, **_k: (fig2, [ax2])
    )
    calls["paga"] = 0
    out2 = await viz_vel._create_velocity_paga_plot(
        adata,
        VisualizationParameters(
            plot_type="velocity", subtype="paga", title="Custom PAGA"
        ),
        context=DummyCtx(),
    )
    assert out2 is fig2
    assert calls["paga"] == 0
    assert ax2.get_title() == "Custom PAGA"
    fig2.clf()


@pytest.mark.asyncio
async def test_paga_auto_selects_first_categorical_cluster_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    calls: dict[str, object] = {"groups": None}
    import scanpy as sc

    def _paga(_adata, groups=None):
        calls["groups"] = groups

    def _plot_paga(*_args, **kwargs):
        kwargs["ax"].scatter([0], [0], c=[1.0])

    monkeypatch.setattr(sc.tl, "paga", _paga)
    monkeypatch.setattr(sc.pl, "paga", _plot_paga)

    fig, ax = plt.subplots()
    monkeypatch.setattr(
        viz_vel, "create_figure_from_params", lambda *_a, **_k: (fig, [ax])
    )

    out = await viz_vel._create_velocity_paga_plot(
        adata,
        VisualizationParameters(plot_type="velocity", subtype="paga"),
        context=DummyCtx(),
    )
    assert out is fig
    assert calls["groups"] == "group"
    fig.clf()
