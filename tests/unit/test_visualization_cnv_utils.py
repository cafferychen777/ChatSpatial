"""Unit tests for CNV visualization routing and branch contracts."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import cnv as viz_cnv
from chatspatial.utils.exceptions import DataNotFoundError, ParameterError


class DummyCtx:
    def __init__(self):
        self.infos: list[str] = []
        self.warnings: list[str] = []

    async def info(self, msg: str):
        self.infos.append(msg)

    async def warning(self, msg: str):
        self.warnings.append(msg)


@pytest.mark.asyncio
async def test_create_cnv_visualization_routes_and_validates_subtype(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    sentinel = object()

    async def _fake(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_cnv, "_create_cnv_heatmap", _fake)
    monkeypatch.setattr(viz_cnv, "_create_spatial_cnv", _fake)

    out_heatmap = await viz_cnv.create_cnv_visualization(
        minimal_spatial_adata,
        VisualizationParameters(plot_type="cnv", subtype="heatmap"),
    )
    out_spatial = await viz_cnv.create_cnv_visualization(
        minimal_spatial_adata,
        VisualizationParameters(plot_type="cnv", subtype="spatial"),
    )
    assert out_heatmap is sentinel
    assert out_spatial is sentinel

    with pytest.raises(ParameterError, match="Unknown CNV visualization subtype"):
        await viz_cnv.create_cnv_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="cnv", subtype="unknown"),
        )


@pytest.mark.asyncio
async def test_spatial_cnv_requires_cnv_features(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(viz_cnv, "require_spatial_coords", lambda _a: _a.obsm["spatial"])
    ctx = DummyCtx()

    with pytest.raises(DataNotFoundError, match="No CNV features found"):
        await viz_cnv._create_spatial_cnv(
            adata,
            VisualizationParameters(plot_type="cnv", subtype="spatial"),
            context=ctx,
        )
    assert any("No CNV features found" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_spatial_cnv_auto_detect_priority_and_categorical_colormap(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["numbat_clone"] = pd.Categorical(["C1"] * adata.n_obs)
    adata.obs["cnv_score"] = np.linspace(0.0, 1.0, adata.n_obs)
    adata.obs["numbat_p_cnv"] = np.linspace(0.0, 1.0, adata.n_obs)

    monkeypatch.setattr(viz_cnv, "require_spatial_coords", lambda _a: _a.obsm["spatial"])
    fig, ax = plt.subplots()
    monkeypatch.setattr(viz_cnv, "create_figure_from_params", lambda *_a, **_k: (fig, [ax]))

    captured: dict[str, Any] = {}

    def _plot_spatial(_adata, _ax, feature, params):
        captured["feature"] = feature
        captured["colormap"] = params.colormap
        _ax.scatter([0], [0], c=[1.0])

    monkeypatch.setattr(viz_cnv, "plot_spatial_feature", _plot_spatial)

    params = VisualizationParameters(plot_type="cnv", subtype="spatial", colormap="")
    ctx = DummyCtx()
    out = await viz_cnv._create_spatial_cnv(adata, params, context=ctx)

    assert out is fig
    assert captured["feature"] == "numbat_clone"
    # Colormap should be applied to the local copy, not mutate the original
    assert captured["colormap"] == "tab20"
    assert params.colormap == ""  # original params unchanged
    assert any("using 'numbat_clone'" in msg for msg in ctx.infos)
    fig.clf()


@pytest.mark.asyncio
async def test_spatial_cnv_numeric_feature_sets_rdbu(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["cnv_score"] = np.linspace(0.0, 1.0, adata.n_obs)

    monkeypatch.setattr(viz_cnv, "require_spatial_coords", lambda _a: _a.obsm["spatial"])
    fig, ax = plt.subplots()
    monkeypatch.setattr(viz_cnv, "create_figure_from_params", lambda *_a, **_k: (fig, [ax]))
    captured_cmap: dict[str, str] = {}
    monkeypatch.setattr(
        viz_cnv,
        "plot_spatial_feature",
        lambda _adata, _ax, feature, params: (
            captured_cmap.__setitem__("colormap", params.colormap),
            _ax.scatter([0], [0], c=[1.0]),
        ),
    )

    params = VisualizationParameters(plot_type="cnv", subtype="spatial", colormap="")
    await viz_cnv._create_spatial_cnv(adata, params, context=DummyCtx())
    # Colormap applied to local copy, not mutating original
    assert captured_cmap["colormap"] == "RdBu_r"
    assert params.colormap == ""  # original unchanged
    fig.clf()


@pytest.mark.asyncio
async def test_spatial_cnv_feature_list_uses_first_entry(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["cnv_score"] = np.linspace(0.0, 1.0, adata.n_obs)
    adata.obs["numbat_p_cnv"] = np.linspace(1.0, 2.0, adata.n_obs)

    monkeypatch.setattr(viz_cnv, "require_spatial_coords", lambda _a: _a.obsm["spatial"])
    fig, ax = plt.subplots()
    monkeypatch.setattr(viz_cnv, "create_figure_from_params", lambda *_a, **_k: (fig, [ax]))
    captured: dict[str, str] = {}
    monkeypatch.setattr(
        viz_cnv,
        "plot_spatial_feature",
        lambda _adata, _ax, feature, params: captured.__setitem__("feature", feature),
    )

    await viz_cnv._create_spatial_cnv(
        adata,
        VisualizationParameters(
            plot_type="cnv",
            subtype="spatial",
            feature=["cnv_score", "numbat_p_cnv"],
            colormap="",
        ),
        context=DummyCtx(),
    )
    assert captured["feature"] == "cnv_score"
    fig.clf()


@pytest.mark.asyncio
async def test_spatial_cnv_auto_detects_numbat_probability(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["numbat_p_cnv"] = np.linspace(0.0, 1.0, adata.n_obs)

    monkeypatch.setattr(viz_cnv, "require_spatial_coords", lambda _a: _a.obsm["spatial"])
    fig, ax = plt.subplots()
    monkeypatch.setattr(viz_cnv, "create_figure_from_params", lambda *_a, **_k: (fig, [ax]))
    captured: dict[str, str] = {}
    monkeypatch.setattr(
        viz_cnv,
        "plot_spatial_feature",
        lambda _adata, _ax, feature, params: captured.__setitem__("feature", feature),
    )
    ctx = DummyCtx()

    await viz_cnv._create_spatial_cnv(
        adata,
        VisualizationParameters(plot_type="cnv", subtype="spatial", colormap=""),
        context=ctx,
    )
    assert captured["feature"] == "numbat_p_cnv"
    assert any("using 'numbat_p_cnv'" in msg for msg in ctx.infos)
    fig.clf()


@pytest.mark.asyncio
async def test_cnv_heatmap_requires_cnv_data(minimal_spatial_adata):
    ctx = DummyCtx()
    with pytest.raises(DataNotFoundError, match="CNV data not found in obsm"):
        await viz_cnv._create_cnv_heatmap(
            minimal_spatial_adata.copy(),
            VisualizationParameters(plot_type="cnv", subtype="heatmap"),
            context=ctx,
        )
    assert any("CNV data not found" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_cnv_heatmap_requires_uns_cnv_metadata(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_cnv"] = np.zeros((adata.n_obs, 10), dtype=float)
    monkeypatch.setattr(viz_cnv, "require", lambda *_a, **_k: None)
    ctx = DummyCtx()

    with pytest.raises(DataNotFoundError, match="CNV metadata not found"):
        await viz_cnv._create_cnv_heatmap(
            adata,
            VisualizationParameters(plot_type="cnv", subtype="heatmap"),
            context=ctx,
        )
    assert any("CNV metadata not found" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_cnv_heatmap_numbat_aggregated_branch(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_cnv_numbat"] = np.random.default_rng(1).normal(size=(adata.n_obs, 50))
    adata.obs["clone"] = ["A"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)

    monkeypatch.setattr(viz_cnv, "require", lambda *_a, **_k: None)
    ctx = DummyCtx()

    fig = await viz_cnv._create_cnv_heatmap(
        adata,
        VisualizationParameters(plot_type="cnv", subtype="heatmap", feature="clone"),
        context=ctx,
    )

    # After fix: original adata should NOT be mutated (copy is used internally)
    assert "X_cnv" not in adata.obsm
    assert "cnv" not in adata.uns
    assert "clone" in fig.axes[0].get_ylabel()
    assert any("Converting Numbat CNV data" in msg for msg in ctx.infos)
    fig.clf()


@pytest.mark.asyncio
async def test_cnv_heatmap_numbat_ungrouped_branch(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_cnv_numbat"] = np.random.default_rng(2).normal(size=(adata.n_obs, 30))

    monkeypatch.setattr(viz_cnv, "require", lambda *_a, **_k: None)
    fig = await viz_cnv._create_cnv_heatmap(
        adata,
        VisualizationParameters(plot_type="cnv", subtype="heatmap"),
        context=DummyCtx(),
    )
    assert "All cells (ungrouped)" in fig.axes[0].get_title()
    fig.clf()


@pytest.mark.asyncio
async def test_cnv_heatmap_infercnvpy_chromosome_branch(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_cnv"] = np.random.default_rng(3).normal(size=(adata.n_obs, 12))
    adata.uns["cnv"] = {"genomic_positions": True}
    adata.var["chromosome"] = ["chr1"] * adata.n_vars

    monkeypatch.setattr(viz_cnv, "require", lambda *_a, **_k: None)
    captured: dict[str, object] = {}

    fake_infercnvpy = ModuleType("infercnvpy")

    def _chrom_heatmap(*_args, **kwargs):
        captured["groupby"] = kwargs.get("groupby")
        captured["dendrogram"] = kwargs.get("dendrogram")
        plt.figure()

    fake_infercnvpy.pl = SimpleNamespace(chromosome_heatmap=_chrom_heatmap)
    monkeypatch.setitem(sys.modules, "infercnvpy", fake_infercnvpy)

    fig = await viz_cnv._create_cnv_heatmap(
        adata,
        VisualizationParameters(plot_type="cnv", subtype="heatmap", cluster_key="group"),
        context=DummyCtx(),
    )
    assert captured["groupby"] == "group"
    assert captured["dendrogram"] is True
    assert fig is not None
    fig.clf()


@pytest.mark.asyncio
async def test_cnv_heatmap_infercnvpy_auto_selects_groupby(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_cnv"] = np.random.default_rng(4).normal(size=(adata.n_obs, 12))
    adata.uns["cnv"] = {"genomic_positions": True}
    adata.var["chromosome"] = ["chr1"] * adata.n_vars
    adata.obs["cell_type"] = ["A"] * (adata.n_obs // 2) + [
        "B"
    ] * (adata.n_obs - adata.n_obs // 2)

    monkeypatch.setattr(viz_cnv, "require", lambda *_a, **_k: None)
    captured: dict[str, object] = {}

    fake_infercnvpy = ModuleType("infercnvpy")

    def _chrom_heatmap(*_args, **kwargs):
        captured["groupby"] = kwargs.get("groupby")
        plt.figure()

    fake_infercnvpy.pl = SimpleNamespace(chromosome_heatmap=_chrom_heatmap)
    monkeypatch.setitem(sys.modules, "infercnvpy", fake_infercnvpy)
    ctx = DummyCtx()

    fig = await viz_cnv._create_cnv_heatmap(
        adata,
        VisualizationParameters(plot_type="cnv", subtype="heatmap"),
        context=ctx,
    )
    assert captured["groupby"] == "cell_type"
    assert any("Grouping CNV heatmap by 'cell_type'" in msg for msg in ctx.infos)
    fig.clf()


@pytest.mark.asyncio
async def test_cnv_heatmap_infercnvpy_requires_groupby_when_none_available(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    del adata.obs["group"]
    adata.obsm["X_cnv"] = np.random.default_rng(5).normal(size=(adata.n_obs, 12))
    adata.uns["cnv"] = {"genomic_positions": True}
    adata.var["chromosome"] = ["chr1"] * adata.n_vars

    monkeypatch.setattr(viz_cnv, "require", lambda *_a, **_k: None)

    with pytest.raises(ParameterError, match="requires a grouping column"):
        await viz_cnv._create_cnv_heatmap(
            adata,
            VisualizationParameters(plot_type="cnv", subtype="heatmap"),
            context=DummyCtx(),
        )
