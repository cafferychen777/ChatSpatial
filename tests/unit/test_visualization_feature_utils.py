"""Unit tests for feature visualization contracts."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import feature as viz_feature
from chatspatial.utils.exceptions import DataError, DataNotFoundError, ParameterError


class DummyCtx:
    def __init__(self):
        self.infos: list[str] = []

    async def info(self, msg: str):
        self.infos.append(msg)

    async def warning(self, msg: str):
        self.infos.append(msg)


def test_parse_lr_pairs_from_features_contract():
    # Only ^ is recognized as LR separator
    regular, lr_pairs = viz_feature._parse_lr_pairs_from_features(
        ["CCL5^CCR5", "CXCL12^CXCR4"]
    )
    assert regular == []
    assert ("CCL5", "CCR5") in lr_pairs
    assert ("CXCL12", "CXCR4") in lr_pairs


def test_parse_lr_pairs_underscore_treated_as_gene():
    """Underscored names like HLA_DRA are regular genes, not LR pairs."""
    regular, lr_pairs = viz_feature._parse_lr_pairs_from_features(
        ["HLA_DRA", "CD3D_CD3E", "_temp"]
    )
    assert regular == ["HLA_DRA", "CD3D_CD3E", "_temp"]
    assert lr_pairs == []


def test_parse_lr_pairs_mixed_raises():
    """Mixing LR pairs with regular features raises ParameterError."""
    with pytest.raises(ParameterError, match="Cannot mix"):
        viz_feature._parse_lr_pairs_from_features(["geneA", "CCL5^CCR5"])


@pytest.mark.asyncio
async def test_create_feature_visualization_invalid_basis(minimal_spatial_adata):
    with pytest.raises(ParameterError, match="Invalid basis"):
        await viz_feature.create_feature_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="feature", basis="bad", feature="gene_0"),
        )


@pytest.mark.asyncio
async def test_create_feature_visualization_requires_features_when_no_cluster(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(viz_feature, "get_cluster_key", lambda *_args, **_kwargs: None)
    with pytest.raises(ParameterError, match="No features specified"):
        await viz_feature.create_feature_visualization(
            adata,
            VisualizationParameters(plot_type="feature", basis="spatial", feature=None),
        )


@pytest.mark.asyncio
async def test_create_feature_visualization_umap_fallback_computation(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx()

    def _ensure_umap(a):
        a.obsm["X_umap"] = np.ones((a.n_obs, 2))
        return True

    monkeypatch.setattr(viz_feature, "ensure_umap", _ensure_umap)

    async def _validated(*_args, **_kwargs):
        return ["gene_0"]

    monkeypatch.setattr(viz_feature, "get_validated_features", _validated)

    async def _single(*_args, **_kwargs):
        return plt.figure()

    monkeypatch.setattr(viz_feature, "_create_single_feature_plot", _single)

    fig = await viz_feature.create_feature_visualization(
        adata,
        VisualizationParameters(plot_type="feature", basis="umap", feature="gene_0"),
        context=ctx,
    )
    assert fig is not None
    assert any("Computed UMAP embedding" in m for m in ctx.infos)


@pytest.mark.asyncio
async def test_create_feature_visualization_routes_lr_pairs(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    sentinel = object()

    async def _lr(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_feature, "_create_lr_pairs_visualization", _lr)

    out = await viz_feature.create_feature_visualization(
        adata,
        VisualizationParameters(
            plot_type="feature",
            basis="spatial",
            feature=["CCL5^CCR5"],
        ),
        context=DummyCtx(),
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_create_single_feature_plot_feature_not_found(minimal_spatial_adata):
    with pytest.raises(DataNotFoundError, match="Feature 'missing' not found"):
        await viz_feature._create_single_feature_plot(
            minimal_spatial_adata,
            VisualizationParameters(
                plot_type="feature", basis="spatial", feature="missing"
            ),
            "missing",
            "spatial",
            minimal_spatial_adata.obsm["spatial"],
        )


@pytest.mark.asyncio
async def test_create_single_feature_plot_gene_branch(minimal_spatial_adata):
    fig = await viz_feature._create_single_feature_plot(
        minimal_spatial_adata,
        VisualizationParameters(plot_type="feature", basis="spatial", feature="gene_0"),
        "gene_0",
        "spatial",
        minimal_spatial_adata.obsm["spatial"],
    )
    assert fig is not None
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_single_feature_plot_categorical_obs_branch(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["cluster"] = pd.Categorical(
        ["A"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )

    fig = await viz_feature._create_single_feature_plot(
        adata,
        VisualizationParameters(
            plot_type="feature", basis="spatial", feature="cluster"
        ),
        "cluster",
        "spatial",
        adata.obsm["spatial"],
    )
    assert fig is not None
    expected_colors = viz_feature.get_category_colors(2)
    actual_colors = [
        collection.get_facecolors()[0] for collection in fig.axes[0].collections
    ]
    np.testing.assert_allclose(
        [color[:3] for color in actual_colors],
        [color[:3] for color in expected_colors],
    )
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_single_feature_plot_does_not_mutate_obs_dtype(
    minimal_spatial_adata,
):
    """Regression: viz must not convert obs column dtype as side effect."""
    adata = minimal_spatial_adata.copy()
    labels = ["A", "B"] * (adata.n_obs // 2) + ["A"] * (adata.n_obs % 2)
    adata.obs["label"] = pd.Series(labels, index=adata.obs_names, dtype="string")
    original_dtype = adata.obs["label"].dtype
    assert isinstance(original_dtype, pd.StringDtype)

    fig = await viz_feature._create_single_feature_plot(
        adata,
        VisualizationParameters(plot_type="feature", basis="spatial", feature="label"),
        "label",
        "spatial",
        adata.obsm["spatial"],
    )
    assert fig is not None
    # Must NOT have been mutated to Categorical
    assert adata.obs["label"].dtype == original_dtype
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_multi_feature_plot_cleans_temp_column(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    fig = await viz_feature._create_multi_feature_plot(
        adata,
        VisualizationParameters(
            plot_type="feature", basis="spatial", feature=["gene_0", "gene_1"]
        ),
        context=None,
        features=["gene_0", "gene_1"],
        basis="spatial",
        coords=adata.obsm["spatial"],
    )
    assert fig is not None
    assert "_feature_viz_temp_99" not in adata.obs.columns
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_feature_visualization_requires_spatial_coordinates(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]

    with pytest.raises((DataNotFoundError, DataError), match="[Ss]patial coordinates"):
        await viz_feature.create_feature_visualization(
            adata,
            VisualizationParameters(
                plot_type="feature", basis="spatial", feature="gene_0"
            ),
        )


@pytest.mark.asyncio
async def test_create_feature_visualization_works_with_alternative_spatial_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """Visualization should succeed when coords are in X_spatial instead of spatial."""
    adata = minimal_spatial_adata.copy()
    coords = adata.obsm.pop("spatial")
    adata.obsm["X_spatial"] = coords

    monkeypatch.setattr(viz_feature, "auto_spot_size", lambda _a, s, **kw: s or 1.0)

    result = await viz_feature.create_feature_visualization(
        adata,
        VisualizationParameters(plot_type="feature", basis="spatial", feature="gene_0"),
    )
    assert result is not None


@pytest.mark.asyncio
async def test_create_feature_visualization_umap_missing_when_not_computable(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm.pop("X_umap", None)
    monkeypatch.setattr(viz_feature, "ensure_umap", lambda _adata: False)

    with pytest.raises(DataNotFoundError, match="UMAP embedding not found"):
        await viz_feature.create_feature_visualization(
            adata,
            VisualizationParameters(
                plot_type="feature", basis="umap", feature="gene_0"
            ),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_create_feature_visualization_requires_pca_embedding(
    minimal_spatial_adata,
):
    with pytest.raises(DataNotFoundError, match="PCA embedding not found"):
        await viz_feature.create_feature_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="feature", basis="pca", feature="gene_0"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_create_feature_visualization_routes_multi_feature_path(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    sentinel = object()

    async def _validated(*_args, **_kwargs):
        return ["gene_0", "gene_1"]

    async def _multi(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_feature, "get_validated_features", _validated)
    monkeypatch.setattr(viz_feature, "_create_multi_feature_plot", _multi)

    out = await viz_feature.create_feature_visualization(
        minimal_spatial_adata.copy(),
        VisualizationParameters(
            plot_type="feature",
            basis="spatial",
            feature=["gene_0", "gene_1"],
        ),
        context=DummyCtx(),
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_create_feature_visualization_pca_uses_first_two_components_and_default_cluster(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_pca"] = np.column_stack(
        [np.arange(adata.n_obs), np.arange(adata.n_obs), np.arange(adata.n_obs)]
    )
    sentinel = plt.figure()

    monkeypatch.setattr(
        viz_feature, "get_cluster_key", lambda *_args, **_kwargs: "group"
    )

    async def _validated(*_args, **_kwargs):
        return ["group"]

    async def _single(_adata, _params, _feature, basis, coords):
        assert basis == "pca"
        assert coords.shape[1] == 2
        return sentinel

    monkeypatch.setattr(viz_feature, "get_validated_features", _validated)
    monkeypatch.setattr(viz_feature, "_create_single_feature_plot", _single)

    out = await viz_feature.create_feature_visualization(
        adata,
        VisualizationParameters(plot_type="feature", basis="pca", feature=None),
        context=DummyCtx(),
    )

    assert out is sentinel
    plt.close(sentinel)


@pytest.mark.asyncio
async def test_create_single_feature_plot_numeric_obs_branch_hides_axes(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["score"] = np.linspace(0, 1, adata.n_obs)
    adata.obsm["X_pca"] = np.column_stack(
        [np.arange(adata.n_obs), np.arange(adata.n_obs)]
    )

    fig = await viz_feature._create_single_feature_plot(
        adata,
        VisualizationParameters(
            plot_type="feature",
            basis="pca",
            feature="score",
            show_colorbar=False,
            show_axes=False,
        ),
        "score",
        "pca",
        adata.obsm["X_pca"],
    )
    assert fig is not None
    assert fig.axes[0].get_title() == "score (pca)"
    assert not fig.axes[0].axison
    plt.close(fig)


@pytest.mark.asyncio
@pytest.mark.parametrize("color_scale", ["log", "sqrt"])
async def test_create_single_feature_plot_applies_gene_scaling_and_umap_axis_labels(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, color_scale: str
):
    adata = minimal_spatial_adata.copy()
    coords = np.column_stack([np.arange(adata.n_obs), np.arange(adata.n_obs)])
    called = {"log": 0, "sqrt": 0}
    orig_log1p = viz_feature.np.log1p
    orig_sqrt = viz_feature.np.sqrt

    def _log1p(arr):
        called["log"] += 1
        return orig_log1p(arr)

    def _sqrt(arr):
        called["sqrt"] += 1
        return orig_sqrt(arr)

    monkeypatch.setattr(viz_feature.np, "log1p", _log1p)
    monkeypatch.setattr(viz_feature.np, "sqrt", _sqrt)

    fig = await viz_feature._create_single_feature_plot(
        adata,
        VisualizationParameters(
            plot_type="feature",
            basis="umap",
            feature="gene_0",
            color_scale=color_scale,
            show_colorbar=False,
        ),
        "gene_0",
        "umap",
        coords,
    )
    assert fig.axes[0].get_xlabel() == "UMAP1"
    assert fig.axes[0].get_ylabel() == "UMAP2"
    if color_scale == "log":
        assert called["log"] >= 1
    else:
        assert called["sqrt"] >= 1
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_single_feature_plot_numeric_obs_branch_adds_colorbar(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["score"] = np.linspace(0, 1, adata.n_obs)
    adata.obsm["X_pca"] = np.column_stack(
        [np.arange(adata.n_obs), np.arange(adata.n_obs)]
    )
    calls = {"count": 0}

    monkeypatch.setattr(
        viz_feature,
        "add_colorbar",
        lambda *_args, **_kwargs: calls.__setitem__("count", calls["count"] + 1),
    )

    fig = await viz_feature._create_single_feature_plot(
        adata,
        VisualizationParameters(
            plot_type="feature",
            basis="pca",
            feature="score",
            show_colorbar=True,
        ),
        "score",
        "pca",
        adata.obsm["X_pca"],
    )
    assert calls["count"] == 1
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_multi_feature_plot_umap_handles_gene_and_numeric_obs(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_umap"] = np.column_stack(
        [np.arange(adata.n_obs), np.arange(adata.n_obs)]
    )
    adata.obs["score"] = np.linspace(0, 10, adata.n_obs)

    fig = await viz_feature._create_multi_feature_plot(
        adata,
        VisualizationParameters(
            plot_type="feature",
            basis="umap",
            feature=["gene_0", "score"],
            show_colorbar=False,
            show_axes=False,
            color_scale="log",
        ),
        context=None,
        features=["gene_0", "score"],
        basis="umap",
        coords=adata.obsm["X_umap"],
    )

    assert fig is not None
    assert "_feature_viz_temp_99" not in adata.obs.columns
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_multi_feature_plot_breaks_when_axes_shorter_than_features(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    fig, ax = plt.subplots()
    monkeypatch.setattr(
        viz_feature,
        "setup_multi_panel_figure",
        lambda **_kwargs: (fig, np.array([ax])),
    )

    out = await viz_feature._create_multi_feature_plot(
        adata,
        VisualizationParameters(
            plot_type="feature", basis="spatial", show_colorbar=False
        ),
        context=None,
        features=["gene_0", "gene_1"],
        basis="spatial",
        coords=adata.obsm["spatial"],
    )
    assert out is fig
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_multi_feature_plot_umap_sqrt_with_colorbars(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_umap"] = np.column_stack(
        [np.arange(adata.n_obs), np.arange(adata.n_obs)]
    )
    adata.obs["score"] = np.linspace(0, 10, adata.n_obs)
    calls = {"count": 0}
    monkeypatch.setattr(
        viz_feature.plt,
        "colorbar",
        lambda *_args, **_kwargs: calls.__setitem__("count", calls["count"] + 1),
    )

    fig = await viz_feature._create_multi_feature_plot(
        adata,
        VisualizationParameters(
            plot_type="feature",
            basis="umap",
            feature=["gene_0", "score"],
            show_colorbar=True,
            color_scale="sqrt",
        ),
        context=None,
        features=["gene_0", "score"],
        basis="umap",
        coords=adata.obsm["X_umap"],
    )
    assert calls["count"] >= 2
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_multi_feature_plot_spatial_categorical_inverts_y(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cluster"] = pd.Categorical(
        ["A"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )

    fig = await viz_feature._create_multi_feature_plot(
        adata,
        VisualizationParameters(
            plot_type="feature",
            basis="spatial",
            feature=["cluster"],
            show_colorbar=False,
        ),
        context=None,
        features=["cluster"],
        basis="spatial",
        coords=adata.obsm["spatial"],
    )
    assert fig.axes[0].yaxis_inverted()
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_lr_pairs_visualization_requires_available_pairs(
    minimal_spatial_adata,
):
    with pytest.raises(DataNotFoundError, match="None of the specified LR pairs found"):
        await viz_feature._create_lr_pairs_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="feature", basis="spatial"),
            context=DummyCtx(),
            lr_pairs=[("LIG", "REC")],
            basis="spatial",
            coords=minimal_spatial_adata.obsm["spatial"],
        )


@pytest.mark.asyncio
async def test_create_lr_pairs_visualization_limits_pairs_and_reports_titles(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_umap"] = np.column_stack(
        [np.arange(adata.n_obs), np.arange(adata.n_obs)]
    )
    ctx = DummyCtx()
    lr_pairs = [
        ("gene_0", "gene_1"),
        ("gene_2", "gene_3"),
        ("gene_4", "gene_5"),
        ("gene_6", "gene_7"),
        ("gene_8", "gene_9"),
    ]

    fig = await viz_feature._create_lr_pairs_visualization(
        adata,
        VisualizationParameters(
            plot_type="feature",
            basis="umap",
            show_colorbar=False,
            show_correlation_stats=False,
            correlation_method="kendall",
            color_scale="sqrt",
        ),
        context=ctx,
        lr_pairs=lr_pairs,
        basis="umap",
        coords=adata.obsm["X_umap"],
    )

    titles = [ax.get_title() for ax in fig.axes]
    assert any("Too many LR pairs" in msg for msg in ctx.infos)
    assert any("Visualizing 4 LR pairs" in msg for msg in ctx.infos)
    assert any(" vs " in title for title in titles)
    assert "_lr_viz_temp_99" not in adata.obs.columns
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_lr_pairs_visualization_spatial_log_pearson_with_colorbars(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    calls = {"count": 0}
    monkeypatch.setattr(
        viz_feature.plt,
        "colorbar",
        lambda *_args, **_kwargs: calls.__setitem__("count", calls["count"] + 1),
    )
    monkeypatch.setattr(viz_feature, "pearsonr", lambda *_args, **_kwargs: (0.42, 0.01))

    fig = await viz_feature._create_lr_pairs_visualization(
        adata,
        VisualizationParameters(
            plot_type="feature",
            basis="spatial",
            show_colorbar=True,
            show_correlation_stats=True,
            correlation_method="pearson",
            color_scale="log",
        ),
        context=DummyCtx(),
        lr_pairs=[("gene_0", "gene_1")],
        basis="spatial",
        coords=adata.obsm["spatial"],
    )

    assert calls["count"] >= 2
    assert any("Correlation:" in ax.get_title() for ax in fig.axes)
    assert "_lr_viz_temp_99" not in adata.obs.columns
    plt.close(fig)


@pytest.mark.asyncio
async def test_create_lr_pairs_visualization_umap_spearman_with_colorbars(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_umap"] = np.column_stack(
        [np.arange(adata.n_obs), np.arange(adata.n_obs)]
    )
    calls = {"count": 0}
    monkeypatch.setattr(
        viz_feature.plt,
        "colorbar",
        lambda *_args, **_kwargs: calls.__setitem__("count", calls["count"] + 1),
    )
    monkeypatch.setattr(
        viz_feature, "spearmanr", lambda *_args, **_kwargs: (0.33, 0.02)
    )

    fig = await viz_feature._create_lr_pairs_visualization(
        adata,
        VisualizationParameters(
            plot_type="feature",
            basis="umap",
            show_colorbar=True,
            show_correlation_stats=True,
            correlation_method="spearman",
        ),
        context=DummyCtx(),
        lr_pairs=[("gene_2", "gene_3")],
        basis="umap",
        coords=adata.obsm["X_umap"],
    )

    assert calls["count"] >= 2
    assert any("Correlation:" in ax.get_title() for ax in fig.axes)
    assert "_lr_viz_temp_99" not in adata.obs.columns
    plt.close(fig)
