"""Unit tests for shared visualization core utilities."""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import seaborn as sns

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import core as viz_core
from chatspatial.utils.exceptions import DataNotFoundError, ParameterError


class _Ctx:
    def __init__(self):
        self.warnings: list[str] = []

    async def warning(self, message: str):
        self.warnings.append(message)


def test_resolve_figure_size_prioritizes_user_size_and_defaults():
    user_params = VisualizationParameters(plot_type="feature", figure_size=(7, 5))
    auto_params = VisualizationParameters(plot_type="feature", figure_size=None)

    assert viz_core.resolve_figure_size(user_params, "spatial") == (7, 5)
    assert viz_core.resolve_figure_size(auto_params, n_panels=4) == (15, 8)
    assert viz_core.resolve_figure_size(auto_params, "unknown") == (10, 8)


def test_temporary_obs_column_never_overwrites_existing_data(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    original = pd.Series(
        np.arange(adata.n_obs), index=adata.obs_names, name="_temporary"
    )
    adata.obs["_temporary"] = original

    with viz_core.temporary_obs_column(
        adata, "_temporary", np.repeat(-1, adata.n_obs)
    ) as key:
        assert key == "_temporary_1"
        np.testing.assert_array_equal(adata.obs[key], -1)
        pd.testing.assert_series_equal(adata.obs["_temporary"], original)

    assert "_temporary_1" not in adata.obs
    pd.testing.assert_series_equal(adata.obs["_temporary"], original)


def test_create_figure_from_params_normalizes_axes_shape():
    params = VisualizationParameters(plot_type="feature", dpi=123)

    fig_single, axes_single = viz_core.create_figure_from_params(
        params, n_rows=1, n_cols=1, squeeze=True
    )
    fig_row, axes_row = viz_core.create_figure_from_params(
        params, n_rows=1, n_cols=3, squeeze=True
    )

    assert isinstance(axes_single, np.ndarray)
    assert len(axes_single) == 1
    assert isinstance(axes_row, np.ndarray)
    assert len(axes_row) == 3
    assert int(fig_single.dpi) == 123

    plt.close(fig_single)
    plt.close(fig_row)


def test_setup_multi_panel_figure_sets_title_and_hides_unused_axes():
    params = VisualizationParameters(
        plot_type="feature",
        panel_layout=(1, 3),
        figure_size=(9, 3),
        title=None,
    )

    fig, axes = viz_core.setup_multi_panel_figure(
        n_panels=2, params=params, default_title="Default title", use_tight_layout=False
    )

    assert len(axes) == 3
    assert fig._suptitle is not None
    assert fig._suptitle.get_text() == "Default title"
    assert not axes[2].axison

    plt.close(fig)


def test_setup_multi_panel_figure_supports_tight_layout_mode():
    params = VisualizationParameters(plot_type="feature", panel_layout=(1, 2))
    fig, axes = viz_core.setup_multi_panel_figure(
        n_panels=2, params=params, default_title="", use_tight_layout=True
    )
    assert len(axes) == 2
    plt.close(fig)


def test_safe_tight_layout_supports_pyplot_and_figure_targets():
    fig, _ = plt.subplots()
    viz_core.safe_tight_layout()
    viz_core.safe_tight_layout(fig=fig)
    plt.close(fig)


@pytest.mark.asyncio
async def test_get_validated_features_supports_genes_obs_obsm_and_truncation(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cluster"] = pd.Categorical(["A"] * adata.n_obs)
    ctx = _Ctx()
    params = VisualizationParameters(
        plot_type="feature",
        feature=["gene_0", "gene_1", "gene_2", "cluster", "missing_feature"],
    )

    validated = await viz_core.get_validated_features(
        adata, params, context=ctx, max_features=3
    )

    # gene_0, gene_1, gene_2, cluster pass; missing_feature rejected; truncated to 3
    assert validated == ["gene_0", "gene_1", "gene_2"]
    assert any("not found in genes or obs" in w for w in ctx.warnings)
    assert any("Too many features" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_get_validated_features_genes_only_warns_for_non_genes(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cluster"] = ["A"] * adata.n_obs
    ctx = _Ctx()
    params = VisualizationParameters(plot_type="feature", feature=["gene_0", "cluster"])

    validated = await viz_core.get_validated_features(
        adata, params, context=ctx, genes_only=True
    )

    assert validated == ["gene_0"]
    assert any("not found in var_names" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_get_validated_features_returns_empty_when_feature_is_none(
    minimal_spatial_adata,
):
    params = VisualizationParameters(plot_type="feature", feature=None)
    validated = await viz_core.get_validated_features(minimal_spatial_adata, params)
    assert validated == []


def test_validate_and_prepare_feature_handles_gene_obs_and_missing(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cluster"] = pd.Categorical(["A"] * adata.n_obs)

    gene_values, gene_name, gene_is_cat = viz_core.validate_and_prepare_feature(
        adata, "gene_0"
    )
    obs_values, obs_name, obs_is_cat = viz_core.validate_and_prepare_feature(
        adata, "cluster"
    )

    assert gene_name == "gene_0"
    assert gene_is_cat is False
    assert isinstance(gene_values, np.ndarray)
    assert obs_name == "cluster"
    assert obs_is_cat is True
    assert len(obs_values) == adata.n_obs

    with pytest.raises(DataNotFoundError, match="Feature 'missing' not found"):
        viz_core.validate_and_prepare_feature(adata, "missing")


def test_categorical_colormap_selection_and_fallback():
    assert viz_core.get_categorical_cmap(8, user_cmap="Set2") == "Set2"
    assert viz_core.get_categorical_cmap(8, user_cmap="coolwarm") == "tab10"
    assert viz_core.get_categorical_cmap(15) == "tab20"
    assert viz_core.get_categorical_cmap(100) == "husl"


def test_category_colors_and_colormap_variants():
    auto_colors = viz_core.get_category_colors(3)
    fallback_colors = viz_core.get_category_colors(4, cmap_name="viridis")
    explicit_colors = viz_core.get_category_colors(4, cmap_name="Set2")

    indexed_colors = viz_core.get_colormap("tab10", n_colors=3)
    seaborn_palette = viz_core.get_colormap("Set2")
    mpl_cmap = viz_core.get_colormap("viridis")

    assert len(auto_colors) == 3
    assert fallback_colors == viz_core.get_category_colors(4)
    assert explicit_colors == sns.color_palette("Set2", n_colors=4)
    assert len(indexed_colors) == 3
    assert len(seaborn_palette) > 0
    assert hasattr(mpl_cmap, "__call__")
    assert viz_core.get_diverging_colormap() == "RdBu_r"


def test_tab20_orders_distinct_hues_before_lighter_companions():
    colors = viz_core.get_category_colors(15, cmap_name="tab20")
    tab20 = sns.color_palette("tab20", n_colors=20)

    assert colors[:10] == tab20[::2]
    assert colors[10:] == tab20[1::2][:5]


def test_auto_spot_size_prioritizes_user_metadata_and_adaptive_fallback(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()

    from_user = viz_core.auto_spot_size(adata, user_spot_size=77.0)
    from_metadata = viz_core.auto_spot_size(adata, user_spot_size=None, basis="spatial")

    tiny = viz_core.auto_spot_size(
        SimpleNamespace(n_obs=100_000, uns={}),
        user_spot_size=None,
        basis="umap",
    )
    huge = viz_core.auto_spot_size(
        SimpleNamespace(n_obs=10, uns={}),
        user_spot_size=None,
        basis="umap",
    )

    assert from_user == 77.0
    assert from_metadata > 100.0
    assert tiny == 5.0
    assert huge == 200.0


def test_plot_spatial_feature_supports_continuous_and_categorical_values(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    params = VisualizationParameters(
        plot_type="feature",
        feature="gene_0",
        show_axes=False,
    )
    fig, ax = plt.subplots()
    mappable = viz_core.plot_spatial_feature(
        adata,
        ax=ax,
        feature="gene_0",
        params=params,
        title="Gene panel",
    )
    assert mappable is not None
    assert ax.get_title() == "Gene panel"
    assert not ax.axison
    plt.close(fig)

    adata.obs["cluster"] = pd.Categorical(
        ["A"] * (adata.n_obs // 2) + ["B"] * (adata.n_obs - adata.n_obs // 2)
    )
    fig2, ax2 = plt.subplots()
    categorical_values = adata.obs["cluster"].values
    mappable2 = viz_core.plot_spatial_feature(
        adata,
        ax=ax2,
        values=categorical_values,
        params=VisualizationParameters(plot_type="feature", show_legend=True),
    )
    assert mappable2 is None
    assert ax2.get_legend() is not None
    plt.close(fig2)


def test_plot_spatial_feature_validates_required_inputs_and_feature_presence(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    fig, ax = plt.subplots()
    with pytest.raises(DataNotFoundError, match="Feature 'missing' not found"):
        viz_core.plot_spatial_feature(adata, ax=ax, feature="missing")

    with pytest.raises(
        ParameterError, match="Either feature or values must be provided"
    ):
        viz_core.plot_spatial_feature(adata, ax=ax, feature=None, values=None)
    plt.close(fig)


def test_plot_spatial_feature_handles_object_string_obs_column(
    minimal_spatial_adata,
):
    """Regression: object/string obs columns must be treated as categorical."""
    adata = minimal_spatial_adata.copy()
    labels = ["A", "B"] * (adata.n_obs // 2) + ["A"] * (adata.n_obs % 2)
    adata.obs["label"] = pd.Series(labels, index=adata.obs_names, dtype="string")
    assert isinstance(adata.obs["label"].dtype, pd.StringDtype)

    fig, ax = plt.subplots()
    mappable = viz_core.plot_spatial_feature(
        adata, ax=ax, feature="label", params=VisualizationParameters()
    )
    assert mappable is None  # categorical → no mappable
    plt.close(fig)


def test_plot_spatial_feature_handles_categorical_with_nan(
    minimal_spatial_adata,
):
    """Regression: categorical obs with NaN must not KeyError."""
    adata = minimal_spatial_adata.copy()
    labels = ["A", None, "B"] * (adata.n_obs // 3)
    labels += ["A"] * (adata.n_obs - len(labels))
    adata.obs["label"] = pd.Categorical(labels)

    fig, ax = plt.subplots()
    mappable = viz_core.plot_spatial_feature(
        adata, ax=ax, feature="label", params=VisualizationParameters()
    )
    assert mappable is None
    plt.close(fig)


def test_get_categorical_columns_and_infer_basis(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["cluster"] = pd.Categorical(["A"] * adata.n_obs)
    adata.obs["sample"] = pd.Series(
        ["s1"] * adata.n_obs,
        index=adata.obs_names,
        dtype="string",
    )
    adata.obsm["X_tsne"] = np.zeros((adata.n_obs, 2))

    cols = viz_core.get_categorical_columns(adata)
    assert "cluster" in cols
    assert "sample" in cols
    assert viz_core.get_categorical_columns(adata, limit=1) == cols[:1]

    assert viz_core.infer_basis(adata, preferred="tsne") == "tsne"
    assert viz_core.infer_basis(adata, priority=["pca", "umap", "spatial"]) == "spatial"

    adata_no_spatial = adata.copy()
    del adata_no_spatial.obsm["spatial"]
    del adata_no_spatial.obsm["X_tsne"]
    adata_no_spatial.obsm["X_custom"] = np.ones((adata.n_obs, 2))
    assert viz_core.infer_basis(adata_no_spatial) == "custom"

    adata_empty = adata.copy()
    adata_empty.obsm.clear()
    assert viz_core.infer_basis(adata_empty) is None
