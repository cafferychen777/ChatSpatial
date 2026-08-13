"""Stable end-to-end workflows for default test runs."""

from __future__ import annotations

import pandas as pd
import pytest

from chatspatial.models.data import (
    AnnotationParameters,
    DifferentialExpressionParameters,
    PreprocessingParameters,
    VisualizationParameters,
)
from chatspatial.server import (
    annotate_cell_types,
    find_markers,
    preprocess_data,
    visualize_data,
)
from chatspatial.tools import annotation as annotation_module
from tests.fixtures.helpers import extract_saved_path, load_generic_dataset


@pytest.mark.e2e
@pytest.mark.asyncio
async def test_e2e_load_preprocess_find_markers(
    spatial_dataset_path, reset_data_manager, e2e_trace
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="e2e_markers")
    e2e_trace.record(
        step="load_data",
        data_id=dataset.id,
        params={
            "data_path": str(spatial_dataset_path),
            "data_type": "generic",
            "name": "e2e_markers",
        },
    )

    preprocess_params = PreprocessingParameters(
        filter_genes_min_cells=1,
        filter_cells_min_genes=1,
        n_hvgs=20,
        subsample_genes=20,
        normalization="log",
    )
    e2e_trace.record(
        step="preprocess_data", data_id=dataset.id, params=preprocess_params
    )
    prep = await preprocess_data(
        dataset.id,
        params=preprocess_params,
    )

    marker_params = DifferentialExpressionParameters(
        group_key="group", method="wilcoxon", n_top_genes=5, min_cells=3
    )
    e2e_trace.record(step="find_markers", data_id=dataset.id, params=marker_params)
    result = await find_markers(dataset.id, params=marker_params)

    ctx = (
        f"data_id={dataset.id}, "
        "preprocess_params={filter_genes_min_cells=1,filter_cells_min_genes=1,n_hvgs=20,subsample_genes=20,normalization=log}, "
        "marker_params={group_key=group,method=wilcoxon,n_top_genes=5,min_cells=3}"
    )
    assert prep.n_cells > 0, ctx
    assert prep.n_genes > 0, ctx
    assert result.n_genes >= 0, ctx
    assert result.top_genes_by_group, ctx
    assert all(result.top_genes_by_group.values()), ctx


@pytest.mark.e2e
@pytest.mark.asyncio
async def test_e2e_load_annotate_with_mock(
    spatial_dataset_path, monkeypatch, reset_data_manager, e2e_trace
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="e2e_annotation")
    e2e_trace.record(
        step="load_data",
        data_id=dataset.id,
        params={
            "data_path": str(spatial_dataset_path),
            "data_type": "generic",
            "name": "e2e_annotation",
        },
    )

    async def _fake_sctype(adata, params, ctx, output_key, confidence_key):
        labels = [
            "mock_type_a" if i % 2 == 0 else "mock_type_b" for i in range(adata.n_obs)
        ]
        adata.obs[output_key] = pd.Categorical(labels)
        adata.obs[confidence_key] = [0.9] * adata.n_obs
        counts = adata.obs[output_key].value_counts().to_dict()
        return annotation_module.AnnotationMethodOutput(
            cell_types=sorted(set(labels)),
            counts=counts,
            confidence={"mock_type_a": 0.9, "mock_type_b": 0.9},
        )

    monkeypatch.setattr(annotation_module, "_annotate_with_sctype", _fake_sctype)

    annotation_params = AnnotationParameters(method="sctype", sctype_tissue="Brain")
    e2e_trace.record(
        step="annotate_cell_types", data_id=dataset.id, params=annotation_params
    )
    result = await annotate_cell_types(
        dataset.id,
        params=annotation_params,
    )

    ctx = (
        f"data_id={dataset.id}, "
        "annotation_params={method=sctype,sctype_tissue=Brain,mock=true}"
    )
    assert result.method == "sctype", ctx
    assert set(result.cell_types) == {"mock_type_a", "mock_type_b"}, ctx
    assert result.output_key == "cell_type_sctype", ctx


@pytest.mark.e2e
@pytest.mark.asyncio
async def test_e2e_load_visualize_feature_png(
    spatial_dataset_path, tmp_path, reset_data_manager, e2e_trace
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="e2e_viz")
    e2e_trace.record(
        step="load_data",
        data_id=dataset.id,
        params={
            "data_path": str(spatial_dataset_path),
            "data_type": "generic",
            "name": "e2e_viz",
        },
    )

    output_file = tmp_path / "feature_plot.png"
    viz_params = VisualizationParameters(
        plot_type="feature",
        feature="gene_0",
        basis="spatial",
        output_path=str(output_file),
        output_format="png",
        dpi=72,
    )
    e2e_trace.record(step="visualize_data", data_id=dataset.id, params=viz_params)
    viz_result = await visualize_data(
        dataset.id,
        params=viz_params,
    )

    saved_path = extract_saved_path(viz_result)
    ctx = (
        f"data_id={dataset.id}, "
        f"visualize_params={VisualizationParameters(plot_type='feature',feature='gene_0',basis='spatial',output_path='{output_file}',output_format='png',dpi=72)}"
    )
    assert saved_path.exists(), ctx
    assert saved_path.suffix == ".png", ctx
    assert saved_path.stat().st_size > 0, ctx
