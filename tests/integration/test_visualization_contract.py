"""Integration contract tests for visualize_data output behavior."""

from __future__ import annotations

import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.server import visualize_data
from tests.fixtures.helpers import extract_saved_path, load_generic_dataset


@pytest.mark.integration
@pytest.mark.asyncio
async def test_visualize_data_appends_suffix_for_custom_output_path(
    spatial_dataset_path, tmp_path, reset_data_manager
):
    dataset = await load_generic_dataset(
        spatial_dataset_path, name="viz_contract_suffix"
    )

    out_no_suffix = tmp_path / "figures" / "feature_plot"
    params = VisualizationParameters(
        plot_type="feature",
        feature="gene_0",
        basis="spatial",
        output_path=str(out_no_suffix),
        output_format="pdf",
        dpi=90,
    )
    result = await visualize_data(dataset.id, params=params)

    saved = extract_saved_path(result)
    assert saved.suffix == ".pdf"
    assert saved.exists()
    assert saved.stat().st_size > 0
    assert "Type: feature" in result
    assert "Format: PDF" in result


@pytest.mark.integration
@pytest.mark.asyncio
async def test_visualize_data_default_output_directory_contract(
    spatial_dataset_path, tmp_path, monkeypatch, reset_data_manager
):
    monkeypatch.chdir(tmp_path)

    dataset = await load_generic_dataset(
        spatial_dataset_path, name="viz_contract_default"
    )

    params = VisualizationParameters(
        plot_type="feature",
        feature="gene_1",
        basis="spatial",
        output_format="png",
        dpi=72,
    )
    result = await visualize_data(dataset.id, params=params)

    saved = extract_saved_path(result)
    assert saved.exists()
    assert saved.suffix == ".png"
    assert saved.parent == (tmp_path / "visualizations")
    assert saved.stat().st_size > 0
    assert "Format: PNG" in result
