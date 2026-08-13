"""Integration contract tests for find_markers behavior."""

from __future__ import annotations

import pytest

from chatspatial.models.data import DifferentialExpressionParameters
from chatspatial.server import data_manager, find_markers
from chatspatial.utils.exceptions import ParameterError
from tests.fixtures.helpers import load_generic_dataset


@pytest.mark.integration
@pytest.mark.asyncio
async def test_find_markers_all_groups_contract(
    spatial_dataset_path, reset_data_manager
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="de_contract")

    params = DifferentialExpressionParameters(
        group_key="group",
        method="wilcoxon",
        n_top_genes=5,
        min_cells=3,
    )
    result = await find_markers(data_id=dataset.id, params=params)

    assert result.data_id == dataset.id
    assert result.comparison == "All groups in group"
    # Every group keeps its own ranking; a flat list could only carry one.
    assert len(result.top_genes_by_group) > 1
    assert all(0 < len(g) <= 5 for g in result.top_genes_by_group.values())
    assert result.n_genes == len(
        {gene for genes in result.top_genes_by_group.values() for gene in genes}
    )

    stored = await data_manager.get_dataset(dataset.id)
    assert "rank_genes_groups" in stored["adata"].uns


@pytest.mark.integration
@pytest.mark.asyncio
async def test_find_markers_group_vs_rest_contract(
    spatial_dataset_path, reset_data_manager
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="de_contract")

    params = DifferentialExpressionParameters(
        group_key="group",
        group1="A",
        method="wilcoxon",
        n_top_genes=6,
        min_cells=3,
    )
    result = await find_markers(data_id=dataset.id, params=params)

    assert result.comparison == "A vs rest"
    assert list(result.top_genes_by_group) == ["A"]
    assert result.n_genes == len(result.top_genes_by_group["A"])
    assert 0 < result.n_genes <= 6


@pytest.mark.integration
@pytest.mark.asyncio
async def test_find_markers_specific_groups_contract(
    spatial_dataset_path, reset_data_manager
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="de_contract")

    params = DifferentialExpressionParameters(
        group_key="group",
        group1="A",
        group2="B",
        method="wilcoxon",
        n_top_genes=4,
        min_cells=3,
    )
    result = await find_markers(data_id=dataset.id, params=params)

    assert result.comparison == "A vs B"
    assert list(result.top_genes_by_group) == ["A"]
    assert result.n_genes == len(result.top_genes_by_group["A"])
    assert 0 < result.n_genes <= 4


@pytest.mark.integration
@pytest.mark.asyncio
async def test_find_markers_invalid_group_raises_parameter_error(
    spatial_dataset_path, reset_data_manager
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="de_contract")

    params = DifferentialExpressionParameters(
        group_key="group",
        group1="Z",
        method="wilcoxon",
        n_top_genes=5,
        min_cells=3,
    )
    with pytest.raises(ParameterError, match="Group 'Z' not found"):
        await find_markers(data_id=dataset.id, params=params)
