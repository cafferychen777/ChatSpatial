"""Integration tests for DefaultSpatialDataManager behavior."""

import pytest

from chatspatial.spatial_mcp_adapter import DefaultSpatialDataManager
from chatspatial.utils.exceptions import DataNotFoundError, ParameterError


@pytest.fixture
def manager() -> DefaultSpatialDataManager:
    """Fresh manager per test to avoid shared state."""
    return DefaultSpatialDataManager()


@pytest.mark.integration
@pytest.mark.asyncio
async def test_data_manager_create_update_get(manager, minimal_spatial_adata):
    data_id = await manager.create_dataset(
        minimal_spatial_adata, prefix="custom", name="demo"
    )
    ds = await manager.get_dataset(data_id)
    assert data_id.startswith("custom_")
    assert ds["name"] == "demo"
    assert ds["adata"].n_obs == minimal_spatial_adata.n_obs
    assert ds["n_cells"] == minimal_spatial_adata.n_obs
    assert ds["n_genes"] == minimal_spatial_adata.n_vars

    subset = minimal_spatial_adata[:20, :10].copy()
    await manager.update_adata(data_id, subset)
    ds2 = await manager.get_dataset(data_id)
    assert ds2["adata"].shape == (20, 10)
    assert ds2["n_cells"] == 20
    assert ds2["n_genes"] == 10

    listed = await manager.list_datasets()
    assert listed == [
        {
            "id": data_id,
            "name": "demo",
            "type": "unknown",
            "n_cells": 20,
            "n_genes": 10,
        }
    ]


@pytest.mark.integration
@pytest.mark.asyncio
async def test_data_manager_auto_generated_ids_are_unique(
    manager, minimal_spatial_adata
):
    id1 = await manager.create_dataset(minimal_spatial_adata, prefix="dup")
    id2 = await manager.create_dataset(minimal_spatial_adata, prefix="dup")

    assert id1 != id2
    assert id1.startswith("dup_")
    assert id2.startswith("dup_")

    with pytest.raises(DataNotFoundError):
        await manager.get_dataset("does_not_exist")


@pytest.mark.integration
@pytest.mark.asyncio
async def test_data_manager_can_publish_a_reserved_id(manager, minimal_spatial_adata):
    reserved = manager.reserve_dataset_id("integrated")

    assert reserved not in manager.data_store
    created = await manager.create_dataset(
        minimal_spatial_adata,
        prefix="integrated",
        data_id=reserved,
    )

    assert created == reserved
    assert (await manager.get_dataset(reserved))["adata"] is minimal_spatial_adata

    with pytest.raises(ParameterError, match="already exists"):
        await manager.create_dataset(
            minimal_spatial_adata,
            prefix="integrated",
            data_id=reserved,
        )


@pytest.mark.integration
@pytest.mark.asyncio
async def test_data_manager_publish_dataset_detaches_metadata(
    manager, minimal_spatial_adata
):
    data_id = manager.reserve_dataset_id("staged")
    staged = {
        "adata": minimal_spatial_adata,
        "name": "staged-data",
        "type": "generic",
    }

    await manager.publish_dataset(data_id, staged)
    staged["name"] = "mutated-after-publish"

    stored = await manager.get_dataset(data_id)
    assert stored["name"] == "staged-data"
    assert stored["n_cells"] == minimal_spatial_adata.n_obs


@pytest.mark.integration
@pytest.mark.asyncio
async def test_data_manager_reservation_skips_existing_generated_id(
    manager, minimal_spatial_adata
):
    manager.data_store["data_1"] = {"adata": minimal_spatial_adata}

    assert manager.reserve_dataset_id("data") == "data_2"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_data_manager_multi_dataset_update_is_atomic(
    manager,
    minimal_spatial_adata,
):
    first_id = await manager.create_dataset(minimal_spatial_adata, prefix="slice")
    original = (await manager.get_dataset(first_id))["adata"]
    candidate = original[:20, :10].copy()

    with pytest.raises(DataNotFoundError, match="missing"):
        await manager.update_adatas(
            {
                first_id: candidate,
                "missing": candidate,
            }
        )

    stored = await manager.get_dataset(first_id)
    assert stored["adata"] is original
    assert stored["n_cells"] == original.n_obs
    assert stored["n_genes"] == original.n_vars
