"""Unit contracts for spatial MCP adapter primitives."""

from __future__ import annotations

from unittest.mock import Mock

import pytest

from chatspatial import spatial_mcp_adapter as adapter
from chatspatial.utils.exceptions import DataNotFoundError


class _FakeMCPContext:
    def __init__(self):
        self.infos: list[str] = []
        self.warnings: list[str] = []
        self.errors: list[str] = []

    async def info(self, msg: str) -> None:
        self.infos.append(msg)

    async def warning(self, msg: str) -> None:
        self.warnings.append(msg)

    async def error(self, msg: str) -> None:
        self.errors.append(msg)


@pytest.mark.asyncio
async def test_data_manager_list_defaults_and_dataset_exists():
    manager = adapter.DefaultSpatialDataManager()
    manager.data_store["d1"] = {"adata": object()}

    listed = await manager.list_datasets()
    assert listed == [
        {
            "id": "d1",
            "name": "Dataset d1",
            "type": "unknown",
            "n_cells": 0,
            "n_genes": 0,
        }
    ]
    assert manager.dataset_exists("d1")
    assert not manager.dataset_exists("missing")


@pytest.mark.asyncio
async def test_data_manager_save_and_update_missing_dataset_raise(
    minimal_spatial_adata,
):
    manager = adapter.DefaultSpatialDataManager()

    with pytest.raises(DataNotFoundError, match="Dataset missing not found"):
        await manager.save_result("missing", "demo", {"ok": True})

    with pytest.raises(DataNotFoundError, match="Dataset missing not found"):
        await manager.update_adata("missing", minimal_spatial_adata.copy())


@pytest.mark.asyncio
async def test_data_manager_create_dataset_filters_reserved_metadata(
    minimal_spatial_adata,
):
    manager = adapter.DefaultSpatialDataManager()
    data_id = await manager.create_dataset(
        minimal_spatial_adata,
        prefix="derived",
        name="derived-data",
        metadata={
            "source": "integration",
            "adata": "should_be_dropped",
            "name": "should_be_dropped",
            "results": {"bad": True},
        },
    )

    stored = await manager.get_dataset(data_id)
    assert data_id.startswith("derived_")
    assert stored["adata"] is minimal_spatial_adata
    assert stored["name"] == "derived-data"
    assert stored["source"] == "integration"
    assert "results" not in stored
    assert stored["type"] == "unknown"
    assert stored["n_cells"] == minimal_spatial_adata.n_obs
    assert stored["n_genes"] == minimal_spatial_adata.n_vars

    listed = await manager.list_datasets()
    assert listed == [
        {
            "id": data_id,
            "name": "derived-data",
            "type": "unknown",
            "n_cells": minimal_spatial_adata.n_obs,
            "n_genes": minimal_spatial_adata.n_vars,
        }
    ]


def test_tool_context_debug_and_log_config_delegate_to_logger():
    logger = Mock()
    ctx = adapter.ToolContext(
        _data_manager=adapter.DefaultSpatialDataManager(), _logger=logger
    )

    ctx.debug("hello")
    ctx.log_config("TestConfig", {"alpha": 1, "beta": "x"})

    assert logger.debug.call_count >= 3
    assert any("hello" in str(call.args[0]) for call in logger.debug.call_args_list)
    assert any(
        "TestConfig" in str(call.args[0]) for call in logger.debug.call_args_list
    )


def test_data_manager_metadata_detects_alternative_spatial_key(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["coordinates"] = adata.obsm.pop("spatial")

    metadata = adapter.DefaultSpatialDataManager._extract_adata_metadata(adata)

    assert metadata["spatial_coordinates_available"] is True
    assert metadata["obsm_keys"] == ["coordinates"]


def test_data_manager_metadata_requires_materialized_tissue_image(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.uns["spatial"] = {"sample": {"images": {"hires": None}}}

    metadata = adapter.DefaultSpatialDataManager._extract_adata_metadata(adata)

    assert metadata["tissue_image_available"] is False


@pytest.mark.asyncio
async def test_tool_context_data_access_and_add_dataset(minimal_spatial_adata):
    manager = adapter.DefaultSpatialDataManager()
    first_id = await manager.create_dataset(minimal_spatial_adata.copy(), prefix="ctx")

    ctx = adapter.ToolContext(_data_manager=manager)
    info = await ctx.get_dataset_info(first_id)
    assert info["adata"].n_obs == minimal_spatial_adata.n_obs

    subset = minimal_spatial_adata[:6, :5].copy()
    await ctx.set_adata(first_id, subset)
    out = await ctx.get_adata(first_id)
    assert out.shape == (6, 5)

    second_id = await ctx.add_dataset(
        minimal_spatial_adata.copy(),
        prefix="derived",
        name="derived_ctx",
        metadata={"source": "ctx"},
    )
    second = await manager.get_dataset(second_id)
    assert second_id.startswith("derived_")
    assert second["name"] == "derived_ctx"
    assert second["source"] == "ctx"


@pytest.mark.asyncio
async def test_tool_context_info_warning_error_use_mcp_context():
    ctx = adapter.ToolContext(
        _data_manager=adapter.DefaultSpatialDataManager(),
        _mcp_context=_FakeMCPContext(),
    )

    await ctx.info("i")
    await ctx.warning("w")
    await ctx.error("e")

    assert ctx._mcp_context is not None
    assert ctx._mcp_context.infos == ["i"]
    assert ctx._mcp_context.warnings == ["w"]
    assert ctx._mcp_context.errors == ["e"]
