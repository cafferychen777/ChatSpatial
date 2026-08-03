"""Unit contracts for spatial MCP adapter primitives."""

from __future__ import annotations

import asyncio
from unittest.mock import Mock

import pytest

from chatspatial import spatial_mcp_adapter as adapter
from chatspatial.models.analysis import PreprocessingResult
from chatspatial.utils.exceptions import DataNotFoundError


class _FakeMCPContext:
    def __init__(self):
        self.progress_updates: list[tuple[float, float | None, str | None]] = []

    async def report_progress(
        self,
        progress: float,
        total: float | None = None,
        message: str | None = None,
    ) -> None:
        self.progress_updates.append((progress, total, message))


class _FailingMCPContext:
    async def report_progress(self, *args, **kwargs) -> None:
        raise RuntimeError("transport unavailable")


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
async def test_data_manager_update_missing_dataset_raises(
    minimal_spatial_adata,
):
    manager = adapter.DefaultSpatialDataManager()

    with pytest.raises(DataNotFoundError, match="Dataset missing not found"):
        await manager.update_adata("missing", minimal_spatial_adata.copy())


@pytest.mark.asyncio
async def test_data_manager_create_dataset_protects_ownership_metadata(
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
    assert stored["results"] == {"bad": True}
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


@pytest.mark.asyncio
async def test_progress_delivery_failure_does_not_fail_tool_context():
    logger = Mock()
    ctx = adapter.ToolContext(
        _data_manager=adapter.DefaultSpatialDataManager(),
        _mcp_context=_FailingMCPContext(),
        _logger=logger,
    )

    await ctx.info("still valid")

    logger.warning.assert_called_once()
    assert ctx._progress_step == 1


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
async def test_tool_context_uses_progress_and_structured_warnings():
    logger = Mock()
    ctx = adapter.ToolContext(
        _data_manager=adapter.DefaultSpatialDataManager(),
        _mcp_context=_FakeMCPContext(),
        _logger=logger,
    )

    await ctx.info("i")
    await ctx.warning("w")
    await ctx.error("e")

    assert ctx._mcp_context is not None
    assert ctx._mcp_context.progress_updates == [
        (1.0, None, "i"),
        (2.0, None, "Warning: w"),
        (3.0, None, "Error: e"),
    ]
    assert ctx.warnings == ("w",)
    assert ctx.finalize({"ok": True}) == {"ok": True, "warnings": ["w"]}
    model_result = PreprocessingResult(
        data_id="d1",
        n_cells=10,
        n_genes=5,
        n_hvgs=3,
        clusters=2,
    )
    finalized_model = ctx.finalize(model_result)
    assert finalized_model is not model_result
    assert finalized_model.warnings == ["w"]
    assert logger.info.called
    assert logger.warning.called
    assert logger.error.called


@pytest.mark.asyncio
async def test_serialize_dataset_access_prevents_overlapping_calls():
    manager = adapter.DefaultSpatialDataManager()
    active_calls = 0
    maximum_active_calls = 0

    @adapter.serialize_dataset_access(manager, "data_id")
    async def operation(data_id: str) -> None:
        nonlocal active_calls, maximum_active_calls
        active_calls += 1
        maximum_active_calls = max(maximum_active_calls, active_calls)
        await asyncio.sleep(0)
        active_calls -= 1

    await asyncio.gather(operation("d1"), operation("d1"))

    assert maximum_active_calls == 1
    assert manager._dataset_locks == {}


@pytest.mark.asyncio
async def test_serialize_dataset_access_releases_locks_after_failure():
    manager = adapter.DefaultSpatialDataManager()

    @adapter.serialize_dataset_access(manager, "data_id")
    async def operation(data_id: str) -> None:
        await manager.get_dataset(data_id)

    for index in range(100):
        with pytest.raises(DataNotFoundError):
            await operation(f"missing_{index}")

    assert manager._dataset_locks == {}


@pytest.mark.asyncio
async def test_serialize_dataset_access_releases_locks_after_cancellation():
    manager = adapter.DefaultSpatialDataManager()
    started = asyncio.Event()
    release = asyncio.Event()

    @adapter.serialize_dataset_access(manager, "data_id")
    async def operation(data_id: str) -> None:
        started.set()
        await release.wait()

    task = asyncio.create_task(operation("d1"))
    await started.wait()
    task.cancel()

    with pytest.raises(asyncio.CancelledError):
        await task

    assert manager._dataset_locks == {}


def test_server_factory_preserves_positional_data_manager_contract():
    manager = adapter.DefaultSpatialDataManager()

    server, spatial_adapter = adapter.create_spatial_mcp_server("Custom", manager)

    assert spatial_adapter.data_manager is manager
    assert server.version == adapter.__version__


def test_server_factory_accepts_explicit_version_metadata():
    server, _ = adapter.create_spatial_mcp_server(
        "Custom",
        server_version="9.8.7",
    )

    assert server.version == "9.8.7"


def test_data_manager_reset_clears_data_locks_and_ids(minimal_spatial_adata):
    manager = adapter.DefaultSpatialDataManager()
    manager.data_store["d1"] = {"adata": minimal_spatial_adata}
    manager._dataset_locks["d1"] = Mock()

    manager.reset()

    assert manager.data_store == {}
    assert manager._dataset_locks == {}
    assert next(manager._id_counter) == 1
