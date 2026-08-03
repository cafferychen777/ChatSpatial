"""End-to-end contracts for modern and legacy MCP protocol eras."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest
from mcp.client import Client
from mcp.server import MCPServer

from chatspatial import __version__
from chatspatial.server import mcp
from chatspatial.utils.exceptions import ProcessingError
from chatspatial.utils.mcp_utils import mcp_tool_error_handler


@pytest.mark.integration
@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("mode", "expected_protocol"),
    [("auto", "2026-07-28"), ("legacy", "2025-11-25")],
)
async def test_server_discovery_and_tool_listing_across_protocol_eras(
    mode: str,
    expected_protocol: str,
):
    async with Client(mcp, mode=mode) as client:
        result = await client.list_tools()

        assert client.protocol_version == expected_protocol
        assert client.server_info is not None
        assert client.server_info.name == "ChatSpatial"
        assert client.server_info.title == "ChatSpatial"
        assert client.server_info.version == __version__
        assert result.result_type == "complete"
        assert len(result.tools) == 20
        assert [tool.name for tool in result.tools] == [
            tool.name for tool in await mcp.list_tools()
        ]

        if mode == "auto":
            assert result.ttl_ms == 3_600_000
            assert result.cache_scope == "public"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_modern_protocol_tool_call_returns_structured_content_and_progress(
    spatial_dataset_path,
    reset_data_manager,
):
    progress_updates: list[tuple[float, float | None, str | None]] = []

    async def capture_progress(
        progress: float,
        total: float | None,
        message: str | None,
    ) -> None:
        progress_updates.append((progress, total, message))

    async with Client(mcp) as client:
        result = await client.call_tool(
            "load_data",
            {
                "data_path": str(spatial_dataset_path),
                "data_type": "generic",
                "name": "protocol-v2-contract",
            },
            progress_callback=capture_progress,
        )

    assert result.is_error is False
    assert result.result_type == "complete"
    assert result.structured_content is not None
    assert result.structured_content["name"] == "protocol-v2-contract"
    assert result.structured_content["data_type"] == "generic"
    assert result.structured_content["n_cells"] > 0
    assert result.structured_content["n_genes"] > 0
    assert [update[0] for update in progress_updates] == [1.0, 2.0]
    assert all(update[1] is None for update in progress_updates)
    assert "Loading data" in (progress_updates[0][2] or "")
    assert "Successfully loaded" in (progress_updates[1][2] or "")


@pytest.mark.integration
@pytest.mark.asyncio
async def test_load_data_keeps_non_finite_column_profiles_schema_valid(
    tmp_path,
    reset_data_manager,
):
    adata = ad.AnnData(np.ones((4, 3)))
    adata.obs["all_nan"] = np.nan
    adata.obsm["spatial"] = np.arange(8, dtype=float).reshape(4, 2)
    data_path = tmp_path / "all_nan_metadata.h5ad"
    adata.write_h5ad(data_path)

    async with Client(mcp) as client:
        result = await client.call_tool(
            "load_data",
            {
                "data_path": str(data_path),
                "data_type": "generic",
                "name": "all-nan-metadata",
            },
        )

    assert result.is_error is False
    column = next(
        item
        for item in result.structured_content["obs_columns"]
        if item["name"] == "all_nan"
    )
    assert column["range"] is None
    assert column["n_unique"] == 0


@pytest.mark.integration
@pytest.mark.asyncio
async def test_unexpected_tool_failure_does_not_expose_traceback_over_mcp():
    test_server = MCPServer("error-boundary-contract")

    @test_server.tool()
    @mcp_tool_error_handler()
    async def fail_processing() -> dict[str, str]:
        raise ProcessingError("controlled wire failure")

    async with Client(test_server) as client:
        result = await client.call_tool("fail_processing", {})

    assert result.is_error is True
    assert len(result.content) == 1
    response_text = result.content[0].text
    assert "controlled wire failure" in response_text
    assert "Traceback:" not in response_text
    assert "/Users/" not in response_text


@pytest.mark.integration
@pytest.mark.asyncio
async def test_obsolete_preprocessing_parameter_is_rejected_over_mcp():
    async with Client(mcp) as client:
        result = await client.call_tool(
            "preprocess_data",
            {"data_id": "unused", "params": {"n_pcs": 30}},
        )

    assert result.is_error is True
    response_text = result.content[0].text
    assert "params.n_pcs" in response_text
    assert "Extra inputs are not permitted" in response_text


@pytest.mark.integration
@pytest.mark.asyncio
async def test_unknown_top_level_parameter_is_rejected_over_mcp():
    async with Client(mcp) as client:
        result = await client.call_tool(
            "export_data",
            {"data_id": "unused", "data_idx": "misspelled"},
        )

    assert result.is_error is True
    response_text = result.content[0].text
    assert "data_idx" in response_text
    assert "Extra inputs are not permitted" in response_text
