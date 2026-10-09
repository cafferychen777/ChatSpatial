"""Client-visible error contract over the real MCP protocol path.

MCP SDK 2.1.0 started replacing the text of every tool exception that is not a
``ToolError`` with ``Error executing tool <name>``. These tests drive an
in-process client against ChatSpatial servers and pin the contract for every
installed ``mcp`` 2.x release: anticipated ChatSpatial errors reach the client
with their type and full message, and unexpected failures never expose a
traceback.
"""

from __future__ import annotations

import sys
import types

import pytest
from mcp.client import Client
from mcp.server.mcpserver import exceptions as mcp_exceptions

from chatspatial.server import mcp
from chatspatial.spatial_mcp_adapter import create_spatial_mcp_server
from chatspatial.utils.exceptions import (
    DataCompatibilityError,
    DataError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)

RAW_COUNTS_MESSAGE = (
    "No raw integer counts found. Sources tried: ['current (normalized)']. "
    "Data appears to be normalized (has_negatives=False, has_decimals=True). "
    "Deconvolution and velocity methods require raw integer counts."
)

# SDK releases that withhold the text of unanticipated tool exceptions.
SDK_MASKS_UNEXPECTED_ERRORS = hasattr(mcp_exceptions, "UnexpectedToolError")


def _error_text(result) -> str:
    assert result.is_error is True
    assert len(result.content) == 1
    return result.content[0].text


@pytest.mark.integration
@pytest.mark.asyncio
async def test_data_error_from_real_tool_reaches_client_verbatim(
    monkeypatch: pytest.MonkeyPatch,
):
    """Reproduce the deconvolve_data report: the DataError text must survive."""

    async def raise_raw_counts_error(*_args, **_kwargs):
        raise DataError(RAW_COUNTS_MESSAGE)

    # The tool imports its backend lazily; a stub keeps this test free of the
    # deep-learning stack while the real tool wrapper and server stay in place.
    stub = types.ModuleType("chatspatial.tools.deconvolution")
    stub.deconvolve_spatial_data = raise_raw_counts_error
    monkeypatch.setitem(sys.modules, "chatspatial.tools.deconvolution", stub)

    async with Client(mcp) as client:
        result = await client.call_tool(
            "deconvolve_data",
            {
                "data_id": "spatial",
                "params": {
                    "method": "flashdeconv",
                    "reference_data_id": "reference",
                    "cell_type_key": "cell_type",
                },
            },
        )

    text = _error_text(result)
    assert f"DataError: {RAW_COUNTS_MESSAGE}" in text
    assert "Traceback" not in text


@pytest.mark.integration
@pytest.mark.asyncio
async def test_missing_dataset_error_reaches_client_verbatim(reset_data_manager):
    async with Client(mcp) as client:
        result = await client.call_tool(
            "preprocess_data",
            {"data_id": "absent_dataset"},
        )

    text = _error_text(result)
    assert "DataNotFoundError: Dataset absent_dataset not found" in text
    assert "Traceback" not in text


@pytest.mark.integration
@pytest.mark.asyncio
@pytest.mark.parametrize(
    "error",
    [
        DataError("adata.obs lacks column 'cell_type'"),
        DataNotFoundError("Deconvolution results not found"),
        DataCompatibilityError("Species mismatch: mouse vs human"),
        ParameterError("n_clusters must be > 0"),
        ProcessingError("Leiden clustering failed to converge"),
        DependencyError("scvi-tools required: pip install scvi-tools"),
        ValueError("Method requires raw count data (integers)"),
        FileNotFoundError("Data path not found: missing.h5ad"),
        PermissionError("Cannot write to output directory: results"),
    ],
    ids=lambda error: type(error).__name__,
)
async def test_anticipated_errors_reach_client_with_type_and_message(
    error: Exception,
):
    server, _adapter = create_spatial_mcp_server("error-visibility")

    @server.tool()
    async def fail() -> str:
        raise error

    async with Client(server) as client:
        result = await client.call_tool("fail", {})

    text = _error_text(result)
    assert f"{type(error).__name__}: {error}" in text
    assert "Traceback" not in text


@pytest.mark.integration
@pytest.mark.asyncio
async def test_direct_calls_keep_chatspatial_exception_types():
    server, _adapter = create_spatial_mcp_server("error-visibility")

    @server.tool()
    async def fail() -> str:
        raise DataError("direct callers see the domain exception")

    with pytest.raises(DataError, match="direct callers see the domain exception"):
        await fail()


@pytest.mark.integration
@pytest.mark.asyncio
async def test_unexpected_error_does_not_leak_traceback_or_private_detail():
    server, _adapter = create_spatial_mcp_server("error-visibility")

    @server.tool()
    async def crash() -> str:
        raise RuntimeError("internal state at /Users/someone/private/cache")

    async with Client(server) as client:
        result = await client.call_tool("crash", {})

    text = _error_text(result)
    assert text.startswith("Error executing tool crash")
    assert "Traceback" not in text
    if SDK_MASKS_UNEXPECTED_ERRORS:
        assert "/Users/" not in text
        assert "RuntimeError" not in text
