"""Unit tests for MCP utility error handling behavior."""

from __future__ import annotations

import pytest

from chatspatial.models.analysis import PreprocessingResult
from chatspatial.utils.exceptions import ParameterError, ProcessingError
from chatspatial.utils.mcp_utils import mcp_tool_error_handler, suppress_output


@pytest.mark.asyncio
async def test_error_handler_for_str_return_type_reraises_user_error():
    @mcp_tool_error_handler()
    async def tool() -> str:
        raise ParameterError("invalid arg")

    with pytest.raises(ParameterError, match="invalid arg"):
        await tool()


@pytest.mark.asyncio
async def test_error_handler_for_simple_type_user_error_no_traceback():
    @mcp_tool_error_handler()
    async def tool() -> dict:
        raise ParameterError("bad input")

    with pytest.raises(ParameterError, match="bad input") as exc:
        await tool()
    assert "Traceback:" not in str(exc.value)


@pytest.mark.asyncio
async def test_error_handler_for_simple_type_non_user_error_has_traceback():
    @mcp_tool_error_handler()
    async def tool() -> dict:
        raise ProcessingError("compute failed")

    with pytest.raises(ProcessingError, match="compute failed") as exc:
        await tool()
    assert "Traceback:" in str(exc.value)


@pytest.mark.asyncio
async def test_error_handler_for_basemodel_reraises():
    @mcp_tool_error_handler()
    async def tool() -> PreprocessingResult:
        raise ProcessingError("must bubble up")

    with pytest.raises(ProcessingError, match="must bubble up"):
        await tool()


def test_suppress_output_discards_large_output(capsys: pytest.CaptureFixture[str]):
    with suppress_output():
        print("x" * 1_000_000)

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
