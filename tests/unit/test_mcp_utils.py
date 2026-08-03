"""Unit tests for MCP utility error handling behavior."""

from __future__ import annotations

import asyncio
import logging
import sys

import pytest

from chatspatial.models.analysis import PreprocessingResult
from chatspatial.utils.exceptions import (
    DependencyError,
    ParameterError,
    ProcessingError,
)
from chatspatial.utils.mcp_utils import (
    mcp_tool_error_handler,
    suppress_output,
    suppress_output_async,
)


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
async def test_error_handler_preserves_dependency_error_without_traceback():
    @mcp_tool_error_handler()
    async def tool() -> dict:
        raise DependencyError("install optional backend")

    with pytest.raises(DependencyError, match="install optional backend") as exc:
        await tool()
    assert "Traceback:" not in str(exc.value)


@pytest.mark.asyncio
async def test_error_handler_logs_non_user_error_without_exposing_traceback(
    caplog: pytest.LogCaptureFixture,
):
    @mcp_tool_error_handler()
    async def tool() -> dict:
        raise ProcessingError("compute failed")

    with (
        caplog.at_level(logging.ERROR),
        pytest.raises(ProcessingError, match="compute failed") as exc,
    ):
        await tool()

    assert "Traceback:" not in str(exc.value)
    assert "MCP tool tool failed" in caplog.text
    assert "ProcessingError: compute failed" in caplog.text


@pytest.mark.asyncio
async def test_error_handler_can_expose_traceback_for_local_debugging():
    @mcp_tool_error_handler(include_traceback=True)
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


def test_suppress_output_restores_streams_and_preserves_logging_level(
    capsys: pytest.CaptureFixture[str],
):
    root_logger = logging.getLogger()
    original_level = root_logger.level
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    with pytest.raises(RuntimeError, match="failed"), suppress_output():
        print("hidden stdout")
        print("hidden stderr", file=sys.stderr)
        raise RuntimeError("failed")

    assert root_logger.level == original_level
    assert sys.stdout is original_stdout
    assert sys.stderr is original_stderr
    captured = capsys.readouterr()
    assert "hidden" not in captured.out
    assert "hidden" not in captured.err


@pytest.mark.asyncio
async def test_async_suppression_is_task_local(
    capsys: pytest.CaptureFixture[str],
) -> None:
    suppression_active = asyncio.Event()
    visible_output_written = asyncio.Event()

    async def hidden_task() -> None:
        async with suppress_output_async():
            print("hidden before await")
            suppression_active.set()
            await visible_output_written.wait()
            print("hidden after await", file=sys.stderr)

    async def visible_task() -> None:
        await suppression_active.wait()
        print("visible stdout")
        print("visible stderr", file=sys.stderr)
        visible_output_written.set()

    await asyncio.gather(hidden_task(), visible_task())

    captured = capsys.readouterr()
    assert captured.out == "visible stdout\n"
    assert captured.err == "visible stderr\n"


@pytest.mark.asyncio
async def test_async_suppression_propagates_through_to_thread(
    capsys: pytest.CaptureFixture[str],
) -> None:
    async with suppress_output_async():
        await asyncio.to_thread(print, "hidden worker output")

    assert "hidden worker output" not in capsys.readouterr().out
