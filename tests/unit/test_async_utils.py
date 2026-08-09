"""Tests for shared async boundaries around blocking scientific code."""

from __future__ import annotations

import asyncio
import contextvars
import time
from pathlib import Path

import pytest

from chatspatial.utils.async_utils import run_sync, run_sync_with_timeout


@pytest.mark.asyncio
async def test_run_sync_keeps_event_loop_responsive_and_propagates_context():
    marker = contextvars.ContextVar("marker", default="missing")
    marker.set("request-context")
    event_loop_progressed = asyncio.Event()

    async def mark_progress() -> None:
        await asyncio.sleep(0.01)
        event_loop_progressed.set()

    progress_task = asyncio.create_task(mark_progress())
    observed = await run_sync(lambda: (time.sleep(0.03), marker.get())[1])
    await progress_task

    assert event_loop_progressed.is_set()
    assert observed == "request-context"


@pytest.mark.asyncio
async def test_run_sync_with_timeout_terminates_worker_on_timeout(tmp_path: Path):
    completion_marker = tmp_path / "completed"

    def slow_operation() -> None:
        time.sleep(0.5)
        completion_marker.write_text("completed")

    started = time.perf_counter()
    with pytest.raises(TimeoutError):
        await run_sync_with_timeout(slow_operation, timeout=0.05)
    elapsed = time.perf_counter() - started

    assert elapsed < 0.4
    await asyncio.sleep(0.55)
    assert not completion_marker.exists()


@pytest.mark.asyncio
async def test_run_sync_with_timeout_returns_result_and_validates_timeout():
    assert await run_sync_with_timeout(lambda: 42, timeout=1) == 42

    with pytest.raises(ValueError, match="greater than zero"):
        await run_sync_with_timeout(lambda: None, timeout=0)
