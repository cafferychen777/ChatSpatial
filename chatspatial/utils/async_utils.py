"""Async boundaries for synchronous scientific-computing operations."""

from __future__ import annotations

import asyncio
from collections.abc import Callable
from contextlib import suppress
from multiprocessing.connection import Connection
from typing import ParamSpec, TypeVar

import cloudpickle

P = ParamSpec("P")
T = TypeVar("T")


async def run_sync(function: Callable[P, T], /, *args: P.args, **kwargs: P.kwargs) -> T:
    """Run blocking work without occupying the MCP event-loop thread."""
    return await asyncio.to_thread(function, *args, **kwargs)


def _run_pickled_callable(serialized: bytes, sender: Connection) -> None:
    """Execute one cloudpickled callable and return a serialized outcome."""
    try:
        function = cloudpickle.loads(serialized)
        payload = cloudpickle.dumps((True, function(), None))
    except BaseException as exc:  # Child must report backend failures to parent.
        import traceback

        remote_traceback = traceback.format_exc()
        try:
            payload = cloudpickle.dumps((False, exc, remote_traceback))
        except Exception:
            payload = cloudpickle.dumps(
                (
                    False,
                    RuntimeError(f"Worker failed with {type(exc).__name__}: {exc}"),
                    remote_traceback,
                )
            )

    try:
        sender.send_bytes(payload)
    finally:
        sender.close()


async def run_sync_with_timeout(
    function: Callable[[], T],
    /,
    *,
    timeout: float,
    process_name: str = "chatspatial-worker",
) -> T:
    """Run blocking work in a disposable process with a hard timeout."""
    if timeout <= 0:
        raise ValueError("timeout must be greater than zero")

    import multiprocessing

    serialized = await asyncio.to_thread(cloudpickle.dumps, function)
    context = multiprocessing.get_context("spawn")
    receiver, sender = context.Pipe(duplex=False)
    process = context.Process(
        target=_run_pickled_callable,
        args=(serialized, sender),
        name=process_name,
    )
    try:
        process.start()
    except BaseException:
        sender.close()
        receiver.close()
        with suppress(ValueError):
            process.close()
        raise
    sender.close()
    receive_task = asyncio.create_task(asyncio.to_thread(receiver.recv_bytes))

    try:
        payload = await asyncio.wait_for(receive_task, timeout=timeout)
        succeeded, value, remote_traceback = cloudpickle.loads(payload)
        if succeeded:
            return value
        if isinstance(value, BaseException):
            if remote_traceback and hasattr(value, "add_note"):
                value.add_note(f"Worker traceback:\n{remote_traceback}")
            raise value
        raise RuntimeError(f"Worker returned an invalid error payload: {value!r}")
    except (TimeoutError, asyncio.CancelledError):
        if process.is_alive():
            process.terminate()
        raise
    finally:
        receiver.close()
        if not receive_task.done():
            receive_task.cancel()
        with suppress(asyncio.CancelledError, EOFError, OSError):
            await receive_task
        await asyncio.to_thread(process.join, 1.0)
        if process.is_alive():
            process.kill()
            await asyncio.to_thread(process.join)
        process.close()
