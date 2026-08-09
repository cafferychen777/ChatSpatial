"""Unit tests for compute device selection utilities."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chatspatial.utils import device_utils as du


class _FakeCtx:
    def __init__(self) -> None:
        self.messages: list[str] = []

    async def warning(self, msg: str) -> None:
        self.messages.append(msg)


@pytest.fixture
def reset_mps_state(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(du, "_mps_configured", False)
    monkeypatch.delenv("PYTORCH_MPS_HIGH_WATERMARK_RATIO", raising=False)


def _install_fake_torch(
    monkeypatch: pytest.MonkeyPatch,
    *,
    cuda_available: bool = False,
    mps_available: bool | None = None,
) -> None:
    backends = SimpleNamespace()
    if mps_available is not None:
        backends.mps = SimpleNamespace(is_available=lambda: mps_available)

    fake_torch = SimpleNamespace(
        cuda=SimpleNamespace(is_available=lambda: cuda_available),
        backends=backends,
    )
    monkeypatch.setattr(
        du,
        "get",
        lambda name: fake_torch if name == "torch" else None,
    )


def test_cuda_available_returns_true_when_torch_cuda_available(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_torch(monkeypatch, cuda_available=True)
    assert du.cuda_available() is True


def test_cuda_and_mps_available_return_false_when_torch_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(du, "get", lambda _name: None)

    assert du.cuda_available() is False
    assert du.mps_available() is False


def test_mps_available_requires_mps_backend(monkeypatch: pytest.MonkeyPatch) -> None:
    _install_fake_torch(monkeypatch, mps_available=True)
    assert du.mps_available() is True

    _install_fake_torch(monkeypatch, mps_available=None)
    assert du.mps_available() is False


def test_get_device_prefers_cuda_without_mps_configuration(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(du, "cuda_available", lambda: True)
    monkeypatch.setattr(du, "mps_available", lambda: True)

    called = {"mps_config": 0}

    def _configure_once() -> None:
        called["mps_config"] += 1

    monkeypatch.setattr(du, "_configure_mps", _configure_once)

    assert du.get_device(prefer_gpu=True) == "cuda:0"
    assert called["mps_config"] == 0


def test_get_device_uses_mps_fallback_and_can_disable_it(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(du, "cuda_available", lambda: False)
    monkeypatch.setattr(du, "mps_available", lambda: True)

    called = {"mps_config": 0}

    def _configure_once() -> None:
        called["mps_config"] += 1

    monkeypatch.setattr(du, "_configure_mps", _configure_once)

    assert du.get_device(prefer_gpu=True, allow_mps=True) == "mps"
    assert called["mps_config"] == 1
    assert du.get_device(prefer_gpu=True, allow_mps=False) == "cpu"


def test_configure_mps_is_idempotent_and_preserves_existing_env(
    monkeypatch: pytest.MonkeyPatch,
    reset_mps_state: None,
) -> None:
    du._configure_mps()
    assert du._mps_configured is True
    assert os_environ_value("PYTORCH_MPS_HIGH_WATERMARK_RATIO") == "0.0"

    monkeypatch.setattr(du, "_mps_configured", False)
    monkeypatch.setenv("PYTORCH_MPS_HIGH_WATERMARK_RATIO", "0.7")
    du._configure_mps()
    assert os_environ_value("PYTORCH_MPS_HIGH_WATERMARK_RATIO") == "0.7"


def test_configure_mps_returns_early_when_already_configured(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(du, "_mps_configured", True)
    monkeypatch.setenv("PYTORCH_MPS_HIGH_WATERMARK_RATIO", "0.9")

    du._configure_mps()

    assert os_environ_value("PYTORCH_MPS_HIGH_WATERMARK_RATIO") == "0.9"


def os_environ_value(key: str) -> str | None:
    import os

    return os.environ.get(key)


def test_get_accelerator_returns_gpu_for_mps_when_requested(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(du, "cuda_available", lambda: False)
    monkeypatch.setattr(du, "mps_available", lambda: True)

    called = {"mps_config": 0}

    def _configure_once() -> None:
        called["mps_config"] += 1

    monkeypatch.setattr(du, "_configure_mps", _configure_once)

    assert du.get_accelerator(prefer_gpu=True) == "gpu"
    assert called["mps_config"] == 1
    assert du.get_accelerator(prefer_gpu=False) == "cpu"


def test_get_accelerator_prefers_cuda_without_mps_configuration(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(du, "cuda_available", lambda: True)
    monkeypatch.setattr(du, "mps_available", lambda: True)
    called = {"mps_config": 0}

    def _configure_once() -> None:
        called["mps_config"] += 1

    monkeypatch.setattr(du, "_configure_mps", _configure_once)
    assert du.get_accelerator(prefer_gpu=True) == "gpu"
    assert called["mps_config"] == 0


@pytest.mark.asyncio
async def test_resolve_device_async_warns_on_cpu_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ctx = _FakeCtx()
    monkeypatch.setattr(du, "get_device", lambda **_kwargs: "cpu")

    device = await du.resolve_device_async(prefer_gpu=True, ctx=ctx)

    assert device == "cpu"
    assert ctx.messages == ["GPU requested but not available - using CPU"]


@pytest.mark.asyncio
async def test_resolve_device_async_can_suppress_warning(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ctx = _FakeCtx()
    monkeypatch.setattr(du, "get_device", lambda **_kwargs: "cpu")

    device = await du.resolve_device_async(
        prefer_gpu=True,
        ctx=ctx,
        warn_on_fallback=False,
    )

    assert device == "cpu"
    assert ctx.messages == []


def test_get_ot_backend_selects_torch_backend_when_cuda_available(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _FakeBackend:
        class TorchBackend:
            def __init__(self) -> None:
                self.name = "torch"

        class NumpyBackend:
            def __init__(self) -> None:
                self.name = "numpy"

    fake_ot = SimpleNamespace(backend=_FakeBackend)
    monkeypatch.setattr(
        du,
        "require",
        lambda name, **_kwargs: (
            fake_ot if name == "ot" else pytest.fail(f"Unexpected dependency: {name}")
        ),
    )
    monkeypatch.setattr(du, "cuda_available", lambda: True)

    backend = du.get_ot_backend(use_gpu=True)
    assert backend.name == "torch"

    monkeypatch.setattr(du, "cuda_available", lambda: False)
    backend = du.get_ot_backend(use_gpu=True)
    assert backend.name == "numpy"


def test_get_ot_backend_uses_numpy_when_gpu_not_requested(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _FakeBackend:
        class TorchBackend:
            def __init__(self) -> None:
                self.name = "torch"

        class NumpyBackend:
            def __init__(self) -> None:
                self.name = "numpy"

    fake_ot = SimpleNamespace(backend=_FakeBackend)
    monkeypatch.setattr(
        du,
        "require",
        lambda name, **_kwargs: (
            fake_ot if name == "ot" else pytest.fail(f"Unexpected dependency: {name}")
        ),
    )
    monkeypatch.setattr(du, "cuda_available", lambda: True)

    backend = du.get_ot_backend(use_gpu=False)
    assert backend.name == "numpy"
