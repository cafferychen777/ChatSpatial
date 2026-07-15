"""Unit tests for runtime configuration utilities."""

from __future__ import annotations

import builtins
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from chatspatial import config as cfg


def test_is_inside_package_dir_with_explicit_and_cwd_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    package_root = tmp_path / "pkg"
    inside = package_root / "sub" / "work"
    outside = tmp_path / "elsewhere"
    inside.mkdir(parents=True, exist_ok=True)
    outside.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(cfg, "PACKAGE_ROOT", package_root)

    assert cfg.is_inside_package_dir(inside) is True
    assert cfg.is_inside_package_dir(outside) is False

    monkeypatch.chdir(inside)
    assert cfg.is_inside_package_dir() is True


def test_is_writable_dir_create_and_permission_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    writable = tmp_path / "writable"
    assert cfg._is_writable_dir(writable, create=True) is True
    assert cfg._is_writable_dir(tmp_path / "missing", create=False) is False

    def _deny_touch(self: Path, *args, **kwargs):
        raise PermissionError("denied")

    monkeypatch.setattr(Path, "touch", _deny_touch)
    assert cfg._is_writable_dir(writable) is False


def test_get_default_output_dir_prefers_safe_env_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    env_dir = tmp_path / "env_out"
    env_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv("CHATSPATIAL_OUTPUT_DIR", str(env_dir))

    monkeypatch.setattr(cfg, "is_inside_package_dir", lambda p=None: False)
    monkeypatch.setattr(cfg, "_is_writable_dir", lambda path, create=False: True)

    out = cfg.get_default_output_dir()
    assert out == env_dir.resolve()


def test_get_default_output_dir_falls_back_from_unsafe_env_to_cwd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    env_dir = tmp_path / "pkg" / "unsafe"
    env_dir.mkdir(parents=True, exist_ok=True)
    work = tmp_path / "work"
    work.mkdir(parents=True, exist_ok=True)

    monkeypatch.setenv("CHATSPATIAL_OUTPUT_DIR", str(env_dir))
    monkeypatch.chdir(work)

    monkeypatch.setattr(
        cfg,
        "is_inside_package_dir",
        lambda p=None: str((Path.cwd() if p is None else Path(p)).resolve()).startswith(
            str((tmp_path / "pkg").resolve())
        ),
    )
    monkeypatch.setattr(cfg, "_is_writable_dir", lambda path, create=False: True)

    out = cfg.get_default_output_dir()
    assert out == work


def test_get_default_output_dir_uses_tmp_as_last_resort(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("CHATSPATIAL_OUTPUT_DIR", raising=False)
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(cfg, "is_inside_package_dir", lambda p=None: True)

    def _writable(path: Path, create: bool = False) -> bool:
        return False

    monkeypatch.setattr(cfg, "_is_writable_dir", _writable)

    out = cfg.get_default_output_dir()
    assert out == Path("/tmp/chatspatial/outputs")
    assert out.exists()


def test_get_default_output_dir_uses_home_when_cwd_unusable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("CHATSPATIAL_OUTPUT_DIR", raising=False)
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(cfg, "is_inside_package_dir", lambda p=None: True)
    monkeypatch.setattr(
        cfg,
        "_is_writable_dir",
        lambda path, create=False: path == cfg.DEFAULT_OUTPUT_DIR and create,
    )

    out = cfg.get_default_output_dir()
    assert out == cfg.DEFAULT_OUTPUT_DIR


def test_configure_environment_sets_required_vars(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("TQDM_DISABLE", raising=False)
    monkeypatch.setenv("DASK_DATAFRAME__QUERY_PLANNING", "False")

    cfg._configure_environment()

    assert os_environ("TQDM_DISABLE") == "1"
    assert os_environ("DASK_DATAFRAME__QUERY_PLANNING") == "False"


def test_configure_warnings_registers_known_filters(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple] = []

    def _capture(*args, **kwargs):
        calls.append((args, kwargs))

    monkeypatch.setattr(cfg.warnings, "filterwarnings", _capture)

    cfg._configure_warnings()

    broad_filters = [
        kwargs
        for _, kwargs in calls
        if kwargs.get("category") in {UserWarning, FutureWarning}
        and "message" not in kwargs
    ]
    assert broad_filters == []
    assert any(
        kwargs.get("message") == "The legacy Dask DataFrame implementation is deprecated"
        for _, kwargs in calls
    )


def test_configure_libraries_handles_missing_dask_and_sets_scanpy(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _ScanpySettings:
        verbosity = 99
        n_jobs = 1

    fake_scanpy = SimpleNamespace(settings=_ScanpySettings())
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.delitem(sys.modules, "dask", raising=False)

    original_import = builtins.__import__

    def _import(name, *args, **kwargs):
        if name == "dask":
            raise ImportError("dask missing")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _import)

    cfg._configure_libraries()

    assert fake_scanpy.settings.verbosity == 0
    assert fake_scanpy.settings.n_jobs == -1


def test_configure_libraries_configures_dask_when_available(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dask_calls: list[dict] = []

    class _FakeDaskConfig:
        @staticmethod
        def set(payload: dict) -> None:
            dask_calls.append(payload)

    fake_dask = SimpleNamespace(config=_FakeDaskConfig)
    fake_scanpy = SimpleNamespace(settings=SimpleNamespace(verbosity=5, n_jobs=2))

    monkeypatch.setitem(sys.modules, "dask", fake_dask)
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    cfg._configure_libraries()

    assert dask_calls == [{"dataframe.query-planning": True}]
    assert fake_scanpy.settings.verbosity == 0
    assert fake_scanpy.settings.n_jobs == -1


def test_init_runtime_is_idempotent_but_verbose_always_prints(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    calls = {"env": 0, "warn": 0, "lib": 0}

    monkeypatch.setattr(cfg, "_initialized", False)
    monkeypatch.setattr(cfg, "_configure_environment", lambda: calls.__setitem__("env", calls["env"] + 1))
    monkeypatch.setattr(cfg, "_configure_warnings", lambda: calls.__setitem__("warn", calls["warn"] + 1))
    monkeypatch.setattr(cfg, "_configure_libraries", lambda: calls.__setitem__("lib", calls["lib"] + 1))
    monkeypatch.setattr(cfg, "get_default_output_dir", lambda: Path("/tmp/chatspatial/outputs"))

    cfg.init_runtime(verbose=True)
    first = capsys.readouterr().err
    cfg.init_runtime(verbose=True)
    second = capsys.readouterr().err

    assert calls == {"env": 1, "warn": 1, "lib": 1}
    assert "ChatSpatial initialized:" in first
    assert "ChatSpatial initialized:" in second


def os_environ(key: str) -> str | None:
    import os

    return os.environ.get(key)
