"""Unit tests for safe output path resolution."""

from __future__ import annotations

from pathlib import Path

import pytest

from chatspatial.utils import path_utils


@pytest.fixture
def fake_default_dir(tmp_path: Path) -> Path:
    d = tmp_path / "safe_outputs"
    d.mkdir(parents=True, exist_ok=True)
    return d


def test_relative_path_outside_package_resolves_to_cwd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    cwd = tmp_path / "work"
    cwd.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(cwd)

    monkeypatch.setattr(path_utils, "is_inside_package_dir", lambda p=None: False)

    out = path_utils.get_safe_output_path("./viz")
    assert out == cwd / "viz"
    assert out.exists()


def test_absolute_path_inside_package_redirects_to_safe_default(
    tmp_path: Path, fake_default_dir: Path, monkeypatch: pytest.MonkeyPatch
):
    package_root = tmp_path / "pkg"
    package_root.mkdir(parents=True, exist_ok=True)
    inside = package_root / "subdir" / "plots"

    monkeypatch.setattr(path_utils, "PACKAGE_ROOT", package_root)
    monkeypatch.setattr(path_utils, "get_default_output_dir", lambda: fake_default_dir)
    monkeypatch.setattr(
        path_utils,
        "is_inside_package_dir",
        lambda p=None: (Path.cwd() if p is None else Path(p))
        .resolve()
        .is_relative_to(package_root.resolve()),
    )

    with pytest.warns(UserWarning, match="inside package directory"):
        out = path_utils.get_safe_output_path(str(inside))

    assert out == fake_default_dir / "subdir" / "plots"
    assert out.exists()


def test_permission_failure_fallback_to_tmp(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    cwd = tmp_path / "work"
    cwd.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(cwd)
    monkeypatch.setattr(path_utils, "is_inside_package_dir", lambda p=None: False)

    fallback = tmp_path / "platform-temp" / "chatspatial" / "outputs"
    fallback.mkdir(parents=True)
    monkeypatch.setattr(path_utils, "is_writable_dir", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(path_utils, "get_temp_output_dir", lambda: fallback)

    with pytest.warns(UserWarning, match="Falling back to"):
        out = path_utils.get_safe_output_path("./cannot_write")

    assert out == fallback
    assert out.exists()


def test_permission_failure_without_fallback_raises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    cwd = tmp_path / "work"
    cwd.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(cwd)
    monkeypatch.setattr(path_utils, "is_inside_package_dir", lambda p=None: False)

    monkeypatch.setattr(path_utils, "is_writable_dir", lambda *_args, **_kwargs: False)

    with pytest.raises(PermissionError, match="Cannot write to output directory"):
        path_utils.get_safe_output_path("./cannot_write", fallback_to_tmp=False)
