"""Branch-focused contracts for safe output path edge handling."""

from __future__ import annotations

from pathlib import Path

import pytest

from chatspatial.utils import path_utils


def test_absolute_inside_package_with_relative_to_fallback_name(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    package_root = tmp_path / "pkg_root"
    package_root.mkdir()
    safe_default = tmp_path / "safe"
    safe_default.mkdir()
    user_abs = tmp_path / "outside" / "plots"

    monkeypatch.setattr(path_utils, "PACKAGE_ROOT", package_root)
    monkeypatch.setattr(path_utils, "get_default_output_dir", lambda: safe_default)
    # Force "inside package" decision even though relative_to() will fail.
    monkeypatch.setattr(path_utils, "is_inside_package_dir", lambda p=None: True)

    with pytest.warns(UserWarning, match="inside package directory"):
        out = path_utils.get_safe_output_path(str(user_abs))

    assert out == safe_default / "plots"
    assert out.exists()


def test_relative_path_inside_package_uses_safe_default_base(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    safe_default = tmp_path / "safe"
    safe_default.mkdir()
    package_root = tmp_path / "pkg"
    cwd_inside = package_root / "src"
    cwd_inside.mkdir(parents=True)

    monkeypatch.setattr(path_utils, "PACKAGE_ROOT", package_root)
    monkeypatch.setattr(path_utils, "get_default_output_dir", lambda: safe_default)
    monkeypatch.setattr(
        path_utils,
        "is_inside_package_dir",
        lambda p=None: (Path.cwd() if p is None else Path(p))
        .resolve()
        .is_relative_to(package_root.resolve()),
    )
    monkeypatch.chdir(cwd_inside)

    out = path_utils.get_safe_output_path("./a/plots")
    assert out == safe_default / "a" / "plots"
    assert out.exists()


def test_double_check_branch_redirects_target_inside_package(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    safe_default = tmp_path / "safe"
    safe_default.mkdir()
    package_root = tmp_path / "pkg_root"
    package_root.mkdir()
    cwd = tmp_path / "work"
    cwd.mkdir()
    monkeypatch.chdir(cwd)

    monkeypatch.setattr(path_utils, "PACKAGE_ROOT", package_root)
    monkeypatch.setattr(path_utils, "get_default_output_dir", lambda: safe_default)

    def _inside(p=None):
        # First call (cwd) => False, second call (target_path) => True.
        candidate = Path.cwd() if p is None else Path(p)
        return candidate.name == "force_inside"

    monkeypatch.setattr(path_utils, "is_inside_package_dir", _inside)

    out = path_utils.get_safe_output_path("./force_inside")
    # relative_to(package_root) fails -> fallback to safe_default / target_path.name
    assert out == safe_default / "force_inside"
    assert out.exists()


def test_absolute_path_branch_then_double_check_relative_to_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    safe_default = tmp_path / "safe"
    safe_default.mkdir()
    package_root = tmp_path / "pkg_root"
    package_root.mkdir()
    user_abs = package_root / "nested" / "plots"

    monkeypatch.setattr(path_utils, "PACKAGE_ROOT", package_root)
    monkeypatch.setattr(path_utils, "get_default_output_dir", lambda: safe_default)

    calls = {"n": 0}

    def _inside(p=None):
        # First check on absolute user path returns False -> line 89 path.
        # Second check on resolved target returns True -> line 108 path.
        calls["n"] += 1
        if calls["n"] == 1:
            return False
        candidate = Path.cwd() if p is None else Path(p)
        return candidate.resolve().is_relative_to(package_root.resolve())

    monkeypatch.setattr(path_utils, "is_inside_package_dir", _inside)

    out = path_utils.get_safe_output_path(str(user_abs))
    assert out == safe_default / "nested" / "plots"
    assert out.exists()
