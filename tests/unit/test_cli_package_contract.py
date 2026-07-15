"""Unit contract for chatspatial.cli package metadata."""

from __future__ import annotations

import importlib
from pathlib import Path

import chatspatial


def test_cli_package_import_contract():
    mod = importlib.import_module("chatspatial.cli")
    assert mod.__all__ == []


def test_package_exposes_pep561_marker():
    package_dir = Path(chatspatial.__file__).parent
    assert (package_dir / "py.typed").is_file()
