#!/usr/bin/env python3
"""Shared path discovery for ChatSpatial reproducibility scripts."""

from __future__ import annotations

import os
from pathlib import Path

REPRODUCIBILITY_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = REPRODUCIBILITY_ROOT.parent


def _env_path(name: str) -> list[Path]:
    value = os.environ.get(name, "").strip()
    if not value:
        return []
    return [Path(value).expanduser().resolve()]


def _unique(paths: list[Path]) -> list[Path]:
    seen = set()
    out = []
    for path in paths:
        key = str(path)
        if key not in seen:
            seen.add(key)
            out.append(path)
    return out


def find_chatspatial_code_dir(required: bool = False) -> Path | None:
    """Return the ChatSpatial package source directory, if discoverable.

    Set CHATSPATIAL_CODE_DIR when the package repository is not in a common
    workspace location.
    """
    candidates = _unique(
        _env_path("CHATSPATIAL_CODE_DIR")
        + [
            PROJECT_ROOT,
            REPRODUCIBILITY_ROOT / "code",
            PROJECT_ROOT.parent / "ChatSpatial" / "code",
            PROJECT_ROOT.parent / "ChatSpatial",
        ]
    )
    for candidate in candidates:
        if (candidate / "chatspatial").is_dir():
            return candidate
        nested_code = candidate / "code"
        if (nested_code / "chatspatial").is_dir():
            return nested_code
    if required:
        raise FileNotFoundError(
            "Could not locate the ChatSpatial source tree. Set "
            "CHATSPATIAL_CODE_DIR=/path/to/ChatSpatial/code."
        )
    return None


def find_benchmarks_dir(required: bool = False) -> Path | None:
    """Return the benchmarks directory containing STAgent/SpatialAgent."""
    candidates = _unique(
        _env_path("CHATSPATIAL_BENCHMARKS_DIR")
        + [
            PROJECT_ROOT / "benchmarks",
            PROJECT_ROOT.parent / "benchmarks",
            REPRODUCIBILITY_ROOT / "benchmarks",
        ]
    )
    for candidate in candidates:
        if candidate.is_dir() and (
            (candidate / "STAgent").exists() or (candidate / "SpatialAgent").exists()
        ):
            return candidate
    if required:
        raise FileNotFoundError(
            "Could not locate benchmarks/STAgent and benchmarks/SpatialAgent. "
            "Set CHATSPATIAL_BENCHMARKS_DIR=/path/to/benchmarks."
        )
    return None


def find_competitor_dir(
    dirname: str,
    env_var: str,
    required: bool = False,
) -> Path | None:
    """Return an installed competitor framework directory."""
    candidates = _env_path(env_var)
    benchmarks_dir = find_benchmarks_dir(required=False)
    if benchmarks_dir is not None:
        candidates.append(benchmarks_dir / dirname)
    candidates.extend(
        [
            PROJECT_ROOT.parent / ".competitor_analysis" / dirname,
            PROJECT_ROOT.parent / "benchmarks" / dirname,
            PROJECT_ROOT / "benchmarks" / dirname,
        ]
    )
    for candidate in _unique(candidates):
        if candidate.is_dir():
            return candidate
    if required:
        raise FileNotFoundError(
            f"Could not locate {dirname}. Set {env_var}=/path/to/{dirname}."
        )
    return None


def find_env_file() -> Path | None:
    """Return the first workspace .env file found."""
    candidates = _unique(
        _env_path("CHATSPATIAL_ENV_FILE")
        + [
            PROJECT_ROOT / ".env",
            PROJECT_ROOT.parent / ".env",
            REPRODUCIBILITY_ROOT / ".env",
        ]
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def load_env_file(path: Path | None = None) -> None:
    """Load simple KEY=VALUE pairs without requiring python-dotenv."""
    env_path = path or find_env_file()
    if env_path is None:
        return
    with open(env_path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            os.environ.setdefault(key.strip(), value.strip().strip("\"'"))
