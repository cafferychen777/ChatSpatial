#!/usr/bin/env python3
"""Audit that release archives contain only intended package files."""

from __future__ import annotations

import argparse
import stat
import sys
import tarfile
import tomllib
import zipfile
from pathlib import Path, PurePosixPath

MAX_ARCHIVE_BYTES = 50 * 1024 * 1024
MAX_MEMBER_BYTES = 20 * 1024 * 1024
MAX_TOTAL_BYTES = 100 * 1024 * 1024
MAX_MEMBERS = 5_000

FORBIDDEN_PARTS = frozenset(
    {
        ".git",
        ".github",
        ".mypy_cache",
        ".pytest_cache",
        ".ruff_cache",
        ".tox",
        ".venv",
        "__pycache__",
        "build",
        "data",
        "dist",
        "docs",
        "env",
        "scripts",
        "tests",
        "third_party",
        "venv",
    }
)
FORBIDDEN_SUFFIXES = (
    ".csv",
    ".csv.gz",
    ".h5",
    ".h5ad",
    ".hdf5",
    ".log",
    ".pyc",
    ".pyo",
    ".temp",
    ".tmp",
)


class DistributionError(ValueError):
    """Raised when a distribution violates the release contract."""


def _project_version(project_root: Path) -> str:
    with (project_root / "pyproject.toml").open("rb") as handle:
        project = tomllib.load(handle)["project"]
    if project["name"] != "chatspatial":
        raise DistributionError("unexpected project name in pyproject.toml")
    return str(project["version"])


def _validated_parts(name: str) -> tuple[str, ...]:
    if "\\" in name:
        raise DistributionError(f"archive member uses a backslash: {name!r}")

    path = PurePosixPath(name)
    parts = path.parts
    if path.is_absolute() or not parts or ".." in parts:
        raise DistributionError(f"unsafe archive member path: {name!r}")

    normalized_parts = tuple(part.casefold() for part in parts)
    forbidden = FORBIDDEN_PARTS.intersection(normalized_parts)
    if forbidden:
        found = ", ".join(sorted(forbidden))
        raise DistributionError(f"forbidden path component ({found}): {name!r}")

    lowered_name = name.casefold()
    if lowered_name.endswith(FORBIDDEN_SUFFIXES):
        raise DistributionError(f"forbidden file type: {name!r}")
    return parts


def _validate_member_budget(sizes: list[int], archive_name: str) -> None:
    if len(sizes) > MAX_MEMBERS:
        raise DistributionError(
            f"{archive_name} contains {len(sizes)} members; limit is {MAX_MEMBERS}"
        )
    if any(size < 0 or size > MAX_MEMBER_BYTES for size in sizes):
        raise DistributionError(
            f"{archive_name} contains a member larger than {MAX_MEMBER_BYTES} bytes"
        )
    if sum(sizes) > MAX_TOTAL_BYTES:
        raise DistributionError(
            f"{archive_name} expands beyond {MAX_TOTAL_BYTES} bytes"
        )


def _validate_wheel(wheel_path: Path, version: str) -> None:
    dist_info = f"chatspatial-{version}.dist-info"
    required = {
        "chatspatial/__init__.py",
        "chatspatial/py.typed",
        f"{dist_info}/METADATA",
        f"{dist_info}/WHEEL",
        f"{dist_info}/RECORD",
        f"{dist_info}/entry_points.txt",
    }

    with zipfile.ZipFile(wheel_path) as archive:
        members = archive.infolist()
        names = [member.filename for member in members]
        if len(names) != len(set(names)):
            raise DistributionError(f"{wheel_path.name} contains duplicate members")

        for member in members:
            parts = _validated_parts(member.filename)
            if parts[0] not in {"chatspatial", dist_info}:
                raise DistributionError(
                    f"unexpected wheel root {parts[0]!r}: {member.filename!r}"
                )
            mode = (member.external_attr >> 16) & 0o170000
            if mode == stat.S_IFLNK or member.flag_bits & 0x1:
                raise DistributionError(
                    f"link or encrypted wheel member: {member.filename!r}"
                )

        _validate_member_budget(
            [member.file_size for member in members], wheel_path.name
        )
        missing = required.difference(names)
        if missing:
            raise DistributionError(
                f"{wheel_path.name} is missing: {', '.join(sorted(missing))}"
            )

        metadata = archive.read(f"{dist_info}/METADATA").decode("utf-8")
        if "Name: chatspatial\n" not in metadata:
            raise DistributionError("wheel metadata has the wrong project name")
        if f"Version: {version}\n" not in metadata:
            raise DistributionError("wheel metadata has the wrong version")


def _validate_sdist(sdist_path: Path, version: str) -> None:
    root = f"chatspatial-{version}"
    required = {
        f"{root}/LICENSE",
        f"{root}/PKG-INFO",
        f"{root}/README.md",
        f"{root}/chatspatial/__init__.py",
        f"{root}/chatspatial/py.typed",
        f"{root}/pyproject.toml",
    }

    with tarfile.open(sdist_path, mode="r:gz") as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        if len(names) != len(set(names)):
            raise DistributionError(f"{sdist_path.name} contains duplicate members")

        for member in members:
            parts = _validated_parts(member.name)
            if parts[0] != root:
                raise DistributionError(
                    f"unexpected sdist root {parts[0]!r}: {member.name!r}"
                )
            if not (member.isfile() or member.isdir()):
                raise DistributionError(f"non-regular sdist member: {member.name!r}")

        _validate_member_budget([member.size for member in members], sdist_path.name)
        missing = required.difference(names)
        if missing:
            raise DistributionError(
                f"{sdist_path.name} is missing: {', '.join(sorted(missing))}"
            )

        package_info = archive.extractfile(f"{root}/PKG-INFO")
        if package_info is None:
            raise DistributionError("cannot read sdist package metadata")
        with package_info:
            metadata = package_info.read().decode("utf-8")
        if "Name: chatspatial\n" not in metadata:
            raise DistributionError("sdist metadata has the wrong project name")
        if f"Version: {version}\n" not in metadata:
            raise DistributionError("sdist metadata has the wrong version")


def validate_distribution(project_root: Path, distribution_dir: Path) -> None:
    """Validate filenames, metadata, paths, types, and size limits."""
    version = _project_version(project_root)
    expected_names = {
        f"chatspatial-{version}-py3-none-any.whl",
        f"chatspatial-{version}.tar.gz",
    }

    if not distribution_dir.is_dir():
        raise DistributionError(
            f"distribution directory does not exist: {distribution_dir}"
        )
    artifacts = {path.name: path for path in distribution_dir.iterdir()}
    if set(artifacts) != expected_names:
        expected = ", ".join(sorted(expected_names))
        actual = ", ".join(sorted(artifacts)) or "<empty>"
        raise DistributionError(f"expected exactly [{expected}], found [{actual}]")

    for artifact in artifacts.values():
        if not artifact.is_file() or artifact.is_symlink():
            raise DistributionError(f"artifact is not a regular file: {artifact}")
        if artifact.stat().st_size <= 0 or artifact.stat().st_size > MAX_ARCHIVE_BYTES:
            raise DistributionError(f"artifact has an invalid size: {artifact}")

    _validate_wheel(artifacts[f"chatspatial-{version}-py3-none-any.whl"], version)
    _validate_sdist(artifacts[f"chatspatial-{version}.tar.gz"], version)
    print(f"[release] Distribution audit passed for chatspatial {version}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "distribution_dir",
        type=Path,
        help="Directory containing one wheel and one source distribution",
    )
    args = parser.parse_args()
    project_root = Path(__file__).resolve().parents[2]

    try:
        validate_distribution(project_root, args.distribution_dir.resolve())
    except (DistributionError, OSError, tarfile.TarError, zipfile.BadZipFile) as exc:
        print(f"[release] Distribution audit failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
