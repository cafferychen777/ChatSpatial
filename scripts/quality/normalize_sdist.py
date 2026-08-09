#!/usr/bin/env python3
"""Rewrite an sdist with deterministic tar and gzip metadata."""

from __future__ import annotations

import argparse
import gzip
import os
import tarfile
import tempfile
from pathlib import Path


def _copy_member(
    source: tarfile.TarFile,
    target: tarfile.TarFile,
    member: tarfile.TarInfo,
    epoch: int,
) -> None:
    if not (member.isfile() or member.isdir()):
        raise ValueError(f"unsupported sdist member type: {member.name!r}")

    normalized = tarfile.TarInfo(member.name)
    normalized.type = tarfile.DIRTYPE if member.isdir() else tarfile.REGTYPE
    normalized.mode = member.mode
    normalized.uid = 0
    normalized.gid = 0
    normalized.uname = ""
    normalized.gname = ""
    normalized.mtime = epoch
    normalized.size = member.size if member.isfile() else 0

    if member.isfile():
        content = source.extractfile(member)
        if content is None:
            raise ValueError(f"cannot read sdist member: {member.name!r}")
        try:
            target.addfile(normalized, content)
        finally:
            content.close()
    else:
        target.addfile(normalized)


def normalize_sdist(sdist_path: Path, epoch: int) -> None:
    """Atomically normalize member order, ownership, timestamps, and gzip header."""
    if epoch < 0:
        raise ValueError("SOURCE_DATE_EPOCH must be non-negative")

    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=sdist_path.parent,
            prefix=f".{sdist_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            with gzip.GzipFile(
                filename="",
                mode="wb",
                compresslevel=9,
                fileobj=temporary,
                mtime=epoch,
            ) as compressed:
                with tarfile.open(
                    fileobj=compressed,
                    mode="w",
                    format=tarfile.PAX_FORMAT,
                ) as target:
                    with tarfile.open(sdist_path, mode="r:gz") as source:
                        for member in sorted(
                            source.getmembers(), key=lambda item: item.name
                        ):
                            _copy_member(source, target, member, epoch)
        os.replace(temporary_path, sdist_path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdist", type=Path)
    parser.add_argument("epoch", type=int)
    args = parser.parse_args()
    normalize_sdist(args.sdist.resolve(), args.epoch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
