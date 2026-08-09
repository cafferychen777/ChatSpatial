#!/usr/bin/env python3
"""Fetch datasets from a registry into an external data directory."""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import zipfile
from pathlib import Path
from urllib.request import urlopen


def _load_registry(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, out_file: Path) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url) as resp, out_file.open("wb") as f:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)


def _extract_if_needed(file_path: Path, target_dir: Path) -> None:
    suffixes = file_path.suffixes
    if suffixes[-1:] == [".zip"]:
        with zipfile.ZipFile(file_path, "r") as zf:
            zf.extractall(target_dir)
        return

    tar_like = suffixes[-2:] in ([".tar", ".gz"], [".tar", ".bz2"]) or suffixes[
        -1:
    ] == [".tgz"]
    if tar_like:
        with tarfile.open(file_path, "r:*") as tf:
            tf.extractall(target_dir)


def _iter_enabled_datasets(registry: dict):
    for ds in registry.get("datasets", []):
        if ds.get("enabled", True):
            yield ds


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch external datasets from registry"
    )
    parser.add_argument("--registry", default="scripts/data/dataset_registry.json")
    parser.add_argument("--list", action="store_true", help="List enabled datasets")
    parser.add_argument(
        "--dataset", action="append", help="Dataset id to fetch (repeatable)"
    )
    parser.add_argument(
        "--dest", default=str(Path.home() / ".cache" / "chatspatial" / "datasets")
    )
    parser.add_argument(
        "--force", action="store_true", help="Re-download even if file exists"
    )
    args = parser.parse_args()

    registry = _load_registry(Path(args.registry))
    datasets = list(_iter_enabled_datasets(registry))

    if args.list:
        for ds in datasets:
            print(f"{ds['id']}: {ds.get('description', '')}")
        return

    selected = set(args.dataset or [])
    if selected:
        datasets = [ds for ds in datasets if ds["id"] in selected]

    if not datasets:
        raise SystemExit(
            "No enabled datasets selected. Check --dataset and registry enabled flags."
        )

    dest = Path(args.dest)
    dest.mkdir(parents=True, exist_ok=True)

    for ds in datasets:
        ds_dir = dest / ds["id"]
        ds_dir.mkdir(parents=True, exist_ok=True)

        filename = ds.get("filename") or Path(ds["url"]).name
        out_file = ds_dir / filename

        if out_file.exists() and not args.force:
            print(f"[skip] {ds['id']} already exists: {out_file}")
        else:
            print(f"[download] {ds['id']} -> {out_file}")
            _download(ds["url"], out_file)

        expected_sha = ds.get("sha256", "").strip()
        if expected_sha:
            got_sha = _sha256(out_file)
            if got_sha.lower() != expected_sha.lower():
                raise ValueError(
                    f"Checksum mismatch for {ds['id']}: {got_sha} != {expected_sha}"
                )

        if ds.get("extract", False):
            print(f"[extract] {ds['id']}")
            _extract_if_needed(out_file, ds_dir)

    print(f"Done. Datasets available under: {dest}")


if __name__ == "__main__":
    main()
