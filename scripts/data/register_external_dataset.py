#!/usr/bin/env python3
"""Register local external dataset paths for reproducible local development."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Register a local dataset path")
    parser.add_argument("dataset_id", help="Logical dataset id")
    parser.add_argument("path", help="Absolute path to local dataset file or directory")
    parser.add_argument("--manifest", default="data/datasets.local.json")
    args = parser.parse_args()

    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)

    payload = {"datasets": {}}
    if manifest.exists():
        with manifest.open("r", encoding="utf-8") as f:
            payload = json.load(f)

    payload.setdefault("datasets", {})[args.dataset_id] = str(
        Path(args.path).expanduser().resolve()
    )

    with manifest.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)

    print(f"Registered {args.dataset_id} -> {payload['datasets'][args.dataset_id]}")


if __name__ == "__main__":
    main()
