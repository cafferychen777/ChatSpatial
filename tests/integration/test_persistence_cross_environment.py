"""Bidirectional interchange between independently installed AnnData versions."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.integration
def test_h5ad_roundtrip_across_environments(tmp_path):
    other_python = os.environ.get("CHATSPATIAL_COMPAT_PYTHON")
    if not other_python:
        pytest.skip("Set CHATSPATIAL_COMPAT_PYTHON to the second runtime")
    worker = Path(__file__).with_name("persistence_roundtrip_worker.py")
    root = Path(__file__).resolve().parents[2]
    env = {**os.environ, "PYTHONPATH": str(root)}
    for label, first, second in [
        ("forward", sys.executable, other_python),
        ("reverse", other_python, sys.executable),
    ]:
        paths = [tmp_path / f"{label}-{i}.h5ad" for i in range(3)]
        for executable, source, destination in [
            (first, "create", paths[0]),
            (second, paths[0], paths[1]),
            (first, paths[1], paths[2]),
        ]:
            subprocess.run(
                [executable, str(worker), str(source), str(destination)],
                env=env,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
        snapshots = [json.loads(Path(str(p) + ".json").read_text()) for p in paths]
        assert snapshots[0] == snapshots[1] == snapshots[2]
