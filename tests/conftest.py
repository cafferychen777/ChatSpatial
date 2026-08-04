"""Shared pytest configuration and fixtures for ChatSpatial tests."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pytest

try:
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover - matplotlib is an optional dependency
    plt = None


try:
    import anndata as ad
except ImportError:  # pragma: no cover - handled by skip
    ad = None


def _to_jsonable(value: Any) -> Any:
    """Convert nested objects to JSON-serializable structures."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _to_jsonable(v) for k, v in value.items()}
    if hasattr(value, "model_dump"):  # pydantic models
        return _to_jsonable(value.model_dump())
    return repr(value)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Expose per-phase test report on test item for fixture finalizers."""
    outcome = yield
    report = outcome.get_result()
    setattr(item, f"rep_{report.when}", report)


@dataclass
class E2ETrace:
    """Per-test trace collector for failure diagnostics."""

    test_nodeid: str
    events: list[dict[str, Any]] = field(default_factory=list)

    def record(
        self,
        *,
        step: str,
        data_id: str | None = None,
        params: Any = None,
        notes: str | None = None,
    ) -> None:
        self.events.append(
            {
                "step": step,
                "data_id": data_id,
                "params": _to_jsonable(params),
                "notes": notes,
            }
        )


@pytest.fixture
def e2e_trace(request: pytest.FixtureRequest):
    """Collect execution trace and dump JSON artifact on e2e test failure."""
    trace = E2ETrace(test_nodeid=request.node.nodeid)
    yield trace

    if request.node.get_closest_marker("e2e") is None:
        return

    rep_call = getattr(request.node, "rep_call", None)
    if rep_call is None or not rep_call.failed:
        return

    artifact_dir = Path("tests/artifacts/e2e_failures")
    artifact_dir.mkdir(parents=True, exist_ok=True)

    safe_name = request.node.nodeid.replace("/", "__").replace("::", "__")
    output_file = artifact_dir / f"{safe_name}.json"
    payload = {
        "test": request.node.nodeid,
        "outcome": "failed",
        "events": trace.events,
    }
    output_file.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )


@pytest.fixture
def rng() -> np.random.Generator:
    """Seeded RNG for reproducible test data."""
    return np.random.default_rng(42)


@pytest.fixture(autouse=True)
def close_matplotlib_figures():
    """Prevent cross-test figure leakage in visualization-heavy test runs."""
    yield
    if plt is not None:
        plt.close("all")


@pytest.fixture
def minimal_spatial_adata(rng: np.random.Generator):
    """Create a small AnnData object with spatial coordinates and group labels."""
    if ad is None:
        pytest.skip("anndata is required for dataset fixtures")

    n_cells = 60
    n_genes = 24
    X = rng.poisson(4, size=(n_cells, n_genes)).astype(np.float32)

    adata = ad.AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]

    # Two balanced groups for DE tests
    adata.obs["group"] = ["A"] * (n_cells // 2) + ["B"] * (n_cells - n_cells // 2)

    # Minimal spatial metadata for visualization tests
    adata.obsm["spatial"] = rng.uniform(0, 100, size=(n_cells, 2))
    adata.uns["spatial"] = {
        "library_id": {
            "images": {},
            "scalefactors": {"spot_diameter_fullres": 35.0},
        }
    }

    return adata


@pytest.fixture
def write_h5ad_dataset(tmp_path: Path):
    """Return helper that writes AnnData to temp h5ad and returns its path."""

    def _writer(adata_obj, filename: str = "dataset.h5ad") -> Path:
        out = tmp_path / filename
        adata_obj.write_h5ad(out)
        return out

    return _writer


@pytest.fixture
def spatial_dataset_path(minimal_spatial_adata, write_h5ad_dataset) -> Path:
    """Write and return a reusable spatial test dataset path."""
    return write_h5ad_dataset(minimal_spatial_adata, "spatial_test.h5ad")


@pytest.fixture
def reset_data_manager():
    """Reset global data manager state before and after each test."""
    from chatspatial.server import data_manager

    data_manager.reset()
    yield data_manager
    data_manager.reset()


@pytest.fixture
def deepspotm_stack(monkeypatch):
    """Install stand-in deepspotm and torch modules for one test.

    Keeps DeepSpot-M coverage free of the real package, its weights, and any
    network access, which is what lets CI run these tests.
    """
    import sys

    from chatspatial.utils import dependency_manager as dm
    from tests.fixtures.helpers import (
        fake_deepspotm_module,
        fake_torch_module,
    )

    recorder: dict[str, Any] = {}
    monkeypatch.setitem(sys.modules, "deepspotm", fake_deepspotm_module(recorder))
    monkeypatch.setitem(sys.modules, "torch", fake_torch_module())
    monkeypatch.setattr(
        "chatspatial.tools.histology._resolve_checkpoint_revision",
        lambda params: params.model_revision or "resolved-main-sha",
    )
    dm._try_import.cache_clear()
    dm._check_spec.cache_clear()
    yield recorder
    dm._try_import.cache_clear()
    dm._check_spec.cache_clear()
