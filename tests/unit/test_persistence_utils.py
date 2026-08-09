"""Unit tests for persistence utilities."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from chatspatial.utils import persistence
from chatspatial.utils.exceptions import ParameterError


def _make_adata(n_obs: int = 8, n_vars: int = 5):
    X = np.arange(n_obs * n_vars, dtype=np.float32).reshape(n_obs, n_vars)
    adata = ad.AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    return adata


def test_get_active_path_uses_home(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    p = persistence.get_active_path("demo")
    assert p == tmp_path / ".chatspatial" / "active" / "demo.h5ad"


def test_get_active_path_rejects_path_traversal(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    with pytest.raises(ParameterError, match="filesystem-safe"):
        persistence.get_active_path("../outside")

    assert not (tmp_path / ".chatspatial" / "outside.h5ad").exists()


def test_export_and_load_roundtrip_with_custom_path(tmp_path: Path):
    adata = _make_adata()
    out = tmp_path / "exported" / "roundtrip.h5ad"

    exported = persistence.export_adata("d1", adata, out)
    loaded = persistence.load_adata_from_active("d1", out)

    assert exported == out
    assert exported.exists()
    assert loaded.shape == adata.shape


def test_export_sanitizes_uns_dataframe_object_columns(tmp_path: Path):
    adata = _make_adata()
    results = pd.DataFrame(
        {
            "gene_a": ["CXCL12", None, np.nan, 123],
            "gene_b": ["CXCR4", "CCR7", None, "CD74"],
            "T|T": [0.1, 0.2, 0.3, 0.4],
        },
        index=[10, 11, 12, 13],
    )
    adata.uns["ccc"] = {
        "method": "cellphonedb",
        "results": results,
        "method_data": {"significant_means": results.copy()},
    }
    adata.uns["ccc_cellphonedb"] = {
        "method": "cellphonedb",
        "results": results.copy(),
    }
    adata.uns["gsea_results"] = pd.DataFrame(
        {
            "Term": ["Pathway A", "Pathway B"],
            "ES": [None, {"unexpected": "object"}],
            "NES": [1.2, -0.4],
            "FDR q-val": pd.Series(["0.05", "0.20"], dtype=object),
        }
    )
    out = tmp_path / "exported" / "uns_roundtrip.h5ad"

    exported = persistence.export_adata("d1", adata, out)
    loaded = persistence.load_adata_from_active("d1", exported)

    loaded_results = loaded.uns["ccc"]["results"]
    assert "ccc_cellphonedb" in loaded.uns
    assert list(loaded_results.columns) == ["gene_a", "gene_b", "T|T"]
    assert list(loaded_results.index) == ["10", "11", "12", "13"]
    assert loaded_results["gene_a"].tolist() == ["CXCL12", "", "", "123"]
    assert np.issubdtype(loaded_results["T|T"].dtype, np.number)
    loaded_gsea = loaded.uns["gsea_results"]
    assert loaded_gsea["ES"].tolist() == ["", "{'unexpected': 'object'}"]
    assert np.issubdtype(loaded_gsea["NES"].dtype, np.number)
    assert loaded_gsea["FDR q-val"].tolist() == [0.05, 0.2]
    assert np.issubdtype(loaded_gsea["FDR q-val"].dtype, np.number)
    assert adata.uns["ccc"]["results"].index.tolist() == [10, 11, 12, 13]


def test_load_adata_missing_file_raises_file_not_found(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="Data file not found"):
        persistence.load_adata_from_active("missing", tmp_path / "missing.h5ad")


def test_export_uses_active_path_when_custom_path_not_provided(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    adata = _make_adata()

    exported = persistence.export_adata("active_ds", adata)

    assert exported == tmp_path / ".chatspatial" / "active" / "active_ds.h5ad"
    assert exported.exists()


def test_export_wraps_underlying_write_errors(tmp_path: Path):
    class _BrokenAnnData:
        def write_h5ad(self, *_args, **_kwargs):
            raise RuntimeError("disk full")

    with pytest.raises(IOError, match="Failed to export data"):
        persistence.export_adata(
            "broken", _BrokenAnnData(), tmp_path / "x" / "broken.h5ad"
        )


def test_export_failure_preserves_previous_complete_file(tmp_path: Path) -> None:
    class _PartiallyWritingAnnData:
        uns: dict[str, object] = {}

        def copy(self):
            return self

        def write_h5ad(self, path: Path, **_kwargs) -> None:
            Path(path).write_bytes(b"partial")
            raise RuntimeError("disk full")

    destination = tmp_path / "active" / "dataset.h5ad"
    destination.parent.mkdir()
    destination.write_bytes(b"previous complete dataset")

    with pytest.raises(IOError, match="Failed to export data"):
        persistence.export_adata("dataset", _PartiallyWritingAnnData(), destination)

    assert destination.read_bytes() == b"previous complete dataset"
    assert list(destination.parent.iterdir()) == [destination]


def test_load_adata_uses_active_path_when_custom_path_not_provided(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    adata = _make_adata()
    persistence.export_adata("active_load", adata)

    loaded = persistence.load_adata_from_active("active_load")

    assert loaded.shape == adata.shape


def test_load_adata_wraps_underlying_read_errors(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    corrupt_path = tmp_path / "broken.h5ad"
    corrupt_path.write_bytes(b"corrupt")

    import anndata

    monkeypatch.setattr(
        anndata,
        "read_h5ad",
        lambda _path: (_ for _ in ()).throw(RuntimeError("bad format")),
    )

    with pytest.raises(IOError, match="Failed to load data"):
        persistence.load_adata_from_active("unused", corrupt_path)
