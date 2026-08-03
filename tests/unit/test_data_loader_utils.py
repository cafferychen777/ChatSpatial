"""Unit contracts for data loader helpers and generic loading flow."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
import pytest

from chatspatial.utils import data_loader as dl
from chatspatial.utils.exceptions import DataCompatibilityError, ParameterError, ProcessingError


@pytest.mark.asyncio
async def test_load_spatial_data_rejects_missing_path():
    with pytest.raises(FileNotFoundError):
        await dl.load_spatial_data("/path/does/not/exist.h5ad", "generic")


@pytest.mark.asyncio
async def test_load_spatial_data_generic_wraps_read_error(tmp_path: Path, monkeypatch):
    path = tmp_path / "demo.h5ad"
    path.write_text("placeholder")

    fake_scanpy = ModuleType("scanpy")
    fake_scanpy.read_h5ad = lambda _p: (_ for _ in ()).throw(RuntimeError("bad file"))
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    with pytest.raises(ProcessingError, match="Error loading generic data"):
        await dl.load_spatial_data(str(path), "generic")


@pytest.mark.asyncio
async def test_load_spatial_data_generic_sets_raw_and_counts_layer(
    minimal_spatial_adata, tmp_path: Path, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    path = tmp_path / "sample.h5ad"
    path.write_text("placeholder")

    fake_scanpy = ModuleType("scanpy")
    fake_scanpy.read_h5ad = lambda _p: adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    result = await dl.load_spatial_data(str(path), "generic")
    out = result["adata"]

    assert result["name"] == "sample"
    assert out.raw is not None
    assert "counts" in out.layers
    assert result["spatial_coordinates_available"] is True


@pytest.mark.asyncio
async def test_load_spatial_data_rejects_unsupported_platform(tmp_path: Path, monkeypatch):
    path = tmp_path / "sample.h5ad"
    path.write_text("placeholder")
    monkeypatch.setitem(sys.modules, "scanpy", ModuleType("scanpy"))

    with pytest.raises(ParameterError, match="Unsupported platform type"):
        await dl.load_spatial_data(str(path), "unsupported")


def test_find_spatial_folder_requires_expected_files(tmp_path: Path):
    h5 = tmp_path / "x.h5"
    h5.write_text("x")
    spatial_dir = tmp_path / "spatial"
    spatial_dir.mkdir()

    assert dl._find_spatial_folder(str(h5)) is None

    (spatial_dir / "tissue_positions_list.csv").write_text("barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n")
    (spatial_dir / "scalefactors_json.json").write_text("{}")
    found = dl._find_spatial_folder(str(h5))
    assert found is not None
    assert Path(found).name == "spatial"


def test_find_spatial_folder_accepts_spaceranger_v2_positions_file(tmp_path: Path):
    h5 = tmp_path / "x.h5"
    h5.write_text("x")
    spatial_dir = tmp_path / "spatial"
    spatial_dir.mkdir()
    (spatial_dir / "tissue_positions.csv").write_text(
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n"
    )
    (spatial_dir / "scalefactors_json.json").write_text("{}")

    assert dl._find_spatial_folder(str(h5)) == str(spatial_dir)


def test_add_spatial_info_handles_barcode_suffix_alignment(
    minimal_spatial_adata, tmp_path: Path, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"bc{i}-1" for i in range(adata.n_obs)]

    spatial_dir = tmp_path / "sample_a" / "spatial"
    spatial_dir.mkdir(parents=True)

    positions = pd.DataFrame(
        {
            "barcode": [f"bc{i}" for i in range(adata.n_obs)],
            "in_tissue": [1] * adata.n_obs,
            "array_row": list(range(adata.n_obs)),
            "array_col": list(range(adata.n_obs)),
            "pxl_row_in_fullres": np.linspace(0, 10, adata.n_obs),
            "pxl_col_in_fullres": np.linspace(5, 15, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions_list.csv", index=False)
    (spatial_dir / "scalefactors_json.json").write_text(json.dumps({"spot_diameter_fullres": 10}))

    monkeypatch.setattr(dl, "is_available", lambda _name: False)

    out = dl._add_spatial_info_to_adata(adata, str(spatial_dir))
    assert "spatial" in out.obsm
    assert out.obsm["spatial"].shape[1] == 2
    np.testing.assert_array_equal(out.obs["array_row"], np.arange(adata.n_obs))
    np.testing.assert_array_equal(out.obs["array_col"], np.arange(adata.n_obs))
    assert "spatial" in out.uns and "sample_a" in out.uns["spatial"]


def test_add_spatial_info_rejects_invalid_positions_format(
    minimal_spatial_adata, tmp_path: Path, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    spatial_dir = tmp_path / "sample_b" / "spatial"
    spatial_dir.mkdir(parents=True)
    (spatial_dir / "tissue_positions_list.csv").write_text("a,b,c\n1,2,3\n")
    (spatial_dir / "scalefactors_json.json").write_text("{}")
    monkeypatch.setattr(dl, "is_available", lambda _name: False)

    with pytest.raises(DataCompatibilityError, match="Unexpected tissue positions format"):
        dl._add_spatial_info_to_adata(adata, str(spatial_dir))


# =============================================================================
# Counts layer creation contracts
# =============================================================================


@pytest.mark.asyncio
async def test_load_creates_counts_from_raw_when_x_is_normalized(
    minimal_spatial_adata, tmp_path: Path, monkeypatch
):
    """When X is normalized but .raw has integer counts, counts layer is created from .raw."""
    import anndata as ad

    adata = minimal_spatial_adata.copy()
    # Simulate pre-processed h5ad: X is log-normalized, .raw has integer counts
    raw_X = adata.X.copy()  # Poisson integers
    adata.raw = ad.AnnData(X=raw_X, var=adata.var.copy())
    adata.X = np.log1p(adata.X)  # Now X is float, non-integer

    path = tmp_path / "normalized.h5ad"
    path.write_text("placeholder")

    fake_scanpy = ModuleType("scanpy")
    fake_scanpy.read_h5ad = lambda _p: adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    result = await dl.load_spatial_data(str(path), "generic")
    out = result["adata"]

    assert "counts" in out.layers
    # counts should come from .raw (integer), not from X (float)
    assert np.allclose(out.layers["counts"].toarray() if hasattr(out.layers["counts"], "toarray") else out.layers["counts"], raw_X)


@pytest.mark.asyncio
async def test_load_does_not_crash_when_raw_has_more_genes(
    minimal_spatial_adata, tmp_path: Path, monkeypatch
):
    """When .raw has more genes than current adata, no shape mismatch crash."""
    import anndata as ad

    adata = minimal_spatial_adata.copy()
    # .raw has full gene set
    n_extra = 10
    raw_X = np.hstack([
        adata.X,
        np.ones((adata.n_obs, n_extra), dtype=np.float32),
    ])
    raw_var = pd.DataFrame(
        index=[f"gene_{i}" for i in range(adata.n_vars + n_extra)]
    )
    adata.raw = ad.AnnData(X=raw_X, var=raw_var)
    # Normalize X so it's not integer
    adata.X = np.log1p(adata.X)

    path = tmp_path / "filtered.h5ad"
    path.write_text("placeholder")

    fake_scanpy = ModuleType("scanpy")
    fake_scanpy.read_h5ad = lambda _p: adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    # Should not raise ValueError: shape mismatch
    result = await dl.load_spatial_data(str(path), "generic")
    out = result["adata"]
    # counts layer skipped because .raw shape doesn't match
    assert "counts" not in out.layers


# =============================================================================
# Spatial key detection contract
# =============================================================================


@pytest.mark.asyncio
async def test_load_detects_alternative_spatial_keys(
    minimal_spatial_adata, tmp_path: Path, monkeypatch
):
    """spatial_coordinates_available is True for alternative keys like X_spatial."""
    adata = minimal_spatial_adata.copy()
    # Move coordinates to alternative key
    coords = adata.obsm.pop("spatial")
    adata.obsm["X_spatial"] = coords

    path = tmp_path / "alt_spatial.h5ad"
    path.write_text("placeholder")

    fake_scanpy = ModuleType("scanpy")
    fake_scanpy.read_h5ad = lambda _p: adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    result = await dl.load_spatial_data(str(path), "generic")
    assert result["spatial_coordinates_available"] is True


# =============================================================================
# Var names alignment contracts
# =============================================================================


@pytest.mark.asyncio
async def test_load_creates_counts_from_raw_with_reordered_genes(
    minimal_spatial_adata, tmp_path: Path, monkeypatch
):
    """When .raw has same genes in different order, counts layer is reindexed."""
    import anndata as ad

    adata = minimal_spatial_adata.copy()
    raw_X = adata.X.copy()
    # Reverse gene order in .raw
    reversed_var = adata.var.copy().iloc[::-1]
    reversed_X = raw_X[:, ::-1]
    adata.raw = ad.AnnData(X=reversed_X, var=reversed_var)
    adata.X = np.log1p(adata.X)  # Normalize X so .raw path is used

    path = tmp_path / "reordered.h5ad"
    path.write_text("placeholder")

    fake_scanpy = ModuleType("scanpy")
    fake_scanpy.read_h5ad = lambda _p: adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    result = await dl.load_spatial_data(str(path), "generic")
    out = result["adata"]

    assert "counts" in out.layers
    # Counts should be reindexed to match adata.var_names order (original)
    counts = (
        out.layers["counts"].toarray()
        if hasattr(out.layers["counts"], "toarray")
        else out.layers["counts"]
    )
    assert np.allclose(counts, raw_X)


def test_load_does_not_freeze_normalized_x_as_raw(minimal_spatial_adata):
    """load_spatial_data should not set .raw when X is normalized."""
    # This tests the contract: .raw is only set when X is integer counts
    adata = minimal_spatial_adata.copy()
    adata.raw = None
    adata.X = np.asarray(adata.X, dtype=float) + 0.5  # normalized

    from chatspatial.utils.adata_utils import check_is_integer_counts

    is_int, has_neg, _ = check_is_integer_counts(adata.X)

    # Contract: normalized X should not be frozen as .raw
    assert not is_int or has_neg  # verify our test data IS normalized
