"""Unit tests for visium/xenium loading branches and spatial helper paths."""

from __future__ import annotations

import gzip
import json
import os
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
import pytest

from chatspatial.utils import data_loader as dl
from chatspatial.utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)


class _FakeScanpy(ModuleType):
    def __init__(self) -> None:
        super().__init__("scanpy")
        self.calls: dict[str, list] = {
            "read_visium": [],
            "read_10x_mtx": [],
            "read_10x_h5": [],
            "read_h5ad": [],
        }
        self._read_visium_ret = None
        self._read_10x_mtx_ret = None
        self._read_10x_h5_ret = None
        self._read_h5ad_ret = None
        self._read_10x_h5_exc: Exception | None = None

    def read_visium(self, path):
        self.calls["read_visium"].append(path)
        return self._read_visium_ret

    def read_10x_mtx(self, path, var_names="gene_symbols", cache=False):
        self.calls["read_10x_mtx"].append((path, var_names, cache))
        return self._read_10x_mtx_ret

    def read_10x_h5(self, path):
        self.calls["read_10x_h5"].append(path)
        if self._read_10x_h5_exc is not None:
            raise self._read_10x_h5_exc
        return self._read_10x_h5_ret

    def read_h5ad(self, path):
        self.calls["read_h5ad"].append(path)
        return self._read_h5ad_ret


@pytest.mark.asyncio
async def test_load_visium_directory_h5_uses_read_visium(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    visium_dir = tmp_path / "sample_visium"
    visium_dir.mkdir()
    (visium_dir / "filtered_feature_bc_matrix.h5").write_text("h5")

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_visium_ret = minimal_spatial_adata.copy()
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    out = await dl.load_spatial_data(str(visium_dir), "visium")

    assert fake_scanpy.calls["read_visium"] == [str(visium_dir)]
    assert out["type"] == "visium"
    assert out["n_cells"] == minimal_spatial_adata.n_obs


@pytest.mark.asyncio
async def test_load_visium_mtx_directory_adds_spatial_coordinates(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    visium_dir = tmp_path / "sample_mtx"
    mtx_dir = visium_dir / "filtered_feature_bc_matrix"
    spatial_dir = visium_dir / "spatial"
    mtx_dir.mkdir(parents=True)
    spatial_dir.mkdir(parents=True)

    (mtx_dir / "matrix.mtx").write_text("mtx")

    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"bc{i}" for i in range(adata.n_obs)]
    adata.obsm.clear()
    adata.uns.clear()

    positions = pd.DataFrame(
        {
            0: adata.obs_names,
            1: np.ones(adata.n_obs, dtype=int),
            2: np.arange(adata.n_obs),
            3: np.arange(adata.n_obs),
            4: np.linspace(0, 100, adata.n_obs),
            5: np.linspace(20, 120, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions_list.csv", index=False, header=False)
    (spatial_dir / "scalefactors_json.json").write_text(json.dumps({"spot_diameter_fullres": 10}))

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_mtx_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    out = await dl.load_spatial_data(str(visium_dir), "visium")
    out_adata = out["adata"]

    assert fake_scanpy.calls["read_10x_mtx"]
    assert "spatial" in out_adata.obsm
    assert out_adata.obsm["spatial"].shape == (adata.n_obs, 2)


@pytest.mark.asyncio
async def test_load_visium_mtx_directory_rejects_missing_matrix_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    visium_dir = tmp_path / "sample_mtx_missing_matrix"
    visium_dir.mkdir(parents=True)

    monkeypatch.setitem(sys.modules, "scanpy", _FakeScanpy())

    with pytest.raises(
        DataCompatibilityError,
        match="does not have the expected 10x Visium structure",
    ):
        await dl.load_spatial_data(str(visium_dir), "visium")


@pytest.mark.asyncio
async def test_load_visium_filtered_matrix_dir_requires_matrix_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    visium_dir = tmp_path / "sample_missing_mtx_file"
    mtx_dir = visium_dir / "filtered_feature_bc_matrix"
    mtx_dir.mkdir(parents=True)

    monkeypatch.setitem(sys.modules, "scanpy", _FakeScanpy())

    with pytest.raises(
        DataCompatibilityError,
        match="missing matrix.mtx or matrix.mtx.gz",
    ):
        await dl.load_spatial_data(str(visium_dir), "visium")


@pytest.mark.asyncio
async def test_load_visium_mtx_directory_parses_header_positions_file(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    visium_dir = tmp_path / "sample_mtx_header"
    mtx_dir = visium_dir / "filtered_feature_bc_matrix"
    spatial_dir = visium_dir / "spatial"
    mtx_dir.mkdir(parents=True)
    spatial_dir.mkdir(parents=True)

    (mtx_dir / "matrix.mtx").write_text("mtx")

    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"bc{i}" for i in range(adata.n_obs)]
    adata.obsm.clear()
    adata.uns.clear()

    positions = pd.DataFrame(
        {
            "barcode": adata.obs_names,
            "in_tissue": np.ones(adata.n_obs, dtype=int),
            "array_row": np.arange(adata.n_obs),
            "array_col": np.arange(adata.n_obs),
            "pxl_row_in_fullres": np.linspace(0, 100, adata.n_obs),
            "pxl_col_in_fullres": np.linspace(10, 110, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions_list.csv", index=False)
    (spatial_dir / "scalefactors_json.json").write_text("{}")

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_mtx_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    out = await dl.load_spatial_data(str(visium_dir), "visium")
    assert out["adata"].obsm["spatial"].shape == (adata.n_obs, 2)


@pytest.mark.asyncio
async def test_load_visium_mtx_directory_supports_spaceranger_v2_positions(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    visium_dir = tmp_path / "sample_mtx_v2"
    mtx_dir = visium_dir / "filtered_feature_bc_matrix"
    spatial_dir = visium_dir / "spatial"
    mtx_dir.mkdir(parents=True)
    spatial_dir.mkdir(parents=True)
    (mtx_dir / "matrix.mtx").write_text("mtx")

    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"bc{i}" for i in range(adata.n_obs)]
    adata.obsm.clear()
    adata.uns.clear()
    positions = pd.DataFrame(
        {
            "barcode": adata.obs_names,
            "in_tissue": np.ones(adata.n_obs, dtype=int),
            "array_row": np.arange(adata.n_obs),
            "array_col": np.arange(adata.n_obs),
            "pxl_row_in_fullres": np.linspace(0, 100, adata.n_obs),
            "pxl_col_in_fullres": np.linspace(10, 110, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions.csv", index=False)
    (spatial_dir / "scalefactors_json.json").write_text(
        json.dumps({"spot_diameter_fullres": 10})
    )

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_mtx_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setattr(dl, "is_available", lambda _name: False)

    out = await dl.load_spatial_data(str(visium_dir), "visium")
    out_adata = out["adata"]

    assert out_adata.obsm["spatial"].shape == (adata.n_obs, 2)
    assert "sample_mtx_v2" in out_adata.uns["spatial"]
    assert out_adata.uns["spatial"]["sample_mtx_v2"]["scalefactors"] == {
        "spot_diameter_fullres": 10
    }


@pytest.mark.asyncio
async def test_load_visium_mtx_directory_errors_when_spatial_alignment_fails(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """Visium loader must fail-fast when spatial info cannot be loaded.

    Silently degrading to expression-only data would violate the user's
    explicit platform='visium' intent and cause confusing downstream errors.
    """
    visium_dir = tmp_path / "sample_mtx_bad_spatial"
    mtx_dir = visium_dir / "filtered_feature_bc_matrix"
    spatial_dir = visium_dir / "spatial"
    mtx_dir.mkdir(parents=True)
    spatial_dir.mkdir(parents=True)

    (mtx_dir / "matrix.mtx").write_text("mtx")

    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"obs{i}" for i in range(adata.n_obs)]
    adata.obsm.clear()
    adata.uns.clear()

    positions = pd.DataFrame(
        {
            "barcode": [f"other{i}" for i in range(adata.n_obs)],
            "in_tissue": np.ones(adata.n_obs, dtype=int),
            "array_row": np.arange(adata.n_obs),
            "array_col": np.arange(adata.n_obs),
            "pxl_row_in_fullres": np.linspace(0, 100, adata.n_obs),
            "pxl_col_in_fullres": np.linspace(10, 110, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions_list.csv", index=False)

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_mtx_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    with pytest.raises(DataCompatibilityError, match="No matching barcodes"):
        await dl.load_spatial_data(str(visium_dir), "visium")


@pytest.mark.asyncio
async def test_load_visium_h5_path_calls_spatial_helpers(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    h5_path = tmp_path / "sample.h5"
    h5_path.write_text("h5")

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_ret = minimal_spatial_adata.copy()
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    calls: dict[str, int] = {"find": 0, "add": 0}

    def _find(_path: str) -> str:
        calls["find"] += 1
        return "/tmp/spatial"

    def _add(adata, _spatial_path: str):
        calls["add"] += 1
        return adata

    monkeypatch.setattr(dl, "_find_spatial_folder", _find)
    monkeypatch.setattr(dl, "_add_spatial_info_to_adata", _add)

    await dl.load_spatial_data(str(h5_path), "visium")

    assert calls == {"find": 1, "add": 1}


@pytest.mark.asyncio
async def test_load_visium_h5_path_requires_spatial_metadata(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    h5_path = tmp_path / "sample_without_spatial.h5"
    h5_path.write_text("h5")

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_ret = minimal_spatial_adata.copy()
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setattr(dl, "_find_spatial_folder", lambda _path: None)

    with pytest.raises(DataNotFoundError, match="Visium spatial metadata not found"):
        await dl.load_spatial_data(str(h5_path), "visium")


@pytest.mark.asyncio
async def test_load_visium_h5_path_errors_when_spatial_helper_fails(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """Visium H5 loader must fail-fast when spatial helper raises."""
    h5_path = tmp_path / "sample_with_bad_spatial.h5"
    h5_path.write_text("h5")

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_ret = minimal_spatial_adata.copy()
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setattr(dl, "_find_spatial_folder", lambda _path: "/tmp/spatial")
    monkeypatch.setattr(
        dl,
        "_add_spatial_info_to_adata",
        lambda _adata, _spatial_path: (_ for _ in ()).throw(RuntimeError("bad spatial")),
    )

    with pytest.raises(DataCompatibilityError, match="Visium spatial information failed"):
        await dl.load_spatial_data(str(h5_path), "visium")


@pytest.mark.asyncio
async def test_load_visium_h5ad_without_spatial_emits_warning(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
):
    h5ad_path = tmp_path / "sample.h5ad"
    h5ad_path.write_text("h5ad")

    adata = minimal_spatial_adata.copy()
    adata.obsm.clear()
    adata.uns.clear()

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_h5ad_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    with caplog.at_level("WARNING"):
        out = await dl.load_spatial_data(str(h5ad_path), "visium")

    assert out["spatial_coordinates_available"] is False
    assert "does not contain spatial information" in caplog.text


@pytest.mark.asyncio
async def test_load_visium_read_10x_h5_error_includes_helpful_suggestions(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    h5_path = tmp_path / "sample.h5"
    h5_path.write_text("h5")

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_exc = RuntimeError("No matching barcodes while read_10x_h5")
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    with pytest.raises(ProcessingError) as exc_info:
        await dl.load_spatial_data(str(h5_path), "visium")

    msg = str(exc_info.value)
    assert "Possible solutions" in msg
    assert "barcode format" in msg


@pytest.mark.asyncio
async def test_load_visium_read_10x_h5_error_reports_invalid_h5_guidance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    h5_path = tmp_path / "sample_bad.h5"
    h5_path.write_text("h5")

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_exc = RuntimeError("read_10x_h5 parser failure")
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    with pytest.raises(ProcessingError) as exc_info:
        await dl.load_spatial_data(str(h5_path), "visium")

    assert "valid 10x H5 file" in str(exc_info.value)


@pytest.mark.asyncio
async def test_load_visium_spatial_issue_error_adds_spatial_guidance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    h5_path = tmp_path / "sample.h5"
    h5_path.write_text("h5")

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_exc = RuntimeError("spatial metadata missing")
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    with pytest.raises(ProcessingError) as exc_info:
        await dl.load_spatial_data(str(h5_path), "visium")

    assert "Spatial data issue detected" in str(exc_info.value)


@pytest.mark.asyncio
async def test_load_visium_rejects_unsupported_visium_input_format(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    txt_path = tmp_path / "sample.txt"
    txt_path.write_text("plain")
    monkeypatch.setitem(sys.modules, "scanpy", _FakeScanpy())

    with pytest.raises(ParameterError, match="Unsupported file format for visium"):
        await dl.load_spatial_data(str(txt_path), "visium")


@pytest.mark.asyncio
async def test_load_visium_wraps_file_not_found_from_scanpy(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    visium_dir = tmp_path / "sample_visium_missing"
    visium_dir.mkdir()
    (visium_dir / "filtered_feature_bc_matrix.h5").write_text("h5")

    fake_scanpy = _FakeScanpy()
    fake_scanpy.read_visium = lambda _path: (_ for _ in ()).throw(
        FileNotFoundError("matrix missing")
    )
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    with pytest.raises(DataNotFoundError, match="File not found"):
        await dl.load_spatial_data(str(visium_dir), "visium")


@pytest.mark.asyncio
async def test_load_xenium_standard_h5_filters_to_cells_with_metadata(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    xen_dir = tmp_path / "xenium_std"
    xen_dir.mkdir()
    (xen_dir / "cell_feature_matrix.h5").write_text("h5")

    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]

    cells = pd.DataFrame(
        {
            "cell_id": [f"cell_{i}" for i in range(adata.n_obs // 2)],
            "x_centroid": np.linspace(0, 10, adata.n_obs // 2),
            "y_centroid": np.linspace(5, 15, adata.n_obs // 2),
            "transcript_counts": np.arange(adata.n_obs // 2),
        }
    )
    with gzip.open(xen_dir / "cells.csv.gz", "wt") as f:
        cells.to_csv(f, index=False)

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    out = await dl.load_spatial_data(str(xen_dir), "xenium")
    out_adata = out["adata"]

    assert out_adata.n_obs == adata.n_obs // 2
    assert "spatial" in out_adata.obsm
    assert "transcript_counts" in out_adata.obs


@pytest.mark.asyncio
async def test_load_xenium_zarr_format_branch_uses_zarr_loader(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    xen_dir = tmp_path / "xenium_zarr"
    xen_dir.mkdir()
    (xen_dir / "cell_feature_matrix.zarr.zip").write_text("zip")
    (xen_dir / "cells.zarr.zip").write_text("zip")

    fake_scanpy = _FakeScanpy()
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setattr(dl, "_load_xenium_zarr", lambda _path: minimal_spatial_adata.copy())

    out = await dl.load_spatial_data(str(xen_dir), "xenium")
    assert out["type"] == "xenium"
    assert out["n_cells"] == minimal_spatial_adata.n_obs


@pytest.mark.asyncio
async def test_load_xenium_standard_matrix_dir_uses_read_10x_mtx(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    xen_dir = tmp_path / "xenium_mtx"
    matrix_dir = xen_dir / "cell_feature_matrix"
    xen_dir.mkdir()
    matrix_dir.mkdir()
    with gzip.open(xen_dir / "cells.csv.gz", "wt") as f:
        pd.DataFrame(
            {
                "cell_id": [f"cell_{i}" for i in range(minimal_spatial_adata.n_obs)],
                "x_centroid": np.linspace(0, 10, minimal_spatial_adata.n_obs),
                "y_centroid": np.linspace(5, 15, minimal_spatial_adata.n_obs),
            }
        ).to_csv(f, index=False)

    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_mtx_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    out = await dl.load_spatial_data(str(xen_dir), "xenium")
    assert fake_scanpy.calls["read_10x_mtx"]
    assert out["adata"].n_obs == adata.n_obs


@pytest.mark.asyncio
async def test_load_xenium_standard_reads_parquet_when_available(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    xen_dir = tmp_path / "xenium_parquet"
    xen_dir.mkdir()
    (xen_dir / "cell_feature_matrix.h5").write_text("h5")
    (xen_dir / "cells.parquet").write_text("placeholder")

    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    monkeypatch.setattr(
        pd,
        "read_parquet",
        lambda _p: pd.DataFrame(
            {
                "cell_id": [f"cell_{i}" for i in range(adata.n_obs)],
                "x_centroid": np.linspace(0, 10, adata.n_obs),
                "y_centroid": np.linspace(5, 15, adata.n_obs),
            }
        ),
    )

    out = await dl.load_spatial_data(str(xen_dir), "xenium")
    assert out["adata"].n_obs == adata.n_obs


@pytest.mark.asyncio
async def test_load_xenium_standard_raises_when_centroid_columns_missing(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    xen_dir = tmp_path / "xenium_missing_centroids"
    xen_dir.mkdir()
    (xen_dir / "cell_feature_matrix.h5").write_text("h5")
    with gzip.open(xen_dir / "cells.csv.gz", "wt") as f:
        pd.DataFrame(
            {
                "cell_id": [f"cell_{i}" for i in range(minimal_spatial_adata.n_obs)],
                "some_other_col": np.arange(minimal_spatial_adata.n_obs),
            }
        ).to_csv(f, index=False)

    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    with pytest.raises(DataCompatibilityError, match="missing x_centroid/y_centroid"):
        await dl.load_spatial_data(str(xen_dir), "xenium")


@pytest.mark.asyncio
async def test_load_xenium_standard_raises_when_no_matching_cell_ids(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    xen_dir = tmp_path / "xenium_nomatch"
    xen_dir.mkdir()
    (xen_dir / "cell_feature_matrix.h5").write_text("h5")

    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"obs_{i}" for i in range(adata.n_obs)]

    cells = pd.DataFrame(
        {
            "cell_id": [f"other_{i}" for i in range(adata.n_obs)],
            "x_centroid": np.linspace(0, 10, adata.n_obs),
            "y_centroid": np.linspace(5, 15, adata.n_obs),
        }
    )
    with gzip.open(xen_dir / "cells.csv.gz", "wt") as f:
        cells.to_csv(f, index=False)

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    with pytest.raises(DataCompatibilityError, match="No matching cell IDs"):
        await dl.load_spatial_data(str(xen_dir), "xenium")


@pytest.mark.asyncio
async def test_load_xenium_raises_not_found_for_invalid_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    xen_dir = tmp_path / "xenium_empty"
    xen_dir.mkdir()

    fake_scanpy = _FakeScanpy()
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    with pytest.raises(DataNotFoundError, match="No valid Xenium data found"):
        await dl.load_spatial_data(str(xen_dir), "xenium")


@pytest.mark.asyncio
async def test_load_xenium_wraps_unexpected_errors_as_processing_error(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    xen_dir = tmp_path / "xenium_runtime_error"
    xen_dir.mkdir()
    (xen_dir / "cell_feature_matrix.h5").write_text("h5")
    with gzip.open(xen_dir / "cells.csv.gz", "wt") as f:
        pd.DataFrame({"cell_id": ["cell_0"]}).to_csv(f, index=False)

    adata = minimal_spatial_adata.copy()
    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_10x_h5_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setattr(pd, "read_csv", lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("csv boom")))

    with pytest.raises(ProcessingError, match="Error loading Xenium data"):
        await dl.load_spatial_data(str(xen_dir), "xenium")


def test_find_spatial_folder_checks_parent_directory_candidate(tmp_path: Path):
    sample_dir = tmp_path / "sample"
    sample_dir.mkdir()
    h5_path = sample_dir / "matrix.h5"
    h5_path.write_text("h5")

    parent_spatial = tmp_path / "spatial"
    parent_spatial.mkdir()
    (parent_spatial / "tissue_positions_list.csv").write_text("barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n")
    (parent_spatial / "scalefactors_json.json").write_text("{}")

    found = dl._find_spatial_folder(str(h5_path))
    assert found is not None
    assert os.path.normpath(found) == os.path.normpath(str(parent_spatial))


def test_add_spatial_info_supports_5_column_positions_and_suffix_removal(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"bc{i}" for i in range(adata.n_obs)]

    spatial_dir = tmp_path / "sample_c" / "spatial"
    spatial_dir.mkdir(parents=True)

    positions = pd.DataFrame(
        {
            0: [f"bc{i}-1" for i in range(adata.n_obs)],
            1: np.arange(adata.n_obs),
            2: np.arange(adata.n_obs),
            3: np.linspace(0, 50, adata.n_obs),
            4: np.linspace(10, 60, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions_list.csv", header=False, index=False)
    (spatial_dir / "scalefactors_json.json").write_text("{}")

    monkeypatch.setattr(dl, "is_available", lambda _name: False)

    out = dl._add_spatial_info_to_adata(adata, str(spatial_dir))

    assert out.n_obs == adata.n_obs
    assert "in_tissue" in out.obs
    assert "sample_c" in out.uns["spatial"]


def test_add_spatial_info_loads_images_when_pillow_available(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"bc{i}" for i in range(adata.n_obs)]

    spatial_dir = tmp_path / "sample_img" / "spatial"
    spatial_dir.mkdir(parents=True)

    positions = pd.DataFrame(
        {
            "barcode": adata.obs_names,
            "in_tissue": [1] * adata.n_obs,
            "array_row": np.arange(adata.n_obs),
            "array_col": np.arange(adata.n_obs),
            "pxl_row_in_fullres": np.linspace(0, 10, adata.n_obs),
            "pxl_col_in_fullres": np.linspace(5, 15, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions_list.csv", index=False)
    (spatial_dir / "scalefactors_json.json").write_text("{}")
    (spatial_dir / "tissue_hires_image.png").write_bytes(b"x")

    fake_pil = ModuleType("PIL")

    class _FakeImage:
        @staticmethod
        def open(_path):
            return np.zeros((4, 4, 3), dtype=np.uint8)

    fake_pil.Image = _FakeImage
    monkeypatch.setitem(sys.modules, "PIL", fake_pil)
    monkeypatch.setattr(dl, "is_available", lambda name: name == "Pillow")

    out = dl._add_spatial_info_to_adata(adata, str(spatial_dir))

    imgs = out.uns["spatial"]["sample_img"]["images"]
    assert "hires" in imgs
    assert imgs["hires"].shape == (4, 4, 3)


def test_add_spatial_info_supports_6_column_positions_without_header(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"bc{i}" for i in range(adata.n_obs)]

    spatial_dir = tmp_path / "sample_6col" / "spatial"
    spatial_dir.mkdir(parents=True)

    positions = pd.DataFrame(
        {
            0: adata.obs_names,
            1: np.ones(adata.n_obs, dtype=int),
            2: np.arange(adata.n_obs),
            3: np.arange(adata.n_obs),
            4: np.linspace(0, 10, adata.n_obs),
            5: np.linspace(5, 15, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions_list.csv", index=False, header=False)
    (spatial_dir / "scalefactors_json.json").write_text("{}")
    monkeypatch.setattr(dl, "is_available", lambda _name: False)

    out = dl._add_spatial_info_to_adata(adata, str(spatial_dir))
    assert out.obsm["spatial"].shape == (adata.n_obs, 2)


def test_add_spatial_info_uses_fallback_library_id_when_parent_missing(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"bc{i}" for i in range(adata.n_obs)]

    spatial_dir = tmp_path / "spatial"
    spatial_dir.mkdir()
    positions = pd.DataFrame(
        {
            "barcode": adata.obs_names,
            "in_tissue": [1] * adata.n_obs,
            "array_row": np.arange(adata.n_obs),
            "array_col": np.arange(adata.n_obs),
            "pxl_row_in_fullres": np.linspace(0, 10, adata.n_obs),
            "pxl_col_in_fullres": np.linspace(5, 15, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions_list.csv", index=False)
    (spatial_dir / "scalefactors_json.json").write_text("{}")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(dl, "is_available", lambda _name: False)

    out = dl._add_spatial_info_to_adata(adata, "spatial")
    assert "sample_1" in out.uns["spatial"]


def test_add_spatial_info_logs_warning_when_image_loading_fails(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
):
    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"bc{i}" for i in range(adata.n_obs)]

    spatial_dir = tmp_path / "sample_bad_image" / "spatial"
    spatial_dir.mkdir(parents=True)
    positions = pd.DataFrame(
        {
            "barcode": adata.obs_names,
            "in_tissue": [1] * adata.n_obs,
            "array_row": np.arange(adata.n_obs),
            "array_col": np.arange(adata.n_obs),
            "pxl_row_in_fullres": np.linspace(0, 10, adata.n_obs),
            "pxl_col_in_fullres": np.linspace(5, 15, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions_list.csv", index=False)
    (spatial_dir / "scalefactors_json.json").write_text("{}")
    (spatial_dir / "tissue_lowres_image.png").write_bytes(b"x")

    fake_pil = ModuleType("PIL")

    class _BadImage:
        @staticmethod
        def open(_path):
            raise RuntimeError("broken image")

    fake_pil.Image = _BadImage
    monkeypatch.setitem(sys.modules, "PIL", fake_pil)
    monkeypatch.setattr(dl, "is_available", lambda name: name == "Pillow")

    with caplog.at_level("WARNING"):
        out = dl._add_spatial_info_to_adata(adata, str(spatial_dir))

    assert "sample_bad_image" in out.uns["spatial"]
    assert "Could not load image tissue_lowres_image.png" in caplog.text


def test_add_spatial_info_raises_when_barcodes_cannot_be_aligned(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs_names = [f"obs{i}" for i in range(adata.n_obs)]

    spatial_dir = tmp_path / "sample_d" / "spatial"
    spatial_dir.mkdir(parents=True)

    positions = pd.DataFrame(
        {
            "barcode": [f"other{i}" for i in range(adata.n_obs)],
            "in_tissue": [1] * adata.n_obs,
            "array_row": np.arange(adata.n_obs),
            "array_col": np.arange(adata.n_obs),
            "pxl_row_in_fullres": np.linspace(0, 10, adata.n_obs),
            "pxl_col_in_fullres": np.linspace(5, 15, adata.n_obs),
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions_list.csv", index=False)
    (spatial_dir / "scalefactors_json.json").write_text("{}")

    monkeypatch.setattr(dl, "is_available", lambda _name: False)

    with pytest.raises(DataCompatibilityError, match="No matching barcodes"):
        dl._add_spatial_info_to_adata(adata, str(spatial_dir))


@pytest.mark.asyncio
async def test_load_generic_sets_tissue_image_available_when_nested_images_exist(
    minimal_spatial_adata,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    path = tmp_path / "generic.h5ad"
    path.write_text("h5ad")

    adata = minimal_spatial_adata.copy()
    adata.uns["spatial"] = {
        "sample1": {
            "images": {"hires": np.zeros((2, 2, 3), dtype=np.uint8)},
            "scalefactors": {},
        }
    }

    fake_scanpy = _FakeScanpy()
    fake_scanpy._read_h5ad_ret = adata
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    out = await dl.load_spatial_data(str(path), "generic")

    assert out["tissue_image_available"] is True
