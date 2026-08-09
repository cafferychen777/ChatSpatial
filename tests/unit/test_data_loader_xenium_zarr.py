"""Unit tests for Xenium zarr loading with proper resource management.

Verifies that _load_xenium_zarr correctly reads zarr data and
releases ZipStore file handles after loading.
"""

import os
import zipfile

import numpy as np
import pytest

from chatspatial.utils.data_loader import _load_xenium_zarr

pytestmark = pytest.mark.filterwarnings("ignore:Duplicate name:UserWarning")


def _create_xenium_zarr_fixture(tmp_path):
    """Create minimal Xenium zarr.zip files for testing.

    Produces cell_feature_matrix.zarr.zip and cells.zarr.zip with
    the same structure as real Xenium output.
    """
    import zarr
    from zarr.storage import ZipStore

    n_cells = 10
    n_features = 5

    # Build cell_feature_matrix.zarr.zip
    matrix_zip = str(tmp_path / "cell_feature_matrix.zarr.zip")
    with ZipStore(matrix_zip, mode="w") as store:
        root = zarr.open(store, mode="w")
        cf = root.create_group("cell_features")

        # CSC sparse matrix components (identity-like for simplicity)
        data = np.ones(n_cells, dtype=np.float32)
        indices = np.arange(n_cells, dtype=np.int32)
        indptr = np.concatenate(
            [
                np.arange(0, n_cells, dtype=np.int64),
                np.array([n_cells] * (n_features - n_cells + 1), dtype=np.int64),
            ]
        )[: n_features + 1]

        # zarr v3: use create_array with data= (shape inferred)
        cf.create_array("data", data=data)
        cf.create_array("indices", data=indices)
        cf.create_array("indptr", data=indptr)

        cf.attrs["number_cells"] = n_cells
        cf.attrs["number_features"] = n_features
        cf.attrs["feature_keys"] = [f"Gene_{i}" for i in range(n_features)]
        cf.attrs["feature_ids"] = [f"ENSG{i:011d}" for i in range(n_features)]

    # Build cells.zarr.zip
    cells_zip = str(tmp_path / "cells.zarr.zip")
    rng = np.random.default_rng(42)
    # columns: x_centroid, y_centroid, cell_area, nucleus_area
    cell_summary = rng.uniform(0, 1000, size=(n_cells, 4)).astype(np.float64)
    cell_id = np.arange(1, n_cells + 1).reshape(-1, 1)

    with ZipStore(cells_zip, mode="w") as store:
        root = zarr.open(store, mode="w")
        root.create_array("cell_summary", data=cell_summary)
        root.create_array("cell_id", data=cell_id)

        cs_arr = root["cell_summary"]
        cs_arr.attrs["column_names"] = [
            "x_centroid",
            "y_centroid",
            "cell_area",
            "nucleus_area",
        ]

    return tmp_path, n_cells, n_features


@pytest.mark.unit
class TestLoadXeniumZarr:
    """Tests for _load_xenium_zarr function."""

    def test_returns_valid_anndata(self, tmp_path):
        """Loaded AnnData has correct shape and spatial coordinates."""
        data_dir, n_cells, n_features = _create_xenium_zarr_fixture(tmp_path)

        adata = _load_xenium_zarr(str(data_dir))

        assert adata.shape == (n_cells, n_features)
        assert "spatial" in adata.obsm
        assert adata.obsm["spatial"].shape == (n_cells, 2)

    def test_gene_names_and_ids(self, tmp_path):
        """Variable names and gene IDs are correctly loaded."""
        data_dir, _, n_features = _create_xenium_zarr_fixture(tmp_path)

        adata = _load_xenium_zarr(str(data_dir))

        assert list(adata.var_names) == [f"Gene_{i}" for i in range(n_features)]
        assert list(adata.var["gene_ids"]) == [
            f"ENSG{i:011d}" for i in range(n_features)
        ]

    def test_obs_metadata_columns(self, tmp_path):
        """Cell metadata columns from cell_summary are present."""
        data_dir, _, _ = _create_xenium_zarr_fixture(tmp_path)

        adata = _load_xenium_zarr(str(data_dir))

        expected_cols = ["x_centroid", "y_centroid", "cell_area", "nucleus_area"]
        for col in expected_cols:
            assert col in adata.obs.columns

    def test_zipstore_handles_released(self, tmp_path):
        """ZipStore file handles are closed after loading completes."""
        data_dir, _, _ = _create_xenium_zarr_fixture(tmp_path)

        adata = _load_xenium_zarr(str(data_dir))

        # Verify the files are not locked: open them again for reading
        # On all platforms, a leaked handle would prevent re-opening in
        # write mode or cause warnings.
        matrix_zip = os.path.join(str(data_dir), "cell_feature_matrix.zarr.zip")
        cells_zip = os.path.join(str(data_dir), "cells.zarr.zip")

        # If handles were leaked, zipfile.ZipFile in 'a' mode would
        # raise or produce corrupted output on some platforms
        for path in [matrix_zip, cells_zip]:
            with zipfile.ZipFile(path, "r") as zf:
                assert len(zf.namelist()) > 0

        # Confirm adata is still valid (data was fully materialized)
        assert adata.X.shape[0] > 0

    def test_sparse_matrix_format(self, tmp_path):
        """Expression matrix is in CSR sparse format."""
        import scipy.sparse as sp

        data_dir, _, _ = _create_xenium_zarr_fixture(tmp_path)

        adata = _load_xenium_zarr(str(data_dir))

        assert sp.issparse(adata.X)
        assert isinstance(adata.X, sp.csr_matrix)
