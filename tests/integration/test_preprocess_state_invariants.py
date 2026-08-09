"""Integration tests for preprocessing state invariants."""

from __future__ import annotations

import pytest

from chatspatial.models.data import PreprocessingParameters
from chatspatial.server import data_manager, load_data, preprocess_data


@pytest.mark.integration
@pytest.mark.asyncio
async def test_preprocess_creates_expected_state_artifacts(
    spatial_dataset_path, reset_data_manager
):
    dataset = await load_data(
        str(spatial_dataset_path), "generic", name="prep_invariants"
    )

    params = PreprocessingParameters(
        filter_genes_min_cells=1,
        filter_cells_min_genes=1,
        normalization="log",
        n_hvgs=15,
        subsample_genes=15,
        remove_mito_genes=False,
    )
    result = await preprocess_data(dataset.id, params=params)

    stored = await data_manager.get_dataset(dataset.id)
    adata = stored["adata"]

    assert result.n_cells == adata.n_obs
    assert result.n_genes == adata.n_vars

    # Invariants: raw + counts layer + preprocessing metadata
    assert adata.raw is not None
    assert adata.raw.n_vars >= adata.n_vars
    assert "counts" in adata.layers
    assert adata.layers["counts"].shape == adata.X.shape

    assert "preprocessing" in adata.uns
    prep_meta = adata.uns["preprocessing"]
    assert prep_meta["completed"] is True
    assert prep_meta["raw_preserved"] is True
    assert prep_meta["counts_layer"] is True


@pytest.mark.integration
@pytest.mark.asyncio
async def test_preprocess_does_not_overwrite_existing_raw(
    minimal_spatial_adata, write_h5ad_dataset, reset_data_manager
):
    adata = minimal_spatial_adata.copy()

    # Pre-existing raw with sentinel metadata to detect unintended overwrite.
    # Keep same shape as X to satisfy data loader's counts-layer expectations.
    raw_adata = adata.copy()
    raw_adata.var["raw_marker"] = "keep_me"
    adata.raw = raw_adata
    dataset_path = write_h5ad_dataset(adata, "with_existing_raw.h5ad")

    dataset = await load_data(str(dataset_path), "generic", name="existing_raw")

    params = PreprocessingParameters(
        filter_genes_min_cells=1,
        filter_cells_min_genes=1,
        normalization="log",
        n_hvgs=12,
        subsample_genes=12,
        remove_mito_genes=False,
    )
    await preprocess_data(dataset.id, params=params)

    stored = await data_manager.get_dataset(dataset.id)
    adata_after = stored["adata"]

    assert adata_after.raw is not None
    assert adata_after.raw.n_vars == 24
    assert "raw_marker" in adata_after.raw.var.columns
    assert set(adata_after.raw.var["raw_marker"]) == {"keep_me"}
    assert "counts" in adata_after.layers
