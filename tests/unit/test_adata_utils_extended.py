"""Extended unit tests for high-impact AnnData utility paths."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from chatspatial.utils import adata_utils as au
from chatspatial.utils.exceptions import DataError, DataNotFoundError, ParameterError


class _WarnCtx:
    def __init__(self) -> None:
        self.messages: list[str] = []

    async def warning(self, msg: str) -> None:
        self.messages.append(msg)


def test_discovery_helpers_find_cell_cluster_and_batch_keys(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["T"] * adata.n_obs
    adata.obs["leiden"] = ["0"] * adata.n_obs
    adata.obs["sample_id"] = ["s1"] * adata.n_obs

    assert au.get_cell_type_key(adata) == "cell_type"
    assert au.get_cluster_key(adata) == "leiden"
    assert au.get_batch_key(adata) == "sample_id"

    empty = minimal_spatial_adata.copy()
    for col in ["group"]:
        if col in empty.obs:
            del empty.obs[col]
    assert au.get_cell_type_key(empty) is None
    assert au.get_cluster_key(empty) is None
    assert au.get_batch_key(empty) is None


def test_discovery_helpers_use_deterministic_priority_when_multiple_keys_exist(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["celltype"] = ["A"] * adata.n_obs
    adata.obs["cell_type"] = ["B"] * adata.n_obs
    adata.obs["cluster"] = ["x"] * adata.n_obs
    adata.obs["leiden"] = ["0"] * adata.n_obs
    adata.obs["sample"] = ["s"] * adata.n_obs
    adata.obs["batch"] = ["b"] * adata.n_obs
    adata.obsm["coords"] = adata.obsm["spatial"].copy()

    # Priority order should be stable and explicit.
    assert au.get_cell_type_key(adata) == "cell_type"
    assert au.get_cluster_key(adata) == "leiden"
    assert au.get_batch_key(adata) == "batch"
    assert au.get_spatial_key(adata) == "spatial"


def test_sample_expression_values_handles_dense_sparse_and_empty_sparse(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()

    dense_sample = au.sample_expression_values(adata, n_samples=5)
    assert dense_sample.shape[0] == 5

    adata_sparse = adata.copy()
    adata_sparse.X = sparse.csr_matrix(np.asarray(adata_sparse.X))
    sparse_sample = au.sample_expression_values(adata_sparse, n_samples=7)
    assert sparse_sample.shape[0] <= 7

    adata_empty = adata.copy()
    adata_empty.X = sparse.csr_matrix((adata_empty.n_obs, adata_empty.n_vars), dtype=float)
    empty_sample = au.sample_expression_values(adata_empty, n_samples=3)
    assert empty_sample.shape[0] == 3
    assert np.all(empty_sample == 0)


def test_require_spatial_coords_validates_nan_inf_and_missing_key(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obsm["spatial"][0, 0] = np.nan
    with pytest.raises(DataError, match="contain NaN"):
        au.require_spatial_coords(adata)

    adata = minimal_spatial_adata.copy()
    adata.obsm["spatial"][0, 0] = np.inf
    with pytest.raises(DataError, match="infinite"):
        au.require_spatial_coords(adata)

    with pytest.raises(DataError, match="not found in adata.obsm"):
        au.require_spatial_coords(minimal_spatial_adata.copy(), spatial_key="not_here")


def test_require_spatial_coords_obs_xy_nan_and_missing_sources(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    adata.obs["x"] = np.arange(adata.n_obs, dtype=float)
    adata.obs["y"] = np.arange(adata.n_obs, dtype=float)
    adata.obs.loc[0, "x"] = np.nan
    with pytest.raises(DataError, match="contain NaN"):
        au.require_spatial_coords(adata)

    adata2 = minimal_spatial_adata.copy()
    del adata2.obsm["spatial"]
    if "x" in adata2.obs:
        del adata2.obs["x"]
    if "y" in adata2.obs:
        del adata2.obs["y"]
    with pytest.raises(DataError, match="No spatial coordinates found"):
        au.require_spatial_coords(adata2)


def test_require_spatial_coords_rejects_less_than_two_dimensions(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obsm["spatial"] = np.ones((adata.n_obs, 1), dtype=float)
    with pytest.raises(DataError, match="at least 2 dimensions"):
        au.require_spatial_coords(adata)


def test_validate_obs_and_var_column_raise_useful_errors(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()

    with pytest.raises(DataError, match="Cell type not found"):
        au.validate_obs_column(adata, "cell_type", "Cell type")

    with pytest.raises(DataError, match="Marker not found"):
        au.validate_var_column(adata, "marker", "Marker")


def test_validate_adata_basics_check_empty_ratio_for_dense_and_sparse(minimal_spatial_adata):
    dense = minimal_spatial_adata.copy()
    dense.X = np.asarray(dense.X)
    dense.X[:20, :] = 0.0
    with pytest.raises(DataError, match="cells .* zero expression"):
        au.validate_adata_basics(dense, check_empty_ratio=True, max_empty_obs_ratio=0.1)

    sparse_adata = minimal_spatial_adata.copy()
    sparse_adata.X = sparse.csr_matrix(np.asarray(sparse_adata.X))
    sparse_matrix = sparse_adata.X.tolil()
    sparse_matrix[:, :18] = 0.0
    sparse_adata.X = sparse_matrix.tocsr()
    with pytest.raises(DataError, match="genes .* zero expression"):
        au.validate_adata_basics(sparse_adata, check_empty_ratio=True, max_empty_vars_ratio=0.5)


def test_validate_adata_basics_rejects_none_and_too_small_dataset(minimal_spatial_adata):
    with pytest.raises(DataError, match="cannot be None"):
        au.validate_adata_basics(None)

    tiny_obs = minimal_spatial_adata[:1, :].copy()
    with pytest.raises(DataError, match="observations"):
        au.validate_adata_basics(tiny_obs, min_obs=2)

    tiny_var = minimal_spatial_adata[:, :1].copy()
    with pytest.raises(DataError, match="variables"):
        au.validate_adata_basics(tiny_var, min_vars=2)


def test_ensure_categorical_converts_and_rejects_missing_column(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["cluster"] = ["A", "B"] * (adata.n_obs // 2)

    au.ensure_categorical(adata, "cluster")
    assert pd.api.types.is_categorical_dtype(adata.obs["cluster"])

    with pytest.raises(DataError, match="Column 'missing' not found"):
        au.ensure_categorical(adata, "missing")


def test_shallow_copy_adata_shares_expression_but_copies_metadata(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.layers["counts"] = np.asarray(adata.X).copy()
    adata.varm["loadings"] = np.ones((adata.n_vars, 2), dtype=float)
    adata.uns["meta"] = {"v": 1}
    adata.raw = adata.copy()

    shallow = au.shallow_copy_adata(adata)

    assert shallow is not adata
    assert shallow.X is adata.X
    assert shallow.layers["counts"] is adata.layers["counts"]
    assert shallow.varm["loadings"] is adata.varm["loadings"]
    assert shallow.raw is adata.raw
    assert shallow.obs is not adata.obs
    assert shallow.var is not adata.var
    assert shallow.uns is not adata.uns


def _make_velovi_workspace(adata, *, velocity_key: str = "velocity_velovi"):
    vel = adata[:, :6].copy()
    vel.layers[velocity_key] = np.arange(
        vel.n_obs * vel.n_vars, dtype=np.float32
    ).reshape(vel.n_obs, vel.n_vars)
    vel.layers["Ms"] = np.full((vel.n_obs, vel.n_vars), 1.0, dtype=np.float32)
    vel.layers["Mu"] = np.full((vel.n_obs, vel.n_vars), 2.0, dtype=np.float32)
    diagonal = np.arange(1, vel.n_obs + 1, dtype=np.float32)
    vel.obsp["connectivities"] = sparse.diags(diagonal, format="csr")
    vel.obsp["distances"] = sparse.diags(diagonal * 2, format="csr")
    return vel


def test_store_and_reconstruct_velovi_essential_data_roundtrip(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    vel = _make_velovi_workspace(adata)

    au.store_velovi_essential_data(adata, vel)
    assert au.has_velovi_essential_data(adata) is True
    assert adata.uns["velovi_obs_names"] == adata.obs_names.tolist()

    reconstructed = au.reconstruct_velovi_adata(adata)
    assert reconstructed.n_obs == adata.n_obs
    assert reconstructed.n_vars == vel.n_vars
    assert "velocity" in reconstructed.layers
    assert "Ms" in reconstructed.layers
    assert "Mu" in reconstructed.layers
    assert "connectivities" in reconstructed.obsp
    assert "distances" in reconstructed.obsp
    assert "neighbors" in reconstructed.uns
    assert "spatial" in reconstructed.obsm


def test_store_velovi_essential_data_uses_velocity_fallback(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    vel = _make_velovi_workspace(adata, velocity_key="velocity")
    adata.uns["velovi_Mu"] = np.array([[99.0]])
    adata.uns["velovi_distances"] = sparse.csr_matrix([[99.0]])
    del vel.layers["Mu"]
    del vel.obsp["distances"]
    au.store_velovi_essential_data(adata, vel)
    assert "velovi_velocity" in adata.uns
    assert np.array_equal(adata.uns["velovi_velocity"], vel.layers["velocity"])
    assert "velovi_Mu" not in adata.uns
    assert "velovi_distances" not in adata.uns


def test_store_velovi_essential_data_rejects_reordered_workspace(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    vel = _make_velovi_workspace(adata)[::-1].copy()
    with pytest.raises(DataError, match="do not match the target dataset order"):
        au.store_velovi_essential_data(adata, vel)
    assert au.has_velovi_essential_data(adata) is False


def test_reconstruct_velovi_adata_aligns_reordered_cell_subset(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    vel = _make_velovi_workspace(adata)
    au.store_velovi_essential_data(adata, vel)
    positions = np.array([7, 2, 11, 0])
    subset = adata[positions].copy()

    reconstructed = au.reconstruct_velovi_adata(subset)

    assert reconstructed.obs_names.equals(subset.obs_names)
    np.testing.assert_array_equal(
        reconstructed.layers["velocity"],
        vel.layers["velocity_velovi"][positions, :],
    )
    expected_graph = vel.obsp["connectivities"][positions, :][:, positions]
    np.testing.assert_array_equal(
        reconstructed.obsp["connectivities"].toarray(),
        expected_graph.toarray(),
    )


def test_reconstruct_velovi_adata_rejects_unknown_current_cells(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    vel = _make_velovi_workspace(adata)
    au.store_velovi_essential_data(adata, vel)
    adata.obs_names = ["new-cell", *adata.obs_names[1:]]

    with pytest.raises(DataError, match="not present in the VELOVI analysis"):
        au.reconstruct_velovi_adata(adata)


def test_reconstruct_velovi_adata_rejects_corrupted_matrix_shape(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    au.store_velovi_essential_data(adata, _make_velovi_workspace(adata))
    adata.uns["velovi_velocity"] = adata.uns["velovi_velocity"][:-1, :]

    with pytest.raises(DataError, match="velovi_velocity has shape"):
        au.reconstruct_velovi_adata(adata)


def test_velovi_payload_survives_h5ad_roundtrip(
    minimal_spatial_adata, tmp_path
):
    adata = minimal_spatial_adata.copy()
    vel = _make_velovi_workspace(adata)
    au.store_velovi_essential_data(adata, vel)
    adata.obsm["velovi_latent_time"] = pd.DataFrame(
        np.arange(adata.n_obs * vel.n_vars, dtype=np.float32).reshape(
            adata.n_obs, vel.n_vars
        ),
        index=adata.obs_names,
        columns=vel.var_names,
    )
    output_path = tmp_path / "velovi_roundtrip.h5ad"
    adata.write_h5ad(output_path)

    loaded = ad.read_h5ad(output_path)
    positions = np.array([9, 3, 0])
    subset = loaded[positions].copy()
    reconstructed = au.reconstruct_velovi_adata(subset)

    assert isinstance(subset.obsm["velovi_latent_time"], pd.DataFrame)
    assert subset.obsm["velovi_latent_time"].columns.equals(vel.var_names)
    np.testing.assert_array_equal(
        reconstructed.layers["velocity"],
        vel.layers["velocity_velovi"][positions, :],
    )


def test_reconstruct_velovi_adata_requires_essential_keys(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    assert au.has_velovi_essential_data(adata) is False
    with pytest.raises(DataNotFoundError, match="Missing VELOVI reconstruction data"):
        au.reconstruct_velovi_adata(adata)

    adata.uns["velovi_obs_names"] = adata.obs_names.tolist()
    adata.uns["velovi_gene_names"] = ["g1", "g2"]
    with pytest.raises(DataNotFoundError, match="velovi_velocity"):
        au.reconstruct_velovi_adata(adata)


def test_standardize_adata_moves_spatial_and_casts_known_categories(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obsm["coords"] = adata.obsm["spatial"].copy()
    del adata.obsm["spatial"]
    adata.obs["cell_type"] = ["T", "B"] * (adata.n_obs // 2)

    out = au.standardize_adata(adata, copy=True)

    assert "spatial" in out.obsm
    assert pd.api.types.is_categorical_dtype(out.obs["cell_type"])
    assert "spatial" not in adata.obsm


def test_standardize_adata_ignores_obs_xy_conversion_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    adata.obs["x"] = ["bad"] * adata.n_obs
    adata.obs["y"] = ["bad"] * adata.n_obs

    monkeypatch.setattr(au.pd, "to_numeric", lambda *_a, **_k: (_ for _ in ()).throw(ValueError("boom")))
    out = au.standardize_adata(adata, copy=True)
    assert "spatial" not in out.obsm


def test_standardize_adata_moves_numeric_obs_xy_to_spatial(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obsm.pop("spatial", None)
    adata.obs["x"] = [str(i) for i in range(adata.n_obs)]
    adata.obs["y"] = [str(i + 1) for i in range(adata.n_obs)]

    out = au.standardize_adata(adata, copy=True)
    assert "spatial" in out.obsm
    assert out.obsm["spatial"].dtype == np.float64
    assert out.obsm["spatial"].shape == (adata.n_obs, 2)


def test_validate_adata_required_keys_and_integrity_checks(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()

    with pytest.raises(DataError, match="Missing required keys"):
        au.validate_adata(adata, required_keys={"obs": ["missing_col"]})

    # Required keys pass, but detailed checks should fail with aggregated issues.
    adata.layers.pop("spliced", None)
    adata.layers.pop("unspliced", None)
    adata.obsm["spatial"] = np.ones((adata.n_obs, 2), dtype=float)
    with pytest.raises(DataError, match="Validation failed"):
        au.validate_adata(
            adata,
            required_keys={},
            check_spatial=True,
            check_velocity=True,
        )


def test_validate_adata_collects_missing_from_none_dict_and_plain_attrs(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    with pytest.raises(DataError, match="Missing required keys"):
        au.validate_adata(
            adata,
            required_keys={
                "missing_category": "k1",
                "obsm": ["not_here"],
                "X": ["not_a_dict"],
            },
        )


def test_validate_adata_spatial_and_velocity_issue_branches(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    with pytest.raises(DataError, match="Validation failed"):
        au.validate_adata(adata, required_keys={}, check_spatial=True)

    adata2 = minimal_spatial_adata.copy()
    adata2.obsm["spatial"] = np.full((adata2.n_obs, 1), np.nan, dtype=float)
    with pytest.raises(DataError, match="Validation failed"):
        au.validate_adata(adata2, required_keys={}, check_spatial=True)

    adata3 = minimal_spatial_adata.copy()
    adata3.layers["spliced"] = sparse.csr_matrix((adata3.n_obs, adata3.n_vars), dtype=float)
    unspliced = sparse.csr_matrix(np.zeros((adata3.n_obs, adata3.n_vars), dtype=float))
    unspliced = unspliced.tolil()
    unspliced[0, 0] = np.nan
    adata3.layers["unspliced"] = unspliced.tocsr()
    with pytest.raises(DataError, match="Validation failed"):
        au.validate_adata(adata3, required_keys={}, check_velocity=True)


def test_validate_adata_velocity_dense_issue_branches(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.layers["spliced"] = np.zeros((adata.n_obs, adata.n_vars), dtype=float)
    unspliced = np.ones((adata.n_obs, adata.n_vars), dtype=float)
    unspliced[0, 0] = np.nan
    adata.layers["unspliced"] = unspliced

    with pytest.raises(DataError, match="spliced.*all zeros"):
        au.validate_adata(adata, required_keys={}, check_velocity=True)
    with pytest.raises(DataError, match="unspliced.*NaN"):
        au.validate_adata(adata, required_keys={}, check_velocity=True)


def test_store_analysis_metadata_and_get_parameter(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()

    au.store_analysis_metadata(
        adata,
        analysis_name="spatial_stats_neighborhood",
        method="squidpy",
        parameters={"cluster_key": "leiden"},
        results_keys={"uns": ["neighborhood_enrichment"]},
        statistics={"n_clusters": 3},
        species="human",
        database="consensus",
        reference_info={"source": "ref_atlas_v1"},
    )

    assert "spatial_stats_neighborhood_metadata" in adata.uns
    assert adata.uns["spatial_stats_neighborhood_metadata"]["reference_info"] == {
        "source": "ref_atlas_v1"
    }
    assert (
        au.get_analysis_parameter(
            adata,
            analysis_name="spatial_stats_neighborhood",
            parameter_name="cluster_key",
            default=None,
        )
        == "leiden"
    )
    assert (
        au.get_analysis_parameter(
            adata,
            analysis_name="missing_analysis",
            parameter_name="cluster_key",
            default="fallback",
        )
        == "fallback"
    )
    adata.uns["no_params_metadata"] = {"method": "x"}
    assert (
        au.get_analysis_parameter(
            adata,
            analysis_name="no_params",
            parameter_name="cluster_key",
            default="fallback2",
        )
        == "fallback2"
    )


def test_hvg_and_gene_selection_paths(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.var["highly_variable"] = False
    adata.var.loc[adata.var_names[:5], "highly_variable"] = True

    hvgs = au.get_highly_variable_genes(adata, max_genes=3)
    assert len(hvgs) == 3

    selected = au.select_genes_for_analysis(adata, genes=["gene_0", "missing"], n_genes=2)
    assert selected == ["gene_0"]

    with pytest.raises(DataError, match="Did you mean"):
        au.select_genes_for_analysis(adata, genes=["gene_00"], n_genes=2)

    adata_no_hvg = minimal_spatial_adata.copy()
    out = au.select_genes_for_analysis(
        adata_no_hvg,
        genes=None,
        n_genes=2,
        require_hvg=False,
    )
    assert out == []


def test_select_genes_for_analysis_hvg_and_require_hvg_error(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.var["highly_variable"] = False
    adata.var.loc[adata.var_names[:4], "highly_variable"] = True

    hvgs = au.select_genes_for_analysis(adata, genes=None, n_genes=3, require_hvg=True)
    assert hvgs == adata.var_names[:3].tolist()

    adata_no_hvg = minimal_spatial_adata.copy()
    with pytest.raises(DataError, match="Highly variable genes \\(HVG\\) required"):
        au.select_genes_for_analysis(
            adata_no_hvg,
            genes=None,
            n_genes=3,
            require_hvg=True,
        )


def test_hvg_fallback_to_variance_for_sparse_and_dense(minimal_spatial_adata):
    dense = minimal_spatial_adata.copy()
    dense.var = dense.var.drop(columns=[c for c in dense.var.columns if c == "highly_variable"], errors="ignore")
    dense_hvgs = au.get_highly_variable_genes(dense, max_genes=4, fallback_to_variance=True)
    assert len(dense_hvgs) == 4

    sparse_adata = minimal_spatial_adata.copy()
    sparse_adata.var = sparse_adata.var.drop(columns=[c for c in sparse_adata.var.columns if c == "highly_variable"], errors="ignore")
    sparse_adata.X = sparse.csr_matrix(np.asarray(sparse_adata.X))
    sparse_hvgs = au.get_highly_variable_genes(sparse_adata, max_genes=4, fallback_to_variance=True)
    assert len(sparse_hvgs) == 4


def test_hvg_fallback_to_variance_returns_descending_ranked_genes():
    import anndata as ad

    X = np.array(
        [
            [1.0, 10.0, 3.0, 0.0],
            [1.0, 0.0, 3.0, 5.0],
            [1.0, 10.0, 3.0, 0.0],
            [1.0, 0.0, 3.0, 5.0],
        ]
    )
    adata = ad.AnnData(X)
    adata.var_names = ["g0", "g1", "g2", "g3"]

    out = au.get_highly_variable_genes(adata, max_genes=3, fallback_to_variance=True)
    assert out == ["g1", "g3", "g2"]


def test_make_unique_names_and_var_name_fixers(minimal_spatial_adata):
    assert au.make_unique_names(["A", "B", "A", "A"]) == ["A", "B", "A-1", "A-2"]

    adata = minimal_spatial_adata.copy()
    adata.var_names = [f"gene_{i // 2}" for i in range(adata.n_vars)]
    n_fixed = au.ensure_unique_var_names(adata)
    assert n_fixed > 0
    assert adata.var_names.is_unique


@pytest.mark.asyncio
async def test_ensure_unique_var_names_async_warns(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.var_names = ["dup"] * adata.n_vars
    ctx = _WarnCtx()

    n_fixed = await au.ensure_unique_var_names_async(adata, ctx, label="reference")

    assert n_fixed == adata.n_vars - 1
    assert any("duplicate gene names" in m for m in ctx.messages)


def test_check_is_integer_counts_and_ensure_counts_layer(minimal_spatial_adata):
    dense = np.array([[1.0, 2.0], [3.0, 4.0]])
    assert au.check_is_integer_counts(dense) == (True, False, False)

    non_int = np.array([[1.2, 2.0], [3.0, -4.0]])
    is_int, has_neg, has_dec = au.check_is_integer_counts(non_int)
    assert is_int is False
    assert has_neg is True
    assert has_dec is True
    assert au.check_is_integer_counts(sparse.csr_matrix(dense)) == (True, False, False)
    assert au.check_is_integer_counts(sparse.csr_matrix((3, 3), dtype=float)) == (
        True,
        False,
        False,
    )

    sparse_non_int = sparse.csr_matrix(np.array([[0.0, 1.5], [0.0, 2.0]]))
    is_int_sparse, has_neg_sparse, has_dec_sparse = au.check_is_integer_counts(
        sparse_non_int
    )
    assert is_int_sparse is False
    assert has_neg_sparse is False
    assert has_dec_sparse is True

    sparse_negative = sparse.csr_matrix(np.array([[0.0, -1.0], [0.0, 2.0]]))
    is_int_neg, has_neg_neg, has_dec_neg = au.check_is_integer_counts(sparse_negative)
    assert is_int_neg is False
    assert has_neg_neg is True
    assert has_dec_neg is False

    adata = minimal_spatial_adata.copy()
    adata.layers["counts"] = np.asarray(adata.X).copy()
    created = au.ensure_counts_layer(adata)
    assert created is False

    adata2 = minimal_spatial_adata.copy()
    adata2.raw = adata2.copy()
    adata2 = adata2[:, :5].copy()
    created = au.ensure_counts_layer(adata2)
    assert created is True
    assert adata2.layers["counts"].shape == (adata2.n_obs, adata2.n_vars)

    adata3 = minimal_spatial_adata.copy()
    adata3.raw = None
    adata3.layers.clear()
    # X contains integer counts (Poisson) → should create layer from X
    created = au.ensure_counts_layer(adata3)
    assert created is True
    assert "counts" in adata3.layers


def test_get_raw_data_source_skips_normalized_raw_and_falls_back_to_counts(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    raw_non_integer = adata.copy()
    raw_non_integer.X = np.asarray(raw_non_integer.X, dtype=float) + 0.25
    adata.raw = raw_non_integer
    adata.layers["counts"] = np.asarray(adata.X, dtype=np.int64)

    out = au.get_raw_data_source(
        adata, prefer_complete_genes=True, require_integer_counts=True
    )
    assert out.source == "counts_layer"


def test_get_gene_and_genes_expression_paths(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.layers["counts"] = np.asarray(adata.X).copy()

    one = au.get_gene_expression(adata, "gene_0")
    assert one.shape == (adata.n_obs,)

    one_counts = au.get_gene_expression(adata, "gene_0", layer="counts")
    assert np.array_equal(one_counts, one)

    many = au.get_genes_expression(adata, ["gene_0", "gene_1"], layer="counts")
    assert many.shape == (adata.n_obs, 2)

    with pytest.raises(DataError, match="Gene 'missing' not found"):
        au.get_gene_expression(adata, "missing")

    with pytest.raises(DataError, match="Layer 'bad' not found"):
        au.get_gene_expression(adata, "gene_0", layer="bad")

    with pytest.raises(DataError, match="Genes not found"):
        au.get_genes_expression(adata, ["gene_0", "missing"])


def test_get_genes_expression_additional_branches(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    out = au.get_genes_expression(adata, ["gene_0", "gene_1"])
    assert out.shape == (adata.n_obs, 2)

    with pytest.raises(DataError, match="Layer 'bad' not found"):
        au.get_genes_expression(adata, ["gene_0"], layer="bad")


def test_profiles_and_overlap_helpers(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["num_col"] = np.arange(adata.n_obs)
    adata.obs["cat_col"] = [f"v{i}" for i in range(adata.n_obs)]

    obs_profile = au.get_column_profile(adata, layer="obs")
    names = {item["name"]: item for item in obs_profile}
    assert names["num_col"]["dtype"] == "numerical"
    assert names["cat_col"]["dtype"] == "categorical"
    assert len(names["cat_col"]["sample_values"]) == 5

    profile = au.get_adata_profile(adata)
    assert "obs_columns" in profile
    assert "top_expressed_genes" in profile

    common = au.find_common_genes(["A", "B", "C"], ["B", "C"], ["C", "D"])
    assert common == ["C"]

    ordered_common = au.find_common_genes(
        ["G3", "G1", "G2", "G1"],
        ["G1", "G2", "G3"],
        ["G2", "G3"],
    )
    assert ordered_common == ["G3", "G2"]

    with pytest.raises(ParameterError, match="at least 2 gene collections"):
        au.find_common_genes(["A"])

    with pytest.raises(DataError, match="Insufficient gene overlap"):
        au.validate_gene_overlap(["A", "B"], source_n_genes=200, target_n_genes=150, min_genes=10)


def test_profiles_large_cardinality_and_gene_profile_fallback(minimal_spatial_adata):
    import anndata as ad

    big = ad.AnnData(np.ones((5, 120), dtype=float))
    big.var_names = [f"g{i}" for i in range(big.n_vars)]
    big.var["huge_cat"] = [f"cat_{i}" for i in range(big.n_vars)]
    var_profile = au.get_column_profile(big, layer="var")
    huge = next(item for item in var_profile if item["name"] == "huge_cat")
    assert len(huge["sample_values"]) == 3

    class _BrokenMeanArray(np.ndarray):
        def mean(self, *args, **kwargs):  # type: ignore[override]
            raise RuntimeError("boom")

    adata = minimal_spatial_adata.copy()
    adata.X = np.asarray(adata.X).view(_BrokenMeanArray)
    _top_hvg, top_expr = au.get_gene_profile(adata)
    assert len(top_expr) == min(10, adata.n_vars)


def test_get_gene_profile_returns_top_expression_in_descending_order():
    import anndata as ad

    expr = np.array(
        [
            [1.0, 10.0, 3.0, 4.0],
            [1.0, 10.0, 3.0, 4.0],
        ],
        dtype=float,
    )
    adata = ad.AnnData(expr)
    adata.var_names = ["g0", "g1", "g2", "g3"]

    _top_hvg, top_expr = au.get_gene_profile(adata)
    assert top_expr[:4] == ["g1", "g3", "g2", "g0"]


def test_get_raw_data_source_raw_error_and_counts_skip_path(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.raw = adata.copy()
    adata.layers["counts"] = np.asarray(adata.X, dtype=float) + 0.25  # non-integer, skip
    adata.X = np.asarray(adata.X, dtype=float)  # integer-valued floats -> valid current

    monkeypatch.setattr(adata.raw, "to_adata", lambda: (_ for _ in ()).throw(RuntimeError("raw bad")))
    out = au.get_raw_data_source(adata, prefer_complete_genes=True, require_integer_counts=True)
    assert out.source == "current"


def test_ensure_counts_layer_rejects_normalized_raw_falls_back_to_x(
    minimal_spatial_adata,
):
    """When .raw has normalized data but X has integer counts, use X fallback."""
    adata = minimal_spatial_adata.copy()
    # Create a .raw with normalized (non-integer) data
    raw_norm = adata.copy()
    raw_norm.X = np.asarray(raw_norm.X, dtype=float) + 0.5  # non-integer
    adata.raw = raw_norm
    # Remove any existing counts layer
    if "counts" in adata.layers:
        del adata.layers["counts"]

    # X still has integer counts (Poisson) → fallback should succeed
    created = au.ensure_counts_layer(adata)
    assert created is True
    assert "counts" in adata.layers
    np.testing.assert_array_equal(
        np.asarray(adata.layers["counts"]),
        np.asarray(adata.X),
    )


def test_ensure_counts_layer_accepts_integer_raw(minimal_spatial_adata):
    """ensure_counts_layer should accept .raw that contains valid integer counts."""
    adata = minimal_spatial_adata.copy()
    raw_int = adata.copy()
    raw_int.X = np.abs(np.round(np.asarray(raw_int.X))).astype(float)
    adata.raw = raw_int
    if "counts" in adata.layers:
        del adata.layers["counts"]

    created = au.ensure_counts_layer(adata)
    assert created is True
    assert "counts" in adata.layers


def test_ensure_counts_layer_raises_when_no_valid_source(
    minimal_spatial_adata,
):
    """ensure_counts_layer raises when no raw, no counts, and X is normalized."""
    adata = minimal_spatial_adata.copy()
    adata.raw = None
    adata.layers.clear()
    adata.X = np.asarray(adata.X, dtype=float) + 0.5  # normalized
    with pytest.raises(DataNotFoundError):
        au.ensure_counts_layer(adata)


def test_ensure_counts_layer_falls_back_to_x_when_raw_absent(
    minimal_spatial_adata,
):
    """ensure_counts_layer should use X when .raw is None and X is integer counts."""
    adata = minimal_spatial_adata.copy()
    adata.raw = None
    if "counts" in adata.layers:
        del adata.layers["counts"]
    # Ensure X is integer
    adata.X = np.abs(np.round(np.asarray(adata.X))).astype(float)

    created = au.ensure_counts_layer(adata)
    assert created is True
    assert "counts" in adata.layers
    np.testing.assert_array_equal(
        np.asarray(adata.layers["counts"]),
        np.asarray(adata.X),
    )


def test_ensure_counts_layer_raises_when_raw_and_x_both_normalized(
    minimal_spatial_adata,
):
    """ensure_counts_layer raises when both .raw and X contain normalized data."""
    adata = minimal_spatial_adata.copy()
    raw_norm = adata.copy()
    raw_norm.X = np.asarray(raw_norm.X, dtype=float) + 0.5
    adata.raw = raw_norm
    adata.X = np.asarray(adata.X, dtype=float) + 0.5  # also normalized
    if "counts" in adata.layers:
        del adata.layers["counts"]
    with pytest.raises(DataNotFoundError):
        au.ensure_counts_layer(adata)
