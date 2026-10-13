"""Numerical mitochondrial QC and shared gene-selection regression tests."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from chatspatial.models.data import (
    PreprocessingParameters,
    SpatialVariableGenesParameters,
)
from chatspatial.tools.preprocessing import (
    _calculate_qc_metrics,
    _filter_and_subsample,
    _select_and_subsample_genes,
)
from chatspatial.tools.spatial_genes import _select_testable_genes
from chatspatial.utils.exceptions import DataError
from chatspatial.utils.gene_annotations import mitochondrial_mask


class Context:
    def __init__(self):
        self.warnings = []

    async def warning(self, message):
        self.warnings.append(message)


def make_data(symbol, sparse_matrix=False):
    counts = np.array([[80.0, 20.0], [5.0, 95.0], [0.0, 100.0], [0.0, 0.0]])
    adata = AnnData(sparse.csr_matrix(counts) if sparse_matrix else counts)
    adata.var_names = [symbol, "GENE1"]
    return adata


@pytest.mark.parametrize("symbol", ["MT-CO2", "MT.CO2", "mt-Co2", "mt.Co2"])
@pytest.mark.parametrize("sparse_matrix", [False, True])
async def test_real_qc_and_filtering(symbol, sparse_matrix):
    adata = make_data(symbol, sparse_matrix)
    column, metrics = _calculate_qc_metrics(adata)
    np.testing.assert_allclose(adata.obs[column], [80, 5, 0, np.nan], equal_nan=True)
    assert metrics["median_mito_pct"] == 5
    assert metrics["max_mito_pct"] == 80
    ctx = Context()
    result = await _filter_and_subsample(
        adata,
        PreprocessingParameters(filter_genes_min_cells=None, filter_cells_min_genes=1),
        ctx,
        metrics,
        column,
    )
    assert list(result.obs_names) == ["1", "2"]
    assert metrics["mito_filter_applied"] is True
    assert not ctx.warnings


async def test_unknown_gene_set_is_unavailable_and_filtering_is_skipped():
    adata = make_data("ENSG00000198712")
    column, metrics = _calculate_qc_metrics(adata)
    assert column is None
    assert metrics["n_mt_genes"] == 0
    assert metrics["mito_qc_available"] is False
    assert "median_mito_pct" not in metrics
    assert adata.obs.pct_counts_mt.isna().all()
    ctx = Context()
    result = await _filter_and_subsample(
        adata,
        PreprocessingParameters(filter_genes_min_cells=None, filter_cells_min_genes=1),
        ctx,
        metrics,
        column,
    )
    assert result.n_obs == 3
    assert metrics["mito_filter_applied"] is False
    assert any("Filtering was skipped" in warning for warning in ctx.warnings)


async def test_annotation_preserved_and_shared_with_hvg_and_svg():
    adata = make_data("ENSG00000198712")
    adata.var["mt"] = [True, False]
    column, metrics = _calculate_qc_metrics(adata)
    assert metrics["mito_annotation_source"] == "var.mt"
    np.testing.assert_allclose(adata.obs[column].iloc[:3], [80, 5, 0])
    params = PreprocessingParameters(normalization="none", n_hvgs=2)
    selected = await _select_and_subsample_genes(adata, params, Context())
    assert selected.var.highly_variable.tolist() == [False, True]
    mask = _select_testable_genes(
        selected,
        list(selected.var_names),
        SpatialVariableGenesParameters(method="flashs", test_only_hvg=False),
    )
    assert mask.tolist() == [False, True]


def test_raw_superset_and_reordered_genes_use_identifiers():
    adata = AnnData(np.ones((2, 3)))
    adata.var_names = ["ENSG_MT", "MT.CO2", "GENE1"]
    adata.var["mt"] = [True, False, False]
    adata.raw = adata.copy()
    adata = adata[:, ["GENE1", "MT.CO2"]].copy()
    names = ["GENE1", "ENSG_MT", "MT.CO2", "mt.Nd1", "MTCO2P12"]
    mask, _ = mitochondrial_mask(adata, names)
    assert mask.tolist() == [False, True, False, True, False]
    mask = _select_testable_genes(
        adata,
        names,
        SpatialVariableGenesParameters(method="flashs", test_only_hvg=False),
    )
    assert mask.tolist() == [True, False, True, False, True]


@pytest.mark.parametrize(
    "annotation", [["True", "False"], pd.array([True, None], dtype="boolean")]
)
def test_invalid_annotation_does_not_silently_cast_to_boolean(annotation):
    adata = make_data("MT.CO2")
    adata.var["mt"] = annotation
    with pytest.raises(DataError, match="boolean values"):
        mitochondrial_mask(adata)
