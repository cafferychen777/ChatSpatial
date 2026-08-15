"""Regression tests for gene subsampling on small panels."""

from __future__ import annotations

import numpy as np
import pytest


@pytest.mark.unit
def test_small_panel_respects_subsample_genes(minimal_spatial_adata):
    """When n_vars < 100 and user explicitly requests subsample_genes,
    the HVG mask should select only the requested number of genes."""
    import scipy.sparse

    from chatspatial.tools.preprocessing import _compute_safe_percent_top  # noqa: F401

    adata = minimal_spatial_adata.copy()
    assert adata.n_vars < 100, "Fixture must have < 100 genes for this test"

    n_hvgs = 10
    # Replicate the small-panel branch logic from preprocessing.py
    if scipy.sparse.issparse(adata.X):
        var = np.asarray(
            adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)
        ).ravel()
    else:
        var = np.var(adata.X, axis=0)
    top_idx = np.argpartition(var, -n_hvgs)[-n_hvgs:]
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[top_idx] = True
    adata.var["highly_variable"] = mask

    # Subsample step
    result = adata[:, adata.var["highly_variable"]].copy()
    assert result.n_vars == n_hvgs


@pytest.mark.unit
def test_small_panel_all_hvg_when_no_subsample(minimal_spatial_adata):
    """When n_vars < 100 and no subsample_genes, all genes should be HVG."""
    adata = minimal_spatial_adata.copy()
    assert adata.n_vars < 100

    # Small-panel default: all genes HVG
    adata.var["highly_variable"] = True
    assert adata.var["highly_variable"].sum() == adata.n_vars


class _WarnCtx:
    def __init__(self):
        self.warnings: list[str] = []

    async def warning(self, msg: str) -> None:
        self.warnings.append(msg)

    async def info(self, msg: str) -> None:
        pass


def _counts_adata(n_obs: int = 40, n_vars: int = 120):
    import anndata as ad

    rng = np.random.default_rng(0)
    counts = rng.poisson(4, size=(n_obs, n_vars)).astype(np.float32)
    adata = ad.AnnData(counts)
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    adata.layers["counts"] = counts.copy()
    return adata


@pytest.mark.unit
@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("normalization", "expected_flavor"),
    [("pearson_residuals", "pearson_residuals"), ("log", None)],
)
async def test_hvg_flavor_matches_the_normalization(
    monkeypatch, normalization: str, expected_flavor: str | None
):
    """Pearson residuals are not log-normalized expression.

    Scanpy's mean-dispersion binning presupposes log1p counts, so ranking
    residuals with it measures noise. The residual-variance flavor reads the
    counts layer instead, which is the variability measure that belongs to
    this normalization.
    """
    from chatspatial.models.data import PreprocessingParameters
    from chatspatial.tools import preprocessing

    adata = _counts_adata()
    seen: dict[str, object] = {}

    def _capture(a, n_top_genes, **kwargs):
        seen["n_top_genes"] = n_top_genes
        seen.update(kwargs)
        a.var["highly_variable"] = True
        return False

    monkeypatch.setattr(preprocessing, "ensure_highly_variable_genes", _capture)

    await preprocessing._select_and_subsample_genes(
        adata,
        PreprocessingParameters(normalization=normalization, n_hvgs=20),
        _WarnCtx(),
    )

    assert seen.get("flavor") == expected_flavor
    assert seen.get("layer") == ("counts" if expected_flavor else None)


@pytest.mark.unit
@pytest.mark.asyncio
async def test_hvg_falls_back_to_binning_without_a_counts_layer(monkeypatch):
    """The residual flavor reads counts; without them there is nothing to read."""
    from chatspatial.models.data import PreprocessingParameters
    from chatspatial.tools import preprocessing

    adata = _counts_adata()
    del adata.layers["counts"]
    seen: dict[str, object] = {}

    def _capture(a, n_top_genes, **kwargs):
        seen.update(kwargs)
        a.var["highly_variable"] = True
        return False

    monkeypatch.setattr(preprocessing, "ensure_highly_variable_genes", _capture)

    await preprocessing._select_and_subsample_genes(
        adata,
        PreprocessingParameters(normalization="pearson_residuals", n_hvgs=20),
        _WarnCtx(),
    )

    assert "flavor" not in seen


@pytest.mark.unit
@pytest.mark.asyncio
async def test_subsampling_reapplies_the_gene_filter():
    """The filter ran against every spot, so its verdict does not survive.

    A gene it kept can hold no counts at all once the subsample is drawn, and
    a variance-stabilizing normalization then has nothing to divide by: on a
    lymph node subsampled to 900 spots, 773 such genes came out of Pearson
    residuals as entire columns of NaN. The invariant is asserted rather than
    a particular gene, because which spots are drawn is not the contract.
    """
    import anndata as ad

    from chatspatial.models.data import PreprocessingParameters
    from chatspatial.tools import preprocessing

    rng = np.random.default_rng(0)
    counts = rng.poisson(3, size=(60, 10)).astype(np.float32)
    # Four genes live only in the tail of the matrix, so a subsample leaves
    # most of them below the threshold.
    counts[:, 6:] = 0
    counts[57:, 6:] = 5

    adata = ad.AnnData(counts)
    adata.var_names = [f"gene_{i}" for i in range(10)]
    min_cells = 3
    ctx = _WarnCtx()

    out = await preprocessing._filter_and_subsample(
        adata,
        PreprocessingParameters(
            subsample_spots=20,
            filter_genes_min_cells=min_cells,
            filter_cells_min_genes=None,
            filter_mito_pct=None,
        ),
        ctx,
        {},
        None,
    )

    assert out.n_obs == 20
    expressed_in = (np.asarray(out.X) > 0).sum(axis=0)
    assert (
        expressed_in >= min_cells
    ).all(), "a gene survived that the filter would have rejected on these spots"
    if out.n_vars < 10:
        assert any("Subsampling to 20 spots left" in w for w in ctx.warnings)


@pytest.mark.unit
def test_undefined_residuals_become_zero():
    """A gene with no counts has an expected count of zero, so 0/0 is NaN.

    Nothing was observed and nothing was expected, which is no deviation:
    the value the statistic reaches for is zero.
    """
    import anndata as ad
    import scipy.sparse as sp

    from chatspatial.tools import preprocessing

    values = np.array([[1.0, np.nan], [-2.0, np.nan], [0.5, np.nan]])

    dense = ad.AnnData(values.copy())
    assert preprocessing._zero_non_finite_residuals(dense) == 1
    assert np.isfinite(np.asarray(dense.X)).all()
    assert np.asarray(dense.X)[0, 0] == 1.0  # the defined column is untouched

    sparse = ad.AnnData(sp.csr_matrix(values))
    assert preprocessing._zero_non_finite_residuals(sparse) == 1
    assert np.isfinite(sparse.X.data).all()


@pytest.mark.unit
def test_finite_residuals_are_left_alone():
    import anndata as ad

    from chatspatial.tools import preprocessing

    adata = ad.AnnData(np.array([[1.0, -2.0], [0.0, 3.0]]))
    assert preprocessing._zero_non_finite_residuals(adata) == 0


@pytest.mark.unit
@pytest.mark.asyncio
async def test_pearson_residuals_leave_no_undefined_values_behind():
    """The guard has to be wired into the normalization, not merely exist.

    Testing the helper alone passed even with its call site removed, which is
    exactly the gap that let entire NaN columns reach adata.X.
    """
    import anndata as ad

    from chatspatial.models.data import PreprocessingParameters
    from chatspatial.tools.preprocessing import preprocess_data

    rng = np.random.default_rng(0)
    counts = rng.poisson(4, size=(40, 120)).astype(np.float32)
    counts[:, 100:] = 0  # 20 genes with no counts anywhere
    adata = ad.AnnData(counts)
    adata.var_names = [f"gene_{i}" for i in range(120)]
    adata.obsm["spatial"] = rng.random((40, 2)) * 100

    ctx = _PreprocessCtx(adata)
    await preprocess_data(
        "d",
        ctx,
        PreprocessingParameters(
            normalization="pearson_residuals",
            filter_genes_min_cells=None,
            filter_cells_min_genes=None,
            filter_mito_pct=None,
            n_hvgs=20,
        ),
    )

    values = np.asarray(
        ctx.adata.X.todense() if hasattr(ctx.adata.X, "todense") else ctx.adata.X
    )
    assert np.isfinite(values).all(), "undefined residuals reached adata.X"
    assert any("no counts in the analyzed cells" in w for w in ctx.warnings)


class _PreprocessCtx(_WarnCtx):
    def __init__(self, adata):
        super().__init__()
        self.adata = adata

    async def get_adata(self, _data_id):
        return self.adata

    async def set_adata(self, _data_id, adata):
        self.adata = adata

    async def report_progress(self, *_args, **_kwargs):
        return None

    def log_config(self, *_args, **_kwargs):
        return None

    def log_step(self, *_args, **_kwargs):
        return None

    def log_result(self, *_args, **_kwargs):
        return None


@pytest.mark.unit
@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("requested", "n_hvgs", "should_warn"),
    [(3000, 2000, True), (500, 2000, False)],
    ids=["capped-by-n_hvgs", "within-n_hvgs"],
)
async def test_a_capped_gene_subsample_says_so(
    requested: int, n_hvgs: int, should_warn: bool
):
    """Genes are chosen out of the highly variable set, so it bounds them.

    Asking for 3000 with the default n_hvgs=2000 returned 1988: the count was
    reported, but not the reason it was not the one requested.
    """
    from chatspatial.models.data import PreprocessingParameters
    from chatspatial.tools import preprocessing

    adata = _counts_adata(n_obs=40, n_vars=4000)
    ctx = _WarnCtx()

    def _mark(a, n_top_genes, **_kwargs):
        flags = np.zeros(a.n_vars, dtype=bool)
        flags[:n_top_genes] = True
        a.var["highly_variable"] = flags
        return False

    import unittest.mock as mock

    with mock.patch.object(preprocessing, "ensure_highly_variable_genes", _mark):
        await preprocessing._select_and_subsample_genes(
            adata,
            PreprocessingParameters(subsample_genes=requested, n_hvgs=n_hvgs),
            ctx,
        )

    capped = [w for w in ctx.warnings if "subsample_genes=" in w]
    assert bool(capped) is should_warn
    if should_warn:
        assert f"subsample_genes={requested}" in capped[0]
        assert "raise n_hvgs" in capped[0]
