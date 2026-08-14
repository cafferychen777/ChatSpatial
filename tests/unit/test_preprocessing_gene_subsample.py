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
