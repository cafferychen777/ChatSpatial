"""SpaGCN keeps spots that have no expression in the genes it receives.

identify_spatial_domains passes SpaGCN the highly variable genes. After UMI
downsampling, a spot can have counts only outside that gene set, so its row is
all zero. SpaGCN builds its spatial graph from every spot, and the domain labels
must still cover every spot. This test trains SpaGCN (torch), so it is marked
slow and runs where spagcn-modern is installed.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest
import scanpy as sc

from chatspatial.models.data import SpatialDomainParameters
from chatspatial.tools import spatial_domains as domains_module
from chatspatial.tools.spatial_domains import identify_spatial_domains


class DummyCtx:
    def __init__(self, adata):
        self._adata = adata
        self.warnings: list[str] = []

    async def get_adata(self, data_id: str):
        return self._adata

    async def set_adata(self, data_id: str, adata):
        self._adata = adata

    async def warning(self, msg: str):
        self.warnings.append(msg)

    async def info(self, msg: str):
        return None


def _banded_visium_like_adata(zero_spot: int) -> ad.AnnData:
    """Three horizontal bands with band-specific genes on a 20 x 15 grid."""
    rng = np.random.default_rng(7)
    rows, cols = np.meshgrid(np.arange(20), np.arange(15), indexing="ij")
    rows, cols = rows.ravel(), cols.ravel()
    band = np.minimum(rows // 7, 2)
    n_hvg, n_other = 60, 20
    rates = np.full((3, n_hvg + n_other), 1.0)
    for b in range(3):
        rates[b, b * 20 : (b + 1) * 20] = 3.0
    counts = rng.poisson(rates[band]).astype(np.float64)
    counts[zero_spot, :n_hvg] = 0.0  # counts only outside the selected genes
    counts[zero_spot, n_hvg] = 2.0

    adata = ad.AnnData(counts)
    adata.obs_names = [f"spot_{i}" for i in range(adata.n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(adata.n_vars)]
    adata.obsm["spatial"] = np.column_stack([cols * 100.0, rows * 100.0])
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.var["highly_variable"] = np.arange(adata.n_vars) < n_hvg
    return adata


@pytest.mark.slow
@pytest.mark.integration
@pytest.mark.asyncio
async def test_spagcn_labels_every_spot_when_one_has_no_selected_gene_counts(
    monkeypatch: pytest.MonkeyPatch,
):
    pytest.importorskip("torch")
    pytest.importorskip("SpaGCN")
    monkeypatch.setattr(domains_module, "export_analysis_result", lambda *a, **k: None)

    zero_spot = 37
    adata = _banded_visium_like_adata(zero_spot)
    hvg = adata.var["highly_variable"].to_numpy()
    assert adata[zero_spot, hvg].X.sum() == 0  # the condition that used to fail

    ctx = DummyCtx(adata)
    params = SpatialDomainParameters(
        method="spagcn", n_domains=3, spagcn_use_histology=False
    )
    result = await identify_spatial_domains("d", ctx, params)

    labels = ctx._adata.obs[result.domain_key]
    assert len(labels) == adata.n_obs
    assert labels.notna().all()
    assert sum(result.domain_counts.values()) == adata.n_obs
