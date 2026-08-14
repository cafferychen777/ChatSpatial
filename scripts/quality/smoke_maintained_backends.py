"""Exercise maintained optional backends from an installed ChatSpatial wheel."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy.sparse import csr_matrix


def _adata(n_obs: int = 24, n_vars: int = 12) -> ad.AnnData:
    rng = np.random.default_rng(2026)
    result = ad.AnnData(rng.poisson(3.0, size=(n_obs, n_vars)).astype(np.float32))
    result.obs_names = [f"spot_{index}" for index in range(n_obs)]
    result.var_names = [f"gene_{index}" for index in range(n_vars)]
    result.var["highly_variable"] = True
    result.obsm["spatial"] = np.column_stack(
        (np.arange(n_obs) % 6, np.arange(n_obs) // 6)
    ).astype(float)
    return result


def smoke_spatialde() -> None:
    import SpatialDE

    rng = np.random.default_rng(7)
    means = np.linspace(1, 20, 12)
    counts = pd.DataFrame(
        [rng.negative_binomial(2.0, 2.0 / (2.0 + mean), size=40) for mean in means],
        index=[f"gene_{index}" for index in range(12)],
    )
    stabilized = SpatialDE.stabilize(counts)
    sample_info = pd.DataFrame(
        {"total_counts": counts.sum(axis=0)}, index=stabilized.columns
    )
    residuals = SpatialDE.regress_out(sample_info, stabilized, "np.log(total_counts)")
    assert residuals.shape == counts.shape
    assert np.isfinite(residuals.to_numpy()).all()


def smoke_graphst() -> None:
    from GraphST.GraphST import GraphST

    result = GraphST(_adata(), epochs=2, dim_output=4, random_seed=7).train()
    assert result.obsm["emb"].shape == result.shape
    assert np.isfinite(result.obsm["emb"]).all()


def smoke_stagate() -> None:
    import STAGATE_pyG

    data = _adata()
    STAGATE_pyG.Cal_Spatial_Net(data, rad_cutoff=1.1, verbose=False)
    STAGATE_pyG.train_STAGATE(
        data,
        hidden_dims=[8, 4],
        n_epochs=2,
        random_seed=7,
        verbose=False,
        device="cpu",
    )
    assert data.obsm["STAGATE"].shape == (data.n_obs, 4)
    assert np.isfinite(data.obsm["STAGATE"]).all()


def smoke_stalign() -> None:
    from STalign import STalign

    grid = [torch.linspace(0, 7, 8), torch.linspace(0, 7, 8)]
    rows, columns = np.mgrid[:8, :8]
    image = np.exp(-((columns - 3) ** 2 + (rows - 3) ** 2) / 4)[None].astype(np.float32)
    result = STalign.LDDMM(
        grid,
        image,
        grid,
        image,
        a=3,
        nt=2,
        niter=1,
        epV=1,
        device="cpu",
        dtype=torch.float32,
    )
    assert result["A"].shape == (3, 3)
    assert torch.isfinite(result["v"]).all()


def smoke_paste() -> None:
    import paste

    first = _adata(n_obs=12, n_vars=8)
    second = first.copy()
    second.obsm["spatial"] = second.obsm["spatial"] + 0.25
    plan = paste.pairwise_align(
        first,
        second,
        numItermax=20,
        gpu_verbose=False,
    )
    stacked = paste.stack_slices_pairwise([first, second], [plan])
    assert plan.shape == (12, 12)
    assert np.isfinite(plan).all()
    assert len(stacked) == 2


def smoke_cellrank() -> None:
    import cellrank as cr

    rng = np.random.default_rng(9)
    transition = rng.random((100, 100))
    transition /= transition.sum(axis=1, keepdims=True)
    data = ad.AnnData(rng.normal(size=(100, 5)))
    kernel = cr.kernels.PrecomputedKernel(csr_matrix(transition), data)
    kernel.compute_transition_matrix()
    estimator = cr.estimators.GPCCA(kernel)
    estimator.compute_eigendecomposition()
    estimator.compute_macrostates(n_states=2)
    assert estimator.macrostates_memberships.shape == (100, 2)
    assert estimator.coarse_T.shape == (2, 2)


def main() -> None:
    smoke_spatialde()
    smoke_graphst()
    smoke_stagate()
    smoke_stalign()
    smoke_paste()
    smoke_cellrank()
    print("Maintained optional backend smoke tests passed")


if __name__ == "__main__":
    main()
