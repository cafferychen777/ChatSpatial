"""Regression tests for trajectory auto-root selection."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from chatspatial.tools.trajectory import (
    compute_dpt_trajectory,
    infer_pseudotime_palantir,
)


class TestPalantirAutoRoot:
    """Issue 5a: Palantir auto root should use abs() for sign invariance."""

    def test_auto_root_uses_abs_of_first_dc(self, minimal_spatial_adata, caplog):
        """Auto-selected root cell should be based on abs(DC1),
        making it invariant to eigenvector sign flips."""
        adata = minimal_spatial_adata.copy()

        # Pre-compute PCA so Palantir doesn't fail
        import scanpy as sc

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata)

        # Mock palantir to capture which start_cell is passed
        captured = {}

        class FakeDmRes(dict):
            pass

        def fake_run_diffusion_maps(pca_df, n_components=10):
            # Create fake eigenvectors where the largest absolute value
            # is negative (cell_2 has -0.9, cell_0 has +0.5)
            n = len(pca_df)
            evecs = np.zeros((n, n_components))
            # Make cell at index 2 have the largest absolute value
            # but with negative sign
            evecs[2, 0] = -0.9
            evecs[0, 0] = 0.5
            evecs[1, 0] = 0.3
            res = FakeDmRes()
            res["EigenVectors"] = evecs
            return res

        class FakePalantirResult:
            def __init__(self, n):
                self.pseudotime = pd.Series(np.linspace(0, 1, n), index=adata.obs_names)
                self.branch_probs = pd.DataFrame(np.ones((n, 1)), index=adata.obs_names)

        def fake_run_palantir(ms_data, start_cell, num_waypoints=500):
            captured["start_cell"] = start_cell
            return FakePalantirResult(len(ms_data))

        import types

        fake_palantir = types.ModuleType("palantir")
        fake_palantir.utils = types.SimpleNamespace(
            run_diffusion_maps=fake_run_diffusion_maps,
            determine_multiscale_space=lambda result: pd.DataFrame(
                result["EigenVectors"], index=adata.obs_names
            ),
        )
        fake_palantir.core = types.SimpleNamespace(run_palantir=fake_run_palantir)
        with caplog.at_level(logging.WARNING):
            infer_pseudotime_palantir(
                adata,
                root_cells=None,
                palantir_module=fake_palantir,
            )

        # With abs(), cell_2 should be selected (|-0.9| > |0.5|)
        assert (
            captured["start_cell"] == "cell_2"
        ), f"Expected cell_2 (largest abs DC1), got {captured['start_cell']}"
        # Warning should be emitted
        assert any("No root cell specified" in msg for msg in caplog.messages)

    def test_explicit_root_skips_auto_selection(self, minimal_spatial_adata):
        """When root_cells is provided, auto-selection is skipped."""
        adata = minimal_spatial_adata.copy()

        import scanpy as sc

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata)

        captured = {}

        class FakeDmRes(dict):
            pass

        def fake_run_diffusion_maps(pca_df, n_components=10):
            n = len(pca_df)
            evecs = np.zeros((n, n_components))
            evecs[0, 0] = 1.0
            res = FakeDmRes()
            res["EigenVectors"] = evecs
            return res

        class FakePalantirResult:
            def __init__(self, n):
                self.pseudotime = pd.Series(np.linspace(0, 1, n), index=adata.obs_names)
                self.branch_probs = pd.DataFrame(np.ones((n, 1)), index=adata.obs_names)

        def fake_run_palantir(ms_data, start_cell, num_waypoints=500):
            captured["start_cell"] = start_cell
            return FakePalantirResult(len(ms_data))

        import types

        fake_palantir = types.ModuleType("palantir")
        fake_palantir.utils = types.SimpleNamespace(
            run_diffusion_maps=fake_run_diffusion_maps,
            determine_multiscale_space=lambda result: pd.DataFrame(
                result["EigenVectors"], index=adata.obs_names
            ),
        )
        fake_palantir.core = types.SimpleNamespace(run_palantir=fake_run_palantir)
        infer_pseudotime_palantir(
            adata,
            root_cells=["cell_5"],
            palantir_module=fake_palantir,
        )
        assert captured["start_cell"] == "cell_5"


class TestDPTAutoRoot:
    """Issue 5b: DPT auto root should use diffusion map, not index 0."""

    def test_auto_root_uses_diffmap(self, minimal_spatial_adata, caplog):
        """DPT auto-root should pick cell with max abs(DC1)
        from X_diffmap, not just index 0."""
        adata = minimal_spatial_adata.copy()

        import scanpy as sc

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.diffmap(adata)

        # Manually set diffmap so cell at index 3 has largest abs value
        n_comps = adata.obsm["X_diffmap"].shape[1]
        adata.obsm["X_diffmap"] = np.zeros((adata.n_obs, n_comps), dtype=np.float64)
        adata.obsm["X_diffmap"][3, 0] = -2.5  # largest absolute
        adata.obsm["X_diffmap"][0, 0] = 1.0  # smaller absolute

        with caplog.at_level(logging.WARNING):
            compute_dpt_trajectory(adata, root_cells=None)

        assert (
            adata.uns["iroot"] == 3
        ), f"Expected iroot=3 (largest abs DC1), got {adata.uns['iroot']}"
        assert any("No root cell specified" in msg for msg in caplog.messages)

    def test_auto_root_fallback_without_diffmap(
        self, minimal_spatial_adata, monkeypatch, caplog
    ):
        """If X_diffmap is somehow missing, fall back to index 0."""
        adata = minimal_spatial_adata.copy()

        import scanpy as sc

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.diffmap(adata)

        # Remove diffmap to test fallback (but it will be recomputed
        # by ensure_diffmap; so we monkeypatch ensure_diffmap to not
        # recompute, simulating the edge case)
        from chatspatial.tools import trajectory as traj_mod

        monkeypatch.setattr(traj_mod, "ensure_diffmap", lambda adata: None)
        if "X_diffmap" in adata.obsm:
            del adata.obsm["X_diffmap"]

        with caplog.at_level(logging.WARNING):
            compute_dpt_trajectory(adata, root_cells=None)

        assert adata.uns["iroot"] == 0
        assert any("No root cell specified" in msg for msg in caplog.messages)
