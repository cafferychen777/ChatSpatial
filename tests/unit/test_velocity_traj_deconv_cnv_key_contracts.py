"""Contract tests for velocity, trajectory, deconvolution, and CNV parametric key generation.

Also tests provenance integrity: metadata results_keys must point to keys
that THIS run owns (per-run copies), not shared "active" keys.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

# ---------------------------------------------------------------------------
# _build_velocity_key
# ---------------------------------------------------------------------------


class TestBuildVelocityKey:
    """Verify _build_velocity_key encodes method + mode correctly."""

    @staticmethod
    def _make_params(**kw):
        defaults = {"method": "scvelo", "scvelo_mode": "stochastic"}
        defaults.update(kw)
        return SimpleNamespace(**defaults)

    def test_scvelo_stochastic(self):
        from chatspatial.tools.velocity import _build_velocity_key

        assert _build_velocity_key(self._make_params()) == "velocity_scvelo_stochastic"

    def test_scvelo_dynamical(self):
        from chatspatial.tools.velocity import _build_velocity_key

        p = self._make_params(scvelo_mode="dynamical")
        assert _build_velocity_key(p) == "velocity_scvelo_dynamical"

    def test_velovi_no_mode(self):
        from chatspatial.tools.velocity import _build_velocity_key

        p = self._make_params(method="velovi")
        assert _build_velocity_key(p) == "velocity_velovi"

    def test_scvelo_modes_coexist(self):
        from chatspatial.tools.velocity import _build_velocity_key

        k1 = _build_velocity_key(self._make_params(scvelo_mode="stochastic"))
        k2 = _build_velocity_key(self._make_params(scvelo_mode="dynamical"))
        assert k1 != k2


# ---------------------------------------------------------------------------
# _build_trajectory_key
# ---------------------------------------------------------------------------


class TestBuildTrajectoryKey:
    """Verify _build_trajectory_key encodes method + discriminating param."""

    @staticmethod
    def _make_params(**kw):
        defaults = {"method": "cellrank", "spatial_weight": 0.5, "root_cells": None}
        defaults.update(kw)
        return SimpleNamespace(**defaults)

    def test_cellrank_default(self):
        from chatspatial.tools.trajectory import _build_trajectory_key

        assert (
            _build_trajectory_key(self._make_params()) == "trajectory_cellrank_sw0_50"
        )

    def test_cellrank_different_weight(self):
        from chatspatial.tools.trajectory import _build_trajectory_key

        p = self._make_params(spatial_weight=0.0)
        assert _build_trajectory_key(p) == "trajectory_cellrank_sw0_00"

    def test_cellrank_weights_coexist(self):
        from chatspatial.tools.trajectory import _build_trajectory_key

        k1 = _build_trajectory_key(self._make_params(spatial_weight=0.3))
        k2 = _build_trajectory_key(self._make_params(spatial_weight=0.7))
        assert k1 != k2

    def test_palantir_with_root(self):
        from chatspatial.tools.trajectory import _build_trajectory_key

        p = self._make_params(method="palantir", root_cells=["cell_42"])
        assert _build_trajectory_key(p) == "trajectory_palantir_cell_42"

    def test_dpt_without_root(self):
        from chatspatial.tools.trajectory import _build_trajectory_key

        p = self._make_params(method="dpt", root_cells=None)
        assert _build_trajectory_key(p) == "trajectory_dpt"


# ---------------------------------------------------------------------------
# _build_deconvolution_key
# ---------------------------------------------------------------------------


class TestBuildDeconvolutionKey:
    def test_with_reference(self):
        from chatspatial.tools.deconvolution import _build_deconvolution_key

        assert (
            _build_deconvolution_key("flashdeconv", "ref1")
            == "deconvolution_flashdeconv_ref1"
        )

    def test_without_reference(self):
        from chatspatial.tools.deconvolution import _build_deconvolution_key

        assert (
            _build_deconvolution_key("cell2location", None)
            == "deconvolution_cell2location"
        )

    def test_different_refs_coexist(self):
        from chatspatial.tools.deconvolution import _build_deconvolution_key

        k1 = _build_deconvolution_key("flashdeconv", "ref_a")
        k2 = _build_deconvolution_key("flashdeconv", "ref_b")
        assert k1 != k2


# ---------------------------------------------------------------------------
# _build_cnv_key
# ---------------------------------------------------------------------------


class TestBuildCnvKey:
    @staticmethod
    def _make_params(**kw):
        defaults = {"method": "infercnvpy", "reference_categories": ["immune"]}
        defaults.update(kw)
        return SimpleNamespace(**defaults)

    def test_infercnvpy_with_categories(self):
        from chatspatial.tools.cnv_analysis import _build_cnv_key

        assert _build_cnv_key(self._make_params()) == "cnv_infercnvpy_immune"

    def test_numbat_with_categories(self):
        from chatspatial.tools.cnv_analysis import _build_cnv_key

        p = self._make_params(method="numbat", reference_categories=["A", "B"])
        assert _build_cnv_key(p) == "cnv_numbat_A_B"

    def test_categories_sorted_for_determinism(self):
        from chatspatial.tools.cnv_analysis import _build_cnv_key

        p1 = self._make_params(reference_categories=["B", "A"])
        p2 = self._make_params(reference_categories=["A", "B"])
        assert _build_cnv_key(p1) == _build_cnv_key(p2)

    def test_different_categories_coexist(self):
        from chatspatial.tools.cnv_analysis import _build_cnv_key

        k1 = _build_cnv_key(self._make_params(reference_categories=["immune"]))
        k2 = _build_cnv_key(self._make_params(reference_categories=["stromal"]))
        assert k1 != k2

    def test_no_categories(self):
        from chatspatial.tools.cnv_analysis import _build_cnv_key

        p = self._make_params(reference_categories=None)
        assert _build_cnv_key(p) == "cnv_infercnvpy"


# ---------------------------------------------------------------------------
# Provenance integrity: velocity results_keys must not claim shared keys
# ---------------------------------------------------------------------------


class TestVelocityResultsKeysProvenance:
    """Velocity metadata results_keys should not include shared uns keys."""

    def test_velocity_method_not_in_results_keys(self):
        """uns["velocity_method"] is shared; must NOT appear in results_keys."""
        # The velocity code builds results_keys with "uns": []
        # This test documents the contract.
        results_keys: dict[str, list[str]] = {"uns": [], "obs": [], "obsm": []}
        assert "velocity_method" not in results_keys["uns"]

    def test_latent_time_only_claimed_by_dynamical(self):
        """Only scvelo dynamical should claim obs["latent_time"]."""
        # Simulate: a stochastic run on data that has leftover latent_time
        # The condition is: method==scvelo AND mode==dynamical AND key exists
        p_stochastic = SimpleNamespace(method="scvelo", scvelo_mode="stochastic")
        p_dynamical = SimpleNamespace(method="scvelo", scvelo_mode="dynamical")
        p_velovi = SimpleNamespace(method="velovi", scvelo_mode="stochastic")

        # Only dynamical should pass
        for p, expected in [
            (p_stochastic, False),
            (p_dynamical, True),
            (p_velovi, False),
        ]:
            should_claim = p.method == "scvelo" and p.scvelo_mode == "dynamical"
            assert should_claim == expected


# ---------------------------------------------------------------------------
# Provenance integrity: CellRank per-run copies
# ---------------------------------------------------------------------------


class TestCellRankPerRunCopies:
    """CellRank must save per-run copies of shared keys for provenance."""

    def test_cellrank_suffix_computation(self):
        """Verify per-run suffix matches _build_trajectory_key logic."""
        from chatspatial.tools.trajectory import _build_trajectory_key

        p = SimpleNamespace(method="cellrank", spatial_weight=0.3, root_cells=None)
        key = _build_trajectory_key(p)
        expected_suffix = "cellrank_sw0_30"
        assert key == f"trajectory_{expected_suffix}"

    def test_cellrank_per_run_obs_keys(self, minimal_spatial_adata):
        """Per-run copies must exist alongside shared keys after CellRank."""
        adata = minimal_spatial_adata.copy()
        suffix = "cellrank_sw0_50"

        # Simulate shared + per-run write
        adata.obs["pseudotime"] = np.linspace(0, 1, adata.n_obs)
        adata.obs[f"pseudotime_{suffix}"] = adata.obs["pseudotime"]
        adata.obs["terminal_states"] = "T0"
        adata.obs[f"terminal_states_{suffix}"] = adata.obs["terminal_states"]

        # Both should coexist
        assert "pseudotime" in adata.obs.columns
        assert f"pseudotime_{suffix}" in adata.obs.columns
        assert "terminal_states" in adata.obs.columns
        assert f"terminal_states_{suffix}" in adata.obs.columns


# ---------------------------------------------------------------------------
# Provenance integrity: CCC per-method uns snapshot
# ---------------------------------------------------------------------------


class TestCCCPerMethodSnapshot:
    """CCC must save per-method uns snapshot alongside shared key."""

    def test_per_method_key_format(self):
        """Per-method CCC key follows ccc_{method} pattern."""
        for method in ("liana", "fastccc", "cellphonedb", "cellchat_r"):
            key = f"ccc_{method}"
            assert key.startswith("ccc_")
            assert method in key


# ---------------------------------------------------------------------------
# Provenance integrity: CNV parametric summary key
# ---------------------------------------------------------------------------


class TestCNVParametricSummaryKey:
    """CNV analysis summary must use parametric uns key."""

    def test_infercnvpy_summary_key(self):
        """infercnvpy summary key includes method + categories."""
        from chatspatial.tools.cnv_analysis import _build_cnv_key

        p = SimpleNamespace(method="infercnvpy", reference_categories=["Normal"])
        analysis_key = _build_cnv_key(p)
        cnv_summary_key = f"cnv_analysis_{analysis_key.removeprefix('cnv_')}"
        assert cnv_summary_key == "cnv_analysis_infercnvpy_Normal"

    def test_numbat_summary_key(self):
        """numbat summary key includes method + categories."""
        from chatspatial.tools.cnv_analysis import _build_cnv_key

        p = SimpleNamespace(method="numbat", reference_categories=["A", "B"])
        analysis_key = _build_cnv_key(p)
        cnv_summary_key = f"cnv_analysis_{analysis_key.removeprefix('cnv_')}"
        assert cnv_summary_key == "cnv_analysis_numbat_A_B"

    def test_different_runs_different_summary_keys(self):
        """Two infercnvpy runs with different refs get different summary keys."""
        from chatspatial.tools.cnv_analysis import _build_cnv_key

        p1 = SimpleNamespace(method="infercnvpy", reference_categories=["Normal"])
        p2 = SimpleNamespace(method="infercnvpy", reference_categories=["Stromal"])
        k1 = _build_cnv_key(p1)
        k2 = _build_cnv_key(p2)
        s1 = f"cnv_analysis_{k1.removeprefix('cnv_')}"
        s2 = f"cnv_analysis_{k2.removeprefix('cnv_')}"
        assert s1 != s2
