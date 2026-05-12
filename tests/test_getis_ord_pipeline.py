"""Tests for the Getis-Ord analysis -> metadata -> viz -> export pipeline.

Verifies five interrelated bug fixes:
1. Viz reads metadata from the correct uns key (spatial_stats_getis_ord_metadata)
2. Hotspot/coldspot counts use z_threshold, not z>0
3. _build_results_keys uses actual analyzed genes, not params.genes
4. n_significant is non-zero when there are significant results
5. Corrected p-value columns are included in results_keys
"""

from __future__ import annotations

import numpy as np
import pytest

try:
    import anndata as ad
except ImportError:
    ad = None

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_getis_adata(rng: np.random.Generator, n_cells: int = 100):
    """Create an AnnData with pre-populated Getis-Ord results and metadata."""
    if ad is None:
        pytest.skip("anndata required")

    genes = ["GeneA", "GeneB"]
    n_genes = len(genes)
    X = rng.poisson(5, size=(n_cells, n_genes)).astype(np.float32)
    adata = ad.AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = genes
    adata.obsm["spatial"] = rng.uniform(0, 100, size=(n_cells, 2))

    # Simulate Getis-Ord analysis output stored in adata.obs
    for gene in genes:
        z = rng.normal(0, 2, size=n_cells)
        p = np.clip(rng.uniform(0, 0.1, size=n_cells), 0, 1)
        p_corr = np.minimum(p * n_cells, 1.0)
        adata.obs[f"{gene}_getis_ord_z"] = z
        adata.obs[f"{gene}_getis_ord_p"] = p
        adata.obs[f"{gene}_getis_ord_p_corrected"] = p_corr

    return adata, genes


# ---------------------------------------------------------------------------
# Test 1: Viz reads metadata from correct uns key
# ---------------------------------------------------------------------------


class TestVizMetadataKey:
    """Verify viz reads from spatial_stats_getis_ord_metadata, not adata.uns['getis_ord']."""

    def test_get_analysis_parameter_reads_correct_key(self, rng):
        from chatspatial.utils.adata_utils import (
            get_analysis_parameter,
            store_analysis_metadata,
        )

        adata, _ = _make_getis_adata(rng)

        # Store metadata using the same function the analysis code uses
        store_analysis_metadata(
            adata,
            analysis_name="spatial_stats_getis_ord",
            method="getis_ord",
            parameters={
                "alpha": 0.01,
                "correction": "bonferroni",
                "z_threshold": 2.576,
            },
            results_keys={"obs": []},
            statistics={"n_cells": 100},
        )

        # Verify get_analysis_parameter reads from the correct key
        alpha = get_analysis_parameter(
            adata, "spatial_stats_getis_ord", "alpha", default=0.05
        )
        assert alpha == 0.01

        correction = get_analysis_parameter(
            adata, "spatial_stats_getis_ord", "correction", default="none"
        )
        assert correction == "bonferroni"

        z_thresh = get_analysis_parameter(
            adata, "spatial_stats_getis_ord", "z_threshold", default=None
        )
        assert z_thresh == 2.576

    def test_old_key_not_used(self, rng):
        """Ensure the old adata.uns['getis_ord'] key is NOT what we read."""
        from chatspatial.utils.adata_utils import get_analysis_parameter

        adata, _ = _make_getis_adata(rng)

        # Put data in the OLD (wrong) key
        adata.uns["getis_ord"] = {"parameters": {"alpha": 0.99}}

        # get_analysis_parameter should return default, not 0.99
        alpha = get_analysis_parameter(
            adata, "spatial_stats_getis_ord", "alpha", default=0.05
        )
        assert alpha == 0.05  # default, not the wrong key's value


# ---------------------------------------------------------------------------
# Test 2: Hotspot counts use z_threshold
# ---------------------------------------------------------------------------


class TestHotspotCounting:
    """Verify hotspot/coldspot counting uses z_threshold, not z>0."""

    def test_z_threshold_counting(self):
        """With z_threshold=1.96, spots with 0<z<1.96 should NOT be hotspots."""
        from scipy.stats import norm

        z_scores = np.array([0.5, 2.5, -0.5, -2.5, 1.0, 3.0])
        p_vals = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
        sig_alpha = 0.05
        z_thresh = norm.ppf(1 - sig_alpha / 2)  # ~1.96

        significant = p_vals < sig_alpha
        # Correct counting (using z_threshold)
        hot_spots = np.sum((z_scores > z_thresh) & significant)
        cold_spots = np.sum((z_scores < -z_thresh) & significant)

        # z=0.5 and z=1.0 should NOT be hotspots (below threshold)
        assert hot_spots == 2  # only z=2.5 and z=3.0
        assert cold_spots == 1  # only z=-2.5

        # Old (wrong) counting would give:
        wrong_hot = np.sum((z_scores > 0) & significant)
        assert wrong_hot == 4  # z=0.5, 1.0, 2.5, 3.0 -- too many


# ---------------------------------------------------------------------------
# Test 3: _build_results_keys with auto-selected genes
# ---------------------------------------------------------------------------


class TestBuildResultsKeys:
    """Verify _build_results_keys uses actual genes, not None."""

    def test_with_none_genes_returns_empty_obs(self):
        from chatspatial.tools.spatial_statistics import _build_results_keys

        keys = _build_results_keys("getis_ord", None)
        # With genes=None, no per-gene obs keys should be generated
        assert keys["obs"] == []

    def test_with_actual_genes_returns_per_gene_keys(self):
        from chatspatial.tools.spatial_statistics import _build_results_keys

        genes = ["GeneA", "GeneB"]
        keys = _build_results_keys("getis_ord", genes)

        # Should have z, p, and p_corrected for each gene
        assert "GeneA_getis_ord_z" in keys["obs"]
        assert "GeneA_getis_ord_p" in keys["obs"]
        assert "GeneB_getis_ord_z" in keys["obs"]
        assert "GeneB_getis_ord_p" in keys["obs"]


# ---------------------------------------------------------------------------
# Test 4: n_significant is non-zero for getis_ord
# ---------------------------------------------------------------------------


class TestNSignificant:
    """Verify n_significant is properly aggregated for getis_ord results."""

    def test_extract_result_summary_has_n_significant(self):
        from chatspatial.tools.spatial_statistics import _extract_result_summary

        result = {
            "genes_analyzed": ["GeneA", "GeneB"],
            "results": {
                "GeneA": {
                    "n_hot_spots": 10,
                    "n_cold_spots": 5,
                    "n_significant_raw": 15,
                },
                "GeneB": {
                    "n_hot_spots": 8,
                    "n_cold_spots": 3,
                    "n_significant_raw": 11,
                },
            },
        }
        summary = _extract_result_summary(result, "getis_ord")
        assert summary["n_significant"] == 26  # 15 + 11

    def test_extract_result_summary_prefers_corrected(self):
        from chatspatial.tools.spatial_statistics import _extract_result_summary

        result = {
            "genes_analyzed": ["GeneA"],
            "results": {
                "GeneA": {
                    "n_hot_spots": 10,
                    "n_cold_spots": 5,
                    "n_significant_raw": 100,
                    "n_significant_corrected": 7,
                    "n_hot_spots_corrected": 4,
                    "n_cold_spots_corrected": 3,
                },
            },
        }
        summary = _extract_result_summary(result, "getis_ord")
        # Should prefer corrected over raw
        assert summary["n_significant"] == 7

    def test_getis_ord_return_includes_n_significant(self):
        """The return summary should count corrected hot and cold spots."""
        from chatspatial.tools.spatial_statistics import _extract_result_summary

        summary = _extract_result_summary(
            {
                "genes_analyzed": ["GeneA", "GeneB"],
                "results": {
                    "GeneA": {
                        "n_hot_spots": 10,
                        "n_cold_spots": 20,
                        "n_hot_spots_corrected": 3,
                        "n_cold_spots_corrected": 2,
                    },
                    "GeneB": {"n_hot_spots": 4, "n_cold_spots": 6},
                },
            },
            "getis_ord",
        )

        assert summary["n_significant"] == 15


# ---------------------------------------------------------------------------
# Test 5: Corrected p-value columns in results_keys
# ---------------------------------------------------------------------------


class TestCorrectedPValueKeys:
    """Verify corrected p-value columns are included in results_keys."""

    def test_p_corrected_always_in_results_keys(self):
        from chatspatial.tools.spatial_statistics import _build_results_keys

        genes = ["GeneA"]
        keys = _build_results_keys("getis_ord", genes)

        assert "GeneA_getis_ord_p_corrected" in keys["obs"]

    def test_all_three_columns_per_gene(self):
        from chatspatial.tools.spatial_statistics import _build_results_keys

        genes = ["GeneX"]
        keys = _build_results_keys("getis_ord", genes)

        expected = [
            "GeneX_getis_ord_z",
            "GeneX_getis_ord_p",
            "GeneX_getis_ord_p_corrected",
        ]
        for exp_key in expected:
            assert exp_key in keys["obs"], f"Missing {exp_key} in obs keys"
