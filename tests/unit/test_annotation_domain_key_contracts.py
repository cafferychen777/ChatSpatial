"""Contract tests for annotation and spatial domain parametric key generation."""

from __future__ import annotations

from chatspatial.models.data import SpatialDomainParameters
from chatspatial.tools.annotation import _build_annotation_suffix
from chatspatial.tools.spatial_domains import (
    _build_domain_suffix,
    _domain_count_control,
)

# ---------------------------------------------------------------------------
# _build_annotation_suffix
# ---------------------------------------------------------------------------


class TestBuildAnnotationSuffix:
    def test_non_reference_method_uses_method_only(self):
        assert _build_annotation_suffix("sctype", None) == "sctype"
        assert _build_annotation_suffix("cellassign", None) == "cellassign"
        assert _build_annotation_suffix("mllmcelltype", None) == "mllmcelltype"

    def test_reference_method_without_ref_uses_method_only(self):
        """tangram/scanvi/singler without reference_data_id fall back to method."""
        assert _build_annotation_suffix("tangram", None) == "tangram"
        assert _build_annotation_suffix("scanvi", None) == "scanvi"
        assert _build_annotation_suffix("singler", None) == "singler"

    def test_reference_method_with_ref_encodes_ref_id(self):
        assert _build_annotation_suffix("tangram", "ref_pbmc") == "tangram_ref_pbmc"
        assert _build_annotation_suffix("scanvi", "atlas_v2") == "scanvi_atlas_v2"
        assert _build_annotation_suffix("singler", "hpca") == "singler_hpca"

    def test_different_refs_produce_different_suffixes(self):
        s1 = _build_annotation_suffix("scanvi", "refA")
        s2 = _build_annotation_suffix("scanvi", "refB")
        assert s1 != s2

    def test_non_reference_method_ignores_ref_id(self):
        """sctype with a ref_id (shouldn't happen, but be defensive)."""
        assert _build_annotation_suffix("sctype", "some_ref") == "sctype"


# ---------------------------------------------------------------------------
# _build_domain_suffix
# ---------------------------------------------------------------------------


def _params(**kwargs) -> SpatialDomainParameters:
    return SpatialDomainParameters(**kwargs)


class TestDomainCountControl:
    """Which parameter a backend actually obeys drives both keys and warnings."""

    def test_resolution_driven_backends_are_named(self):
        assert _domain_count_control(_params(method="leiden")) == "resolution"
        assert _domain_count_control(_params(method="louvain")) == "resolution"
        assert (
            _domain_count_control(_params(method="banksy"))
            == "banksy_cluster_resolution"
        )

    def test_count_driven_backends_report_none(self):
        for method in ("spagcn", "stagate", "graphst"):
            assert _domain_count_control(_params(method=method)) is None

    def test_aestetik_depends_on_its_clustering_method(self):
        assert (
            _domain_count_control(
                _params(method="aestetik", aestetik_clustering_method="kmeans")
            )
            is None
        )
        assert (
            _domain_count_control(
                _params(method="aestetik", aestetik_clustering_method="leiden")
            )
            == "resolution"
        )


class TestBuildDomainSuffix:
    def test_leiden_encodes_resolution(self):
        assert (
            _build_domain_suffix(_params(method="leiden", resolution=0.5))
            == "leiden_res0_50"
        )
        assert (
            _build_domain_suffix(_params(method="leiden", resolution=1.0))
            == "leiden_res1_00"
        )

    def test_louvain_encodes_resolution(self):
        assert (
            _build_domain_suffix(_params(method="louvain", resolution=0.8))
            == "louvain_res0_80"
        )

    def test_count_driven_backends_encode_n_domains(self):
        assert (
            _build_domain_suffix(_params(method="spagcn", n_domains=7)) == "spagcn_n7"
        )
        assert (
            _build_domain_suffix(_params(method="stagate", n_domains=10))
            == "stagate_n10"
        )

    def test_banksy_encodes_its_own_resolution_not_n_domains(self):
        """Regression: the key advertised an n_domains banksy never honoured."""
        assert (
            _build_domain_suffix(_params(method="banksy", n_domains=12))
            == "banksy_res0_50"
        )
        assert (
            _build_domain_suffix(
                _params(method="banksy", banksy_cluster_resolution=1.2)
            )
            == "banksy_res1_20"
        )

    def test_suffixes_separate_runs_that_differ(self):
        assert _build_domain_suffix(
            _params(method="leiden", resolution=0.5)
        ) != _build_domain_suffix(_params(method="leiden", resolution=1.0))
        assert _build_domain_suffix(
            _params(method="spagcn", n_domains=5)
        ) != _build_domain_suffix(_params(method="spagcn", n_domains=10))
        assert _build_domain_suffix(
            _params(method="leiden", resolution=0.5)
        ) != _build_domain_suffix(_params(method="louvain", resolution=0.5))
