"""Contract tests for differential expression parametric key generation."""

from __future__ import annotations

from chatspatial.tools.differential import _build_de_key


# ---------------------------------------------------------------------------
# _build_de_key unit tests
# ---------------------------------------------------------------------------


def test_build_de_key_all_groups():
    assert _build_de_key("wilcoxon", None, None) == "de_wilcoxon_all_groups"


def test_build_de_key_pairwise():
    assert _build_de_key("wilcoxon", "T", "B") == "de_wilcoxon_T_vs_B"


def test_build_de_key_pairwise_vs_rest():
    assert _build_de_key("t-test", "T", None) == "de_t-test_T_vs_rest"
    assert _build_de_key("t-test", "T", "rest") == "de_t-test_T_vs_rest"


def test_build_de_key_pydeseq2():
    assert _build_de_key("pydeseq2", "treated", "control") == "de_pydeseq2_treated_vs_control"


# ---------------------------------------------------------------------------
# Coexistence: different methods/comparisons produce distinct keys
# ---------------------------------------------------------------------------


def test_different_methods_produce_different_keys():
    """wilcoxon all_groups and t-test all_groups must not collide."""
    k1 = _build_de_key("wilcoxon", None, None)
    k2 = _build_de_key("t-test", None, None)
    assert k1 != k2


def test_different_comparisons_same_method_produce_different_keys():
    """T vs B and T vs rest using the same method must not collide."""
    k1 = _build_de_key("wilcoxon", "T", "B")
    k2 = _build_de_key("wilcoxon", "T", None)
    assert k1 != k2


def test_all_groups_vs_pairwise_produce_different_keys():
    """all_groups and a specific pairwise comparison must not collide."""
    k1 = _build_de_key("wilcoxon", None, None)
    k2 = _build_de_key("wilcoxon", "T", "B")
    assert k1 != k2
