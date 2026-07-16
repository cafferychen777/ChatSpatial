"""Unit tests for CNV analysis routing and infercnvpy contracts."""

from __future__ import annotations

from types import ModuleType, SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from chatspatial.models.data import CNVParameters
from chatspatial.tools import cnv_analysis as cnv
from chatspatial.tools.cnv_analysis import _expand_gene_aligned_layer
from chatspatial.utils.dependency_manager import REnvironment
from chatspatial.utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
)


class DummyCtx:
    def __init__(self, adata):
        self.adata = adata
        self.warnings: list[str] = []
        self.committed = None

    async def get_adata(self, _data_id: str):
        return self.adata

    async def set_adata(self, _data_id: str, adata):
        self.adata = adata
        self.committed = adata

    async def warning(self, msg: str):
        self.warnings.append(msg)


def _add_gene_positions(adata, chromosomes: list[str] | None = None):
    if chromosomes is None:
        chromosomes = ["chr1"] * (adata.n_vars // 2) + ["chr2"] * (
            adata.n_vars - adata.n_vars // 2
        )
    adata.var["chromosome"] = chromosomes
    adata.var["start"] = np.arange(adata.n_vars) * 1000
    adata.var["end"] = adata.var["start"] + 999


def _required_dependency(name: str, *_args, **_kwargs):
    return __import__("sys").modules[name]


def _numbat_allele_data(cell: str) -> pd.DataFrame:
    """Return one valid phased SNP row for Numbat boundary tests."""
    return pd.DataFrame(
        {
            "cell": [cell],
            "snp_id": ["chr1_100_A_G"],
            "CHROM": ["chr1"],
            "POS": [100],
            "cM": [0.1],
            "REF": ["A"],
            "ALT": ["G"],
            "AD": [3],
            "DP": [10],
            "GT": ["0|1"],
        }
    )


@pytest.mark.asyncio
async def test_infer_cnv_rejects_unknown_method_via_runtime_guard(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    params = CNVParameters(
        method="infercnvpy",
        reference_key="cell_type",
        reference_categories=["A"],
    ).model_copy(update={"method": "unknown"})

    with pytest.raises(ParameterError, match="Unknown CNV method"):
        await cnv.infer_cnv("d1", DummyCtx(adata), params)


@pytest.mark.asyncio
async def test_infer_cnv_infercnvpy_success_sparse_stats_and_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    _add_gene_positions(adata)
    captured: dict[str, object] = {}

    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **_kwargs):
        adata_obj.obsm["X_cnv"] = sparse.csr_matrix(
            np.tile(np.array([0.0, 1.0, 0.0, 2.0]), (adata_obj.n_obs, 1))
        )
        adata_obj.uns["cnv"] = {"ok": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(
        cnv,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_args, **_kwargs: [])

    ctx = DummyCtx(adata)
    out = await cnv.infer_cnv(
        "d1",
        ctx,
        CNVParameters(
            method="infercnvpy",
            reference_key="cell_type",
            reference_categories=["A"],
            cluster_cells=False,
            dendrogram=False,
        ),
    )

    assert out.method == "infercnvpy"
    assert out.cnv_score_key == "X_cnv"
    assert out.n_chromosomes == 2
    assert "mean_cnv" in out.statistics
    assert "std_cnv" in out.statistics
    assert "median_cnv" in out.statistics
    assert ctx.committed is not adata
    assert "cnv_analysis_infercnvpy_A" not in adata.uns
    assert "cnv_analysis_infercnvpy_A" in ctx.committed.uns
    assert captured["analysis_name"] == "cnv_infercnvpy_A"
    assert captured["results_keys"]["uns"] == ["cnv", "cnv_analysis_infercnvpy_A"]
    assert captured["results_keys"]["obsm"] == ["X_cnv"]


@pytest.mark.asyncio
async def test_infer_cnv_late_failure_does_not_publish_partial_results(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    _add_gene_positions(adata)
    ctx = DummyCtx(adata)

    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **_kwargs):
        adata_obj.obsm["X_cnv"] = np.ones((adata_obj.n_obs, 2), dtype=float)
        adata_obj.uns["cnv"] = {"partial": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(
        cnv,
        "export_analysis_result",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("export failed")),
    )

    with pytest.raises(RuntimeError, match="export failed"):
        await cnv.infer_cnv(
            "d1",
            ctx,
            CNVParameters(
                method="infercnvpy",
                reference_key="cell_type",
                reference_categories=["A"],
                cluster_cells=False,
                dendrogram=False,
            ),
        )

    assert ctx.committed is None
    assert "X_cnv" not in adata.obsm
    assert "cnv" not in adata.uns


@pytest.mark.asyncio
async def test_infer_cnv_infercnvpy_workspace_isolation_avoids_leaking_temp_mutations(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    _add_gene_positions(adata)
    adata.obsm["keep"] = np.ones((adata.n_obs, 2), dtype=float)

    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **_kwargs):
        adata_obj.obs["_tmp_obs"] = "x"
        adata_obj.var["_tmp_var"] = "y"
        adata_obj.uns["_tmp_uns"] = {"z": 1}
        adata_obj.obsm["X_cnv"] = sparse.csr_matrix(
            np.tile(np.array([0.0, 1.0]), (adata_obj.n_obs, 1))
        )
        adata_obj.uns["cnv"] = {"ok": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_args, **_kwargs: None)

    ctx = DummyCtx(adata)
    await cnv.infer_cnv(
        "d1",
        ctx,
        CNVParameters(
            method="infercnvpy",
            reference_key="cell_type",
            reference_categories=["A"],
            cluster_cells=False,
            dendrogram=False,
        ),
    )

    assert "_tmp_obs" not in adata.obs.columns
    assert "_tmp_var" not in adata.var.columns
    assert "_tmp_uns" not in adata.uns
    assert "keep" in adata.obsm
    assert "_tmp_obs" not in ctx.committed.obs.columns
    assert "_tmp_var" not in ctx.committed.var.columns
    assert "_tmp_uns" not in ctx.committed.uns
    assert "keep" in ctx.committed.obsm
    assert "X_cnv" in ctx.committed.obsm


@pytest.mark.asyncio
async def test_infer_cnv_infercnvpy_rejects_missing_gene_positions(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30

    fake_infercnvpy = ModuleType("infercnvpy")
    fake_infercnvpy.tl = SimpleNamespace(
        infercnv=lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("boom"))
    )
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)

    with pytest.raises(DataCompatibilityError, match="Missing columns"):
        await cnv.infer_cnv(
            "d1",
            DummyCtx(adata),
            CNVParameters(
                method="infercnvpy",
                reference_key="cell_type",
                reference_categories=["A"],
            ),
        )


@pytest.mark.asyncio
async def test_infer_cnv_rejects_missing_reference_categories(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30

    with pytest.raises(ParameterError, match="Reference categories"):
        await cnv.infer_cnv(
            "d2",
            DummyCtx(adata),
            CNVParameters(
                method="infercnvpy",
                reference_key="cell_type",
                reference_categories=["MISSING"],
            ),
        )


def test_infer_cnv_numbat_requires_allele_dataframe(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30

    monkeypatch.setattr(
        cnv,
        "validate_r_environment",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("R should load only after allele-data validation")
        ),
    )

    with pytest.raises(ParameterError, match="numbat_allele_data_raw"):
        cnv._infer_cnv_numbat(
            "d3",
            adata,
            CNVParameters(
                method="numbat",
                reference_key="cell_type",
                reference_categories=["A"],
            ),
            DummyCtx(adata),
        )


def test_infer_cnv_numbat_dependency_error_when_rpy2_missing(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    adata.uns["numbat_allele_data_raw"] = _numbat_allele_data(adata.obs_names[0])
    monkeypatch.setattr(
        cnv,
        "validate_r_environment",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("rpy2 is required")
        ),
    )

    with pytest.raises(DependencyError, match="rpy2 is required"):
        cnv._infer_cnv_numbat(
            "d4",
            adata,
            CNVParameters(
                method="numbat",
                reference_key="cell_type",
                reference_categories=["A"],
            ),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_infer_cnv_infercnvpy_without_cnv_matrix_returns_non_visual_result(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    _add_gene_positions(adata)

    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **_kwargs):
        adata_obj.uns["cnv"] = {"ok": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_args, **_kwargs: None)

    out = await cnv.infer_cnv(
        "d5",
        DummyCtx(adata),
        CNVParameters(
            method="infercnvpy",
            reference_key="cell_type",
            reference_categories=["A"],
            cluster_cells=False,
            dendrogram=False,
        ),
    )

    assert out.cnv_score_key is None
    assert out.visualization_available is False
    assert out.statistics["n_reference_cells"] == 30


def _install_fake_rpy2_stack(monkeypatch: pytest.MonkeyPatch):
    class _CM:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class _Converter:
        def __add__(self, _other):
            return self

    fake_anndata2ri = ModuleType("anndata2ri")
    fake_anndata2ri.converter = _Converter()
    monkeypatch.setitem(__import__("sys").modules, "anndata2ri", fake_anndata2ri)

    fake_openrlib_mod = ModuleType("rpy2.rinterface_lib")
    fake_openrlib_mod.openrlib = SimpleNamespace(rlock=_CM())
    monkeypatch.setitem(
        __import__("sys").modules, "rpy2.rinterface_lib", fake_openrlib_mod
    )

    fake_robj = ModuleType("rpy2.robjects")
    fake_robj.r = lambda *_a, **_k: None
    fake_robj.globalenv = {}
    fake_robj.conversion = SimpleNamespace(localconverter=lambda *_a, **_k: _CM())
    fake_robj.default_converter = _Converter()
    fake_robj.numpy2ri = SimpleNamespace(
        converter=_Converter(), deactivate=lambda: None
    )
    fake_robj.pandas2ri = SimpleNamespace(
        converter=_Converter(), deactivate=lambda: None
    )
    monkeypatch.setitem(__import__("sys").modules, "rpy2.robjects", fake_robj)

    # Also patch top-level package so `import rpy2.robjects` resolves to fake
    fake_rpy2_pkg = ModuleType("rpy2")
    fake_rpy2_pkg.robjects = fake_robj
    fake_rpy2_pkg.rinterface_lib = fake_openrlib_mod
    monkeypatch.setitem(__import__("sys").modules, "rpy2", fake_rpy2_pkg)

    environment = REnvironment(
        robjects=fake_robj,
        pandas2ri=fake_robj.pandas2ri,
        numpy2ri=fake_robj.numpy2ri,
        packages=SimpleNamespace(importr=lambda _package: object()),
        conversion=fake_robj.conversion,
        openrlib=fake_openrlib_mod.openrlib,
        anndata2ri=fake_anndata2ri,
    )
    monkeypatch.setattr(
        cnv,
        "validate_r_environment",
        lambda *_args, **_kwargs: environment,
    )
    return environment


@pytest.mark.asyncio
async def test_infer_cnv_infercnvpy_excludes_chromosomes_before_inference(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    _add_gene_positions(adata, ["chr1"] * 10 + ["chr2"] * 10 + ["chrM"] * 4)

    seen = {"n_vars": None, "exclude_chromosomes": "unset"}
    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **_kwargs):
        seen["n_vars"] = adata_obj.n_vars
        seen["exclude_chromosomes"] = _kwargs.get("exclude_chromosomes")
        adata_obj.obsm["X_cnv"] = np.ones((adata_obj.n_obs, 3), dtype=float)
        adata_obj.uns["cnv"] = {"ok": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_a, **_k: [])
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_a, **_k: None)

    out = await cnv.infer_cnv(
        "d6",
        DummyCtx(adata),
        CNVParameters(
            method="infercnvpy",
            reference_key="cell_type",
            reference_categories=["A"],
            exclude_chromosomes=["chrM"],
            cluster_cells=False,
            dendrogram=False,
        ),
    )

    assert seen["n_vars"] == 20
    assert seen["exclude_chromosomes"] is None
    assert out.n_genes_analyzed == 20


@pytest.mark.asyncio
async def test_infer_cnv_none_exclusion_keeps_sex_chromosomes(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    _add_gene_positions(adata, ["chr1"] * 20 + ["chrX"] * 2 + ["chrY"] * 2)
    captured: dict[str, object] = {}

    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **kwargs):
        captured["genes"] = adata_obj.var_names.tolist()
        captured["exclude_chromosomes"] = kwargs["exclude_chromosomes"]
        adata_obj.obsm["X_cnv"] = np.ones((adata_obj.n_obs, 2), dtype=float)
        adata_obj.uns["cnv"] = {"ok": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_a, **_k: [])
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_a, **_k: None)

    out = await cnv.infer_cnv(
        "d_keep_sex",
        DummyCtx(adata),
        CNVParameters(
            method="infercnvpy",
            reference_key="cell_type",
            reference_categories=["A"],
            exclude_chromosomes=None,
            cluster_cells=False,
            dendrogram=False,
        ),
    )

    assert captured["genes"] == adata.var_names.tolist()
    assert captured["exclude_chromosomes"] is None
    assert out.n_genes_analyzed == adata.n_vars


@pytest.mark.parametrize("sparse_input", [False, True])
def test_expand_gene_aligned_layer_maps_columns_by_gene_name(sparse_input):
    source_genes = pd.Index(["g3", "g1"])
    target_genes = pd.Index(["g1", "g2", "g3"])
    values = np.array([[30.0, 10.0], [31.0, 11.0]], dtype=np.float32)
    layer = sparse.csr_matrix(values) if sparse_input else values

    expanded = _expand_gene_aligned_layer(layer, source_genes, target_genes, n_obs=2)
    if sparse.issparse(expanded):
        expanded = expanded.toarray()

    np.testing.assert_array_equal(
        expanded,
        np.array([[10.0, 0.0, 30.0], [11.0, 0.0, 31.0]], dtype=np.float32),
    )


@pytest.mark.asyncio
async def test_infer_cnv_infercnvpy_cluster_and_dendrogram_failures_emit_warnings(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    _add_gene_positions(adata)
    ctx = DummyCtx(adata)

    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **_kwargs):
        adata_obj.obsm["X_cnv"] = np.ones((adata_obj.n_obs, 2), dtype=float)
        adata_obj.uns["cnv"] = {"ok": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)

    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_a, **_k: [])
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_a, **_k: None)

    monkeypatch.setattr(
        cnv.sc.pp,
        "neighbors",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("nb fail")),
    )
    monkeypatch.setattr(
        cnv.sc.tl,
        "dendrogram",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("dendro fail")),
    )

    out = await cnv.infer_cnv(
        "d7",
        ctx,
        CNVParameters(
            method="infercnvpy",
            reference_key="cell_type",
            reference_categories=["A"],
            cluster_cells=True,
            dendrogram=True,
        ),
    )

    assert out.method == "infercnvpy"
    assert any("Failed to cluster cells by CNV" in w for w in ctx.warnings)
    assert any("Failed to compute dendrogram" in w for w in ctx.warnings)


def test_infer_cnv_numbat_rejects_allele_dataframe_with_missing_required_columns(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    adata.uns["numbat_allele_data_raw"] = pd.DataFrame(
        {
            "cell": [adata.obs_names[0]],
            "CHROM": ["chr1"],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["G"],
            "AD": [3],
        }
    )
    _install_fake_rpy2_stack(monkeypatch)

    with pytest.raises(ParameterError, match="missing required columns"):
        cnv._infer_cnv_numbat(
            "d8",
            adata,
            CNVParameters(
                method="numbat",
                reference_key="cell_type",
                reference_categories=["A"],
            ),
            DummyCtx(adata),
        )


@pytest.mark.parametrize("invalid", [None, [], pd.DataFrame()])
def test_validate_numbat_allele_data_rejects_invalid_containers(invalid):
    with pytest.raises(ParameterError):
        cnv._validate_numbat_allele_data(invalid)


def test_validate_and_align_numbat_outputs_preserves_cell_identity():
    obs_names = pd.Index(["s1", "s2", "s3"])
    clone_post = pd.DataFrame(
        {
            "cell": ["s3", "s1"],
            "clone_opt": ["c3", "c1"],
            "p_cnv": [0.8, 0.2],
            "compartment_opt": ["tumor", "normal"],
        }
    )
    geno = pd.DataFrame(
        {
            "cell": ["s3", "s1"],
            "seg1": [0.9, 0.1],
            "seg2": [0.7, 0.3],
        }
    )

    normalized, aligned_clones, aligned_genotypes = (
        cnv._validate_and_align_numbat_outputs(clone_post, geno, obs_names)
    )

    assert normalized["cell"].tolist() == ["s3", "s1"]
    assert aligned_clones.index.tolist() == obs_names.tolist()
    assert aligned_clones["clone_opt"].tolist() == ["c1", "unassigned", "c3"]
    np.testing.assert_allclose(aligned_genotypes[[0, 2]], [[0.1, 0.3], [0.9, 0.7]])
    assert np.isnan(aligned_genotypes[1]).all()


@pytest.mark.parametrize(
    ("table_name", "mutator", "message"),
    [
        (
            "clone",
            lambda table: pd.concat([table, table.iloc[[0]]], ignore_index=True),
            "duplicate cell identifiers",
        ),
        (
            "geno",
            lambda table: pd.concat([table, table.iloc[[0]]], ignore_index=True),
            "duplicate cell identifiers",
        ),
        (
            "clone",
            lambda table: table.drop(columns="p_cnv"),
            "missing required columns",
        ),
        ("geno", lambda table: table[["cell"]], "no CNV segment columns"),
        (
            "geno",
            lambda table: table.assign(cell=["unknown", "s2"]),
            "cells absent from the dataset",
        ),
        (
            "clone",
            lambda table: table.assign(p_cnv=1.5),
            r"probabilities in \[0, 1\]",
        ),
    ],
)
def test_validate_and_align_numbat_outputs_rejects_ambiguous_axes(
    table_name, mutator, message
):
    obs_names = pd.Index(["s1", "s2"])
    clone_post = pd.DataFrame(
        {
            "cell": obs_names,
            "clone_opt": ["c1", "c2"],
            "p_cnv": [0.2, 0.8],
            "compartment_opt": ["normal", "tumor"],
        }
    )
    geno = pd.DataFrame({"cell": obs_names, "seg1": [0.1, 0.9]})
    if table_name == "clone":
        clone_post = mutator(clone_post)
    else:
        geno = mutator(geno)

    with pytest.raises(DataCompatibilityError, match=message):
        cnv._validate_and_align_numbat_outputs(clone_post, geno, obs_names)


def test_infer_cnv_numbat_requires_nonempty_reference_cells(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    adata.uns["numbat_allele_data_raw"] = _numbat_allele_data(adata.obs_names[0])
    _install_fake_rpy2_stack(monkeypatch)

    with pytest.raises(ParameterError, match="No reference cells found"):
        cnv._infer_cnv_numbat(
            "d9",
            adata,
            CNVParameters(
                method="numbat",
                reference_key="cell_type",
                reference_categories=["MISSING"],
            ),
            DummyCtx(adata),
        )


def _install_fake_rpy2_stack_with_runner(
    monkeypatch: pytest.MonkeyPatch,
    runner,
):
    class _CM:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class _Converter:
        def __add__(self, _other):
            return self

    fake_anndata2ri = ModuleType("anndata2ri")
    fake_anndata2ri.converter = _Converter()
    monkeypatch.setitem(__import__("sys").modules, "anndata2ri", fake_anndata2ri)

    fake_openrlib_mod = ModuleType("rpy2.rinterface_lib")
    fake_openrlib_mod.openrlib = SimpleNamespace(rlock=_CM())
    monkeypatch.setitem(
        __import__("sys").modules, "rpy2.rinterface_lib", fake_openrlib_mod
    )

    fake_robj = ModuleType("rpy2.robjects")
    fake_robj.globalenv = {}

    def _fake_r(code: str):
        if "run_numbat" in code:
            runner(fake_robj.globalenv)
        return None

    fake_robj.r = _fake_r
    fake_robj.conversion = SimpleNamespace(localconverter=lambda *_a, **_k: _CM())
    fake_robj.default_converter = _Converter()
    fake_robj.numpy2ri = SimpleNamespace(
        converter=_Converter(), deactivate=lambda: None
    )
    fake_robj.pandas2ri = SimpleNamespace(
        converter=_Converter(), deactivate=lambda: None
    )
    monkeypatch.setitem(__import__("sys").modules, "rpy2.robjects", fake_robj)

    # Also patch top-level package so `import rpy2.robjects` resolves to fake
    fake_rpy2_pkg = ModuleType("rpy2")
    fake_rpy2_pkg.robjects = fake_robj
    fake_rpy2_pkg.rinterface_lib = fake_openrlib_mod
    monkeypatch.setitem(__import__("sys").modules, "rpy2", fake_rpy2_pkg)

    environment = REnvironment(
        robjects=fake_robj,
        pandas2ri=fake_robj.pandas2ri,
        numpy2ri=fake_robj.numpy2ri,
        packages=SimpleNamespace(importr=lambda _package: object()),
        conversion=fake_robj.conversion,
        openrlib=fake_openrlib_mod.openrlib,
        anndata2ri=fake_anndata2ri,
    )
    monkeypatch.setattr(
        cnv,
        "validate_r_environment",
        lambda *_args, **_kwargs: environment,
    )
    return environment


@pytest.mark.asyncio
async def test_infer_cnv_numbat_success_parses_outputs_and_writes_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, tmp_path
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    adata.uns["numbat_allele_data_raw"] = _numbat_allele_data(adata.obs_names[0])

    def _runner(env: dict):
        out_dir = env["out_dir"]
        cell_barcodes = list(env["cell_barcodes"])

        analyzed_cells = cell_barcodes[::2]
        clone_post = pd.DataFrame(
            {
                "cell": analyzed_cells,
                "clone_opt": ["c1"] * len(analyzed_cells),
                "p_cnv": [0.7] * len(analyzed_cells),
                "compartment_opt": ["tumor"] * len(analyzed_cells),
            }
        )
        clone_post.to_csv(f"{out_dir}/clone_post_2.tsv", sep="\t", index=False)

        geno = pd.DataFrame(
            {
                "cell": analyzed_cells,
                "seg1": np.ones(len(analyzed_cells)),
                "seg2": np.zeros(len(analyzed_cells)),
            }
        )
        geno.to_csv(f"{out_dir}/geno_2.tsv", sep="\t", index=False)

        segs = pd.DataFrame({"segment": ["s1"], "chr": ["chr1"], "note": [None]})
        segs.to_csv(
            f"{out_dir}/segs_consensus_2.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

        with open(f"{out_dir}/tree_final_2.rds", "wb") as f:
            f.write(b"tree")

    _install_fake_rpy2_stack_with_runner(monkeypatch, _runner)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_a, **_k: [])
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_a, **_k: None)

    def _mkdtemp_success(prefix, dir):
        _ = prefix, dir
        p = tmp_path / "numbat_out"
        p.mkdir(parents=True, exist_ok=True)
        return str(p)

    monkeypatch.setattr(__import__("tempfile"), "mkdtemp", _mkdtemp_success)

    out = cnv._infer_cnv_numbat(
        "d10",
        adata,
        CNVParameters(
            method="numbat",
            reference_key="cell_type",
            reference_categories=["A"],
        ),
        DummyCtx(adata),
    )

    assert out.method == "numbat"
    assert out.cnv_score_key == "X_cnv_numbat"
    assert out.statistics["n_segments"] == 2
    assert out.statistics["n_clones"] == 1
    assert out.statistics["n_cells_analyzed"] == adata.n_obs // 2
    assert out.statistics["n_cells_missing"] == adata.n_obs // 2
    assert np.isfinite(out.statistics["mean_cnv"])
    assert adata.obsm["X_cnv_numbat"].shape == (adata.n_obs, 2)
    assert np.isnan(adata.obsm["X_cnv_numbat"][1]).all()
    assert adata.obs["numbat_clone"].iloc[1] == "unassigned"
    assert "numbat_clone" in adata.obs
    assert "numbat_segments" in adata.uns
    assert "numbat_phylogeny" in adata.uns
    phylogeny = adata.uns["numbat_phylogeny"]
    assert phylogeny == {
        "generated": True,
        "retained": False,
        "source_filename": "tree_final_2.rds",
        "tree_type": "tbl_graph",
    }
    assert "tree_file" not in phylogeny
    assert not (tmp_path / "numbat_out").exists()


def test_infer_cnv_numbat_missing_output_files_preserves_data_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, tmp_path
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    adata.uns["numbat_allele_data_raw"] = _numbat_allele_data(adata.obs_names[0])

    def _runner(_env: dict):
        return None

    _install_fake_rpy2_stack_with_runner(monkeypatch, _runner)

    def _mkdtemp_missing(prefix, dir):
        _ = prefix, dir
        p = tmp_path / "numbat_out_missing"
        p.mkdir(parents=True, exist_ok=True)
        return str(p)

    monkeypatch.setattr(__import__("tempfile"), "mkdtemp", _mkdtemp_missing)

    with pytest.raises(DataNotFoundError, match="No Numbat clone_post output found"):
        cnv._infer_cnv_numbat(
            "d11",
            adata,
            CNVParameters(
                method="numbat",
                reference_key="cell_type",
                reference_categories=["A"],
            ),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_infer_cnv_infercnvpy_cluster_and_dendrogram_success_copies_outputs(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    _add_gene_positions(adata)
    captured: dict[str, object] = {}

    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **_kwargs):
        # 75% zeros so sparse median branch uses exact zero
        arr = np.tile(np.array([0.0, 0.0, 0.0, 2.0]), (adata_obj.n_obs, 1))
        adata_obj.obsm["X_cnv"] = sparse.csr_matrix(arr)
        adata_obj.uns["cnv"] = {"ok": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_a, **_k: [])
    monkeypatch.setattr(
        cnv, "store_analysis_metadata", lambda _adata, **kwargs: captured.update(kwargs)
    )
    monkeypatch.setattr(cnv.sc.pp, "neighbors", lambda *_a, **_k: None)

    def _fake_leiden(adata_obj, key_added="cnv_clusters"):
        adata_obj.obs[key_added] = pd.Categorical(
            ["c0"] * (adata_obj.n_obs // 2)
            + ["c1"] * (adata_obj.n_obs - adata_obj.n_obs // 2)
        )

    monkeypatch.setattr(cnv.sc.tl, "leiden", _fake_leiden)
    monkeypatch.setattr(
        cnv.sc.tl,
        "dendrogram",
        lambda adata_obj, groupby="cnv_clusters": adata_obj.uns.__setitem__(
            f"dendrogram_{groupby}", {"linkage": "ok"}
        ),
    )

    ctx = DummyCtx(adata)
    out = await cnv.infer_cnv(
        "d12",
        ctx,
        CNVParameters(
            method="infercnvpy",
            reference_key="cell_type",
            reference_categories=["A"],
            cluster_cells=True,
            dendrogram=True,
        ),
    )

    assert out.statistics["median_cnv"] == 0.0
    assert "cnv_clusters" not in adata.obs
    assert "dendrogram_cnv_clusters" not in adata.uns
    assert "cnv_clusters" in ctx.committed.obs
    assert "dendrogram_cnv_clusters" in ctx.committed.uns
    assert "obs" in captured["results_keys"]
    assert "cnv_clusters" in captured["results_keys"]["obs"]
    assert "dendrogram_cnv_clusters" in captured["results_keys"]["uns"]


@pytest.mark.asyncio
async def test_infer_cnv_infercnvpy_uses_cnv_layer_when_obsm_missing(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    _add_gene_positions(adata)

    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **_kwargs):
        adata_obj.layers["cnv"] = np.ones(
            (adata_obj.n_obs, adata_obj.n_vars), dtype=float
        )
        adata_obj.uns["cnv"] = {"ok": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_a, **_k: [])
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_a, **_k: None)

    out = await cnv.infer_cnv(
        "d13",
        DummyCtx(adata),
        CNVParameters(
            method="infercnvpy",
            reference_key="cell_type",
            reference_categories=["A"],
            cluster_cells=False,
            dendrogram=False,
        ),
    )

    assert out.cnv_score_key == "cnv"
    assert out.n_chromosomes == 2


def test_infer_cnv_numbat_dependency_error_when_r_package_unavailable(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    adata.uns["numbat_allele_data_raw"] = _numbat_allele_data(adata.obs_names[0])
    monkeypatch.setattr(
        cnv,
        "validate_r_environment",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("R package 'numbat' not installed")
        ),
    )

    with pytest.raises(DependencyError, match="R package 'numbat' not installed"):
        cnv._infer_cnv_numbat(
            "d14",
            adata,
            CNVParameters(
                method="numbat",
                reference_key="cell_type",
                reference_categories=["A"],
            ),
            DummyCtx(adata),
        )


def test_infer_cnv_numbat_missing_geno_file_preserves_data_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, tmp_path
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    adata.uns["numbat_allele_data_raw"] = _numbat_allele_data(adata.obs_names[0])

    def _runner(env: dict):
        out_dir = env["out_dir"]
        cell_barcodes = list(env["cell_barcodes"])
        clone_post = pd.DataFrame(
            {
                "cell": cell_barcodes,
                "clone_opt": ["c1"] * len(cell_barcodes),
                "p_cnv": [0.7] * len(cell_barcodes),
                "compartment_opt": ["tumor"] * len(cell_barcodes),
            }
        )
        clone_post.to_csv(f"{out_dir}/clone_post_2.tsv", sep="\t", index=False)
        # Intentionally do not write geno_2.tsv

    _install_fake_rpy2_stack_with_runner(monkeypatch, _runner)

    def _mkdtemp(prefix, dir):
        _ = prefix, dir
        p = tmp_path / "numbat_out_geno_missing"
        p.mkdir(parents=True, exist_ok=True)
        return str(p)

    monkeypatch.setattr(__import__("tempfile"), "mkdtemp", _mkdtemp)

    with pytest.raises(DataNotFoundError, match="geno_2.tsv"):
        cnv._infer_cnv_numbat(
            "d15",
            adata,
            CNVParameters(
                method="numbat",
                reference_key="cell_type",
                reference_categories=["A"],
            ),
            DummyCtx(adata),
        )


def test_infer_cnv_numbat_cell_mismatch_preserves_compatibility_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, tmp_path
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    adata.uns["numbat_allele_data_raw"] = _numbat_allele_data(adata.obs_names[0])

    def _runner(env: dict):
        out_dir = env["out_dir"]
        cell_barcodes = list(env["cell_barcodes"])
        clone_post = pd.DataFrame(
            {
                "cell": cell_barcodes,
                "clone_opt": ["c1"] * len(cell_barcodes),
                "p_cnv": [0.7] * len(cell_barcodes),
                "compartment_opt": ["tumor"] * len(cell_barcodes),
            }
        )
        clone_post.to_csv(f"{out_dir}/clone_post_2.tsv", sep="\t", index=False)

        bad_cells = cell_barcodes[:-1] + ["UNKNOWN_CELL"]
        geno = pd.DataFrame({"cell": bad_cells, "seg1": np.ones(len(bad_cells))})
        geno.to_csv(f"{out_dir}/geno_2.tsv", sep="\t", index=False)

    _install_fake_rpy2_stack_with_runner(monkeypatch, _runner)

    def _mkdtemp(prefix, dir):
        _ = prefix, dir
        p = tmp_path / "numbat_out_bad_cells"
        p.mkdir(parents=True, exist_ok=True)
        return str(p)

    monkeypatch.setattr(__import__("tempfile"), "mkdtemp", _mkdtemp)

    with pytest.raises(
        DataCompatibilityError, match="geno output contains cells absent"
    ):
        cnv._infer_cnv_numbat(
            "d16",
            adata,
            CNVParameters(
                method="numbat",
                reference_key="cell_type",
                reference_categories=["A"],
            ),
            DummyCtx(adata),
        )


def test_infer_cnv_numbat_cleanup_failure_is_swallowed(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, tmp_path
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    adata.uns["numbat_allele_data_raw"] = _numbat_allele_data(adata.obs_names[0])

    def _runner(_env: dict):
        return None  # Missing outputs -> ProcessingError path

    _install_fake_rpy2_stack_with_runner(monkeypatch, _runner)

    def _mkdtemp(prefix, dir):
        _ = prefix, dir
        p = tmp_path / "numbat_out_cleanup_fail"
        p.mkdir(parents=True, exist_ok=True)
        return str(p)

    monkeypatch.setattr(__import__("tempfile"), "mkdtemp", _mkdtemp)
    monkeypatch.setattr(
        __import__("shutil"),
        "rmtree",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("rm fail")),
    )

    with pytest.raises(DataNotFoundError, match="No Numbat clone_post output found"):
        cnv._infer_cnv_numbat(
            "d17",
            adata,
            CNVParameters(
                method="numbat",
                reference_key="cell_type",
                reference_categories=["A"],
            ),
            DummyCtx(adata),
        )


# =============================================================================
# Issue 1 regression: layers["cnv"] shape mismatch after exclude_chromosomes
# =============================================================================


@pytest.mark.asyncio
async def test_infer_cnv_layers_cnv_padded_after_exclude_chromosomes(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """When exclude_chromosomes filters genes, layers['cnv'] must be padded to
    match original adata shape instead of crashing with ValueError."""
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    # 20 genes on chr1/chr2, 4 on chrM (will be excluded)
    _add_gene_positions(adata, ["chr1"] * 10 + ["chr2"] * 10 + ["chrM"] * 4)
    original_n_vars = adata.n_vars  # 24

    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **_kwargs):
        # infercnvpy puts results in layers["cnv"] (gene-level, not window)
        adata_obj.layers["cnv"] = np.ones(
            (adata_obj.n_obs, adata_obj.n_vars), dtype=np.float32
        )
        adata_obj.uns["cnv"] = {"ok": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_a, **_k: [])
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_a, **_k: None)

    # This used to raise ValueError: incorrect shape
    ctx = DummyCtx(adata)
    out = await cnv.infer_cnv(
        "d_pad",
        ctx,
        CNVParameters(
            method="infercnvpy",
            reference_key="cell_type",
            reference_categories=["A"],
            exclude_chromosomes=["chrM"],
            cluster_cells=False,
            dendrogram=False,
        ),
    )

    assert out.cnv_score_key == "cnv"
    assert "cnv" not in adata.layers
    result_adata = ctx.committed
    assert "cnv" in result_adata.layers
    # Layer must match original shape
    assert result_adata.layers["cnv"].shape == (adata.n_obs, original_n_vars)
    # chrM columns (indices 20-23) should be zero-padded
    cnv_layer = result_adata.layers["cnv"]
    if hasattr(cnv_layer, "toarray"):
        cnv_layer = cnv_layer.toarray()
    # Analyzed genes (chr1+chr2) should have value 1.0
    assert cnv_layer[:, :20].sum() > 0
    # Excluded genes (chrM) should be zero
    np.testing.assert_array_equal(cnv_layer[:, 20:], 0.0)


@pytest.mark.asyncio
async def test_infer_cnv_layers_cnv_sparse_padded_after_exclude(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """Same padding must work when infercnvpy returns sparse layers['cnv']."""
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    _add_gene_positions(adata, ["chr1"] * 10 + ["chr2"] * 10 + ["chrM"] * 4)

    fake_infercnvpy = ModuleType("infercnvpy")

    def _fake_infercnv(adata_obj, **_kwargs):
        adata_obj.layers["cnv"] = sparse.csr_matrix(
            np.ones((adata_obj.n_obs, adata_obj.n_vars), dtype=np.float32)
        )
        adata_obj.uns["cnv"] = {"ok": True}

    fake_infercnvpy.tl = SimpleNamespace(infercnv=_fake_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_a, **_k: [])
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_a, **_k: None)

    ctx = DummyCtx(adata)
    out = await cnv.infer_cnv(
        "d_sparse_pad",
        ctx,
        CNVParameters(
            method="infercnvpy",
            reference_key="cell_type",
            reference_categories=["A"],
            exclude_chromosomes=["chrM"],
            cluster_cells=False,
            dendrogram=False,
        ),
    )

    assert out.cnv_score_key == "cnv"
    assert "cnv" not in adata.layers
    cnv_layer = ctx.committed.layers["cnv"]
    assert cnv_layer.shape == (adata.n_obs, adata.n_vars)
    assert sparse.issparse(cnv_layer)
    dense = cnv_layer.toarray()
    np.testing.assert_array_equal(dense[:, 20:], 0.0)


# =============================================================================
# Issue 2 regression: missing gene positions fail before inference
# =============================================================================


@pytest.mark.asyncio
async def test_infer_cnv_rejects_exclude_chromosomes_without_gene_positions(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    if "chromosome" in adata.var.columns:
        del adata.var["chromosome"]

    fake_infercnvpy = ModuleType("infercnvpy")
    fake_infercnvpy.tl = SimpleNamespace(
        infercnv=lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("boom"))
    )
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", fake_infercnvpy)
    monkeypatch.setattr(cnv, "require", _required_dependency)

    with pytest.raises(DataCompatibilityError, match="Missing columns"):
        await cnv.infer_cnv(
            "d_warn",
            DummyCtx(adata),
            CNVParameters(
                method="infercnvpy",
                reference_key="cell_type",
                reference_categories=["A"],
                exclude_chromosomes=["chrM"],
                cluster_cells=False,
                dendrogram=False,
            ),
        )
