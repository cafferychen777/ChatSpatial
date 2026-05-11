"""Additional integration contracts for server wrappers around analysis tools."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest

from chatspatial.models.analysis import (
    AnnotationResult,
    CNVResult,
    CellCommunicationResult,
    ConditionComparisonResult,
    DeconvolutionResult,
    PreprocessingResult,
    RNAVelocityResult,
    SpatialDomainResult,
    SpatialVariableGenesResult,
    TrajectoryResult,
)
from chatspatial.models.data import (
    AnnotationParameters,
    CNVParameters,
    CellCommunicationParameters,
    ConditionComparisonParameters,
    DeconvolutionParameters,
    PreprocessingParameters,
    RNAVelocityParameters,
    SpatialDomainParameters,
    SpatialVariableGenesParameters,
    TrajectoryParameters,
    VisualizationParameters,
)
from chatspatial.tools.embeddings import EmbeddingParameters
from chatspatial.server import (
    analyze_cell_communication,
    analyze_cnv,
    annotate_cell_types,
    preprocess_data,
    analyze_trajectory_data,
    analyze_velocity_data,
    compare_conditions,
    compute_embeddings,
    data_manager,
    deconvolve_data,
    find_spatial_genes,
    identify_spatial_domains,
    visualize_data,
)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_compute_embeddings_returns_model_dump_and_passes_params(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    calls: dict[str, object] = {}

    class FakeEmbeddingResult:
        def model_dump(self):
            return {"ok": True, "computed": ["pca", "neighbors"]}

    async def fake_compute(data_id, ctx, params):
        calls["data_id"] = data_id
        calls["ctx"] = ctx
        calls["params"] = params
        return FakeEmbeddingResult()

    fake_module = SimpleNamespace(
        EmbeddingParameters=EmbeddingParameters,
        compute_embeddings=fake_compute,
    )
    monkeypatch.setitem(sys.modules, "chatspatial.tools.embeddings", fake_module)

    saved: dict[str, object] = {}

    async def fake_save_result(data_id: str, result_type: str, result):
        saved["data_id"] = data_id
        saved["result_type"] = result_type
        saved["result"] = result

    monkeypatch.setattr(data_manager, "save_result", fake_save_result)

    embed_params = EmbeddingParameters(
        clustering_method="louvain",
        n_pcs=22,
        n_neighbors=9,
        force=True,
    )
    result = await compute_embeddings("d1", params=embed_params)

    assert result == {"ok": True, "computed": ["pca", "neighbors"]}
    assert calls["data_id"] == "d1"
    params = calls["params"]
    assert isinstance(params, EmbeddingParameters)
    assert params.clustering_method == "louvain"
    assert params.n_pcs == 22
    assert params.n_neighbors == 9
    assert params.force is True
    # Verify embeddings result is cached
    assert saved["data_id"] == "d1"
    assert saved["result_type"] == "embeddings"
    assert saved["result"] == result


@pytest.mark.integration
@pytest.mark.asyncio
async def test_visualize_data_returns_fallback_message_when_tool_returns_none(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    async def fake_visualize(data_id, ctx, params):
        assert isinstance(params, VisualizationParameters)
        return None

    import chatspatial.tools.visualization.main as visualization_main

    monkeypatch.setattr(visualization_main, "visualize_data", fake_visualize)

    result = await visualize_data("d1")
    assert result.startswith("Visualization generation failed")


@pytest.mark.integration
@pytest.mark.asyncio
async def test_compare_conditions_saves_result_with_expected_key(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    calls: dict[str, object] = {}

    async def fake_compare(data_id, ctx, params):
        calls["data_id"] = data_id
        calls["params"] = params
        return ConditionComparisonResult(
            data_id=data_id,
            method="pseudobulk",
            comparison="treated vs control",
            condition_key="condition",
            condition1="treated",
            condition2="control",
            sample_key="sample",
            n_samples_condition1=3,
            n_samples_condition2=3,
            global_n_significant=5,
            results_key="condition_results",
            statistics={"n_genes_tested": 120},
        )

    fake_module = SimpleNamespace(compare_conditions=fake_compare)
    monkeypatch.setitem(
        sys.modules, "chatspatial.tools.condition_comparison", fake_module
    )

    saved: dict[str, object] = {}

    async def fake_save_result(data_id: str, result_type: str, result):
        saved["data_id"] = data_id
        saved["result_type"] = result_type
        saved["result"] = result

    monkeypatch.setattr(data_manager, "save_result", fake_save_result)

    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
    )
    result = await compare_conditions(data_id="d2", params=params)

    assert isinstance(result, ConditionComparisonResult)
    assert calls["data_id"] == "d2"
    assert saved["data_id"] == "d2"
    assert saved["result_type"] == "condition_comparison_treated_vs_control"
    assert saved["result"] is result


@pytest.mark.integration
@pytest.mark.asyncio
async def test_analyze_cnv_saves_result_with_expected_key(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    async def fake_cnv(*, data_id, ctx, params):
        return CNVResult(
            data_id=data_id,
            method=params.method,
            reference_key=params.reference_key,
            reference_categories=params.reference_categories,
            n_chromosomes=22,
            n_genes_analyzed=100,
        )

    from chatspatial.tools.cnv_analysis import _build_cnv_key

    fake_module = SimpleNamespace(infer_cnv=fake_cnv, _build_cnv_key=_build_cnv_key)
    monkeypatch.setitem(sys.modules, "chatspatial.tools.cnv_analysis", fake_module)

    saved: dict[str, object] = {}

    async def fake_save_result(data_id: str, result_type: str, result):
        saved["data_id"] = data_id
        saved["result_type"] = result_type

    monkeypatch.setattr(data_manager, "save_result", fake_save_result)

    cnv_params = CNVParameters(
        reference_key="cell_type", reference_categories=["immune"]
    )
    result = await analyze_cnv("d3", params=cnv_params)
    assert isinstance(result, CNVResult)
    assert saved["data_id"] == "d3"
    assert saved["result_type"] == "cnv_infercnvpy_immune"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_velocity_and_trajectory_wrappers_save_expected_result_keys(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    async def fake_velocity(data_id, ctx, params):
        assert isinstance(params, RNAVelocityParameters)
        return RNAVelocityResult(
            data_id=data_id, velocity_computed=True, mode="stochastic"
        )

    async def fake_trajectory(data_id, ctx, params):
        assert isinstance(params, TrajectoryParameters)
        return TrajectoryResult(
            data_id=data_id,
            pseudotime_computed=True,
            velocity_computed=True,
            pseudotime_key="dpt_pseudotime",
            method="dpt",
            spatial_weight=0.0,
        )

    from chatspatial.tools.trajectory import _build_trajectory_key
    from chatspatial.tools.velocity import _build_velocity_key

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.velocity",
        SimpleNamespace(
            analyze_rna_velocity=fake_velocity, _build_velocity_key=_build_velocity_key
        ),
    )
    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.trajectory",
        SimpleNamespace(
            analyze_trajectory=fake_trajectory,
            _build_trajectory_key=_build_trajectory_key,
        ),
    )

    saved_calls: list[tuple[str, str]] = []

    async def fake_save_result(data_id: str, result_type: str, result):
        saved_calls.append((data_id, result_type))

    monkeypatch.setattr(data_manager, "save_result", fake_save_result)

    velocity_result = await analyze_velocity_data("d4")
    trajectory_result = await analyze_trajectory_data("d4")

    assert isinstance(velocity_result, RNAVelocityResult)
    assert isinstance(trajectory_result, TrajectoryResult)
    assert ("d4", "velocity_scvelo_stochastic") in saved_calls
    assert ("d4", "trajectory_cellrank_sw0_50") in saved_calls


@pytest.mark.integration
@pytest.mark.asyncio
async def test_preprocess_and_annotation_wrappers_materialize_default_params(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    async def fake_preprocess(data_id, ctx, params):
        assert isinstance(params, PreprocessingParameters)
        return PreprocessingResult(
            data_id=data_id,
            n_cells=12,
            n_genes=50,
            n_hvgs=20,
            clusters=0,
            qc_metrics={"ok": True},
        )

    async def fake_annotate(data_id, ctx, params):
        assert isinstance(params, AnnotationParameters)
        return AnnotationResult(
            data_id=data_id,
            method=params.method,
            output_key="cell_type_tangram",
            confidence_key="confidence_tangram",
            cell_types=["T", "B"],
            counts={"T": 1, "B": 1},
            confidence_scores={"T": 0.8, "B": 0.7},
        )

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.preprocessing",
        SimpleNamespace(preprocess_data=fake_preprocess),
    )
    from chatspatial.tools.annotation import _build_annotation_suffix

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.annotation",
        SimpleNamespace(
            annotate_cell_types=fake_annotate,
            _build_annotation_suffix=_build_annotation_suffix,
        ),
    )

    saved_calls: list[tuple[str, str]] = []

    async def fake_save_result(data_id: str, result_type: str, result):
        saved_calls.append((data_id, result_type))

    monkeypatch.setattr(data_manager, "save_result", fake_save_result)

    prep_result = await preprocess_data("d_default")
    ann_result = await annotate_cell_types("d_default")

    assert isinstance(prep_result, PreprocessingResult)
    assert isinstance(ann_result, AnnotationResult)
    assert ("d_default", "preprocessing") in saved_calls
    # Default tangram without reference → suffix is just "tangram"
    assert ("d_default", "annotation_tangram") in saved_calls


@pytest.mark.integration
@pytest.mark.asyncio
async def test_domains_and_spatial_genes_wrappers_materialize_default_params(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    async def fake_domains(data_id, ctx, params):
        assert isinstance(params, SpatialDomainParameters)
        return SpatialDomainResult(
            data_id=data_id,
            method=params.method,
            n_domains=2,
            domain_key="spatial_domains",
            domain_counts={"0": 5, "1": 7},
        )

    async def fake_spatial_genes(data_id, ctx, params):
        assert isinstance(params, SpatialVariableGenesParameters)
        return SpatialVariableGenesResult(
            data_id=data_id,
            method=params.method,
            n_genes_analyzed=20,
            n_significant_genes=4,
            spatial_genes=["gene_1", "gene_2"],
            results_key=f"spatial_genes_{params.method}",
        )

    from chatspatial.tools.spatial_domains import _build_domain_suffix

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.spatial_domains",
        SimpleNamespace(
            identify_spatial_domains=fake_domains,
            _build_domain_suffix=_build_domain_suffix,
        ),
    )
    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.spatial_genes",
        SimpleNamespace(identify_spatial_genes=fake_spatial_genes),
    )

    saved_calls: list[tuple[str, str]] = []

    async def fake_save_result(data_id: str, result_type: str, result):
        saved_calls.append((data_id, result_type))

    monkeypatch.setattr(data_manager, "save_result", fake_save_result)

    domain_result = await identify_spatial_domains("d_defaults")
    spatial_gene_result = await find_spatial_genes("d_defaults")

    assert isinstance(domain_result, SpatialDomainResult)
    assert isinstance(spatial_gene_result, SpatialVariableGenesResult)
    # Default spagcn with n_domains=7 → suffix is "spagcn_n7"
    assert ("d_defaults", "spatial_domains_spagcn_n7") in saved_calls
    assert ("d_defaults", "spatial_genes_flashs") in saved_calls


@pytest.mark.integration
@pytest.mark.asyncio
async def test_deconvolution_domains_comm_and_spatial_genes_save_expected_keys(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    async def fake_deconv(data_id, ctx, params):
        return DeconvolutionResult(
            data_id=data_id,
            method=params.method,
            dominant_type_key="dominant_cell_type",
            cell_types=["T", "B"],
            n_cell_types=2,
            proportions_key="cell_type_proportions",
        )

    async def fake_domains(data_id, ctx, params):
        return SpatialDomainResult(
            data_id=data_id,
            method=params.method,
            n_domains=3,
            domain_key="spatial_domains",
            domain_counts={"0": 12, "1": 11, "2": 10},
        )

    async def fake_comm(data_id, ctx, params):
        return CellCommunicationResult(
            data_id=data_id,
            method=params.method,
            species=params.species,
            database="consensus",
            analysis_type="cluster",
            n_lr_pairs=20,
            n_significant_pairs=5,
        )

    async def fake_spatial_genes(data_id, ctx, params):
        return SpatialVariableGenesResult(
            data_id=data_id,
            method=params.method,
            n_genes_analyzed=50,
            n_significant_genes=7,
            spatial_genes=["gene_1", "gene_2"],
            results_key=f"spatial_genes_{params.method}",
        )

    from chatspatial.tools.spatial_domains import _build_domain_suffix

    from chatspatial.tools.deconvolution import _build_deconvolution_key

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.deconvolution",
        SimpleNamespace(
            deconvolve_spatial_data=fake_deconv,
            _build_deconvolution_key=_build_deconvolution_key,
        ),
    )
    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.spatial_domains",
        SimpleNamespace(
            identify_spatial_domains=fake_domains,
            _build_domain_suffix=_build_domain_suffix,
        ),
    )
    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.cell_communication",
        SimpleNamespace(analyze_cell_communication=fake_comm),
    )
    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.spatial_genes",
        SimpleNamespace(identify_spatial_genes=fake_spatial_genes),
    )

    saved_calls: list[tuple[str, str]] = []

    async def fake_save_result(data_id: str, result_type: str, result):
        saved_calls.append((data_id, result_type))

    monkeypatch.setattr(data_manager, "save_result", fake_save_result)

    deconv_params = DeconvolutionParameters(
        method="flashdeconv",
        reference_data_id="ref_1",
        cell_type_key="cell_type",
    )
    domain_params = SpatialDomainParameters(method="leiden")
    comm_params = CellCommunicationParameters(
        method="liana",
        species="human",
        cell_type_key="cell_type",
    )
    spatial_gene_params = SpatialVariableGenesParameters(method="flashs")

    deconv_result = await deconvolve_data("d5", params=deconv_params)
    domain_result = await identify_spatial_domains("d5", params=domain_params)
    comm_result = await analyze_cell_communication("d5", params=comm_params)
    spatial_gene_result = await find_spatial_genes("d5", params=spatial_gene_params)

    assert isinstance(deconv_result, DeconvolutionResult)
    assert isinstance(domain_result, SpatialDomainResult)
    assert isinstance(comm_result, CellCommunicationResult)
    assert isinstance(spatial_gene_result, SpatialVariableGenesResult)
    assert ("d5", "deconvolution_flashdeconv_ref_1") in saved_calls
    # leiden with default resolution=0.5 → suffix "leiden_res0_50"
    assert ("d5", "spatial_domains_leiden_res0_50") in saved_calls
    assert ("d5", "cell_communication_liana") in saved_calls
    assert ("d5", "spatial_genes_flashs") in saved_calls
