"""Unit tests for Pydantic parameter model validation contracts."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from chatspatial.models.analysis import (
    CellTypeComparisonResult,
    ConditionComparisonResult,
    DEGene,
    EmbeddingResult,
    EnrichmentResult,
    PreprocessingResult,
    SpatialRegistrationResult,
)
from chatspatial.models.data import (
    AnnotationParameters,
    CellCommunicationParameters,
    CNVParameters,
    ColumnInfo,
    DeconvolutionParameters,
    DifferentialExpressionParameters,
    EnrichmentParameters,
    IntegrationParameters,
    PreprocessingParameters,
    RegistrationParameters,
    RNAVelocityParameters,
    SpatialDataset,
    SpatialDomainParameters,
    SpatialStatisticsParameters,
    SpatialVariableGenesParameters,
    TrajectoryParameters,
    VisualizationParameters,
)


def _de_gene() -> DEGene:
    return DEGene(gene="G", log2fc=1.0, pvalue=0.01, padj=0.02)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_public_result_models_reject_non_finite_numbers(value: float):
    with pytest.raises(ValidationError, match="finite number"):
        DEGene(gene="G", log2fc=value, pvalue=0.01, padj=0.02)

    with pytest.raises(ValidationError, match="finite number"):
        ColumnInfo(
            name="score",
            dtype="numerical",
            n_unique=1,
            range=[value, value],
        )

    with pytest.raises(ValidationError, match="finite number"):
        PreprocessingResult(
            data_id="d1",
            n_cells=1,
            n_genes=1,
            n_hvgs=1,
            clusters=0,
            qc_metrics={"median_counts": value},
        )


@pytest.mark.parametrize(
    ("model_cls", "kwargs"),
    [
        (
            EmbeddingResult,
            {"data_id": "d1", "computed": [], "skipped": [], "unknown": True},
        ),
        (
            SpatialRegistrationResult,
            {
                "method": "paste",
                "source_id": "source",
                "target_id": "target",
                "n_source_spots": 1,
                "n_target_spots": 1,
                "spatial_key_registered": "spatial_registered",
                "unknown": True,
            },
        ),
        (
            SpatialDataset,
            {
                "id": "d1",
                "name": "demo",
                "data_type": "generic",
                "unknown": True,
            },
        ),
    ],
)
def test_public_result_models_reject_unknown_fields(model_cls, kwargs):
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        model_cls(**kwargs)


def test_excluded_enrichment_details_preserve_infinite_odds_ratios():
    result = EnrichmentResult(
        method="ora",
        n_gene_sets=1,
        n_significant=1,
        top_gene_sets=["perfect_separation"],
        top_depleted_sets=[],
        enrichment_scores={"perfect_separation": float("inf")},
    )

    assert result.enrichment_scores["perfect_separation"] == float("inf")
    assert "enrichment_scores" not in result.model_dump()
    assert "Infinity" not in result.model_dump_json()


def test_nested_analysis_values_do_not_publish_warning_fields():
    cell_type_result = CellTypeComparisonResult(
        cell_type="T",
        n_cells_condition1=10,
        n_cells_condition2=10,
        n_samples_condition1=2,
        n_samples_condition2=2,
        n_significant_genes=1,
        top_upregulated=[_de_gene()],
        top_downregulated=[],
    )

    dumped = cell_type_result.model_dump()
    schema = ConditionComparisonResult.model_json_schema()

    assert "warnings" not in dumped
    assert "warnings" not in dumped["top_upregulated"][0]
    assert "warnings" not in schema["$defs"]["CellTypeComparisonResult"]["properties"]
    assert "warnings" not in schema["$defs"]["DEGene"]["properties"]


def test_top_level_analysis_result_only_serializes_nonempty_warnings():
    result = ConditionComparisonResult(
        data_id="d1",
        method="pseudobulk",
        comparison="treated vs control",
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
        n_samples_condition1=2,
        n_samples_condition2=2,
        global_n_significant=1,
        global_top_upregulated=[_de_gene()],
        results_key="condition_comparison",
        statistics={},
    )

    assert "warnings" not in result.model_dump()
    assert result.model_copy(update={"warnings": ["low sample count"]}).model_dump()[
        "warnings"
    ] == ["low sample count"]


def test_differential_expression_requires_group_key():
    with pytest.raises(ValidationError):
        DifferentialExpressionParameters()


def test_differential_expression_rejects_invalid_method():
    with pytest.raises(ValidationError):
        DifferentialExpressionParameters(group_key="cluster", method="unknown")


def test_preprocessing_rejects_invalid_normalization():
    with pytest.raises(ValidationError):
        PreprocessingParameters(normalization="invalid")


def test_preprocessing_rejects_out_of_range_scrublet_rate():
    with pytest.raises(ValidationError):
        PreprocessingParameters(scrublet_expected_doublet_rate=0.9)


def test_visualization_preprocesses_string_input():
    params = VisualizationParameters.model_validate("gene:CCL21")
    assert params.plot_type == "feature"
    assert params.feature == "CCL21"


def test_visualization_preprocess_params_handles_none_alias_and_passthrough():
    assert VisualizationParameters.preprocess_params(None) == {}
    assert VisualizationParameters.preprocess_params("CXCL12") == {
        "feature": "CXCL12",
        "plot_type": "feature",
    }
    assert VisualizationParameters.preprocess_params({"features": ["g1", "g2"]}) == {
        "feature": ["g1", "g2"]
    }
    sentinel = object()
    assert VisualizationParameters.preprocess_params(sentinel) is sentinel


def test_visualization_statistics_requires_subtype():
    with pytest.raises(
        ValidationError, match="subtype is required when plot_type='statistics'"
    ):
        VisualizationParameters(plot_type="statistics")


@pytest.mark.parametrize(
    ("plot_type", "expected_subtype"),
    [
        ("cnv", "heatmap"),
        ("velocity", "stream"),
        ("enrichment", "barplot"),
        ("trajectory", "pseudotime"),
    ],
)
def test_visualization_defaults_subtype_by_plot_type(
    plot_type: str, expected_subtype: str
):
    params = VisualizationParameters(plot_type=plot_type)
    assert params.subtype == expected_subtype


def test_visualization_communication_subtype_defaults_at_runtime():
    params = VisualizationParameters(plot_type="communication")
    assert params.subtype is None


def test_annotation_rejects_invalid_method_literal():
    with pytest.raises(ValidationError):
        AnnotationParameters(method="not_supported")


def test_annotation_sctype_remote_options_default_to_false():
    params = AnnotationParameters(method="sctype")
    assert params.sctype_allow_remote is False
    assert params.sctype_allow_runtime_r_install is False


def test_preprocessing_uses_canonical_scrublet_field_name():
    params = PreprocessingParameters(use_scrublet=True)
    assert params.use_scrublet is True
    assert "use_scrublet" in PreprocessingParameters.model_fields


def test_spatial_domain_includes_banksy_literal():
    params = SpatialDomainParameters(method="banksy")
    assert params.method == "banksy"


def test_spatial_domain_includes_aestetik_literal():
    params = SpatialDomainParameters(method="aestetik")
    assert params.method == "aestetik"
    assert params.aestetik_transcriptomics_key == "X_pca"
    assert params.aestetik_morphology_key == "X_pca_morphology"


def test_spatial_domain_rejects_even_aestetik_window_size():
    with pytest.raises(ValidationError):
        SpatialDomainParameters(method="aestetik", aestetik_window_size=4)


def test_spatial_domain_rejects_aestetik_morphology_weight_above_total_weight():
    with pytest.raises(ValidationError):
        SpatialDomainParameters(method="aestetik", aestetik_morphology_weight=3.1)


def test_spatial_variable_genes_default_is_flashs():
    params = SpatialVariableGenesParameters()
    assert params.method == "flashs"


def test_enrichment_default_is_spatial_enrichmap():
    params = EnrichmentParameters(species="human")
    assert params.method == "spatial_enrichmap"


def test_cell_communication_uses_canonical_cellchat_literal():
    params = CellCommunicationParameters(
        method="cellchat_r",
        species="mouse",
        cell_type_key="cell_type",
    )
    assert params.method == "cellchat_r"

    with pytest.raises(ValidationError):
        CellCommunicationParameters(
            method="cellchat",
            species="human",
            cell_type_key="cell_type",
        )


def test_cellphonedb_debug_seed_is_explicit_and_nonnegative():
    params = CellCommunicationParameters(
        method="cellphonedb",
        species="human",
        cell_type_key="cell_type",
    )
    assert params.cellphonedb_debug_seed is None

    with pytest.raises(ValidationError):
        CellCommunicationParameters(
            method="cellphonedb",
            species="human",
            cell_type_key="cell_type",
            cellphonedb_debug_seed=-1,
        )


def test_numbat_genome_matches_supported_backend_references():
    for genome in ["hg38", "hg19", "mm10"]:
        params = CNVParameters(
            method="numbat",
            reference_key="cell_type",
            reference_categories=["normal"],
            numbat_genome=genome,
        )
        assert params.numbat_genome == genome

    with pytest.raises(ValidationError):
        CNVParameters(
            method="numbat",
            reference_key="cell_type",
            reference_categories=["normal"],
            numbat_genome="mm39",
        )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("velovi_n_hidden", 0),
        ("velovi_n_latent", 0),
        ("velovi_n_layers", 0),
        ("velovi_n_epochs", 0),
        ("velovi_dropout_rate", -0.1),
        ("velovi_dropout_rate", 1.0),
        ("velovi_learning_rate", 0.0),
    ],
)
def test_rna_velocity_rejects_invalid_velovi_hyperparameters(field: str, value: float):
    with pytest.raises(ValidationError):
        RNAVelocityParameters.model_validate({field: value})


def test_spatial_statistics_accepts_extended_analysis_literals():
    for analysis_type in [
        "local_join_count",
        "network_properties",
        "spatial_centrality",
    ]:
        params = SpatialStatisticsParameters(
            analysis_type=analysis_type, cluster_key="cluster"
        )
        assert params.analysis_type == analysis_type


@pytest.mark.parametrize(
    ("plot_type", "subtype"),
    [
        ("deconvolution", "imputation"),
        ("trajectory", "circular"),
        ("trajectory", "fate_heatmap"),
        ("trajectory", "palantir"),
        ("velocity", "proportions"),
        ("velocity", "heatmap"),
        ("statistics", "getis_ord"),
        ("statistics", "centrality"),
        ("integration", "highlight"),
    ],
)
def test_visualization_accepts_documented_subtypes(plot_type: str, subtype: str):
    params = VisualizationParameters(plot_type=plot_type, subtype=subtype)
    assert params.subtype == subtype


def _collect_schema_array_issues(
    schema: object, path: tuple[str, ...] = ()
) -> list[str]:
    issues: list[str] = []
    if isinstance(schema, dict):
        if schema.get("type") == "array":
            if "items" not in schema:
                issues.append(f"{'/'.join(path)} is an array without items")
            if "prefixItems" in schema:
                issues.append(f"{'/'.join(path)} uses prefixItems")
        for key, value in schema.items():
            issues.extend(_collect_schema_array_issues(value, (*path, str(key))))
    elif isinstance(schema, list):
        for index, value in enumerate(schema):
            issues.extend(_collect_schema_array_issues(value, (*path, str(index))))
    return issues


@pytest.mark.parametrize(
    "model_cls",
    [
        ColumnInfo,
        SpatialDataset,
        VisualizationParameters,
        SpatialStatisticsParameters,
        TrajectoryParameters,
        RegistrationParameters,
    ],
)
def test_json_bound_pair_models_emit_anthropic_compatible_array_schemas(model_cls):
    issues = _collect_schema_array_issues(model_cls.model_json_schema())
    assert not issues


def test_pair_parameters_accept_tuple_like_inputs_as_json_arrays():
    assert VisualizationParameters(panel_layout=(2, 3)).panel_layout == [2, 3]
    assert VisualizationParameters(figure_size=(6, 4)).figure_size == [6, 4]
    assert VisualizationParameters(lr_pairs=[("L1", "R1")]).lr_pairs == [["L1", "R1"]]
    assert SpatialStatisticsParameters(gene_pairs=[("G1", "G2")]).gene_pairs == [
        ["G1", "G2"]
    ]
    assert TrajectoryParameters(
        cellrank_kernel_weights=(2.0, 1.0)
    ).cellrank_kernel_weights == pytest.approx([2 / 3, 1 / 3])
    assert RegistrationParameters(stalign_image_size=(256, 128)).stalign_image_size == [
        256,
        128,
    ]
    assert ColumnInfo(
        name="score", dtype="numerical", n_unique=3, range=(0.0, 1.0)
    ).range == [
        0.0,
        1.0,
    ]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"dpi": 0},
        {"alpha": -0.1},
        {"alpha": 1.1},
        {"figure_size": [0, 5]},
        {"panel_layout": [2, -1]},
        {"spot_size": 0},
        {"colorbar_size": "wide"},
        {"colorbar_size": "0%"},
    ],
)
def test_visualization_parameters_reject_invalid_rendering_bounds(kwargs):
    with pytest.raises(ValidationError):
        VisualizationParameters(**kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"n_epochs": 0},
        {"scvi_n_hidden": 0},
        {"scvi_n_latent": -1},
        {"scvi_n_layers": 0},
        {"scvi_dropout_rate": 1.0},
    ],
)
def test_integration_parameters_reject_invalid_training_bounds(kwargs):
    with pytest.raises(ValidationError):
        IntegrationParameters(**kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"destvi_n_hidden": 0},
        {"destvi_n_latent": 0},
        {"destvi_n_layers": 0},
        {"destvi_dropout_rate": 1.0},
        {"destvi_learning_rate": 0},
        {"stereoscope_n_epochs": 0},
        {"stereoscope_learning_rate": -0.1},
        {"stereoscope_batch_size": 0},
    ],
)
def test_deconvolution_parameters_reject_invalid_training_bounds(kwargs):
    with pytest.raises(ValidationError):
        DeconvolutionParameters(cell_type_key="cell_type", **kwargs)


def test_spatial_domain_timeout_must_be_positive():
    with pytest.raises(ValidationError):
        SpatialDomainParameters(timeout=0)


@pytest.mark.parametrize(
    ("model_cls", "kwargs"),
    [
        (VisualizationParameters, {"panel_layout": [2]}),
        (VisualizationParameters, {"figure_size": [6, 4, 2]}),
        (VisualizationParameters, {"lr_pairs": [["L1"]]}),
        (SpatialStatisticsParameters, {"gene_pairs": [["G1", "G2", "G3"]]}),
        (TrajectoryParameters, {"cellrank_kernel_weights": [1.0]}),
        (RegistrationParameters, {"stalign_image_size": [128, 128, 64]}),
        (
            ColumnInfo,
            {"name": "score", "dtype": "numerical", "n_unique": 3, "range": [0.0]},
        ),
    ],
)
def test_pair_parameters_reject_invalid_lengths(model_cls, kwargs):
    with pytest.raises(ValidationError):
        model_cls(**kwargs)
