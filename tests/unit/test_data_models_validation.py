"""Unit tests for Pydantic parameter model validation contracts."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from chatspatial.models.data import (
    AnnotationParameters,
    CellCommunicationParameters,
    ColumnInfo,
    DifferentialExpressionParameters,
    EnrichmentParameters,
    PreprocessingParameters,
    RegistrationParameters,
    SpatialDataset,
    SpatialDomainParameters,
    SpatialStatisticsParameters,
    SpatialVariableGenesParameters,
    TrajectoryParameters,
    VisualizationParameters,
)


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
