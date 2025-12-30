"""
Data models for spatial transcriptomics analysis.
"""

# Import result models from analysis module
from .analysis import (
    AnnotationResult,
    CellCommunicationResult,
    CNVResult,
    DeconvolutionResult,
    DifferentialExpressionResult,
    EnrichmentResult,
    IntegrationResult,
    PreprocessingResult,
    RNAVelocityResult,
    SpatialDomainResult,
    SpatialStatisticsResult,
    SpatialVariableGenesResult,
    TrajectoryResult,
)

# Import parameter models from data module
from .data import (
    AnnotationParameters,
    CellCommunicationParameters,
    CNVParameters,
    ColumnInfo,
    DeconvolutionParameters,
    EnrichmentParameters,
    IntegrationParameters,
    PreprocessingParameters,
    RNAVelocityParameters,
    SpatialDataset,
    SpatialDomainParameters,
    SpatialStatisticsParameters,
    SpatialVariableGenesParameters,
    TrajectoryParameters,
    VisualizationParameters,
)

__all__ = [
    # Result models
    "AnnotationResult",
    "CellCommunicationResult",
    "CNVResult",
    "DeconvolutionResult",
    "DifferentialExpressionResult",
    "EnrichmentResult",
    "IntegrationResult",
    "PreprocessingResult",
    "RNAVelocityResult",
    "SpatialDomainResult",
    "SpatialStatisticsResult",
    "SpatialVariableGenesResult",
    "TrajectoryResult",
    # Parameter models
    "AnnotationParameters",
    "CellCommunicationParameters",
    "CNVParameters",
    "ColumnInfo",
    "DeconvolutionParameters",
    "EnrichmentParameters",
    "IntegrationParameters",
    "PreprocessingParameters",
    "RNAVelocityParameters",
    "SpatialDataset",
    "SpatialDomainParameters",
    "SpatialStatisticsParameters",
    "SpatialVariableGenesParameters",
    "TrajectoryParameters",
    "VisualizationParameters",
]
