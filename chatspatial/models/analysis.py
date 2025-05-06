"""
Analysis result models for spatial transcriptomics data.
"""

from typing import Dict, List, Optional, Any
from pydantic import BaseModel
from mcp.server.fastmcp.utilities.types import Image


class PreprocessingResult(BaseModel):
    """Result of data preprocessing"""
    data_id: str
    n_cells: int
    n_genes: int
    n_hvgs: int
    clusters: int
    qc_metrics: Optional[Dict[str, Any]] = None


class DifferentialExpressionResult(BaseModel):
    """Result of differential expression analysis"""
    data_id: str
    comparison: str
    n_genes: int
    top_genes: List[str]
    statistics: Dict[str, Any]


class AnnotationResult(BaseModel):
    """Result of cell type annotation"""
    data_id: str
    method: str
    cell_types: List[str]
    counts: Dict[str, int]
    confidence_scores: Optional[Dict[str, float]] = None


class SpatialAnalysisResult(BaseModel):
    """Result of spatial analysis"""
    data_id: str
    analysis_type: str
    statistics: Optional[Dict[str, Any]] = None
    result_image: Optional[str] = None  # Base64 encoded image


class RNAVelocityResult(BaseModel):
    """Result of RNA velocity analysis"""
    data_id: str
    velocity_computed: bool
    visualization: Optional[Image] = None  # Image object

    class Config:
        arbitrary_types_allowed = True


class TrajectoryResult(BaseModel):
    """Result of trajectory analysis"""
    data_id: str
    pseudotime_computed: bool
    velocity_computed: bool
    pseudotime_key: str
    pseudotime_visualization: Image  # Image object
    velocity_visualization: Optional[Image] = None  # Image object

    class Config:
        arbitrary_types_allowed = True


class IntegrationResult(BaseModel):
    """Result of sample integration"""
    data_id: str
    n_samples: int
    integration_method: str
    umap_visualization: Image  # UMAP visualization image
    spatial_visualization: Image  # Spatial coordinates visualization image

    class Config:
        arbitrary_types_allowed = True


class DeconvolutionResult(BaseModel):
    """Result of spatial deconvolution"""
    data_id: str
    method: str
    cell_types: List[str]
    n_cell_types: int
    proportions_key: str  # Key in adata.obsm where cell type proportions are stored
    visualization: Image  # Visualization of cell type proportions
    statistics: Dict[str, Any]  # Statistics about the deconvolution results

    class Config:
        arbitrary_types_allowed = True