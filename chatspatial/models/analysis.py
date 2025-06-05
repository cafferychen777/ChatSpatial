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
    tangram_mapping_score: Optional[float] = None  # For Tangram method - mapping score
    visualization: Optional[Image] = None  # For Tangram method - visualization of cell type mapping

    class Config:
        arbitrary_types_allowed = True


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
    visualization_params: Optional[Dict[str, Any]] = None  # Parameters for visualizing the results

    class Config:
        arbitrary_types_allowed = True


class SpatialDomainResult(BaseModel):
    """Result of spatial domain identification"""
    data_id: str
    method: str
    n_domains: int
    domain_key: str  # Key in adata.obs where domain labels are stored
    domain_counts: Dict[str, int]  # Number of spots in each domain
    refined_domain_key: Optional[str] = None  # Key for refined domains if refinement was applied
    visualization: Optional[Image] = None  # Visualization of spatial domains
    statistics: Dict[str, Any]  # Statistics about the domain identification
    embeddings_key: Optional[str] = None  # Key in adata.obsm where embeddings are stored (for STAGATE)

    class Config:
        arbitrary_types_allowed = True


class SpatialVariableGenesResult(BaseModel):
    """Result of spatial variable genes identification"""
    data_id: str
    method: str
    n_significant_genes: int
    n_tested_genes: int
    significance_threshold: float
    top_genes: List[str]  # List of top spatially variable genes
    results_key: str  # Key in adata.var where results are stored

    # SpatialDE specific results
    spatialde_results: Optional[Dict[str, Any]] = None  # Full SpatialDE results table

    # AEH results (if performed)
    aeh_performed: bool = False
    aeh_patterns_key: Optional[str] = None  # Key in adata.obsm where spatial patterns are stored
    aeh_membership_key: Optional[str] = None  # Key in adata.var where pattern membership is stored
    n_patterns: Optional[int] = None

    # Visualization
    visualization: Optional[Image] = None  # Visualization of top spatial genes
    patterns_visualization: Optional[Image] = None  # Visualization of spatial patterns (if AEH performed)

    # Statistics
    statistics: Dict[str, Any]  # General statistics about the analysis

    class Config:
        arbitrary_types_allowed = True


class CellCommunicationResult(BaseModel):
    """Result of cell-cell communication analysis"""
    data_id: str
    method: str
    species: str
    database: str
    n_lr_pairs: int  # Total number of LR pairs tested
    n_significant_pairs: int  # Number of significant LR pairs

    # Global analysis results
    global_results_key: Optional[str] = None  # Key in adata.uns where global results are stored
    top_lr_pairs: List[str]  # List of top significant LR pairs

    # Local analysis results (if performed)
    local_analysis_performed: bool = False
    local_results_key: Optional[str] = None  # Key in adata.uns where local results are stored
    communication_matrices_key: Optional[str] = None  # Key in adata.obsp where communication matrices are stored

    # COMMOT specific results
    commot_sender_key: Optional[str] = None  # Key in adata.obsm for sender signals
    commot_receiver_key: Optional[str] = None  # Key in adata.obsm for receiver signals

    # SpatialDM specific results
    spatialdm_selected_spots_key: Optional[str] = None  # Key in adata.uns for selected spots
    spatialdm_weight_matrix_key: Optional[str] = None  # Key in adata.obsp for weight matrix

    # Communication patterns (if identified)
    patterns_identified: bool = False
    n_patterns: Optional[int] = None
    patterns_key: Optional[str] = None  # Key in adata.obs where communication patterns are stored

    # Visualization
    visualization: Optional[Image] = None  # Visualization of top LR pairs
    network_visualization: Optional[Image] = None  # Communication network visualization

    # Statistics
    statistics: Dict[str, Any]  # General statistics about the communication analysis

    class Config:
        arbitrary_types_allowed = True