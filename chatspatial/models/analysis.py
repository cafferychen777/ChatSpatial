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
    velocity_graph_key: Optional[str] = None  # Key for velocity graph in adata.uns
    mode: str  # RNA velocity computation mode
    visualization: Optional[Image] = None  # Image object

    class Config:
        arbitrary_types_allowed = True


class TrajectoryResult(BaseModel):
    """Result of trajectory analysis"""
    data_id: str
    pseudotime_computed: bool
    velocity_computed: bool
    pseudotime_key: str
    method: str  # Trajectory analysis method used
    spatial_weight: float  # Spatial information weight
    pseudotime_visualization: Optional[Image] = None  # Image object - now optional
    velocity_visualization: Optional[Image] = None  # Image object

    class Config:
        arbitrary_types_allowed = True


class IntegrationResult(BaseModel):
    """Result of sample integration"""
    data_id: str
    n_samples: int
    integration_method: str
    umap_visualization: Optional[Image] = None  # UMAP visualization image
    spatial_visualization: Optional[Image] = None  # Spatial coordinates visualization image

    class Config:
        arbitrary_types_allowed = True


class DeconvolutionResult(BaseModel):
    """Result of spatial deconvolution"""
    data_id: str
    method: str
    cell_types: List[str]
    n_cell_types: int
    proportions_key: str  # Key in adata.obsm where cell type proportions are stored
    visualization: Optional[Image] = None  # Visualization of cell type proportions
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
    statistics: Dict[str, Any]  # Statistics about the domain identification
    embeddings_key: Optional[str] = None  # Key in adata.obsm where embeddings are stored

    class Config:
        arbitrary_types_allowed = True




class SpatialVariableGenesResult(BaseModel):
    """Result of GASTON spatial variable genes identification"""
    data_id: str
    preprocessing_method: str
    n_components: int
    n_epochs_trained: int
    final_loss: float

    # Model architecture info
    spatial_hidden_layers: List[int]
    expression_hidden_layers: List[int]

    # Spatial domains and patterns
    n_spatial_domains: int
    spatial_domains_key: str  # Key in adata.obs where spatial domain assignments are stored
    isodepth_key: str  # Key in adata.obs where isodepth values are stored

    # Gene classification results
    continuous_gradient_genes: Dict[str, List[int]]  # Gene -> list of domains with continuous gradients
    discontinuous_genes: Dict[str, List[int]]  # Gene -> list of domain boundaries with discontinuities
    n_continuous_genes: int
    n_discontinuous_genes: int

    # Model outputs stored in adata
    model_predictions_key: str  # Key in adata.obsm where model predictions are stored
    spatial_embedding_key: str  # Key in adata.obsm where spatial embeddings are stored

    # Visualizations (use visualize_data tool for GASTON visualizations)
    isodepth_map_visualization: Optional[Image] = None  # Deprecated: use visualize_data with plot_type="gaston_isodepth"
    spatial_domains_visualization: Optional[Image] = None  # Deprecated: use visualize_data with plot_type="gaston_domains"
    top_genes_visualization: Optional[Image] = None  # Deprecated: use visualize_data with plot_type="gaston_genes"

    # Statistics and metrics
    model_performance: Dict[str, Any]  # Model performance metrics
    spatial_autocorrelation: Dict[str, float]  # Spatial autocorrelation metrics

    # File paths for saved model
    model_checkpoint_path: Optional[str] = None  # Path to saved model checkpoint

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



    # LIANA+ specific results
    liana_results_key: Optional[str] = None  # Key in adata.uns for LIANA cluster results
    liana_spatial_results_key: Optional[str] = None  # Key in adata.uns for LIANA spatial results
    liana_spatial_scores_key: Optional[str] = None  # Key in adata.obsm for spatial scores
    analysis_type: Optional[str] = None  # Type of LIANA analysis: 'cluster' or 'spatial'

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