"""
Analysis result models for spatial transcriptomics data.
"""

from typing import Any, Dict, List, Optional

from pydantic import BaseModel

try:
    from mcp.server.fastmcp.utilities.types import Image
except ImportError:
    # Fallback for when MCP is not available
    Image = Any


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
    """Result of cell type annotation

    Attributes:
        data_id: Dataset identifier
        method: Annotation method used
        cell_types: List of unique cell types identified
        counts: Number of cells per cell type
        confidence_scores: Confidence scores per cell type (when available).
                          Empty dict or None indicates no confidence data available.
                          Only contains real statistical measures, never arbitrary values.
        tangram_mapping_score: For Tangram method - overall mapping quality score
    """

    data_id: str
    method: str
    cell_types: List[str]
    counts: Dict[str, int]
    confidence_scores: Optional[Dict[str, float]] = None
    tangram_mapping_score: Optional[float] = None  # For Tangram method - mapping score

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
    method_fallback_used: Optional[bool] = (
        False  # Whether fallback to simpler method was used
    )

    class Config:
        arbitrary_types_allowed = True


class IntegrationResult(BaseModel):
    """Result of sample integration"""

    data_id: str
    n_samples: int
    integration_method: str

    class Config:
        arbitrary_types_allowed = True


class DeconvolutionResult(BaseModel):
    """Result of spatial deconvolution"""

    data_id: str
    method: str
    cell_types: List[str]
    n_cell_types: int
    proportions_key: str  # Key in adata.obsm where cell type proportions are stored
    statistics: Dict[str, Any]  # Statistics about the deconvolution results

    class Config:
        arbitrary_types_allowed = True


class SpatialDomainResult(BaseModel):
    """Result of spatial domain identification"""

    data_id: str
    method: str
    n_domains: int
    domain_key: str  # Key in adata.obs where domain labels are stored
    domain_counts: Dict[str, int]  # Number of spots in each domain
    refined_domain_key: Optional[str] = (
        None  # Key for refined domains if refinement was applied
    )
    statistics: Dict[str, Any]  # Statistics about the domain identification
    embeddings_key: Optional[str] = (
        None  # Key in adata.obsm where embeddings are stored
    )

    class Config:
        arbitrary_types_allowed = True


class SpatialVariableGenesResult(BaseModel):
    """Result of spatial variable genes identification"""

    data_id: str
    method: str  # Method used for analysis

    # Common results for all methods
    n_genes_analyzed: int  # Total number of genes analyzed
    n_significant_genes: int  # Number of significant spatial genes found
    spatial_genes: List[str]  # List of significant spatial variable gene names

    # Statistical results (available for all methods)
    gene_statistics: Dict[str, float]  # Gene name -> primary statistic value
    p_values: Dict[str, float]  # Gene name -> p-value
    q_values: Dict[str, float]  # Gene name -> FDR-corrected p-value

    # Storage keys for results in adata
    results_key: str  # Base key for storing results in adata

    # Method-specific results (optional, only populated for respective methods)
    gaston_results: Optional[Dict[str, Any]] = None  # GASTON-specific results
    spatialde_results: Optional[Dict[str, Any]] = None  # SpatialDE-specific results
    sparkx_results: Optional[Dict[str, Any]] = None  # SPARK-X specific results
    somde_results: Optional[Dict[str, Any]] = None  # SOMDE-specific results

    # Visualization hints (optional)
    isodepth_visualization: Optional[Dict[str, Any]] = None  # For GASTON isodepth plots
    spatial_domains_visualization: Optional[Dict[str, Any]] = (
        None  # For spatial domain plots
    )
    top_genes_visualization: Optional[Dict[str, Any]] = (
        None  # For top genes visualization
    )

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
    global_results_key: Optional[str] = (
        None  # Key in adata.uns where global results are stored
    )
    top_lr_pairs: List[str]  # List of top significant LR pairs

    # Local analysis results (if performed)
    local_analysis_performed: bool = False
    local_results_key: Optional[str] = (
        None  # Key in adata.uns where local results are stored
    )
    communication_matrices_key: Optional[str] = (
        None  # Key in adata.obsp where communication matrices are stored
    )

    # LIANA+ specific results
    liana_results_key: Optional[str] = (
        None  # Key in adata.uns for LIANA cluster results
    )
    liana_spatial_results_key: Optional[str] = (
        None  # Key in adata.uns for LIANA spatial results
    )
    liana_spatial_scores_key: Optional[str] = (
        None  # Key in adata.obsm for spatial scores
    )
    analysis_type: Optional[str] = (
        None  # Type of LIANA analysis: 'cluster' or 'spatial'
    )

    # Communication patterns (if identified)
    patterns_identified: bool = False
    n_patterns: Optional[int] = None
    patterns_key: Optional[str] = (
        None  # Key in adata.obs where communication patterns are stored
    )

    # Statistics
    statistics: Dict[str, Any]  # General statistics about the communication analysis

    class Config:
        arbitrary_types_allowed = True


class EnrichmentResult(BaseModel):
    """Result from gene set enrichment analysis"""

    # Basic information
    method: str  # Method used (pathway_gsea, pathway_ora, pathway_enrichr, pathway_ssgsea, spatial_enrichmap)
    n_gene_sets: int  # Number of gene sets analyzed
    n_significant: int  # Number of significant gene sets

    # Enrichment scores and statistics
    enrichment_scores: Dict[str, float]  # Enrichment scores for each gene set
    pvalues: Dict[str, float]  # Raw p-values
    adjusted_pvalues: Dict[str, float]  # Adjusted p-values
    gene_set_statistics: Dict[str, Dict[str, Any]]  # Additional statistics per gene set

    # Spatial metrics (for enrichmap)
    spatial_metrics: Optional[Dict[str, Any]] = None  # Spatial autocorrelation metrics
    spatial_scores_key: Optional[str] = (
        None  # Key in adata.obsm for spatial enrichment scores
    )

    # Gene set information
    gene_sets_used: Dict[str, List[str]]  # Gene sets that were actually used
    genes_found: Dict[str, List[str]]  # Genes found in data for each gene set

    # Top results
    top_gene_sets: List[str]  # Top enriched gene sets
    top_depleted_sets: List[str]  # Top depleted gene sets

    # Additional metadata
    parameters_used: Dict[str, Any]  # Parameters used for analysis
    computation_time: float  # Time taken for computation

    class Config:
        arbitrary_types_allowed = True


class SpatialStatisticsIntegrationResult(BaseModel):
    """Result of spatial statistics integration and batch analysis"""

    data_id: str
    analysis_name: str
    timestamp: str
    integration_key: str  # Key prefix used for storing integrated results

    # Dataset information
    n_observations: int
    n_variables: int

    # Analysis methods used
    methods_used: List[
        str
    ]  # e.g., ['local_moran', 'getis_ord', 'neighborhood_enrichment']
    total_features_analyzed: int

    # Integration metadata
    local_statistics: Dict[str, Any]  # Information about local spatial statistics
    global_statistics: Dict[str, Any]  # Information about global spatial statistics
    result_keys: Dict[str, str]  # Mapping of method names to their result keys in adata

    # Performance metrics (if from batch analysis)
    analysis_times: Optional[Dict[str, float]] = None  # Time taken for each analysis
    total_analysis_time: Optional[float] = None  # Total time for batch analysis

    # Visualization readiness
    visualization_features: Optional[Dict[str, List[str]]] = (
        None  # Available features for visualization
    )

    class Config:
        arbitrary_types_allowed = True


class SpatialVisualizationData(BaseModel):
    """Formatted spatial statistics data ready for visualization"""

    result_key: str
    analysis_metadata: Dict[str, Any]  # Integration metadata

    # Spatial information
    spatial_coordinates: List[List[float]]  # x, y coordinates
    n_observations: int

    # Feature data organized by analysis method
    local_moran_features: Optional[Dict[str, Dict[str, Any]]] = (
        None  # Gene -> feature data
    )
    getis_ord_features: Optional[Dict[str, Dict[str, Any]]] = (
        None  # Gene -> feature data
    )
    bivariate_moran_features: Optional[Dict[str, Dict[str, Any]]] = (
        None  # Pair -> feature data
    )
    neighborhood_enrichment: Optional[Dict[str, Any]] = None  # Enrichment matrix data

    # Visualization configuration
    color_maps: Dict[str, str]  # Method -> colormap name
    data_types: Dict[str, str]  # Method -> data type (continuous, categorical, matrix)
    plot_configs: Dict[str, Dict[str, Any]]  # Method -> plot configuration

    # Export format
    output_format: str  # 'dict', 'dataframe', 'plotly'

    class Config:
        arbitrary_types_allowed = True


class SpatialAnalysisSummary(BaseModel):
    """Comprehensive summary of spatial analysis results"""

    result_key: str
    analysis_name: str
    timestamp: str

    # Overview information
    dataset_info: Dict[str, int]  # n_observations, n_variables
    methods_used: List[str]
    total_features_analyzed: int

    # Method-specific summaries
    method_summaries: Dict[str, Dict[str, Any]]  # Method -> summary statistics

    # Significant findings
    significant_findings: Dict[str, Dict[str, Any]]  # Method -> significant results

    # Interpretation guide
    interpretation_guide: Optional[Dict[str, str]] = (
        None  # Method -> interpretation text
    )

    # Output format
    output_format: str  # 'dict', 'markdown', 'html', 'json'
    formatted_content: Optional[str] = None  # Formatted content if not dict

    class Config:
        arbitrary_types_allowed = True
