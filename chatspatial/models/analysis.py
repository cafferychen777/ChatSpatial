"""
Analysis result models for spatial transcriptomics data.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional

from pydantic import BaseModel, ConfigDict, Field

if TYPE_CHECKING:
    from mcp.types import ImageContent
else:
    try:
        from mcp.types import ImageContent
    except ImportError:
        # Fallback for when MCP is not available
        ImageContent = Any  # type: ignore[misc,assignment]


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
        output_key: Column name in adata.obs where cell types are stored (e.g., "cell_type_tangram")
        confidence_key: Column name in adata.obs where confidence scores are stored (e.g., "confidence_tangram")
        cell_types: List of unique cell types identified
        counts: Number of cells per cell type
        confidence_scores: Confidence scores per cell type (when available).
                          Empty dict or None indicates no confidence data available.
                          Only contains real statistical measures, never arbitrary values.
        tangram_mapping_score: For Tangram method - overall mapping quality score
    """

    data_id: str
    method: str
    output_key: str  # Column name where cell types are stored
    confidence_key: Optional[str] = (
        None  # Column name where confidence scores are stored
    )
    cell_types: List[str]
    counts: Dict[str, int]
    confidence_scores: Optional[Dict[str, float]] = None
    tangram_mapping_score: Optional[float] = None  # For Tangram method - mapping score

    model_config = ConfigDict(arbitrary_types_allowed=True)


class SpatialStatisticsResult(BaseModel):
    """Result of spatial analysis

    Note: Visualization is handled separately via the visualize_data tool.
    This model only contains statistical results and metadata.
    """

    data_id: str
    analysis_type: str
    statistics: Optional[Dict[str, Any]] = None


class RNAVelocityResult(BaseModel):
    """Result of RNA velocity analysis"""

    data_id: str
    velocity_computed: bool
    velocity_graph_key: Optional[str] = None  # Key for velocity graph in adata.uns
    mode: str  # RNA velocity computation mode

    model_config = ConfigDict(arbitrary_types_allowed=True)


class TrajectoryResult(BaseModel):
    """Result of trajectory analysis"""

    data_id: str
    pseudotime_computed: bool
    velocity_computed: bool
    pseudotime_key: str
    method: str  # Trajectory analysis method used
    spatial_weight: float  # Spatial information weight

    model_config = ConfigDict(arbitrary_types_allowed=True)


class IntegrationResult(BaseModel):
    """Result of sample integration"""

    data_id: str
    n_samples: int
    integration_method: str

    model_config = ConfigDict(arbitrary_types_allowed=True)


class DeconvolutionResult(BaseModel):
    """Result of spatial deconvolution

    Attributes:
        data_id: Dataset identifier
        method: Deconvolution method used
        dominant_type_key: Column name in adata.obs where dominant cell type is stored (e.g., "dominant_celltype_cell2location")
        cell_types: List of cell types identified
        n_cell_types: Number of cell types
        proportions_key: Key in adata.obsm where cell type proportions are stored
        statistics: Statistics about the deconvolution results
    """

    data_id: str
    method: str
    dominant_type_key: str  # Column name where dominant cell type is stored
    cell_types: List[str]
    n_cell_types: int
    proportions_key: str  # Key in adata.obsm where cell type proportions are stored
    statistics: Dict[str, Any]  # Statistics about the deconvolution results

    model_config = ConfigDict(arbitrary_types_allowed=True)


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

    model_config = ConfigDict(arbitrary_types_allowed=True)


class SpatialVariableGenesResult(BaseModel):
    """Result of spatial variable genes identification"""

    data_id: str
    method: str  # Method used for analysis

    # Common results for all methods
    n_genes_analyzed: int  # Total number of genes analyzed (input)
    n_significant_genes: int  # Total number of significant genes found (q < 0.05)
    n_returned_genes: int  # Number of genes actually returned
    spatial_genes: List[str]  # List of returned gene names (length = n_returned_genes)

    # Statistical results (available for all methods)
    gene_statistics: Dict[str, float]  # Gene name -> primary statistic value
    p_values: Dict[str, float]  # Gene name -> p-value
    q_values: Dict[str, float]  # Gene name -> FDR-corrected p-value

    # Storage keys for results in adata
    results_key: str  # Base key for storing results in adata

    # Method-specific results (optional, only populated for respective methods)
    spatialde_results: Optional[Dict[str, Any]] = None  # SpatialDE-specific results
    sparkx_results: Optional[Dict[str, Any]] = None  # SPARK-X specific results

    model_config = ConfigDict(arbitrary_types_allowed=True)


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

    model_config = ConfigDict(arbitrary_types_allowed=True)


class EnrichmentResult(BaseModel):
    """Result from gene set enrichment analysis

    Note on serialization:
        To minimize MCP response size (~12k tokens -> ~0.5k tokens), large
        dictionaries are excluded from JSON serialization using Field(exclude=True).
        These fields are still stored in the Python object and saved to adata.uns
        for downstream visualization.

        Fields included in MCP response (sent to LLM):
        - method, n_gene_sets, n_significant (basic info)
        - top_gene_sets, top_depleted_sets (top 10 pathway names)
        - spatial_scores_key (for spatial methods)

        Fields excluded from MCP response (stored in adata.uns):
        - enrichment_scores, pvalues, adjusted_pvalues (full dicts)
        - gene_set_statistics (detailed stats per pathway)
        - spatial_metrics (spatial autocorrelation data)
    """

    # Basic information - always included in MCP response
    method: str  # Method used (pathway_gsea, pathway_ora, etc.)
    n_gene_sets: int  # Number of gene sets analyzed
    n_significant: int  # Number of significant gene sets

    # Top results - always included (compact, just pathway names)
    top_gene_sets: List[str]  # Top enriched gene sets (max 10)
    top_depleted_sets: List[str]  # Top depleted gene sets (max 10)

    # Spatial info key - included
    spatial_scores_key: Optional[str] = None  # Key in adata.obsm

    # ============================================================
    # EXCLUDED FROM MCP RESPONSE - stored in adata.uns for viz
    # Full data available via visualize_data() tool
    # ============================================================
    enrichment_scores: Dict[str, float] = Field(
        default_factory=dict,
        exclude=True,  # Exclude from JSON serialization to LLM
    )
    pvalues: Optional[Dict[str, float]] = Field(
        default=None,
        exclude=True,
    )
    adjusted_pvalues: Optional[Dict[str, float]] = Field(
        default=None,
        exclude=True,
    )
    gene_set_statistics: Dict[str, Dict[str, Any]] = Field(
        default_factory=dict,
        exclude=True,
    )
    spatial_metrics: Optional[Dict[str, Any]] = Field(
        default=None,
        exclude=True,
    )

    model_config = ConfigDict(arbitrary_types_allowed=True)


class CNVResult(BaseModel):
    """Result of Copy Number Variation (CNV) analysis

    Attributes:
        data_id: Dataset identifier
        method: CNV inference method used (e.g., "infercnvpy")
        reference_key: Column name in adata.obs used for reference cell types
        reference_categories: List of cell types/clusters used as normal reference
        n_chromosomes: Number of chromosomes analyzed
        n_genes_analyzed: Number of genes included in CNV analysis
        cnv_score_key: Key in adata.obsm where CNV scores are stored (e.g., "X_cnv")
        statistics: Statistics about CNV detection (e.g., mean CNV, cell counts)
        visualization_available: Whether CNV heatmap visualization data is available
    """

    data_id: str
    method: str  # Method used (e.g., "infercnvpy")
    reference_key: str  # Column used for reference cells
    reference_categories: List[str]  # Categories used as reference
    n_chromosomes: int  # Number of chromosomes analyzed
    n_genes_analyzed: int  # Number of genes analyzed
    cnv_score_key: Optional[str] = None  # Key in adata.obsm (e.g., "X_cnv")
    statistics: Optional[Dict[str, Any]] = None  # CNV statistics
    visualization_available: bool = False  # Whether visualization is available

    model_config = ConfigDict(arbitrary_types_allowed=True)
