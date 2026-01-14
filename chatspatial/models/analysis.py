"""
Analysis result models for spatial transcriptomics data.
"""

from typing import TYPE_CHECKING, Any, Optional

from pydantic import BaseModel, ConfigDict, Field

if TYPE_CHECKING:
    from mcp.types import ImageContent
else:
    try:
        from mcp.types import ImageContent
    except ImportError:
        # Fallback for when MCP is not available
        ImageContent = Any  # type: ignore[misc,assignment]


class BaseAnalysisResult(BaseModel):
    """Base class for all analysis results.

    Provides common configuration and optional shared fields.
    All analysis result models should inherit from this class.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)


class PreprocessingResult(BaseAnalysisResult):
    """Result of data preprocessing"""

    data_id: str
    n_cells: int
    n_genes: int
    n_hvgs: int
    clusters: int
    qc_metrics: Optional[dict[str, Any]] = None


class DifferentialExpressionResult(BaseAnalysisResult):
    """Result of differential expression analysis

    Note on serialization:
        For consistency with other result models, the statistics dict is excluded
        from JSON serialization. Key summary info is in explicit fields.

        Fields included in MCP response:
        - data_id, comparison (basic info)
        - n_genes (count)
        - top_genes (top differentially expressed genes)

        Fields excluded from MCP response:
        - statistics (detailed DE metrics per group)
    """

    data_id: str
    comparison: str
    n_genes: int
    top_genes: list[str] = Field(default_factory=list)

    # Detailed statistics - excluded from MCP response
    statistics: dict[str, Any] = Field(
        default_factory=dict,
        exclude=True,  # Exclude from JSON serialization to LLM
    )


class AnnotationResult(BaseAnalysisResult):
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
    cell_types: list[str]
    counts: dict[str, int]
    confidence_scores: Optional[dict[str, float]] = None
    tangram_mapping_score: Optional[float] = None  # For Tangram method - mapping score


class SpatialStatisticsResult(BaseAnalysisResult):
    """Result of spatial analysis

    Note on serialization:
        To minimize MCP response size, detailed per-gene/per-spot statistics are
        excluded from JSON serialization using Field(exclude=True). Summary fields
        are always included.

        Fields included in MCP response:
        - data_id, analysis_type (basic info)
        - n_features_analyzed, n_significant (summary counts)
        - top_features (top significant genes/clusters)
        - summary_metrics (compact key metrics)
        - results_key (for accessing full results)

        Fields excluded from MCP response (stored in adata):
        - statistics (full detailed results dict)

        Visualization is handled separately via the visualize_data tool.
    """

    data_id: str
    analysis_type: str

    # Summary fields - always included in MCP response
    n_features_analyzed: int = 0
    n_significant: int = 0
    top_features: list[str] = Field(default_factory=list)
    summary_metrics: dict[str, float] = Field(default_factory=dict)
    results_key: Optional[str] = None  # Key in adata.uns for full results

    # Detailed statistics - excluded from MCP response
    statistics: Optional[dict[str, Any]] = Field(
        default=None,
        exclude=True,  # Exclude from JSON serialization to LLM
    )


class RNAVelocityResult(BaseAnalysisResult):
    """Result of RNA velocity analysis"""

    data_id: str
    velocity_computed: bool
    velocity_graph_key: Optional[str] = None  # Key for velocity graph in adata.uns
    mode: str  # RNA velocity computation mode


class TrajectoryResult(BaseAnalysisResult):
    """Result of trajectory analysis"""

    data_id: str
    pseudotime_computed: bool
    velocity_computed: bool
    pseudotime_key: str
    method: str  # Trajectory analysis method used
    spatial_weight: float  # Spatial information weight


class IntegrationResult(BaseAnalysisResult):
    """Result of sample integration"""

    data_id: str
    n_samples: int
    integration_method: str


class DeconvolutionResult(BaseAnalysisResult):
    """Result of spatial deconvolution

    Note on serialization:
        To minimize MCP response size, detailed per-cell-type statistics are
        excluded from JSON serialization using Field(exclude=True).

        Fields included in MCP response:
        - data_id, method, n_cell_types, cell_types (basic info)
        - n_spots, genes_used (summary counts)
        - dominant_type_key, proportions_key (storage keys)

        Fields excluded from MCP response (stored in adata):
        - statistics (includes mean_proportions, dominant_types dicts)
    """

    data_id: str
    method: str
    dominant_type_key: str  # Column name where dominant cell type is stored
    cell_types: list[str]
    n_cell_types: int
    proportions_key: str  # Key in adata.obsm where cell type proportions are stored

    # Summary fields - always included
    n_spots: int = 0
    genes_used: int = 0

    # Detailed statistics - excluded from MCP response
    statistics: dict[str, Any] = Field(
        default_factory=dict,
        exclude=True,  # Exclude from JSON serialization to LLM
    )


class SpatialDomainResult(BaseAnalysisResult):
    """Result of spatial domain identification

    Note on serialization:
        For consistency with other result models, the detailed statistics dict
        is excluded from JSON serialization. Key summary info is in explicit fields.

        Fields included in MCP response:
        - data_id, method, n_domains (basic info)
        - domain_key, refined_domain_key, embeddings_key (storage keys)
        - domain_counts (number of spots per domain - typically compact)

        Fields excluded from MCP response:
        - statistics (method parameters, stored in adata.uns)
    """

    data_id: str
    method: str
    n_domains: int
    domain_key: str  # Key in adata.obs where domain labels are stored
    domain_counts: dict[str, int]  # Number of spots in each domain
    refined_domain_key: Optional[str] = (
        None  # Key for refined domains if refinement was applied
    )
    embeddings_key: Optional[str] = (
        None  # Key in adata.obsm where embeddings are stored
    )

    # Detailed statistics - excluded from MCP response
    statistics: dict[str, Any] = Field(
        default_factory=dict,
        exclude=True,  # Exclude from JSON serialization to LLM
    )


class SpatialVariableGenesResult(BaseAnalysisResult):
    """Result of spatial variable genes identification.

    Note on serialization:
        To minimize MCP response size, detailed statistics are excluded from
        JSON serialization using Field(exclude=True). These fields are still
        stored in the Python object and saved to adata.var for downstream
        visualization and export.

        Access complete statistics via:
        - adata.var['spatialde_pval'], adata.var['spatialde_qval'] (SpatialDE)
        - adata.var['sparkx_pval'], adata.var['sparkx_qval'] (SPARK-X)
    """

    data_id: str
    method: str  # Method used for analysis

    # Summary statistics - always returned to LLM
    n_genes_analyzed: int  # Total number of genes analyzed
    n_significant_genes: int  # Total significant genes found (q < 0.05)

    # Top spatial genes - returned to LLM (truncated for token efficiency)
    spatial_genes: list[str]

    # Storage key for accessing full results in adata
    results_key: str

    # ============================================================
    # Fields excluded from MCP response (stored in adata.var)
    # ============================================================
    gene_statistics: dict[str, float] = Field(
        default_factory=dict,
        exclude=True,  # Exclude from JSON serialization to LLM
    )
    p_values: dict[str, float] = Field(
        default_factory=dict,
        exclude=True,
    )
    q_values: dict[str, float] = Field(
        default_factory=dict,
        exclude=True,
    )
    spatialde_results: Optional[dict[str, Any]] = Field(
        default=None,
        exclude=True,
    )
    sparkx_results: Optional[dict[str, Any]] = Field(
        default=None,
        exclude=True,
    )


class CellCommunicationResult(BaseAnalysisResult):
    """Result of cell-cell communication analysis

    Note on serialization:
        To minimize MCP response size, detailed statistics are excluded from
        JSON serialization. Key summary info is in explicit fields.

        Fields included in MCP response:
        - data_id, method, species, database (basic info)
        - n_lr_pairs, n_significant_pairs, top_lr_pairs (summary)
        - Various *_key fields (storage keys for accessing full results)

        Fields excluded from MCP response:
        - statistics (detailed analysis metrics)
    """

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
    top_lr_pairs: list[str] = Field(default_factory=list)  # Top significant LR pairs

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

    # Detailed statistics - excluded from MCP response
    statistics: dict[str, Any] = Field(
        default_factory=dict,
        exclude=True,  # Exclude from JSON serialization to LLM
    )


class EnrichmentResult(BaseAnalysisResult):
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
    top_gene_sets: list[str]  # Top enriched gene sets (max 10)
    top_depleted_sets: list[str]  # Top depleted gene sets (max 10)

    # Spatial info key - included
    spatial_scores_key: Optional[str] = None  # Key in adata.obsm

    # ============================================================
    # EXCLUDED FROM MCP RESPONSE - stored in adata.uns for viz
    # Full data available via visualize_data() tool
    # ============================================================
    enrichment_scores: dict[str, float] = Field(
        default_factory=dict,
        exclude=True,  # Exclude from JSON serialization to LLM
    )
    pvalues: Optional[dict[str, float]] = Field(
        default=None,
        exclude=True,
    )
    adjusted_pvalues: Optional[dict[str, float]] = Field(
        default=None,
        exclude=True,
    )
    gene_set_statistics: dict[str, dict[str, Any]] = Field(
        default_factory=dict,
        exclude=True,
    )
    spatial_metrics: Optional[dict[str, Any]] = Field(
        default=None,
        exclude=True,
    )


class CNVResult(BaseAnalysisResult):
    """Result of Copy Number Variation (CNV) analysis

    Note on serialization:
        For consistency with other result models, the statistics dict is excluded
        from JSON serialization. Key summary info is in explicit fields.

        Fields included in MCP response:
        - data_id, method, reference_key, reference_categories (basic info)
        - n_chromosomes, n_genes_analyzed (summary counts)
        - cnv_score_key (storage key)
        - visualization_available (status flag)

        Fields excluded from MCP response:
        - statistics (detailed CNV metrics)
    """

    data_id: str
    method: str  # Method used (e.g., "infercnvpy")
    reference_key: str  # Column used for reference cells
    reference_categories: list[str]  # Categories used as reference
    n_chromosomes: int  # Number of chromosomes analyzed
    n_genes_analyzed: int  # Number of genes analyzed
    cnv_score_key: Optional[str] = None  # Key in adata.obsm (e.g., "X_cnv")
    visualization_available: bool = False  # Whether visualization is available

    # Detailed statistics - excluded from MCP response
    statistics: Optional[dict[str, Any]] = Field(
        default=None,
        exclude=True,  # Exclude from JSON serialization to LLM
    )


class DEGene(BaseAnalysisResult):
    """A single differentially expressed gene with statistics"""

    gene: str
    log2fc: float
    pvalue: float
    padj: float
    mean_expr_condition1: Optional[float] = None
    mean_expr_condition2: Optional[float] = None


class CellTypeComparisonResult(BaseAnalysisResult):
    """Differential expression result for a single cell type"""

    cell_type: str
    n_cells_condition1: int
    n_cells_condition2: int
    n_samples_condition1: int
    n_samples_condition2: int
    n_significant_genes: int
    top_upregulated: list[DEGene]  # Upregulated in condition1
    top_downregulated: list[DEGene]  # Downregulated in condition1
    all_de_genes: list[DEGene] = Field(
        default_factory=list,
        exclude=True,  # Exclude from MCP response to reduce size
    )


class ConditionComparisonResult(BaseAnalysisResult):
    """Result of multi-sample condition comparison analysis.

    Attributes:
        data_id: Dataset identifier
        method: Method used for differential expression
        comparison: Human-readable comparison string (e.g., "Treatment vs Control")
        condition_key: Column used for condition grouping
        condition1: First condition (experimental group)
        condition2: Second condition (reference group)
        sample_key: Column used for sample identification
        cell_type_key: Column used for cell type stratification (if provided)
        n_samples_condition1: Number of samples in condition1
        n_samples_condition2: Number of samples in condition2
        global_results: Results when no cell type stratification (cell_type_key=None)
        cell_type_results: Results stratified by cell type (when cell_type_key provided)
        results_key: Key in adata.uns where full results are stored
        statistics: Overall statistics about the comparison
    """

    data_id: str
    method: str
    comparison: str
    condition_key: str
    condition1: str
    condition2: str
    sample_key: str
    cell_type_key: Optional[str] = None

    # Sample counts
    n_samples_condition1: int
    n_samples_condition2: int

    # Global results (when cell_type_key is None)
    global_n_significant: Optional[int] = None
    global_top_upregulated: Optional[list[DEGene]] = None
    global_top_downregulated: Optional[list[DEGene]] = None

    # Cell type stratified results (when cell_type_key is provided)
    cell_type_results: Optional[list[CellTypeComparisonResult]] = None

    # Storage keys
    results_key: str  # Key in adata.uns for full results

    # Summary statistics
    statistics: dict[str, Any]
