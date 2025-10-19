"""
Data models for spatial transcriptomics analysis.
"""

from __future__ import annotations

from typing import Annotated, Any, Dict, List, Literal, Optional, Tuple, Union

from pydantic import BaseModel, Field, model_validator
from typing_extensions import Self


class ColumnInfo(BaseModel):
    """Metadata column information for dataset profiling"""

    name: str
    dtype: Literal["categorical", "numerical"]
    n_unique: int
    sample_values: Optional[List[str]] = None  # Sample values for categorical
    range: Optional[Tuple[float, float]] = None  # Value range for numerical


class SpatialDataset(BaseModel):
    """Spatial transcriptomics dataset model with comprehensive metadata profile"""

    id: str
    name: str
    data_type: Literal[
        "10x_visium", "slide_seq", "merfish", "seqfish", "other", "h5ad", "auto"
    ]
    description: Optional[str] = None

    # Basic statistics
    n_cells: int = 0
    n_genes: int = 0
    spatial_coordinates_available: bool = False
    tissue_image_available: bool = False

    # Metadata profiles - let LLM interpret the structure
    obs_columns: Optional[List[ColumnInfo]] = None  # Cell-level metadata
    var_columns: Optional[List[ColumnInfo]] = None  # Gene-level metadata
    obsm_keys: Optional[List[str]] = None  # Multi-dimensional data keys
    uns_keys: Optional[List[str]] = None  # Unstructured data keys

    # Gene expression profiles
    top_highly_variable_genes: Optional[List[str]] = None
    top_expressed_genes: Optional[List[str]] = None


class AnalysisParameters(BaseModel):
    """Analysis parameters model"""

    # Data filtering and subsampling parameters (user controlled)
    filter_genes_min_cells: Optional[Annotated[int, Field(gt=0)]] = (
        3  # Filter genes expressed in < N cells
    )
    filter_cells_min_genes: Optional[Annotated[int, Field(gt=0)]] = (
        30  # Filter cells expressing < N genes
    )
    subsample_spots: Optional[Annotated[int, Field(gt=0, le=50000)]] = (
        None  # Subsample to N spots (None = no subsampling)
    )
    subsample_genes: Optional[Annotated[int, Field(gt=0, le=50000)]] = (
        None  # Keep top N variable genes (None = keep all filtered genes)
    )
    subsample_random_seed: int = 42  # Random seed for subsampling

    # Normalization and scaling parameters
    normalization: Literal["log", "sct", "pearson_residuals", "none", "scvi"] = Field(
        default="log",
        description=(
            "Normalization method for gene expression data.\n\n"
            "AVAILABLE OPTIONS:\n"
            "• 'log' (default): Standard log(x+1) normalization after library size correction. "
            "Robust and widely used for most analyses.\n"
            "• 'pearson_residuals': GLM-based variance stabilization for sparse UMI data. "
            "Requires raw integer counts and scanpy>=1.9.0. Best for single-cell resolution data.\n"
            "• 'none': Skip normalization. Use when data is already pre-normalized.\n\n"
            "NOT IMPLEMENTED (will raise error):\n"
            "• 'sct': SCTransform normalization - not implemented. Will raise NotImplementedError.\n"
            "• 'scvi': scVI normalization - use use_scvi_preprocessing=True instead. Will raise NotImplementedError.\n\n"
            "REQUIREMENTS:\n"
            "• pearson_residuals: Raw count data (integers only), scanpy>=1.9.0, sufficient memory for dense operations\n"
            "• none: Data should already be normalized (will warn if raw counts detected)\n\n"
            "ERROR HANDLING:\n"
            "If the requested normalization fails, an error will be raised with specific reasons "
            "and suggested alternatives. No silent fallbacks will occur.\n\n"
            "RECOMMENDATIONS:\n"
            "• For raw Visium/Xenium/MERFISH data: 'pearson_residuals' or 'log'\n"
            "• For pre-processed data: 'none'\n"
            "• For batch effect correction: use 'log' with batch correction or use_scvi_preprocessing=True"
        ),
    )
    scale: bool = True
    n_hvgs: Annotated[int, Field(gt=0, le=5000)] = 2000
    n_pcs: Annotated[int, Field(gt=0, le=100)] = 30

    # ========== Normalization Control Parameters ==========
    normalize_target_sum: Optional[float] = Field(
        default=None,  # Adaptive default - uses median counts
        ge=1.0,  # Must be positive if specified
        le=1e8,  # Reasonable upper bound
        description=(
            "Target sum for total count normalization per cell/spot. "
            "Controls the library size after normalization. "
            "\n"
            "RECOMMENDED VALUES BY TECHNOLOGY:\n"
            "• None (default): Uses median of total counts - most adaptive, recommended for unknown data\n"
            "• 1e4 (10,000): Standard for 10x Visium spatial transcriptomics\n"
            "• 1e6 (1,000,000): CPM normalization, standard for MERFISH/CosMx/Xenium\n"
            "• Custom value: Match to your expected counts per cell/spot\n"
            "\n"
            "DECISION GUIDE:\n"
            "- Multi-cellular spots (Visium): Use 1e4\n"
            "- Single-cell imaging (MERFISH, Xenium, CosMx): Use 1e6\n"
            "- High-depth sequencing: Consider 1e5 or higher\n"
            "- Low-depth/targeted panels: Consider 1e3-1e4\n"
            "- Cross-sample integration: Use same value for all samples\n"
            "- Spatial domain analysis: Consider skipping normalization (None)\n"
            "\n"
            "SCIENTIFIC RATIONALE:\n"
            "This parameter scales all cells/spots to have the same total count, "
            "removing technical variation due to sequencing depth or capture efficiency. "
            "The choice affects the magnitude of normalized expression values and "
            "can influence downstream analyses like HVG selection and clustering."
        ),
    )

    scale_max_value: Optional[float] = Field(
        default=10.0,
        ge=1.0,  # Must be positive if specified
        le=100.0,  # Reasonable upper bound
        description=(
            "Maximum value for clipping after scaling to unit variance (in standard deviations). "
            "Prevents extreme outliers from dominating downstream analyses. "
            "\n"
            "RECOMMENDED VALUES:\n"
            "• 10.0 (default): Standard in single-cell field, balances outlier control with data preservation\n"
            "• None: No clipping - preserves all variation, use for high-quality data\n"
            "• 5.0-8.0: More aggressive clipping for noisy data\n"
            "• 15.0-20.0: Less aggressive for clean imaging data\n"
            "\n"
            "DECISION GUIDE BY DATA TYPE:\n"
            "- Standard scRNA-seq or Visium: 10.0\n"
            "- High-quality imaging (MERFISH/Xenium): 15.0 or None\n"
            "- Noisy/low-quality data: 5.0-8.0\n"
            "- Exploratory analysis: Start with 10.0\n"
            "- Final analysis: Consider None to preserve all variation\n"
            "\n"
            "TECHNICAL DETAILS:\n"
            "After scaling each gene to zero mean and unit variance, "
            "values exceeding ±max_value standard deviations are clipped. "
            "This prevents a few extreme values from dominating PCA and clustering. "
            "Lower values increase robustness but may remove biological signal."
        ),
    )

    # scVI preprocessing parameters
    use_scvi_preprocessing: bool = False  # Whether to use scVI for preprocessing
    scvi_n_hidden: int = 128
    scvi_n_latent: int = 10
    scvi_n_layers: int = 1
    scvi_dropout_rate: float = 0.1
    scvi_gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb"

    # ResolVI preprocessing parameters
    use_resolvi_preprocessing: bool = Field(
        default=False,
        description="Use ResolVI for spatial denoising (high-resolution spatial data only)",
    )
    resolvi_params: Optional["ResolVIPreprocessingParameters"] = Field(
        default=None, description="ResolVI-specific parameters (uses defaults if None)"
    )

    # Key naming parameters (configurable hard-coded keys)
    cluster_key: str = Field(
        "leiden", alias="clustering_key"
    )  # Key name for storing clustering results
    spatial_key: Optional[str] = Field(
        default=None,
        description="Spatial coordinate key in obsm (auto-detected if None)",
    )  # Changed from hardcoded "spatial" to allow auto-detection
    batch_key: str = "batch"  # Key name for batch information in obs

    # User-controllable parameters (scientifically-informed defaults)
    n_neighbors: Annotated[int, Field(gt=2, le=100)] = Field(
        default=15,
        description=(
            "Number of neighbors for k-NN graph construction. "
            "Default 15 aligns with Scanpy industry standard and UMAP developer recommendations (10-15 range). "
            "Larger values (20-50) preserve more global structure, smaller values (5-10) emphasize local patterns. "
            "For spatial transcriptomics: 15 captures meaningful tissue neighborhoods in both Visium (55μm) and Visium HD (2μm) data."
        ),
    )
    clustering_resolution: Annotated[float, Field(gt=0.1, le=2.0)] = Field(
        default=1.0,
        description=(
            "Leiden clustering resolution parameter controlling clustering coarseness. "
            "Higher values (1.5-2.0) produce more numerous, smaller clusters; "
            "lower values (0.2-0.5) produce fewer, broader clusters. "
            "Common values: 0.25, 0.5, 1.0. Default 1.0 matches scanpy standard and works well for most spatial datasets."
        ),
    )

    # Advanced preprocessing options
    enable_rna_velocity: bool = False  # Whether to include RNA velocity preprocessing
    velocity_params: Optional["RNAVelocityParameters"] = (
        None  # RNA velocity parameters. If None and enable_rna_velocity=True, defaults will be used
    )

    enable_trajectory_analysis: bool = (
        False  # Whether to include trajectory analysis preprocessing
    )
    dpt_root_cell: Optional[str] = (
        None  # Root cell for diffusion pseudotime (cell barcode)
    )
    enable_spatial_domains: bool = (
        False  # Whether to include spatial domain-specific preprocessing
    )


class VisualizationParameters(BaseModel):
    """Visualization parameters model"""

    model_config = {"extra": "forbid"}  # Strict validation after preprocessing

    feature: Optional[Union[str, List[str]]] = Field(
        None,
        description="Single feature or list of features (accepts both 'feature' and 'features')",
    )  # Single feature or list of features

    @model_validator(mode="before")
    @classmethod
    def handle_features_alias(cls, data):
        """Handle 'features' parameter as alias for 'feature'"""
        if isinstance(data, dict) and "features" in data and "feature" not in data:
            data = data.copy()
            data["feature"] = data.pop("features")
        return data

    plot_type: Literal[
        "spatial",
        "heatmap",
        "violin",
        "umap",
        "spatial_domains",
        "cell_communication",
        "deconvolution",
        "trajectory",
        "rna_velocity",
        "spatial_analysis",
        "multi_gene",
        "lr_pairs",
        "gene_correlation",
        "pathway_enrichment",
        "spatial_enrichment",  # Clear enrichment types
        "spatial_interaction",
        "batch_integration",  # Batch integration quality assessment
        "cnv_heatmap",  # CNV analysis heatmap
        "spatial_cnv",  # CNV spatial projection
        "card_imputation",  # CARD imputation high-resolution results
    ] = "spatial"
    colormap: str = "viridis"

    # Spatial analysis visualization parameters
    analysis_type: Optional[
        Literal[
            "neighborhood",
            "co_occurrence",
            "ripley",
            "moran",
            "centrality",
            "getis_ord",
        ]
    ] = Field(
        None,
        description="Analysis type for spatial_analysis plot type. Determines which spatial analysis result to visualize.",
    )
    cluster_key: Optional[str] = Field(
        None,
        description=(
            "Column name in adata.obs for grouping (e.g., 'leiden', 'louvain', 'cell_type'). "
            "REQUIRED for plot_type='heatmap'. "
            "Optional for other plot types."
        ),
    )

    # Multi-gene visualization parameters
    multi_panel: bool = False  # Whether to create multi-panel plots
    panel_layout: Optional[Tuple[int, int]] = (
        None  # (rows, cols) - auto-determined if None
    )

    # Ligand-receptor pair parameters
    lr_pairs: Optional[List[Tuple[str, str]]] = None  # List of (ligand, receptor) pairs
    lr_database: str = "cellchat"  # Database for LR pairs

    # Gene correlation parameters
    correlation_method: Literal["pearson", "spearman", "kendall"] = "pearson"
    show_correlation_stats: bool = True

    # Figure parameters
    figure_size: Optional[Tuple[int, int]] = (
        None  # (width, height) - auto-determined if None
    )
    dpi: int = 300  # Publication quality (Nature/Cell standard)
    alpha: float = 0.9  # Spot transparency (higher = more opaque)
    spot_size: Optional[float] = None  # Auto-determined if None (recommended)
    alpha_img: float = Field(
        0.5,
        ge=0.0,
        le=1.0,
        description="Background tissue image transparency (lower = dimmer, helps spots stand out)",
    )

    # Color parameters
    vmin: Optional[float] = None  # Minimum value for color scale
    vmax: Optional[float] = None  # Maximum value for color scale
    color_scale: Literal["linear", "log", "sqrt"] = "linear"  # Color scaling

    # Display parameters
    title: Optional[str] = None
    show_legend: bool = True
    show_colorbar: bool = True
    show_axes: bool = True
    add_gene_labels: bool = True  # Whether to add gene names as labels

    # Trajectory visualization parameters
    basis: Optional[str] = (
        None  # Basis for trajectory visualization (e.g., 'spatial', 'umap', 'pca')
    )

    # GSEA visualization parameters
    gsea_results_key: str = "gsea_results"  # Key in adata.uns for GSEA results
    gsea_plot_type: Literal["enrichment_plot", "barplot", "dotplot", "spatial"] = (
        "barplot"
    )
    n_top_pathways: int = 10  # Number of top pathways to show in barplot

    # NEW: Spatial plot enhancement parameters
    add_outline: bool = Field(
        False, description="Add cluster outline/contour overlay to spatial plots"
    )
    outline_color: str = Field("black", description="Color for cluster outlines")
    outline_width: float = Field(
        0.4, description="Line width for cluster outlines (Nature/Cell standard)"
    )
    outline_cluster_key: Optional[str] = Field(
        None, description="Cluster key for outlines (e.g., 'leiden')"
    )

    # NEW: UMAP enhancement parameters
    size_by: Optional[str] = Field(
        None,
        description="Feature for point size encoding in UMAP (dual color+size encoding)",
    )
    show_velocity: bool = Field(
        False, description="Overlay RNA velocity vectors on UMAP"
    )
    show_trajectory: bool = Field(
        False, description="Show trajectory connections on UMAP (PAGA)"
    )
    velocity_scale: float = Field(1.0, description="Scaling factor for velocity arrows")

    # NEW: Heatmap enhancement parameters
    obs_annotation: Optional[List[str]] = Field(
        None, description="List of obs keys to show as column annotations"
    )
    var_annotation: Optional[List[str]] = Field(
        None, description="List of var keys to show as row annotations"
    )
    annotation_colors: Optional[Dict[str, str]] = Field(
        None, description="Custom colors for annotations"
    )

    # NEW: Integration assessment parameters
    batch_key: str = Field(
        "batch", description="Key in adata.obs for batch/sample identifier"
    )
    integration_method: Optional[str] = Field(
        None, description="Integration method used (for display)"
    )

    # NEW: Network visualization parameters (for neighborhood analysis)
    show_network: bool = Field(
        False,
        description="Show network-style visualization for neighborhood enrichment",
    )
    network_threshold: float = Field(
        2.0, description="Z-score threshold for network edges"
    )
    network_layout: Literal[
        "spring", "circular", "kamada_kawai", "fruchterman_reingold"
    ] = Field("spring", description="Network layout algorithm")

    # Deconvolution visualization parameters
    n_cell_types: Annotated[
        int,
        Field(
            gt=0,
            le=10,
            description="Number of top cell types to show in deconvolution visualization. Must be between 1-10. Default: 4",
        ),
    ] = 4

    @model_validator(mode="after")
    def validate_conditional_parameters(self) -> Self:
        """Validate parameter dependencies and provide helpful error messages."""

        # Spatial analysis validation
        if self.plot_type == "spatial_analysis":
            if not self.analysis_type or (
                isinstance(self.analysis_type, str) and not self.analysis_type.strip()
            ):
                available_types = [
                    "neighborhood",
                    "co_occurrence",
                    "ripley",
                    "moran",
                    "centrality",
                    "getis_ord",
                ]
                raise ValueError(
                    f"Parameter dependency error: analysis_type is required when plot_type='spatial_analysis'.\n"
                    f"Available analysis types: {', '.join(available_types)}\n"
                    f"Example usage: VisualizationParameters(plot_type='spatial_analysis', analysis_type='neighborhood')\n"
                    f"For more details, see spatial analysis documentation."
                )

        # Future: Add other conditional validations here
        # if self.plot_type == "cell_communication" and not self.method:
        #     raise ValueError("method required for cell_communication plot_type")

        return self


class AnnotationParameters(BaseModel):
    """Cell type annotation parameters model"""

    method: Literal[
        "tangram",
        "scanvi",
        "cellassign",
        "mllmcelltype",
        "sctype",
        "singler",
    ] = "tangram"
    marker_genes: Optional[Dict[str, List[str]]] = None
    reference_data: Optional[str] = None
    reference_data_id: Optional[str] = (
        None  # For Tangram method - ID of reference single-cell dataset
    )
    training_genes: Optional[List[str]] = (
        None  # For Tangram method - genes to use for mapping
    )
    num_epochs: int = (
        100  # For Tangram/ScanVI methods - number of training epochs (reduced for faster training)
    )
    mode: Literal["cells", "clusters"] = "cells"  # For Tangram method - mapping mode
    cluster_label: Optional[str] = (
        None  # For mLLMCellType method - cluster label in spatial data. Only required when method='mllmcelltype'
    )
    cell_type_key: Optional[str] = Field(
        default=None,
        description=(
            "Column name for cell types in REFERENCE data. "
            "\n\n"
            "REQUIRED FOR METHODS USING REFERENCE DATA:\n"
            "  • tangram: REQUIRED - maps spatial data to reference using cell type labels\n"
            "  • scanvi: REQUIRED - transfers labels from reference to query data\n"
            "  • singler: REQUIRED - correlates expression with reference cell types\n"
            "\n"
            "NOT REQUIRED FOR METHODS WITHOUT REFERENCE:\n"
            "  • cellassign: Not needed - uses marker_genes parameter instead\n"
            "  • sctype: Not needed - uses built-in database or custom markers\n"
            "  • mllmcelltype: Not needed - uses LLM for annotation\n"
            "\n"
            "Common column names in reference data: 'cell_type', 'cell_types', 'celltype', 'annotation', 'label', 'cell_type_original'\n"
            "\n"
            "The LLM will auto-detect from metadata if not specified, but explicit specification is recommended."
        ),
    )

    # Tangram-specific parameters (aligned with official API)
    tangram_density_prior: Literal["rna_count_based", "uniform"] = (
        "rna_count_based"  # Density prior for mapping
    )
    tangram_device: str = "cpu"  # Device for computation ('cpu' or 'cuda:0')
    tangram_lambda_r: Optional[float] = None  # Regularization parameter
    tangram_lambda_neighborhood: Optional[float] = None  # Spatial regularization
    tangram_learning_rate: float = 0.1  # Learning rate for optimization
    tangram_compute_validation: bool = False  # Whether to compute validation metrics
    tangram_project_genes: bool = False  # Whether to project gene expression

    # General parameters for batch effect and data handling
    batch_key: Optional[str] = None  # For batch effect correction
    layer: Optional[str] = None  # Which layer to use for analysis

    # scANVI parameters (scvi-tools semi-supervised label transfer)
    scanvi_n_hidden: int = Field(
        default=128,
        description="Number of hidden units per layer. Official default: 128",
    )
    scanvi_n_latent: int = Field(
        default=10,
        description=(
            "Dimensionality of latent space. Official default: 10\n"
            "scvi-tools recommendation for large integration: 30\n"
            "⚠️  Empirical (not official): Small datasets may need 3-5 to avoid NaN"
        ),
    )
    scanvi_n_layers: int = Field(
        default=1,
        description=(
            "Number of hidden layers. Official default: 1\n"
            "scvi-tools recommendation for large integration: 2"
        ),
    )
    scanvi_dropout_rate: float = Field(
        default=0.1,
        description=(
            "Dropout rate for regularization. Official default: 0.1\n"
            "⚠️  Empirical (not official): 0.2-0.3 may help small datasets"
        ),
    )
    scanvi_unlabeled_category: str = Field(
        default="Unknown",
        description="Label for unlabeled cells in semi-supervised learning",
    )

    # SCVI pretraining parameters (official best practice)
    scanvi_use_scvi_pretrain: bool = Field(
        default=True,
        description=(
            "Whether to pretrain with SCVI before SCANVI training. Default: True\n"
            "Official scvi-tools best practice: SCVI pretraining improves stability\n"
            "⚠️  For small datasets: Set to False if encountering NaN errors"
        ),
    )
    scanvi_scvi_epochs: int = Field(
        default=200, description="Number of epochs for SCVI pretraining. Default: 200"
    )
    scanvi_n_samples_per_label: int = Field(
        default=100,
        description="Number of samples per label for semi-supervised training",
    )

    # Query training parameters
    scanvi_query_epochs: int = Field(
        default=100,
        description=(
            "Number of epochs for training on query data. Default: 100\n"
            "⚠️  For small datasets: Recommend 50 to prevent overfitting"
        ),
    )
    scanvi_check_val_every_n_epoch: int = Field(
        default=10, description="Validation check frequency during training"
    )

    # CellAssign parameters
    cellassign_n_hidden: int = 100
    cellassign_learning_rate: float = 0.001
    cellassign_max_iter: int = 200

    # mLLMCellType parameters
    mllm_n_marker_genes: Annotated[int, Field(gt=0, le=50)] = (
        20  # Number of marker genes per cluster
    )
    mllm_species: Literal["human", "mouse"] = "human"  # Species
    mllm_tissue: Optional[str] = None  # Tissue type (e.g., "brain", "liver")
    mllm_provider: Literal[
        "openai",
        "anthropic",
        "gemini",
        "deepseek",
        "qwen",
        "zhipu",
        "stepfun",
        "minimax",
        "grok",
        "openrouter",
    ] = "openai"  # LLM provider (use 'gemini' not 'google')
    mllm_model: Optional[str] = (
        None  # Model name. Defaults: openai="gpt-5", anthropic="claude-sonnet-4-20250514", gemini="gemini-2.5-pro-preview-03-25"
        # Examples: "gpt-5", "claude-sonnet-4-5-20250929", "claude-opus-4-1-20250805", "gemini-2.5-pro", "qwen-max-2025-01-25"
    )
    mllm_api_key: Optional[str] = None  # API key for the LLM provider
    mllm_additional_context: Optional[str] = None  # Additional context for annotation
    mllm_use_cache: bool = True  # Whether to use caching for API calls
    mllm_base_urls: Optional[Union[str, Dict[str, str]]] = None  # Custom API endpoints
    mllm_verbose: bool = False  # Whether to print detailed logs
    mllm_force_rerun: bool = False  # Force reanalysis bypassing cache

    # Multi-model consensus parameters (interactive_consensus_annotation)
    mllm_use_consensus: bool = False  # Whether to use multi-model consensus
    mllm_models: Optional[List[Union[str, Dict[str, str]]]] = (
        None  # List of models for consensus
    )
    mllm_api_keys: Optional[Dict[str, str]] = None  # Dict mapping provider to API key
    mllm_consensus_threshold: float = 0.7  # Agreement threshold for consensus
    mllm_entropy_threshold: float = 1.0  # Entropy threshold for controversy detection
    mllm_max_discussion_rounds: int = 3  # Maximum discussion rounds
    mllm_consensus_model: Optional[Union[str, Dict[str, str]]] = (
        None  # Model for consensus checking
    )
    mllm_clusters_to_analyze: Optional[List[str]] = None  # Specific clusters to analyze

    # ScType parameters
    sctype_tissue: Optional[str] = (
        None  # Tissue type (supported: "Adrenal", "Brain", "Eye", "Heart", "Hippocampus", "Immune system", "Intestine", "Kidney", "Liver", "Lung", "Muscle", "Pancreas", "Placenta", "Spleen", "Stomach", "Thymus")
    )
    sctype_db_: Optional[str] = (
        None  # Custom database path (if None, uses default ScTypeDB)
    )
    sctype_scaled: bool = True  # Whether input data is scaled
    sctype_custom_markers: Optional[Dict[str, Dict[str, List[str]]]] = (
        None  # Custom markers: {"CellType": {"positive": [...], "negative": [...]}}
    )
    sctype_use_cache: bool = True  # Whether to cache results to avoid repeated R calls

    # SingleR parameters (for enhanced marker_genes method)
    singler_reference: Optional[str] = (
        None  # Reference name from celldex (e.g., "blueprint_encode", "dice", "hpca")
    )
    singler_integrated: bool = (
        False  # Whether to use integrated annotation with multiple references
    )
    num_threads: int = 4  # Number of threads for parallel processing


class SpatialAnalysisParameters(BaseModel):
    """Spatial analysis parameters model"""

    analysis_type: Literal[
        "neighborhood",
        "co_occurrence",
        "ripley",
        "moran",
        "local_moran",  # Added: Local Moran's I (LISA)
        "geary",
        "centrality",
        "getis_ord",
        "bivariate_moran",
        "join_count",  # Traditional Join Count for binary data (2 categories)
        "local_join_count",  # Local Join Count for multi-category data (>2 categories)
        "network_properties",
        "spatial_centrality",
    ] = "neighborhood"
    cluster_key: Optional[str] = Field(
        default=None,
        description=(
            "Column name for cluster/cell type labels in adata.obs. "
            "\n\n"
            "REQUIRED FOR GROUP-BASED ANALYSES:\n"
            "  • neighborhood: REQUIRED - analyzes enrichment between cell type groups\n"
            "  • co_occurrence: REQUIRED - measures spatial co-occurrence of groups\n"
            "  • ripley: REQUIRED - analyzes spatial point patterns by group\n"
            "  • join_count: REQUIRED - for BINARY categorical data (2 categories)\n"
            "  • local_join_count: REQUIRED - for MULTI-CATEGORY data (>2 categories)\n"
            "\n"
            "OPTIONAL/NOT REQUIRED FOR GENE-BASED ANALYSES:\n"
            "  • moran: Not required - analyzes gene expression spatial patterns\n"
            "  • local_moran: Not required - identifies local spatial clusters for genes\n"
            "  • geary: Not required - measures gene expression spatial autocorrelation\n"
            "  • getis_ord: Not required - detects hot/cold spots for gene expression\n"
            "  • bivariate_moran: Not required - analyzes gene pair spatial correlation\n"
            "  • centrality: Not required - computes spatial network centrality\n"
            "  • network_properties: Not required - analyzes spatial network structure\n"
            "  • spatial_centrality: Not required - measures spatial importance\n"
            "\n"
            "Common column names: 'leiden', 'louvain', 'cell_type', 'cell_type_tangram', 'seurat_clusters', 'clusters'\n"
            "\n"
            "The LLM will auto-detect from metadata if not specified for required analyses."
        ),
    )
    n_neighbors: Annotated[int, Field(gt=0)] = Field(
        8,
        description=(
            "Number of nearest neighbors for spatial graph construction. "
            "Default: 8 (recommended by ArcGIS for Getis-Ord analysis). "
            "Adjust based on dataset density and spatial scale."
        ),
    )

    # Unified gene selection parameter (NEW)
    genes: Optional[List[str]] = Field(
        None,
        description="Specific genes to analyze. If None, uses HVG or defaults based on analysis type",
    )
    n_top_genes: Annotated[int, Field(gt=0, le=500)] = Field(
        20,
        description="Number of top HVGs to analyze (default 20, up to 500 for comprehensive analysis)",
    )

    # Parallel processing parameters
    n_jobs: Optional[int] = Field(
        1,
        description="Number of parallel jobs. 1 = no parallelization (recommended for small datasets), None = auto-detect, -1 = all cores",
    )
    backend: Literal["loky", "threading", "multiprocessing"] = Field(
        "threading",
        description="Parallelization backend (threading is safer than loky)",
    )

    # Moran's I specific parameters
    moran_n_perms: Annotated[int, Field(gt=0, le=10000)] = Field(
        10,
        description="Number of permutations (default 10 for speed, use 100+ for publication)",
    )
    moran_two_tailed: bool = Field(False, description="Use two-tailed test")

    # Getis-Ord Gi* specific parameters
    getis_ord_correction: Literal["bonferroni", "fdr_bh", "none"] = Field(
        "fdr_bh",
        description=(
            "Multiple testing correction method for Getis-Ord analysis. "
            "Options: 'fdr_bh' (Benjamini-Hochberg FDR, recommended for multi-gene), "
            "'bonferroni' (conservative), 'none' (no correction)"
        ),
    )
    getis_ord_alpha: Annotated[float, Field(gt=0.0, le=1.0)] = Field(
        0.05,
        description=(
            "Significance level (alpha) for Getis-Ord hotspot detection. "
            "Determines Z-score threshold via norm.ppf(1 - alpha/2). "
            "Common values: 0.05 (z=1.96), 0.01 (z=2.576), 0.10 (z=1.645)"
        ),
    )

    # Bivariate Moran's I specific parameters
    gene_pairs: Optional[List[Tuple[str, str]]] = Field(
        None, description="Gene pairs for bivariate analysis"
    )


class RNAVelocityParameters(BaseModel):
    """RNA velocity analysis parameters model"""

    model_config = {
        "extra": "forbid"
    }  # Strict validation - no extra parameters allowed

    # Velocity computation method selection
    method: Literal["scvelo", "velovi"] = "scvelo"

    # scVelo specific parameters
    mode: Literal["deterministic", "stochastic", "dynamical"] = "stochastic"
    n_pcs: Annotated[int, Field(gt=0, le=100)] = 30
    basis: str = "spatial"

    # Preprocessing parameters for velocity computation
    min_shared_counts: Annotated[int, Field(gt=0)] = (
        30  # Minimum shared counts for filtering
    )
    n_top_genes: Annotated[int, Field(gt=0)] = 2000  # Number of top genes to retain
    n_neighbors: Annotated[int, Field(gt=0)] = (
        30  # Number of neighbors for moments computation
    )

    # VELOVI specific parameters
    velovi_n_hidden: int = 128
    velovi_n_latent: int = 10
    velovi_n_layers: int = 1
    velovi_n_epochs: int = 1000
    velovi_dropout_rate: float = 0.1
    velovi_learning_rate: float = 1e-3
    velovi_use_gpu: bool = False


class TrajectoryParameters(BaseModel):
    """Trajectory analysis parameters model"""

    method: Literal["cellrank", "palantir", "dpt"] = "cellrank"
    spatial_weight: Annotated[float, Field(ge=0.0, le=1.0)] = (
        0.5  # Spatial information weight
    )
    root_cells: Optional[List[str]] = None  # For Palantir method

    # CellRank specific parameters
    cellrank_kernel_weights: Tuple[float, float] = (
        0.8,
        0.2,
    )  # (velocity_weight, connectivity_weight)
    cellrank_n_states: Annotated[int, Field(gt=0, le=20)] = (
        5  # Number of macrostates for CellRank
    )

    # Palantir specific parameters
    palantir_n_diffusion_components: Annotated[int, Field(gt=0, le=50)] = (
        10  # Number of diffusion components
    )
    palantir_num_waypoints: Annotated[int, Field(gt=0)] = (
        500  # Number of waypoints for Palantir
    )

    # Fallback control
    # Removed: allow_fallback_to_dpt - No longer doing automatic fallbacks
    # LLMs should explicitly choose which method to use


class IntegrationParameters(BaseModel):
    """Sample integration parameters model"""

    method: Literal["harmony", "bbknn", "scanorama", "scvi"] = "harmony"
    batch_key: str = "batch"  # Batch information key
    n_pcs: Annotated[int, Field(gt=0, le=100)] = (
        30  # Number of principal components for integration
    )
    align_spatial: bool = True  # Whether to align spatial coordinates
    reference_batch: Optional[str] = None  # Reference batch for spatial alignment

    # Common scvi-tools parameters
    use_gpu: bool = False  # Whether to use GPU acceleration for scvi-tools methods
    n_epochs: Optional[int] = None  # Number of training epochs (None = auto-determine)

    # scVI integration parameters
    scvi_n_hidden: int = 128
    scvi_n_latent: int = 10
    scvi_n_layers: int = 1
    scvi_dropout_rate: float = 0.1
    scvi_gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb"


class DeconvolutionParameters(BaseModel):
    """Spatial deconvolution parameters model"""

    method: Literal[
        "cell2location", "rctd", "destvi", "stereoscope", "spotlight", "tangram", "card"
    ] = "cell2location"
    reference_data_id: Optional[str] = (
        None  # Reference single-cell data for deconvolution
    )
    cell_type_key: str  # REQUIRED: Key in reference data for cell type information. LLM will infer from metadata. Common values: 'cell_type', 'celltype', 'annotation', 'label'
    n_top_genes: Annotated[int, Field(gt=0, le=5000)] = (
        2000  # Number of top genes to use
    )
    use_gpu: bool = False  # Whether to use GPU for cell2location
    ref_model_epochs: Annotated[int, Field(gt=0)] = (
        250  # Number of epochs for reference model training (NB regression). Official recommendation: 250
    )
    n_epochs: Annotated[int, Field(gt=0)] = (
        30000  # Number of epochs for Cell2location spatial mapping model training. Official recommendation: 30000
    )
    n_cells_per_spot: int = (
        30  # Expected number of cells per spatial location (tissue-dependent). Official recommendation: 30
    )
    reference_profiles: Optional[Dict[str, List[float]]] = (
        None  # Reference expression profiles
    )

    # Cell2location specific parameters
    detection_alpha: Annotated[float, Field(gt=0)] = (
        200.0  # RNA detection sensitivity parameter for cell2location. Higher values = less sensitivity correction. Official recommendation: 200
    )

    # SPOTlight specific parameters
    hvg: Optional[int] = None  # Number of highly variable genes to use (None = use all)

    # DestVI parameters
    destvi_n_hidden: int = 128
    destvi_n_latent: int = 10
    destvi_n_layers: int = 1
    destvi_dropout_rate: float = 0.1
    destvi_learning_rate: float = 1e-3

    # Stereoscope parameters
    stereoscope_n_epochs: int = 150000
    stereoscope_learning_rate: float = 0.01
    stereoscope_batch_size: int = 128

    # RCTD specific parameters
    rctd_mode: Literal["full", "doublet", "multi"] = "full"  # RCTD mode
    max_cores: Annotated[int, Field(gt=0, le=16)] = 4  # Maximum number of cores to use
    rctd_confidence_threshold: Annotated[float, Field(gt=0)] = (
        10.0  # Confidence threshold for cell type assignment
    )
    rctd_doublet_threshold: Annotated[float, Field(gt=0)] = (
        25.0  # Threshold for doublet detection
    )

    # CARD specific parameters
    card_minCountGene: Annotated[int, Field(gt=0)] = (
        100  # Minimum gene counts for CARD QC
    )
    card_minCountSpot: Annotated[int, Field(gt=0)] = (
        5  # Minimum spots per gene for CARD QC
    )
    card_sample_key: Optional[str] = (
        None  # Optional sample/batch column in reference data for CARD
    )
    card_imputation: bool = (
        False  # Whether to perform CARD spatial imputation for higher resolution
    )
    card_NumGrids: Annotated[int, Field(gt=0)] = (
        2000  # Number of grids for CARD imputation (default: 2000, higher = finer resolution)
    )
    card_ineibor: Annotated[int, Field(gt=0)] = (
        10  # Number of neighbors for CARD imputation (default: 10)
    )


class SpatialDomainParameters(BaseModel):
    """Spatial domain identification parameters model"""

    method: Literal["spagcn", "leiden", "louvain", "stagate", "graphst"] = "spagcn"
    n_domains: Annotated[int, Field(gt=0, le=50)] = (
        7  # Number of spatial domains to identify
    )

    # SpaGCN specific parameters
    spagcn_s: Annotated[float, Field(gt=0.0)] = (
        1.0  # Weight given to histology in SpaGCN
    )
    spagcn_b: Annotated[int, Field(gt=0)] = (
        49  # Area of each spot when extracting color intensity
    )
    spagcn_p: Annotated[float, Field(ge=0.0, le=1.0)] = (
        0.5  # Percentage of total expression contributed by neighborhoods
    )
    spagcn_use_histology: bool = True  # Whether to use histology image in SpaGCN
    spagcn_random_seed: int = 100  # Random seed for SpaGCN

    # General clustering parameters
    resolution: float = 0.5  # Resolution for leiden/louvain clustering
    use_highly_variable: bool = True  # Whether to use highly variable genes only
    refine_domains: bool = (
        True  # Whether to refine spatial domains using spatial smoothing
    )
    refinement_threshold: Annotated[float, Field(ge=0.0, le=1.0)] = (
        0.5  # Threshold for refinement: only relabel if >=threshold of neighbors differ (0.5 = 50%, following SpaGCN)
    )

    # Clustering-specific parameters for leiden/louvain methods
    cluster_n_neighbors: Optional[Annotated[int, Field(gt=0)]] = (
        None  # Number of neighbors for clustering (default: 15)
    )
    cluster_spatial_weight: Optional[Annotated[float, Field(ge=0.0, le=1.0)]] = (
        None  # Weight for spatial information (default: 0.3)
    )
    cluster_resolution: Optional[float] = None  # Resolution parameter for clustering

    # STAGATE specific parameters
    stagate_rad_cutoff: Optional[float] = (
        None  # Radius cutoff for spatial neighbors (default: 150)
    )
    stagate_learning_rate: Optional[float] = None  # Learning rate (default: 0.001)
    stagate_weight_decay: Optional[float] = None  # Weight decay (default: 0.0001)
    stagate_epochs: Optional[int] = None  # Number of training epochs (default: 1000)
    stagate_dim_output: Optional[int] = (
        None  # Dimension of output representation (default: 15)
    )
    stagate_random_seed: Optional[int] = None  # Random seed (default: 42)

    # GraphST specific parameters
    graphst_use_gpu: bool = False  # Whether to use GPU acceleration
    graphst_clustering_method: Literal["mclust", "leiden", "louvain"] = (
        "leiden"  # Clustering method for GraphST
    )
    graphst_refinement: bool = True  # Whether to refine domains using spatial info
    graphst_radius: int = 50  # Radius for spatial refinement
    graphst_random_seed: int = 42  # Random seed for GraphST
    graphst_n_clusters: Optional[int] = (
        None  # Number of clusters (if None, uses n_domains)
    )

    # Simple timeout configuration
    timeout: Optional[int] = None  # Timeout in seconds (default: 600)


class SpatialVariableGenesParameters(BaseModel):
    """Spatial variable genes identification parameters model"""

    # Method selection
    method: Literal["spatialde", "sparkx"] = (
        "sparkx"  # Default to SPARK-X (best accuracy)
    )

    # Common parameters for all methods
    n_top_genes: Optional[Annotated[int, Field(gt=0, le=5000)]] = (
        None  # Number of top spatial variable genes to return (None = all significant)
    )
    spatial_key: str = "spatial"  # Key in obsm containing spatial coordinates

    # SpatialDE-specific parameters
    spatialde_normalized: bool = True  # Whether data is already normalized
    spatialde_kernel: str = "SE"  # Kernel function type for SpatialDE

    # SPARK-X specific parameters
    sparkx_percentage: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.1  # Percentage of total expression for filtering
    )
    sparkx_min_total_counts: Annotated[int, Field(gt=0)] = (
        10  # Minimum total counts per gene
    )
    sparkx_num_core: Annotated[int, Field(gt=0, le=16)] = (
        1  # Number of cores for parallel processing
    )
    sparkx_option: Literal["single", "mixture"] = (
        "mixture"  # Kernel testing: "single" (faster) or "mixture" (11 kernels)
    )
    sparkx_verbose: bool = False  # Whether to print detailed R output

    # Gene filtering parameters
    filter_mt_genes: bool = (
        True  # Filter mitochondrial genes (MT-*) - standard practice
    )
    filter_ribo_genes: bool = (
        False  # Filter ribosomal genes (RPS*, RPL*) - optional, may remove housekeeping
    )
    test_only_hvg: bool = (
        True  # Test only highly variable genes - 2024 best practice for reducing housekeeping dominance
        # Requires preprocessing with HVG detection first; set to False to test all genes (not recommended)
    )
    warn_housekeeping: bool = True  # Warn if >30% of top genes are housekeeping genes


class CellCommunicationParameters(BaseModel):
    """Cell-cell communication analysis parameters model with explicit user control"""

    # ========== Basic Method Selection ==========
    method: Literal["liana", "cellphonedb", "cellchat_liana"] = "liana"

    # ========== Species and Resource Control ==========
    species: Literal["human", "mouse", "zebrafish"] = (
        "human"  # Species for ligand-receptor database
    )

    # Enhanced resource selection including mouseconsensus
    liana_resource: Literal[
        "consensus",
        "mouseconsensus",
        "cellchat",
        "cellphonedb",
        "connectome",
        "omnipath",
        "celltalkdb",
        "icellnet",
    ] = "consensus"  # LR database resource

    # Species validation control
    validate_species_match: bool = (
        True  # Whether to validate species matches gene patterns
    )
    species_validation_strictness: Literal["strict", "warning", "disabled"] = "strict"
    # strict: Error on mismatch, warning: Warn but continue, disabled: Skip validation

    # ========== Data Usage Control ==========
    # Explicit data source control
    data_source: Literal["current", "raw"] = "current"
    # current: Use adata.X, raw: Use adata.raw (if exists)

    # ========== Spatial Analysis Control ==========
    perform_spatial_analysis: bool = (
        True  # Whether to perform spatial bivariate analysis
    )

    # Spatial connectivity control
    spatial_connectivity_handling: Literal[
        "require_existing", "compute_with_params", "skip"
    ] = "require_existing"
    # require_existing: Must exist or error, compute_with_params: Compute with provided params, skip: Skip spatial analysis

    # Spatial connectivity parameters (required when compute_with_params)
    spatial_neighbors_kwargs: Optional[Dict[str, Any]] = None
    # Example: {"coord_type": "generic", "n_neighs": 6, "radius": 150}

    # ========== Cell Type Control ==========
    # Cell type column control
    cell_type_column: str  # REQUIRED: Which column to use for cell types. LLM will infer from metadata. Common values: 'cell_type', 'celltype', 'leiden', 'louvain', 'seurat_clusters'
    cell_type_handling: Literal["require", "create_from_column", "skip_validation"] = (
        "require"
    )
    # require: Must exist, create_from_column: Create from clustering, skip_validation: Skip (dangerous)

    cell_type_source_column: Optional[str] = (
        None  # Source column when create_from_column
    )
    # Example: "leiden", "louvain", "seurat_clusters"

    # ========== LIANA Specific Parameters ==========
    liana_local_metric: Literal["cosine", "pearson", "spearman", "jaccard"] = (
        "cosine"  # Local spatial metric
    )
    liana_global_metric: Literal["morans", "lee"] = "morans"  # Global spatial metric
    liana_n_perms: Annotated[int, Field(gt=0)] = 100  # Number of permutations for LIANA
    liana_nz_prop: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.2  # Minimum expression proportion
    )
    liana_bandwidth: Optional[int] = None  # Bandwidth for spatial connectivity
    liana_cutoff: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.1  # Cutoff for spatial connectivity
    )

    # ========== Advanced Control ==========
    # Failure handling
    on_validation_failure: Literal["error", "warning", "ignore"] = "error"

    # ========== Expression Filtering Parameters ==========
    min_cells: Annotated[int, Field(ge=0)] = (
        3  # Minimum cells expressing ligand or receptor (required by LIANA for statistical validity)
    )

    # ========== Result Control ==========
    plot_top_pairs: Annotated[int, Field(gt=0, le=20)] = (
        6  # Number of top LR pairs to include in results
    )

    # ========== CellPhoneDB Specific Parameters ==========
    cellphonedb_threshold: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.1  # Expression threshold
    )
    cellphonedb_iterations: Annotated[int, Field(gt=0, le=10000)] = (
        1000  # Statistical permutations
    )
    cellphonedb_result_precision: Annotated[int, Field(gt=0, le=5)] = (
        3  # Result decimal precision
    )
    cellphonedb_pvalue: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.05  # P-value significance threshold
    )
    cellphonedb_use_microenvironments: bool = (
        True  # Whether to use spatial microenvironments
    )
    cellphonedb_spatial_radius: Optional[Annotated[float, Field(gt=0.0)]] = (
        None  # Spatial radius for microenvironments
    )
    cellphonedb_debug_seed: Optional[int] = None  # Random seed for reproducible results

    # ========== Gene Filtering for CellPhoneDB Compatibility (ULTRATHINK) ==========
    enable_gene_filtering: bool = (
        True  # Whether to enable automatic gene filtering for CellPhoneDB compatibility
    )
    gene_filtering_strategy: Literal[
        "none", "conservative", "moderate", "aggressive", "adaptive"
    ] = "conservative"
    # Gene filtering strategies:
    # - none: No filtering applied
    # - conservative: Only remove core problematic genes (ICAM3, ITGAD, CD11D) [DEFAULT]
    # - moderate: Remove core + ITGB2 for added safety
    # - aggressive: Remove all potentially related genes
    # - adaptive: Automatically adjust based on dataset characteristics


class EnrichmentParameters(BaseModel):
    """Parameters for gene set enrichment analysis"""

    # REQUIRED: Species specification (no default value)
    species: Literal["human", "mouse", "zebrafish"]
    # Must explicitly specify the species for gene set matching:
    # - "human": For human data (genes like CD5L, PTPRC - all uppercase)
    # - "mouse": For mouse data (genes like Cd5l, Ptprc - capitalize format)
    # - "zebrafish": For zebrafish data

    # Method selection
    method: Literal[
        "spatial_enrichmap",
        "pathway_gsea",
        "pathway_ora",
        "pathway_enrichr",
        "pathway_ssgsea",
    ] = "spatial_enrichmap"  # Enrichment method

    # Gene sets
    gene_sets: Optional[Union[List[str], Dict[str, List[str]]]] = (
        None  # Gene sets to analyze
    )
    score_keys: Optional[Union[str, List[str]]] = None  # Names for gene signatures

    # Gene set database - choose species-appropriate option
    gene_set_database: Optional[
        Literal[
            "GO_Biological_Process",  # Default (auto-adapts to species)
            "GO_Molecular_Function",  # GO molecular function terms
            "GO_Cellular_Component",  # GO cellular component terms
            "KEGG_Pathways",  # KEGG pathways (species-specific: human=2021, mouse=2019)
            "Reactome_Pathways",  # Reactome pathway database (2022 version)
            "MSigDB_Hallmark",  # MSigDB hallmark gene sets (2020 version)
            "Cell_Type_Markers",  # Cell type marker genes
        ]
    ] = "GO_Biological_Process"

    # Spatial parameters (for spatial_enrichmap)
    spatial_key: str = "spatial"  # Key for spatial coordinates
    n_neighbors: Annotated[int, Field(gt=0)] = 6  # Number of spatial neighbors
    smoothing: bool = True  # Whether to perform spatial smoothing
    correct_spatial_covariates: bool = True  # Whether to correct for spatial covariates

    # Analysis parameters
    batch_key: Optional[str] = None  # Column for batch-wise normalization
    gene_weights: Optional[Dict[str, Dict[str, float]]] = (
        None  # Pre-computed gene weights
    )
    min_genes: Annotated[int, Field(gt=0)] = 10  # Minimum genes in gene set
    max_genes: Annotated[int, Field(gt=0)] = 500  # Maximum genes in gene set

    # Statistical parameters
    pvalue_cutoff: Annotated[float, Field(gt=0.0, lt=1.0)] = 0.05  # P-value cutoff
    adjust_method: Literal["bonferroni", "fdr", "none"] = (
        "fdr"  # Multiple testing correction
    )
    n_permutations: Annotated[int, Field(gt=0)] = (
        1000  # Number of permutations for GSEA
    )

    # Result filtering parameters
    plot_top_terms: Annotated[int, Field(gt=0, le=30)] = (
        10  # Number of top terms to include in results
    )


class CNVParameters(BaseModel):
    """Copy Number Variation (CNV) analysis parameters model"""

    # Method selection
    method: Literal["infercnvpy", "numbat"] = Field(
        "infercnvpy",
        description=(
            "CNV analysis method. 'infercnvpy': expression-based (default), "
            "'numbat': haplotype-aware (requires allele data)"
        ),
    )

    # Reference cell specification
    reference_key: str = Field(
        ...,
        description=(
            "Column name in adata.obs containing cell type or cluster labels "
            "for identifying reference (normal) cells. Common values: "
            "'cell_type', 'leiden', 'louvain', 'seurat_clusters'"
        ),
    )
    reference_categories: List[str] = Field(
        ...,
        description=(
            "List of cell types/clusters to use as reference (normal) cells. "
            "These should be non-malignant cells like immune cells, fibroblasts, etc. "
            "Example: ['T cells', 'B cells', 'Macrophages']"
        ),
    )

    # infercnvpy parameters
    window_size: Annotated[int, Field(gt=0, le=500)] = Field(
        100, description="Number of genes for CNV averaging window (default: 100)"
    )
    step: Annotated[int, Field(gt=0, le=100)] = Field(
        10, description="Step size for sliding window (default: 10)"
    )

    # Analysis options
    exclude_chromosomes: Optional[List[str]] = Field(
        None,
        description=(
            "Chromosomes to exclude from analysis (e.g., ['chrX', 'chrY', 'chrM'])"
        ),
    )
    dynamic_threshold: Optional[float] = Field(
        1.5,
        gt=0.0,
        description="Threshold for dynamic CNV calling (default: 1.5)",
    )

    # Clustering and visualization options (infercnvpy)
    cluster_cells: bool = Field(
        False, description="Whether to cluster cells by CNV pattern"
    )
    dendrogram: bool = Field(
        False, description="Whether to compute hierarchical clustering dendrogram"
    )

    # Numbat-specific parameters
    numbat_genome: Literal["hg38", "hg19", "mm10", "mm39"] = Field(
        "hg38", description="Reference genome for Numbat (default: hg38)"
    )
    numbat_allele_data_key: str = Field(
        "allele_counts",
        description="Layer name in adata containing allele count data",
    )
    numbat_t: Annotated[float, Field(gt=0.0, le=1.0)] = Field(
        0.15, description="Transition probability threshold (default: 0.15)"
    )
    numbat_max_entropy: Annotated[float, Field(gt=0.0, le=1.0)] = Field(
        0.8,
        description=(
            "Maximum entropy threshold. Use 0.8 for spatial data, "
            "0.5 for scRNA-seq (default: 0.8)"
        ),
    )
    numbat_min_cells: Annotated[int, Field(gt=0)] = Field(
        10, description="Minimum cells per CNV event (default: 10)"
    )
    numbat_ncores: Annotated[int, Field(gt=0, le=16)] = Field(
        1, description="Number of cores for parallel processing (default: 1)"
    )
    numbat_skip_nj: bool = Field(
        False, description="Skip neighbor-joining tree reconstruction (default: False)"
    )


class ResolVIPreprocessingParameters(BaseModel):
    """ResolVI preprocessing parameters for high-resolution spatial transcriptomics

    ResolVI is a deep learning method for denoising and correcting molecular
    misassignment in cellular-resolved spatial transcriptomics data.

    Requirements:
    - High-resolution cellular-resolved spatial data (Xenium, MERFISH, CosMx, etc.)
    - Real spatial coordinates
    - Raw count data (not log-transformed)
    - Minimum ~1000 cells recommended
    """

    # Core model parameters
    n_latent: int = Field(
        default=30,
        gt=0,
        le=100,
        description="Number of latent dimensions for ResolVI representation",
    )
    n_hidden: int = Field(
        default=128,
        gt=0,
        le=512,
        description="Number of hidden units in neural network layers",
    )
    n_epochs: int = Field(
        default=100,
        gt=0,
        le=1000,
        description="Number of training epochs (100 recommended by official docs)",
    )

    # Data specifications
    spatial_key: Optional[str] = Field(
        default=None,
        description="Spatial coordinate key in obsm (auto-detected if None)",
    )
    count_layer: Optional[str] = Field(
        default=None, description="Layer containing count data (uses X if None)"
    )

    # Semi-supervised options
    labels_key: Optional[str] = Field(
        default=None, description="Cell type labels column for semi-supervised mode"
    )
    batch_key: Optional[str] = Field(
        default=None, description="Batch key for multi-sample integration"
    )
    semisupervised: bool = Field(
        default=False, description="Enable semi-supervised mode with cell type labels"
    )

    # Output options
    compute_corrected_counts: bool = Field(
        default=True, description="Compute and store molecularly corrected counts"
    )
    compute_differential_abundance: bool = Field(
        default=False, description="Enable differential abundance computation"
    )

    # Hardware configuration
    use_gpu: bool = Field(
        default=False, description="Use GPU acceleration (highly recommended for speed)"
    )

    # Data validation
    require_spatial: bool = Field(
        default=True,
        description="Require real spatial coordinates (should always be True for ResolVI)",
    )
    min_cells: int = Field(
        default=1000,
        gt=0,
        description="Minimum cells required for meaningful ResolVI analysis",
    )
