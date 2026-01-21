"""
Data models for spatial transcriptomics analysis.
"""

from __future__ import annotations

from typing import Annotated, Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field, model_validator
from typing_extensions import Self

# =============================================================================
# Spatial Data Type Definitions (Single Source of Truth)
# =============================================================================

# Actual platform types for storage - these are the only valid final states
# - visium: 10x Visium spatial transcriptomics (spot-based, ~55μm resolution)
# - xenium: 10x Xenium in situ (single-cell resolution, imaging-based)
# - slide_seq: Slide-seq/Slide-seqV2 (bead-based, ~10μm resolution)
# - merfish: MERFISH imaging-based (single-cell resolution)
# - seqfish: seqFISH/seqFISH+ (single-cell resolution)
# - generic: General spatial data (h5ad files, unknown platforms)
SpatialPlatform = Literal[
    "visium", "xenium", "slide_seq", "merfish", "seqfish", "generic"
]


class ColumnInfo(BaseModel):
    """Metadata column information for dataset profiling"""

    name: str
    dtype: Literal["categorical", "numerical"]
    n_unique: int
    sample_values: Optional[list[str]] = None  # Sample values for categorical
    range: Optional[tuple[float, float]] = None  # Value range for numerical


class SpatialDataset(BaseModel):
    """Spatial transcriptomics dataset model with comprehensive metadata profile"""

    id: str
    name: str
    data_type: SpatialPlatform  # Only valid platform types, never "auto" or "h5ad"
    description: Optional[str] = None

    # Basic statistics
    n_cells: int = 0
    n_genes: int = 0
    spatial_coordinates_available: bool = False
    tissue_image_available: bool = False

    # Metadata profiles - let LLM interpret the structure
    obs_columns: Optional[list[ColumnInfo]] = None  # Cell-level metadata
    var_columns: Optional[list[ColumnInfo]] = None  # Gene-level metadata
    obsm_keys: Optional[list[str]] = None  # Multi-dimensional data keys
    uns_keys: Optional[list[str]] = None  # Unstructured data keys

    # Gene expression profiles
    top_highly_variable_genes: Optional[list[str]] = None
    top_expressed_genes: Optional[list[str]] = None


class PreprocessingParameters(BaseModel):
    """Preprocessing parameters model"""

    # Data filtering and subsampling parameters (user controlled)
    filter_genes_min_cells: Optional[int] = Field(
        default=3, gt=0, description="Filter genes expressed in fewer than N cells."
    )
    filter_cells_min_genes: Optional[int] = Field(
        default=30, gt=0, description="Filter cells expressing fewer than N genes."
    )
    subsample_spots: Optional[int] = Field(
        default=None,
        gt=0,
        le=50000,
        description="Subsample to N spots. None = no subsampling.",
    )
    subsample_genes: Optional[int] = Field(
        default=None,
        gt=0,
        le=50000,
        description="Keep top N variable genes. None = keep all.",
    )
    subsample_random_seed: int = 42

    # ========== Mitochondrial and Ribosomal Gene Filtering ==========
    filter_mito_pct: Optional[float] = Field(
        default=20.0,
        ge=0.0,
        le=100.0,
        description="Filter spots with mitochondrial percentage above threshold. Use 30-50 for muscle/neurons.",
    )
    remove_mito_genes: bool = Field(
        default=True,
        description="Exclude mitochondrial genes from HVG selection. Set False if studying mito biology.",
    )
    remove_ribo_genes: bool = Field(
        default=False,
        description="Exclude ribosomal genes from HVG selection. Usually not needed.",
    )

    # Normalization and scaling parameters
    normalization: Literal["log", "sct", "pearson_residuals", "none", "scvi"] = Field(
        default="log",
        description="Normalization method. 'sct' requires R, 'scvi' for batch correction.",
    )
    scale: bool = Field(
        default=False,
        description="Scale to unit variance before PCA. Usually not needed for spatial data.",
    )
    n_hvgs: int = Field(
        default=2000,
        gt=0,
        le=5000,
        description="Number of highly variable genes to select.",
    )
    n_pcs: int = Field(
        default=30, gt=0, le=100, description="Number of principal components for PCA."
    )

    # ========== Normalization Control Parameters ==========
    normalize_target_sum: Optional[float] = Field(
        default=None,
        ge=1.0,
        le=1e8,
        description="Target counts per cell. None=median, 1e4=Visium, 1e6=MERFISH/Xenium.",
    )

    scale_max_value: Optional[float] = Field(
        default=10.0,
        ge=1.0,
        le=100.0,
        description="Max value for clipping after scaling (in std devs). None=no clipping.",
    )

    # SCTransform preprocessing parameters (requires R + sctransform package via rpy2)
    # Installation: R -e 'install.packages("sctransform")' && pip install rpy2
    sct_var_features_n: int = Field(
        default=3000,
        ge=100,
        le=10000,
        description="Number of variable features for SCTransform.",
    )
    sct_method: Literal["offset", "fix-slope"] = Field(
        default="fix-slope",
        description="SCTransform version. 'fix-slope' (v2, recommended) or 'offset' (v1).",
    )
    sct_exclude_poisson: bool = Field(
        default=True,
        description="Exclude Poisson genes from SCTransform regularization.",
    )
    sct_n_cells: Optional[int] = Field(
        default=5000,
        ge=100,
        description="Cells for parameter estimation. None=all cells (better for small datasets).",
    )

    # scVI preprocessing parameters - architecture
    use_scvi_preprocessing: bool = False  # Whether to use scVI for preprocessing
    scvi_n_hidden: int = 128
    scvi_n_latent: int = 10
    scvi_n_layers: int = 1
    scvi_dropout_rate: float = 0.1
    scvi_gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb"

    # scVI preprocessing parameters - training (user-configurable)
    scvi_max_epochs: Annotated[int, Field(gt=0, le=2000)] = Field(
        default=400,
        description="Training epochs for scVI. Increase to 600-800 for large datasets without early stopping.",
    )
    scvi_early_stopping: bool = Field(
        default=True,
        description="Enable early stopping. Set False for debugging or exact epoch control.",
    )
    scvi_early_stopping_patience: Annotated[int, Field(gt=0, le=100)] = Field(
        default=20,
        description="Epochs to wait before early stopping. Increase (30-50) for noisy data, decrease (10-15) for speed.",
    )
    scvi_train_size: Annotated[float, Field(gt=0.5, le=1.0)] = Field(
        default=0.9,
        description="Training data fraction (rest for validation). 1.0 disables validation and early stopping.",
    )

    # Key naming parameters (configurable hard-coded keys)
    cluster_key: str = Field(
        default="leiden",
        alias="clustering_key",
        description="Key name for clustering results in obs.",
    )
    spatial_key: Optional[str] = Field(
        default=None,
        description="Spatial coordinate key in obsm (auto-detected if None)",
    )
    batch_key: str = Field(
        default="batch", description="Key name for batch information in obs."
    )

    # User-controllable parameters (scientifically-informed defaults)
    n_neighbors: Annotated[int, Field(gt=2, le=100)] = Field(
        default=15,
        description="Neighbors for k-NN graph. Larger (20-50) for global structure, smaller (5-10) for local patterns.",
    )
    clustering_resolution: Annotated[float, Field(gt=0.1, le=2.0)] = Field(
        default=1.0,
        description="Leiden clustering resolution. Higher (1.5-2.0) for more clusters, lower (0.2-0.5) for fewer.",
    )


class DifferentialExpressionParameters(BaseModel):
    """Differential expression analysis parameters model.

    This model encapsulates all parameters for differential expression analysis,
    following the unified (data_id, ctx, params) signature pattern.
    """

    group_key: str = Field(
        ...,
        description="Grouping column in adata.obs. Common: 'leiden', 'cell_type', 'seurat_clusters'.",
    )

    group1: Optional[str] = Field(
        None,
        description="First group for comparison. None=find markers for all groups (one-vs-rest).",
    )

    group2: Optional[str] = Field(
        None,
        description="Second group. None/'rest'=compare against all others. Requires group1.",
    )

    method: Literal[
        "wilcoxon", "t-test", "t-test_overestim_var", "logreg", "pydeseq2"
    ] = Field(
        "wilcoxon",
        description="Statistical method. 'pydeseq2' requires sample_key for pseudobulk analysis.",
    )

    sample_key: Optional[str] = Field(
        None,
        description="Sample column for pseudobulk aggregation. Required for 'pydeseq2'. Common: 'sample', 'patient_id', 'batch'.",
    )

    n_top_genes: Annotated[int, Field(gt=0, le=500)] = Field(
        50,
        description="Top DE genes to return per group.",
    )

    pseudocount: Annotated[float, Field(gt=0, le=100)] = Field(
        1.0,
        description="Pseudocount for log2FC. Lower (0.1-0.5) for low-expression sensitivity, higher (1-10) for sparse data.",
    )

    min_cells: Annotated[int, Field(gt=0, le=1000)] = Field(
        3,
        description="Minimum cells per group. Increase to 10-30 for more robust results.",
    )


class VisualizationParameters(BaseModel):
    """Visualization parameters model"""

    model_config = ConfigDict(extra="forbid")  # Strict validation after preprocessing

    feature: Optional[Union[str, list[str]]] = Field(
        None,
        description="Single feature or list of features (accepts both 'feature' and 'features')",
    )  # Single feature or list of features

    @model_validator(mode="before")
    @classmethod
    def preprocess_params(cls, data):
        """
        Preprocess visualization parameters to handle different input formats.

        Handles:
        - None: Returns empty dict
        - str: Converts to feature parameter (supports "gene:CCL21" and "CCL21" formats)
        - dict: Normalizes features/feature naming
        """
        # Handle None input
        if data is None:
            return {}

        # Handle string format parameters (shorthand for feature)
        if isinstance(data, str):
            if data.startswith("gene:"):
                feature = data.split(":", 1)[1]
                return {"feature": feature, "plot_type": "feature"}
            else:
                return {"feature": data, "plot_type": "feature"}

        # Handle dict format - normalize features/feature naming
        if isinstance(data, dict):
            data_copy = data.copy()
            # Handle 'features' as alias for 'feature'
            if "features" in data_copy and "feature" not in data_copy:
                data_copy["feature"] = data_copy.pop("features")
            return data_copy

        # For other types (e.g., VisualizationParameters instances), return as-is
        return data

    # Refactored plot_type: 10 unified types (from original 19)
    # - feature: spatial/UMAP feature visualization (basis='spatial'|'umap')
    # - expression: heatmap/violin/dotplot/correlation (subtype='heatmap'|'violin'|'dotplot'|'correlation')
    # - deconvolution: cell type proportions (subtype='spatial_multi'|'pie'|'dominant'|'stacked_bar'|'imputation')
    # - communication: cell-cell communication (subtype='heatmap'|'chord'|'network'|'dotplot')
    # - interaction: spatial ligand-receptor pairs
    # - trajectory: pseudotime/fate analysis (subtype='pseudotime'|'fate'|'gene_trends')
    # - velocity: RNA velocity (subtype='stream'|'grid'|'arrow')
    # - statistics: spatial statistics (subtype='neighborhood'|'co_occurrence'|'ripley'|'moran'|'centrality'|'getis_ord')
    # - enrichment: pathway enrichment (subtype='barplot'|'dotplot'|'heatmap'|'network')
    # - cnv: copy number variation (subtype='heatmap'|'spatial')
    # - integration: batch integration quality
    plot_type: Literal[
        "feature",  # Unified spatial/UMAP feature visualization
        "expression",  # Heatmap, violin, dotplot, gene correlation
        "deconvolution",  # Cell type proportions with subtypes
        "communication",  # Cell-cell communication patterns
        "interaction",  # Spatial ligand-receptor pairs
        "trajectory",  # Pseudotime and fate analysis
        "velocity",  # RNA velocity visualization
        "statistics",  # Spatial statistics (Moran's I, etc.)
        "enrichment",  # Pathway/gene set enrichment
        "cnv",  # Copy number variation
        "integration",  # Batch integration quality
    ] = "feature"
    colormap: str = "coolwarm"

    # Unified subtype parameter for plot_types with multiple visualization modes
    subtype: Optional[str] = Field(
        None,
        description=(
            "Visualization subtype. Options by plot_type:\n"
            "- expression: 'heatmap'|'violin'|'dotplot'|'correlation'\n"
            "- deconvolution: 'spatial_multi'|'pie'|'dominant'|'stacked_bar'|'imputation'\n"
            "- communication: 'heatmap'|'chord'|'network'|'dotplot'\n"
            "- trajectory: 'pseudotime'|'fate'|'gene_trends'\n"
            "- velocity: 'stream'|'grid'|'arrow'\n"
            "- statistics: 'neighborhood'|'co_occurrence'|'ripley'|'moran'|'centrality'|'getis_ord' (required)\n"
            "- enrichment: 'barplot'|'dotplot'|'heatmap'|'network'\n"
            "- cnv: 'heatmap'|'spatial'"
        ),
    )
    cluster_key: Optional[str] = Field(
        None,
        description="Cluster/cell type column. Required for heatmap, violin, dotplot.",
    )

    # Multi-gene visualization parameters
    multi_panel: bool = False  # Whether to create multi-panel plots
    panel_layout: Optional[tuple[int, int]] = (
        None  # (rows, cols) - auto-determined if None
    )

    # GridSpec subplot spacing parameters (for multi-panel plots)
    subplot_wspace: float = Field(
        0.0,
        ge=-0.3,
        le=1.0,
        description="Horizontal subplot spacing. Lower for tighter, higher for looser.",
    )
    subplot_hspace: float = Field(
        0.3,
        ge=0.0,
        le=1.0,
        description="Vertical subplot spacing. Lower for tighter, higher for looser.",
    )

    # Colorbar parameters (for spatial plots with make_axes_locatable)
    colorbar_pad: float = Field(
        0.02,
        ge=0.0,
        le=0.2,
        description="Distance between subplot and colorbar.",
    )
    colorbar_size: str = Field(
        "3%",
        description="Colorbar width as percentage (e.g., '3%', '5%').",
    )

    # Ligand-receptor pair parameters
    lr_pairs: Optional[list[tuple[str, str]]] = None  # List of (ligand, receptor) pairs
    lr_database: str = "cellchat"  # Database for LR pairs
    plot_top_pairs: int = Field(
        6,
        gt=0,
        le=100,
        description="Top LR pairs to display. Use higher values (e.g., 50) for chord diagrams.",
    )

    # Gene correlation parameters
    correlation_method: Literal["pearson", "spearman", "kendall"] = "pearson"
    show_correlation_stats: bool = True

    # Figure parameters
    figure_size: Optional[tuple[int, int]] = (
        None  # (width, height) - auto-determined if None
    )
    dpi: int = 300  # Publication quality (Nature/Cell standard)
    alpha: float = 0.9  # Spot transparency (higher = more opaque)
    spot_size: Optional[float] = Field(
        None,
        description="Spot size in pixels. None=auto. Set manually: 50 for dense, 150 for sparse.",
    )
    alpha_img: float = Field(
        0.3,
        ge=0.0,
        le=1.0,
        description="Background image transparency. Lower=dimmer, higher to show tissue structure.",
    )
    show_tissue_image: bool = Field(
        True,
        description="Show tissue histology image as background. False for clean coordinate system.",
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

    # Coordinate basis for feature visualization
    basis: Optional[str] = Field(
        "spatial",
        description=(
            "Coordinate basis for feature/interaction visualization:\n"
            "- 'spatial': Spatial coordinates (default)\n"
            "- 'umap': UMAP coordinates\n"
            "- 'pca': PCA coordinates (trajectory only)"
        ),
    )

    # GSEA visualization parameters
    gsea_results_key: str = "gsea_results"  # Key in adata.uns for GSEA results
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
    velocity_scale: float = Field(1.0, description="Scaling factor for velocity arrows")

    # NEW: Heatmap enhancement parameters
    obs_annotation: Optional[list[str]] = Field(
        None, description="List of obs keys to show as column annotations"
    )
    var_annotation: Optional[list[str]] = Field(
        None, description="List of var keys to show as row annotations"
    )
    annotation_colors: Optional[dict[str, str]] = Field(
        None, description="Custom colors for annotations"
    )

    # NEW: Integration assessment parameters
    batch_key: str = Field(
        "batch", description="Key in adata.obs for batch/sample identifier"
    )
    integration_method: Optional[str] = Field(
        None, description="Integration method used (for display)"
    )

    # Dotplot visualization parameters
    dotplot_dendrogram: bool = Field(
        False,
        description="Whether to show dendrogram for gene clustering in dotplot",
    )
    dotplot_swap_axes: bool = Field(
        False,
        description="Swap axes to show genes on x-axis and groups on y-axis",
    )
    dotplot_standard_scale: Optional[Literal["var", "group"]] = Field(
        None,
        description="Standardize expression. 'var'=per gene, 'group'=per cluster.",
    )
    dotplot_dot_max: Optional[float] = Field(
        None,
        ge=0.0,
        le=1.0,
        description="Maximum dot size fraction. None=use observed max.",
    )
    dotplot_dot_min: Optional[float] = Field(
        None,
        ge=0.0,
        le=1.0,
        description="Minimum dot size fraction. None=use observed min.",
    )
    dotplot_smallest_dot: float = Field(
        0.0,
        ge=0.0,
        le=50.0,
        description="Dot size when fraction=0. Default 0 hides unexpressed genes.",
    )
    dotplot_var_groups: Optional[dict[str, list[str]]] = Field(
        None,
        description="Group genes by category. Example: {'T cells': ['CD3D', 'CD4']}.",
    )
    dotplot_categories_order: Optional[list[str]] = Field(
        None,
        description="Custom order for groups (clusters/cell types) on the axis",
    )

    # Deconvolution visualization parameters
    n_cell_types: Annotated[
        int,
        Field(
            gt=0,
            le=10,
            description="Top cell types to show in deconvolution.",
        ),
    ] = 4
    deconv_method: Optional[str] = Field(
        None,
        description="Deconvolution method (e.g., 'cell2location', 'rctd'). Auto-selected if only one exists.",
    )
    min_proportion_threshold: float = Field(
        0.3,
        ge=0.0,
        le=1.0,
        description="Threshold for 'pure' vs 'mixed' spots in dominant_type plot.",
    )
    show_mixed_spots: bool = Field(
        True,
        description="Mark mixed spots in dominant_type visualization.",
    )
    pie_scale: float = Field(
        0.4,
        gt=0.0,
        le=2.0,
        description="Scale factor for pie charts in scatterpie.",
    )
    scatterpie_alpha: float = Field(
        1.0,
        ge=0.0,
        le=1.0,
        description="Pie chart transparency (0=transparent, 1=opaque).",
    )
    max_spots: int = Field(
        100,
        gt=0,
        le=1000,
        description="Max spots to show in stacked_bar.",
    )
    sort_by: Literal["dominant_type", "spatial", "cluster"] = Field(
        "dominant_type",
        description="Sorting for stacked_bar: by dominant type, spatial order, or cluster.",
    )

    # Export options (unified in visualize_data, no separate save tool needed)
    output_path: Optional[str] = Field(
        None,
        description="Save path. If provided, saves to this path instead of default ./visualizations/.",
    )
    output_format: Literal["png", "pdf", "svg", "eps", "tiff", "jpg"] = Field(
        "png",
        description="Output format. 'pdf'/'svg' for vector graphics (publication).",
    )

    @model_validator(mode="after")
    def validate_conditional_parameters(self) -> Self:
        """Validate parameter dependencies and provide helpful error messages."""

        # Statistics validation - subtype is required
        if self.plot_type == "statistics" and (
            not self.subtype
            or (isinstance(self.subtype, str) and not self.subtype.strip())
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
                f"Parameter dependency error: subtype is required when plot_type='statistics'.\n"
                f"Available subtypes: {', '.join(available_types)}\n"
                f"Example usage: VisualizationParameters(plot_type='statistics', subtype='neighborhood')"
            )

        # Deconvolution - set default subtype if not provided
        if self.plot_type == "deconvolution" and not self.subtype:
            self.subtype = "spatial_multi"

        # Expression - set default subtype if not provided
        if self.plot_type == "expression" and not self.subtype:
            self.subtype = "heatmap"

        # CNV - set default subtype if not provided
        if self.plot_type == "cnv" and not self.subtype:
            self.subtype = "heatmap"

        # Velocity - set default subtype if not provided
        if self.plot_type == "velocity" and not self.subtype:
            self.subtype = "stream"

        # Enrichment - set default subtype if not provided
        if self.plot_type == "enrichment" and not self.subtype:
            self.subtype = "barplot"

        # Communication - set default subtype if not provided
        if self.plot_type == "communication" and not self.subtype:
            self.subtype = "heatmap"

        # Trajectory - set default subtype if not provided
        if self.plot_type == "trajectory" and not self.subtype:
            self.subtype = "pseudotime"

        return self


class AnnotationParameters(BaseModel):
    """Cell type annotation parameters model"""

    method: Literal[
        "tangram", "scanvi", "cellassign", "mllmcelltype", "sctype", "singler"
    ] = Field(
        default="tangram",
        description="tangram/scanvi: require reference_data_id. cellassign: requires marker_genes. singler: uses celldex reference.",
    )
    marker_genes: Optional[dict[str, list[str]]] = None
    reference_data: Optional[str] = None
    reference_data_id: Optional[str] = (
        None  # For Tangram method - ID of reference single-cell dataset
    )
    training_genes: Optional[list[str]] = (
        None  # For Tangram method - genes to use for mapping
    )
    num_epochs: int = (
        100  # For Tangram/ScanVI methods - number of training epochs (reduced for faster training)
    )
    tangram_mode: Literal["cells", "clusters"] = (
        "cells"  # Tangram mapping mode: 'cells' (cell-level) or 'clusters' (cluster-level)
    )
    cluster_label: Optional[str] = (
        None  # For mLLMCellType method - cluster label in spatial data. Only required when method='mllmcelltype'
    )
    cell_type_key: Optional[str] = Field(
        default=None,
        description="Cell type column in reference data. Required for tangram, scanvi, singler.",
    )

    # Tangram-specific parameters (aligned with scvi.external.Tangram API)
    tangram_density_prior: Literal["rna_count_based", "uniform"] = (
        "rna_count_based"  # Density prior for mapping
    )
    tangram_device: str = "cpu"  # Device for computation ('cpu' or 'cuda:0')
    tangram_learning_rate: float = 0.1  # Learning rate for optimization
    tangram_compute_validation: bool = False  # Whether to compute validation metrics
    tangram_project_genes: bool = False  # Whether to project gene expression

    # Tangram regularization parameters (optional)
    tangram_lambda_r: Optional[float] = (
        None  # Regularization parameter for entropy term in Tangram loss
    )
    tangram_lambda_neighborhood: Optional[float] = (
        None  # Neighborhood regularization parameter for spatial smoothness
    )

    # General parameters for batch effect and data handling
    batch_key: Optional[str] = None  # For batch effect correction
    layer: Optional[str] = None  # Which layer to use for analysis

    # scANVI parameters (scvi-tools semi-supervised label transfer)
    scanvi_n_hidden: int = Field(
        default=128,
        description="Hidden units per layer.",
    )
    scanvi_n_latent: int = Field(
        default=10,
        description="Latent space dimensions. Use 3-5 for small datasets to avoid NaN.",
    )
    scanvi_n_layers: int = Field(
        default=1,
        description="Hidden layers. Use 2 for large integration.",
    )
    scanvi_dropout_rate: float = Field(
        default=0.1,
        description="Dropout rate. Use 0.2-0.3 for small datasets.",
    )
    scanvi_unlabeled_category: str = Field(
        default="Unknown",
        description="Label for unlabeled cells.",
    )

    # SCVI pretraining parameters
    scanvi_use_scvi_pretrain: bool = Field(
        default=True,
        description="SCVI pretraining. Set False if NaN errors on small datasets.",
    )
    scanvi_scvi_epochs: int = Field(default=200, description="SCVI pretraining epochs.")
    scanvi_scanvi_epochs: int = Field(
        default=20,
        description="SCANVI training epochs after pretraining. Increase to 50-100 for complex datasets.",
    )
    scanvi_n_samples_per_label: int = Field(
        default=100,
        description="Samples per label for semi-supervised training.",
    )

    # Query training parameters
    scanvi_query_epochs: int = Field(
        default=100,
        description="Query training epochs. Use 50 for small datasets.",
    )
    scanvi_check_val_every_n_epoch: int = Field(
        default=10, description="Validation check frequency."
    )

    # CellAssign parameters
    cellassign_n_hidden: int = 100
    cellassign_learning_rate: float = 0.001
    cellassign_max_iter: int = 200

    # mLLMCellType parameters
    mllm_n_marker_genes: int = Field(
        default=20,
        gt=0,
        le=50,
        description="Number of marker genes per cluster for LLM annotation.",
    )
    mllm_species: Literal["human", "mouse"] = "human"
    mllm_tissue: Optional[str] = None
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
    ] = "openai"
    mllm_model: Optional[str] = None
    mllm_api_key: Optional[str] = None
    mllm_additional_context: Optional[str] = None
    mllm_use_cache: bool = True
    mllm_base_urls: Optional[Union[str, dict[str, str]]] = None
    mllm_verbose: bool = False
    mllm_force_rerun: bool = False

    # Multi-model consensus parameters (interactive_consensus_annotation)
    mllm_use_consensus: bool = False  # Whether to use multi-model consensus
    mllm_models: Optional[list[Union[str, dict[str, str]]]] = (
        None  # List of models for consensus
    )
    mllm_api_keys: Optional[dict[str, str]] = None  # Dict mapping provider to API key
    mllm_consensus_threshold: float = 0.7  # Agreement threshold for consensus
    mllm_entropy_threshold: float = 1.0  # Entropy threshold for controversy detection
    mllm_max_discussion_rounds: int = 3  # Maximum discussion rounds
    mllm_consensus_model: Optional[Union[str, dict[str, str]]] = (
        None  # Model for consensus checking
    )
    mllm_clusters_to_analyze: Optional[list[str]] = None  # Specific clusters to analyze

    # ScType parameters
    sctype_tissue: Optional[str] = (
        None  # Tissue type (supported: "Adrenal", "Brain", "Eye", "Heart", "Hippocampus", "Immune system", "Intestine", "Kidney", "Liver", "Lung", "Muscle", "Pancreas", "Placenta", "Spleen", "Stomach", "Thymus")
    )
    sctype_db_: Optional[str] = (
        None  # Custom database path (if None, uses default ScTypeDB)
    )
    sctype_scaled: bool = True  # Whether input data is scaled
    sctype_custom_markers: Optional[dict[str, dict[str, list[str]]]] = (
        None  # Custom markers: {"CellType": {"positive": [...], "negative": [...]}}
    )
    sctype_use_cache: bool = True  # Whether to cache results to avoid repeated R calls

    # SingleR parameters (for enhanced marker_genes method)
    singler_reference: Optional[str] = Field(
        default=None,
        description="Celldex reference. Human: 'hpca' (default). Mouse: 'immgen' (default).",
    )
    singler_integrated: bool = Field(
        default=False,
        description="Whether to use integrated annotation with multiple references",
    )
    singler_fine_tune: bool = Field(
        default=True,
        description="Whether to perform fine-tuning step in SingleR annotation (refines labels based on marker genes)",
    )
    num_threads: int = 4  # Number of threads for parallel processing


class SpatialStatisticsParameters(BaseModel):
    """Spatial statistics parameters model"""

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
        description="Cluster column. Required for neighborhood, co_occurrence, ripley, join_count analyses.",
    )
    n_neighbors: Annotated[int, Field(gt=0)] = Field(
        8,
        description="Nearest neighbors for spatial graph.",
    )

    # Unified gene selection parameter (NEW)
    genes: Optional[list[str]] = Field(
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
        description="Parallel jobs. 1=none, None=auto, -1=all cores.",
    )
    backend: Literal["loky", "threading", "multiprocessing"] = Field(
        "threading",
        description="Parallelization backend.",
    )

    # Moran's I specific parameters
    moran_n_perms: Annotated[int, Field(gt=0, le=10000)] = Field(
        10,
        description="Permutations. Use 100+ for publication.",
    )
    moran_two_tailed: bool = Field(False, description="Use two-tailed test.")

    # Local Moran's I (LISA) specific parameters
    local_moran_permutations: Annotated[int, Field(gt=0, le=9999)] = Field(
        999,
        description="Permutations for p-value. Use 9999 for publication quality.",
    )
    local_moran_alpha: Annotated[float, Field(gt=0.0, lt=1.0)] = Field(
        0.05,
        description="Significance level. 0.01 (conservative), 0.05 (standard), 0.10 (exploratory).",
    )
    local_moran_fdr_correction: bool = Field(
        True,
        description="Apply FDR correction. Set False for exploratory analysis.",
    )

    # Getis-Ord Gi* specific parameters
    getis_ord_correction: Literal["bonferroni", "fdr_bh", "none"] = Field(
        "fdr_bh",
        description="Multiple testing correction. 'fdr_bh' (recommended), 'bonferroni', 'none'.",
    )
    getis_ord_alpha: Annotated[float, Field(gt=0.0, le=1.0)] = Field(
        0.05,
        description="Significance level for hotspot detection.",
    )

    # Co-occurrence specific parameters
    co_occurrence_interval: Optional[int] = Field(
        50,
        description="Distance intervals for co-occurrence. Lower=faster, higher=finer resolution.",
    )

    # Bivariate Moran's I specific parameters
    gene_pairs: Optional[list[tuple[str, str]]] = Field(
        None, description="Gene pairs for bivariate analysis"
    )


class RNAVelocityParameters(BaseModel):
    """RNA velocity analysis parameters model"""

    model_config = ConfigDict(
        extra="forbid"
    )  # Strict validation - no extra parameters allowed

    # Velocity computation method selection
    method: Literal["scvelo", "velovi"] = Field(
        default="scvelo",
        description="REQUIRES 'spliced' and 'unspliced' layers (from velocyto/kallisto/STARsolo).",
    )

    # scVelo specific parameters
    scvelo_mode: Literal["deterministic", "stochastic", "dynamical"] = Field(
        default="stochastic",
        description="'dynamical' mode REQUIRED for CellRank gene_trends visualization.",
    )
    n_pcs: int = Field(
        default=30,
        gt=0,
        le=100,
        description="Principal components for velocity computation.",
    )
    basis: str = "spatial"

    # Preprocessing parameters for velocity computation
    min_shared_counts: int = Field(
        default=30, gt=0, description="Minimum shared counts for gene filtering."
    )
    n_top_genes: int = Field(
        default=2000, gt=0, description="Number of top variable genes to retain."
    )
    n_neighbors: int = Field(
        default=30, gt=0, description="Neighbors for moments computation."
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

    method: Literal["cellrank", "palantir", "dpt"] = Field(
        default="cellrank",
        description="'cellrank' requires velocity data. 'palantir'/'dpt' work without velocity.",
    )
    spatial_weight: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Weight for spatial kernel in CellRank (0=no spatial, 1=spatial only). Only used with method='cellrank'.",
    )
    root_cells: Optional[list[str]] = Field(
        default=None,
        description="Starting cells for trajectory inference (palantir/dpt).",
    )

    # CellRank specific parameters
    cellrank_kernel_weights: tuple[float, float] = (0.8, 0.2)
    cellrank_n_states: int = Field(
        default=5, gt=0, le=20, description="Number of macrostates for CellRank."
    )

    # Palantir specific parameters
    palantir_n_diffusion_components: int = Field(
        default=10, gt=0, le=50, description="Diffusion components for Palantir."
    )
    palantir_num_waypoints: int = Field(
        default=500, gt=0, description="Number of waypoints for Palantir."
    )

    # Fallback control
    # Removed: allow_fallback_to_dpt - No longer doing automatic fallbacks
    # LLMs should explicitly choose which method to use


class IntegrationParameters(BaseModel):
    """Sample integration parameters model"""

    method: Literal["harmony", "bbknn", "scanorama", "scvi"] = Field(
        default="harmony",
        description="Integration method. 'scvi' requires scvi-tools and supports GPU.",
    )
    batch_key: str = "batch"
    n_pcs: int = Field(
        default=30, gt=0, le=100, description="Principal components for integration."
    )
    align_spatial: bool = True
    reference_batch: Optional[str] = None

    # Common scvi-tools parameters
    use_gpu: bool = False
    n_epochs: Optional[int] = None

    # scVI integration parameters
    scvi_n_hidden: int = 128
    scvi_n_latent: int = 10
    scvi_n_layers: int = 1
    scvi_dropout_rate: float = 0.1
    scvi_gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb"


class DeconvolutionParameters(BaseModel):
    """Spatial deconvolution parameters model"""

    method: Literal[
        "flashdeconv",
        "cell2location",
        "rctd",
        "destvi",
        "stereoscope",
        "spotlight",
        "tangram",
        "card",
    ] = Field(
        default="flashdeconv",
        description="All methods require reference_data_id and cell_type_key.",
    )

    reference_data_id: Optional[str] = Field(
        default=None,
        description="REQUIRED: ID of loaded single-cell reference dataset.",
    )

    cell_type_key: str = Field(
        description="REQUIRED: Column in reference data with cell type labels.",
    )

    # Universal GPU parameter
    use_gpu: bool = Field(
        False,
        description="GPU acceleration. Supported by cell2location, destvi, stereoscope, tangram.",
    )

    # Cell2location specific parameters
    cell2location_ref_model_epochs: Annotated[int, Field(gt=0)] = Field(
        250,
        description="Reference model training epochs. Cell2location only.",
    )
    cell2location_n_epochs: Annotated[int, Field(gt=0)] = Field(
        30000,
        description="Spatial mapping model epochs. Cell2location only.",
    )
    cell2location_n_cells_per_spot: Annotated[int, Field(gt=0)] = Field(
        30,
        description="Expected cells per spot. 30 for Visium, 5-10 for MERFISH. Cell2location only.",
    )
    cell2location_detection_alpha: Annotated[float, Field(gt=0)] = Field(
        20.0,
        description="RNA detection sensitivity. 20 for high variability, 200 for low. Cell2location only.",
    )

    # Batch and covariate correction for cell2location
    cell2location_batch_key: Optional[str] = Field(
        None,
        description="Batch column for batch effect correction. Cell2location only.",
    )
    cell2location_categorical_covariate_keys: Optional[list[str]] = Field(
        None,
        description="Categorical covariate columns for batch correction. Cell2location only.",
    )

    # Gene filtering parameters (Cell2location-specific preprocessing)
    cell2location_apply_gene_filtering: bool = Field(
        True,
        description="Apply permissive gene filtering (different from HVG). Cell2location only.",
    )
    cell2location_gene_filter_cell_count_cutoff: int = Field(
        5,
        description="Minimum cells expressing a gene. Cell2location only.",
    )
    cell2location_gene_filter_cell_percentage_cutoff2: float = Field(
        0.03,
        description="Minimum cell percentage for gene filtering. Cell2location only.",
    )
    cell2location_gene_filter_nonz_mean_cutoff: float = Field(
        1.12,
        description="Minimum non-zero mean expression. Cell2location only.",
    )

    # Phase 2: Training enhancement parameters (Cell2location)
    cell2location_ref_model_lr: Annotated[float, Field(gt=0)] = Field(
        0.002,
        description="Reference model learning rate. Cell2location only.",
    )
    cell2location_lr: Annotated[float, Field(gt=0)] = Field(
        0.005,
        description="Model learning rate. Cell2location only.",
    )
    cell2location_ref_model_train_size: Annotated[float, Field(gt=0, le=1)] = Field(
        1.0,
        description="Fraction of reference data for training. Cell2location only.",
    )
    cell2location_train_size: Annotated[float, Field(gt=0, le=1)] = Field(
        1.0,
        description="Fraction of spatial data for training. Cell2location only.",
    )
    cell2location_enable_qc_plots: bool = Field(
        False,
        description="Generate QC diagnostic plots. Cell2location only.",
    )
    cell2location_qc_output_dir: Optional[str] = Field(
        None,
        description="Output directory for QC plots. Cell2location only.",
    )

    # Phase 3: Runtime optimization parameters (Cell2location)
    cell2location_early_stopping: bool = Field(
        False,
        description="Enable early stopping. Cell2location only.",
    )
    cell2location_early_stopping_patience: Annotated[int, Field(gt=0)] = Field(
        45,
        description="Epochs to wait before stopping. Cell2location only.",
    )
    cell2location_early_stopping_threshold: Annotated[float, Field(gt=0)] = Field(
        0.0,
        description="Minimum relative change for improvement. Cell2location only.",
    )
    cell2location_use_aggressive_training: bool = Field(
        False,
        description="Use aggressive training for >50k locations. Cell2location only.",
    )
    cell2location_validation_size: Annotated[float, Field(gt=0, lt=1)] = Field(
        0.1,
        description="Validation set fraction. Required if early_stopping=True. Cell2location only.",
    )

    # SPOTlight specific parameters
    spotlight_n_top_genes: Annotated[int, Field(gt=0, le=5000)] = Field(
        2000,
        description="Number of HVGs for deconvolution. SPOTlight only.",
    )
    spotlight_nmf_model: Literal["ns"] = Field(
        "ns",
        description="NMF model type (only 'ns' supported). SPOTlight only.",
    )
    spotlight_min_prop: Annotated[float, Field(ge=0, le=1)] = Field(
        0.01,
        description="Minimum cell type proportion threshold. SPOTlight only.",
    )
    spotlight_scale: bool = Field(
        True,
        description="Scale/normalize data. SPOTlight only.",
    )
    spotlight_weight_id: str = Field(
        "mean.AUC",
        description="Marker gene weight column. SPOTlight only.",
    )

    # DestVI parameters
    destvi_n_epochs: Annotated[int, Field(gt=0)] = Field(
        2000,
        description="Training epochs. DestVI only.",
    )
    destvi_n_hidden: int = 128
    destvi_n_latent: int = 10
    destvi_n_layers: int = 1
    destvi_dropout_rate: float = 0.1
    destvi_learning_rate: float = 1e-3

    # DestVI advanced parameters (official scvi-tools defaults)
    destvi_train_size: Annotated[float, Field(gt=0.0, le=1.0)] = Field(
        default=0.9,
        description="Training data fraction. DestVI only.",
    )
    destvi_vamp_prior_p: Annotated[int, Field(ge=1)] = Field(
        default=15,
        description="Number of VampPrior components. DestVI only.",
    )
    destvi_l1_reg: Annotated[float, Field(ge=0.0)] = Field(
        default=10.0,
        description="L1 regularization for sparsity. DestVI only.",
    )

    # Stereoscope parameters
    stereoscope_n_epochs: int = 150000
    stereoscope_learning_rate: float = 0.01
    stereoscope_batch_size: int = 128

    # RCTD specific parameters
    rctd_mode: Literal["full", "doublet", "multi"] = Field(
        "full",
        description="'doublet' for high-res (Slide-seq, MERFISH), 'full' for low-res (Visium), 'multi' for constrained mixing.",
    )
    max_cores: int = Field(
        default=4, gt=0, le=16, description="Maximum CPU cores for R-based methods."
    )
    rctd_confidence_threshold: float = Field(
        default=10.0,
        gt=0,
        description="Cell type assignment confidence (higher = more stringent).",
    )
    rctd_doublet_threshold: float = Field(
        default=25.0,
        gt=0,
        description="Doublet detection threshold (doublet/multi modes).",
    )
    rctd_max_multi_types: Annotated[int, Field(ge=2, le=10)] = Field(
        4,
        description="Max cell types per spot in multi mode. 4-6 for Visium, 2-3 for higher resolution.",
    )

    # CARD specific parameters
    card_minCountGene: Annotated[int, Field(gt=0)] = Field(
        100,
        description="Minimum total counts per gene. CARD only.",
    )
    card_minCountSpot: Annotated[int, Field(gt=0)] = Field(
        5,
        description="Minimum spots expressing a gene. CARD only.",
    )
    card_sample_key: Optional[str] = Field(
        None,
        description="Sample/batch column for multi-sample analysis. CARD only.",
    )
    card_imputation: bool = Field(
        False,
        description="Enable spatial imputation for higher resolution. CARD only.",
    )
    card_NumGrids: Annotated[int, Field(gt=0)] = Field(
        2000,
        description="Grid points for imputation. 2000 (standard), 5000 (high-res), 10000 (ultra). CARD only.",
    )
    card_ineibor: Annotated[int, Field(gt=0)] = Field(
        10,
        description="Neighbors for imputation smoothing. Higher = smoother. CARD only.",
    )

    # Tangram specific parameters
    tangram_n_epochs: Annotated[int, Field(gt=0)] = Field(
        1000,
        description="Spatial mapping epochs. Tangram only.",
    )
    tangram_mode: Literal["cells", "clusters", "constrained"] = Field(
        "cells",
        description="'cells' for cell-level, 'clusters' requires cluster_label. Tangram only.",
    )
    tangram_learning_rate: Annotated[float, Field(gt=0)] = Field(
        0.1,
        description="Optimizer learning rate. Tangram only.",
    )
    tangram_density_prior: Literal["rna_count_based", "uniform"] = Field(
        "rna_count_based",
        description="'rna_count_based' weights by RNA counts, 'uniform' for equal weights. Tangram only.",
    )

    # FlashDeconv specific parameters (DEFAULT METHOD - ultra-fast, atlas-scale)
    flashdeconv_sketch_dim: Annotated[int, Field(gt=0, le=2048)] = Field(
        512,
        description="Sketched space dimension. FlashDeconv only.",
    )
    flashdeconv_lambda_spatial: Annotated[float, Field(gt=0)] = Field(
        5000.0,
        description="Spatial regularization. 5000 for Visium, 50000+ for high-res. FlashDeconv only.",
    )
    flashdeconv_n_hvg: Annotated[int, Field(gt=0, le=5000)] = Field(
        2000,
        description="Number of HVGs. FlashDeconv only.",
    )
    flashdeconv_n_markers_per_type: Annotated[int, Field(gt=0, le=500)] = Field(
        50,
        description="Marker genes per cell type. FlashDeconv only.",
    )


class SpatialDomainParameters(BaseModel):
    """Spatial domain identification parameters model"""

    method: Literal["spagcn", "leiden", "louvain", "stagate", "graphst"] = Field(
        default="spagcn",
        description="'spagcn' uses histology image. 'stagate'/'graphst' for data without images.",
    )
    n_domains: int = Field(
        default=7, gt=0, le=50, description="Number of spatial domains to identify."
    )

    # SpaGCN specific parameters
    spagcn_s: float = Field(
        default=1.0, gt=0.0, description="Histology weight in SpaGCN."
    )
    spagcn_b: int = Field(
        default=49, gt=0, description="Spot area for color intensity extraction."
    )
    spagcn_p: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Neighborhood expression contribution fraction.",
    )
    spagcn_use_histology: bool = True
    spagcn_random_seed: int = 100

    # General clustering parameters
    resolution: float = 0.5
    use_highly_variable: bool = True
    refine_domains: bool = True
    refinement_threshold: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Neighbor agreement threshold for domain refinement.",
    )

    # Clustering-specific parameters for leiden/louvain methods
    cluster_n_neighbors: Optional[int] = Field(
        default=None, gt=0, description="Neighbors for clustering (leiden/louvain)."
    )
    cluster_spatial_weight: Optional[float] = Field(
        default=None,
        ge=0.0,
        le=1.0,
        description="Spatial information weight (leiden/louvain).",
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
    method: Literal["spatialde", "sparkx"] = Field(
        default="sparkx",
        description="'sparkx' is faster and recommended. 'spatialde' uses Gaussian process.",
    )

    # Common parameters for all methods
    n_top_genes: Optional[int] = Field(
        default=None,
        gt=0,
        le=5000,
        description="Top spatial variable genes to return. None = all significant.",
    )
    spatial_key: str = "spatial"

    # SpatialDE-specific parameters
    spatialde_normalized: bool = True
    spatialde_kernel: str = "SE"
    spatialde_pi0: Optional[float] = Field(
        default=None,
        gt=0.0,
        le=1.0,
        description="Null hypothesis prior for FDR. None = adaptive (recommended). Higher = more conservative. SpatialDE only.",
    )

    # SPARK-X specific parameters
    sparkx_percentage: float = Field(
        default=0.1, gt=0.0, le=1.0, description="Expression filtering percentage."
    )
    sparkx_min_total_counts: int = Field(
        default=10, gt=0, description="Minimum total counts per gene."
    )
    sparkx_num_core: int = Field(
        default=1, gt=0, le=16, description="CPU cores for parallel processing."
    )
    sparkx_option: Literal["single", "mixture"] = "mixture"
    sparkx_verbose: bool = False

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
    method: Literal["liana", "cellphonedb", "cellchat_r", "fastccc"] = Field(
        default="fastccc",
        description="fastccc/cellphonedb: human only. liana/cellchat_r: human, mouse, zebrafish.",
    )

    # ========== Species and Resource Control ==========
    species: Literal["human", "mouse", "zebrafish"] = Field(
        description="REQUIRED. fastccc/cellphonedb: human only. mouse/zebrafish: use liana or cellchat_r.",
    )

    # LIANA resource selection (matches actual LIANA+ supported resources)
    liana_resource: Literal[
        "consensus",  # Default: consensus of multiple databases (recommended)
        "mouseconsensus",  # Mouse consensus database
        "baccin2019",  # Baccin et al. 2019 resource
        "cellcall",  # CellCall database
        "cellchatdb",  # CellChat database
        "cellinker",  # CellLinker database
        "cellphonedb",  # CellPhoneDB database (curated, stringent)
        "celltalkdb",  # CellTalkDB database (large)
        "connectomedb2020",  # Connectome database 2020
        "embrace",  # EMBRACE database
        "guide2pharma",  # Guide to Pharmacology
        "hpmr",  # Human Plasma Membrane Receptome
        "icellnet",  # iCellNet database (immune focus)
        "italk",  # iTALK database
        "kirouac2010",  # Kirouac et al. 2010
        "lrdb",  # LRdb database
        "ramilowski2015",  # Ramilowski et al. 2015
    ] = "consensus"  # LR database resource

    # ========== Spatial Analysis Control ==========
    perform_spatial_analysis: bool = (
        True  # Whether to perform spatial bivariate analysis
    )

    # ========== Cell Type Control ==========
    # Cell type key (unified naming with other tools)
    cell_type_key: str  # REQUIRED: Which column to use for cell types. LLM will infer from metadata. Common values: 'cell_type', 'celltype', 'leiden', 'louvain', 'seurat_clusters'

    # ========== LIANA Specific Parameters ==========
    liana_local_metric: Literal["cosine", "pearson", "spearman", "jaccard"] = "cosine"
    liana_global_metric: Literal["morans", "lee"] = "morans"
    liana_n_perms: int = Field(
        default=1000, gt=0, description="Permutations for p-value calculation."
    )
    liana_nz_prop: float = Field(
        default=0.2, gt=0.0, le=1.0, description="Minimum expression proportion."
    )
    liana_bandwidth: Optional[int] = None
    liana_cutoff: float = Field(
        default=0.1, gt=0.0, le=1.0, description="Spatial connectivity cutoff."
    )
    liana_significance_alpha: float = Field(
        default=0.05,
        gt=0.0,
        lt=1.0,
        description="FDR significance threshold. 0.01 for stringent, 0.10 for exploratory. LIANA only.",
    )

    # ========== Expression Filtering Parameters ==========
    min_cells: int = Field(
        default=3, ge=0, description="Minimum cells expressing ligand or receptor."
    )

    # ========== Result Control ==========
    plot_top_pairs: int = Field(
        default=6,
        gt=0,
        le=100,
        description="Top LR pairs to display. Use higher values (e.g., 50) for chord diagrams.",
    )

    # ========== CellPhoneDB Specific Parameters ==========
    cellphonedb_threshold: float = Field(
        default=0.1, gt=0.0, le=1.0, description="Expression threshold."
    )
    cellphonedb_iterations: int = Field(
        default=1000, gt=0, le=10000, description="Statistical permutations."
    )
    cellphonedb_result_precision: int = Field(
        default=3, gt=0, le=5, description="Result decimal precision."
    )
    cellphonedb_pvalue: float = Field(
        default=0.05, gt=0.0, le=1.0, description="P-value significance threshold."
    )
    cellphonedb_use_microenvironments: bool = True
    cellphonedb_spatial_radius: Optional[float] = Field(
        default=None, gt=0.0, description="Spatial radius for microenvironments."
    )
    cellphonedb_debug_seed: Optional[int] = None

    # Multiple testing correction for CellPhoneDB
    # When using minimum p-value across multiple cell type pairs, correction is needed
    # to control false positive rate (e.g., 7 clusters = 49 pairs → FPR 91.9% without correction)
    cellphonedb_correction_method: Literal["fdr_bh", "bonferroni", "sidak", "none"] = (
        "fdr_bh"  # Multiple testing correction method (default: Benjamini-Hochberg FDR)
    )
    # Options:
    # - "fdr_bh": Benjamini-Hochberg FDR (recommended, balances sensitivity & specificity)
    # - "bonferroni": Bonferroni correction (most conservative, controls FWER)
    # - "sidak": Šidák correction (similar to Bonferroni but more accurate for independent tests)
    # - "none": No correction (NOT recommended, leads to ~92% FPR with 7 clusters)

    # ========== CellChat R Specific Parameters ==========
    # These parameters are only used when method="cellchat_r"
    cellchat_db_category: Literal[
        "Secreted Signaling",
        "ECM-Receptor",
        "Cell-Cell Contact",
        "All",
    ] = "All"
    # CellChatDB category to use:
    # - "Secreted Signaling": Ligand-receptor pairs for secreted signaling
    # - "ECM-Receptor": Extracellular matrix-receptor interactions
    # - "Cell-Cell Contact": Direct cell-cell contact interactions
    # - "All": Use all categories (default)

    cellchat_type: Literal["triMean", "truncatedMean", "thresholdedMean", "median"] = (
        "triMean"
    )
    # CellChat expression aggregation method:
    # - "trimean": Tukey's trimean (robust, default, produces fewer interactions)
    # - "truncatedMean": Truncated mean (more interactions, use with trim parameter)

    cellchat_trim: float = Field(
        default=0.1,
        ge=0.0,
        le=0.5,
        description="Trim proportion for truncatedMean method.",
    )

    cellchat_population_size: bool = True

    cellchat_min_cells: int = Field(
        default=10, ge=1, description="Minimum cells per group for filterCommunication."
    )

    cellchat_distance_use: bool = True

    cellchat_interaction_range: float = Field(
        default=250.0, gt=0.0, description="Max ligand diffusion range in microns."
    )

    cellchat_scale_distance: float = Field(
        default=0.01, gt=0.0, description="Scale factor for distance calculation."
    )

    cellchat_contact_knn_k: int = Field(
        default=6,
        ge=1,
        description="Nearest neighbors for contact-dependent signaling.",
    )

    cellchat_contact_range: Optional[float] = Field(
        default=None,
        gt=0.0,
        description="Distance threshold for contact signaling. None = use contact_knn_k.",
    )

    # CellChat spatial conversion factors (platform-specific)
    cellchat_pixel_ratio: Annotated[float, Field(gt=0.0)] = Field(
        default=0.5,
        description="Pixel to micrometer ratio. 0.5 for Visium, 0.18 for CosMx. CellChat only.",
    )

    cellchat_spatial_tol: Annotated[float, Field(gt=0.0)] = Field(
        default=27.5,
        description="Spatial tolerance (half spot diameter). 27.5 for Visium, 5-10 for single-cell. CellChat only.",
    )

    # ========== FastCCC Specific Parameters ==========
    # FastCCC is a permutation-free framework using FFT-based convolution
    # Reference: Nature Communications 2025 (https://github.com/Svvord/FastCCC)
    # Key advantage: Ultra-fast (16M cells in minutes vs hours for permutation methods)

    fastccc_single_unit_summary: Literal["Mean", "Median", "Q3", "Quantile_0.9"] = (
        Field(
            default="Mean",
            description="Gene expression aggregation. 'Median' for outlier robustness. FastCCC only.",
        )
    )

    fastccc_complex_aggregation: Literal["Minimum", "Average"] = Field(
        default="Minimum",
        description="Complex subunit aggregation. 'Minimum' ensures all subunits present. FastCCC only.",
    )

    fastccc_lr_combination: Literal["Arithmetic", "Geometric"] = Field(
        default="Arithmetic",
        description="Ligand-receptor score combination. 'Geometric' is more conservative. FastCCC only.",
    )

    fastccc_min_percentile: Annotated[float, Field(ge=0.0, le=1.0)] = Field(
        default=0.1,
        description="Minimum expression percentile. FastCCC only.",
    )

    fastccc_use_cauchy: bool = Field(
        default=True,
        description="Use Cauchy combination for robust p-values (slower). FastCCC only.",
    )

    fastccc_pvalue_threshold: Annotated[float, Field(gt=0.0, le=1.0)] = Field(
        default=0.05,
        description="P-value threshold for significance. FastCCC only.",
    )

    fastccc_use_deg: bool = Field(
        default=False,
        description="Filter to DEGs only. FastCCC only.",
    )


class EnrichmentParameters(BaseModel):
    """Parameters for gene set enrichment analysis"""

    model_config = ConfigDict(extra="forbid")

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
    ] = Field(
        default="spatial_enrichmap",
        description="'spatial_enrichmap' for spatial patterns. 'pathway_gsea'/'pathway_ora' for standard enrichment.",
    )

    # Gene sets
    gene_sets: Optional[Union[list[str], dict[str, list[str]]]] = (
        None  # Gene sets to analyze
    )
    score_keys: Optional[Union[str, list[str]]] = None  # Names for gene signatures

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
    spatial_key: str = "spatial"
    n_neighbors: int = Field(
        default=6, gt=0, description="Spatial neighbors for enrichment mapping."
    )
    smoothing: bool = True
    correct_spatial_covariates: bool = True

    # Analysis parameters
    batch_key: Optional[str] = None
    min_genes: int = Field(default=10, gt=0, description="Minimum genes in gene set.")
    max_genes: int = Field(default=500, gt=0, description="Maximum genes in gene set.")

    # Statistical parameters
    pvalue_cutoff: float = Field(
        default=0.05, gt=0.0, lt=1.0, description="P-value significance cutoff."
    )
    adjust_method: Literal["bonferroni", "fdr", "none"] = "fdr"
    n_permutations: int = Field(
        default=1000, gt=0, description="Permutations for GSEA."
    )


class CNVParameters(BaseModel):
    """Copy Number Variation (CNV) analysis parameters model"""

    # Method selection
    method: Literal["infercnvpy", "numbat"] = Field(
        "infercnvpy",
        description="'infercnvpy' for expression-based, 'numbat' requires allele data.",
    )

    # Reference cell specification
    reference_key: str = Field(
        ...,
        description="Column with cell type labels for reference cells.",
    )
    reference_categories: list[str] = Field(
        ...,
        description="Cell types to use as reference (normal) cells.",
    )

    # infercnvpy parameters
    window_size: Annotated[int, Field(gt=0, le=500)] = Field(
        100, description="Genes per CNV window."
    )
    step: Annotated[int, Field(gt=0, le=100)] = Field(
        10, description="Sliding window step size."
    )

    # Analysis options
    exclude_chromosomes: Optional[list[str]] = Field(
        None,
        description="Chromosomes to exclude (e.g., ['chrX', 'chrY']).",
    )
    dynamic_threshold: Optional[float] = Field(
        1.5,
        gt=0.0,
        description="Threshold for CNV calling.",
    )

    # Clustering and visualization options (infercnvpy)
    cluster_cells: bool = Field(False, description="Cluster cells by CNV pattern.")
    dendrogram: bool = Field(False, description="Compute hierarchical clustering.")

    # Numbat-specific parameters
    numbat_genome: Literal["hg38", "hg19", "mm10", "mm39"] = Field(
        "hg38", description="Reference genome. Numbat only."
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
        description="Max entropy threshold. 0.8 for spatial, 0.5 for scRNA-seq. Numbat only.",
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


class RegistrationParameters(BaseModel):
    """Spatial registration parameters for aligning multiple tissue slices."""

    method: Literal["paste", "stalign"] = Field(
        "paste",
        description="'paste' for optimal transport, 'stalign' for complex deformations.",
    )
    reference_idx: Optional[int] = Field(
        None,
        ge=0,
        description="Reference slice index (0-indexed). None = first slice.",
    )

    # PASTE-specific parameters
    paste_alpha: Annotated[float, Field(gt=0, le=1)] = Field(
        0.1,
        description="Spatial vs expression weight. Higher = more spatial. PASTE only.",
    )
    paste_n_components: Annotated[int, Field(gt=0, le=100)] = Field(
        30,
        description="Number of PCA components for PASTE center alignment (default: 30).",
    )
    paste_numItermax: Annotated[int, Field(gt=0, le=1000)] = Field(
        200,
        description="Maximum iterations for optimal transport solver (default: 200).",
    )

    # STalign-specific parameters
    stalign_image_size: tuple[int, int] = Field(
        (128, 128),
        description="Image size for STalign rasterization (height, width).",
    )
    stalign_niter: Annotated[int, Field(gt=0, le=500)] = Field(
        50,
        description="Number of LDDMM iterations for STalign (default: 50).",
    )
    stalign_a: Annotated[float, Field(gt=0)] = Field(
        500.0,
        description="Regularization parameter 'a' for STalign (default: 500).",
    )
    stalign_use_expression: bool = Field(
        True,
        description="Use gene expression for STalign intensity (vs uniform).",
    )

    # Common parameters
    use_gpu: bool = Field(
        False,
        description="Use GPU acceleration (PASTE with PyTorch backend, STalign).",
    )


class ConditionComparisonParameters(BaseModel):
    """Parameters for multi-sample condition comparison analysis.

    This tool compares gene expression between experimental conditions (e.g., Treatment vs Control)
    across multiple biological samples, using proper statistical methods that account for
    sample-level variation.

    Key difference from find_markers:
    - find_markers: Compares cell types/clusters WITHIN a dataset (e.g., T cell vs B cell)
    - compare_conditions: Compares CONDITIONS ACROSS samples (e.g., Treatment vs Control)
    """

    # Required parameters
    condition_key: str = Field(
        ...,
        description="Column with experimental conditions (e.g., 'treatment').",
    )

    condition1: str = Field(
        ...,
        description="First condition (typically experimental/treatment group).",
    )

    condition2: str = Field(
        ...,
        description="Second condition (typically control/reference group).",
    )

    sample_key: str = Field(
        ...,
        description="Column with sample/replicate IDs. Critical for statistical inference.",
    )

    # Optional parameters
    cell_type_key: Optional[str] = Field(
        None,
        description="Cell type column. If provided, analysis is cell type-stratified.",
    )

    method: Literal["pseudobulk"] = Field(
        "pseudobulk",
        description="Pseudobulk aggregation with DESeq2. Requires 2+ samples per condition.",
    )

    n_top_genes: Annotated[int, Field(gt=0, le=500)] = Field(
        50,
        description="Top DEGs to return per comparison.",
    )

    min_cells_per_sample: Annotated[int, Field(gt=0)] = Field(
        10,
        description="Minimum cells per sample for inclusion.",
    )

    min_samples_per_condition: Annotated[int, Field(gt=0)] = Field(
        2,
        description="Minimum samples per condition for DESeq2.",
    )

    padj_threshold: Annotated[float, Field(gt=0, lt=1)] = Field(
        0.05,
        description="Adjusted p-value threshold.",
    )

    log2fc_threshold: Annotated[float, Field(ge=0)] = Field(
        0.0,
        description="Minimum absolute log2 fold change. 0 = no filtering.",
    )
