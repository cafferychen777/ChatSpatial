"""
Data models for spatial transcriptomics analysis.
"""

from __future__ import annotations

from typing import Annotated, Dict, List, Literal, Optional, Tuple, Union

from pydantic import BaseModel, ConfigDict, Field, model_validator
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


class PreprocessingParameters(BaseModel):
    """Preprocessing parameters model"""

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

    # ========== Mitochondrial and Ribosomal Gene Filtering ==========
    filter_mito_pct: Optional[float] = Field(
        default=20.0,
        ge=0.0,
        le=100.0,
        description=(
            "Filter spots/cells with mitochondrial percentage above this threshold.\n\n"
            "DEFAULT: 20.0 (remove spots with >20% mitochondrial reads)\n\n"
            "RATIONALE:\n"
            "High mitochondrial content often indicates cell stress, damage, or apoptosis.\n"
            "Damaged cells release cytoplasmic mRNA while retaining mitochondrial transcripts.\n\n"
            "RECOMMENDED VALUES:\n"
            "• 20.0 (default): Standard threshold for most tissues\n"
            "• 5-10: Stringent filtering for high-quality data\n"
            "• 30-50: Relaxed for tissues with naturally high mito (muscle, neurons)\n"
            "• None: Disable filtering (not recommended)\n\n"
            "TISSUE-SPECIFIC CONSIDERATIONS:\n"
            "• Brain: White matter naturally has higher mito% than gray matter\n"
            "• Muscle/Heart: High mito% is biologically normal\n"
            "• Tumor samples: May have elevated mito% due to metabolic changes\n\n"
            "REFERENCE:\n"
            "OSTA Book: lmweber.org/OSTA/pages/seq-quality-control.html"
        ),
    )
    remove_mito_genes: bool = Field(
        default=True,
        description=(
            "Remove mitochondrial genes (MT-*, mt-*) before HVG selection.\n\n"
            "DEFAULT: True (recommended for most analyses)\n\n"
            "RATIONALE:\n"
            "Mitochondrial genes can dominate HVG selection due to high expression\n"
            "and technical variation, masking biologically relevant genes.\n\n"
            "WHEN TO ENABLE (True):\n"
            "• Standard spatial transcriptomics analysis\n"
            "• Clustering and cell type identification\n"
            "• Trajectory analysis\n\n"
            "WHEN TO DISABLE (False):\n"
            "• Studying mitochondrial biology or metabolism\n"
            "• Analyzing mitochondrial heteroplasmy\n"
            "• When mito genes are biologically relevant to your question\n\n"
            "NOTE: Genes are only excluded from HVG selection, not removed from data.\n"
            "They remain available in adata.raw for downstream analyses."
        ),
    )
    remove_ribo_genes: bool = Field(
        default=False,
        description=(
            "Remove ribosomal genes (RPS*, RPL*, Rps*, Rpl*) before HVG selection.\n\n"
            "DEFAULT: False (ribosomal genes often carry biological signal)\n\n"
            "RATIONALE:\n"
            "Ribosomal genes are highly expressed housekeeping genes. While they\n"
            "add noise in some analyses, they can be informative for cell state.\n\n"
            "WHEN TO ENABLE (True):\n"
            "• When ribosomal genes dominate your HVG list\n"
            "• For cleaner clustering focused on cell type markers\n"
            "• Following certain published pipelines that recommend it\n\n"
            "WHEN TO KEEP DISABLED (False):\n"
            "• Standard analyses (ribosomal content varies by cell type)\n"
            "• Studying translation or ribosome biogenesis\n"
            "• When unsure - ribosomal genes rarely cause problems"
        ),
    )

    # Normalization and scaling parameters
    normalization: Literal["log", "sct", "pearson_residuals", "none", "scvi"] = Field(
        default="log",
        description=(
            "Normalization method for gene expression data.\n\n"
            "AVAILABLE OPTIONS:\n"
            "• 'log' (default): Standard log(x+1) normalization after library size correction. "
            "Robust and widely used for most analyses.\n"
            "• 'sct': SCTransform v2 variance-stabilizing normalization via R's sctransform package. "
            "Best for raw UMI counts from 10x platforms. Based on regularized negative binomial regression.\n"
            "• 'pearson_residuals': Analytic Pearson residuals (scanpy built-in, similar to SCTransform). "
            "Requires raw integer counts and scanpy>=1.9.0. Faster than SCTransform with similar results.\n"
            "• 'none': Skip normalization. Use when data is already pre-normalized.\n"
            "• 'scvi': Deep learning-based normalization using scVI variational autoencoder. "
            "Learns a latent representation (X_scvi) that replaces PCA. Best for batch correction and denoising.\n\n"
            "REQUIREMENTS:\n"
            "• sct: R with sctransform package (R -e 'install.packages(\"sctransform\")') + rpy2\n"
            "• pearson_residuals: Raw count data (integers only), scanpy>=1.9.0\n"
            "• scvi: scvi-tools package (pip install scvi-tools), raw count data\n"
            "• none: Data should already be normalized (will warn if raw counts detected)\n\n"
            "RECOMMENDATIONS:\n"
            "• For raw Visium/Xenium/MERFISH data: 'sct', 'pearson_residuals', or 'log'\n"
            "• For Seurat workflow compatibility: 'sct' (SCTransform v2)\n"
            "• For speed with similar results: 'pearson_residuals'\n"
            "• For pre-processed data: 'none'\n"
            "• For batch effect correction and denoising: 'scvi' (deep learning-based)"
        ),
    )
    scale: bool = Field(
        default=False,
        description=(
            "Scale gene expression to unit variance before PCA.\n\n"
            "DEFAULT: False (following Scanpy spatial transcriptomics best practices)\n\n"
            "RATIONALE:\n"
            "The standard Scanpy spatial transcriptomics tutorials do NOT include scaling:\n"
            "  normalize_total → log1p → HVG selection → PCA\n"
            "Scaling is omitted because log-normalization already stabilizes variance.\n\n"
            "WHEN TO ENABLE (scale=True):\n"
            "• Using methods that explicitly require scaled input (e.g., GraphST)\n"
            "• When gene expression magnitudes vary dramatically\n"
            "• For compatibility with Seurat's ScaleData() workflow\n\n"
            "WHEN TO KEEP DISABLED (scale=False):\n"
            "• Standard Visium/spatial analysis with Scanpy/Squidpy\n"
            "• Using SCTransform normalization (already variance-stabilized)\n"
            "• Using Pearson residuals normalization\n\n"
            "REFERENCE:\n"
            "Scanpy spatial tutorial: scanpy-tutorials.readthedocs.io/en/latest/spatial/"
        ),
    )
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

    # SCTransform preprocessing parameters (requires R + sctransform package via rpy2)
    # Installation: R -e 'install.packages("sctransform")' && pip install rpy2
    sct_var_features_n: int = Field(
        default=3000,
        ge=100,
        le=10000,
        description="Number of highly variable features for SCTransform (default: 3000)",
    )
    sct_method: Literal["offset", "fix-slope"] = Field(
        default="fix-slope",
        description=(
            "SCTransform regularization method:\n"
            "• 'fix-slope' (default, v2): Fixed slope regularization, more robust and recommended.\n"
            "• 'offset': Original offset model from v1."
        ),
    )
    sct_exclude_poisson: bool = Field(
        default=True,
        description="Exclude Poisson genes from regularization (v2 default: True). "
        "Improves robustness by excluding genes where variance ≤ mean.",
    )
    sct_n_cells: Optional[int] = Field(
        default=5000,
        ge=100,
        description="Number of cells to subsample for parameter estimation (default: 5000). "
        "Set to None to use all cells (slower but may be more accurate for small datasets).",
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
        description=(
            "Maximum number of training epochs for scVI. "
            "Default 400 is sufficient for most datasets with early stopping enabled. "
            "Increase to 600-800 for large/complex datasets without early stopping."
        ),
    )
    scvi_early_stopping: bool = Field(
        default=True,
        description=(
            "Whether to enable early stopping based on validation ELBO. "
            "STRONGLY RECOMMENDED: Prevents overfitting and reduces training time. "
            "Set to False only for debugging or when you need exact epoch control."
        ),
    )
    scvi_early_stopping_patience: Annotated[int, Field(gt=0, le=100)] = Field(
        default=20,
        description=(
            "Number of epochs to wait for validation improvement before stopping. "
            "Default 20 balances convergence detection with training stability. "
            "Increase to 30-50 for noisy data, decrease to 10-15 for faster training."
        ),
    )
    scvi_train_size: Annotated[float, Field(gt=0.5, le=1.0)] = Field(
        default=0.9,
        description=(
            "Fraction of data used for training (rest for validation). "
            "Default 0.9 (90% train, 10% validation) is standard practice. "
            "Use 1.0 to disable validation (NOT RECOMMENDED - no early stopping)."
        ),
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


class VisualizationParameters(BaseModel):
    """Visualization parameters model"""

    model_config = ConfigDict(extra="forbid")  # Strict validation after preprocessing

    feature: Optional[Union[str, List[str]]] = Field(
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
                return {"feature": feature, "plot_type": "spatial"}
            else:
                return {"feature": data, "plot_type": "spatial"}

        # Handle dict format - normalize features/feature naming
        if isinstance(data, dict):
            data_copy = data.copy()
            # Handle 'features' as alias for 'feature'
            if "features" in data_copy and "feature" not in data_copy:
                data_copy["feature"] = data_copy.pop("features")
            return data_copy

        # For other types (e.g., VisualizationParameters instances), return as-is
        return data

    plot_type: Literal[
        "spatial",
        "heatmap",
        "violin",
        "umap",
        "dotplot",  # Marker gene expression dotplot
        "cell_communication",
        "deconvolution",
        "trajectory",
        "rna_velocity",
        "spatial_statistics",
        "multi_gene",
        "lr_pairs",
        "gene_correlation",
        "pathway_enrichment",
        "spatial_interaction",
        "batch_integration",  # Batch integration quality assessment
        "cnv_heatmap",  # CNV analysis heatmap
        "spatial_cnv",  # CNV spatial projection
        "card_imputation",  # CARD imputation high-resolution results
    ] = "spatial"
    colormap: str = "coolwarm"

    # Unified subtype parameter for all visualization types with subtypes
    subtype: Optional[str] = Field(
        None,
        description=(
            "Unified subtype parameter for visualization variants. "
            "Usage depends on plot_type:\n"
            "- rna_velocity: 'stream' (default, velocity embedding stream), "
            "'phase' (spliced vs unspliced phase plot), 'proportions' (pie chart of spliced/unspliced ratios), "
            "'heatmap' (gene expression by latent_time), 'paga' (PAGA with velocity arrows)\n"
            "- trajectory: 'pseudotime' (default, pseudotime on embedding), "
            "'circular' (CellRank circular projection), 'fate_map' (aggregated fate probabilities), "
            "'gene_trends' (gene expression along lineages), 'fate_heatmap' (smoothed expression heatmap), "
            "'palantir' (Palantir comprehensive results)\n"
            "- pathway_enrichment: 'barplot', 'dotplot' (traditional ORA/GSEA), "
            "'spatial_score', 'spatial_correlogram', 'spatial_variogram', 'spatial_cross_correlation' (spatial EnrichMap)\n"
            "- deconvolution: 'spatial_multi', 'dominant_type', 'diversity', 'stacked_bar', 'scatterpie', 'umap'\n"
            "- spatial_statistics: 'neighborhood', 'co_occurrence', 'ripley', 'moran', 'centrality', 'getis_ord'\n"
            "- Other plot types may not require this parameter"
        ),
    )
    cluster_key: Optional[str] = Field(
        None,
        description=(
            "Column name in adata.obs containing cluster or cell type labels "
            "(e.g., 'leiden', 'louvain', 'cell_type'). "
            "REQUIRED for plot_type='heatmap' and 'violin'. "
            "NOTE: ChatSpatial uses 'cluster_key' (not 'groupby' as in Scanpy) "
            "for consistency with Squidpy spatial analysis functions."
        ),
    )

    # Multi-gene visualization parameters
    multi_panel: bool = False  # Whether to create multi-panel plots
    panel_layout: Optional[Tuple[int, int]] = (
        None  # (rows, cols) - auto-determined if None
    )

    # GridSpec subplot spacing parameters (for multi-panel plots)
    subplot_wspace: float = Field(
        0.0,
        ge=-0.3,  # Allow larger negative values for extreme tight spacing
        le=1.0,
        description=(
            "Horizontal spacing between subplots (GridSpec wspace parameter). "
            "Fraction of average subplot width. "
            "Default 0.0 provides tight spacing for spatial plots with colorbars. "
            "Common values: 0.0 (tight), 0.05 (compact), 0.1 (normal), 0.2 (loose). "
            "Negative values (-0.1 to -0.2) create overlapping spacing for extreme compactness."
        ),
    )
    subplot_hspace: float = Field(
        0.3,
        ge=0.0,
        le=1.0,
        description=(
            "Vertical spacing between subplots (GridSpec hspace parameter). "
            "Fraction of average subplot height. "
            "Default 0.3 provides comfortable spacing. "
            "Common values: 0.2 (tight), 0.3 (normal), 0.4 (loose)."
        ),
    )

    # Colorbar parameters (for spatial plots with make_axes_locatable)
    colorbar_pad: float = Field(
        0.02,
        ge=0.0,
        le=0.2,
        description=(
            "Distance between subplot and colorbar (as fraction of subplot width). "
            "Default 0.02 provides tight spacing. "
            "Common values: 0.02 (tight), 0.03 (compact), 0.05 (normal)."
        ),
    )
    colorbar_size: str = Field(
        "3%",
        description=(
            "Width of colorbar as percentage of subplot width. "
            "Default '3%' provides narrow colorbar to save space. "
            "Common values: '3%' (narrow), '4%' (compact), '5%' (normal)."
        ),
    )

    # Ligand-receptor pair parameters
    lr_pairs: Optional[List[Tuple[str, str]]] = None  # List of (ligand, receptor) pairs
    lr_database: str = "cellchat"  # Database for LR pairs
    plot_top_pairs: int = Field(
        6,
        gt=0,
        le=100,
        description="Number of top LR pairs to display in cell communication visualization. Default: 6. For chord diagrams, use higher values (e.g., 50) to show more interactions.",
    )

    # Gene correlation parameters
    correlation_method: Literal["pearson", "spearman", "kendall"] = "pearson"
    show_correlation_stats: bool = True

    # Figure parameters
    figure_size: Optional[Tuple[int, int]] = (
        None  # (width, height) - auto-determined if None
    )
    dpi: int = 300  # Publication quality (Nature/Cell standard)
    alpha: float = 0.9  # Spot transparency (higher = more opaque)
    spot_size: Optional[float] = Field(
        150.0,
        description=(
            "Size of spots in spatial plots (in pixels). "
            "Default 150 provides good balance for most 10x Visium data. "
            "Adjust based on data density: "
            "dense (>3000 spots): 100-150, "
            "sparse (<2000 spots): 150-200. "
            "Set to None for scanpy auto-sizing (not recommended - usually too small)."
        ),
    )
    alpha_img: float = Field(
        0.3,
        ge=0.0,
        le=1.0,
        description=(
            "Background tissue image transparency (lower = dimmer, helps spots stand out). "
            "Default 0.3 provides good contrast. "
            "Increase to 0.4-0.5 to emphasize tissue structure."
        ),
    )
    show_tissue_image: bool = Field(
        True,
        description=(
            "Whether to show tissue histology image in spatial plots. "
            "If False, only plot spots on coordinates without background image. "
            "This option only applies when tissue image is available. "
            "When False, spots are plotted on a clean coordinate system for clearer visualization. "
            "Default: True"
        ),
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
        description=(
            "Standardize expression values for dotplot. "
            "'var' = standardize per gene (row), "
            "'group' = standardize per group (column)"
        ),
    )
    dotplot_dot_max: Optional[float] = Field(
        None,
        ge=0.0,
        le=1.0,
        description=(
            "Maximum dot size as fraction (0-1). "
            "If None, maximum observed fraction is used"
        ),
    )
    dotplot_dot_min: Optional[float] = Field(
        None,
        ge=0.0,
        le=1.0,
        description=(
            "Minimum dot size as fraction (0-1). "
            "If None, minimum observed fraction is used"
        ),
    )
    dotplot_smallest_dot: float = Field(
        0.0,
        ge=0.0,
        le=50.0,
        description=(
            "Size of dot when expression fraction is 0. "
            "Default 0 hides genes with no expression in a group"
        ),
    )
    dotplot_var_groups: Optional[Dict[str, List[str]]] = Field(
        None,
        description=(
            "Group genes by category for organized display. "
            "Example: {'T cell markers': ['CD3D', 'CD4'], 'B cell markers': ['CD19', 'MS4A1']}"
        ),
    )
    dotplot_categories_order: Optional[List[str]] = Field(
        None,
        description="Custom order for groups (clusters/cell types) on the axis",
    )

    # Deconvolution visualization parameters
    n_cell_types: Annotated[
        int,
        Field(
            gt=0,
            le=10,
            description="Number of top cell types to show in deconvolution visualization. Must be between 1-10. Default: 4",
        ),
    ] = 4
    deconv_method: Optional[str] = Field(
        None,
        description=(
            "Deconvolution method name (e.g., 'cell2location', 'rctd'). "
            "If None and only one result exists, auto-select and notify. "
            "If None and multiple results exist, raise error requiring explicit specification. "
            "This ensures you visualize the intended analysis for scientific reproducibility."
        ),
    )
    min_proportion_threshold: float = Field(
        0.3,
        ge=0.0,
        le=1.0,
        description="Minimum proportion threshold for marking spots as 'pure' vs 'mixed' (dominant_type visualization). Default: 0.3",
    )
    show_mixed_spots: bool = Field(
        True,
        description="Whether to mark mixed/heterogeneous spots in dominant_type visualization. Default: True",
    )
    pie_scale: float = Field(
        0.4,
        gt=0.0,
        le=2.0,
        description="Size scale factor for pie charts in scatterpie visualization. Default: 0.4",
    )
    scatterpie_alpha: float = Field(
        1.0,
        ge=0.0,
        le=1.0,
        description="Transparency of pie charts in scatterpie visualization (0=transparent, 1=opaque). Default: 1.0",
    )
    max_spots: int = Field(
        100,
        gt=0,
        le=1000,
        description="Maximum number of spots to show in stacked_bar visualization. Default: 100",
    )
    sort_by: Literal["dominant_type", "spatial", "cluster"] = Field(
        "dominant_type",
        description="Sorting method for stacked_bar visualization. Options: dominant_type (group by dominant cell type), spatial (spatial order), cluster (cluster order). Default: dominant_type",
    )

    @model_validator(mode="after")
    def validate_conditional_parameters(self) -> Self:
        """Validate parameter dependencies and provide helpful error messages."""

        # Spatial statistics validation
        if self.plot_type == "spatial_statistics":
            if not self.subtype or (
                isinstance(self.subtype, str) and not self.subtype.strip()
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
                    f"Parameter dependency error: subtype is required when plot_type='spatial_statistics'.\n"
                    f"Available subtypes: {', '.join(available_types)}\n"
                    f"Example usage: VisualizationParameters(plot_type='spatial_statistics', subtype='neighborhood')\n"
                    f"For more details, see spatial statistics documentation."
                )

        # Deconvolution validation - set default subtype if not provided
        if self.plot_type == "deconvolution":
            if not self.subtype:
                self.subtype = (
                    "spatial_multi"  # Default deconvolution visualization type
                )

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
    tangram_mode: Literal["cells", "clusters"] = (
        "cells"  # Tangram mapping mode: 'cells' (cell-level) or 'clusters' (cluster-level)
    )
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

    # Tangram-specific parameters (aligned with scvi.external.Tangram API)
    tangram_density_prior: Literal["rna_count_based", "uniform"] = (
        "rna_count_based"  # Density prior for mapping
    )
    tangram_device: str = "cpu"  # Device for computation ('cpu' or 'cuda:0')
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
            "WARNING:Empirical (not official): Small datasets may need 3-5 to avoid NaN"
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
            "WARNING:Empirical (not official): 0.2-0.3 may help small datasets"
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
            "WARNING:For small datasets: Set to False if encountering NaN errors"
        ),
    )
    scanvi_scvi_epochs: int = Field(
        default=200, description="Number of epochs for SCVI pretraining. Default: 200"
    )
    scanvi_scanvi_epochs: int = Field(
        default=20,
        description=(
            "Number of epochs for SCANVI model training after SCVI pretraining. Default: 20\n"
            "This is the second stage training that fine-tunes the model for label transfer.\n"
            "Official scvi-tools recommendation: 20 epochs is usually sufficient after pretraining.\n"
            "Increase to 50-100 for complex datasets or if label transfer accuracy is low."
        ),
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
            "WARNING:For small datasets: Recommend 50 to prevent overfitting"
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
    singler_reference: Optional[str] = Field(
        default=None,
        description=(
            "Reference dataset name from celldex package (Python naming convention).\n\n"
            "Valid references:\n"
            "  Human: 'hpca' (Human Primary Cell Atlas, recommended), 'blueprint_encode', "
            "'dice', 'monaco_immune', 'novershtern_hematopoietic'\n"
            "  Mouse: 'immgen' (ImmGen, recommended), 'mouse_rnaseq'\n\n"
            "Common mistakes:\n"
            "  'HumanPrimaryCellAtlasData' - WRONG, use 'hpca'\n"
            "  'ImmGenData' - WRONG, use 'immgen'\n\n"
            "If None, uses species-appropriate default ('hpca' for human, 'immgen' for mouse)."
        ),
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

    # Local Moran's I (LISA) specific parameters
    local_moran_permutations: Annotated[int, Field(gt=0, le=9999)] = Field(
        999,
        description=(
            "Number of permutations for pseudo p-value calculation in Local Moran's I. "
            "Higher values increase precision: 99 -> precision 0.01, 999 -> precision 0.001. "
            "Default 999 is standard practice. Use 9999 for publication-quality results."
        ),
    )
    local_moran_alpha: Annotated[float, Field(gt=0.0, lt=1.0)] = Field(
        0.05,
        description=(
            "Significance level (alpha) for Local Moran's I hotspot/coldspot detection. "
            "Used with FDR correction to determine significant spatial clusters. "
            "Common values: 0.05 (standard), 0.01 (conservative), 0.10 (exploratory)."
        ),
    )
    local_moran_fdr_correction: bool = Field(
        True,
        description=(
            "Whether to apply FDR (False Discovery Rate) correction for multiple testing. "
            "STRONGLY RECOMMENDED: Each location is tested separately, creating a multiple "
            "testing problem. FDR correction controls the expected proportion of false positives. "
            "Set to False only for exploratory analysis."
        ),
    )

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

    model_config = ConfigDict(
        extra="forbid"
    )  # Strict validation - no extra parameters allowed

    # Velocity computation method selection
    method: Literal["scvelo", "velovi"] = "scvelo"

    # scVelo specific parameters
    scvelo_mode: Literal["deterministic", "stochastic", "dynamical"] = "stochastic"
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
        "flashdeconv",
        "cell2location",
        "rctd",
        "destvi",
        "stereoscope",
        "spotlight",
        "tangram",
        "card",
    ] = "flashdeconv"
    reference_data_id: Optional[str] = (
        None  # Reference single-cell data for deconvolution
    )
    cell_type_key: str  # REQUIRED: Key in reference data for cell type information. LLM will infer from metadata. Common values: 'cell_type', 'celltype', 'annotation', 'label'

    # Universal GPU parameter
    use_gpu: bool = Field(
        False,
        description=(
            "Whether to use GPU acceleration for training. "
            "Supported by: Cell2location, DestVI, Stereoscope, Tangram. "
            "Not supported by: RCTD, SPOTlight, CARD (R-based methods). "
            "Requires CUDA-compatible GPU and proper PyTorch installation."
        ),
    )

    # Cell2location specific parameters
    cell2location_ref_model_epochs: Annotated[int, Field(gt=0)] = Field(
        250,
        description=(
            "Number of epochs for Cell2location reference model training (NB regression). "
            "This is the first stage training for estimating reference cell type signatures. "
            "Official recommendation: 250. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_n_epochs: Annotated[int, Field(gt=0)] = Field(
        30000,
        description=(
            "Number of epochs for Cell2location spatial mapping model training. "
            "Official recommendation: 30000. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_n_cells_per_spot: Annotated[int, Field(gt=0)] = Field(
        30,
        description=(
            "Expected number of cells per spatial location for Cell2location. "
            "This is tissue-dependent (e.g., 30 for Visium, 5-10 for MERFISH). "
            "Official recommendation: 30 for Visium data. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_detection_alpha: Annotated[float, Field(gt=0)] = Field(
        20.0,
        description=(
            "RNA detection sensitivity parameter for Cell2location. "
            "NEW DEFAULT (2024): 20 for high technical variability, 200 for low variability. "
            "Recommendation: test both values on your data. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )

    # Batch and covariate correction for cell2location
    cell2location_batch_key: Optional[str] = Field(
        None,
        description=(
            "Column name in adata.obs for batch information (e.g., 'sample_id', 'batch'). "
            "Used for batch effect correction in Cell2location. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_categorical_covariate_keys: Optional[List[str]] = Field(
        None,
        description=(
            "List of column names in adata.obs for categorical technical covariates "
            "(e.g., ['platform', 'donor_id']) for Cell2location. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )

    # Gene filtering parameters (Cell2location-specific preprocessing)
    cell2location_apply_gene_filtering: bool = Field(
        True,
        description=(
            "Apply Cell2location's recommended permissive gene filtering before training. "
            "ONLY USED BY CELL2LOCATION. This is NOT the same as HVG selection:\n"
            "• Cell2location uses permissive filtering to keep rare cell type markers\n"
            "• Yields ~10k-16k genes (more than typical 2k HVGs)\n"
            "• Official recommendation: avoid further gene selection for robust results\n"
            "Other methods use different strategies (see spotlight_n_top_genes parameter)."
        ),
    )
    cell2location_gene_filter_cell_count_cutoff: int = Field(
        5,
        description=(
            "Minimum cells expressing a gene for Cell2location filtering (official default: 5). "
            "Low cutoff preserves rare cell type markers. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_gene_filter_cell_percentage_cutoff2: float = Field(
        0.03,
        description=(
            "Minimum percentage of cells expressing for Cell2location (official default: 0.03 = 3%). "
            "Genes detected in ≥3% of cells are always included. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_gene_filter_nonz_mean_cutoff: float = Field(
        1.12,
        description=(
            "Minimum non-zero mean expression for Cell2location (official default: 1.12). "
            "For genes between cutoffs, only keep if avg expression in non-zero cells > 1.12. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )

    # Phase 2: Training enhancement parameters (Cell2location)
    cell2location_ref_model_lr: Annotated[float, Field(gt=0)] = Field(
        0.002,
        description=(
            "Reference model learning rate for Cell2location (official default: 0.002 with ClippedAdam optimizer). "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_lr: Annotated[float, Field(gt=0)] = Field(
        0.005,
        description=(
            "Cell2location model learning rate (official default: 0.005). "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_ref_model_train_size: Annotated[float, Field(gt=0, le=1)] = Field(
        1.0,
        description=(
            "Fraction of reference data for training in Cell2location. "
            "DEFAULT: 1.0 (official tutorial recommendation - use all data). "
            "IMPORTANT: RegressionModel validation is not yet implemented, so train_size=1 is standard practice. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_train_size: Annotated[float, Field(gt=0, le=1)] = Field(
        1.0,
        description=(
            "Fraction of spatial data for training in Cell2location. "
            "DEFAULT: 1.0 (official tutorial: 'we need to estimate cell abundance at all locations'). "
            "Using train_size=1 ensures all spatial locations are included in training. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_enable_qc_plots: bool = Field(
        False,
        description=(
            "Generate QC diagnostic plots for Cell2location (ELBO history, convergence diagnostics). "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_qc_output_dir: Optional[str] = Field(
        None,
        description=(
            "Output directory for Cell2location QC plots (None = plots not saved to disk). "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )

    # Phase 3: Runtime optimization parameters (Cell2location)
    cell2location_early_stopping: bool = Field(
        False,
        description=(
            "Enable early stopping to reduce Cell2location training time. "
            "DEFAULT: False (following official tutorial best practice). "
            "IMPORTANT: RegressionModel does not support validation, so early stopping is not recommended. "
            "Official tutorial uses train_size=1 without early stopping. "
            "Only enable if you have specific convergence monitoring needs. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_early_stopping_patience: Annotated[int, Field(gt=0)] = Field(
        45,
        description=(
            "Epochs to wait before stopping if no improvement for Cell2location (official default: 45). "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_early_stopping_threshold: Annotated[float, Field(gt=0)] = Field(
        0.0,
        description=(
            "Minimum relative change to qualify as improvement for Cell2location (0 = any improvement). "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_use_aggressive_training: bool = Field(
        False,
        description=(
            "Use train_aggressive() method for large-scale datasets in Cell2location. "
            "DEFAULT: False (standard train() method, following official tutorial). "
            "WHEN TO USE: Only for datasets with >50k locations that require mini-batch training due to GPU memory constraints. "
            "Standard Visium datasets (<50k locations) should use train_size=1 with batch_size=None (official best practice). "
            "Aggressive training implements amortised inference for scalability to 100k-1M+ locations. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )
    cell2location_validation_size: Annotated[float, Field(gt=0, lt=1)] = Field(
        0.1,
        description=(
            "Fraction of data for validation set in Cell2location (required if early_stopping=True). "
            "NOTE: Official tutorial uses train_size=1 (no validation split) for standard workflows. "
            "ONLY USED BY CELL2LOCATION METHOD."
        ),
    )

    # SPOTlight specific parameters
    spotlight_n_top_genes: Annotated[int, Field(gt=0, le=5000)] = Field(
        2000,
        description=(
            "Number of top highly variable genes (HVGs) to use for SPOTlight deconvolution. "
            "ONLY USED BY SPOTLIGHT METHOD. Other methods use different gene selection strategies:\n"
            "• Cell2location: Uses permissive gene filtering (apply_gene_filtering parameter)\n"
            "• RCTD/DestVI/Stereoscope/CARD/Tangram: Use all common genes between datasets\n"
            "Default: 2000. Recommended range: 1000-3000 for standard Visium data."
        ),
    )
    spotlight_nmf_model: Literal["ns"] = Field(
        "ns",
        description=(
            "NMF model type for SPOTlight. ONLY USED BY SPOTLIGHT METHOD.\n\n"
            "Currently only 'ns' (non-smooth NMF) is supported. This method produces "
            "sparser, more interpretable deconvolution results.\n\n"
            "NOTE: SPOTlight documentation mentions 'std' (standard NMF) as an option, "
            "but it is currently broken in SPOTlight (internally creates 'stdNMF' algorithm "
            "which doesn't exist in the NMF package registry). We only expose working parameters.\n\n"
            "Reference: Elosua-Bayes et al. (2021) Nucleic Acids Research."
        ),
    )
    spotlight_min_prop: Annotated[float, Field(ge=0, le=1)] = Field(
        0.01,
        description=(
            "Minimum cell type proportion threshold for SPOTlight. "
            "Cell types contributing less than this value are filtered out as noise. "
            "Official default: 0.01 (1%). "
            "Lower values = keep more cell types but more noise. "
            "Higher values = stricter filtering but may lose rare cell types. "
            "ONLY USED BY SPOTLIGHT METHOD."
        ),
    )
    spotlight_scale: bool = Field(
        True,
        description=(
            "Whether to scale/normalize data in SPOTlight. "
            "Affects gene expression scale handling. "
            "Default: True (recommended). "
            "ONLY USED BY SPOTLIGHT METHOD."
        ),
    )
    spotlight_weight_id: str = Field(
        "mean.AUC",
        description=(
            "Column name for marker gene weights in SPOTlight. "
            "Specifies which metric to use for weighting marker genes. "
            "Common values: 'mean.AUC' (default), 'median.AUC'. "
            "ONLY USED BY SPOTLIGHT METHOD."
        ),
    )

    # DestVI parameters
    destvi_n_epochs: Annotated[int, Field(gt=0)] = Field(
        2000,
        description=(
            "Number of epochs for DestVI training. "
            "Official recommendation: 2000 (minimum 1000). "
            "ONLY USED BY DESTVI METHOD."
        ),
    )
    destvi_n_hidden: int = 128
    destvi_n_latent: int = 10
    destvi_n_layers: int = 1
    destvi_dropout_rate: float = 0.1
    destvi_learning_rate: float = 1e-3

    # DestVI advanced parameters (official scvi-tools defaults)
    destvi_train_size: Annotated[float, Field(gt=0.0, le=1.0)] = Field(
        default=0.9,
        description=(
            "Fraction of data to use for training DestVI (rest for validation). "
            "Official scvi-tools default: 0.9. "
            "Lower values (0.8) provide more robust validation but less training data. "
            "ONLY USED BY DESTVI METHOD."
        ),
    )
    destvi_vamp_prior_p: Annotated[int, Field(ge=1)] = Field(
        default=15,
        description=(
            "Number of VampPrior components for DestVI. "
            "Official scvi-tools default: 15. "
            "Higher values may improve modeling of complex cell type distributions. "
            "ONLY USED BY DESTVI METHOD."
        ),
    )
    destvi_l1_reg: Annotated[float, Field(ge=0.0)] = Field(
        default=10.0,
        description=(
            "L1 regularization strength for DestVI to encourage sparsity. "
            "Official scvi-tools default: 10.0. "
            "Higher values encourage sparser cell type assignments per spot. "
            "ONLY USED BY DESTVI METHOD."
        ),
    )

    # Stereoscope parameters
    stereoscope_n_epochs: int = 150000
    stereoscope_learning_rate: float = 0.01
    stereoscope_batch_size: int = 128

    # RCTD specific parameters
    rctd_mode: Literal["full", "doublet", "multi"] = Field(
        "full",
        description=(
            "RCTD deconvolution mode (Cable et al. 2022):\n"
            "• 'doublet': Assigns 1-2 cell types per spot, classifies each as 'singlet' or 'doublet'. "
            "Recommended for HIGH-RESOLUTION spatial data (Slide-seq ~10μm, MERFISH, Visium HD)\n"
            "• 'full' (default): Assigns any number of cell types per spot. "
            "Recommended for LOW-RESOLUTION data (standard Visium 55μm spots, 100μm spacing)\n"
            "• 'multi': Extension of doublet mode using greedy algorithm to add multiple cell types. "
            "Alternative to 'full' with more constraints on cell type mixing"
        ),
    )
    max_cores: Annotated[int, Field(gt=0, le=16)] = 4  # Maximum number of cores to use
    rctd_confidence_threshold: Annotated[float, Field(gt=0)] = (
        10.0  # Confidence threshold for cell type assignment (higher = more stringent)
    )
    rctd_doublet_threshold: Annotated[float, Field(gt=0)] = (
        25.0  # Threshold for doublet detection (used in doublet/multi modes)
    )
    rctd_max_multi_types: Annotated[int, Field(ge=2, le=10)] = Field(
        4,
        description=(
            "Maximum number of cell types per spot in RCTD multi mode. "
            "Recommended: 4-6 for Visium (100μm spots), 2-3 for higher resolution. "
            "Must be less than total number of cell types in reference data."
        ),
    )

    # CARD specific parameters
    card_minCountGene: Annotated[int, Field(gt=0)] = Field(
        100,
        description="Minimum total counts per gene across all spots for CARD quality control filtering",
    )
    card_minCountSpot: Annotated[int, Field(gt=0)] = Field(
        5,
        description="Minimum number of spots where a gene must be expressed for CARD quality control",
    )
    card_sample_key: Optional[str] = Field(
        None,
        description="Optional sample/batch column name in reference data for multi-sample CARD analysis",
    )
    card_imputation: bool = Field(
        False,
        description=(
            "Enable CARD spatial imputation to create enhanced high-resolution spatial maps. "
            "CARD's unique CAR (Conditional AutoRegressive) model allows imputation at unmeasured locations, "
            "constructing refined tissue maps with arbitrarily higher resolution than the original measurement. "
            "Extremely fast: 0.4s for all genes (5816x faster than BayesSpace). "
            "Use for: Enhancing Visium to near-cellular resolution, filling tissue gaps, smoothing artifacts"
        ),
    )
    card_NumGrids: Annotated[int, Field(gt=0)] = Field(
        2000,
        description=(
            "Number of spatial grid points for CARD imputation (default: 2000). "
            "Higher values = finer spatial resolution but increased computation. "
            "Typical values: 2000 (standard), 5000 (high-res), 10000 (ultra high-res). "
            "The imputed map will have ~NumGrids locations covering the tissue area"
        ),
    )
    card_ineibor: Annotated[int, Field(gt=0)] = Field(
        10,
        description=(
            "Number of nearest neighbors for CARD spatial imputation (default: 10). "
            "Controls the spatial smoothness of imputed results. "
            "Higher values = smoother maps, lower values = preserve local variation"
        ),
    )

    # Tangram specific parameters
    tangram_n_epochs: Annotated[int, Field(gt=0)] = Field(
        1000,
        description=(
            "Number of epochs for Tangram spatial mapping. "
            "Official recommendation: 1000. "
            "ONLY USED BY TANGRAM METHOD."
        ),
    )
    tangram_mode: Literal["cells", "clusters", "constrained"] = Field(
        "cells",
        description=(
            "Tangram mapping mode. "
            "'cells': Cell-level mapping (default). "
            "'clusters': Cluster-level mapping (requires cluster_label). "
            "'constrained': Constrained optimization with target_count. "
            "Official recommendation: 'cells' for most applications. "
            "ONLY USED BY TANGRAM METHOD."
        ),
    )
    tangram_learning_rate: Annotated[float, Field(gt=0)] = Field(
        0.1,
        description=(
            "Learning rate for Tangram optimizer. "
            "Official default: 0.1. "
            "Higher values = faster convergence but less stable. "
            "Lower values = more stable but slower. "
            "ONLY USED BY TANGRAM METHOD."
        ),
    )
    tangram_density_prior: Literal["rna_count_based", "uniform"] = Field(
        "rna_count_based",
        description=(
            "Spatial density prior for Tangram. "
            "'rna_count_based': Weight by RNA counts (default, recommended). "
            "'uniform': Equal weight for all spots. "
            "Official recommendation: 'rna_count_based' for better biological interpretation. "
            "ONLY USED BY TANGRAM METHOD."
        ),
    )

    # FlashDeconv specific parameters (DEFAULT METHOD - ultra-fast, atlas-scale)
    flashdeconv_sketch_dim: Annotated[int, Field(gt=0, le=2048)] = Field(
        512,
        description=(
            "Dimension of the sketched space for FlashDeconv. "
            "Higher values preserve more information but increase computation. "
            "Default: 512 (recommended for most datasets). "
            "ONLY USED BY FLASHDECONV METHOD."
        ),
    )
    flashdeconv_lambda_spatial: Annotated[float, Field(gt=0)] = Field(
        5000.0,
        description=(
            "Spatial regularization strength for FlashDeconv. "
            "Higher values encourage smoother spatial patterns. "
            "Recommended values by platform:\n"
            "• Standard Visium (55μm): 1000-10000 (default: 5000)\n"
            "• Visium HD (16μm): 5000-20000\n"
            "• Visium HD (8μm): 10000-50000\n"
            "• Visium HD (2μm): 50000-100000\n"
            "• Stereo-seq/Seq-Scope: 50000-200000\n"
            "Use 'auto' for automatic tuning (may underestimate for real data). "
            "ONLY USED BY FLASHDECONV METHOD."
        ),
    )
    flashdeconv_n_hvg: Annotated[int, Field(gt=0, le=5000)] = Field(
        2000,
        description=(
            "Number of highly variable genes to select for FlashDeconv. "
            "Default: 2000. "
            "ONLY USED BY FLASHDECONV METHOD."
        ),
    )
    flashdeconv_n_markers_per_type: Annotated[int, Field(gt=0, le=500)] = Field(
        50,
        description=(
            "Number of marker genes per cell type for FlashDeconv. "
            "Default: 50. "
            "ONLY USED BY FLASHDECONV METHOD."
        ),
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
    spatialde_pi0: Optional[float] = Field(
        default=None,
        gt=0.0,
        le=1.0,
        description=(
            "Prior probability of null hypothesis for SpatialDE q-value estimation. "
            "This represents the expected proportion of genes WITHOUT spatial patterns. "
            "\n\n"
            "VALUES:\n"
            "- None (default, RECOMMENDED): Uses adaptive pi0 estimation from SpatialDE\n"
            "- 0.9: Assumes 10% of genes have spatial patterns (conservative)\n"
            "- 0.5: Assumes 50% of genes have spatial patterns (moderate)\n"
            "- 0.1: Assumes 90% of genes have spatial patterns (aggressive, may increase false positives)\n"
            "\n"
            "SCIENTIFIC NOTE:\n"
            "The pi0 parameter directly affects the stringency of FDR correction. "
            "Lower pi0 values assume more genes are truly spatial, leading to more "
            "liberal q-value estimates and potentially more false positives. "
            "The default adaptive estimation (None) is recommended for most analyses "
            "as it learns pi0 from the data distribution."
        ),
    )

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
    method: Literal["liana", "cellphonedb", "cellchat_r"] = "liana"
    # Methods:
    # - "liana": LIANA+ framework (Python, supports multiple resources)
    # - "cellphonedb": CellPhoneDB v5 (Python)
    # - "cellchat_r": Native R CellChat (full features with mediator proteins & pathways)

    # ========== Species and Resource Control ==========
    species: Literal["human", "mouse", "zebrafish"]
    # REQUIRED: Must explicitly specify species for ligand-receptor database
    # - "human": For human data (genes like ACTB, GAPDH - all uppercase)
    # - "mouse": For mouse data (genes like Actb, Gapdh - capitalized)
    # - "zebrafish": For zebrafish data

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
    liana_local_metric: Literal["cosine", "pearson", "spearman", "jaccard"] = (
        "cosine"  # Local spatial metric
    )
    liana_global_metric: Literal["morans", "lee"] = "morans"  # Global spatial metric
    liana_n_perms: Annotated[int, Field(gt=0)] = (
        1000  # Number of permutations for LIANA (1000 minimum for publication-quality p-values)
    )
    liana_nz_prop: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.2  # Minimum expression proportion
    )
    liana_bandwidth: Optional[int] = None  # Bandwidth for spatial connectivity
    liana_cutoff: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.1  # Cutoff for spatial connectivity
    )
    liana_significance_alpha: Annotated[float, Field(gt=0.0, lt=1.0)] = Field(
        default=0.05,
        description=(
            "Significance threshold (alpha) for FDR-corrected p-values in LIANA analysis.\n"
            "Default: 0.05 (standard statistical threshold).\n"
            "Use 0.01 for more stringent filtering, 0.10 for exploratory analysis.\n"
            "This controls both cluster-level (magnitude_rank) and spatial (FDR-corrected) significance."
        ),
    )

    # ========== Expression Filtering Parameters ==========
    min_cells: Annotated[int, Field(ge=0)] = (
        3  # Minimum cells expressing ligand or receptor (required by LIANA for statistical validity)
    )

    # ========== Result Control ==========
    plot_top_pairs: Annotated[int, Field(gt=0, le=100)] = (
        6  # Number of top LR pairs to include in results (chord diagrams may use 50+)
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

    cellchat_trim: Annotated[float, Field(ge=0.0, le=0.5)] = 0.1
    # Trim proportion for truncatedMean method (0.1 = 10% truncated mean)

    cellchat_population_size: bool = True
    # Whether to consider cell population size effect in communication probability

    cellchat_min_cells: Annotated[int, Field(ge=1)] = 10
    # Minimum number of cells required in each cell group for filterCommunication

    cellchat_distance_use: bool = True
    # Whether to use spatial distance constraints (for spatial data)

    cellchat_interaction_range: Annotated[float, Field(gt=0.0)] = 250.0
    # Maximum interaction/diffusion range of ligands in microns (for spatial data)

    cellchat_scale_distance: Annotated[float, Field(gt=0.0)] = 0.01
    # Scale factor for distance calculation (adjust based on imaging technology)

    cellchat_contact_knn_k: Annotated[int, Field(ge=1)] = 6
    # Number of nearest neighbors for defining contact-dependent signaling
    # Used for spatial data to determine which cells are in contact range

    cellchat_contact_range: Optional[Annotated[float, Field(gt=0.0)]] = None
    # Alternative to contact_knn_k: explicit distance threshold for contact signaling
    # If None, uses contact_knn_k instead (recommended for most spatial data)

    # CellChat spatial conversion factors (platform-specific)
    cellchat_pixel_ratio: Annotated[float, Field(gt=0.0)] = Field(
        default=0.5,
        description=(
            "Conversion factor from image pixels to micrometers (um).\n"
            "Platform-specific defaults:\n"
            "  - Visium (10x): 0.5 (1 pixel ≈ 0.5 um at full resolution)\n"
            "  - MERFISH: Varies by imaging setup, typically 0.1-1.0\n"
            "  - Slide-seq: ~0.5 (10 um beads)\n"
            "  - CosMx: 0.18 (imaging resolution)\n"
            "Used in CellChat's spatial.factors for coordinate conversion."
        ),
    )

    cellchat_spatial_tol: Annotated[float, Field(gt=0.0)] = Field(
        default=27.5,
        description=(
            "Spatial tolerance (half of spot/cell diameter) in micrometers.\n"
            "Platform-specific defaults:\n"
            "  - Visium (10x): 27.5 um (spot diameter ~55um, half is ~27.5)\n"
            "  - MERFISH: 5-10 um (single cell resolution)\n"
            "  - Slide-seq: 5 um (10 um bead diameter / 2)\n"
            "  - CosMx: 5-10 um (single cell resolution)\n"
            "Used in CellChat's spatial.factors.tol for defining spatial proximity."
        ),
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


class RegistrationParameters(BaseModel):
    """Spatial registration parameters for aligning multiple tissue slices."""

    method: Literal["paste", "stalign"] = Field(
        "paste",
        description=(
            "Registration method. 'paste': Probabilistic Alignment of ST Experiments "
            "(optimal transport-based, recommended). 'stalign': STalign diffeomorphic "
            "mapping (LDDMM-based, for complex deformations)."
        ),
    )
    reference_idx: Optional[int] = Field(
        None,
        ge=0,
        description="Index of reference slice (0-indexed). If None, uses first slice.",
    )

    # PASTE-specific parameters
    paste_alpha: Annotated[float, Field(gt=0, le=1)] = Field(
        0.1,
        description=(
            "Spatial regularization parameter for PASTE (0-1). "
            "Higher values give more weight to spatial coordinates vs expression. "
            "Default: 0.1 (expression-dominated alignment)."
        ),
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
    stalign_image_size: Tuple[int, int] = Field(
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
