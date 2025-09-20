"""
Data models for spatial transcriptomics analysis.
"""

from __future__ import annotations

from typing import Annotated, Any, Dict, List, Literal, Optional, Tuple, Union

from pydantic import BaseModel, Field, model_validator
from typing_extensions import Self


class SpatialDataset(BaseModel):
    """Spatial transcriptomics dataset model"""

    id: str
    name: str
    data_type: Literal[
        "10x_visium", "slide_seq", "merfish", "seqfish", "other", "h5ad", "auto"
    ]
    description: Optional[str] = None


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
    normalization: Literal["log", "sct", "pearson_residuals", "none", "scvi"] = "log"
    scale: bool = True
    n_hvgs: Annotated[int, Field(gt=0, le=5000)] = 2000
    n_pcs: Annotated[int, Field(gt=0, le=100)] = 30

    # scVI preprocessing parameters
    use_scvi_preprocessing: bool = False  # Whether to use scVI for preprocessing
    scvi_n_hidden: int = 128
    scvi_n_latent: int = 10
    scvi_n_layers: int = 1
    scvi_dropout_rate: float = 0.1
    scvi_gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb"

    # Key naming parameters (configurable hard-coded keys)
    cluster_key: str = Field(
        "leiden", alias="clustering_key"
    )  # Key name for storing clustering results
    spatial_key: str = "spatial"  # Key name for spatial coordinates in obsm
    batch_key: str = "batch"  # Key name for batch information in obs

    # User-controllable adaptive parameters (None = use automatic detection)
    n_neighbors: Optional[Annotated[int, Field(gt=2, le=50)]] = (
        None  # Number of neighbors for graph construction (None = adaptive: 3-10 based on dataset size)
    )
    clustering_resolution: Optional[Annotated[float, Field(gt=0.1, le=2.0)]] = (
        None  # Leiden clustering resolution (None = adaptive: 0.4-0.8 based on dataset size)
    )

    # Advanced preprocessing options
    enable_rna_velocity: bool = False  # Whether to include RNA velocity preprocessing
    velocity_params: Optional["RNAVelocityParameters"] = (
        None  # Embedded velocity parameters (if None, uses defaults)
    )

    # Deprecated: Individual velocity fields (kept for backward compatibility)
    velocity_mode: Optional[Literal["stochastic", "deterministic", "dynamical"]] = (
        None  # Deprecated: use velocity_params.mode
    )
    velocity_min_shared_counts: Optional[int] = (
        None  # Deprecated: use velocity_params.min_shared_counts
    )
    velocity_n_top_genes: Optional[int] = (
        None  # Deprecated: use velocity_params.n_top_genes
    )
    velocity_n_pcs: Optional[int] = None  # Deprecated: use velocity_params.n_pcs
    velocity_n_neighbors: Optional[int] = (
        None  # Deprecated: use velocity_params.n_neighbors
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
        "gaston_isodepth",
        "gaston_domains",
        "gaston_genes",
        "pathway_enrichment",
        "spatial_enrichment",  # Clear enrichment types
        "spatial_interaction",
        "integration_check",  # NEW: Enhanced visualization types
    ] = "spatial"
    colormap: str = "viridis"

    # Spatial analysis visualization parameters
    analysis_sub_type: Optional[
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
        description="Sub-type for spatial_analysis plot type. Determines which spatial analysis result to visualize.",
    )
    cluster_key: Optional[str] = Field(
        None,
        description="Cluster key for spatial analysis visualization (e.g., 'leiden')",
    )  # For spatial analysis visualization

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
    dpi: int = 100
    alpha: float = 0.8
    spot_size: Optional[float] = None  # Auto-determined if None

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
    outline_width: float = Field(1.0, description="Line width for cluster outlines")
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

    # Legacy parameters (for backward compatibility)
    show_deconvolution: bool = False  # Whether to show deconvolution results
    n_cell_types: Annotated[int, Field(gt=0, le=10)] = (
        4  # Number of top cell types to show
    )

    @model_validator(mode="after")
    def validate_conditional_parameters(self) -> Self:
        """Validate parameter dependencies and provide helpful error messages."""

        # Spatial analysis validation
        if self.plot_type == "spatial_analysis":
            if not self.analysis_sub_type or (
                isinstance(self.analysis_sub_type, str)
                and not self.analysis_sub_type.strip()
            ):
                available_subtypes = [
                    "neighborhood",
                    "co_occurrence",
                    "ripley",
                    "moran",
                    "centrality",
                    "getis_ord",
                ]
                raise ValueError(
                    f"Parameter dependency error: analysis_sub_type is required when plot_type='spatial_analysis'.\n"
                    f"Available analysis types: {', '.join(available_subtypes)}\n"
                    f"Example usage: VisualizationParameters(plot_type='spatial_analysis', analysis_sub_type='neighborhood')\n"
                    f"For more details, see spatial analysis documentation."
                )

        # Future: Add other conditional validations here
        # if self.plot_type == "cell_communication" and not self.method:
        #     raise ValueError("method required for cell_communication plot_type")

        return self


class AnnotationParameters(BaseModel):
    """Cell type annotation parameters model"""

    method: Literal[
        "marker_genes",
        "tangram",
        "scanvi",
        "cellassign",
        "mllmcelltype",
        "sctype",
        "singler",
    ] = "marker_genes"
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
        None  # For Tangram method - cluster label in reference data
    )
    cell_type_key: Optional[str] = (
        None  # IMPORTANT: Column name for cell types in REFERENCE data. Common values: 'cell_type', 'cell_types', 'celltype'. Leave as None to auto-detect. Only set if you know the exact column name!
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

    # scANVI parameters
    scanvi_n_hidden: int = 128
    scanvi_n_latent: int = 10
    scanvi_n_layers: int = 1
    scanvi_dropout_rate: float = 0.1
    scanvi_unlabeled_category: str = "Unknown"

    # SCVI pretraining parameters (official best practice - enabled by default)
    scanvi_use_scvi_pretrain: bool = True  # Default: True (official best practice)
    scanvi_scvi_epochs: int = (
        200  # SCVI pretraining epochs (reduced for faster training)
    )
    scanvi_n_samples_per_label: int = 100  # For semi-supervised training

    # Query training parameters
    scanvi_query_epochs: int = 100  # Epochs for query data (official recommendation)
    scanvi_check_val_every_n_epoch: int = 10  # Validation frequency

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
        "openai", "anthropic", "deepseek", "kimi", "glm", "qianfan", "ollama"
    ] = "openai"  # LLM provider
    mllm_model: Optional[str] = (
        None  # Model name (e.g., "gpt-5", "claude-sonnet-4-20250514")
    )
    mllm_api_key: Optional[str] = None  # API key for the LLM provider

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
        "geary",
        "centrality",
        "getis_ord",
        "bivariate_moran",
        "join_count",
        "network_properties",
        "spatial_centrality",
        # "scviva"  # TODO: Under development - produces NaN values, needs fixing
    ] = "neighborhood"
    cluster_key: str = "leiden"
    n_neighbors: Annotated[int, Field(gt=0)] = 15

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
    moran_genes: Optional[List[str]] = Field(
        None, description="Specific genes for Moran's I (None = use HVG)"
    )
    moran_n_genes: Annotated[int, Field(gt=0, le=100)] = Field(
        10, description="Number of HVG for Moran's I (default 10 for speed)"
    )
    moran_n_perms: Annotated[int, Field(gt=0, le=10000)] = Field(
        10,
        description="Number of permutations (default 10 for speed, use 100+ for publication)",
    )
    moran_two_tailed: bool = Field(False, description="Use two-tailed test")

    # Getis-Ord Gi* specific parameters
    getis_ord_genes: Optional[List[str]] = (
        None  # Specific genes to analyze (if None, use highly variable genes)
    )
    getis_ord_n_genes: Annotated[int, Field(gt=0, le=100)] = (
        20  # Number of top highly variable genes to analyze
    )
    getis_ord_correction: Literal["bonferroni", "fdr_bh", "none"] = (
        "fdr_bh"  # Multiple testing correction
    )
    getis_ord_alpha: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.05  # Significance threshold
    )

    # Bivariate Moran's I specific parameters
    gene_pairs: Optional[List[Tuple[str, str]]] = Field(
        None, description="Gene pairs for bivariate analysis"
    )

    # SCVIVA deep learning parameters
    scviva_n_epochs: int = Field(
        1000, description="Number of training epochs for SCVIVA"
    )
    scviva_n_hidden: int = Field(128, description="Number of hidden units in SCVIVA")
    scviva_n_latent: int = Field(
        10, description="Number of latent dimensions in SCVIVA"
    )
    scviva_use_gpu: bool = Field(False, description="Use GPU for SCVIVA training")


class RNAVelocityParameters(BaseModel):
    """RNA velocity analysis parameters model"""

    model_config = {
        "extra": "forbid"
    }  # Strict validation - no extra parameters allowed

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


class TrajectoryParameters(BaseModel):
    """Trajectory analysis parameters model"""

    method: Literal["cellrank", "palantir", "velovi", "dpt"] = "cellrank"
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

    # VeloVI parameters
    velovi_n_hidden: int = 128
    velovi_n_latent: int = 10
    velovi_n_layers: int = 1
    velovi_dropout_rate: float = 0.1
    velovi_learning_rate: float = 1e-3

    # Fallback control
    # Removed: allow_fallback_to_dpt - No longer doing automatic fallbacks
    # LLMs should explicitly choose which method to use


class IntegrationParameters(BaseModel):
    """Sample integration parameters model"""

    method: Literal[
        "harmony", "bbknn", "scanorama", "mnn", "scvi", "multivi", "totalvi"
    ] = "harmony"
    batch_key: str = "batch"  # Batch information key
    n_pcs: Annotated[int, Field(gt=0, le=100)] = (
        30  # Number of principal components for integration
    )
    align_spatial: bool = True  # Whether to align spatial coordinates
    reference_batch: Optional[str] = None  # Reference batch for spatial alignment

    # scVI integration parameters
    scvi_n_hidden: int = 128
    scvi_n_latent: int = 10
    scvi_n_layers: int = 1
    scvi_dropout_rate: float = 0.1
    scvi_gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb"

    # MultiVI parameters
    multivi_n_hidden: int = 128
    multivi_n_latent: int = 10
    multivi_n_layers: int = 1
    multivi_dropout_rate: float = 0.1

    # TotalVI parameters
    totalvi_n_hidden: int = 128
    totalvi_n_latent: int = 10
    totalvi_n_layers: int = 1
    totalvi_dropout_rate: float = 0.1
    totalvi_protein_dispersion: Literal[
        "gene", "gene-batch", "gene-label", "gene-cell"
    ] = "gene"


class DeconvolutionParameters(BaseModel):
    """Spatial deconvolution parameters model"""

    method: Literal[
        "cell2location", "rctd", "destvi", "stereoscope", "spotlight", "tangram"
    ] = "cell2location"
    reference_data_id: Optional[str] = (
        None  # Reference single-cell data for deconvolution
    )
    cell_type_key: str = "cell_type"  # Key in reference data for cell type information
    n_top_genes: Annotated[int, Field(gt=0, le=5000)] = (
        2000  # Number of top genes to use
    )
    use_gpu: bool = False  # Whether to use GPU for cell2location
    n_epochs: Annotated[int, Field(gt=0)] = 10000  # Number of epochs for cell2location
    n_cells_per_spot: Optional[int] = None  # Expected number of cells per spot
    reference_profiles: Optional[Dict[str, List[float]]] = (
        None  # Reference expression profiles
    )

    # DestVI parameters
    destvi_n_hidden: int = 128
    destvi_n_latent: int = 10
    destvi_n_layers: int = 1
    destvi_dropout_rate: float = 0.1
    destvi_learning_rate: float = 1e-3

    # Stereoscope parameters
    stereoscope_n_epochs: int = 10000
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


class SpatialDomainParameters(BaseModel):
    """Spatial domain identification parameters model"""

    method: Literal["spagcn", "leiden", "louvain", "stagate", "banksy"] = "spagcn"
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

    # BANKSY specific parameters
    banksy_n_neighbors: Optional[int] = (
        None  # Number of spatial neighbors (default: 15)
    )
    banksy_lambda: Optional[float] = (
        None  # Lambda parameter for spatial weight (default: 0.2)
    )
    banksy_max_m: Optional[int] = None  # Maximum order of neighbors (default: 1)
    banksy_decay_type: Optional[
        Literal["uniform", "reciprocal", "gaussian", "scaled_gaussian"]
    ] = None  # Decay type (default: "scaled_gaussian")
    banksy_n_pcs: Optional[int] = None  # Number of principal components (default: 20)
    
    # Simple timeout configuration
    timeout: Optional[int] = None  # Timeout in seconds (default: 600)


class SpatialVariableGenesParameters(BaseModel):
    """Spatial variable genes identification parameters model"""

    # Method selection
    method: Literal["gaston", "spatialde", "sparkx"] = (
        "sparkx"  # Default to SPARK-X (fastest, no installation needed)
    )

    # Common parameters for all methods
    n_top_genes: Optional[Annotated[int, Field(gt=0, le=5000)]] = (
        None  # Number of top spatial variable genes to return (None = all significant)
    )
    spatial_key: str = "spatial"  # Key in obsm containing spatial coordinates

    # GASTON-specific parameters
    # Preprocessing parameters
    preprocessing_method: Literal["glmpca", "pearson_residuals"] = (
        "glmpca"  # Official default: GLM-PCA
    )
    n_components: Annotated[int, Field(gt=0, le=50)] = (
        10  # Number of components for dimensionality reduction
    )

    # Neural network architecture parameters
    spatial_hidden_layers: List[Annotated[int, Field(gt=0, le=1000)]] = [
        20,
        20,
    ]  # Architecture for spatial embedding f_S
    expression_hidden_layers: List[Annotated[int, Field(gt=0, le=1000)]] = [
        20,
        20,
    ]  # Architecture for expression function f_A

    # Training parameters
    epochs: Annotated[int, Field(gt=0, le=50000)] = 10000  # Number of training epochs
    learning_rate: Annotated[float, Field(gt=0.0, le=1.0)] = 0.001  # Learning rate
    optimizer: Literal["adam", "sgd", "adagrad"] = "adam"  # Optimizer type
    batch_size: Optional[Annotated[int, Field(gt=0, le=10000)]] = (
        None  # Batch size (None for full batch)
    )

    # Positional encoding parameters
    use_positional_encoding: bool = False  # Whether to use positional encoding
    embedding_size: Annotated[int, Field(gt=0, le=20)] = (
        4  # Positional encoding embedding size
    )
    sigma: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.2  # Positional encoding sigma parameter
    )

    # Model saving parameters
    checkpoint_interval: Annotated[int, Field(gt=0, le=5000)] = (
        500  # Save model every N epochs (official tutorial: 500)
    )
    random_seed: Annotated[int, Field(ge=0, le=1000)] = (
        0  # Random seed for reproducibility
    )

    # Analysis parameters
    n_domains: Annotated[int, Field(gt=0, le=20)] = (
        5  # Number of spatial domains to identify
    )
    num_bins: Annotated[int, Field(gt=10, le=200)] = (
        70  # Number of bins for isodepth binning
    )
    umi_threshold: Annotated[int, Field(gt=0)] = (
        500  # Minimum UMI count threshold for genes
    )

    # Gene classification parameters
    continuous_quantile: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.9  # Quantile threshold for continuous genes
    )
    discontinuous_quantile: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.9  # Quantile threshold for discontinuous genes
    )
    pvalue_threshold: Annotated[float, Field(gt=0.0, le=1.0)] = (
        0.1  # P-value threshold for slope significance
    )
    zero_fit_threshold: Annotated[int, Field(ge=0)] = (
        0  # Minimum non-zero count per domain
    )

    # Poisson regression parameters
    isodepth_mult_factor: Annotated[float, Field(gt=0.0)] = (
        1.0  # Scaling factor for isodepth values
    )
    regularization: Annotated[float, Field(ge=0.0)] = 0.0  # Regularization parameter

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

    # Gene count thresholds
    min_genes_required: int = 1000  # Minimum genes required for LIANA analysis
    force_analysis_with_few_genes: bool = (
        False  # Force analysis with insufficient genes
    )

    # Preprocessing validation
    validate_preprocessing: bool = True  # Whether to validate data preprocessing
    preprocessing_strictness: Literal["strict", "warning", "disabled"] = "warning"

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
    cell_type_column: str = "cell_type"  # Which column to use for cell types
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

    # Logging control
    verbose_validation: bool = True  # Output detailed validation info
    log_parameter_choices: bool = True  # Log parameter choice reasons

    # Expert mode
    expert_mode: bool = False  # Expert mode skips some safety checks

    # ========== Backward Compatibility Parameters ==========
    min_cells: Annotated[int, Field(ge=0)] = (
        3  # Minimum cells expressing ligand or receptor
    )

    # ========== Result Control ==========
    plot_top_pairs: Annotated[int, Field(gt=0, le=20)] = (
        6  # Number of top LR pairs to include in results
    )

    # Result filtering
    min_lr_score_threshold: Optional[float] = None  # Minimum L-R score threshold
    filter_low_confidence_pairs: bool = True

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

    # ========== CellChat Specific Parameters ==========
    cellchat_type: Literal["triMean", "truncatedMean", "median"] = (
        "triMean"  # Method for computing communication probability
    )
    cellchat_trim: Annotated[float, Field(ge=0.0, le=0.5)] = (
        0.1  # Trimming parameter for truncatedMean
    )
    cellchat_population_size: bool = False  # Whether to consider population size effect

    # ========== Advanced User Options ==========
    custom_lr_pairs: Optional[List[Tuple[str, str]]] = (
        None  # Custom LR pairs as (ligand, receptor) tuples
    )


class EnrichmentParameters(BaseModel):
    """Parameters for gene set enrichment analysis"""

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
    gene_set_database: Optional[str] = (
        "GO_Biological_Process"  # Gene set database for enrichr
    )

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
