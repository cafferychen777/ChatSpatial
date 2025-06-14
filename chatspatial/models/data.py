"""
Data models for spatial transcriptomics analysis.
"""

from typing import Annotated, Dict, List, Literal, Optional, Tuple, Union
from pydantic import BaseModel, Field


class SpatialDataset(BaseModel):
    """Spatial transcriptomics dataset model"""
    id: str
    name: str
    data_type: Literal["10x_visium", "slide_seq", "merfish", "seqfish", "other", "h5ad", "auto"]
    description: Optional[str] = None


class AnalysisParameters(BaseModel):
    """Analysis parameters model"""
    # Data filtering and subsampling parameters (user controlled)
    filter_genes_min_cells: Optional[Annotated[int, Field(gt=0)]] = None  # Filter genes expressed in < N cells (None = auto: 3 for regular data, adaptive for MERFISH)
    filter_cells_min_genes: Optional[Annotated[int, Field(gt=0)]] = None  # Filter cells expressing < N genes (None = auto: 200 for regular data, adaptive for MERFISH)
    subsample_spots: Optional[Annotated[int, Field(gt=0, le=50000)]] = None  # Subsample to N spots (None = no subsampling)
    subsample_genes: Optional[Annotated[int, Field(gt=0, le=50000)]] = None  # Keep top N variable genes (None = keep all filtered genes)
    subsample_random_seed: int = 42  # Random seed for subsampling

    # Normalization and scaling parameters
    normalization: Literal["log", "sct", "none"] = "log"
    scale: bool = True
    n_hvgs: Annotated[int, Field(gt=0, le=5000)] = 2000
    n_pcs: Annotated[int, Field(gt=0, le=100)] = 30


class VisualizationParameters(BaseModel):
    """Visualization parameters model"""
    feature: Optional[str] = None
    plot_type: Literal[
        "spatial", "heatmap", "violin", "umap",
        "spatial_domains", "cell_communication", "deconvolution",
        "trajectory", "spatial_analysis", "multi_gene", "lr_pairs", "gene_correlation",
        "gaston_isodepth", "gaston_domains", "gaston_genes"
    ] = "spatial"
    colormap: str = "viridis"

    # Spatial analysis visualization parameters
    analysis_sub_type: Optional[Literal["neighborhood", "co_occurrence", "ripley", "moran", "centrality", "getis_ord"]] = Field(
        None, description="Sub-type for spatial_analysis plot type. Determines which spatial analysis result to visualize."
    )
    cluster_key: Optional[str] = Field(None, description="Cluster key for spatial analysis visualization (e.g., 'leiden')")  # For spatial analysis visualization

    # Multi-gene visualization parameters
    features: Optional[List[str]] = None  # Multiple features for multi-panel plots
    multi_panel: bool = False  # Whether to create multi-panel plots
    panel_layout: Optional[Tuple[int, int]] = None  # (rows, cols) - auto-determined if None

    # Ligand-receptor pair parameters
    lr_pairs: Optional[List[Tuple[str, str]]] = None  # List of (ligand, receptor) pairs
    lr_database: str = "cellchat"  # Database for LR pairs

    # Gene correlation parameters
    correlation_method: Literal["pearson", "spearman", "kendall"] = "pearson"
    show_correlation_stats: bool = True

    # Figure parameters
    figure_size: Optional[Tuple[int, int]] = None  # (width, height) - auto-determined if None
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
    basis: Optional[str] = None  # Basis for trajectory visualization (e.g., 'spatial', 'umap', 'pca')

    # Legacy parameters (for backward compatibility)
    show_deconvolution: bool = False  # Whether to show deconvolution results
    n_cell_types: Annotated[int, Field(gt=0, le=10)] = 4  # Number of top cell types to show


class AnnotationParameters(BaseModel):
    """Cell type annotation parameters model"""
    method: Literal["marker_genes", "correlation", "supervised", "popv", "gptcelltype", "scrgcl", "tangram"] = "marker_genes"
    marker_genes: Optional[Dict[str, List[str]]] = None
    reference_data: Optional[str] = None
    reference_data_id: Optional[str] = None  # For Tangram method - ID of reference single-cell dataset
    training_genes: Optional[List[str]] = None  # For Tangram method - genes to use for mapping
    num_epochs: int = 500  # For Tangram method - number of training epochs
    mode: Literal["cells", "clusters"] = "cells"  # For Tangram method - mapping mode
    cluster_label: Optional[str] = None  # For Tangram method - cluster label in reference data


class SpatialAnalysisParameters(BaseModel):
    """Spatial analysis parameters model"""
    analysis_type: Literal["neighborhood", "co_occurrence", "ripley", "moran", "centrality", "getis_ord"] = "neighborhood"
    cluster_key: str = "leiden"
    n_neighbors: Annotated[int, Field(gt=0)] = 15

    # Getis-Ord Gi* specific parameters
    getis_ord_genes: Optional[List[str]] = None  # Specific genes to analyze (if None, use highly variable genes)
    getis_ord_n_genes: Annotated[int, Field(gt=0, le=100)] = 20  # Number of top highly variable genes to analyze
    getis_ord_correction: Literal["bonferroni", "fdr_bh", "none"] = "fdr_bh"  # Multiple testing correction
    getis_ord_alpha: Annotated[float, Field(gt=0.0, le=1.0)] = 0.05  # Significance threshold

    include_image: bool = True
    image_dpi: Annotated[int, Field(gt=0, le=300)] = 100
    image_format: Literal["png", "jpg"] = "png"


class RNAVelocityParameters(BaseModel):
    """RNA velocity analysis parameters model"""
    mode: Literal["deterministic", "stochastic", "dynamical"] = "stochastic"
    n_pcs: Annotated[int, Field(gt=0, le=100)] = 30
    basis: str = "spatial"
    color: Optional[str] = None
    reference_data_id: Optional[str] = None  # For SIRV method
    labels: Optional[List[str]] = None  # For SIRV method


class TrajectoryParameters(BaseModel):
    """Trajectory analysis parameters model"""
    method: Literal["cellrank", "palantir"] = "cellrank"
    spatial_weight: Annotated[float, Field(ge=0.0, le=1.0)] = 0.5  # Spatial information weight
    root_cells: Optional[List[str]] = None  # For Palantir method
    
    # CellRank specific parameters
    cellrank_kernel_weights: Tuple[float, float] = (0.8, 0.2)  # (velocity_weight, connectivity_weight)
    cellrank_n_states: Annotated[int, Field(gt=0, le=20)] = 5  # Number of macrostates for CellRank
    
    # Fallback control
    allow_fallback_to_dpt: bool = True  # Whether to fall back to DPT if other methods fail


class IntegrationParameters(BaseModel):
    """Sample integration parameters model"""
    method: Literal["harmony", "bbknn", "scanorama", "mnn"] = "harmony"
    batch_key: str = "batch"  # Batch information key
    n_pcs: Annotated[int, Field(gt=0, le=100)] = 30  # Number of principal components for integration
    align_spatial: bool = True  # Whether to align spatial coordinates
    reference_batch: Optional[str] = None  # Reference batch for spatial alignment


class DeconvolutionParameters(BaseModel):
    """Spatial deconvolution parameters model"""
    method: Literal["cell2location", "spotiphy", "rctd"] = "cell2location"
    reference_data_id: Optional[str] = None  # Reference single-cell data for deconvolution
    cell_type_key: str = "cell_type"  # Key in reference data for cell type information
    n_top_genes: Annotated[int, Field(gt=0, le=5000)] = 2000  # Number of top genes to use
    use_gpu: bool = False  # Whether to use GPU for cell2location or spotiphy
    n_epochs: Annotated[int, Field(gt=0)] = 10000  # Number of epochs for cell2location or spotiphy
    n_cells_per_spot: Optional[int] = None  # Expected number of cells per spot
    reference_profiles: Optional[Dict[str, List[float]]] = None  # Reference expression profiles

    # Spotiphy specific parameters
    spotiphy_batch_prior: Annotated[float, Field(gt=0)] = 2.0  # Parameter of the prior distribution of the batch effect
    spotiphy_adam_lr: Annotated[float, Field(gt=0)] = 0.003  # Learning rate for Adam optimizer in Spotiphy
    spotiphy_adam_betas: Tuple[float, float] = (0.95, 0.999)  # Beta parameters for Adam optimizer in Spotiphy

    # RCTD specific parameters
    rctd_mode: Literal["full", "doublet", "multi"] = "full"  # RCTD mode
    max_cores: Annotated[int, Field(gt=0, le=16)] = 4  # Maximum number of cores to use
    rctd_confidence_threshold: Annotated[float, Field(gt=0)] = 10.0  # Confidence threshold for cell type assignment
    rctd_doublet_threshold: Annotated[float, Field(gt=0)] = 25.0  # Threshold for doublet detection


class SpatialDomainParameters(BaseModel):
    """Spatial domain identification parameters model"""
    method: Literal["spagcn", "leiden", "louvain"] = "spagcn"
    n_domains: Annotated[int, Field(gt=0, le=50)] = 7  # Number of spatial domains to identify

    # SpaGCN specific parameters
    spagcn_s: Annotated[float, Field(gt=0.0)] = 1.0  # Weight given to histology in SpaGCN
    spagcn_b: Annotated[int, Field(gt=0)] = 49  # Area of each spot when extracting color intensity
    spagcn_p: Annotated[float, Field(ge=0.0, le=1.0)] = 0.5  # Percentage of total expression contributed by neighborhoods
    spagcn_use_histology: bool = True  # Whether to use histology image in SpaGCN
    spagcn_random_seed: int = 100  # Random seed for SpaGCN

    # General clustering parameters
    resolution: float = 0.5  # Resolution for leiden/louvain clustering
    use_highly_variable: bool = True  # Whether to use highly variable genes only
    refine_domains: bool = True  # Whether to refine spatial domains using spatial smoothing
    
    # Clustering-specific parameters for leiden/louvain methods
    cluster_n_neighbors: Optional[Annotated[int, Field(gt=0)]] = None  # Number of neighbors for clustering (default: 15)
    cluster_spatial_weight: Optional[Annotated[float, Field(ge=0.0, le=1.0)]] = None  # Weight for spatial information (default: 0.3)



class SpatialVariableGenesParameters(BaseModel):
    """GASTON spatial variable genes identification parameters model"""
    # Preprocessing parameters
    preprocessing_method: Literal["glmpca", "pearson_residuals"] = "glmpca"
    n_components: Annotated[int, Field(gt=0, le=50)] = 10  # Number of components for dimensionality reduction

    # Neural network architecture parameters
    spatial_hidden_layers: List[Annotated[int, Field(gt=0, le=1000)]] = [50]  # Architecture for spatial embedding f_S
    expression_hidden_layers: List[Annotated[int, Field(gt=0, le=1000)]] = [10]  # Architecture for expression function f_A

    # Training parameters
    epochs: Annotated[int, Field(gt=0, le=50000)] = 10000  # Number of training epochs
    learning_rate: Annotated[float, Field(gt=0.0, le=1.0)] = 0.001  # Learning rate
    optimizer: Literal["adam", "sgd", "adagrad"] = "adam"  # Optimizer type
    batch_size: Optional[Annotated[int, Field(gt=0, le=10000)]] = None  # Batch size (None for full batch)

    # Positional encoding parameters
    use_positional_encoding: bool = False  # Whether to use positional encoding
    embedding_size: Annotated[int, Field(gt=0, le=20)] = 4  # Positional encoding embedding size
    sigma: Annotated[float, Field(gt=0.0, le=1.0)] = 0.2  # Positional encoding sigma parameter

    # Model saving parameters
    checkpoint_interval: Annotated[int, Field(gt=0, le=5000)] = 500  # Save model every N epochs
    random_seed: Annotated[int, Field(ge=0, le=1000)] = 0  # Random seed for reproducibility

    # Analysis parameters
    n_domains: Annotated[int, Field(gt=0, le=20)] = 5  # Number of spatial domains to identify
    num_bins: Annotated[int, Field(gt=10, le=200)] = 70  # Number of bins for isodepth binning
    umi_threshold: Annotated[int, Field(gt=0)] = 500  # Minimum UMI count threshold for genes

    # Gene classification parameters
    continuous_quantile: Annotated[float, Field(gt=0.0, le=1.0)] = 0.9  # Quantile threshold for continuous genes
    discontinuous_quantile: Annotated[float, Field(gt=0.0, le=1.0)] = 0.9  # Quantile threshold for discontinuous genes
    pvalue_threshold: Annotated[float, Field(gt=0.0, le=1.0)] = 0.1  # P-value threshold for slope significance
    zero_fit_threshold: Annotated[int, Field(ge=0)] = 0  # Minimum non-zero count per domain

    # Poisson regression parameters
    isodepth_mult_factor: Annotated[float, Field(gt=0.0)] = 1.0  # Scaling factor for isodepth values
    regularization: Annotated[float, Field(ge=0.0)] = 0.0  # Regularization parameter


class CellCommunicationParameters(BaseModel):
    """Cell-cell communication analysis parameters model"""
    method: Literal["liana"] = "liana"  # Only LIANA+ is supported

    # General parameters
    species: Literal["human", "mouse", "zebrafish"] = "human"  # Species for ligand-receptor database
    min_cells: Annotated[int, Field(ge=0)] = 3  # Minimum cells expressing ligand or receptor

    # LIANA+ specific parameters
    liana_resource: Literal["consensus", "cellchat", "cellphonedb", "connectome", "omnipath"] = "consensus"  # LR database resource
    liana_local_metric: Literal["cosine", "pearson", "spearman", "jaccard"] = "cosine"  # Local spatial metric
    liana_global_metric: Literal["morans", "lee"] = "morans"  # Global spatial metric
    liana_n_perms: Annotated[int, Field(gt=0)] = 100  # Number of permutations for LIANA
    liana_nz_prop: Annotated[float, Field(gt=0.0, le=1.0)] = 0.2  # Minimum expression proportion
    liana_bandwidth: Optional[int] = None  # Bandwidth for spatial connectivity (auto-determined if None)
    liana_cutoff: Annotated[float, Field(gt=0.0, le=1.0)] = 0.1  # Cutoff for spatial connectivity
    perform_spatial_analysis: bool = True  # Whether to perform spatial bivariate analysis (vs cluster-based)

    # Custom ligand-receptor pairs (for advanced users)
    custom_lr_pairs: Optional[List[Tuple[str, str]]] = None  # Custom LR pairs as (ligand, receptor) tuples

    # Visualization parameters
    include_image: bool = True  # Whether to include visualization
    plot_top_pairs: Annotated[int, Field(gt=0, le=20)] = 6  # Number of top LR pairs to visualize
    image_dpi: Annotated[int, Field(gt=0, le=300)] = 100  # DPI for output image
    image_format: Literal["png", "jpg"] = "png"  # Image format