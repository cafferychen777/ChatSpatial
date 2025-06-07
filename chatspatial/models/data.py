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
    plot_type: Literal["spatial", "heatmap", "violin", "umap"] = "spatial"
    colormap: str = "viridis"
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
    analysis_type: Literal["neighborhood", "co_occurrence", "ripley", "moran", "centrality"] = "neighborhood"
    cluster_key: str = "leiden"
    n_neighbors: Annotated[int, Field(gt=0)] = 15
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


class IntegrationParameters(BaseModel):
    """Sample integration parameters model"""
    method: Literal["harmony", "bbknn", "scanorama", "mnn"] = "harmony"
    batch_key: str = "batch"  # Batch information key
    n_pcs: Annotated[int, Field(gt=0, le=100)] = 30  # Number of principal components for integration
    align_spatial: bool = True  # Whether to align spatial coordinates
    reference_batch: Optional[str] = None  # Reference batch for spatial alignment


class DeconvolutionParameters(BaseModel):
    """Spatial deconvolution parameters model"""
    method: Literal["cell2location", "nnls", "spotiphy"] = "nnls"
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


class SpatialDomainParameters(BaseModel):
    """Spatial domain identification parameters model"""
    method: Literal["stagate", "spagcn", "leiden", "louvain"] = "stagate"
    n_domains: Annotated[int, Field(gt=0, le=50)] = 7  # Number of spatial domains to identify

    # STAGATE specific parameters
    stagate_alpha: Annotated[float, Field(ge=0.0, le=1.0)] = 0.0  # Weight of cell type-aware spatial neighbor network
    stagate_hidden_dims: List[int] = [512, 30]  # Hidden dimensions for STAGATE encoder
    stagate_n_epochs: Annotated[int, Field(gt=0)] = 500  # Number of epochs for STAGATE training
    stagate_lr: float = 0.0001  # Learning rate for STAGATE
    stagate_random_seed: int = 2020  # Random seed for STAGATE

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

    # Visualization parameters
    include_image: bool = True  # Whether to include visualization
    image_dpi: Annotated[int, Field(gt=0, le=300)] = 100  # DPI for output image
    image_format: Literal["png", "jpg"] = "png"  # Image format



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
    spatial_domain_threshold: Annotated[float, Field(gt=0.0, le=1.0)] = 0.95  # Quantile threshold for spatial domain identification
    gradient_threshold: Annotated[float, Field(gt=0.0, le=1.0)] = 0.95  # Quantile threshold for gradient detection

    # Visualization parameters
    include_visualization: bool = True  # Whether to generate visualizations
    plot_isodepth_map: bool = True  # Whether to plot isodepth map
    plot_spatial_domains: bool = True  # Whether to plot spatial domains
    plot_top_genes: Annotated[int, Field(gt=0, le=20)] = 6  # Number of top genes to visualize
    image_dpi: Annotated[int, Field(gt=0, le=300)] = 100  # DPI for output images
    image_format: Literal["png", "jpg"] = "png"  # Image format


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