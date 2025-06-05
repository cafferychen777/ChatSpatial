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
    """Spatial variable genes identification parameters model"""
    method: Literal["spatialde", "spark", "trendsceek"] = "spatialde"

    # General parameters
    n_top_genes: Annotated[int, Field(gt=0, le=10000)] = 100  # Number of top spatial genes to return
    significance_threshold: Annotated[float, Field(gt=0.0, le=1.0)] = 0.05  # FDR threshold for significance

    # SpatialDE specific parameters
    spatialde_kernel: Literal["SE", "linear", "periodic"] = "SE"  # Kernel type for SpatialDE
    spatialde_normalize: bool = True  # Whether to normalize and stabilize expression data
    spatialde_regress_out_total_counts: bool = True  # Whether to regress out total counts

    # SPARK specific parameters (if implemented)
    spark_num_core: int = 1  # Number of cores for SPARK
    spark_verbose: bool = False  # Verbose output for SPARK

    # Automatic Expression Histology (AEH) parameters
    perform_aeh: bool = False  # Whether to perform automatic expression histology
    aeh_n_patterns: Annotated[int, Field(gt=0, le=20)] = 5  # Number of spatial patterns for AEH
    aeh_length_scale: Optional[float] = None  # Length scale for AEH (auto-determined if None)

    # Visualization parameters
    include_image: bool = True  # Whether to include visualization
    plot_top_genes: Annotated[int, Field(gt=0, le=20)] = 6  # Number of top genes to visualize
    image_dpi: Annotated[int, Field(gt=0, le=300)] = 100  # DPI for output image
    image_format: Literal["png", "jpg"] = "png"  # Image format


class CellCommunicationParameters(BaseModel):
    """Cell-cell communication analysis parameters model"""
    method: Literal["commot", "spatialdm", "cellphonedb"] = "commot"

    # General parameters
    species: Literal["human", "mouse", "zebrafish"] = "human"  # Species for ligand-receptor database
    database: Literal["cellchat", "cellphonedb", "user"] = "cellchat"  # Ligand-receptor database
    min_cells: Annotated[int, Field(ge=0)] = 3  # Minimum cells expressing ligand or receptor

    # COMMOT specific parameters
    commot_dis_thr: Annotated[float, Field(gt=0)] = 200.0  # Distance threshold for COMMOT
    commot_heteromeric: bool = True  # Whether to consider heteromeric complexes
    commot_n_permutations: Annotated[int, Field(gt=0)] = 100  # Number of permutations for significance testing

    # SpatialDM specific parameters
    spatialdm_l: Optional[float] = None  # RBF kernel parameter (auto-determined if None)
    spatialdm_cutoff: Annotated[float, Field(gt=0.0, le=1.0)] = 0.1  # Weight cutoff for SpatialDM
    spatialdm_n_neighbors: Optional[int] = None  # Number of neighbors (auto-determined if None)
    spatialdm_n_permutations: Annotated[int, Field(gt=0)] = 1000  # Number of permutations
    spatialdm_method: Literal["z-score", "permutation", "both"] = "z-score"  # Statistical method
    spatialdm_fdr: bool = True  # Whether to apply FDR correction
    spatialdm_threshold: Annotated[float, Field(gt=0.0, le=1.0)] = 0.1  # Significance threshold

    # Analysis options
    perform_global_analysis: bool = True  # Whether to perform global LR pair selection
    perform_local_analysis: bool = True  # Whether to perform local spot analysis
    identify_communication_patterns: bool = False  # Whether to identify communication patterns

    # Custom ligand-receptor pairs (if database="user")
    custom_lr_pairs: Optional[List[Tuple[str, str]]] = None  # Custom LR pairs as (ligand, receptor) tuples

    # Visualization parameters
    include_image: bool = True  # Whether to include visualization
    plot_top_pairs: Annotated[int, Field(gt=0, le=20)] = 6  # Number of top LR pairs to visualize
    image_dpi: Annotated[int, Field(gt=0, le=300)] = 100  # DPI for output image
    image_format: Literal["png", "jpg"] = "png"  # Image format