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


class AnnotationParameters(BaseModel):
    """Cell type annotation parameters model"""
    method: Literal["marker_genes", "correlation", "supervised", "popv", "gptcelltype", "scrgcl"] = "marker_genes"
    marker_genes: Optional[Dict[str, List[str]]] = None
    reference_data: Optional[str] = None


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