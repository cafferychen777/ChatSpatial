"""
Unified constants for spatial transcriptomics analysis.

This module centralizes commonly used default parameters to improve
code maintainability and consistency across the codebase.
"""

# Analysis default parameters
DEFAULT_N_NEIGHBORS = 30
DEFAULT_RESOLUTION = 0.5
DEFAULT_RANDOM_STATE = 42
DEFAULT_N_PCS = 50

# Preprocessing parameters
DEFAULT_TARGET_SUM = 1e4
MAX_SCALE_VALUE = 10
MIN_NEIGHBORS = 3
MAX_NEIGHBORS_RATIO = 0.1

# Clustering parameters
MIN_KMEANS_CLUSTERS = 2
CLUSTERING_RESOLUTIONS = {
    'small': 0.4,   # < 100 cells
    'medium': 0.6,  # 100-500 cells  
    'large': 0.8    # > 500 cells
}

# Data type detection thresholds
MERFISH_GENE_THRESHOLD = 200

# Visualization parameters
MAX_TSNE_PCA_COMPONENTS = 50

# Spatial analysis parameters  
DEFAULT_SPATIAL_KEY = 'spatial'
DEFAULT_N_STATES = 5  # For trajectory analysis

# File size and computation limits
MAX_CELLS_FOR_FULL_ANALYSIS = 10000
MAX_GENES_FOR_HVG_SELECTION = 2000

# sc-type integration
SCTYPE_VALID_TISSUES = {
    "Adrenal", "Brain", "Eye", "Heart", "Hippocampus", "Immune system",
    "Intestine", "Kidney", "Liver", "Lung", "Muscle", "Pancreas", 
    "Placenta", "Spleen", "Stomach", "Thymus"
}

# Error handling
DEFAULT_TIMEOUT_SECONDS = 300  # 5 minutes for long-running analyses