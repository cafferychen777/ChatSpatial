"""
Data standards and schema definitions for ChatSpatial.

This module defines the canonical data structure standards for AnnData objects
used across all ChatSpatial tools. Following Linus's principle: "Bad programmers 
worry about the code. Good programmers worry about data structures."

The goal is to eliminate all special cases and have one unified data structure
that works across all 17 tool modules.
"""

# Removed type imports for Python 2.7 compatibility
# import numpy as np


class ChatSpatialDataStandards(object):
    """
    Canonical data structure standards for ChatSpatial AnnData objects.
    
    These standards define the expected field names and formats that all
    tools should use. No exceptions, no special cases.
    """
    
    def __init__(self):
        # === Core field naming standards ===
        self.spatial_key = "spatial"        # adata.obsm['spatial'] - (n_obs, 2) float array of x,y coordinates
        self.cell_type_key = "cell_type"    # adata.obs['cell_type'] - categorical cell type annotations
        self.cluster_key = "leiden"         # adata.obs['leiden'] - categorical cluster assignments  
        self.batch_key = "batch"            # adata.obs['batch'] - categorical batch/sample identifiers
        
        # === Data type requirements ===
        self.spatial_dtype = "float64"      # Required dtype for spatial coordinates
        self.spatial_dims = 2               # Required number of spatial dimensions (x, y)
        
        # === Gene naming standards ===
        self.gene_name_format = "symbol"    # Preferred gene identifier format
        
        # === Alternative field names (for input compatibility) ===
        self.spatial_alternatives = None
        self.cell_type_alternatives = None
        self.cluster_alternatives = None
        self.batch_alternatives = None
        
        self._initialize_alternatives()
    
    def _initialize_alternatives(self):
        """Initialize alternative field name sets."""
        if self.spatial_alternatives is None:
            self.spatial_alternatives = {
                "spatial", "X_spatial", "coordinates", "coords", 
                "spatial_coords", "positions"
            }
        
        if self.cell_type_alternatives is None:
            self.cell_type_alternatives = {
                "cell_type", "celltype", "cell_types", "annotation", 
                "cell_annotation", "predicted_celltype"
            }
            
        if self.cluster_alternatives is None:
            self.cluster_alternatives = {
                "leiden", "louvain", "clusters", "cluster", "clustering",
                "cluster_labels", "spatial_domains"
            }
            
        if self.batch_alternatives is None:
            self.batch_alternatives = {
                "batch", "sample", "dataset", "experiment", "replicate",
                "batch_id", "sample_id"
            }


# Global instance - the single source of truth
CHATSPATIAL_STANDARDS = ChatSpatialDataStandards()


class DataValidationResult(object):
    """Result of data structure validation."""
    
    def __init__(self, passed=True, errors=None, warnings=None, suggestions=None, standardized_fields=None):
        self.passed = passed
        self.errors = errors if errors is not None else []
        self.warnings = warnings if warnings is not None else []
        self.suggestions = suggestions if suggestions is not None else []
        self.standardized_fields = standardized_fields if standardized_fields is not None else {}
    
    def add_error(self, message):
        """Add an error message and mark validation as failed."""
        self.errors.append(message)
        self.passed = False
    
    def add_warning(self, message):
        """Add a warning message."""
        self.warnings.append(message)
    
    def add_suggestion(self, message):
        """Add a suggestion for fixing issues."""
        self.suggestions.append(message)


class RequiredDataFields(object):
    """Defines which data fields are required for different analysis types."""
    
    def __init__(self):
        # Core fields required by most tools
        self.CORE_SPATIAL = None
        self.CORE_EXPRESSION = None
        
        # Tool-specific requirements
        self.CELL_COMMUNICATION = None
        self.DECONVOLUTION = None
        self.SPATIAL_ANALYSIS = None
        self.VISUALIZATION = None
        self.TRAJECTORY = None
        self.INTEGRATION = None
        
        self._initialize_requirements()
    
    def _initialize_requirements(self):
        """Initialize required field sets."""
        if self.CORE_SPATIAL is None:
            self.CORE_SPATIAL = {"spatial_coords", "gene_expression"}
            
        if self.CORE_EXPRESSION is None:
            self.CORE_EXPRESSION = {"X", "var_names", "obs_names"}
            
        if self.CELL_COMMUNICATION is None:
            self.CELL_COMMUNICATION = self.CORE_SPATIAL | {"cell_type", "spatial_neighbors"}
            
        if self.DECONVOLUTION is None:
            self.DECONVOLUTION = self.CORE_SPATIAL | {"reference_data", "common_genes"}
            
        if self.SPATIAL_ANALYSIS is None:
            self.SPATIAL_ANALYSIS = self.CORE_SPATIAL | {"cluster_labels", "spatial_neighbors"}
            
        if self.VISUALIZATION is None:
            self.VISUALIZATION = self.CORE_SPATIAL.copy()  # Most flexible requirements
            
        if self.TRAJECTORY is None:
            self.TRAJECTORY = self.CORE_SPATIAL | {"cluster_labels", "highly_variable_genes"}
            
        if self.INTEGRATION is None:
            self.INTEGRATION = self.CORE_SPATIAL | {"batch_info", "highly_variable_genes"}


# Global instance for required fields
REQUIRED_FIELDS = RequiredDataFields()


def get_field_mapping():
    """
    Get the mapping from alternative field names to standard field names.
    
    This is the single source of truth for field name standardization.
    All tools should use this mapping to convert input data to standard format.
    
    Returns:
        Dictionary mapping alternative names to standard names
    """
    standards = CHATSPATIAL_STANDARDS
    mapping = {}
    
    # Spatial coordinates mapping
    for alt in standards.spatial_alternatives:
        mapping[alt] = standards.spatial_key
    
    # Cell type mapping  
    for alt in standards.cell_type_alternatives:
        mapping[alt] = standards.cell_type_key
        
    # Cluster mapping
    for alt in standards.cluster_alternatives:
        mapping[alt] = standards.cluster_key
        
    # Batch mapping
    for alt in standards.batch_alternatives:
        mapping[alt] = standards.batch_key
    
    return mapping


def get_obsm_standard_keys():
    """Get list of standard keys that should be in adata.obsm."""
    return [CHATSPATIAL_STANDARDS.spatial_key]


def get_obs_standard_keys():  
    """Get list of standard keys that should be in adata.obs."""
    standards = CHATSPATIAL_STANDARDS
    return [
        standards.cell_type_key,
        standards.cluster_key, 
        standards.batch_key
    ]


def get_data_format_schema():
    """
    Get the complete data format schema for ChatSpatial AnnData objects.
    
    This defines the expected structure, types, and constraints for all
    data fields used by ChatSpatial tools.
    
    Returns:
        Schema dictionary with field specifications
    """
    standards = CHATSPATIAL_STANDARDS
    
    schema = {
        "obsm": {
            standards.spatial_key: {
                "type": "ndarray",
                "dtype": standards.spatial_dtype,
                "shape": ("n_obs", standards.spatial_dims),
                "description": "Spatial coordinates (x, y) for each observation",
                "required": True
            }
        },
        "obs": {
            standards.cell_type_key: {
                "type": "categorical", 
                "description": "Cell type annotations",
                "required": False,
                "constraints": ["must_be_categorical", "no_missing_allowed"]
            },
            standards.cluster_key: {
                "type": "categorical",
                "description": "Cluster assignments (e.g., from Leiden clustering)",
                "required": False,
                "constraints": ["must_be_categorical", "no_missing_allowed"]
            },
            standards.batch_key: {
                "type": "categorical",
                "description": "Batch/sample identifiers", 
                "required": False,
                "constraints": ["must_be_categorical"]
            }
        },
        "X": {
            "type": "sparse_or_dense",
            "dtype": "numeric",
            "shape": ("n_obs", "n_vars"),
            "description": "Gene expression matrix",
            "required": True,
            "constraints": ["non_negative", "finite_values"]
        },
        "var": {
            "gene_names": {
                "type": "index",
                "description": "Gene identifiers (symbols preferred)",
                "required": True,
                "constraints": ["unique", "non_empty", "no_whitespace"]
            }
        },
        "uns": {
            "spatial": {
                "type": "dict",
                "description": "Spatial metadata (for Visium data)",
                "required": False
            }
        }
    }
    
    return schema


def describe_standards():
    """
    Return a human-readable description of ChatSpatial data standards.
    
    This is useful for error messages and documentation.
    """
    standards = CHATSPATIAL_STANDARDS
    
    description = """
ChatSpatial Data Structure Standards
===================================

Core Field Names:
- Spatial coordinates: adata.obsm['""" + standards.spatial_key + """'] (required)
- Cell type annotations: adata.obs['""" + standards.cell_type_key + """'] (optional)  
- Cluster assignments: adata.obs['""" + standards.cluster_key + """'] (optional)
- Batch information: adata.obs['""" + standards.batch_key + """'] (optional)

Data Requirements:
- Spatial coordinates must be float64 array of shape (n_obs, 2)
- Categorical fields must be pandas Categorical type
- Gene names should be unique symbols when possible
- Expression matrix must contain finite, non-negative values

Alternative Names (automatically converted):
- Spatial: """ + ', '.join(sorted(standards.spatial_alternatives)) + """
- Cell types: """ + ', '.join(sorted(standards.cell_type_alternatives)) + """
- Clusters: """ + ', '.join(sorted(standards.cluster_alternatives)) + """
- Batches: """ + ', '.join(sorted(standards.batch_alternatives)) + """

All ChatSpatial tools expect data in this standard format.
Use chatspatial.utils.data_adapter.standardize_adata() to convert your data.
"""
    
    return description.strip()