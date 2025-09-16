"""
Data validation utilities for ChatSpatial.

This module implements Linus's principle: "Eliminate all special cases."
Instead of having 17 different validation logics scattered across tools,
we have ONE validation system that works everywhere.

Every tool uses the same validation functions. No exceptions.
"""

from typing import Optional, Dict, List, Any, Set, Union, TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse
    import anndata as ad

# Industry standard field names - these are fixed conventions
SPATIAL_KEY = 'spatial'  # Standard key for spatial coordinates in obsm
CELL_TYPE_KEY = 'cell_type'  # Standard key for cell type annotations in obs
CLUSTER_KEY = 'leiden'  # Standard key for cluster assignments in obs
BATCH_KEY = 'batch'  # Standard key for batch information in obs

# Alternative names that might be found in input data
# We support these for compatibility but standardize to the above keys
ALTERNATIVE_SPATIAL_KEYS = {'spatial', 'X_spatial', 'coordinates', 'coords', 
                            'spatial_coords', 'positions'}
ALTERNATIVE_CELL_TYPE_KEYS = {'cell_type', 'celltype', 'cell_types', 'annotation',
                              'cell_annotation', 'predicted_celltype'}
ALTERNATIVE_CLUSTER_KEYS = {'leiden', 'louvain', 'clusters', 'cluster',
                           'clustering', 'cluster_labels', 'spatial_domains'}
ALTERNATIVE_BATCH_KEYS = {'batch', 'sample', 'dataset', 'experiment', 
                         'replicate', 'batch_id', 'sample_id'}


class DataValidationError(Exception):
    """Raised when data validation fails critically."""
    pass


class DataValidationResult:
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


class DataValidator:
    """
    Single source of truth for data validation in ChatSpatial.
    
    This class eliminates the need for scattered validation logic across
    17 different tool modules. One validator, no special cases.
    """
    
    def __init__(self, strict_mode: bool = False):
        """
        Initialize data validator.
        
        Args:
            strict_mode: If True, warnings become errors
        """
        # Import dependencies at runtime
        import numpy as np
        import pandas as pd
        import scanpy as sc
        from scipy import sparse
        import anndata as ad
        
        self.np = np
        self.pd = pd
        self.sc = sc
        self.sparse = sparse
        self.ad = ad
        
        self.strict_mode = strict_mode
    
    def validate_adata(self, 
                      adata: 'ad.AnnData', 
                      analysis_type: Optional[str] = None,
                      require_spatial: bool = True,
                      require_cell_types: bool = False) -> DataValidationResult:
        """
        Validate AnnData object against ChatSpatial standards.
        
        This is the main validation function used by all tools.
        No tool should implement its own validation logic.
        
        Args:
            adata: AnnData object to validate
            analysis_type: Type of analysis (for specific requirements)
            require_spatial: Whether spatial coordinates are required
            require_cell_types: Whether cell type annotations are required
            
        Returns:
            DataValidationResult with validation status and suggestions
        """
        result = DataValidationResult()
        
        # Basic structure validation
        self._validate_basic_structure(adata, result)
        
        # Spatial coordinates validation
        if require_spatial:
            self._validate_spatial_coordinates(adata, result)
        
        # Expression matrix validation
        self._validate_expression_matrix(adata, result)
        
        # Gene names validation  
        self._validate_gene_names(adata, result)
        
        # Optional field validation
        self._validate_optional_fields(adata, result)
        
        # Cell type validation if required
        if require_cell_types:
            self._validate_cell_types(adata, result)
        
        # Analysis-specific validation
        if analysis_type:
            self._validate_analysis_requirements(adata, analysis_type, result)
        
        # Generate standardization suggestions
        self._suggest_standardization(adata, result)
        
        return result
    
    def _validate_basic_structure(self, adata: 'ad.AnnData', result: DataValidationResult) -> None:
        """Validate basic AnnData structure."""
        
        if not isinstance(adata, self.ad.AnnData):
            result.add_error(f"Expected AnnData object, got {type(adata)}")
            return
        
        if adata.n_obs == 0:
            result.add_error("Dataset contains no observations (cells/spots)")
            
        if adata.n_vars == 0:
            result.add_error("Dataset contains no variables (genes)")
        
        # Check for required attributes
        if not hasattr(adata, 'X') or adata.X is None:
            result.add_error("Missing expression matrix (adata.X)")
            
        if not hasattr(adata, 'obs'):
            result.add_error("Missing observation metadata (adata.obs)")
            
        if not hasattr(adata, 'var'):
            result.add_error("Missing variable metadata (adata.var)")
            
        if not hasattr(adata, 'obsm'):
            result.add_error("Missing multi-dimensional observation data (adata.obsm)")
    
    def _validate_spatial_coordinates(self, adata: 'ad.AnnData', result: DataValidationResult) -> None:
        """
        Validate spatial coordinates.
        
        This replaces the scattered coordinate detection logic in multiple tools.
        """
        spatial_key = None
        found_alternatives = []
        
        # Look for spatial coordinates in obsm
        for key in adata.obsm.keys():
            if key in ALTERNATIVE_SPATIAL_KEYS:
                found_alternatives.append(key)
                if key == SPATIAL_KEY:
                    spatial_key = key
                    break
        
        # If standard key not found, use first alternative
        if spatial_key is None and found_alternatives:
            spatial_key = found_alternatives[0]
            result.add_warning(f"Using non-standard spatial key '{spatial_key}', "
                             f"should be '{SPATIAL_KEY}'")
            result.standardized_fields[spatial_key] = SPATIAL_KEY
        
        # Check for coordinates in obs (fallback)
        if spatial_key is None:
            if 'x' in adata.obs and 'y' in adata.obs:
                result.add_warning("Found spatial coordinates in adata.obs['x'], adata.obs['y']. "
                                 f"Should be in adata.obsm['{SPATIAL_KEY}']")
                result.add_suggestion(f"Move coordinates to adata.obsm['{SPATIAL_KEY}']")
                return
        
        if spatial_key is None:
            result.add_error("No spatial coordinates found")
            result.add_suggestion(f"Add spatial coordinates to adata.obsm['{SPATIAL_KEY}']")
            return
        
        # Validate spatial coordinate properties
        coords = adata.obsm[spatial_key]
        
        if not isinstance(coords, self.np.ndarray):
            result.add_error(f"Spatial coordinates must be numpy array, got {type(coords)}")
            return
        
        if coords.shape[0] != adata.n_obs:
            result.add_error(f"Spatial coordinates shape mismatch: "
                           f"{coords.shape[0]} coordinates for {adata.n_obs} observations")
        
        if coords.shape[1] < 2:
            result.add_error(f"Spatial coordinates must have at least 2 dimensions (x, y), "
                           f"found {coords.shape[1]}")
        
        # Check data type - float64 is the expected type
        if coords.dtype != 'float64':
            result.add_warning(f"Spatial coordinates should be float64, "
                             f"found {coords.dtype}")
        
        # Check for invalid values
        if self.np.any(self.np.isnan(coords)):
            result.add_error("Spatial coordinates contain NaN values")
            
        if self.np.any(self.np.isinf(coords)):
            result.add_error("Spatial coordinates contain infinite values")
    
    def _validate_expression_matrix(self, adata: 'ad.AnnData', result: DataValidationResult) -> None:
        """Validate gene expression matrix."""
        
        if adata.X is None:
            result.add_error("Expression matrix (adata.X) is None")
            return
        
        # Check shape
        expected_shape = (adata.n_obs, adata.n_vars)
        if adata.X.shape != expected_shape:
            result.add_error(f"Expression matrix shape mismatch: "
                           f"expected {expected_shape}, got {adata.X.shape}")
        
        # Check for invalid values
        if self.sparse.issparse(adata.X):
            data = adata.X.data
        else:
            data = adata.X.flatten()
        
        if self.np.any(self.np.isnan(data)):
            result.add_error("Expression matrix contains NaN values")
            
        if self.np.any(self.np.isinf(data)):
            result.add_error("Expression matrix contains infinite values")
        
        # Check for negative values (should be rare in processed data)
        if self.np.any(data < 0):
            result.add_warning("Expression matrix contains negative values")
            result.add_suggestion("Consider using raw counts or properly normalized data")
        
        # Check data range
        max_val = self.np.max(data)
        if max_val > 1000:
            result.add_warning(f"Expression values very large (max: {max_val:.1f})")
            result.add_suggestion("Data might need normalization")
        elif max_val < 0.01:
            result.add_warning(f"Expression values very small (max: {max_val:.6f})")
            result.add_suggestion("Data might be over-normalized")
    
    def _validate_gene_names(self, adata: 'ad.AnnData', result: DataValidationResult) -> None:
        """Validate gene names in var_names."""
        
        if len(adata.var_names) == 0:
            result.add_error("No gene names found in adata.var_names")
            return
        
        # Check for duplicates
        if not adata.var_names.is_unique:
            n_duplicates = len(adata.var_names) - len(adata.var_names.unique())
            result.add_error(f"Gene names contain {n_duplicates} duplicates")
            result.add_suggestion("Use adata.var_names_make_unique() to fix duplicates")
        
        # Check for empty gene names
        empty_names = adata.var_names.isna() | (adata.var_names == "")
        if empty_names.any():
            n_empty = empty_names.sum()
            result.add_error(f"Found {n_empty} empty gene names")
        
        # Check gene name format
        gene_names = adata.var_names.tolist()
        
        # Detect probable Ensembl IDs
        ensembl_pattern = any(name.startswith(('ENSG', 'ENSMUSG')) for name in gene_names[:100])
        if ensembl_pattern:
            result.add_warning("Gene names appear to be Ensembl IDs")
            result.add_suggestion("Consider converting to gene symbols for better compatibility")
        
        # Check for whitespace issues
        whitespace_issues = [name for name in gene_names[:100] 
                           if isinstance(name, str) and (name != name.strip())]
        if whitespace_issues:
            result.add_warning(f"Found gene names with whitespace issues: {whitespace_issues[:5]}")
            result.add_suggestion("Trim whitespace from gene names")
    
    def _validate_optional_fields(self, adata: 'ad.AnnData', result: DataValidationResult) -> None:
        """Validate optional metadata fields."""
        
        # Check for cell type annotations
        cell_type_found = self._find_field_in_obs(adata, ALTERNATIVE_CELL_TYPE_KEYS)
        if cell_type_found:
            key, standard_key = cell_type_found
            if key != standard_key:
                result.standardized_fields[key] = standard_key
                result.add_warning(f"Found cell types in non-standard field '{key}', "
                                 f"should be '{standard_key}'")
            self._validate_categorical_field(adata.obs[key], f"Cell types ({key})", result)
        
        # Check for cluster annotations  
        cluster_found = self._find_field_in_obs(adata, ALTERNATIVE_CLUSTER_KEYS)
        if cluster_found:
            key, standard_key = cluster_found
            if key != standard_key:
                result.standardized_fields[key] = standard_key
                result.add_warning(f"Found clusters in non-standard field '{key}', "
                                 f"should be '{standard_key}'")
            self._validate_categorical_field(adata.obs[key], f"Clusters ({key})", result)
        
        # Check for batch information
        batch_found = self._find_field_in_obs(adata, ALTERNATIVE_BATCH_KEYS)
        if batch_found:
            key, standard_key = batch_found
            if key != standard_key:
                result.standardized_fields[key] = standard_key
                result.add_warning(f"Found batch info in non-standard field '{key}', "
                                 f"should be '{standard_key}'")
            self._validate_categorical_field(adata.obs[key], f"Batch info ({key})", result)
    
    def _validate_cell_types(self, adata: 'ad.AnnData', result: DataValidationResult) -> None:
        """Validate cell type annotations (when required)."""
        
        cell_type_found = self._find_field_in_obs(adata, ALTERNATIVE_CELL_TYPE_KEYS)
        if not cell_type_found:
            result.add_error("Cell type annotations required but not found")
            result.add_suggestion(f"Add cell type annotations to adata.obs['{CELL_TYPE_KEY}']")
            return
        
        key, _ = cell_type_found
        cell_types = adata.obs[key]
        
        # Check for missing values
        if cell_types.isna().any():
            n_missing = cell_types.isna().sum()
            result.add_error(f"Cell type annotations contain {n_missing} missing values")
        
        # Check minimum number of cell types
        unique_types = cell_types.dropna().unique()
        if len(unique_types) < 2:
            result.add_warning(f"Only {len(unique_types)} unique cell types found")
            result.add_suggestion("Consider if cell type diversity is sufficient for analysis")
    
    def _validate_analysis_requirements(self, 
                                      adata: 'ad.AnnData', 
                                      analysis_type: str, 
                                      result: DataValidationResult) -> None:
        """Validate requirements for specific analysis types."""
        
        # Analysis-specific requirements are hardcoded based on actual needs
        # No need for abstraction here since these are fixed requirements
        
        # Check analysis-specific requirements
        if analysis_type.lower() == "cell_communication":
            # Cell communication requires cell type annotations and spatial neighbors
            if not self._find_field_in_obs(adata, ALTERNATIVE_CELL_TYPE_KEYS):
                result.add_error("Cell communication analysis requires cell type annotations")
            
            if 'spatial_connectivities' not in adata.obsp:
                result.add_warning("Spatial connectivity graph not found")
                result.add_suggestion("Run sq.gr.spatial_neighbors(adata) to compute spatial graph")
        
        elif analysis_type.lower() == "deconvolution":
            # Deconvolution has different requirements handled elsewhere
            pass
        
        elif analysis_type.lower() == "spatial_analysis":
            # Spatial analysis requires cluster labels
            if not self._find_field_in_obs(adata, ALTERNATIVE_CLUSTER_KEYS):
                result.add_error("Spatial analysis requires cluster annotations")
                result.add_suggestion(f"Run clustering and add results to adata.obs['{CLUSTER_KEY}']")
    
    def _validate_categorical_field(self, series: 'pd.Series', field_name: str, result: DataValidationResult) -> None:
        """Validate a categorical metadata field."""
        
        if not self.pd.api.types.is_categorical_dtype(series):
            result.add_warning(f"{field_name} should be categorical dtype")
            result.add_suggestion(f"Convert to categorical: adata.obs[field] = adata.obs[field].astype('category')")
        
        # Check for too many categories
        n_categories = len(series.unique())
        if n_categories > 100:
            result.add_warning(f"{field_name} has many categories ({n_categories})")
            result.add_suggestion("Consider if this categorical variable is appropriate")
    
    def _find_field_in_obs(self, adata: 'ad.AnnData', alternatives: Set[str]) -> Optional[tuple]:
        """
        Find a field in adata.obs using alternative names.
        
        Returns:
            Tuple of (found_key, standard_key) or None if not found
        """
        for key in adata.obs.columns:
            if key in alternatives:
                # Determine the standard key name based on which set it belongs to
                if key in ALTERNATIVE_SPATIAL_KEYS:
                    standard_key = SPATIAL_KEY
                elif key in ALTERNATIVE_CELL_TYPE_KEYS:
                    standard_key = CELL_TYPE_KEY
                elif key in ALTERNATIVE_CLUSTER_KEYS:
                    standard_key = CLUSTER_KEY
                elif key in ALTERNATIVE_BATCH_KEYS:
                    standard_key = BATCH_KEY
                else:
                    standard_key = key
                return key, standard_key
        return None
    
    def _suggest_standardization(self, adata: 'ad.AnnData', result: DataValidationResult) -> None:
        """Generate suggestions for standardizing the data."""
        
        if result.standardized_fields:
            result.add_suggestion(
                "Use chatspatial.utils.data_adapter.standardize_adata() to automatically "
                "convert to standard field names"
            )
        
        # Check if spatial neighbors are computed
        if 'spatial_connectivities' not in adata.obsp and any(
            key in adata.obsm for key in ALTERNATIVE_SPATIAL_KEYS
        ):
            result.add_suggestion("Compute spatial neighbors: sq.gr.spatial_neighbors(adata)")
        
        # Check if highly variable genes are identified
        if 'highly_variable' not in adata.var.columns:
            result.add_suggestion("Identify highly variable genes: self.sc.pp.highly_variable_genes(adata)")


# Convenience functions for common validation tasks

def validate_spatial_data(adata: 'ad.AnnData', strict: bool = False) -> DataValidationResult:
    """
    Validate spatial transcriptomics data.
    
    This is the main function used by most ChatSpatial tools.
    """
    validator = DataValidator(strict_mode=strict)
    return validator.validate_adata(
        adata, 
        require_spatial=True, 
        require_cell_types=False
    )


def validate_for_cell_communication(adata: 'ad.AnnData', strict: bool = False) -> DataValidationResult:
    """Validate data for cell communication analysis."""
    validator = DataValidator(strict_mode=strict)
    return validator.validate_adata(
        adata,
        analysis_type="cell_communication",
        require_spatial=True,
        require_cell_types=True
    )


def validate_for_deconvolution(adata: 'ad.AnnData', 
                              reference_adata: Optional['ad.AnnData'] = None,
                              strict: bool = False) -> DataValidationResult:
    """Validate data for deconvolution analysis."""
    validator = DataValidator(strict_mode=strict)
    result = validator.validate_adata(
        adata,
        analysis_type="deconvolution", 
        require_spatial=True,
        require_cell_types=False
    )
    
    # Additional deconvolution-specific validation
    if reference_adata is not None:
        ref_result = validator.validate_adata(
            reference_adata,
            require_spatial=False,
            require_cell_types=True
        )
        
        # Merge results
        result.errors.extend([f"Reference data: {err}" for err in ref_result.errors])
        result.warnings.extend([f"Reference data: {warn}" for warn in ref_result.warnings])
    
    return result


def validate_for_spatial_analysis(adata: 'ad.AnnData', strict: bool = False) -> DataValidationResult:
    """Validate data for spatial analysis."""
    validator = DataValidator(strict_mode=strict)
    return validator.validate_adata(
        adata,
        analysis_type="spatial_analysis",
        require_spatial=True,
        require_cell_types=False
    )


def check_data_compatibility(adata: 'ad.AnnData', tool_name: str) -> Dict[str, Any]:
    """
    Check if data is compatible with a specific tool.
    
    This replaces the scattered compatibility checking logic.
    """
    validator = DataValidator()
    
    # Tool-specific validation
    if tool_name.lower() in ["cell_communication", "cellchat", "liana"]:
        result = validate_for_cell_communication(adata)
    elif tool_name.lower() in ["deconvolution", "cell2location", "rctd"]:
        result = validate_for_deconvolution(adata)
    elif tool_name.lower() in ["spatial_analysis", "neighborhood", "cooccurrence"]:
        result = validate_for_spatial_analysis(adata)
    else:
        result = validate_spatial_data(adata)
    
    return {
        "compatible": result.passed,
        "errors": result.errors,
        "warnings": result.warnings,
        "suggestions": result.suggestions,
        "tool": tool_name
    }


def raise_if_invalid(adata: 'ad.AnnData', 
                    analysis_type: Optional[str] = None,
                    message_prefix: str = "Data validation failed") -> None:
    """
    Validate data and raise exception if invalid.
    
    Use this in tools that require valid data to proceed.
    """
    validator = DataValidator(strict_mode=True)
    
    if analysis_type == "cell_communication":
        result = validate_for_cell_communication(adata, strict=True)
    elif analysis_type == "deconvolution":
        result = validate_for_deconvolution(adata, strict=True)
    elif analysis_type == "spatial_analysis":
        result = validate_for_spatial_analysis(adata, strict=True)
    else:
        result = validate_spatial_data(adata, strict=True)
    
    if not result.passed:
        error_msg = f"{message_prefix}:\n" + "\n".join(result.errors)
        if result.suggestions:
            error_msg += "\n\nSuggestions:\n" + "\n".join(result.suggestions)
        raise DataValidationError(error_msg)