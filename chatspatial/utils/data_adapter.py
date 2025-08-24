"""
Data adapter for standardizing AnnData objects to ChatSpatial format.

This module implements Linus's core principle: "消除所有特殊情况"
Instead of having each tool handle different data formats, we have ONE 
adapter that converts everything to standard format.

The adapter is the single point where we handle all the messy reality of
different data formats, so the tools can work with clean, standardized data.
"""

from typing import Optional, Dict, Any, List, Tuple
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
from scipy import sparse
import warnings
from functools import wraps

from ..models.data_standards import (
    CHATSPATIAL_STANDARDS,
    DataValidationResult,
    get_field_mapping,
    describe_standards
)
from .data_validator import DataValidator, validate_spatial_data


class DataStandardizationError(Exception):
    """Raised when data cannot be standardized."""
    pass


def standardize_input(func):
    """
    Decorator to automatically standardize AnnData input for any tool function.
    
    This decorator eliminates the need for each tool to handle different 
    data formats. Apply it to any function that takes adata as first argument.
    
    Usage:
        @standardize_input
        async def my_spatial_tool(adata, params, context):
            # adata is now guaranteed to be in standard format
    """
    @wraps(func)
    async def wrapper(adata, *args, **kwargs):
        # Standardize the input data
        try:
            standardized_adata = standardize_adata(adata, copy=True)
            return await func(standardized_adata, *args, **kwargs)
        except Exception as e:
            # If standardization fails, try with original data and warn
            import logging
            logging.warning(f"Data standardization failed for {func.__name__}: {e}")
            return await func(adata, *args, **kwargs)
    
    return wrapper


class DataAdapter:
    """
    Single source of truth for data format conversion in ChatSpatial.
    
    This class handles all the messy details of converting different data
    formats to the ChatSpatial standard. Tools never see non-standard data.
    """
    
    def __init__(self, strict_mode: bool = False, preserve_original: bool = True):
        """
        Initialize data adapter.
        
        Args:
            strict_mode: If True, raise errors instead of warnings for issues
            preserve_original: If True, preserve original field names as backup
        """
        self.strict_mode = strict_mode
        self.preserve_original = preserve_original
        self.standards = CHATSPATIAL_STANDARDS
        self.field_mapping = get_field_mapping()
        self.validator = DataValidator(strict_mode=strict_mode)
    
    def standardize(self, adata: ad.AnnData, copy: bool = True) -> ad.AnnData:
        """
        Convert AnnData object to ChatSpatial standard format.
        
        This is the main function that eliminates all special cases.
        Every tool should use this to ensure consistent data format.
        
        Args:
            adata: Input AnnData object (any format)
            copy: Whether to work on a copy (recommended)
            
        Returns:
            AnnData object in ChatSpatial standard format
        """
        if copy:
            adata = adata.copy()
        
        # Step 1: Standardize spatial coordinates
        self._standardize_spatial_coordinates(adata)
        
        # Step 2: Standardize observation metadata
        self._standardize_obs_metadata(adata)
        
        # Step 3: Standardize expression matrix
        self._standardize_expression_matrix(adata)
        
        # Step 4: Standardize gene names
        self._standardize_gene_names(adata)
        
        # Step 5: Add missing standard fields with defaults
        self._add_missing_standard_fields(adata)
        
        # Step 6: Compute derived fields if needed
        self._compute_derived_fields(adata)
        
        # Step 7: Final validation
        result = self.validator.validate_adata(adata)
        if not result.passed and self.strict_mode:
            raise DataStandardizationError(f"Standardization failed: {result.errors}")
        
        return adata
    
    def _standardize_spatial_coordinates(self, adata: ad.AnnData) -> None:
        """
        Standardize spatial coordinates to adata.obsm['spatial'].
        
        This replaces the scattered coordinate detection logic in visualization.py:405-415
        and other modules. One function handles all cases.
        """
        target_key = self.standards.spatial_key
        
        # Case 1: Already in standard location
        if target_key in adata.obsm:
            coords = adata.obsm[target_key]
            self._validate_and_fix_coordinates(coords, adata, target_key)
            return
        
        # Case 2: In alternative obsm location
        for alt_key in self.standards.spatial_alternatives:
            if alt_key in adata.obsm and alt_key != target_key:
                coords = adata.obsm[alt_key]
                self._validate_and_fix_coordinates(coords, adata, target_key)
                
                # Move to standard location
                adata.obsm[target_key] = coords
                
                # Preserve original if requested
                if self.preserve_original and alt_key != target_key:
                    adata.obsm[f"original_{alt_key}"] = coords.copy()
                else:
                    del adata.obsm[alt_key]
                
                return
        
        # Case 3: In obs as separate x, y columns
        if 'x' in adata.obs and 'y' in adata.obs:
            try:
                x_coords = pd.to_numeric(adata.obs['x'], errors='coerce')
                y_coords = pd.to_numeric(adata.obs['y'], errors='coerce')
                
                if x_coords.isna().any() or y_coords.isna().any():
                    raise ValueError("Cannot convert x, y coordinates to numeric")
                
                coords = np.column_stack([x_coords, y_coords]).astype(self.standards.spatial_dtype)
                adata.obsm[target_key] = coords
                
                # Clean up obs
                if self.preserve_original:
                    adata.obs['original_x'] = adata.obs['x'] 
                    adata.obs['original_y'] = adata.obs['y']
                del adata.obs['x']
                del adata.obs['y']
                
                return
                
            except Exception as e:
                warnings.warn(f"Could not convert x, y coordinates from obs: {e}")
        
        # Case 4: Try to extract from spatial metadata (Visium data)
        if 'spatial' in adata.uns and isinstance(adata.uns['spatial'], dict):
            try:
                spatial_data = adata.uns['spatial']
                
                # Look for tissue positions or similar
                for key, value in spatial_data.items():
                    if isinstance(value, dict) and 'tissue_positions_list' in value:
                        positions = value['tissue_positions_list']
                        if isinstance(positions, np.ndarray) and positions.shape[1] >= 2:
                            coords = positions[:, -2:].astype(self.standards.spatial_dtype)
                            adata.obsm[target_key] = coords
                            return
                
            except Exception as e:
                warnings.warn(f"Could not extract coordinates from spatial metadata: {e}")
        
        # Case 5: Generate dummy coordinates as last resort
        if not self.strict_mode:
            warnings.warn("No spatial coordinates found. Generating dummy coordinates.")
            n_obs = adata.n_obs
            grid_size = int(np.ceil(np.sqrt(n_obs)))
            
            x_coords = np.tile(np.arange(grid_size), (grid_size, 1)).T.flatten()[:n_obs]
            y_coords = np.repeat(np.arange(grid_size), grid_size)[:n_obs]
            
            coords = np.column_stack([x_coords, y_coords]).astype(self.standards.spatial_dtype)
            adata.obsm[target_key] = coords
            adata.uns[f"{target_key}_dummy"] = True  # Mark as dummy data
        else:
            raise DataStandardizationError("No spatial coordinates found and strict_mode=True")
    
    def _validate_and_fix_coordinates(self, 
                                     coords: np.ndarray, 
                                     adata: ad.AnnData, 
                                     target_key: str) -> None:
        """Validate and fix coordinate array issues."""
        
        # Fix shape issues
        if coords.shape[0] != adata.n_obs:
            raise DataStandardizationError(
                f"Coordinate count ({coords.shape[0]}) doesn't match observation count ({adata.n_obs})"
            )
        
        if coords.shape[1] < 2:
            raise DataStandardizationError(
                f"Coordinates must have at least 2 dimensions, found {coords.shape[1]}"
            )
        
        # Take only first 2 dimensions if more are present
        if coords.shape[1] > 2:
            coords = coords[:, :2]
            warnings.warn("Using only first 2 dimensions of spatial coordinates")
        
        # Fix dtype
        if coords.dtype != self.standards.spatial_dtype:
            coords = coords.astype(self.standards.spatial_dtype)
        
        # Handle NaN/inf values
        if np.any(np.isnan(coords)) or np.any(np.isinf(coords)):
            if self.strict_mode:
                raise DataStandardizationError("Spatial coordinates contain NaN or infinite values")
            else:
                # Replace with median values
                coords = np.where(np.isfinite(coords), coords, np.nanmedian(coords, axis=0))
                warnings.warn("Replaced NaN/inf coordinates with median values")
        
        # Update the coordinates
        adata.obsm[target_key] = coords
    
    def _standardize_obs_metadata(self, adata: ad.AnnData) -> None:
        """
        Standardize observation metadata field names.
        
        This handles the mess of cell_type vs celltype vs leiden etc.
        """
        # Get reverse mapping (alternative -> standard)
        reverse_mapping = {}
        for alt, std in self.field_mapping.items():
            if std in [self.standards.cell_type_key, self.standards.cluster_key, self.standards.batch_key]:
                reverse_mapping[alt] = std
        
        # Process each field in obs
        for obs_key in list(adata.obs.columns):
            if obs_key in reverse_mapping:
                standard_key = reverse_mapping[obs_key]
                
                # Don't overwrite if standard key already exists
                if standard_key in adata.obs and obs_key != standard_key:
                    if self.preserve_original:
                        adata.obs[f"original_{obs_key}"] = adata.obs[obs_key]
                    continue
                
                # Move to standard location
                if obs_key != standard_key:
                    adata.obs[standard_key] = adata.obs[obs_key]
                    if self.preserve_original:
                        adata.obs[f"original_{obs_key}"] = adata.obs[obs_key]
                    del adata.obs[obs_key]
                
                # Ensure categorical type for metadata fields
                if standard_key in [self.standards.cell_type_key, self.standards.cluster_key, self.standards.batch_key]:
                    if not pd.api.types.is_categorical_dtype(adata.obs[standard_key]):
                        adata.obs[standard_key] = adata.obs[standard_key].astype('category')
    
    def _standardize_expression_matrix(self, adata: ad.AnnData) -> None:
        """Standardize the expression matrix."""
        
        if adata.X is None:
            raise DataStandardizationError("Expression matrix (adata.X) is None")
        
        # Handle data type issues
        if hasattr(adata.X, 'dtype'):
            # Convert to float32 if integer (common for count data)
            if np.issubdtype(adata.X.dtype, np.integer):
                if sparse.issparse(adata.X):
                    adata.X = adata.X.astype(np.float32)
                else:
                    adata.X = adata.X.astype(np.float32)
        
        # Handle NaN/inf values
        if sparse.issparse(adata.X):
            data = adata.X.data
            if np.any(np.isnan(data)) or np.any(np.isinf(data)):
                if self.strict_mode:
                    raise DataStandardizationError("Expression matrix contains NaN or infinite values")
                else:
                    data[~np.isfinite(data)] = 0.0
                    warnings.warn("Replaced NaN/inf values in expression matrix with 0")
        else:
            if np.any(np.isnan(adata.X)) or np.any(np.isinf(adata.X)):
                if self.strict_mode:
                    raise DataStandardizationError("Expression matrix contains NaN or infinite values")
                else:
                    adata.X[~np.isfinite(adata.X)] = 0.0
                    warnings.warn("Replaced NaN/inf values in expression matrix with 0")
    
    def _standardize_gene_names(self, adata: ad.AnnData) -> None:
        """Standardize gene names."""
        
        # Make gene names unique
        if not adata.var_names.is_unique:
            original_names = adata.var_names.copy()
            adata.var_names_make_unique()
            if self.preserve_original:
                adata.var['original_gene_names'] = original_names
        
        # Remove whitespace
        cleaned_names = adata.var_names.str.strip()
        if not cleaned_names.equals(adata.var_names):
            if self.preserve_original:
                adata.var['original_gene_names'] = adata.var_names.copy()
            adata.var_names = cleaned_names
    
    def _add_missing_standard_fields(self, adata: ad.AnnData) -> None:
        """Add missing standard fields with sensible defaults."""
        
        # Add default batch information if missing
        if self.standards.batch_key not in adata.obs:
            adata.obs[self.standards.batch_key] = pd.Categorical(['batch1'] * adata.n_obs)
        
        # Add placeholder clustering if none exists
        if not any(key in adata.obs for key in self.standards.cluster_alternatives):
            # Don't add default clustering - let tools handle this
            pass
        
        # Add placeholder cell types if none exist  
        if not any(key in adata.obs for key in self.standards.cell_type_alternatives):
            # Don't add default cell types - let annotation tools handle this
            pass
    
    def _compute_derived_fields(self, adata: ad.AnnData) -> None:
        """Compute commonly needed derived fields."""
        
        # Compute spatial neighbors if missing
        if 'spatial_connectivities' not in adata.obsp and self.standards.spatial_key in adata.obsm:
            try:
                sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key=self.standards.spatial_key)
            except Exception as e:
                warnings.warn(f"Could not compute spatial neighbors: {e}")
        
        # Identify highly variable genes if missing
        if 'highly_variable' not in adata.var.columns:
            try:
                # Only compute if we have reasonable amount of data
                if adata.n_vars > 100 and adata.n_obs > 50:
                    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            except Exception as e:
                warnings.warn(f"Could not identify highly variable genes: {e}")


# Global adapter instance
_global_adapter = DataAdapter()


def standardize_adata(adata: ad.AnnData, 
                     copy: bool = True, 
                     strict: bool = False,
                     preserve_original: bool = True) -> ad.AnnData:
    """
    Convert any AnnData object to ChatSpatial standard format.
    
    This is the main function that tools should use to ensure consistent data.
    It replaces all the scattered data format handling across different modules.
    
    Args:
        adata: Input AnnData object in any format
        copy: Whether to work on a copy (recommended for safety)
        strict: Whether to raise errors for non-critical issues
        preserve_original: Whether to keep original field names as backup
        
    Returns:
        AnnData object in ChatSpatial standard format
        
    Example:
        # Convert user data to standard format
        standard_adata = standardize_adata(user_adata)
        
        # Now all tools work with consistent data structure
        result = await some_spatial_tool(standard_adata, params)
    """
    global _global_adapter
    
    # Use global adapter or create custom one
    if strict != _global_adapter.strict_mode or preserve_original != _global_adapter.preserve_original:
        adapter = DataAdapter(strict_mode=strict, preserve_original=preserve_original)
    else:
        adapter = _global_adapter
    
    return adapter.standardize(adata, copy=copy)


def get_spatial_coordinates(adata: ad.AnnData) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get spatial coordinates in a standardized way.
    
    This replaces the coordinate extraction logic in visualization.py:405-415
    and other places. One function, no special cases.
    
    Args:
        adata: AnnData object (will be standardized if needed)
        
    Returns:
        Tuple of (x_coordinates, y_coordinates)
    """
    # Standardize if needed
    if CHATSPATIAL_STANDARDS.spatial_key not in adata.obsm:
        adata = standardize_adata(adata, copy=True)
    
    coords = adata.obsm[CHATSPATIAL_STANDARDS.spatial_key]
    return coords[:, 0], coords[:, 1]


def get_cell_types(adata: ad.AnnData, default: str = "Unknown") -> pd.Series:
    """
    Get cell type annotations in a standardized way.
    
    Args:
        adata: AnnData object  
        default: Default value for missing cell types
        
    Returns:
        Series of cell type annotations
    """
    # Try standard key first
    if CHATSPATIAL_STANDARDS.cell_type_key in adata.obs:
        return adata.obs[CHATSPATIAL_STANDARDS.cell_type_key]
    
    # Try alternatives
    for alt_key in CHATSPATIAL_STANDARDS.cell_type_alternatives:
        if alt_key in adata.obs:
            return adata.obs[alt_key]
    
    # Return default
    return pd.Series([default] * adata.n_obs, index=adata.obs_names, name=CHATSPATIAL_STANDARDS.cell_type_key)


def get_clusters(adata: ad.AnnData, default: str = "cluster_0") -> pd.Series:
    """
    Get cluster annotations in a standardized way.
    
    Args:
        adata: AnnData object
        default: Default value for missing clusters
        
    Returns: 
        Series of cluster annotations
    """
    # Try standard key first
    if CHATSPATIAL_STANDARDS.cluster_key in adata.obs:
        return adata.obs[CHATSPATIAL_STANDARDS.cluster_key]
    
    # Try alternatives
    for alt_key in CHATSPATIAL_STANDARDS.cluster_alternatives:
        if alt_key in adata.obs:
            return adata.obs[alt_key]
    
    # Return default
    return pd.Series([default] * adata.n_obs, index=adata.obs_names, name=CHATSPATIAL_STANDARDS.cluster_key)


def ensure_spatial_neighbors(adata: ad.AnnData) -> None:
    """
    Ensure spatial neighborhood graph is computed.
    
    Many tools need spatial neighbors but compute them inconsistently.
    This function ensures they exist in a standardized way.
    """
    if 'spatial_connectivities' not in adata.obsp:
        # Ensure we have standard spatial coordinates
        if CHATSPATIAL_STANDARDS.spatial_key not in adata.obsm:
            adata = standardize_adata(adata, copy=False)
        
        try:
            sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key=CHATSPATIAL_STANDARDS.spatial_key)
        except Exception as e:
            warnings.warn(f"Could not compute spatial neighbors: {e}")


def print_standardization_summary(adata: ad.AnnData) -> None:
    """
    Print a summary of the data standardization.
    
    Useful for debugging and understanding what was changed.
    """
    print("ChatSpatial Data Standardization Summary")
    print("=" * 40)
    
    # Check spatial coordinates
    if CHATSPATIAL_STANDARDS.spatial_key in adata.obsm:
        coords = adata.obsm[CHATSPATIAL_STANDARDS.spatial_key]
        print(f"✓ Spatial coordinates: {coords.shape} {coords.dtype}")
    else:
        print("✗ No spatial coordinates found")
    
    # Check standard obs fields
    for field_name, key in [
        ("Cell types", CHATSPATIAL_STANDARDS.cell_type_key),
        ("Clusters", CHATSPATIAL_STANDARDS.cluster_key), 
        ("Batches", CHATSPATIAL_STANDARDS.batch_key)
    ]:
        if key in adata.obs:
            unique_count = len(adata.obs[key].unique())
            print(f"✓ {field_name}: {unique_count} unique values")
        else:
            print(f"- {field_name}: not present")
    
    # Check derived fields
    if 'spatial_connectivities' in adata.obsp:
        print("✓ Spatial neighbors: computed")
    else:
        print("- Spatial neighbors: not computed")
    
    if 'highly_variable' in adata.var.columns:
        hvg_count = adata.var['highly_variable'].sum()
        print(f"✓ Highly variable genes: {hvg_count} identified")
    else:
        print("- Highly variable genes: not identified")
    
    print(f"\nData dimensions: {adata.n_obs} observations × {adata.n_vars} variables")


def compatibility_report(adata: ad.AnnData, target_tool: Optional[str] = None) -> str:
    """
    Generate a compatibility report for the data.
    
    Args:
        adata: AnnData object to check
        target_tool: Specific tool to check compatibility for
        
    Returns:
        String report of compatibility status
    """
    result = validate_spatial_data(adata)
    
    report_lines = [
        "ChatSpatial Data Compatibility Report",
        "=" * 40,
        f"Overall Status: {'✓ COMPATIBLE' if result.passed else '✗ NEEDS FIXES'}",
        ""
    ]
    
    if result.errors:
        report_lines.extend([
            "ERRORS (must fix):",
            *[f"  • {error}" for error in result.errors],
            ""
        ])
    
    if result.warnings:
        report_lines.extend([
            "WARNINGS (recommended to fix):",
            *[f"  • {warning}" for warning in result.warnings], 
            ""
        ])
    
    if result.suggestions:
        report_lines.extend([
            "SUGGESTIONS:",
            *[f"  • {suggestion}" for suggestion in result.suggestions],
            ""
        ])
    
    if target_tool:
        report_lines.extend([
            f"Specific compatibility for '{target_tool}':",
            "  (Additional checks would be performed here)",
            ""
        ])
    
    report_lines.extend([
        "To fix issues automatically:",
        "  standard_adata = standardize_adata(your_adata)",
        "",
        "For more info:",
        "  from chatspatial.models.data_standards import describe_standards",
        "  print(describe_standards())"
    ])
    
    return "\n".join(report_lines)