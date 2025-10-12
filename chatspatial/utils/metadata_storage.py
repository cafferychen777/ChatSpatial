"""
Metadata storage utilities for scientific provenance tracking.

This module provides utilities to store analysis metadata in AnnData objects,
focusing on scientifically important information for reproducibility and quality control.
"""

from typing import Any, Dict, List, Optional

from anndata import AnnData


def store_analysis_metadata(
    adata: AnnData,
    analysis_name: str,
    method: str,
    parameters: Dict[str, Any],
    results_keys: Dict[str, List[str]],
    statistics: Optional[Dict[str, Any]] = None,
    species: Optional[str] = None,
    database: Optional[str] = None,
    reference_info: Optional[Dict[str, Any]] = None,
) -> None:
    """Store analysis metadata in adata.uns for scientific provenance tracking.

    This function stores ONLY scientifically important metadata:
    - Method name (required for reproducibility)
    - Parameters (required for reproducibility)
    - Results locations (required for data access)
    - Statistics (required for quality assessment)
    - Species/Database (required for biological interpretation)
    - Reference info (required for reference-based methods)

    Non-scientific metadata (timestamp, computation time, user, package versions) are NOT stored.

    Args:
        adata: AnnData object to store metadata in
        analysis_name: Name of the analysis (e.g., "annotation_tangram")
        method: Method name (e.g., "tangram", "liana", "cellrank")
        parameters: Dictionary of analysis parameters
        results_keys: Dictionary mapping storage location to list of keys
            Example: {"obs": ["cell_type_tangram"], "obsm": ["tangram_ct_pred"]}
        statistics: Optional dictionary of quality/summary statistics
        species: Optional species identifier (critical for communication/enrichment)
        database: Optional database/resource name (critical for communication/enrichment)
        reference_info: Optional reference dataset information

    Examples:
        >>> # Annotation metadata
        >>> store_analysis_metadata(
        ...     adata,
        ...     analysis_name="annotation_tangram",
        ...     method="tangram",
        ...     parameters={"device": "cpu", "n_epochs": 1000},
        ...     results_keys={"obs": ["cell_type_tangram"], "obsm": ["tangram_ct_pred"]},
        ...     statistics={"mapping_score": 0.85, "n_cell_types": 10}
        ... )

        >>> # Cell communication metadata
        >>> store_analysis_metadata(
        ...     adata,
        ...     analysis_name="cell_communication_liana",
        ...     method="liana",
        ...     parameters={"n_perms": 100, "nz_prop": 0.2},
        ...     results_keys={"uns": ["liana_res"], "obsm": ["liana_spatial_scores"]},
        ...     statistics={"n_lr_pairs": 500, "n_significant": 50},
        ...     species="mouse",
        ...     database="mouseconsensus"
        ... )

        >>> # Result stored in adata.uns will be:
        >>> # {
        >>> #     'method': 'liana',
        >>> #     'parameters': {'n_perms': 100, 'nz_prop': 0.2},
        >>> #     'results_keys': {'uns': ['liana_res'], 'obsm': ['liana_spatial_scores']},
        >>> #     'statistics': {'n_lr_pairs': 500, 'n_significant': 50},
        >>> #     'species': 'mouse',
        >>> #     'database': 'mouseconsensus'
        >>> # }
    """
    # Build metadata dictionary - only scientifically important information
    metadata = {
        "method": method,
        "parameters": parameters,
        "results_keys": results_keys,
    }

    # Add optional scientific metadata
    if statistics is not None:
        metadata["statistics"] = statistics

    if species is not None:
        metadata["species"] = species

    if database is not None:
        metadata["database"] = database

    if reference_info is not None:
        metadata["reference_info"] = reference_info

    # Store in adata.uns with unique key
    metadata_key = f"{analysis_name}_metadata"
    adata.uns[metadata_key] = metadata


def get_analysis_metadata(
    adata: AnnData, analysis_name: str
) -> Optional[Dict[str, Any]]:
    """Retrieve analysis metadata from adata.uns.

    Args:
        adata: AnnData object
        analysis_name: Name of the analysis

    Returns:
        Dictionary of metadata, or None if not found
    """
    metadata_key = f"{analysis_name}_metadata"
    return adata.uns.get(metadata_key)


def list_analysis_metadata(adata: AnnData) -> List[str]:
    """List all analysis metadata keys stored in adata.uns.

    Args:
        adata: AnnData object

    Returns:
        List of analysis names that have metadata stored
    """
    metadata_keys = [key for key in adata.uns.keys() if key.endswith("_metadata")]
    # Remove "_metadata" suffix to get analysis names
    return [key[: -len("_metadata")] for key in metadata_keys]
