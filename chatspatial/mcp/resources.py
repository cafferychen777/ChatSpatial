"""
MCP Resource Management

This module handles resource discovery and management for spatial datasets,
analysis results, and visualizations.
"""

from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, field
import json
from datetime import datetime


@dataclass
class Resource:
    """MCP Resource representation"""
    uri: str
    name: str
    mimeType: str
    description: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to MCP resource format"""
        result = {
            "uri": self.uri,
            "name": self.name,
            "mimeType": self.mimeType
        }
        if self.description:
            result["description"] = self.description
        if self.metadata:
            result["metadata"] = self.metadata
        return result


class ResourceManager:
    """Manages MCP resources for spatial data"""
    
    def __init__(self, data_store: Dict[str, Any]):
        self.data_store = data_store
    
    def get_resource_list(self) -> List[Resource]:
        """Get list of available resources based on loaded datasets
        
        Returns:
            List of Resource objects representing available data
        """
        resources = []
        
        # Add resources for each loaded dataset
        for data_id, dataset_info in self.data_store.items():
            # Main dataset resource
            resources.append(Resource(
                uri=f"spatial://datasets/{data_id}",
                name=dataset_info.get("name", f"Dataset {data_id}"),
                mimeType="application/x-anndata",
                description=f"Spatial transcriptomics dataset with {dataset_info.get('n_cells', 0)} cells and {dataset_info.get('n_genes', 0)} genes",
                metadata={
                    "n_cells": dataset_info.get("n_cells", 0),
                    "n_genes": dataset_info.get("n_genes", 0),
                    "data_type": dataset_info.get("type", "unknown"),
                    "has_spatial": dataset_info.get("has_spatial", False)
                }
            ))
            
            # Add resources for analysis results if available
            if "adata" in dataset_info and hasattr(dataset_info["adata"], "uns"):
                uns_keys = dataset_info["adata"].uns.keys()
                
                # Spatial domains
                if any("domain" in key for key in uns_keys):
                    resources.append(Resource(
                        uri=f"spatial://results/{data_id}/domains",
                        name=f"Spatial domains for {dataset_info.get('name', data_id)}",
                        mimeType="application/json",
                        description="Identified spatial domains and clustering results"
                    ))
                
                # Differential expression results
                if any("rank_genes" in key for key in uns_keys):
                    resources.append(Resource(
                        uri=f"spatial://results/{data_id}/markers",
                        name=f"Marker genes for {dataset_info.get('name', data_id)}",
                        mimeType="application/json",
                        description="Differential expression and marker gene results"
                    ))
                
                # Cell communication results
                if any("communication" in key or "cellchat" in key for key in uns_keys):
                    resources.append(Resource(
                        uri=f"spatial://results/{data_id}/communication",
                        name=f"Cell communication for {dataset_info.get('name', data_id)}",
                        mimeType="application/json",
                        description="Cell-cell communication analysis results"
                    ))
        
        # Add general resources
        resources.extend([
            Resource(
                uri="spatial://plots/current",
                name="Current visualization",
                mimeType="image/png",
                description="Most recently generated visualization"
            ),
            Resource(
                uri="spatial://logs/session",
                name="Analysis session log",
                mimeType="text/plain",
                description="Log of current analysis session"
            )
        ])
        
        return resources
    
    def read_resource_content(self, uri: str) -> str:
        """Read content of a resource
        
        Args:
            uri: Resource URI (e.g., "spatial://datasets/data1")
            
        Returns:
            JSON string with resource content
            
        Raises:
            ValueError: If resource not found
        """
        parts = uri.replace("spatial://", "").split("/")
        
        if len(parts) < 2:
            raise ValueError(f"Invalid resource URI: {uri}")
        
        resource_type = parts[0]
        
        if resource_type == "datasets":
            data_id = parts[1]
            if data_id not in self.data_store:
                raise ValueError(f"Dataset {data_id} not found")
            
            dataset_info = self.data_store[data_id]
            summary = {
                "id": data_id,
                "name": dataset_info.get("name", "Unknown"),
                "type": dataset_info.get("type", "Unknown"),
                "n_cells": dataset_info.get("n_cells", 0),
                "n_genes": dataset_info.get("n_genes", 0),
                "has_spatial": dataset_info.get("has_spatial", False),
                "preprocessing_done": dataset_info.get("preprocessing_done", False),
                "available_keys": {
                    "obs": list(dataset_info.get("obs_columns", [])),
                    "var": list(dataset_info.get("var_columns", [])),
                    "obsm": list(dataset_info.get("obsm_keys", [])),
                    "uns": list(dataset_info.get("uns_keys", []))
                }
            }
            return json.dumps(summary, indent=2)
        
        elif resource_type == "results":
            data_id = parts[1]
            result_type = parts[2] if len(parts) > 2 else "all"
            
            if data_id not in self.data_store:
                raise ValueError(f"Dataset {data_id} not found")
            
            dataset_info = self.data_store[data_id]
            if "adata" not in dataset_info:
                return json.dumps({"error": "No analysis results available"})
            
            adata = dataset_info["adata"]
            results = {}
            
            if result_type == "domains" or result_type == "all":
                domain_keys = [key for key in adata.uns.keys() if "domain" in key]
                if domain_keys:
                    results["domains"] = {key: str(adata.uns[key]) for key in domain_keys[:3]}  # Limit output
            
            if result_type == "markers" or result_type == "all":
                marker_keys = [key for key in adata.uns.keys() if "rank_genes" in key]
                if marker_keys:
                    results["markers"] = {key: "Available" for key in marker_keys}
            
            if result_type == "communication" or result_type == "all":
                comm_keys = [key for key in adata.uns.keys() if "communication" in key or "cellchat" in key]
                if comm_keys:
                    results["communication"] = {key: "Available" for key in comm_keys}
            
            return json.dumps(results, indent=2)
        
        elif resource_type == "plots":
            return json.dumps({"message": "Visualization available through visualize_data tool"})
        
        elif resource_type == "logs":
            return json.dumps({"message": "Session logs available"})
        
        else:
            raise ValueError(f"Unknown resource type: {resource_type}")