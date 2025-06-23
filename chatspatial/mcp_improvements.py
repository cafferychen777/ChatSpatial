"""
MCP Improvements for ChatSpatial
Phase 1: Resources, Tool Annotations, and Error Handling
"""

from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, field
from enum import Enum
import json

# Tool Annotations based on MCP specification
TOOL_ANNOTATIONS = {
    "load_data": {
        "title": "Load Spatial Data",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": True  # Reads external files
    },
    "preprocess_data": {
        "title": "Preprocess Spatial Data",
        "readOnlyHint": False,  # Modifies data
        "destructiveHint": False,
        "idempotentHint": False,  # Results may vary with parameters
        "openWorldHint": False
    },
    "visualize_data": {
        "title": "Visualize Spatial Data",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False
    },
    "annotate_cells": {
        "title": "Annotate Cell Types",
        "readOnlyHint": False,  # Adds annotations
        "destructiveHint": False,
        "idempotentHint": False,  # Results may vary
        "openWorldHint": True  # May use external resources
    },
    "analyze_spatial_data": {
        "title": "Spatial Pattern Analysis",
        "readOnlyHint": False,  # Stores results
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False
    },
    "find_markers": {
        "title": "Find Marker Genes",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False
    },
    "analyze_velocity_data": {
        "title": "RNA Velocity Analysis",
        "readOnlyHint": False,
        "destructiveHint": False,
        "idempotentHint": False,
        "openWorldHint": False
    },
    "analyze_trajectory_data": {
        "title": "Trajectory Analysis",
        "readOnlyHint": False,
        "destructiveHint": False,
        "idempotentHint": False,
        "openWorldHint": False
    },
    "integrate_samples": {
        "title": "Integrate Multiple Samples",
        "readOnlyHint": False,
        "destructiveHint": False,
        "idempotentHint": False,
        "openWorldHint": False
    },
    "deconvolve_data": {
        "title": "Spatial Deconvolution",
        "readOnlyHint": False,
        "destructiveHint": False,
        "idempotentHint": False,
        "openWorldHint": True  # May use external reference data
    },
    "identify_spatial_domains": {
        "title": "Identify Spatial Domains",
        "readOnlyHint": False,
        "destructiveHint": False,
        "idempotentHint": False,
        "openWorldHint": False
    },
    "analyze_cell_communication": {
        "title": "Cell Communication Analysis",
        "readOnlyHint": False,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": True  # May use external LR databases
    },
    "analyze_enrichment": {
        "title": "Gene Set Enrichment Analysis",
        "readOnlyHint": False,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False
    },
    "find_spatial_genes": {
        "title": "Find Spatial Variable Genes",
        "readOnlyHint": False,
        "destructiveHint": False,
        "idempotentHint": False,  # Neural network training
        "openWorldHint": False
    }
}


@dataclass
class Resource:
    """MCP Resource representation"""
    uri: str
    name: str
    mimeType: str
    description: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Prompt:
    """MCP Prompt representation"""
    name: str
    description: str
    arguments: List[Dict[str, Any]] = field(default_factory=list)


@dataclass
class PromptArgument:
    """MCP Prompt Argument"""
    name: str
    description: str
    required: bool = True


class ErrorType(Enum):
    """MCP Error Types"""
    INVALID_REQUEST = -32600
    METHOD_NOT_FOUND = -32601
    INVALID_PARAMS = -32602
    INTERNAL_ERROR = -32603
    
    # Custom error codes (range -32000 to -32099)
    DATASET_NOT_FOUND = -32000
    INVALID_DATA_FORMAT = -32001
    ANALYSIS_FAILED = -32002
    VISUALIZATION_ERROR = -32003
    REFERENCE_DATA_ERROR = -32004


@dataclass
class MCPError:
    """MCP Error Format"""
    code: int
    message: str
    data: Optional[Dict[str, Any]] = None
    
    def to_dict(self) -> Dict[str, Any]:
        result = {
            "code": self.code,
            "message": self.message
        }
        if self.data:
            result["data"] = self.data
        return result


# Prompt Templates for Spatial Analysis
SPATIAL_PROMPTS = [
    Prompt(
        name="analyze-spatial-expression",
        description="Analyze spatial gene expression patterns",
        arguments=[
            {"name": "genes", "description": "Genes to analyze", "required": True},
            {"name": "method", "description": "Analysis method", "required": False}
        ]
    ),
    Prompt(
        name="find-cell-types",
        description="Identify cell types in spatial data",
        arguments=[
            {"name": "method", "description": "Cell type identification method", "required": False},
            {"name": "reference_data", "description": "Reference dataset path", "required": False}
        ]
    ),
    Prompt(
        name="compare-conditions",
        description="Compare spatial patterns between conditions",
        arguments=[
            {"name": "condition_key", "description": "Column defining conditions", "required": True},
            {"name": "condition1", "description": "First condition", "required": True},
            {"name": "condition2", "description": "Second condition", "required": True}
        ]
    ),
    Prompt(
        name="generate-visualization",
        description="Generate spatial visualization",
        arguments=[
            {"name": "plot_type", "description": "Type of visualization", "required": True},
            {"name": "features", "description": "Features to visualize", "required": False},
            {"name": "save_path", "description": "Path to save figure", "required": False}
        ]
    ),
    Prompt(
        name="quality-control",
        description="Perform quality control on spatial data",
        arguments=[
            {"name": "metrics", "description": "QC metrics to compute", "required": False},
            {"name": "thresholds", "description": "QC thresholds", "required": False}
        ]
    ),
    Prompt(
        name="batch-correction",
        description="Correct batch effects in integrated data",
        arguments=[
            {"name": "batch_key", "description": "Column defining batches", "required": True},
            {"name": "method", "description": "Batch correction method", "required": False}
        ]
    ),
    Prompt(
        name="spatial-clustering",
        description="Perform spatial clustering analysis",
        arguments=[
            {"name": "n_clusters", "description": "Number of clusters", "required": False},
            {"name": "resolution", "description": "Clustering resolution", "required": False}
        ]
    ),
    Prompt(
        name="trajectory-inference",
        description="Infer cellular trajectories",
        arguments=[
            {"name": "start_cell", "description": "Starting cell type", "required": False},
            {"name": "end_cell", "description": "Ending cell type", "required": False}
        ]
    )
]


def format_mcp_error(error_type: ErrorType, message: str, data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Format error according to MCP specification"""
    return MCPError(
        code=error_type.value,
        message=message,
        data=data
    ).to_dict()


def get_resource_list(data_store: Dict[str, Any]) -> List[Resource]:
    """Generate list of available resources based on loaded data"""
    resources = []
    
    # Add resources for each loaded dataset
    for data_id, dataset_info in data_store.items():
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


def read_resource_content(uri: str, data_store: Dict[str, Any]) -> str:
    """Read content of a resource"""
    parts = uri.replace("spatial://", "").split("/")
    
    if len(parts) < 2:
        raise ValueError(f"Invalid resource URI: {uri}")
    
    resource_type = parts[0]
    
    if resource_type == "datasets":
        data_id = parts[1]
        if data_id not in data_store:
            raise ValueError(f"Dataset {data_id} not found")
        
        dataset_info = data_store[data_id]
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
        
        if data_id not in data_store:
            raise ValueError(f"Dataset {data_id} not found")
        
        dataset_info = data_store[data_id]
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


def get_tool_with_annotations(tool_name: str) -> Dict[str, Any]:
    """Get tool information including annotations"""
    if tool_name not in TOOL_ANNOTATIONS:
        return {"name": tool_name}
    
    return {
        "name": tool_name,
        **TOOL_ANNOTATIONS[tool_name]
    }


def get_all_tools_with_annotations() -> List[Dict[str, Any]]:
    """Get all tools with their annotations"""
    return [
        {"name": name, **annotations}
        for name, annotations in TOOL_ANNOTATIONS.items()
    ]


def prompt_to_tool_params(prompt_name: str, arguments: Dict[str, Any], data_store: Dict[str, Any]) -> Dict[str, Any]:
    """Convert prompt arguments to tool parameters"""
    
    # Get the first available dataset if not specified
    if not data_store:
        raise ValueError("No datasets loaded")
    
    data_id = list(data_store.keys())[0]
    
    if prompt_name == "analyze-spatial-expression":
        return {
            "tool": "analyze_spatial_data",
            "params": {
                "data_id": data_id,
                "params": {
                    "analysis_type": "spatial_autocorrelation",
                    "genes": arguments.get("genes", []),
                    "method": arguments.get("method", "moran")
                }
            }
        }
    
    elif prompt_name == "find-cell-types":
        return {
            "tool": "annotate_cells",
            "params": {
                "data_id": data_id,
                "params": {
                    "method": arguments.get("method", "marker_genes"),
                    "reference_data_id": arguments.get("reference_data")
                }
            }
        }
    
    elif prompt_name == "compare-conditions":
        return {
            "tool": "find_markers",
            "params": {
                "data_id": data_id,
                "group_key": arguments["condition_key"],
                "group1": arguments["condition1"],
                "group2": arguments["condition2"]
            }
        }
    
    elif prompt_name == "generate-visualization":
        vis_params = {
            "plot_type": arguments["plot_type"]
        }
        if "features" in arguments:
            vis_params["feature"] = arguments["features"][0] if isinstance(arguments["features"], list) else arguments["features"]
        
        return {
            "tool": "visualize_data",
            "params": {
                "data_id": data_id,
                "params": vis_params
            }
        }
    
    elif prompt_name == "quality-control":
        return {
            "tool": "preprocess_data",
            "params": {
                "data_id": data_id,
                "params": {
                    "filter_genes": True,
                    "filter_cells": True,
                    "normalize": False  # Just QC, no normalization
                }
            }
        }
    
    elif prompt_name == "batch-correction":
        return {
            "tool": "integrate_samples",
            "params": {
                "data_ids": [data_id],  # Would need multiple IDs in practice
                "params": {
                    "batch_key": arguments["batch_key"],
                    "method": arguments.get("method", "harmony")
                }
            }
        }
    
    elif prompt_name == "spatial-clustering":
        return {
            "tool": "identify_spatial_domains",
            "params": {
                "data_id": data_id,
                "params": {
                    "method": "stlearn",
                    "n_clusters": arguments.get("n_clusters"),
                    "resolution": arguments.get("resolution", 1.0)
                }
            }
        }
    
    elif prompt_name == "trajectory-inference":
        return {
            "tool": "analyze_trajectory_data",
            "params": {
                "data_id": data_id,
                "params": {
                    "method": "cellrank",
                    "start_cell_type": arguments.get("start_cell"),
                    "end_cell_type": arguments.get("end_cell")
                }
            }
        }
    
    else:
        raise ValueError(f"Unknown prompt: {prompt_name}")