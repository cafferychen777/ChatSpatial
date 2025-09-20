"""
Tool Annotations for MCP

This module provides metadata annotations for spatial analysis tools,
including hints about their behavior and resource usage.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional


@dataclass
class ToolAnnotation:
    """Tool annotation metadata for MCP"""

    title: str
    readOnlyHint: bool = False
    destructiveHint: bool = False
    idempotentHint: bool = True
    openWorldHint: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Convert to MCP annotation format"""
        return {
            "title": self.title,
            "readOnlyHint": self.readOnlyHint,
            "destructiveHint": self.destructiveHint,
            "idempotentHint": self.idempotentHint,
            "openWorldHint": self.openWorldHint,
        }


# Tool annotations dictionary - updated names to match actual tool names
TOOL_ANNOTATIONS: Dict[str, ToolAnnotation] = {
    "load_spatial_data": ToolAnnotation(
        title="Load Spatial Data",
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=True,  # Reads external files
    ),
    "preprocess_data": ToolAnnotation(
        title="Preprocess Spatial Data",
        readOnlyHint=False,  # Modifies data
        destructiveHint=False,
        idempotentHint=False,  # Results may vary with parameters
        openWorldHint=False,
    ),
    "visualize_data": ToolAnnotation(
        title="Visualize Spatial Data",
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    ),
    "annotate_cell_types": ToolAnnotation(
        title="Annotate Cell Types",
        readOnlyHint=False,  # Adds annotations
        destructiveHint=False,
        idempotentHint=False,  # Results may vary
        openWorldHint=True,  # May use external resources
    ),
    "perform_differential_expression": ToolAnnotation(
        title="Differential Expression Analysis",
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    ),
    "find_spatial_variable_genes": ToolAnnotation(
        title="Find Spatial Variable Genes",
        readOnlyHint=False,
        destructiveHint=False,
        idempotentHint=False,  # Neural network training
        openWorldHint=False,
    ),
    "identify_spatial_domains": ToolAnnotation(
        title="Identify Spatial Domains",
        readOnlyHint=False,
        destructiveHint=False,
        idempotentHint=False,
        openWorldHint=False,
    ),
    "analyze_cell_communication": ToolAnnotation(
        title="Cell Communication Analysis",
        readOnlyHint=False,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=True,  # May use external LR databases
    ),
    "perform_enrichment_analysis": ToolAnnotation(
        title="Gene Set Enrichment Analysis",
        readOnlyHint=False,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=True,  # Uses external databases
    ),
    "analyze_rna_velocity": ToolAnnotation(
        title="RNA Velocity Analysis",
        readOnlyHint=False,
        destructiveHint=False,
        idempotentHint=False,
        openWorldHint=False,
    ),
    "analyze_trajectory": ToolAnnotation(
        title="Trajectory Analysis",
        readOnlyHint=False,
        destructiveHint=False,
        idempotentHint=False,
        openWorldHint=False,
    ),
    "integrate_datasets": ToolAnnotation(
        title="Integrate Multiple Samples",
        readOnlyHint=False,
        destructiveHint=False,
        idempotentHint=False,
        openWorldHint=False,
    ),
    "perform_deconvolution": ToolAnnotation(
        title="Spatial Deconvolution",
        readOnlyHint=False,
        destructiveHint=False,
        idempotentHint=False,
        openWorldHint=True,  # May use external reference data
    ),
    "list_datasets": ToolAnnotation(
        title="List Datasets",
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    ),
}


def get_tool_annotation(tool_name: str) -> Optional[ToolAnnotation]:
    """Get annotation for a specific tool

    Args:
        tool_name: Name of the tool

    Returns:
        ToolAnnotation if found, None otherwise
    """
    return TOOL_ANNOTATIONS.get(tool_name)


def get_tool_with_annotations(tool_info: Dict[str, Any]) -> Dict[str, Any]:
    """Enhance tool info with MCP annotations

    Args:
        tool_info: Basic tool information

    Returns:
        Tool info enhanced with MCP annotations
    """
    tool_name = tool_info.get("name", "")
    annotation = get_tool_annotation(tool_name)

    if annotation:
        tool_info["annotations"] = annotation.to_dict()

    return tool_info


def get_all_tools_with_annotations(tools: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Enhance all tools with MCP annotations

    Args:
        tools: List of tool information dictionaries

    Returns:
        List of tools enhanced with MCP annotations
    """
    return [get_tool_with_annotations(tool) for tool in tools]
