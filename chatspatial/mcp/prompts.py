"""
MCP Prompt System

This module provides prompt templates for common spatial analysis workflows.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field


@dataclass
class PromptArgument:
    """MCP Prompt Argument"""
    name: str
    description: str
    required: bool = True
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to MCP prompt argument format"""
        return {
            "name": self.name,
            "description": self.description,
            "required": self.required
        }


@dataclass
class Prompt:
    """MCP Prompt representation"""
    name: str
    description: str
    arguments: List[PromptArgument] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to MCP prompt format"""
        return {
            "name": self.name,
            "description": self.description,
            "arguments": [arg.to_dict() for arg in self.arguments] if self.arguments else []
        }


class PromptManager:
    """Manages MCP prompts for spatial analysis"""
    
    def __init__(self):
        self.prompts = self._initialize_prompts()
    
    def _initialize_prompts(self) -> List[Prompt]:
        """Initialize spatial analysis prompts"""
        return [
            Prompt(
                name="analyze-spatial-expression",
                description="Analyze spatial gene expression patterns",
                arguments=[
                    PromptArgument("genes", "Genes to analyze", required=True),
                    PromptArgument("method", "Analysis method", required=False)
                ]
            ),
            Prompt(
                name="find-cell-types",
                description="Identify cell types in spatial data",
                arguments=[
                    PromptArgument("method", "Cell type identification method", required=False),
                    PromptArgument("reference_data", "Reference dataset path", required=False)
                ]
            ),
            Prompt(
                name="compare-conditions",
                description="Compare spatial patterns between conditions",
                arguments=[
                    PromptArgument("condition_key", "Column defining conditions", required=True),
                    PromptArgument("condition1", "First condition", required=True),
                    PromptArgument("condition2", "Second condition", required=True)
                ]
            ),
            Prompt(
                name="generate-visualization",
                description="Generate spatial visualization",
                arguments=[
                    PromptArgument("plot_type", "Type of visualization", required=True),
                    PromptArgument("feature", "Feature(s) to visualize (single gene or list of genes)", required=False),
                    PromptArgument("save_path", "Path to save figure", required=False)
                ]
            ),
            Prompt(
                name="quality-control",
                description="Perform quality control on spatial data",
                arguments=[
                    PromptArgument("metrics", "QC metrics to compute", required=False),
                    PromptArgument("thresholds", "QC thresholds", required=False)
                ]
            ),
            Prompt(
                name="batch-correction",
                description="Correct batch effects in integrated data",
                arguments=[
                    PromptArgument("batch_key", "Column defining batches", required=True),
                    PromptArgument("method", "Batch correction method", required=False)
                ]
            ),
            Prompt(
                name="spatial-clustering",
                description="Perform spatial clustering analysis",
                arguments=[
                    PromptArgument("n_clusters", "Number of clusters", required=False),
                    PromptArgument("resolution", "Clustering resolution", required=False)
                ]
            ),
            Prompt(
                name="trajectory-inference",
                description="Infer cellular trajectories",
                arguments=[
                    PromptArgument("start_cell", "Starting cell type", required=False),
                    PromptArgument("end_cell", "Ending cell type", required=False)
                ]
            )
        ]
    
    def get_prompt(self, name: str) -> Optional[Prompt]:
        """Get prompt by name"""
        for prompt in self.prompts:
            if prompt.name == name:
                return prompt
        return None
    
    def list_prompts(self) -> List[Prompt]:
        """List all available prompts"""
        return self.prompts
    
    def prompt_to_tool_params(self, prompt_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Convert prompt arguments to tool parameters
        
        Args:
            prompt_name: Name of the prompt
            arguments: User-provided arguments
            
        Returns:
            Dictionary of tool parameters
        """
        # Map prompt names to tool functions
        prompt_tool_map = {
            "analyze-spatial-expression": "find_spatial_variable_genes",
            "find-cell-types": "annotate_cell_types",
            "compare-conditions": "perform_differential_expression",
            "generate-visualization": "visualize_data",
            "quality-control": "preprocess_data",
            "batch-correction": "integrate_datasets",
            "spatial-clustering": "identify_spatial_domains",
            "trajectory-inference": "analyze_trajectory"
        }
        
        tool_name = prompt_tool_map.get(prompt_name)
        if not tool_name:
            raise ValueError(f"No tool mapping for prompt: {prompt_name}")
        
        # Transform arguments based on prompt type
        tool_params = {"tool": tool_name}
        
        # Handle prompt-specific parameter mappings
        if prompt_name == "analyze-spatial-expression":
            if "genes" in arguments:
                tool_params["genes"] = [g.strip() for g in arguments["genes"].split(",")]
            if "method" in arguments:
                tool_params["method"] = arguments["method"]
        
        elif prompt_name == "find-cell-types":
            if "method" in arguments:
                tool_params["method"] = arguments["method"]
            if "reference_data" in arguments:
                tool_params["reference_data"] = arguments["reference_data"]
        
        elif prompt_name == "compare-conditions":
            tool_params["groupby"] = arguments.get("condition_key")
            tool_params["groups"] = [arguments.get("condition1"), arguments.get("condition2")]
        
        elif prompt_name == "generate-visualization":
            tool_params["plot_type"] = arguments.get("plot_type")
            if "feature" in arguments:
                tool_params["feature"] = arguments["feature"]
            # Support legacy 'features' parameter for backward compatibility
            elif "features" in arguments:
                tool_params["feature"] = arguments["features"]
            if "save_path" in arguments:
                tool_params["save_path"] = arguments["save_path"]
        
        elif prompt_name == "spatial-clustering":
            if "n_clusters" in arguments:
                tool_params["n_clusters"] = int(arguments["n_clusters"])
            if "resolution" in arguments:
                tool_params["resolution"] = float(arguments["resolution"])
        
        elif prompt_name == "trajectory-inference":
            if "start_cell" in arguments:
                tool_params["start_cell"] = arguments["start_cell"]
            if "end_cell" in arguments:
                tool_params["end_cell"] = arguments["end_cell"]
        
        # Add any additional arguments not explicitly mapped
        for key, value in arguments.items():
            if key not in tool_params and key != "prompt":
                tool_params[key] = value
        
        return tool_params