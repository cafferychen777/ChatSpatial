"""
MCP Prompt System

This module provides prompt templates for common spatial analysis workflows.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


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
            "required": self.required,
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
            "arguments": (
                [arg.to_dict() for arg in self.arguments] if self.arguments else []
            ),
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
                    PromptArgument("method", "Analysis method", required=False),
                ],
            ),
            Prompt(
                name="find-cell-types",
                description="Identify cell types in spatial data",
                arguments=[
                    PromptArgument(
                        "method", "Cell type identification method", required=False
                    ),
                    PromptArgument(
                        "reference_data", "Reference dataset path", required=False
                    ),
                ],
            ),
            Prompt(
                name="compare-conditions",
                description="Compare spatial patterns between conditions",
                arguments=[
                    PromptArgument(
                        "condition_key", "Column defining conditions", required=True
                    ),
                    PromptArgument("condition1", "First condition", required=True),
                    PromptArgument("condition2", "Second condition", required=True),
                ],
            ),
            Prompt(
                name="generate-visualization",
                description="Generate spatial visualization",
                arguments=[
                    PromptArgument("plot_type", "Type of visualization", required=True),
                    PromptArgument(
                        "feature",
                        "Feature(s) to visualize (single gene or list of genes)",
                        required=False,
                    ),
                    PromptArgument("save_path", "Path to save figure", required=False),
                ],
            ),
            Prompt(
                name="quality-control",
                description="Perform quality control on spatial data",
                arguments=[
                    PromptArgument("metrics", "QC metrics to compute", required=False),
                    PromptArgument("thresholds", "QC thresholds", required=False),
                ],
            ),
            Prompt(
                name="batch-correction",
                description="Correct batch effects in integrated data",
                arguments=[
                    PromptArgument(
                        "batch_key", "Column defining batches", required=True
                    ),
                    PromptArgument("method", "Batch correction method", required=False),
                ],
            ),
            Prompt(
                name="spatial-clustering",
                description="Perform spatial clustering analysis",
                arguments=[
                    PromptArgument("n_clusters", "Number of clusters", required=False),
                    PromptArgument(
                        "resolution", "Clustering resolution", required=False
                    ),
                ],
            ),
            Prompt(
                name="trajectory-inference",
                description="Infer cellular trajectories",
                arguments=[
                    PromptArgument("start_cell", "Starting cell type", required=False),
                    PromptArgument("end_cell", "Ending cell type", required=False),
                ],
            ),
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
