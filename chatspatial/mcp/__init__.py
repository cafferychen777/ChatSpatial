"""
MCP (Model Context Protocol) utilities for ChatSpatial

This package contains modular components for MCP integration:
- errors: Error handling and MCP error types
- resources: Resource management for datasets and results
- prompts: Prompt templates for common workflows
- annotations: Tool metadata and annotations
"""

from .errors import ErrorType, MCPError, format_mcp_error
from .resources import Resource, ResourceManager
from .prompts import Prompt, PromptArgument, PromptManager
from .annotations import ToolAnnotation, get_tool_annotation, TOOL_ANNOTATIONS

__all__ = [
    # Errors
    'ErrorType', 'MCPError', 'format_mcp_error',
    # Resources
    'Resource', 'ResourceManager',
    # Prompts
    'Prompt', 'PromptArgument', 'PromptManager',
    # Annotations
    'ToolAnnotation', 'get_tool_annotation', 'TOOL_ANNOTATIONS'
]