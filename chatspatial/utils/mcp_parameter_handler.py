"""
MCP Parameter Handler - Handle parameter validation at the MCP level.

This module provides a way to bypass FastMCP's automatic Pydantic validation
and handle it ourselves with better error messages.
"""

from typing import Dict, Any, Union, Optional
from pydantic import BaseModel, ValidationError
import json

from .tool_error_handling import create_error_result


def validate_parameters_manually(
    params_dict: Union[Dict[str, Any], BaseModel, None],
    model_class: type,
    param_name: str = "params"
) -> Union[BaseModel, Dict[str, Any]]:
    """
    Manually validate parameters with user-friendly error messages.
    
    This function takes raw parameter data and validates it against a Pydantic model,
    returning either the validated model or a properly formatted error.
    
    Args:
        params_dict: Raw parameter data (dict, model instance, or None)
        model_class: Pydantic model class to validate against
        param_name: Name of the parameter for error messages
        
    Returns:
        Validated Pydantic model instance
        
    Raises:
        ValueError: With user-friendly error message if validation fails
    """
    # Handle None or empty params
    if params_dict is None:
        try:
            return model_class()
        except ValidationError as e:
            raise ValueError(f"Default parameter validation failed: {format_pydantic_errors(e.errors())}")
    
    # If already a model instance, return as-is
    if isinstance(params_dict, model_class):
        return params_dict
    
    # Handle dict input
    if isinstance(params_dict, dict):
        try:
            return model_class(**params_dict)
        except ValidationError as e:
            error_msg = f"Parameter '{param_name}' validation failed:\n{format_pydantic_errors(e.errors())}"
            raise ValueError(error_msg)
    
    # Handle other types
    try:
        return model_class(params_dict)
    except (ValidationError, TypeError) as e:
        error_msg = f"Parameter '{param_name}' type error: expected dict or {model_class.__name__} instance, got {type(params_dict).__name__}"
        raise ValueError(error_msg)


def format_pydantic_errors(errors: list) -> str:
    """Format Pydantic validation errors into user-friendly messages"""
    formatted_errors = []
    
    for error in errors:
        field_path = ".".join(str(x) for x in error["loc"]) if error["loc"] else "root field"
        error_type = error["type"]
        error_input = error.get("input")
        ctx = error.get("ctx", {})
        
        # Create user-friendly error messages in English
        if error_type == "greater_than":
            min_value = ctx.get("gt", 0)
            message = f"  • {field_path}: must be greater than {min_value}, got: {error_input}"
        elif error_type == "greater_than_equal":
            min_value = ctx.get("ge", 0)  
            message = f"  • {field_path}: must be greater than or equal to {min_value}, got: {error_input}"
        elif error_type == "less_than":
            max_value = ctx.get("lt")
            message = f"  • {field_path}: must be less than {max_value}, got: {error_input}"
        elif error_type == "less_than_equal":
            max_value = ctx.get("le")
            message = f"  • {field_path}: must be less than or equal to {max_value}, got: {error_input}"
        elif error_type == "type_error":
            expected_type = ctx.get("expected_type", "unknown")
            message = f"  • {field_path}: type error, expected: {expected_type}, got: {type(error_input).__name__}"
        elif error_type == "missing":
            message = f"  • {field_path}: required parameter is missing"
        elif error_type == "literal_error":
            expected_values = ctx.get("expected", [])
            message = f"  • {field_path}: invalid value, allowed: {expected_values}, got: {error_input}"
        elif error_type == "string_too_short":
            min_length = ctx.get("min_length", 0)
            message = f"  • {field_path}: string too short, minimum {min_length} characters, got: {len(str(error_input)) if error_input else 0}"
        elif error_type == "string_too_long":
            max_length = ctx.get("max_length", 0)
            message = f"  • {field_path}: string too long, maximum {max_length} characters, got: {len(str(error_input)) if error_input else 0}"
        else:
            # Fallback for unknown error types
            message = f"  • {field_path}: {error.get('msg', 'validation failed')} (value: {error_input})"
        
        formatted_errors.append(message)
    
    return "\n".join(formatted_errors)


# Common parameter validation helpers
def validate_analysis_params(params: Union[Dict[str, Any], BaseModel, None]):
    """Validate AnalysisParameters with friendly error messages"""
    from ..models.data import AnalysisParameters
    return validate_parameters_manually(params, AnalysisParameters, "analysis_params")


def validate_visualization_params(params: Union[Dict[str, Any], BaseModel, None]):
    """Validate VisualizationParameters with friendly error messages and preprocessing"""
    from ..models.data import VisualizationParameters

    # Preprocess parameters to handle different input formats
    preprocessed_params = _preprocess_visualization_params(params)

    return validate_parameters_manually(preprocessed_params, VisualizationParameters, "visualization_params")


def _preprocess_visualization_params(params: Union[Dict[str, Any], BaseModel, str, None]) -> Union[Dict[str, Any], BaseModel, None]:
    """
    Preprocess visualization parameters to handle different input formats.

    Handles:
    - None: Returns empty dict
    - str: Converts to feature parameter (supports "gene:CCL21" and "CCL21" formats)
    - dict: Normalizes features/feature naming
    - VisualizationParameters: Returns as-is
    """
    if params is None:
        return {}

    if isinstance(params, str):
        # Handle string format parameters
        if params.startswith("gene:"):
            feature = params.split(":", 1)[1]
            return {"feature": feature, "plot_type": "spatial"}
        else:
            return {"feature": params, "plot_type": "spatial"}

    if isinstance(params, dict):
        # Handle dictionary format parameters
        result = params.copy()

        # features/feature parameter handled by Pydantic model validator

        return result

    # For VisualizationParameters instances or other types, return as-is
    return params


def validate_spatial_analysis_params(params: Union[Dict[str, Any], BaseModel, None]):
    """Validate SpatialAnalysisParameters with friendly error messages"""
    from ..models.data import SpatialAnalysisParameters
    return validate_parameters_manually(params, SpatialAnalysisParameters, "spatial_analysis_params")


def validate_cell_communication_params(params: Union[Dict[str, Any], BaseModel, None]):
    """Validate CellCommunicationParameters with friendly error messages"""
    from ..models.data import CellCommunicationParameters
    return validate_parameters_manually(params, CellCommunicationParameters, "cell_communication_params")


def validate_annotation_params(params: Union[Dict[str, Any], BaseModel, None]):
    """Validate AnnotationParameters with friendly error messages"""
    from ..models.data import AnnotationParameters
    return validate_parameters_manually(params, AnnotationParameters, "annotation_params")


# Decorator for tools that need manual parameter validation
def manual_parameter_validation(*param_validators):
    """
    Decorator for manual parameter validation.
    
    Usage:
    @manual_parameter_validation(
        ("params", validate_analysis_params)
    )
    async def my_tool(data_id: str, params: Any = None):
        # params is now validated and converted to proper type
        ...
    """
    def decorator(func):
        from functools import wraps
        import inspect
        
        @wraps(func)
        async def async_wrapper(*args, **kwargs):
            # Get function signature
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            
            # Apply parameter validators
            for param_name, validator in param_validators:
                if param_name in bound_args.arguments:
                    try:
                        validated_param = validator(bound_args.arguments[param_name])
                        bound_args.arguments[param_name] = validated_param
                    except ValueError as e:
                        # Re-raise for the outer error handler to catch
                        raise e
            
            # Call the original function with validated parameters
            return await func(**bound_args.arguments)
        
        @wraps(func)
        def sync_wrapper(*args, **kwargs):
            # Get function signature
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            
            # Apply parameter validators
            for param_name, validator in param_validators:
                if param_name in bound_args.arguments:
                    try:
                        validated_param = validator(bound_args.arguments[param_name])
                        bound_args.arguments[param_name] = validated_param
                    except ValueError as e:
                        # Re-raise for the outer error handler to catch
                        raise e
            
            # Call the original function with validated parameters
            return func(**bound_args.arguments)
        
        # Return appropriate wrapper
        if inspect.iscoroutinefunction(func):
            return async_wrapper
        else:
            return sync_wrapper
    
    return decorator