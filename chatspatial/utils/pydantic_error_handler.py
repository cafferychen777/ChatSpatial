"""
Pydantic validation error handler for MCP tools.

This module provides utilities to catch Pydantic validation errors
and convert them to MCP-compliant error format before they reach
the MCP framework.
"""

from typing import Any, Dict, Type, Optional, get_type_hints
from functools import wraps
import inspect
from pydantic import BaseModel, ValidationError

from .tool_error_handling import create_error_result, invalid_parameter_error


def get_pydantic_models_from_signature(func) -> Dict[str, Type[BaseModel]]:
    """Extract Pydantic models from function signature"""
    sig = inspect.signature(func)
    type_hints = get_type_hints(func)
    
    pydantic_params = {}
    for param_name, param in sig.parameters.items():
        param_type = type_hints.get(param_name, param.annotation)
        
        # Check if it's a Pydantic model
        if (inspect.isclass(param_type) and 
            issubclass(param_type, BaseModel)):
            pydantic_params[param_name] = param_type
            
    return pydantic_params


def validate_pydantic_params(func, args, kwargs) -> Dict[str, Any]:
    """Validate and convert Pydantic parameters with better error messages"""
    sig = inspect.signature(func)
    bound_args = sig.bind(*args, **kwargs)
    bound_args.apply_defaults()
    
    pydantic_models = get_pydantic_models_from_signature(func)
    
    for param_name, model_class in pydantic_models.items():
        if param_name in bound_args.arguments:
            param_value = bound_args.arguments[param_name]
            
            # Skip if already a model instance
            if isinstance(param_value, model_class):
                continue
            
            # Validate and convert
            try:
                if isinstance(param_value, dict):
                    # Convert dict to model
                    validated_model = model_class(**param_value)
                    bound_args.arguments[param_name] = validated_model
                elif param_value is None and param_name in sig.parameters:
                    # Use default if available
                    default = sig.parameters[param_name].default
                    if default != inspect.Parameter.empty:
                        if isinstance(default, model_class):
                            bound_args.arguments[param_name] = default
                        else:
                            bound_args.arguments[param_name] = model_class()
                            
            except ValidationError as e:
                # Create user-friendly error message
                error_details = []
                for error in e.errors():
                    field_path = ".".join(str(x) for x in error["loc"])
                    field_name = field_path if field_path else param_name
                    
                    # Map Pydantic error types to user-friendly messages
                    error_msg = format_validation_error(error, field_name, param_value)
                    error_details.append(error_msg)
                
                # Create comprehensive error message
                full_error = f"参数验证失败 - {param_name}:\n" + "\n".join(f"  • {detail}" for detail in error_details)
                
                # Return the error in the proper format
                raise ValueError(full_error)
    
    return bound_args.arguments


def format_validation_error(error: Dict[str, Any], field_name: str, param_value: Any) -> str:
    """Format Pydantic validation error into user-friendly message"""
    error_type = error["type"]
    error_input = error.get("input", param_value)
    
    # Map common validation errors to Chinese messages
    error_messages = {
        "greater_than": f"'{field_name}' 必须大于 {error.get('ctx', {}).get('gt', 0)}，当前值: {error_input}",
        "greater_than_equal": f"'{field_name}' 必须大于等于 {error.get('ctx', {}).get('ge', 0)}，当前值: {error_input}",
        "less_than": f"'{field_name}' 必须小于 {error.get('ctx', {}).get('lt', 0)}，当前值: {error_input}",
        "less_than_equal": f"'{field_name}' 必须小于等于 {error.get('ctx', {}).get('le', 0)}，当前值: {error_input}",
        "type_error": f"'{field_name}' 类型错误，期望: {error.get('ctx', {}).get('expected_type', 'unknown')}，实际: {type(error_input).__name__}",
        "value_error": f"'{field_name}' 值错误: {error.get('msg', 'Invalid value')}",
        "missing": f"缺少必需参数: '{field_name}'",
        "string_too_short": f"'{field_name}' 长度不足，最少需要 {error.get('ctx', {}).get('min_length', 0)} 个字符",
        "string_too_long": f"'{field_name}' 长度过长，最多允许 {error.get('ctx', {}).get('max_length', 0)} 个字符",
        "literal_error": f"'{field_name}' 值无效，允许的值: {error.get('ctx', {}).get('expected', [])}",
    }
    
    return error_messages.get(error_type, f"'{field_name}' 验证失败: {error.get('msg', 'Unknown error')}")


def mcp_pydantic_error_handler():
    """
    Decorator to handle Pydantic validation errors in MCP tools.
    
    This decorator should be used in combination with @mcp_tool_error_handler():
    
    @mcp.tool()
    @mcp_tool_error_handler()
    @mcp_pydantic_error_handler()
    async def my_tool(data_id: str, params: MyParameters = MyParameters()):
        ...
    """
    def decorator(func):
        @wraps(func)
        async def async_wrapper(*args, **kwargs):
            try:
                # Validate Pydantic parameters first
                validated_kwargs = validate_pydantic_params(func, args, kwargs)
                
                # Call the original function with validated parameters
                return await func(**validated_kwargs)
                
            except ValueError as e:
                # Re-raise ValueError for the outer error handler to catch
                raise e
            except ValidationError as e:
                # This shouldn't happen if validate_pydantic_params works correctly,
                # but include as fallback
                error_msg = f"参数验证失败: {str(e)}"
                raise ValueError(error_msg)
        
        @wraps(func)
        def sync_wrapper(*args, **kwargs):
            try:
                # Validate Pydantic parameters first
                validated_kwargs = validate_pydantic_params(func, args, kwargs)
                
                # Call the original function with validated parameters
                return func(**validated_kwargs)
                
            except ValueError as e:
                # Re-raise ValueError for the outer error handler to catch
                raise e
            except ValidationError as e:
                # Fallback for validation errors
                error_msg = f"参数验证失败: {str(e)}"
                raise ValueError(error_msg)
        
        # Return appropriate wrapper based on function type
        if inspect.iscoroutinefunction(func):
            return async_wrapper
        else:
            return sync_wrapper
    
    return decorator


# Alternative approach: Create custom parameter validation function
def validate_analysis_parameters(params_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate analysis parameters with custom logic and user-friendly error messages.
    
    This is an alternative to Pydantic validation that provides more control
    over error messages and validation logic.
    """
    from ..models.data import AnalysisParameters
    
    # Common parameter validation rules
    validation_rules = {
        "n_pcs": {
            "type": int,
            "min": 1,
            "max": 100,
            "description": "主成分数量必须在 1-100 之间"
        },
        "n_hvgs": {
            "type": int, 
            "min": 100,
            "max": 5000,
            "description": "高变基因数量必须在 100-5000 之间"
        },
        "subsample_spots": {
            "type": int,
            "min": 1,
            "max": 50000,
            "description": "spot 子采样数量必须在 1-50000 之间",
            "optional": True
        },
        "subsample_genes": {
            "type": int,
            "min": 1, 
            "max": 50000,
            "description": "基因子采样数量必须在 1-50000 之间",
            "optional": True
        }
    }
    
    errors = []
    
    for param_name, rules in validation_rules.items():
        if param_name in params_dict:
            value = params_dict[param_name]
            
            # Skip None values for optional parameters
            if value is None and rules.get("optional", False):
                continue
                
            # Type validation
            if not isinstance(value, rules["type"]):
                errors.append(f"{param_name}: 类型错误，期望 {rules['type'].__name__}，实际 {type(value).__name__}")
                continue
            
            # Range validation
            if "min" in rules and value < rules["min"]:
                errors.append(f"{param_name}: 值 {value} 小于最小值 {rules['min']}")
            
            if "max" in rules and value > rules["max"]:
                errors.append(f"{param_name}: 值 {value} 大于最大值 {rules['max']}")
    
    if errors:
        error_msg = "参数验证失败:\n" + "\n".join(f"  • {error}" for error in errors)
        raise ValueError(error_msg)
    
    # If validation passes, try to create the Pydantic model
    try:
        return AnalysisParameters(**params_dict)
    except ValidationError as e:
        # Fallback to Pydantic's error message if our validation missed something
        error_msg = f"参数验证失败: {str(e)}"
        raise ValueError(error_msg)