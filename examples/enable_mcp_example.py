"""
Example: How to enable MCP enhanced features (using new mcp module)
"""

# 1. Add resource management in server.py
async def update_resources(mcp, data_store):
    """Update available resource list"""
    from chatspatial.mcp.resources import ResourceManager

    resource_manager = ResourceManager(data_store)
    resources = resource_manager.get_resource_list()

    for resource in resources:
        # Register resource to MCP
        await mcp.add_resource(
            uri=resource.uri,
            name=resource.name,
            mimeType=resource.mimeType,
            description=resource.description
        )

# 2. Call after data loading
async def enhanced_load_data(file_path: str, params: dict, context):
    # Original loading logic
    result = await load_spatial_data(file_path, params)

    # Update resources
    await update_resources(context.mcp, context.data_store)

    return result

# 3. Enable prompt functionality
from chatspatial.mcp.prompts import PromptManager

prompt_manager = PromptManager()

@mcp.list_prompts
async def list_prompts():
    return [prompt.to_dict() for prompt in prompt_manager.list_prompts()]

@mcp.get_prompt
async def get_prompt(name: str, arguments: dict = None):
    prompt = prompt_manager.get_prompt(name)
    if not prompt:
        raise ValueError(f"Prompt '{name}' not found")
    
    if arguments:
        # Convert to tool call
        tool_params = prompt_manager.prompt_to_tool_params(name, arguments)
        return {
            "description": prompt.description,
            "messages": [{
                "role": "user",
                "content": f"Run {tool_params['tool']} with {tool_params}"
            }]
        }
    
    return {
        "description": prompt.description,
        "messages": [{
            "role": "user",
            "content": prompt.description
        }]
    }

# 4. Use tool annotations
from chatspatial.mcp.annotations import get_tool_annotation

def get_enhanced_tool_info(tool_name: str) -> dict:
    """Get tool information with MCP annotations"""
    annotation = get_tool_annotation(tool_name)
    
    tool_info = {
        "name": tool_name,
        # ... other tool information
    }
    
    if annotation:
        tool_info["annotations"] = annotation.to_dict()
    
    return tool_info