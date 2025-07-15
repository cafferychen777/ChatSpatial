"""
示例：如何启用 MCP 增强功能（使用新的 mcp 模块）
"""

# 1. 在 server.py 中添加资源管理
async def update_resources(mcp, data_store):
    """更新可用资源列表"""
    from chatspatial.mcp.resources import ResourceManager
    
    resource_manager = ResourceManager(data_store)
    resources = resource_manager.get_resource_list()
    
    for resource in resources:
        # 注册资源到 MCP
        await mcp.add_resource(
            uri=resource.uri,
            name=resource.name,
            mimeType=resource.mimeType,
            description=resource.description
        )

# 2. 在数据加载后调用
async def enhanced_load_data(file_path: str, params: dict, context):
    # 原有的加载逻辑
    result = await load_spatial_data(file_path, params)
    
    # 更新资源
    await update_resources(context.mcp, context.data_store)
    
    return result

# 3. 启用提示功能
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
        # 转换为工具调用
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

# 4. 使用工具注解
from chatspatial.mcp.annotations import get_tool_annotation

def get_enhanced_tool_info(tool_name: str) -> dict:
    """获取带有 MCP 注解的工具信息"""
    annotation = get_tool_annotation(tool_name)
    
    tool_info = {
        "name": tool_name,
        # ... 其他工具信息
    }
    
    if annotation:
        tool_info["annotations"] = annotation.to_dict()
    
    return tool_info