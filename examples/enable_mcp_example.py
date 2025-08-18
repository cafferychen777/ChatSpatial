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
    
    # Generate structured prompt message based on template and arguments
    message_content = _generate_prompt_message(name, arguments or {})
    
    return {
        "description": prompt.description,
        "messages": [{
            "role": "user",
            "content": message_content
        }]
    }

def _generate_prompt_message(prompt_name: str, arguments: dict) -> str:
    """Generate structured prompt message content based on template and arguments"""
    
    if prompt_name == "analyze-spatial-expression":
        genes = arguments.get("genes", "")
        method = arguments.get("method", "Moran's I")
        if genes:
            return f"Analyze spatial expression patterns for genes: {genes}. Use {method} statistic to identify spatial variability and create visualizations."
        else:
            return "Analyze spatial gene expression patterns using statistical methods to identify spatially variable genes."
    
    elif prompt_name == "find-cell-types":
        method = arguments.get("method", "automated annotation")
        reference = arguments.get("reference_data", "")
        ref_text = f" using reference data: {reference}" if reference else ""
        return f"Identify cell types in spatial transcriptomics data using {method}{ref_text}. Provide detailed cell type annotations with confidence scores."
    
    elif prompt_name == "compare-conditions":
        condition_key = arguments.get("condition_key", "condition")
        cond1 = arguments.get("condition1", "condition A")
        cond2 = arguments.get("condition2", "condition B") 
        return f"Compare spatial patterns between {cond1} and {cond2} using the '{condition_key}' grouping. Identify differentially expressed genes and spatial differences."
    
    elif prompt_name == "generate-visualization":
        plot_type = arguments.get("plot_type", "spatial plot")
        feature = arguments.get("feature", "")
        feature_text = f" for {feature}" if feature else ""
        return f"Generate a {plot_type}{feature_text}. Create high-quality visualizations with proper legends and color scales."
    
    elif prompt_name == "quality-control":
        metrics = arguments.get("metrics", "standard QC metrics")
        return f"Perform quality control analysis using {metrics}. Generate QC plots and filtering recommendations."
    
    elif prompt_name == "batch-correction":
        batch_key = arguments.get("batch_key", "batch")
        method = arguments.get("method", "Harmony")
        return f"Correct batch effects using {method} method with batch information from '{batch_key}' column. Evaluate correction effectiveness."
    
    elif prompt_name == "spatial-clustering":
        n_clusters = arguments.get("n_clusters", "")
        resolution = arguments.get("resolution", "")
        params = []
        if n_clusters: params.append(f"{n_clusters} clusters")
        if resolution: params.append(f"resolution {resolution}")
        param_text = f" with {', '.join(params)}" if params else ""
        return f"Perform spatial clustering analysis{param_text}. Identify spatially coherent domains and visualize results."
    
    elif prompt_name == "trajectory-inference":
        start_cell = arguments.get("start_cell", "")
        end_cell = arguments.get("end_cell", "")
        trajectory_text = ""
        if start_cell and end_cell:
            trajectory_text = f" from {start_cell} to {end_cell}"
        elif start_cell:
            trajectory_text = f" starting from {start_cell}"
        return f"Infer cellular trajectories{trajectory_text}. Compute pseudotime and identify trajectory-associated genes."
    
    else:
        # Fallback for unknown prompts
        return f"Execute {prompt_name} analysis with the provided parameters."

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