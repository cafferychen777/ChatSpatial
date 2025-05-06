# ChatSpatial - Interactive Spatial Transcriptomics Assistant

ChatSpatial is an interactive spatial transcriptomics data analysis assistant based on the Model Context Protocol (MCP), providing a suite of tools for spatial transcriptomics data processing, visualization, and analysis. It helps researchers analyze spatial transcriptomics data through natural language dialogue.

## Features

- **Data Loading**: Support for various spatial transcriptomics data formats (10x Visium, Slide-seq, MERFISH, seqFISH, etc.)
- **Data Preprocessing**: Normalization, batch effect correction, dimensionality reduction, etc.
- **Spatial Visualization**: Spatial distribution of gene expression, visualization of clustering results, etc.
- **Differential Expression Analysis**: Identification of differentially expressed genes between cell populations
- **Cell Type Annotation**: Multiple cell type annotation methods, including marker-based, correlation-based, and supervised classification
- **Spatial Analysis**: Spatial autocorrelation, neighborhood analysis, spatial trajectories, etc.
- **Deconvolution Analysis**: Support for NNLS, Cell2location, and Spotiphy methods with enhanced error handling and user feedback
- **Standardized Image Processing**: Unified image processing module ensuring all visualization functions return standardized Image objects

## Installation

```bash
# Clone the repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial

# Install basic dependencies
pip install -e .

# Install all optional dependencies
pip install -e .[all]

# Or install specific optional dependencies
pip install -e .[deconvolution]  # Install deconvolution-related dependencies
```

## Usage

### Starting the Server

```bash
# Using stdio transport (default)
chatspatial

# Using SSE transport with a specified port
chatspatial --transport sse --port 8000
```

### Client Example

```python
import asyncio
from mcp.client.session import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client

async def main():
    async with stdio_client(
        StdioServerParameters(command="chatspatial")
    ) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()

            # List available tools
            tools = await session.list_tools()
            print(tools)

            # Load data
            result = await session.call_tool("load_data", {
                "data_path": "/path/to/spatial_data.h5ad",
                "data_type": "10x_visium",
                "name": "Mouse Brain"
            })
            print(result)

            # Preprocess data
            result = await session.call_tool("preprocess", {
                "data_id": result["id"],
                "params": {
                    "normalization": "log",
                    "n_hvgs": 2000,
                    "n_pcs": 30
                }
            })
            print(result)

            # Visualize data
            result = await session.call_tool("visualize", {
                "data_id": "data_1",
                "params": {
                    "feature": "Cd8a",
                    "plot_type": "spatial",
                    "colormap": "viridis"
                }
            })
            # Process the returned image...

asyncio.run(main())
```

## Available Tools

The server provides the following tools:

1. `load_data` - Load spatial transcriptomics data
2. `preprocess` - Preprocess data
3. `visualize` - Visualize data
4. `annotate` - Cell type annotation
5. `analyze_spatial` - Spatial analysis
6. `find_markers` - Differential expression analysis
7. `list_datasets` - List loaded datasets
8. `integrate_samples` - Multi-sample integration
9. `analyze_trajectory` - Trajectory analysis
10. `analyze_velocity` - RNA velocity analysis
11. `deconvolve` - Spatial transcriptomics deconvolution

## Resources

The server provides the following resources:

- `dataset://{data_id}` - Get dataset information

## Dependencies

- mcp - Model Context Protocol Python SDK
- numpy, pandas - Data processing
- matplotlib - Visualization
- scanpy - Single-cell data analysis
- squidpy - Spatial transcriptomics analysis
- anndata - AnnData data structure
- scikit-learn - Machine learning algorithms
- cell2location (optional) - Spatial transcriptomics deconvolution
- Spotiphy (optional) - Spatial transcriptomics deconvolution
- torch, pyro-ppl (optional) - For Spotiphy and cell2location
- cellrank (optional) - Trajectory analysis
- scvelo (optional) - RNA velocity analysis

## License

MIT

## Contributing

Pull Requests and Issues are welcome!
