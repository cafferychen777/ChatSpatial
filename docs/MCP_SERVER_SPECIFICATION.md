# ChatSpatial MCP Server Specification

## Overview

ChatSpatial is a **Model Context Protocol (MCP) server** for spatial transcriptomics analysis. It provides AI assistants with access to comprehensive spatial transcriptomics analysis capabilities through a standardized MCP interface.

**MCP Server Information:**
- **Name**: `ChatSpatial`
- **Version**: `0.3.x`
- **Protocol**: MCP v2024-11-05
- **Capabilities**: `tools`, `resources`, `prompts` (planned)

## MCP Architecture

### Server Capabilities

This server implements the MCP protocol to provide:

1. **Tools (32 implemented)**: Interactive analysis functions for spatial transcriptomics
2. **Resources (integrated)**: Access to datasets, analysis results, and visualizations  
3. **Prompts (planned)**: Workflow automation for common analysis patterns

### Transport Options

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/python",
      "args": ["-m", "chatspatial.server"],
      "env": {
        "PYTHONPATH": "/path/to/chatspatial"
      }
    }
  }
}
```

## MCP Tools Specification

All tools follow MCP `tools/call` request format and return standardized `ToolResult` responses.

### 1. Data Management Tools

#### `load_data`
**Description**: Loads spatial transcriptomics data files and creates an MCP resource.

**MCP Tool Schema**:
```json
{
  "name": "load_data",
  "description": "Load spatial transcriptomics data and register as MCP resource",
  "inputSchema": {
    "type": "object",
    "properties": {
      "data_path": {
        "type": "string",
        "description": "Absolute path to data file or directory"
      },
      "data_type": {
        "type": "string", 
        "enum": ["auto", "10x_visium", "slide_seq", "merfish", "seqfish", "h5ad"],
        "default": "auto",
        "description": "Data format type. 'auto' attempts to infer format."
      },
      "name": {
        "type": "string",
        "description": "Optional dataset name. Auto-generated if not provided."
      }
    },
    "required": ["data_path"]
  }
}
```

**MCP Tool Response**:
```json
{
  "content": [
    {
      "type": "text",
      "text": "Successfully loaded dataset 'mouse_brain_visium' with 15000 cells and 20000 genes."
    },
    {
      "type": "resource",
      "resource": {
        "uri": "spatial://datasets/mouse_brain_visium",
        "name": "mouse_brain_visium", 
        "mimeType": "application/x-anndata"
      }
    }
  ]
}
```

#### `preprocess_data` 
**Description**: Preprocess spatial transcriptomics data with quality control and normalization.

**MCP Tool Schema**:
```json
{
  "name": "preprocess_data",
  "description": "Preprocess spatial data with QC, normalization, and dimensionality reduction",
  "inputSchema": {
    "type": "object", 
    "properties": {
      "data_id": {
        "type": "string",
        "description": "Dataset identifier from load_data"
      },
      "params": {
        "type": "object",
        "properties": {
          "min_genes": {"type": "integer", "default": 200},
          "max_genes": {"type": "integer", "default": 5000},
          "min_cells": {"type": "integer", "default": 3},
          "mt_threshold": {"type": "number", "default": 20.0},
          "normalization": {
            "type": "string",
            "enum": ["log", "sct", "scvi"],
            "default": "log"
          },
          "n_pcs": {"type": "integer", "default": 50},
          "clustering_resolution": {"type": "number", "default": 0.5}
        }
      }
    },
    "required": ["data_id"]
  }
}
```

**MCP Tool Response**:
```json
{
  "content": [
    {
      "type": "text", 
      "text": "Preprocessing completed. Filtered to 12500 cells and 18000 genes. Identified 2000 highly variable genes."
    }
  ]
}
```

### 2. Visualization Tools

#### `visualize_data`
**Description**: Generate spatial transcriptomics visualizations as MCP image content.

**MCP Tool Schema**:
```json
{
  "name": "visualize_data", 
  "description": "Create spatial transcriptomics visualizations",
  "inputSchema": {
    "type": "object",
    "properties": {
      "data_id": {
        "type": "string",
        "description": "Dataset identifier"
      },
      "params": {
        "type": "object", 
        "properties": {
          "plot_type": {
            "type": "string",
            "enum": [
              "spatial", "umap", "violin", "heatmap", "dotplot",
              "spatial_domains", "cell_communication", "deconvolution", 
              "trajectory", "gaston_isodepth", "gaston_domains"
            ],
            "description": "Type of visualization to generate"
          },
          "feature": {
            "type": "string", 
            "description": "Gene or feature to visualize"
          },
          "colormap": {"type": "string", "default": "viridis"},
          "figure_size": {
            "type": "array",
            "items": {"type": "integer"},
            "minItems": 2,
            "maxItems": 2
          }
        }
      }
    },
    "required": ["data_id"]
  }
}
```

**MCP Tool Response**:
```json
{
  "content": [
    {
      "type": "image",
      "data": "<base64_encoded_png_data>",
      "mimeType": "image/png"
    },
    {
      "type": "resource",
      "resource": {
        "uri": "spatial://visualizations/mouse_brain_spatial_1673892000", 
        "name": "Spatial plot - Foxp2",
        "mimeType": "image/png"
      }
    }
  ]
}
```

### 3. Analysis Tools

#### `annotate_cells`
**Description**: Perform automated cell type annotation using multiple methods.

**MCP Tool Schema**:
```json
{
  "name": "annotate_cells",
  "description": "Annotate cell types in spatial transcriptomics data", 
  "inputSchema": {
    "type": "object",
    "properties": {
      "data_id": {"type": "string"},
      "params": {
        "type": "object",
        "properties": {
          "method": {
            "type": "string",
            "enum": [
              "marker_genes", "correlation", "tangram", "scanvi", 
              "cellassign", "mllmcelltype", "supervised"
            ],
            "default": "marker_genes"
          },
          "reference_data_id": {
            "type": "string",
            "description": "Reference single-cell dataset (required for tangram, scanvi)"
          },
          "marker_genes": {
            "type": "object",
            "description": "Cell type marker genes dictionary"
          }
        }
      }
    },
    "required": ["data_id"]
  }
}
```

#### `find_spatial_genes`
**Description**: Identify spatially variable genes using deep learning and statistical methods.

**MCP Tool Schema**:
```json
{
  "name": "find_spatial_genes",
  "description": "Identify spatially variable genes using GASTON, SpatialDE, or SPARK",
  "inputSchema": {
    "type": "object",
    "properties": {
      "data_id": {"type": "string"},
      "params": {
        "type": "object", 
        "properties": {
          "method": {
            "type": "string",
            "enum": ["gaston", "spatialde", "spark"],
            "default": "gaston",
            "description": "Spatial gene detection method"
          },
          "n_genes": {
            "type": "integer", 
            "description": "Maximum number of genes to analyze"
          },
          "preprocessing_method": {
            "type": "string",
            "enum": ["glmpca", "pearson_residuals"],
            "default": "pearson_residuals"
          }
        }
      }
    },
    "required": ["data_id"]
  }
}
```

#### `analyze_cell_communication`
**Description**: Analyze cell-cell communication patterns using LIANA+ framework.

**MCP Tool Schema**:
```json
{
  "name": "analyze_cell_communication",
  "description": "Analyze spatial cell-cell communication patterns",
  "inputSchema": {
    "type": "object",
    "properties": {
      "data_id": {"type": "string"},
      "params": {
        "type": "object",
        "properties": {
          "method": {
            "type": "string", 
            "enum": ["liana"],
            "default": "liana"
          },
          "analysis_mode": {
            "type": "string",
            "enum": ["global", "cluster", "spatial_bivariate"],
            "default": "global" 
          },
          "sender_cell_type": {"type": "string"},
          "receiver_cell_type": {"type": "string"}
        }
      }
    },
    "required": ["data_id"]
  }
}
```

### 4. Advanced Analysis Tools

#### `identify_spatial_domains`
**Description**: Identify spatial domains and tissue architecture using graph-based methods.

#### `deconvolve_data` 
**Description**: Deconvolve spatial spots to estimate cell type compositions using scvi-tools methods.

#### `integrate_samples`
**Description**: Integrate multiple spatial transcriptomics samples using Harmony, scVI, or other methods.

#### `analyze_trajectory_data`
**Description**: Infer cellular trajectories and pseudotime using Palantir, CellRank, or DPT.

#### `analyze_enrichment`
**Description**: Perform pathway enrichment analysis using GSEA, ORA, or Enrichr.

## MCP Resources Integration

ChatSpatial creates and manages the following MCP resource types:

### Dataset Resources
```
URI Pattern: spatial://datasets/{data_id}
MIME Type: application/x-anndata
Content: Dataset metadata and summary statistics
```

### Analysis Result Resources  
```
URI Pattern: spatial://results/{data_id}/{analysis_type}
MIME Type: application/json
Content: Structured analysis results
```

### Visualization Resources
```
URI Pattern: spatial://visualizations/{viz_id} 
MIME Type: image/png
Content: Generated plot images
```

## Error Handling

All tools implement standardized MCP error handling:

```json
{
  "isError": true,
  "content": [
    {
      "type": "text",
      "text": "Error: Dataset 'invalid_id' not found. Available datasets: ['mouse_brain', 'human_liver']"
    }
  ]
}
```

**Error Categories**:
- `DataNotFoundError`: Dataset or result not found
- `InvalidParameterError`: Invalid parameter values
- `ProcessingError`: Analysis execution failures  
- `DependencyError`: Missing required packages

## Usage Examples

### Basic Spatial Analysis Workflow

1. **Load Data**:
```json
{
  "method": "tools/call",
  "params": {
    "name": "load_data", 
    "arguments": {
      "data_path": "/path/to/visium_data",
      "data_type": "10x_visium"
    }
  }
}
```

2. **Preprocess**:
```json
{
  "method": "tools/call",
  "params": {
    "name": "preprocess_data",
    "arguments": {
      "data_id": "visium_mouse_brain",
      "params": {
        "min_genes": 200,
        "normalization": "log"
      }
    }
  }
}
```

3. **Visualize**:
```json
{
  "method": "tools/call", 
  "params": {
    "name": "visualize_data",
    "arguments": {
      "data_id": "visium_mouse_brain",
      "params": {
        "plot_type": "spatial",
        "feature": "Foxp2"
      }
    }
  }
}
```

## Integration with AI Assistants

ChatSpatial's MCP interface enables natural language interaction:

**User**: "Load the Visium mouse brain data and show me the spatial expression of Foxp2"

**AI Assistant** → MCP calls:
1. `load_data(data_path="/data/mouse_brain", data_type="10x_visium")`  
2. `visualize_data(data_id="result", params={"plot_type": "spatial", "feature": "Foxp2"})`

The MCP protocol ensures seamless integration with Claude Desktop, Continue, and other MCP-compatible AI tools.

## Performance Considerations

- **Timeout Handling**: All tools respect MCP timeout limits (typically 4 minutes for Claude Desktop)
- **Large Dataset Support**: Automatic chunking and progress reporting for large analyses
- **Resource Caching**: Efficient caching of visualizations and analysis results
- **Memory Management**: Smart memory usage for large spatial datasets

## Dependencies

**Core Requirements**:
- Python 3.10+
- scanpy, anndata, pandas, numpy
- matplotlib, seaborn
- mcp (Model Context Protocol SDK)

**Optional Analysis Methods**:
- scvi-tools (for advanced deep learning methods)
- liana (for cell communication analysis) 
- gaston-spatial (for GASTON spatial gene detection)
- R + rpy2 (for SPARK analysis)

## Contributing

ChatSpatial follows MCP best practices for extensibility. New analysis methods can be added by:

1. Implementing the analysis function in `chatspatial/tools/`
2. Adding MCP tool decorator with proper schema
3. Registering in the server module
4. Adding documentation following this specification format

---

*This specification covers ChatSpatial v0.3.x implementing MCP protocol version 2024-11-05*