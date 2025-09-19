---
layout: default
title: Data Models
parent: API Reference
grand_parent: Reference
nav_order: 1
description: "Data models and schemas used throughout ChatSpatial"
---

# Data Models

This document describes the data models, schemas, and parameter structures used throughout ChatSpatial MCP tools.

## Overview

ChatSpatial uses strongly-typed data models to ensure consistent communication between the MCP client and server. All parameters and return values are validated against JSON schemas.

## Common Data Types

### DataID
A unique identifier for loaded spatial datasets.

```json
{
  "type": "string",
  "description": "Unique identifier for a spatial dataset",
  "pattern": "^[a-zA-Z0-9_-]+$"
}
```

### Coordinates
Spatial coordinate information.

```json
{
  "type": "object",
  "properties": {
    "x": {"type": "number"},
    "y": {"type": "number"},
    "z": {"type": "number", "optional": true}
  },
  "required": ["x", "y"]
}
```

## Tool Parameters

### load_data Parameters

```json
{
  "type": "object",
  "properties": {
    "data_path": {
      "type": "string",
      "description": "Path to the spatial data file (.h5ad, .h5, .csv, etc.)"
    },
    "name": {
      "type": "string",
      "description": "Identifier name for the dataset"
    },
    "data_type": {
      "type": "string",
      "enum": ["visium", "merfish", "slide_seq", "starmap", "seqfish", "xenium"],
      "description": "Type of spatial transcriptomics technology"
    }
  },
  "required": ["data_path", "name"]
}
```

### preprocess_data Parameters

```json
{
  "type": "object",
  "properties": {
    "data_id": {
      "type": "string",
      "description": "Dataset identifier"
    },
    "normalization_method": {
      "type": "string",
      "enum": ["log", "sctransform", "pearson_residuals"],
      "default": "log"
    },
    "filter_genes": {
      "type": "boolean",
      "default": true
    },
    "min_cells": {
      "type": "integer",
      "minimum": 1,
      "default": 3
    },
    "min_genes": {
      "type": "integer",
      "minimum": 1,
      "default": 200
    }
  },
  "required": ["data_id"]
}
```

### annotate_cells Parameters

```json
{
  "type": "object",
  "properties": {
    "data_id": {
      "type": "string",
      "description": "Dataset identifier"
    },
    "method": {
      "type": "string",
      "enum": ["marker_based", "tangram", "sctype", "cell2location", "scanvi", "cellassign", "mllmcelltype"],
      "description": "Cell type annotation method"
    },
    "reference_path": {
      "type": "string",
      "description": "Path to reference dataset (required for some methods)"
    },
    "markers": {
      "type": "object",
      "description": "Cell type markers (for marker-based method)"
    }
  },
  "required": ["data_id", "method"]
}
```

### identify_spatial_domains Parameters

```json
{
  "type": "object",
  "properties": {
    "data_id": {
      "type": "string",
      "description": "Dataset identifier"
    },
    "method": {
      "type": "string",
      "enum": ["spagcn", "stagate", "banksy", "leiden", "louvain"],
      "description": "Spatial domain identification method"
    },
    "n_domains": {
      "type": "integer",
      "minimum": 2,
      "description": "Number of spatial domains to identify"
    },
    "resolution": {
      "type": "number",
      "minimum": 0.1,
      "maximum": 2.0,
      "description": "Clustering resolution parameter"
    }
  },
  "required": ["data_id", "method"]
}
```

### analyze_cell_communication Parameters

```json
{
  "type": "object",
  "properties": {
    "data_id": {
      "type": "string",
      "description": "Dataset identifier"
    },
    "method": {
      "type": "string",
      "enum": ["liana", "cellphonedb", "cellchat"],
      "description": "Cell communication analysis method"
    },
    "cell_type_column": {
      "type": "string",
      "description": "Column name containing cell type annotations"
    },
    "min_cells": {
      "type": "integer",
      "minimum": 3,
      "default": 10,
      "description": "Minimum cells per cell type"
    }
  },
  "required": ["data_id", "method", "cell_type_column"]
}
```

### visualize_data Parameters

```json
{
  "type": "object",
  "properties": {
    "data_id": {
      "type": "string",
      "description": "Dataset identifier"
    },
    "plot_type": {
      "type": "string",
      "enum": [
        "spatial", "umap", "tsne", "pca", "heatmap", "violin", 
        "spatial_domains", "cell_communication", "trajectory",
        "gene_expression", "marker_genes", "qc_metrics"
      ],
      "description": "Type of plot to generate"
    },
    "genes": {
      "type": "array",
      "items": {"type": "string"},
      "description": "Genes to visualize (for expression plots)"
    },
    "color_by": {
      "type": "string",
      "description": "Column to use for coloring"
    },
    "save_path": {
      "type": "string",
      "description": "Path to save the plot"
    }
  },
  "required": ["data_id", "plot_type"]
}
```

## Return Value Schemas

### StandardResult

```json
{
  "type": "object",
  "properties": {
    "success": {
      "type": "boolean",
      "description": "Whether the operation was successful"
    },
    "message": {
      "type": "string",
      "description": "Human-readable result message"
    },
    "data": {
      "type": "object",
      "description": "Operation-specific result data"
    },
    "metadata": {
      "type": "object",
      "description": "Additional metadata about the operation"
    }
  },
  "required": ["success", "message"]
}
```

### LoadDataResult

```json
{
  "type": "object",
  "properties": {
    "success": {"type": "boolean"},
    "message": {"type": "string"},
    "data": {
      "type": "object",
      "properties": {
        "id": {"type": "string"},
        "n_obs": {"type": "integer"},
        "n_vars": {"type": "integer"},
        "technology": {"type": "string"},
        "spatial_coords": {"type": "boolean"}
      }
    }
  }
}
```

### VisualizationResult

```json
{
  "type": "object",
  "properties": {
    "success": {"type": "boolean"},
    "message": {"type": "string"},
    "data": {
      "type": "object",
      "properties": {
        "plot_path": {"type": "string"},
        "plot_type": {"type": "string"},
        "image_data": {"type": "string", "format": "base64"}
      }
    }
  }
}
```

## Data Validation

All input parameters are validated using JSON Schema validation. The server will return detailed validation errors if parameters don't match the expected schema.

### Validation Error Format

```json
{
  "success": false,
  "error": {
    "type": "ValidationError",
    "message": "Parameter validation failed",
    "details": [
      {
        "field": "method",
        "error": "Value 'invalid_method' is not allowed",
        "allowed_values": ["spagcn", "stagate", "banksy"]
      }
    ]
  }
}
```

## Dataset Metadata

ChatSpatial tracks metadata for each loaded dataset:

```json
{
  "id": "my_dataset",
  "name": "Mouse Brain Visium",
  "technology": "visium",
  "n_obs": 2696,
  "n_vars": 32285,
  "has_spatial": true,
  "processed": true,
  "annotations": {
    "cell_types": true,
    "spatial_domains": true
  },
  "loaded_at": "2024-01-15T10:30:00Z"
}
```

This metadata is used to validate tool compatibility and provide contextual information to the LLM agent.