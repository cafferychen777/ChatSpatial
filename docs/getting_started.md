# Getting Started

This guide helps you install and run ChatSpatial quickly.

## Installation

Follow the minimal steps:

```bash
conda create -n chatspatial python=3.10
conda activate chatspatial
pip install -e .
chatspatial --help
```

For more detailed instructions, see [INSTALLATION.md](../INSTALLATION.md).

## Minimal MCP Setup

Add to your MCP client configuration:

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/your/python",
      "args": ["-m", "chatspatial"],
      "env": {}
    }
  }
}
```

## First Task

- Load a Visium dataset (H5AD)
- Run preprocess_data
- Run identify_spatial_domains (SpaGCN)
- Visualize with visualize_data

