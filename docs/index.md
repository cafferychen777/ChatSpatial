# ChatSpatial: Agentic Spatial Transcriptomics Analysis

![ChatSpatial](https://img.shields.io/badge/ChatSpatial-MCP%20Server-blue?style=for-the-badge&logo=python)

**Interactive Spatial Transcriptomics Analysis via Model Context Protocol**

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://python.org)
[![MCP](https://img.shields.io/badge/MCP-Compatible-green.svg)](https://modelcontextprotocol.io)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/cafferychen777/ChatSpatial/blob/main/LICENSE)
[![GitHub](https://img.shields.io/badge/GitHub-Repository-black.svg?logo=github)](https://github.com/cafferychen777/ChatSpatial)

## 🤖 What is Model Context Protocol (MCP)?

**Model Context Protocol (MCP)** is an open standard that enables LLM agents to securely connect with external tools and data sources. MCP functions as a universal plugin system for agentic systems that enables natural language interaction with specialized software.

**📚 [Learn more about MCP →](resources/what_is_mcp.md)**

### 🎥 New to MCP? Watch This!

[![MCP Demo Video](https://img.shields.io/badge/▶️_Watch-MCP_Demo_Video-red?style=for-the-badge&logo=youtube)](https://www.youtube.com/watch?v=sfCBCyNyw7U)

*Simple 5-minute explanation of what MCP is and why it matters*

### Why MCP + ChatSpatial?

Rather than learning complex bioinformatics tools, users can ask:

> *"Load my Visium data and identify spatial domains"*
> *"Which genes show spatial patterns in my dataset?"*
> *"Visualize cell communication networks"*

ChatSpatial handles the technical complexity while you focus on biological insights.

## 🎯 What is ChatSpatial?

ChatSpatial is a production-ready MCP server that provides LLM agents with comprehensive spatial transcriptomics analysis capabilities. It enables natural language interaction with complex spatial data analysis through 16 standardized MCP tools.

### Key Features

- **16+ Analysis Tools**: Complete spatial transcriptomics workflow from preprocessing to visualization
- **Agent-Native**: Designed for seamless integration with Claude, GPT, and other LLM agents
- **Multiple Technologies**: Support for Visium, MERFISH, Slide-seq, and other spatial platforms
- **Rich Visualizations**: 20 plot types with MCP image objects for direct agent display
- **Production Ready**: Robust error handling, validation, and comprehensive testing

## 🚀 Quick Start

### 1. Installation

```bash
# Create environment
conda create -n chatspatial python=3.10
conda activate chatspatial

# Install ChatSpatial
pip install -e .
```

### 2. MCP Configuration

Add to your Claude Desktop configuration:

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

### 3. First Analysis

```python
# Load spatial data
result = load_data(data_path="data.h5ad", name="my_data")

# Preprocess
preprocess_data(data_id=result.id)

# Identify spatial domains
identify_spatial_domains(data_id=result.id, method="spagcn")

# Visualize results
visualize_data(data_id=result.id, plot_type="spatial_domains")
```

## 📚 Documentation Structure

### 🏁 Getting Started

- **[Getting Started Guide](getting-started/)** - Complete setup and installation
- **[Installation](getting-started/installation.md)** - Detailed installation instructions
- **[Quick Start](getting-started/quick-start.md)** - Get up and running quickly

### 📖 Tutorials

Structured learning paths from beginner to advanced:

- **[Core Tutorials](tutorials/core/)** - Essential concepts and basic workflows
  - [Basic Spatial Analysis](tutorials/core/basic_spatial_analysis.md)
  - [Spatial Statistics](tutorials/core/spatial_statistics.md)
  - [Visualization Tutorial](tutorials/core/visualization_tutorial.md)
- **[Analysis Tutorials](tutorials/analysis/)** - Specialized analysis methods
  - [Cell Type Annotation](tutorials/analysis/cell_type_annotation.md)
  - [Cell Communication Analysis](tutorials/analysis/cell_communication_analysis.md)
  - [Spatial Enrichment](tutorials/analysis/spatial_enrichment.md)
- **[Advanced Tutorials](tutorials/advanced/)** - Complex multi-modal analysis
  - [Batch Integration](tutorials/advanced/batch_integration.md)
  - [Spatial Registration](tutorials/advanced/spatial_registration.md)
  - [Trajectory Analysis](tutorials/advanced/trajectory_analysis.md)
- **[Learning Paths](tutorials/learning-paths/)** - Guided skill development tracks

### 🔧 Reference

Complete reference materials:

- **[API Reference](reference/api/)** - Complete MCP tool documentation
- **[Quick Reference](reference/quick-reference/)** - Cheat sheets and quick guides
- **[Configuration](reference/configuration.md)** - Settings and configuration options
- **[Data Formats](reference/data_formats.md)** - Supported formats and schemas
- **[Troubleshooting](reference/troubleshooting/)** - Common issues and solutions
- **[Performance](reference/performance.md)** - Optimization guidelines
- **[Benchmarks](reference/benchmarks.md)** - Performance comparisons

### 🎯 Examples & Use Cases

Real-world examples and workflows:

- **[Example Workflows](examples/workflows/)** - Complete analysis examples
- **[Datasets](examples/datasets/)** - Sample datasets and usage patterns
- **[Integration Guides](examples/guides/)** - Setup and integration examples

## 🛠️ Core Analysis Tools

| Category | Tools | Description |
|----------|-------|-------------|
| **Data Management** | `load_data`, `preprocess_data` | Data loading and preprocessing |
| **Cell Annotation** | `annotate_cells` | 7 annotation methods (Tangram, scType, etc.) |
| **Spatial Analysis** | `analyze_spatial_data`, `identify_spatial_domains` | Pattern analysis and domain identification |
| **Cell Communication** | `analyze_cell_communication` | LIANA+, CellPhoneDB, CellChat |
| **Deconvolution** | `deconvolve_data` | Cell2location, NNLS methods |
| **Gene Analysis** | `identify_spatial_genes`, `find_markers` | Spatial variable genes and markers |
| **Visualization** | `visualize_data` | 20 plot types with agent-friendly outputs |

## 🏗️ MCP Architecture

ChatSpatial implements the Model Context Protocol for seamless agentic integration:

```mermaid
graph LR
    A[LLM Agent] --> B[MCP Client]
    B --> C[ChatSpatial Server]
    C --> D[Analysis Tools]
    D --> E[Results]
    E --> F[Visualizations]
    F --> B
    B --> A
```

### Why MCP?

- **Standardized Interface**: Consistent tool discovery and invocation
- **Type Safety**: JSON Schema validation for all inputs/outputs
- **Error Handling**: Structured error reporting and recovery
- **Streaming Support**: Real-time progress updates
- **Security**: Built-in access controls and validation

## 🔬 Supported Technologies

| Platform | Status | Features |
|----------|--------|----------|
| **10x Visium** | ✅ Full Support | Spatial domains, deconvolution, communication |
| **MERFISH** | ✅ Full Support | High-resolution analysis, trajectory inference |
| **Slide-seq** | ✅ Full Support | Subcellular resolution analysis |
| **STARmap** | ✅ Supported | 3D spatial analysis |
| **seqFISH+** | ✅ Supported | Single-cell resolution |
| **Xenium** | 🔄 In Progress | Latest spatial technology |

## 📊 Analysis Capabilities

### Preprocessing & QC

- Normalization (log, SCTransform, Pearson residuals)
- Quality control metrics
- Batch effect correction
- Spatial coordinate validation

### Cell Type Annotation

- **Marker-based**: Traditional marker gene approach
- **Tangram**: Spatial mapping with reference data
- **scType**: Automated cell type identification
- **Cell2location**: Probabilistic deconvolution
- **scANVI**: Semi-supervised annotation
- **CellAssign**: Probabilistic assignment
- **MLLMCellType**: Multi-modal LLM classifier

### Spatial Domain Identification

- **SpaGCN**: Graph convolutional networks
- **STAGATE**: Spatial-temporal attention
- **BANKSY**: Spatial clustering
- **Leiden/Louvain**: Community detection

### Cell Communication Analysis

- **LIANA+**: Comprehensive ligand-receptor analysis
- **CellPhoneDB**: Statistical interaction testing
- **CellChat**: Systematic communication analysis

## 🎨 Visualization Gallery

ChatSpatial provides rich visualizations optimized for LLM agents:

- **Spatial Plots**: Gene expression, domains, communication
- **UMAP/t-SNE**: Dimensionality reduction with annotations
- **Heatmaps**: Expression patterns and correlations
- **Violin Plots**: Distribution comparisons
- **Trajectory Plots**: Developmental pathways
- **Communication Networks**: Cell-cell interaction graphs

All visualizations are returned as MCP Image objects for direct display in agent interfaces.

## 🤝 Contributing

We welcome contributions! See our [Contributing Guide](../CONTRIBUTING.md) for details.

## 📄 License

ChatSpatial is released under the MIT License. See [LICENSE](../LICENSE) for details.

## 🆘 Support

- **Documentation**: Browse this site for comprehensive guides
- **Issues**: Report bugs on [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
- **Discussions**: Join our [GitHub Discussions](https://github.com/cafferychen777/ChatSpatial/discussions)

