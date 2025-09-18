---
layout: default
title: Home
nav_order: 1
description: "ChatSpatial is an MCP server for spatial transcriptomics analysis through natural language interaction with Claude and other LLMs."
permalink: /
---

# ChatSpatial Documentation
{: .fs-9 }

Spatial transcriptomics analysis through natural language interaction
{: .fs-6 .fw-300 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View it on GitHub](https://github.com/cafferychen777/ChatSpatial){: .btn .fs-5 .mb-4 .mb-md-0 }

---

{: .warning }
> This documentation is currently being migrated to Just the Docs. Some links may be broken during the transition.

## What is ChatSpatial?

ChatSpatial is a Model Context Protocol (MCP) server that enables spatial transcriptomics analysis through natural language interaction. It integrates with Claude and other LLMs to perform complex spatial genomics analysis through simple conversational commands.

### Key Features

- **Natural Language Interface**: Analyze spatial data using conversational commands
- **Comprehensive Analysis Tools**: Cell type annotation, spatial domains, trajectory analysis, and more
- **Multiple Data Formats**: Support for H5AD, 10X Visium, Slide-seq, and other spatial formats
- **Advanced Algorithms**: Integration with state-of-the-art spatial analysis methods
- **Visualization**: Built-in plotting and visualization capabilities

### Supported Analysis Types

- **Cell Type Annotation**: Marker-based, reference-based (Tangram, scANVI), and ML approaches
- **Spatial Domain Identification**: SpaGCN, STAGATE, BANKSY, and clustering methods
- **Cell Communication**: LIANA, CellPhoneDB, and spatial interaction analysis
- **Trajectory Analysis**: RNA velocity, pseudotime, and developmental trajectories
- **Spatial Statistics**: Moran's I, spatial autocorrelation, and neighborhood analysis
- **Deconvolution**: Cell2location, RCTD, Stereoscope for cell type proportions

## Getting Started

### Prerequisites

- Python 3.8 or higher
- Claude Desktop or compatible MCP client
- Git for installation

### Quick Installation

```bash
# Clone the repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial

# Install dependencies
pip install -e .

# Run the MCP server
python -m chatspatial
```

For detailed installation instructions, see the [Installation Guide](getting-started/installation.html).

### Your First Analysis

Once installed, you can start analyzing spatial data immediately:

```
Load the mouse brain Visium dataset and show me a spatial plot of the top variable genes
```

```
Identify spatial domains in my data using SpaGCN and visualize the results
```

```
Perform cell type annotation using marker genes and show the spatial distribution
```

## Documentation Structure

<div class="code-example" markdown="1">

**Getting Started**
- [Installation](getting-started/installation.html)
- [Quick Start Guide](getting-started/quick-start.html)
- [Basic Configuration](getting-started/configuration.html)

**Tutorials**
- [Core Analysis](tutorials/core/)
- [Advanced Methods](tutorials/advanced/)
- [Learning Paths](tutorials/learning-paths/)

**Reference**
- [API Documentation](reference/api/)
- [Data Formats](reference/data-formats.html)
- [Troubleshooting](reference/troubleshooting/)

**Examples**
- [Example Workflows](examples/)
- [Case Studies](examples/workflows/case_studies/)

</div>

## Community and Support

- **GitHub Issues**: [Report bugs and request features](https://github.com/cafferychen777/ChatSpatial/issues)
- **Discussions**: [Community Q&A and discussions](https://github.com/cafferychen777/ChatSpatial/discussions)
- **Documentation**: This site and inline help

## License

ChatSpatial is distributed under the MIT License. See [LICENSE](https://github.com/cafferychen777/ChatSpatial/blob/main/LICENSE) for more information.

---

### About This Documentation

This documentation is built with [Just the Docs](https://just-the-docs.github.io/just-the-docs/), a documentation theme for Jekyll.