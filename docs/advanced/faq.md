# Frequently Asked Questions

Common questions and their answers.

---

## General Questions

### What is ChatSpatial?

ChatSpatial is an LLM agent for spatial transcriptomics analysis via MCP. It integrates 60+ methods from Python and R ecosystems into a unified conversational interface.

### What is MCP?

MCP (Model Context Protocol) is an open standard that enables secure connections between AI applications and external tools. ChatSpatial uses MCP to allow LLMs like Claude to perform spatial analysis through natural language.

### What data formats does ChatSpatial support?

- H5AD (AnnData format)
- 10X Visium spatial folders
- H5 files from 10X Genomics
- MTX matrix files
- Slide-seq, MERFISH, and other spatial formats

### Do I need programming experience?

No! ChatSpatial is designed for natural language conversation. However, understanding basic spatial transcriptomics concepts is helpful.

---

## Installation and Setup

### What Python version do I need?

Python 3.11 or higher (3.12 recommended, 3.13 supported).

### Should I use a virtual environment?

Yes, strongly recommended! Virtual environments prevent dependency conflicts.

### How do I configure ChatSpatial?

See our [Configuration Guide](configuration.md) for MCP client setup instructions.

### Can I use ChatSpatial without Claude?

Yes! ChatSpatial works with any MCP-compatible client.

---

## Analysis Questions

### How do I choose the right method?

See our [Concepts Guide](../concepts.md) for detailed method comparisons. Quick summary:

**Cell Type Annotation:**
| Situation | Recommended Method |
|-----------|-------------------|
| Have reference scRNA-seq | Tangram or scANVI |
| Have marker gene list | CellAssign |
| Want automatic annotation | mLLMCelltype |

**Spatial Domains:**
| Situation | Recommended Method |
|-----------|-------------------|
| Visium with H&E image | SpaGCN |
| High-resolution, no image | STAGATE or GraphST |
| Quick exploration | Leiden |

**Deconvolution:**
| Situation | Recommended Method |
|-----------|-------------------|
| Quick exploration | FlashDeconv (fast) |
| Publication quality | Cell2location (accurate) |
| R environment | RCTD |

### Why does my analysis take so long?

Large datasets and complex methods can be time-consuming. Tips:
- Use GPU acceleration when available (Cell2location, scVI)
- Reduce dataset size for initial testing
- Use faster methods for exploratory analysis
- Check system resources (RAM, CPU)

### How much memory do I need?

Depends on dataset size:
- Small (<5000 cells): 8GB RAM sufficient
- Medium (5000-50000 cells): 16GB RAM recommended
- Large (>50000 cells): 32GB+ RAM needed

---

## Advanced Topics

### Can I use ChatSpatial in research publications?

Yes! ChatSpatial is open-source under MIT license. Please cite our paper if you use it in published research.

### How do I contribute?

See our [Contributing Guide](https://github.com/cafferychen777/ChatSpatial/blob/main/CONTRIBUTING.md).

### Can I add my own analysis methods?

Yes! ChatSpatial's modular architecture makes it easy to add new tools. See developer documentation.

### Is GPU acceleration supported?

Yes, for many methods including Cell2location, scANVI, STAGATE, and VeloVI. Set `use_gpu=True` in parameters to enable.

See [Methods Reference](methods-reference.md) for the full list of GPU-accelerated methods.

---

## Data and Privacy

### Is my data sent to external servers?

No. ChatSpatial runs locally. Your data never leaves your computer. Only natural language commands are sent to the LLM for interpretation.

### Can I use ChatSpatial offline?

Analysis tools work offline, but you need internet to communicate with the LLM (Claude) for command interpretation.

### How is my data stored?

Data is stored in standard H5AD format. You control save locations through the `CHATSPATIAL_DATA_DIR` environment variable.

---

## Still Have Questions?

- Check the [Quick Start](../quickstart.md)
- Read the [Methods Reference](methods-reference.md)
- Browse [GitHub Discussions](https://github.com/cafferychen777/ChatSpatial/discussions)
- [Report issues on GitHub](https://github.com/cafferychen777/ChatSpatial/issues)
