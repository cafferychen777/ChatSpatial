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

Python 3.10 or higher. Python 3.11-3.12 recommended for optimal performance.

### Should I use a virtual environment?

Yes, strongly recommended! Virtual environments prevent dependency conflicts.

### How do I configure ChatSpatial?

See our [Configuration Guide](configuration.md) for MCP client setup instructions.

### Can I use ChatSpatial without Claude?

Yes! ChatSpatial works with any MCP-compatible client.

---

## Analysis Questions

### Which cell type annotation method should I use?

Depends on your data:
- **CellAssign**: When you have known marker genes
- **Tangram/scANVI**: When you have single-cell reference data
- **MLLMCellType**: For automated annotation

### How do I choose the right spatial domain method?

Consider your data type:
- **SpaGCN**: Best for Visium with histology images
- **STAGATE**: Good for high-resolution data without images
- **GraphST**: Advanced graph-based learning
- **Leiden/Louvain**: Simple clustering-based approach

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

## Troubleshooting

### Tool calls are failing. What should I check?

1. Verify data is loaded and preprocessed
2. Check all required parameters are provided
3. Review error messages for specifics
4. See [Troubleshooting Guide](troubleshooting.md)

### How do I update ChatSpatial?

```bash
pip install --upgrade chatspatial
# or for full version:
pip install --upgrade chatspatial[full]
```

### What if I encounter a bug?

[Report it on GitHub](https://github.com/cafferychen777/ChatSpatial/issues) with:
- Error message
- Steps to reproduce
- System information
- Dataset details

---

## Advanced Topics

### Can I use ChatSpatial in research publications?

Yes! ChatSpatial is open-source under MIT license. Please cite our paper if you use it in published research.

### How do I contribute?

See our [Contributing Guide](https://github.com/cafferychen777/ChatSpatial/blob/main/CONTRIBUTING.md).

### Can I add my own analysis methods?

Yes! ChatSpatial's modular architecture makes it easy to add new tools. See developer documentation.

### Is GPU acceleration supported?

Yes, for certain methods:
- Cell2location (deconvolution)
- scVI-based methods (preprocessing, annotation)
- VeloVI (RNA velocity)

Set `use_gpu=True` in parameters to enable.

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
