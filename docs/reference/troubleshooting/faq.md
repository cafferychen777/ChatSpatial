# Frequently Asked Questions (FAQ)

Common questions about ChatSpatial and their answers.

## General Questions

### What is ChatSpatial?

ChatSpatial is an agentic workflow orchestration platform for spatial transcriptomics analysis. It integrates 60+ state-of-the-art methods from Python and R ecosystems into a unified conversational interface via the Model Context Protocol (MCP).

### What is MCP?

MCP (Model Context Protocol) is an open standard that enables secure, controlled connections between AI applications and external data sources and tools. ChatSpatial uses MCP to allow LLMs like Claude to perform spatial transcriptomics analysis through natural language.

### What data formats does ChatSpatial support?

ChatSpatial supports:
- H5AD (AnnData format)
- 10X Visium spatial folders
- H5 files from 10X Genomics
- MTX matrix files
- Slide-seq, MERFISH, and other spatial transcriptomics formats

### Do I need programming experience?

No! ChatSpatial is designed to be used through natural language conversation. However, understanding basic spatial transcriptomics concepts is helpful.

## Installation and Setup

### What Python version do I need?

ChatSpatial requires Python 3.10 or higher. Python 3.11 is recommended for optimal performance.

### Should I use a virtual environment?

Yes, strongly recommended! Virtual environments prevent dependency conflicts and make it easier to manage installations.

### How do I configure ChatSpatial with Claude Desktop?

See our [Configuration Guide](../../getting-started/configuration.md) for detailed instructions on setting up the MCP connection.

### Can I use ChatSpatial without Claude Desktop?

Yes! ChatSpatial can work with any MCP-compatible client. You can also use it programmatically via Python.

## Analysis Questions

### Which cell type annotation method should I use?

It depends on your data and needs:
- **Marker-based (CellAssign)**: When you have known marker genes
- **Reference-based (Tangram, scANVI)**: When you have single-cell reference data
- **Machine learning (MLLMCellType)**: For automated annotation

See [Cell Type Annotation Tutorial](../../tutorials/analysis/cell_type_annotation.md) for details.

### How do I choose the right spatial domain method?

Consider your data type:
- **SpaGCN**: Best for Visium with histology images
- **STAGATE**: Good for high-resolution data without images
- **GraphST**: Advanced method with graph self-supervised learning
- **Leiden/Louvain**: Simple clustering-based approach

### Why does my analysis take so long?

Large datasets and complex methods can be time-consuming. Tips to improve speed:
- Use GPU acceleration when available (Cell2location, etc.)
- Reduce dataset size for initial testing
- Use faster methods for exploratory analysis
- Check system resources (RAM, CPU usage)

### How much memory do I need?

Memory requirements depend on dataset size:
- Small datasets (<5000 cells): 8GB RAM sufficient
- Medium datasets (5000-50000 cells): 16GB RAM recommended
- Large datasets (>50000 cells): 32GB+ RAM may be needed

## Troubleshooting

### My tool calls are failing. What should I check?

1. Verify MCP server is running
2. Check data is properly loaded and preprocessed
3. Confirm all required parameters are provided
4. Review error messages for specific issues

See [Common Issues](common_issues.md) for detailed troubleshooting.

### How do I update ChatSpatial?

```bash
pip install --upgrade chatspatial
```

For development version:
```bash
cd ChatSpatial
git pull
pip install -e ".[full]"
```

### What if I encounter a bug?

Please [report it on GitHub](https://github.com/cafferychen777/ChatSpatial/issues) with:
- Error message
- Steps to reproduce
- System information
- Dataset details (if relevant)

## Advanced Topics

### Can I use ChatSpatial in my research publication?

Yes! ChatSpatial is open-source under the MIT license. Please cite our paper if you use it in published research.

### How do I contribute to ChatSpatial?

We welcome contributions! See our [Contributing Guide](https://github.com/cafferychen777/ChatSpatial/blob/main/CONTRIBUTING.md) for details.

### Can I add my own analysis methods?

Yes! ChatSpatial's modular architecture makes it easy to add new tools. See our developer documentation for guidance.

### Is GPU acceleration supported?

Yes, for certain methods:
- Cell2location (deconvolution)
- scVI-based methods
- VeloVI (RNA velocity)

Set `use_gpu=True` in parameters to enable.

## Data and Privacy

### Is my data sent to external servers?

No. ChatSpatial runs locally on your machine. Your data never leaves your computer. The only external communication is with the LLM (e.g., Claude) for interpreting your commands, but actual data and results stay local.

### Can I use ChatSpatial offline?

The analysis tools work offline, but you need an internet connection to communicate with the LLM (Claude) for natural language interpretation.

### How is my data stored?

Data is stored in standard H5AD format. You control where files are saved through the `CHATSPATIAL_DATA_DIR` environment variable.

## Still Have Questions?

- Check the [full documentation](../../index.md)
- Browse [tutorials](../../tutorials/index.md)
- Search [GitHub discussions](https://github.com/cafferychen777/ChatSpatial/discussions)
- Ask the community
