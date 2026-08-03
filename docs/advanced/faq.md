# Frequently Asked Questions

This page gives **short answers and pointers** to the canonical docs.

---

## General

### What is ChatSpatial?

ChatSpatial is an MCP server for spatial transcriptomics analysis through natural language. See the project overview on the [documentation home](../index.rst).

### What is MCP?

MCP (Model Context Protocol) is the interface layer that lets AI clients call external tools safely. ChatSpatial uses MCP so clients can run spatial analysis tools through structured commands.

### What data formats does ChatSpatial support?

Common formats include H5AD, 10X Visium folders, 10X H5, MTX, Slide-seq, and MERFISH-style spatial data. For a first successful run, see [Quick Start](../quickstart.md).

### Do I need programming experience?

No. ChatSpatial is designed for natural-language use, though basic spatial transcriptomics knowledge still helps.

---

## Setup

### What Python version do I need?

Python 3.11-3.14 is supported, with 3.12 recommended. See [Installation](../installation.md).

### Should I use a virtual environment?

Yes. Use a dedicated environment to avoid dependency conflicts. See [Installation](../installation.md).

### How do I configure ChatSpatial in my client?

Use the [Configuration Guide](configuration.md) for exact client syntax.

### Which clients can I use?

ChatSpatial works with any MCP-compatible client, including Claude Code, Claude
Desktop, Codex, OpenCode, and other MCP-capable tools.

---

## Analysis

### How do I choose the right method?

Use:
- [Concepts](../concepts.md) for method selection guidance
- [Examples](../examples.md) for prompt recipes
- [Methods Reference](methods-reference.md) for exact parameters and defaults

### Why does my analysis take so long?

Large datasets and heavier methods can be slow. Try:
- smaller datasets for testing
- faster baseline methods first
- GPU where supported
- checking memory and CPU limits

If performance is failing rather than merely slow, see [Troubleshooting](troubleshooting.md).

### How much memory do I need?

A rough rule:
- small datasets: 8GB RAM
- medium datasets: 16GB RAM
- large datasets: 32GB+ RAM

See [Troubleshooting](troubleshooting.md) for failure-mode guidance.

---

## Advanced

### Can I use ChatSpatial in publications?

Yes. It is MIT licensed. See citation details in the [documentation home](../index.rst).

### How do I contribute?

See [Contributing](../contributing.md).

### Can I add my own analysis methods?

Yes. The project is modular, but contributor-facing implementation guidance belongs in the codebase and contributor docs rather than the FAQ.

### Is GPU acceleration supported?

Yes, for many methods. See [Methods Reference](methods-reference.md) for exact support.

---

## Data and Privacy

### Is my data sent to external servers?

No. ChatSpatial runs locally; your analysis data stays on your machine.

### Can I use ChatSpatial offline?

The analysis stack can run locally, but you still need an LLM client for natural-language interaction.

### How is my data stored?

Loaded datasets live in the running MCP server memory. Export/reload defaults,
Docker mount behavior, and figure output locations are described in the
[Configuration Guide](configuration.md), with Docker-specific details in
[Docker / GHCR](../docker.md).

---

## Still Have Questions?

- [Installation](../installation.md)
- [Quick Start](../quickstart.md)
- [Configuration Guide](configuration.md)
- [Troubleshooting](troubleshooting.md)
- [Methods Reference](methods-reference.md)
- [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
