<div align="center">

# ChatSpatial

**MCP server for spatial transcriptomics analysis via natural language**

[![Paper](https://img.shields.io/badge/bioRxiv-2026.02.26.708361-b31b1b.svg)](https://doi.org/10.64898/2026.02.26.708361)
[![MLGenX @ ICLR 2026](https://img.shields.io/badge/MLGenX%20@%20ICLR%202026-Accepted-brightgreen.svg)](https://openreview.net/forum?id=xZ814yNaUW)
[![ENAR 2026](https://img.shields.io/badge/ENAR%202026-Oral-blue.svg)](https://www.enar.org/meetings/spring2026/)
[![IBC 2026](https://img.shields.io/badge/IBC%202026-Oral-blue.svg)](https://www.ibc2026.org/home)
[![CI](https://github.com/cafferychen777/ChatSpatial/actions/workflows/ci.yml/badge.svg)](https://github.com/cafferychen777/ChatSpatial/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chatspatial)](https://pypi.org/project/chatspatial/)
[![Python 3.11-3.13](https://img.shields.io/badge/python-3.11--3.13-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docs](https://img.shields.io/badge/docs-available-blue)](https://cafferychen777.github.io/ChatSpatial/)

</div>

<p align="center">
  <img src="assets/images/overview.jpg" alt="ChatSpatial Overview" width="900">
</p>

ChatSpatial replaces ad-hoc LLM code generation with **schema-enforced orchestration**. Instead of generating arbitrary scripts, the LLM selects tools and parameters from a curated registry, making spatial transcriptomics workflows more reproducible across sessions and clients.

It exposes **60+ spatial transcriptomics methods** as MCP tools, so any MCP-compatible client can analyze data through natural language.

---

## Start Here

1. **Install ChatSpatial** — [Installation Guide](INSTALLATION.md)
2. **Configure your MCP client** — [Configuration Guide](docs/advanced/configuration.md)
3. **Run your first analysis** — [Quick Start](docs/quickstart.md)

**Minimal example prompt:**

```text
Load /absolute/path/to/spatial_data.h5ad and show me the tissue structure
```

> ChatSpatial works with **any MCP-compatible client** — Claude Code, Claude Desktop, Codex, OpenCode, and other MCP-capable tools.

---

## Capabilities

60+ methods across 11 categories. Supports 10x Visium, Xenium, Slide-seq v2, MERFISH, seqFISH.

| Category | Methods |
|----------|---------|
| **Spatial Domains** | SpaGCN, STAGATE, GraphST, BANKSY, Leiden, Louvain |
| **Deconvolution** | FlashDeconv, Cell2location, RCTD, DestVI, Stereoscope, SPOTlight, Tangram, CARD |
| **Cell Communication** | LIANA+, CellPhoneDB, CellChat (`cellchat_r`), FastCCC |
| **Cell Type Annotation** | Tangram, scANVI, CellAssign, mLLMCelltype, scType, SingleR |
| **Differential Expression** | Wilcoxon, t-test, Logistic Regression, pyDESeq2 |
| **Trajectory & Velocity** | CellRank, Palantir, DPT, scVelo, VeloVI |
| **Spatial Statistics** | Moran's I, Local Moran, Geary's C, Getis-Ord Gi*, Ripley's K, Co-occurrence, Neighborhood Enrichment, Centrality Scores, Local Join Count, Network Properties |
| **Enrichment** | GSEA, ORA, Enrichr, ssGSEA, Spatial EnrichMap |
| **Spatial Genes** | SpatialDE, SPARK-X, FlashS |
| **Integration** | Harmony, BBKNN, Scanorama, scVI |
| **Other** | CNV Analysis (InferCNVPy, Numbat), Spatial Registration (PASTE, STalign) |

---

## Documentation

| Guide | Owns |
|-------|------|
| [Installation](INSTALLATION.md) | Environment setup, package install, platform notes |
| [Quick Start](docs/quickstart.md) | First successful analysis after setup |
| [Concepts](docs/concepts.md) | Method selection and analysis reasoning |
| [Examples](docs/examples.md) | Prompt recipes and workflow examples |
| [Configuration](docs/advanced/configuration.md) | Exact MCP client configuration syntax |
| [Troubleshooting](docs/advanced/troubleshooting.md) | Symptom → fix guidance |
| [Methods Reference](docs/advanced/methods-reference.md) | Canonical tool parameters and defaults |
| [Full Docs](https://cafferychen777.github.io/ChatSpatial/) | Complete documentation site |

---

## Citation

If you use ChatSpatial in your research, please cite:

```bibtex
@article{Yang2026.02.26.708361,
  author = {Yang, Chen and Zhang, Xianyang and Chen, Jun},
  title = {ChatSpatial: Schema-Enforced Agentic Orchestration for Reproducible and Cross-Platform Spatial Transcriptomics},
  elocation-id = {2026.02.26.708361},
  year = {2026},
  doi = {10.64898/2026.02.26.708361},
  publisher = {Cold Spring Harbor Laboratory},
  URL = {https://www.biorxiv.org/content/early/2026/03/01/2026.02.26.708361},
  journal = {bioRxiv}
}
```

ChatSpatial orchestrates many excellent third-party methods. **Please also cite the original tools your analysis used.**

---

## Contributing

Documentation improvements, bug reports, and new analysis methods are all welcome. See [CONTRIBUTING.md](CONTRIBUTING.md).

<div align="center">

**MIT License** · [GitHub](https://github.com/cafferychen777/ChatSpatial) · [Issues](https://github.com/cafferychen777/ChatSpatial/issues)

</div>

<!-- mcp-name: io.github.cafferychen777/chatspatial -->

## Hosted deployment

A hosted deployment is available on [Fronteir AI](https://fronteir.ai/mcp/cafferychen777-chatspatial).

