# ChatSpatial Documentation

**Chat with your spatial transcriptomics data. No coding required.**

ChatSpatial is a Model Context Protocol (MCP) server that integrates 60+ spatial analysis methods into Claude. Analyze your data through natural language conversation.

---

## 🚀 Getting Started (Choose Your Path)

<table>
<tr>
<td width="50%">

### 🎯 **New Users Start Here**

**I want to analyze data quickly:**

1. [**Quick Start**](quickstart.md) - 5-minute setup
2. [**Examples**](examples.md) - Copy-paste workflows
3. Start chatting with your data!

✅ Perfect for: Researchers, biologists, anyone new to ChatSpatial

</td>
<td width="50%">

### 🔧 **Advanced Users**

**I need detailed documentation:**

- [**Methods Reference**](advanced/methods-reference.md) - All tools and parameters
- [**Installation Guide**](advanced/installation.md) - Detailed setup options
- [**Configuration**](advanced/configuration.md) - Advanced settings
- [**Troubleshooting**](advanced/troubleshooting.md) - Problem solving
- [**FAQ**](advanced/faq.md) - Common questions

✅ Perfect for: Power users, developers, troubleshooting

</td>
</tr>
</table>

---

## 📖 Documentation Structure

This documentation is organized for **quick access**:

```
docs/
├── quickstart.md              # Start here! (5 minutes)
├── examples.md                # Real analysis examples
├── advanced/                  # Detailed documentation
│   ├── methods-reference.md  # Complete tool reference
│   ├── installation.md       # Detailed setup guide
│   ├── configuration.md      # Advanced configuration
│   ├── troubleshooting.md    # Problem solving
│   └── faq.md                # Frequently asked questions
└── examples/                  # Sample datasets and workflows
```

**80% of users only need:** `quickstart.md` + `examples.md`

---

## 💡 Quick Reference

### What Can ChatSpatial Do?

| Analysis Type | Natural Language Example |
|--------------|--------------------------|
| **Load Data** | "Load my Visium dataset" |
| **Spatial Domains** | "Identify tissue regions" |
| **Cell Types** | "Annotate cell types using reference" |
| **Deconvolution** | "Deconvolve spots with Cell2location" |
| **Communication** | "Analyze cell-cell interactions" |
| **Trajectories** | "Find developmental paths" |
| **Enrichment** | "Run pathway analysis" |
| **Visualization** | "Create spatial heatmap" |

**60+ methods across 12 analysis categories** - all through conversation!

---

## 🎓 Learning Paths

### Path 1: First-Time User (30 minutes)
1. [Quick Start](quickstart.md) - Setup and first analysis
2. [Examples - Basic Analysis](examples.md#-basic-spatial-analysis) - Learn fundamentals
3. [Examples - Spatial Domains](examples.md#%EF%B8%8F-spatial-domain-analysis) - Identify tissue regions
4. Try your own data!

### Path 2: Intermediate User (1 hour)
1. [Examples - Cell Type Analysis](examples.md#-cell-type-analysis) - Annotation methods
2. [Examples - Deconvolution](examples.md#-deconvolution-analysis) - Estimate compositions
3. [Examples - Communication](examples.md#-cell-communication-analysis) - Interactions
4. [Complete Workflows](examples.md#-complete-workflows) - Multi-step analysis

### Path 3: Advanced User (2+ hours)
1. [Examples - Advanced Analysis](examples.md#-advanced-analysis) - Trajectories, CNV, etc.
2. [Methods Reference](advanced/methods-reference.md) - Deep dive into all tools
3. [Configuration](advanced/configuration.md) - Optimize for your use case
4. Combine multiple methods for publication-quality analysis

---

## 🆘 Getting Help

### Quick Troubleshooting

**Problem: ChatSpatial not showing up in Claude**
- ✅ Restart Claude after configuration
- ✅ Check Python path is correct
- ✅ See [Troubleshooting Guide](advanced/troubleshooting.md)

**Problem: Analysis failing**
- ✅ Run preprocessing first
- ✅ Use absolute file paths
- ✅ Check data format
- ✅ See [FAQ](advanced/faq.md)

**Problem: Don't know which method to use**
- ✅ Just ask Claude naturally!
- ✅ Claude will pick the right tool
- ✅ See [Examples](examples.md) for ideas

### Need More Help?

- **Search**: Use the search box (top of page)
- **Issues**: [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
- **Discussions**: [GitHub Discussions](https://github.com/cafferychen777/ChatSpatial/discussions)

---

## 🎯 Common Questions

**Q: Do I need to know Python?**
A: No! ChatSpatial understands natural language. Just chat with Claude.

**Q: What data formats are supported?**
A: H5AD, 10x Visium, Slide-seq, MERFISH, seqFISH, and more. [Details →](advanced/methods-reference.md)

**Q: How do I choose between methods?**
A: Just describe what you want. Claude and ChatSpatial will pick the right tool.

**Q: Can I use my own reference data?**
A: Yes! Load your reference and ask ChatSpatial to use it for annotation or deconvolution.

**Q: Is ChatSpatial free?**
A: Yes! MIT license. Free for academic and commercial use.

[**More Questions →**](advanced/faq.md)

---

## 📊 Example Conversation

```text
👤 "Load /path/to/visium_data.h5ad and identify spatial domains"

🤖 ✅ Loaded 3,456 spots, 18,078 genes
    ✅ Identified 7 spatial domains using SpaGCN
    ✅ Generated visualization

👤 "Find marker genes for domain 3 and show me what cell type it is"

🤖 ✅ Found 23 significant markers (adj. p < 0.05)
    Top markers: GFAP, S100B, AQP4
    ✅ Domain 3 shows astrocyte signature
    ✅ Created expression heatmap

👤 "Now analyze cell communication between all domains"

🤖 ✅ Analyzed 142 ligand-receptor pairs using LIANA
    ✅ Found 18 significant interactions (p < 0.05)
    ✅ Generated communication network diagram
```

**That's ChatSpatial. Natural conversation → Scientific results.** 🎉

---

## 🔗 Quick Links

**Essential:**
- [Quick Start](quickstart.md) - Get running in 5 minutes
- [Examples](examples.md) - Real analysis workflows
- [Methods Reference](advanced/methods-reference.md) - All tools

**Resources:**
- [Main README](../README.md) - Project overview
- [GitHub Repository](https://github.com/cafferychen777/ChatSpatial)
- [Sample Datasets](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data)

**Support:**
- [Troubleshooting](advanced/troubleshooting.md)
- [FAQ](advanced/faq.md)
- [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)

---

<div align="center">

**Ready to start?** [→ Quick Start Guide](quickstart.md)

Made with ❤️ for the spatial transcriptomics community

</div>
