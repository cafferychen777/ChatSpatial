
title: Spatial Statistics
description: Statistical analysis of spatial patterns and autocorrelation
---

# Spatial Statistics Methods with ChatSpatial

Identify spatial patterns in your tissue using detailed statistical methods through natural conversation. This tutorial shows you how to perform spatial statistics analysis by asking ChatSpatial what you want to explore - no statistical background required.

## What You'll Learn

By the end of this tutorial, you'll know how to:

1. **Perform Bivariate Moran's I Analysis** - discover spatial correlations between gene pairs
2. **Run Local Indicators of Spatial Association (LISA)** - identify hotspots and coldspots
3. **Execute Neighborhood Enrichment Analysis** - understand cell type co-occurrence patterns
4. **Apply Ripley's K Function** - analyze point pattern distributions
5. **Calculate Spatial Centrality Metrics** - find network hubs in your tissue
6. **Conduct Join Count Analysis** - measure categorical spatial autocorrelation

For each method, you'll learn:
- What biological questions it answers
- How to ask ChatSpatial to run it
- When to use each approach
- How to interpret the results

## Before You Start

### What You Need
- ✅ **Claude Desktop** installed and running
- ✅ **ChatSpatial configured** (see [Getting Started](../../getting-started/index.md))
- ✅ **Your spatial data** already loaded and preprocessed
- ✅ **Basic spatial analysis completed** (see [Basic Spatial Analysis](basic_spatial_analysis.md))
- ✅ **Interest in spatial patterns** in your tissue

### Understanding Spatial Statistics

**Traditional statistics** asks: "Is this gene differentially expressed?"

**Spatial statistics** asks: "Is this biological pattern organized spatially in my tissue?"

This approach identifies:
- 🗺️ **Spatial hotspots** - regions of concentrated activity
- 🔗 **Spatial correlations** - how genes co-vary across space
- 📍 **Point patterns** - how cells are distributed spatially
- 🌐 **Network structures** - connectivity patterns in tissue architecture

## Method 1: Bivariate Moran's I Analysis

### Understanding Bivariate Spatial Correlation

**What it does:** Measures whether two genes show similar spatial patterns - are they both high or low in the same regions?

**Biological applications:**

- Co-regulated gene pairs in development
- Metabolic pathway coordination across tissue
- Stress response gene co-localization
- Transcription factor-target relationships

### 💬 Real Conversation Example

**👤 You:** "I have two genes, CREB1 and BDNF, that I think might be co-regulated in memory formation. Can you check if they show similar spatial patterns in my mouse brain data?"

**🤖 Claude:** "Good hypothesis! CREB1 and BDNF are indeed key memory formation genes. Let me perform bivariate Moran's I analysis to see if they show coordinated spatial patterns in your brain tissue."

*[ChatSpatial analyzes bivariate spatial correlation]*

**🤖 Claude:** "Results show that your memory genes display strong spatial co-regulation:

**🧠 Bivariate Moran's I Results for CREB1 vs BDNF:**

**Global Bivariate Correlation:**
- **Bivariate Moran's I**: 0.67 (strong positive spatial correlation)
- **p-value**: < 0.001 (highly significant)
- **Interpretation**: These genes are spatially co-expressed across your brain tissue

**Local Spatial Co-regulation:**
- **Hippocampus CA1**: Both genes peak together (local I = 0.85)
- **Cortical layers 2-3**: Strong coordinated expression (local I = 0.72)
- **White matter**: Both consistently low (local I = 0.61)

**Biological Interpretation:** CREB1 and BDNF show spatial coordination, consistent with their role as co-regulated memory formation factors. They are highest in brain regions critical for learning and memory."

**👤 You:** "That is informative! Can you check a few more gene pairs involved in memory circuits?"

**🤖 Claude:** "Absolutely! Let me analyze multiple memory-related gene pairs to map your brain's memory circuit architecture."

---

### 🎯 When to Use Bivariate Moran's I

**Used for:**

- 🧬 **Testing co-regulation hypotheses** - do pathway genes coordinate spatially?
- 📊 **Validating gene networks** - are network partners spatially correlated?
- 🔬 **Discovering new interactions** - find unexpected spatial correlations
- 🎯 **Pathway analysis** - understand spatial organization of biological processes

### 💡 Natural Ways to Ask

**Questions you can ask:**

- "Do these two genes show similar spatial patterns?"
- "Are GENE1 and GENE2 spatially correlated?"
- "Can you find genes that co-vary with my favorite gene spatially?"
- "Which gene pairs show the strongest spatial correlation?"
- "Test these pathway genes for spatial coordination"

## Method 2: Local Indicators of Spatial Association (LISA)

### Understanding LISA Analysis

**What it does:** Identifies specific tissue locations where genes show unusual spatial clustering - hotspots where expression is much higher than expected.

**Biological applications:**

- Disease lesion identification
- Stem cell niche detection
- Immune infiltration zones
- Developmental organizing centers

### 💬 Real Conversation Example

**👤 You:** "I want to find hotspots of immune activity in my tumor tissue. Can you identify where immune response genes are spatially clustered?"

**🤖 Claude:** "Good question for understanding tumor immunity! Let me perform LISA analysis to identify immune hotspots in your tumor tissue using local spatial association statistics."

*[ChatSpatial runs LISA analysis on immune genes]*

**🤖 Claude:** "Analysis reveals immune spatial organization:

**🎯 LISA Results for Immune Response Genes:**

**CD68 (Macrophage Marker):**
- **High-High clusters**: 3 distinct hotspots (89 spots total)
- **Hotspot 1**: Tumor edge - dense macrophage infiltration
- **Hotspot 2**: Necrotic core boundary - cleanup macrophages
- **Hotspot 3**: Vascular regions - patrolling macrophages
- **Significance**: p < 0.001 for all hotspots

**IFNG (T-cell Activation):**
- **High-High clusters**: 2 major hotspots (56 spots total)
- **Hotspot 1**: Tumor-stroma interface - active immune response
- **Hotspot 2**: Lymphoid aggregate - T-cell education zone
- **Low-Low regions**: Tumor core - immune exclusion zones

**IL1B (Inflammatory Cytokine):**
- **High-High clusters**: Overlaps with macrophage hotspots
- **Spatial pattern**: Ring-like distribution around tumor areas
- **Biological significance**: Inflammatory barrier formation

**🔍 Key Finding:** Your tumor has 3 distinct immune microenvironments - infiltrated edges, excluded cores, and vascular patrol zones."

**👤 You:** "Can you create a map showing these immune hotspots?"

---

### 🎯 When to Use LISA Analysis

**Used for:**

- 🎯 **Hotspot identification** - where are the action zones?
- 🏥 **Disease analysis** - locate pathological processes
- 🧪 **Quality control** - detect technical artifacts or batch effects
- 🌱 **Development studies** - find organizing centers and gradients

### 💡 Natural Ways to Ask

**Questions you can ask:**

- "Where are the hotspots for this gene?"
- "Can you find spatial clusters of immune activity?"
- "Identify regions where my gene is unusually high"
- "Show me local spatial patterns for this pathway"
- "Find the organizing centers in my developmental data"

## Method 3: Neighborhood Enrichment Analysis

### Understanding Cell Type Co-occurrence

**What it does:** Determines whether specific cell types prefer to be neighbors - are certain cell types more likely to be adjacent than expected by chance?

**Biological applications:**

- Tissue architecture analysis
- Cell-cell communication zones
- Developmental lineage tracking
- Pathological tissue organization

### 💬 Real Conversation Example

**👤 You:** "In my brain organoid data, I want to understand how different cell types organize themselves. Do neurons prefer to be near astrocytes? Can you analyze cell type neighborhoods?"

**🤖 Claude:** "Good question for understanding brain tissue architecture! Let me perform neighborhood enrichment analysis to map how your brain organoid cell types organize spatially."

*[ChatSpatial analyzes cell type neighborhoods]*

**🤖 Claude:** "Analysis reveals cellular architecture patterns:

**🔍 Neighborhood Enrichment Results:**

**Strong Positive Associations (Preferred Neighbors):**
- **Excitatory Neurons ↔ Astrocytes**: Z-score = +4.2 (p < 0.001)
  - *Biology*: Astrocytes provide metabolic support to active neurons
- **Inhibitory Neurons ↔ Excitatory Neurons**: Z-score = +3.8 (p < 0.001)
  - *Biology*: Balanced excitation-inhibition circuits
- **Oligodendrocytes ↔ Neurons**: Z-score = +3.1 (p < 0.01)
  - *Biology*: Myelination of axonal projections

**Neutral Associations (Random Distribution):**
- **Microglia ↔ All Types**: Z-scores near 0
  - *Biology*: Uniform tissue surveillance pattern

**Negative Associations (Spatial Avoidance):**
- **Neural Stem Cells ↔ Mature Neurons**: Z-score = -2.9 (p < 0.01)
  - *Biology*: Stem cells concentrate in neurogenic niches
- **Astrocytes ↔ Astrocytes**: Z-score = -2.1 (p < 0.05)
  - *Biology*: Territorial spacing to cover tissue efficiently

**🎯 Key Finding:** Your brain organoid recapitulates natural tissue architecture - supporting cells cluster with neurons they serve, while maintaining functional territories."

**👤 You:** "Can you show me which regions have the most organized cell type patterns?"

---

### 🎯 When to Use Neighborhood Enrichment

**Used for:**

- 🏗️ **Tissue architecture** - understanding structural organization
- 💬 **Communication analysis** - identifying interaction zones
- 🔬 **Developmental biology** - tracking lineage relationships
- 🏥 **Pathology studies** - how disease disrupts normal organization

### 💡 Natural Ways to Ask

**Questions you can ask:**

- "Which cell types prefer to be neighbors?"
- "How do my cell types organize spatially?"
- "Are there cell-cell communication hotspots?"
- "Do tumor cells avoid immune cells?"
- "Show me cell type territorial patterns"

## Method 4: Ripley's K Function Analysis

### Understanding Point Pattern Analysis

**What it does:** Analyzes whether cells of specific types are randomly distributed, clustered together, or regularly spaced across tissue.

**Biological applications:**

- Stem cell niche organization
- Immune cell patrol patterns
- Blood vessel distribution analysis
- Tumor cell invasion patterns

### 💬 Real Conversation Example

**👤 You:** "I have cancer stem cells marked in my tumor tissue. Are they randomly scattered or do they form organized clusters? Can you analyze their spatial distribution pattern?"

**🤖 Claude:** "Good question about cancer stem cell organization! Ripley's K function analysis will reveal whether your stem cells show clustering, regular spacing, or random distribution - each pattern has different biological implications."

*[ChatSpatial performs Ripley's K analysis]*

**🤖 Claude:** "Analysis reveals cancer stem cell spatial organization:

**📊 Ripley's K Function Results:**

**Cancer Stem Cell Distribution Pattern:**
- **Pattern Type**: Strong clustering (non-random)
- **Clustering Scale**: Most pronounced at 50-150 μm distances
- **Peak Clustering**: At ~80 μm (typical niche size)
- **Significance**: p < 0.001 vs random distribution

**Spatial Characteristics:**
- **Cluster Count**: 7 distinct stem cell clusters identified
- **Average Cluster Size**: 12.3 cells per cluster
- **Inter-cluster Distance**: 280 ± 45 μm (regular spacing)
- **Pattern Interpretation**: Territorial niche organization

**Cross-type Analysis (Stem Cells vs Other Cancer Cells):**
- **Attraction Zone**: 0-60 μm (stem cells recruit other cancer cells)
- **Neutral Zone**: 60-120 μm (random association)
- **Exclusion Zone**: >120 μm (stem cell territorial boundaries)

**🔬 Biological Interpretation:** Your cancer stem cells form organized niches with territorial spacing - they create supportive microenvironments while maintaining optimal resource access through regular spacing."

**👤 You:** "Can you compare this pattern to how normal stem cells are distributed?"

---

### 🎯 When to Use Ripley's K Function

**Used for:**

- 🧬 **Stem cell studies** - understanding niche organization
- 🩸 **Vascular analysis** - blood vessel distribution patterns
- 🦠 **Immune surveillance** - patrol pattern analysis
- 🔬 **Development studies** - how cells organize during morphogenesis

### 💡 Natural Ways to Ask

**Questions you can ask:**

- "Are my cells randomly distributed or clustered?"
- "Do these cell types form organized patterns?"
- "What's the typical spacing between these cells?"
- "Analyze the distribution pattern of my marked cells"
- "Compare clustering patterns between conditions"

## Method 5: Spatial Centrality Analysis

### Understanding Network-Based Spatial Metrics

**What it does:** Identifies the most "important" locations in your tissue based on spatial connectivity - which regions are hubs that connect different tissue areas?

**Biological applications:**

- Identifying tissue organizing centers
- Finding communication bottlenecks
- Discovering critical transition zones
- Understanding information flow patterns

### 💬 Real Conversation Example

**👤 You:** "In my brain tissue, I want to find the most important regions that connect different brain areas. Which locations act as hubs in the spatial network?"

**🤖 Claude:** "Good question about brain connectivity! Spatial centrality analysis will identify the hub regions that are most critical for connecting different parts of your brain tissue - like finding the major intersections in a city."

*[ChatSpatial analyzes spatial network centrality]*

**🤖 Claude:** "Analysis reveals spatial network architecture:

**🎯 Spatial Centrality Analysis Results:**

**Top Network Hubs (Highest Betweenness Centrality):**
- **Hub 1**: Corpus callosum region (centrality = 0.89)
  - *Role*: Bridge between left-right hemispheres
  - *Connection count*: Links 23 different brain regions
- **Hub 2**: Hippocampal-cortical boundary (centrality = 0.76)
  - *Role*: Memory-cortex integration center
  - *Connection count*: Links 18 regions
- **Hub 3**: White matter junction (centrality = 0.68)
  - *Role*: Multi-tract convergence zone
  - *Connection count*: Links 15 regions

**Network Properties:**
- **Most Connected Spots**: Average 12.4 connections per hub
- **Network Efficiency**: 0.73 (highly efficient connectivity)
- **Clustering Coefficient**: 0.58 (modular organization)
- **Small-World Index**: 2.1 (optimal information transfer)

**Critical Pathways:**
- **Information bottlenecks**: 3 spots that, if removed, would disconnect major brain regions
- **Redundant pathways**: Multiple routes between cortical areas (robust network)
- **Modular boundaries**: Clear boundaries between functional brain modules

**🔍 Key Finding:** Your brain tissue shows classic small-world network architecture with efficient hubs that minimize connection distances while maintaining modular organization."

**👤 You:** "Can you identify which genes are most highly expressed in these network hubs?"

---

### 🎯 When to Use Spatial Centrality Analysis

**Used for:**

- 🧠 **Brain connectivity studies** - understanding neural network architecture
- 🩸 **Vascular network analysis** - identifying critical vessel junctions
- 🔬 **Developmental studies** - finding organizing centers
- 🏥 **Pathology research** - understanding how disease disrupts connectivity

### 💡 Natural Ways to Ask

**Questions you can ask:**

- "Which regions are the most important hubs in my tissue?"
- "Find the network centers in my spatial data"
- "What are the critical connection points?"
- "Identify bottlenecks in my tissue architecture"
- "Show me the spatial network structure"

## Method 6: Join Count Analysis

### Understanding Categorical Spatial Autocorrelation

**What it does:** Measures whether similar categories (like cell types or tissue states) tend to cluster together spatially more than expected by chance.

**Biological applications:**

- Tissue zonation analysis
- Cell state transition boundaries
- Pathological region organization
- Developmental compartment formation

### 💬 Real Conversation Example

**👤 You:** "I have liver tissue with healthy, fibrotic, and necrotic regions marked. Do similar tissue states cluster together, or are they randomly mixed? Can you analyze the spatial organization of these tissue states?"

**🤖 Claude:** "Good question for understanding tissue pathology organization! Join count analysis will reveal whether your liver tissue states show spatial clustering - this tells us about disease progression patterns and tissue organization principles."

*[ChatSpatial performs join count analysis]*

**🤖 Claude:** "Analysis reveals tissue organization patterns:

**🔍 Join Count Analysis Results:**

**Healthy-Healthy Joins:**
- **Observed**: 287 adjacent healthy-healthy pairs
- **Expected**: 198 pairs (if randomly distributed)
- **Z-score**: +6.8 (p < 0.001)
- **Interpretation**: Healthy tissue strongly clusters together

**Fibrotic-Fibrotic Joins:**
- **Observed**: 156 adjacent fibrotic-fibrotic pairs
- **Expected**: 89 pairs (if random)
- **Z-score**: +5.2 (p < 0.001)
- **Interpretation**: Fibrotic regions form connected zones

**Necrotic-Necrotic Joins:**
- **Observed**: 78 adjacent necrotic-necrotic pairs
- **Expected**: 34 pairs (if random)
- **Z-score**: +7.1 (p < 0.001)
- **Interpretation**: Necrotic areas cluster in discrete patches

**Cross-State Boundaries:**
- **Healthy-Fibrotic**: Less than expected (Z = -3.4) - sharp boundaries
- **Fibrotic-Necrotic**: More than expected (Z = +2.8) - transition zones
- **Healthy-Necrotic**: Much less than expected (Z = -5.9) - rarely adjacent

**🎯 Key Finding:** Your liver shows organized pathological zonation - healthy regions remain protected, fibrotic areas form transition zones around necrotic cores, and there are sharp boundaries preventing healthy-necrotic contact."

**👤 You:** "Can you quantify how sharp these boundaries are between tissue states?"

---

### 🎯 When to Use Join Count Analysis

**Used for:**

- 🏥 **Pathology studies** - understanding disease organization
- 🧬 **Cell state analysis** - measuring state transition sharpness
- 🌱 **Development research** - compartment boundary formation
- 🔬 **Tissue zonation** - liver zones, kidney regions, brain layers

### 💡 Natural Ways to Ask

**Questions you can ask:**

- "Do similar cell types cluster together spatially?"
- "Are my tissue states randomly mixed or organized?"
- "How sharp are the boundaries between these regions?"
- "Analyze clustering of my categorical data"
- "Do disease states show spatial organization?"

## Specialized Integration: Combining Multiple Methods

### 💬 Multi-Method Analysis Example

**👤 You:** "I want to do a complete spatial analysis of my brain development data. Can you combine multiple spatial statistics methods to give me the full picture?"

**🤖 Claude:** "Good approach! A multi-method spatial analysis will reveal different aspects of your brain development patterns. Let me run a complete spatial statistics battery."

*[ChatSpatial performs integrated analysis]*

**🤖 Claude:** "Complete spatial architecture analysis finished:

**🎯 Integrated Spatial Statistics Results:**

**Method 1 - LISA Hotspots:**
- **NEUROG2**: 4 neurogenesis hotspots in ventricular zones
- **SHH**: Organizing center hotspot in floor plate region
- **PAX6**: Dorsal telencephalon development hotspot

**Method 2 - Bivariate Correlations:**
- **NEUROG2 ↔ ASCL1**: r = 0.73 (proneural gene coordination)
- **SHH ↔ GLI1**: r = 0.84 (signaling pathway integrity)
- **WNT1 ↔ EN1**: r = 0.69 (midbrain patterning coordination)

**Method 3 - Cell Type Neighborhoods:**
- **Neural Progenitors ↔ Radial Glia**: Strong association (Z = +5.2)
- **Differentiating Neurons ↔ Progenitors**: Moderate association (Z = +2.8)
- **Mature Neurons ↔ Progenitors**: Avoidance pattern (Z = -3.1)

**Method 4 - Network Centrality:**
- **Hub 1**: Ventricular zone organizing center (centrality = 0.92)
- **Hub 2**: Cortical plate formation zone (centrality = 0.76)
- **Hub 3**: Boundary region coordinators (centrality = 0.64)

**Method 5 - Ripley's K Patterns:**
- **Neural Progenitors**: Strong clustering (50-100 μm scale)
- **Differentiating Neurons**: Transition from clustered to regular spacing
- **Mature Neurons**: Regular spacing pattern (optimal connectivity)

**🧠 Integrated Biological Interpretation:** Your developing brain shows classic developmental organization - organizing centers create morphogen gradients (LISA), coordinate gene networks (bivariate), establish cell type territories (neighborhoods), form connectivity hubs (centrality), and regulate cell distribution patterns (Ripley's K) in a coordinated developmental program."

---

## Method Selection Guide

### 🎯 Which Method Should I Use?

| Research Question | Best Method | Why This Method? |
|------------------|-------------|------------------|
| "Do these two genes coordinate spatially?" | **Bivariate Moran's I** | Measures spatial correlation between gene pairs |
| "Where are the hotspots for this process?" | **LISA Analysis** | Identifies local clusters and spatial outliers |
| "How do cell types organize spatially?" | **Neighborhood Enrichment** | Quantifies cell-cell spatial associations |
| "Are cells randomly distributed?" | **Ripley's K Function** | Tests point pattern randomness vs clustering |
| "Which regions are network hubs?" | **Spatial Centrality** | Identifies important connectivity nodes |
| "Do similar states cluster together?" | **Join Count Analysis** | Measures categorical spatial autocorrelation |

### 💡 Combining Methods Strategy

**Start with LISA** → Find interesting hotspots
**Add Bivariate Moran's** → Understand gene coordination in hotspots
**Include Neighborhood Analysis** → See how cell types organize around hotspots
**Apply Ripley's K** → Understand cell distribution patterns
**Finish with Centrality** → Identify network importance of regions

## Best Practices and Tips

### ✅ Before Analysis

**Data Preparation:**

- Ensure spatial coordinates are accurate
- Check for batch effects or technical artifacts
- Verify cell type annotations are reliable
- Remove low-quality spots/cells

**Parameter Selection:**

- Start with standard neighbor numbers (6-30)
- Test different spatial scales for your tissue
- Use appropriate significance thresholds
- Consider multiple testing correction

### ✅ During Analysis

**Quality Control:**

- Check that results make biological sense
- Validate with known positive controls
- Compare across different parameter settings
- Look for technical artifacts (edge effects, etc.)

**Interpretation:**

- Consider both statistical and biological significance
- Account for tissue architecture and biology
- Compare with literature and atlases
- Think about multiple scales of organization

### ✅ After Analysis

**Validation:**

- Cross-validate with independent datasets
- Confirm key findings with targeted experiments
- Compare results across analysis methods
- Test robustness to parameter choices

**Integration:**

- Combine with other omics data types
- Link to functional studies
- Connect to biological pathways
- Consider clinical or therapeutic relevance

## Common Pitfalls to Avoid

### 🚨 Statistical Issues

**❌ Multiple Testing Problems:**

- Testing hundreds of genes without correction
- *Solution*: Use FDR correction, focus on strong effects

**❌ Ignoring Spatial Scale:**

- Using inappropriate neighbor distances
- *Solution*: Test multiple scales, consider tissue architecture

**❌ Technical Artifacts:**

- Edge effects, batch effects, imaging artifacts
- *Solution*: Quality control, technical replicates

### 🚨 Interpretation Issues

**❌ Over-interpretation:**

- Assuming correlation implies causation
- *Solution*: Validate with experiments, consider alternative explanations

**❌ Ignoring Biology:**

- Focusing only on statistics, not biological meaning
- *Solution*: Connect results to known biology, literature

**❌ Single-Method Bias:**

- Relying on one analysis approach
- *Solution*: Use multiple complementary methods

## Troubleshooting Guide

### "No significant results found"
**Possible causes:**

- Spatial scale mismatch with biological processes
- Low signal-to-noise ratio in data
- Inappropriate parameter settings

**Solutions to try:**

- Adjust neighbor parameters (try 6, 15, 30 neighbors)
- Focus on highly variable genes
- Check data quality and preprocessing
- Reduce statistical stringency initially

### "Results do not make biological sense"
**Possible causes:**

- Technical artifacts or batch effects
- Incorrect spatial coordinates
- Poor cell type annotations

**Solutions to try:**

- Validate spatial coordinates
- Check for batch effects
- Compare with positive controls
- Examine raw data quality

### "Analysis takes too long"
**Possible causes:**

- Large datasets with many features
- Computationally intensive methods
- Suboptimal parameters

**Solutions to try:**

- Focus on subset of highly variable genes
- Use spatial subsampling for exploration
- Optimize neighbor graph computation
- Run analysis in batches

## Summary: Your Spatial Statistics Toolkit

Congratulations! You now know how to:

🎯 **Perform six spatial statistics methods** through natural conversation
📊 **Choose the right method** for your biological question
🔬 **Interpret results biologically** with confidence
🧬 **Combine multiple approaches** for complete analysis
⭐ **Avoid common pitfalls** and troubleshoot issues

### 🚀 Next Steps

Ready to explore more advanced spatial analysis? Try these tutorials:

- **[Spatial Gene Set Enrichment](../analysis/spatial_enrichment.md)** - Map biological pathways spatially
- **[Cell Communication Analysis](../analysis/cell_communication_analysis.md)** - Discover cell-cell signaling patterns
- **[Trajectory Analysis](../advanced/trajectory_analysis.md)** - Track developmental and dynamic processes

### 📚 Additional Resources

- **[Spatial Statistics Theory](https://rspatial.org/raster/analysis/3-spauto.md)** - Mathematical background
- **[PySAL Documentation](https://pysal.org/index.md)** - Spatial analysis library details
- **[Spatial Transcriptomics Best Practices](../guides/index.md)** - Analysis guidelines

---

**Ready to identify spatial patterns in your tissue?** Start asking ChatSpatial about the spatial statistics methods that interest you most.