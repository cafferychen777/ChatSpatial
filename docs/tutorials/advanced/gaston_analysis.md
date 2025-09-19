---
layout: default
title: GASTON Analysis
parent: Advanced Tutorials
grand_parent: Tutorials
nav_order: 1
description: "Machine learning-powered spatial gene discovery with GASTON"
---

# GASTON: ML-Powered Spatial Gene Discovery

**Machine Learning Approach to Spatial Gene Analysis**

GASTON (Generative Adversarial Spatial Transcriptomics Optimization Network) is a machine learning method for identifying spatial gene patterns in transcriptomics data. This tutorial demonstrates how to use GASTON for spatial gene discovery through ChatSpatial's conversational interface.

## What Makes GASTON Different

### Traditional vs ML-Powered Approaches

**Traditional Methods** (SpatialDE, SPARK):
- Statistical correlation-based detection
- Limited to simple linear relationships
- Binary classifications (spatial vs non-spatial)
- Manual parameter tuning required

**GASTON Machine Learning Approach**:
- **Neural network pattern recognition** that learns complex spatial relationships
- **Isodepth learning** that maps gene expression using topographic methods
- **Gradient classification** that distinguishes continuous from discontinuous patterns
- **Adversarial training** that identifies spatial domains
- **Conversational interface** that provides access to ML methods

### When to Choose GASTON

Use GASTON when you need:

- **Complex pattern detection** beyond simple correlations
- **Continuous gradient analysis** to understand expression transitions
- **Spatial domain discovery** with adversarial neural networks
- **High-resolution pattern mapping** with isodepth visualization
- **Advanced feature engineering** through GLM-PCA preprocessing

Use Traditional Methods (SpatialDE/SPARK) when:

- Simple statistical correlation is sufficient
- Quick screening of large gene lists
- Computational resources are limited
- Established statistical significance testing needed

## Core GASTON Concepts

### Neural Network Architecture

GASTON uses a dual-encoder architecture:

- **Spatial Encoder**: Maps tissue coordinates to latent spatial representations
- **Expression Encoder**: Processes gene expression through advanced feature engineering
- **Adversarial Training**: Competing networks refine spatial domain boundaries

### Isodepth Learning and Visualization

Think of GASTON's isodepth maps as **topographic elevation maps** for gene expression:

- **Isodepth lines**: Connect points of equal expression "elevation"
- **Gradients**: Show direction and rate of expression change
- **Peaks and valleys**: Identify expression hotspots and cold regions
- **Terrain mapping**: Visualize complex spatial expression landscapes

### Gradient Classification System

GASTON classifies genes into distinct spatial patterns:

**Continuous Gradients**:

- Smooth expression transitions across tissue
- No sharp boundaries
- Gradual spatial changes
- Examples: Morphogen gradients, metabolic zones

**Discontinuous Patterns**:

- Sharp expression boundaries
- Distinct spatial domains
- Binary on/off patterns  
- Examples: Cell type markers, tissue boundaries

### GLM-PCA Feature Engineering

GASTON's preprocessing pipeline includes advanced feature engineering:

**Standard Approach**: Basic log-normalization and PCA

**GASTON's GLM-PCA**: Generalized Linear Model PCA that:

- Accounts for count data distribution properties
- Reduces technical noise while preserving biological signal
- Improves features for spatial pattern detection
- Enables more accurate neural network training

## Real Conversation Examples

### Starting a GASTON Analysis

Here's how typical GASTON conversations unfold:

---

**User:** "I want to analyze spatial gene patterns in my mouse brain data using ML methods. Can you use GASTON to find genes with complex spatial patterns?"

**Claude:** "I will use GASTON to analyze complex spatial gene patterns using neural networks. Let me set up the ML analysis for your mouse brain data.

I'll use GASTON's neural network approach with:

- **Isodepth learning** to map expression topology
- **Gradient classification** to distinguish pattern types
- **GLM-PCA feature engineering** for optimal preprocessing
- **Adversarial domain discovery** to find spatial territories

Let me start the GASTON analysis now..."

*[GASTON neural network training begins]*

**Claude:** "**GASTON Analysis Complete**

GASTON has identified spatial patterns in your data:

**📈 Continuous Gradient Genes (78 genes)**:

- *Hpca*: Hippocampal gradient marker (isodepth correlation: 0.92)
- *Camk2a*: Learning/memory pathway gradient  
- *Grin1*: NMDA receptor spatial organization

**🔲 Discontinuous Pattern Genes (134 genes)**:

- *Slc17a7*: Sharp cortical layer boundaries
- *Gad2*: Distinct inhibitory territories
- *Mbp*: White matter domain marker

**Spatial Domains Discovered**: 8 distinct ML-identified territories
**Isodepth Range**: [0.124, 0.887] expression elevation units

The neural network achieved 87% spatial pattern reconstruction accuracy."

---

### Diving Deeper into GASTON Results

**User:** "Can you show me the isodepth visualization for the hippocampal genes? I want to understand the topographic patterns."

**Claude:** "I will create the isodepth visualization for your hippocampal genes. Let me generate the topographic expression maps."

*[Creates GASTON isodepth plots]*

**Claude:** "**Hippocampal Expression Topography**

Here's what GASTON's isodepth analysis reveals:

**Hpca Expression Landscape**:

- **Peak elevation**: Hippocampal CA fields (isodepth: 0.85)
- **Gradient slope**: Smooth transition to cortex (slope: -0.23/μm)
- **Valley regions**: Corpus callosum (isodepth: 0.12)
- **Contour pattern**: Concentric expression rings around CA1/CA3

**Camk2a Terrain**:

- **Ridge formation**: Pyramidal layer enrichment
- **Plateau regions**: Stable cortical expression
- **Drainage patterns**: Expression flows toward subcortical regions

This topographic view shows how gene expression varies across your tissue. GASTON's neural networks have learned to identify these molecular patterns."

---

### Comparing GASTON vs Traditional Methods

**User:** "How does GASTON's analysis compare to traditional methods like SpatialDE? Are we finding different genes?"

**Claude:** "I will run a comparative analysis between GASTON and SpatialDE to show the differences between ML and statistical approaches."

*[Runs both analyses]*

**Claude:** "**GASTON vs SpatialDE Comparison**

**Unique GASTON Discoveries (42 genes)**:

- Complex gradient patterns missed by correlation-based methods
- Multi-modal spatial distributions
- Subtle expression terrains requiring ML pattern recognition

**SpatialDE Exclusive Genes (28 genes)**:

- Strong linear correlations with distance
- Simple radial patterns
- Statistically significant but spatially basic

**Overlapping Findings (156 genes)**:

- Clear spatial markers validated by both approaches
- Strong confidence in these spatial genes

**Key Finding**: GASTON's neural network discovers **64% more complex patterns** that traditional statistics miss, while confirming all major spatial markers. The ML approach identifies additional topographic complexity in the tissue."

---

### Advanced GASTON Parameters

**User:** "Can you tune GASTON for more sensitive pattern detection? I suspect there are subtle gradients I am missing."

**Claude:** "I will tune GASTON's neural architecture for more sensitive pattern detection. I will use this configuration:

- **Deeper spatial encoder**: 5 hidden layers for complex pattern learning
- **Enhanced adversarial training**: More discriminator updates
- **Fine-grained isodepth bins**: 50 elevation levels
- **GLM-PCA preprocessing**: 100 components for maximum signal capture
- **Extended training**: 2000 epochs for pattern convergence

This high-sensitivity configuration will detect more subtle spatial gradients."

*[Advanced GASTON analysis]*

**Claude:** "**High-Sensitivity GASTON Results**

The enhanced ML analysis identified **23 additional subtle gradient genes**:

**Newly Detected Patterns**:

- *Atp1a3*: Ultra-fine cortical lamination (gradient strength: 0.34)
- *Snap25*: Micro-domain organization in CA1
- *Grin2b*: Developmental gradient remnants

**Pattern Details**:

- Average gradient strength improved from 0.45 to 0.62
- Spatial resolution enhanced by 34%
- Domain boundary precision increased to 94%

The neural networks have learned to identify molecular patterns not detected by traditional analysis."

---

## Understanding GASTON Visualizations

### Isodepth Expression Maps

GASTON's signature visualization shows gene expression as topographic elevation:

**Reading the Map**:

- **Warm colors (red/yellow)**: Expression peaks and ridges  
- **Cool colors (blue/green)**: Expression valleys and basins
- **Contour lines**: Connect points of equal expression
- **Arrows**: Show gradient direction and strength
- **Elevation legend**: Maps color to expression level

**Biological Interpretation**:

- **Mountain ranges**: Core expression domains
- **River valleys**: Expression boundaries between regions
- **Plateaus**: Stable, uniform expression zones
- **Cliff faces**: Sharp discontinuous transitions

### Spatial Domain Visualization

GASTON's adversarial training discovers distinct spatial territories:

**Domain Features**:

- **Color-coded regions**: Each ML-discovered domain
- **Boundary precision**: Sharp vs gradual transitions
- **Domain statistics**: Size, gene content, connectivity
- **Hierarchical organization**: Sub-domains within domains

### Gradient Classification Plots

Distinguishes continuous from discontinuous spatial patterns:

**Continuous Genes**:

- Smooth gradient arrows across tissue
- Progressive color transitions
- Flow-like pattern visualization

**Discontinuous Genes**:

- Sharp boundary demarcations  
- Binary on/off regions
- Territory-based coloring

## Integration with Other ChatSpatial Analyses

### Cell Communication + GASTON

**User:** "Can GASTON help me understand how spatial gene gradients relate to cell communication patterns?"

**Claude:** "I will integrate GASTON's spatial gradient analysis with cell communication mapping to show how expression topography relates to cellular crosstalk."

*[Integrates GASTON gradients with ligand-receptor analysis]*

**Claude:** "**Spatial Gradients and Cell Communication**

**Gradient-Driven Communication**:

- **Wnt3a gradient** creates **Fzd signaling streams** across cortical layers
- **Bdnf expression peaks** drive **TrkB receptor hotspots** in hippocampus  
- **Vegfa valleys** correspond to **low angiogenic signaling** zones

**Communication Topology**:

- Ligand gradients create **directional signaling flows**
- Receptor distributions follow **isodepth contours**
- Communication strength correlates with **gradient steepness**

GASTON shows that the tissue's communication network follows the molecular topography."

### Spatial Domains + GASTON

**User:** "How do GASTON's ML-discovered domains compare to traditional clustering methods?"

**Claude:** "I will analyze how GASTON's adversarial domain discovery differs from standard spatial clustering approaches."

*[Compares GASTON domains with Leiden clustering]*

**Claude:** "**ML Domains vs Traditional Clustering**

**GASTON Advantages**:

- **Gradient-aware boundaries**: Domains respect expression flows
- **Multi-scale detection**: Discovers both large territories and micro-domains
- **Pattern-based grouping**: Groups by spatial behavior, not just similarity

**Traditional Clustering Strengths**:

- **Cell type accuracy**: Better for discrete cell populations
- **Computational efficiency**: Faster for large datasets
- **Statistical validation**: Established significance testing

**Best Practice**: Use GASTON for **spatial pattern discovery**, traditional clustering for **cell type identification**, then integrate results for comprehensive tissue understanding."

## Advanced GASTON Workflows

### Multi-Sample GASTON Analysis

**User:** "I have multiple brain sections. Can GASTON find conserved spatial patterns across samples?"

**Claude:** "GASTON can learn conserved spatial motifs across multiple samples using transfer learning approaches."

### Temporal GASTON Analysis  

**User:** "Can GASTON track how spatial patterns change during development?"

**Claude:** "GASTON's temporal extension can map how expression topography evolves across developmental stages."

### GASTON Parameter Optimization

| Parameter | Default | High Sensitivity | Fast Screening |
|-----------|---------|------------------|----------------|
| Spatial Hidden Layers | [128, 64, 32] | [256, 128, 64, 32, 16] | [64, 32] |
| Expression Hidden Layers | [64, 32] | [128, 64, 32, 16] | [32, 16] |  
| Training Epochs | 500 | 2000 | 200 |
| Isodepth Bins | 20 | 50 | 10 |
| Learning Rate | 0.001 | 0.0005 | 0.002 |
| Batch Size | 64 | 32 | 128 |

## Troubleshooting GASTON

### Common Issues and Solutions

#### "Neural network training failed"

```bash
Solution: Reduce batch size or learning rate
Alternative: Use Pearson residuals instead of GLM-PCA
```

#### "No spatial patterns detected"

```bash
Check: Data preprocessing quality
Try: Increase training epochs or reduce regularization
```

#### "Isodepth visualization looks noisy"

```bash
Adjust: Increase number of isodepth bins
Filter: Use only high-confidence gradient genes
```

#### "GASTON finds too many genes"

```bash
Tune: Increase gradient thresholds
Focus: Use discontinuous genes only for cleaner patterns
```

## Machine Learning in Spatial Biology

GASTON represents an example of ML-driven spatial transcriptomics. Future developments may include:

- **Foundation models** pre-trained on large spatial datasets
- **Multi-modal integration** combining imaging, spatial, and temporal data
- **Predictive spatial biology** for expression pattern forecasting
- **Real-time analysis** during tissue collection
- **Personalized spatial medicine** based on individual tissue patterns

ChatSpatial's conversational interface provides access to these ML capabilities for spatial analysis.

## Key Takeaways

1. **GASTON uses ML** to find complex patterns missed by statistics
2. **Isodepth visualization** reveals expression topography like terrain maps
3. **Gradient classification** distinguishes smooth from sharp spatial transitions
4. **Neural networks learn** spatial relationships beyond human intuition
5. **Conversational interface** makes advanced ML accessible to biologists
6. **Integration capabilities** combine GASTON with other spatial analyses

To begin spatial gene discovery with GASTON, start a conversation with ChatSpatial and ask: *"Can you use GASTON to find ML-powered spatial patterns in my data?"*

---

**Next Steps**:

- Try [Cell Communication Analysis](../analysis/cell_communication_analysis.md) with GASTON gradients
- Explore [Spatial Enrichment](../analysis/spatial_enrichment.md) using GASTON domains  
- Create [Advanced Visualizations](../core/visualization_tutorial.md) of isodepth patterns

*Explore the molecular landscape of your tissue through ML-powered conversation.*
