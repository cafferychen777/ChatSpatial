
title: Cell Type Annotation
description: Identifying cell types in spatial transcriptomics data
---

# Identifying Cell Types with ChatSpatial

Discover what types of cells are in your spatial transcriptomics data through natural conversation with ChatSpatial.

## What You'll Learn

By the end of this tutorial, you'll know how to:

1. **Ask ChatSpatial to identify cell types** in your tissue
2. **Understand different annotation methods** and when to use them
3. **Explore cell type distributions** across your tissue
4. **Validate and refine** cell type assignments
5. **Visualize cell types** in spatial context

## Before You Start

### What You Need
- ✅ **Spatial data loaded** in ChatSpatial
- ✅ **Basic preprocessing completed** (see [Getting Started](../../getting-started/index.md))
- ✅ **Curiosity about** what cell types are in your tissue!

### No Cell Biology PhD Required!
This tutorial is designed for anyone curious about their tissue composition. ChatSpatial will guide you through the biology and help interpret the results.

## Step 1: Starting Cell Type Identification

### 💬 The Discovery Begins

**👤 You:** "I have my spatial data loaded and I can see different regions. Now I'm curious - what types of cells are actually in my tissue?"

**🤖 ChatSpatial:** "That's one of the most exciting questions in spatial biology! I can help you identify the cell types in your tissue using several approaches. Let me start with marker gene analysis to give you a quick overview."

*[ChatSpatial analyzes marker gene expression]*

**🤖 ChatSpatial:** "Fascinating! Based on marker gene expression, I can see several distinct cell types in your brain tissue:

**Identified Cell Types:**
- **Neurons** (35% of spots) - High NeuN, MAP2, RBFOX3
- **Astrocytes** (28% of spots) - High GFAP, AQP4, S100B  
- **Oligodendrocytes** (18% of spots) - High MBP, PLP1, MOG
- **Microglia** (12% of spots) - High IBA1, CX3CR1, P2RY12
- **Endothelial cells** (5% of spots) - High PECAM1, VWF, CDH5
- **Pericytes** (2% of spots) - High PDGFRB, ACTA2, RGS5

This looks like a healthy brain tissue with all major cell types represented! Would you like me to show you where these cell types are located spatially, or explore any specific cell type in more detail?"

---

### 🎯 What Just Happened?

ChatSpatial automatically:
- ✅ Analyzed expression of known marker genes
- ✅ Assigned cell type labels based on marker expression
- ✅ Calculated the proportion of each cell type
- ✅ Provided the key genes that define each cell type

## Step 2: Exploring Spatial Distribution

### 💬 Seeing Where Cell Types Live

**👤 You:** "This is amazing! I'd love to see where these different cell types are located in my tissue. Can you show me a spatial map?"

**🤖 ChatSpatial:** "Absolutely! Let me create a spatial cell type map for you."

*[Beautiful spatial visualization appears showing different colored regions]*

**🤖 ChatSpatial:** "Here's your spatial cell type map! 🧠✨

**Spatial Organization:**
- **Neurons** (blue): Concentrated in cortical layers and hippocampus
- **Astrocytes** (green): Distributed throughout, higher in gray matter
- **Oligodendrocytes** (yellow): Enriched in white matter regions
- **Microglia** (red): Scattered throughout with some clustering
- **Endothelial cells** (purple): Following vascular patterns
- **Pericytes** (orange): Co-localized with endothelial cells

**Biological Insights:**
Your tissue shows beautiful anatomical organization! The cell types are distributed exactly as we'd expect in healthy brain tissue:
- Neurons form the functional layers
- Oligodendrocytes cluster in myelinated areas
- Vascular cells outline blood vessels
- Microglia patrol throughout the tissue

This suggests your tissue has preserved its natural architecture!"

**👤 You:** "That's so cool! I notice the neurons seem to be in distinct layers. Can you tell me more about the neuronal populations?"

**🤖 ChatSpatial:** "Excellent observation! Let me dive deeper into the neuronal subtypes."

*[ChatSpatial performs detailed neuronal analysis]*

**🤖 ChatSpatial:** "Great eye for detail! I can identify several neuronal subtypes based on their marker expression:

**Neuronal Subtypes:**
- **Excitatory neurons** (Layer 2/3): High SLC17A7, CAMK2A
- **Excitatory neurons** (Layer 5): High BCL11B, FEZF2
- **Inhibitory neurons** (GABAergic): High GAD1, GAD2, PVALB
- **Pyramidal neurons**: High TBR1, SATB2
- **Interneurons**: High SST, VIP, CCK

**Spatial Pattern:**
- **Layer 2/3**: Superficial cortical neurons
- **Layer 5**: Deep cortical projection neurons  
- **Interneurons**: Scattered throughout cortical layers
- **Hippocampal neurons**: Distinct CA1/CA3 populations

This layered organization is a hallmark of healthy cortical development!"

---

## Step 3: Refining Cell Type Assignments

### 💬 Getting More Precise

**👤 You:** "This is getting really detailed! But I want to make sure these cell type assignments are accurate. Can you use a more sophisticated method to validate them?"

**🤖 ChatSpatial:** "Absolutely! Let me use a more advanced method called scType that uses comprehensive gene signatures for more precise cell type identification."

*[ChatSpatial runs advanced cell type annotation]*

**🤖 ChatSpatial:** "Excellent idea! The scType analysis confirms and refines our initial findings:

**Refined Cell Type Annotations:**
- **Excitatory neurons**: 
  - Layer 2/3 pyramidal (15%)
  - Layer 5 pyramidal (12%)
  - Layer 6 pyramidal (8%)
- **Inhibitory neurons**:
  - Parvalbumin+ interneurons (6%)
  - Somatostatin+ interneurons (4%)
  - VIP+ interneurons (2%)
- **Glial cells**:
  - Protoplasmic astrocytes (20%)
  - Fibrous astrocytes (8%)
  - Mature oligodendrocytes (15%)
  - Oligodendrocyte precursors (3%)
- **Other cells**:
  - Microglia (7%)
  - Endothelial cells (5%)
  - Pericytes (2%)

**Confidence Scores:**
- High confidence (>0.8): 78% of spots
- Medium confidence (0.6-0.8): 18% of spots  
- Low confidence (<0.6): 4% of spots

The high confidence scores suggest these annotations are very reliable!"

**👤 You:** "Wow, that's much more detailed! What about those spots with low confidence? What might they be?"

**🤖 ChatSpatial:** "Great question! Let me investigate those ambiguous spots."

*[ChatSpatial analyzes uncertain regions]*

**🤖 ChatSpatial:** "The low-confidence spots are interesting! They fall into a few categories:

**Ambiguous Regions:**
- **Transition zones** (40%): Between tissue layers with mixed signals
- **Damaged areas** (25%): Lower gene expression, possible artifacts
- **Rare cell types** (20%): Possibly ependymal cells or choroid plexus
- **Mixed spots** (15%): Multiple cell types in single measurement

**Recommendations:**
- **Transition zones**: Normal - these are biologically meaningful boundaries
- **Damaged areas**: Consider filtering if they affect analysis
- **Rare cell types**: Interesting! Might warrant further investigation
- **Mixed spots**: Expected in spatial data due to resolution limits

Overall, 96% confident annotation is excellent for spatial data!"

---

## Step 4: Biological Validation

### 💬 Making Sure It Makes Sense

**👤 You:** "This all looks great, but how do I know these cell type assignments actually make biological sense?"

**🤖 ChatSpatial:** "Excellent scientific thinking! Let me validate these assignments using multiple approaches."

*[ChatSpatial performs biological validation]*

**🤖 ChatSpatial:** "Let me validate these cell types using several biological criteria:

**✅ Validation Results:**

**1. Marker Gene Consistency:**
- Neurons: NeuN+ (98% of assigned spots)
- Astrocytes: GFAP+ (95% of assigned spots)
- Oligodendrocytes: MBP+ (92% of assigned spots)
- Microglia: IBA1+ (94% of assigned spots)

**2. Spatial Organization:**
- Cortical layering preserved ✅
- White matter enriched in oligodendrocytes ✅
- Vascular cells follow blood vessel patterns ✅
- Microglia evenly distributed ✅

**3. Gene Expression Profiles:**
- Each cell type shows expected transcriptional signatures ✅
- No unexpected marker combinations ✅
- Expression levels match literature values ✅

**4. Proportions Check:**
- Cell type ratios match published brain composition ✅
- No impossible cell type combinations ✅

**Biological Validation Score: 96/100** - Excellent!"

**👤 You:** "That's reassuring! Are there any cell types that might be missing or any surprises in the data?"

**🤖 ChatSpatial:** "Good question! Let me check for potentially missing or unusual cell types."

*[ChatSpatial searches for additional cell types]*

**🤖 ChatSpatial:** "Interesting findings from the comprehensive search:

**Potentially Missing Cell Types:**
- **Ependymal cells**: Low abundance, might be present in ventricular areas
- **Choroid plexus**: Possible in a few spots near ventricles
- **Neural stem cells**: Very rare, might be in neurogenic niches

**Unexpected Findings:**
- **Activated microglia**: Small population showing inflammatory markers
- **Reactive astrocytes**: Subset with stress response signatures
- **Immature oligodendrocytes**: More than expected, suggests active myelination

**Biological Interpretation:**
- The **activated immune cells** might indicate normal tissue maintenance
- **Reactive astrocytes** could be response to tissue processing
- **Active myelination** suggests healthy, dynamic tissue

These findings actually add biological richness to your dataset!"

---

## Step 5: Creating Publication-Ready Visualizations

### 💬 Showing Your Results

**👤 You:** "This analysis is fantastic! Can you help me create some nice visualizations for my presentation or paper?"

**🤖 ChatSpatial:** "Absolutely! Let me create several publication-quality visualizations for you."

*[ChatSpatial generates multiple beautiful plots]*

**🤖 ChatSpatial:** "Here are your publication-ready cell type visualizations! 📊✨

**Generated Plots:**
1. **Spatial cell type map** - Overview of tissue organization
2. **Cell type proportions** - Bar chart with percentages
3. **Marker gene heatmap** - Expression patterns by cell type
4. **UMAP embedding** - Dimensionality reduction colored by cell type
5. **Confidence score map** - Reliability of annotations
6. **Layer-specific analysis** - Cortical organization detail

**Figure Suggestions:**
- **Main figure**: Spatial map + proportions chart
- **Supplementary**: Marker heatmap + UMAP + confidence scores
- **Methods validation**: Marker gene expression by cell type

All plots are high-resolution and ready for publication! Would you like me to adjust colors, labels, or create additional specialized plots?"

**👤 You:** "Perfect! One last question - how should I describe this analysis in my methods section?"

**🤖 ChatSpatial:** "Great question! Here's a methods description you can adapt:

**Suggested Methods Text:**
'Cell type annotation was performed using ChatSpatial's integrated pipeline. Initial cell type identification utilized marker gene expression analysis with established brain cell type markers. Annotations were refined using scType with comprehensive gene signature databases. Spatial organization was validated by examining anatomical distribution patterns and marker gene consistency. Cell type assignments achieved >95% confidence for the majority of spots, with validation scores confirming biological plausibility of the annotations.'

**Key Points to Include:**
- Multiple validation approaches used
- High confidence scores achieved
- Spatial organization preserved
- Biologically consistent results

Feel free to modify this based on your specific journal requirements!"

---

## Common Questions and Troubleshooting

### 💬 "What if ChatSpatial can't identify some cell types?"

**🤖 ChatSpatial:** "This can happen for several reasons:
- **Novel cell states**: Your tissue might have unique populations
- **Technical limitations**: Some cell types need higher resolution
- **Rare populations**: Very low-abundance cells are challenging
- **Mixed signals**: Spatial resolution might capture multiple cells

I can help you investigate these ambiguous regions and suggest targeted approaches!"

### 💬 "Can I add my own custom cell type markers?"

**🤖 ChatSpatial:** "Yes! If you have specific markers for your tissue or research focus, just tell me:
- 'Use these genes to identify cell type X'
- 'I think this region contains cell type Y based on gene Z'
- 'Can you check for expression of my custom marker list?'

I can incorporate your expertise into the analysis!"

### 💬 "How do I compare cell types between different samples?"

**🤖 ChatSpatial:** "Great question! I can help you compare:
- Cell type proportions between conditions
- Spatial distribution changes
- Marker gene expression differences
- Novel cell states in disease vs healthy

Just load multiple datasets and ask me to compare them!"

## Next Steps

### 🚀 Continue Your Analysis

**Now that you know your cell types, try:**
- **"How do these cell types communicate with each other?"** → [Cell Communication Tutorial](cell_communication_analysis.md)
- **"What genes are specifically expressed in each cell type?"**
- **"Can you find spatial domains within cell types?"**
- **"Show me how cell types change across tissue regions"**

### 📚 Learn More

- [Basic Spatial Analysis](../core/basic_spatial_analysis.md) - Start here if you're new
- [Cell Communication Analysis](cell_communication_analysis.md) - Next logical step
- [Visualization Tutorial](../core/visualization_tutorial.md) - Make beautiful plots

**Remember**: Cell type identification is the foundation for understanding tissue organization and function. Every cell type has a story to tell about your tissue's biology!