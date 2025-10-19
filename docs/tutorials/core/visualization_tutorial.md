---
layout: default
title: Visualization Tutorial
parent: Core Tutorials
grand_parent: Tutorials
nav_order: 3
description: "Creating beautiful visualizations with ChatSpatial"
---

# Creating Beautiful Visualizations with ChatSpatial

Learn how to create stunning, publication-ready visualizations of your spatial transcriptomics data through simple conversation with ChatSpatial.

## What You'll Create

By the end of this tutorial, you'll know how to:

1. **Ask for any type of plot** you can imagine
2. **Customize colors, styles, and layouts** through conversation
3. **Create publication-ready figures** automatically
4. **Combine multiple visualizations** into comprehensive panels
5. **Export high-quality images** for presentations and papers

## Before You Start

### What You Need
- ✅ **Spatial data loaded** in ChatSpatial
- ✅ **Some analysis completed** (spatial domains, cell types, etc.)
- ✅ **Ideas about** what you want to visualize!

### No Design Skills Required!
ChatSpatial handles all the technical details of creating beautiful plots. You just need to describe what you want to see!

## Step 1: Basic Spatial Visualizations

### 💬 Starting with Simple Requests

---

**👤 You:** "I'd like to see my spatial data. Can you create a nice visualization showing the tissue structure?"

**🤖 ChatSpatial:** "Absolutely! Let me create a beautiful spatial overview of your tissue."

*[Stunning spatial plot appears]*

**🤖 ChatSpatial:** "Here's your tissue overview! 🎨

**What you're seeing:**
- Each dot represents a measurement spot
- Colors show total gene expression (UMI counts)
- Spatial coordinates preserve tissue architecture
- Color scale from blue (low) to red (high) expression

**Visual insights:**
- Your tissue has clear regional differences
- Some areas show much higher gene activity
- The spatial organization looks well-preserved
- No obvious technical artifacts

Would you like me to color this by something specific, like cell types or a particular gene?"

**👤 You:** "Yes! Can you show me the same tissue but colored by the spatial domains we identified earlier?"

**🤖 ChatSpatial:** "Perfect! Let me create a spatial domain visualization."

*[Beautiful domain-colored spatial plot appears]*

**🤖 ChatSpatial:** "Here's your tissue colored by spatial domains! 🌈

**Domain Organization:**
- **Domain 1** (blue): Cortical layers - organized neural tissue
- **Domain 2** (green): White matter - myelinated regions  
- **Domain 3** (red): Hippocampus - memory formation area
- **Domain 4** (yellow): Ventricular zone - stem cell regions
- **Domain 5** (purple): Boundary areas - transition zones

The colors clearly show how your tissue is organized into distinct functional regions. Each domain has its own biological identity!"

---

### 🎯 What Just Happened?

ChatSpatial automatically:
- ✅ Created publication-quality spatial plots
- ✅ Chose appropriate color schemes
- ✅ Added informative legends and labels
- ✅ Provided biological interpretation of visual patterns

## Step 2: Gene Expression Visualizations

### 💬 Exploring Specific Genes

---

**👤 You:** "This looks great! Now I'm curious about specific genes. Can you show me where GFAP is expressed in my tissue?"

**🤖 ChatSpatial:** "Excellent choice! GFAP is a key astrocyte marker. Let me visualize its expression pattern."

*[Gene expression spatial plot appears]*

**🤖 ChatSpatial:** "Here's GFAP expression across your tissue! 🧬

**Expression Pattern:**
- **High expression** (red/yellow): Scattered throughout gray matter
- **Moderate expression** (orange): Transition zones
- **Low expression** (blue): White matter and neuronal layers

**Biological Interpretation:**
This pattern is exactly what we'd expect for GFAP! It's marking astrocytes, which are:
- Abundant in gray matter (supporting neurons)
- Less common in white matter (mainly oligodendrocytes there)
- Present at tissue boundaries (maintaining structure)

The spatial pattern confirms healthy astrocyte distribution!"

**👤 You:** "That's fascinating! Can you show me multiple genes at once? I'd like to compare GFAP with some neuronal markers."

**🤖 ChatSpatial:** "Great idea! Let me create a multi-gene comparison panel."

*[Multi-panel gene expression plot appears]*

**🤖 ChatSpatial:** "Here's your gene expression comparison panel! 🔬

**Four-Gene Panel:**
- **GFAP** (top-left): Astrocyte marker - scattered distribution
- **MAP2** (top-right): Neuronal marker - layer-specific pattern
- **MBP** (bottom-left): Oligodendrocyte marker - white matter enriched
- **IBA1** (bottom-right): Microglia marker - evenly distributed

**Comparative Insights:**
- **Complementary patterns**: Each gene marks different cell types
- **Spatial segregation**: Clear tissue organization visible
- **Expected biology**: All patterns match known cell distributions
- **Tissue health**: Normal cellular architecture preserved

This multi-gene view beautifully shows how different cell types organize your tissue!"

---

## Step 3: Advanced Visualization Techniques

### 💬 Getting Creative with Plots

---

**👤 You:** "These spatial plots are amazing! But can you create other types of visualizations? Maybe something showing the relationships between different regions?"

**🤖 ChatSpatial:** "Absolutely! Let me create some advanced visualizations to show different perspectives on your data."

*[Multiple advanced plots appear]*

**🤖 ChatSpatial:** "Here are several advanced visualizations of your data! 📊✨

**1. UMAP Embedding:**
- Shows relationships between all spots
- Colored by spatial domains
- Reveals transcriptional similarities
- Clusters match spatial organization

**2. Violin Plots:**
- Gene expression distributions by domain
- Shows expression variability
- Compares domains statistically
- Highlights domain-specific genes

**3. Heatmap:**
- Top marker genes for each domain
- Hierarchical clustering of genes and domains
- Color intensity shows expression level
- Reveals gene expression signatures

**4. Network Plot:**
- Connections between similar regions
- Node size = region importance
- Edge thickness = similarity strength
- Shows tissue connectivity patterns

Each plot reveals different aspects of your tissue organization!"

**👤 You:** "Wow, these are all so informative! Can you create a summary figure that combines the most important visualizations for a presentation?"

**🤖 ChatSpatial:** "Perfect! Let me create a comprehensive summary figure."

*[Beautiful multi-panel summary figure appears]*

**🤖 ChatSpatial:** "Here's your presentation-ready summary figure! 🎯

**Panel Layout:**
- **A**: Spatial domain map (main result)
- **B**: Cell type distribution (biological context)
- **C**: Key marker genes (validation)
- **D**: Domain statistics (quantification)

**Figure Features:**
- **Consistent color scheme** across all panels
- **Clear labels and legends** for easy interpretation
- **High resolution** ready for projection or printing
- **Logical flow** telling the story of your tissue

This figure tells the complete story: spatial organization → biological identity → molecular evidence → quantitative summary!"

---

## Step 4: Customizing Your Visualizations

### 💬 Making It Perfect

---

**👤 You:** "This summary figure is great, but I'd like to customize it a bit. Can you change the colors and maybe adjust the layout?"

**🤖 ChatSpatial:** "Of course! I love helping you get the perfect visualization. What specific changes would you like?"

**👤 You:** "Can you use a colorblind-friendly palette? And maybe make the spatial plot larger since that's the main result?"

**🤖 ChatSpatial:** "Excellent suggestions! Let me create a colorblind-friendly version with an improved layout."

*[Redesigned figure with new colors and layout appears]*

**🤖 ChatSpatial:** "Here's your customized figure! ♿🎨

**Improvements Made:**
- **Colorblind-friendly palette**: Using viridis colors that work for everyone
- **Larger spatial plot**: Now takes up 60% of the figure space
- **Improved contrast**: Better visibility for all viewers
- **Cleaner layout**: More space between panels
- **Enhanced labels**: Larger, clearer text

**Accessibility Features:**
- Colors distinguishable for all types of color vision
- High contrast ratios for readability
- Clear symbols and patterns as backup to color
- Large enough text for presentations

This version will work perfectly for diverse audiences!"

**👤 You:** "Perfect! Can you also create a version optimized for publication in a journal?"

**🤖 ChatSpatial:** "Absolutely! Let me create a publication-optimized version."

*[Journal-ready figure appears]*

**🤖 ChatSpatial:** "Here's your journal-ready figure! 📰

**Publication Optimizations:**
- **High DPI** (300+): Crisp printing quality
- **Vector graphics**: Scalable without quality loss
- **Journal color scheme**: Professional, conservative colors
- **Compact layout**: Efficient use of space
- **Standard fonts**: Times/Arial for compatibility
- **Clear scale bars**: Proper size references

**Format Options Available:**
- **PDF**: Vector format for journals
- **PNG**: High-res raster for online
- **SVG**: Editable vector format
- **EPS**: Legacy journal format

**Figure Caption Suggestion:**
'Spatial transcriptomics analysis reveals distinct tissue domains. (A) Spatial organization showing five distinct domains. (B) Cell type distribution across domains. (C) Marker gene expression validation. (D) Quantitative domain statistics.'

Ready for submission!"

---

## Step 5: Specialized Analysis Visualizations

### 💬 Exploring Advanced Plot Types

---

**👤 You:** "I've heard about spatial statistics and spatially variable genes. Can you show me what these advanced analyses look like?"

**🤖 ChatSpatial:** "Absolutely! Spatial statistics provide unique insights into spatial gene patterns and cell type organization. Let me create several specialized visualizations."

*[Spatial statistics plots appear]*

**🤖 ChatSpatial:** "Here are your spatial analysis results! 🧬📊

**Spatial Statistics Suite:**

**1. Moran's I Autocorrelation:**
- Measures spatial clustering of gene expression
- High Moran's I indicates strong spatial patterns
- Identifies genes with coordinated spatial expression
- Statistical significance testing included

**2. Spatially Variable Genes (SPARK-X/SpatialDE):**
- Highlights genes with significant spatial patterns
- Statistical significance color-coded
- Shows spatial correlation coefficients
- Identifies genes driving spatial organization

**3. Neighborhood Enrichment:**
- Analyzes cell type spatial relationships
- Shows which cell types prefer to co-localize
- Identifies spatial avoidance patterns
- Quantifies tissue architecture organization

**Biological Insights:**
- 234 genes show significant spatial variability
- Strong spatial autocorrelation indicates organized tissue architecture
- Cell type distributions show non-random spatial patterns
- Expression gradients follow anatomical boundaries

These patterns suggest well-preserved tissue architecture!"

**👤 You:** "That's fascinating! Can you show me more spatial statistics? I want to understand if my cell types are randomly distributed or if there are specific patterns."

**🤖 ChatSpatial:** "Great question! Spatial statistics can reveal whether cell distributions are random, clustered, or dispersed. Let me run several spatial analysis types."

*[Spatial statistics plots appear]*

**🤖 ChatSpatial:** "Here are your spatial statistics results! 📐🔍

**Spatial Analysis Suite (6 types available):**

**1. Neighborhood Enrichment (`spatial_statistics: neighborhood`):**
- Z-scores showing cell type co-occurrence
- Red = enriched interactions (avoid each other)
- Blue = depleted interactions (cluster together)
- Network-style visualization available

**2. Moran's I Autocorrelation (`spatial_statistics: moran`):**
- Measures spatial clustering of gene expression
- Values near +1 = strong clustering
- Values near -1 = checkerboard pattern  
- Values near 0 = random distribution

**3. Ripley's K Function (`spatial_statistics: ripley`):**
- Tests for clustering vs. dispersion
- K(r) above expected = clustering at distance r
- K(r) below expected = regular spacing
- Multiple distance scales analyzed

**4. Getis-Ord Statistics (`spatial_statistics: getis_ord`):**
- Identifies local hotspots and coldspots
- Red areas = significant hotspots (high-high clusters)
- Blue areas = coldspots (low-low clusters)
- Shows local spatial anomalies

**Key Statistical Findings:**
- **Neurons and astrocytes are co-enriched** (positive Z-score: +3.2)
- **Oligodendrocytes avoid other cell types** (isolation pattern)
- **Strong spatial autocorrelation** for most genes (Moran's I > 0.6)
- **Hotspots concentrated at tissue boundaries** (Getis-Ord)
- **Regular spacing** detected for certain cell types (Ripley's K)

Your tissue shows highly organized, non-random spatial patterns!"

---

## Step 6: Integration and Quality Control Visualizations

### 💬 Ensuring Data Quality

---

**👤 You:** "These analyses look great, but I want to make sure my data integration worked properly. Can you show me quality control visualizations?"

**🤖 ChatSpatial:** "Excellent practice! Quality control is crucial for reliable results. Let me create comprehensive QC visualizations."

*[Integration QC plots appear]*

**🤖 ChatSpatial:** "Here's your integration quality assessment! ✅🔍

**Batch Integration Quality (`batch_integration`):**
- **UMAP by batch**: Shows how well batches mix
- **Gene expression consistency**: Batch-specific artifacts detected
- **Integration metrics**: Quantitative mixing scores
- **Before/after comparison**: Shows integration improvement

**Deconvolution Results (`deconvolution`):**
- **Cell type proportions**: Spatial maps of estimated cell types
- **Confidence scores**: Reliability of predictions
- **Mixed regions**: Areas with multiple cell types
- **Validation metrics**: Agreement with known markers

**Quality Assessment Results:**
- ✅ **Excellent batch mixing**: No obvious batch effects
- ✅ **Consistent gene expression**: No artificial gradients
- ✅ **High deconvolution confidence**: 87% of spots confidently assigned
- ✅ **Biological consistency**: Results match known anatomy
- ⚠️ **Minor edge effects**: Some boundary artifacts detected

**Quality Score: 9.2/10** - Excellent data quality for analysis!"

**👤 You:** "Perfect! Now I want to explore cell communication. Can you show me the different ways to visualize cell-cell interactions?"

**🤖 ChatSpatial:** "Absolutely! Cell communication analysis has several visualization approaches. Let me show you the full suite."

*[Cell communication visualization panel appears]*

**🤖 ChatSpatial:** "Here's your comprehensive cell communication analysis! 🤝💬

**Communication Visualization Types:**

**1. Cell Communication Networks (`cell_communication`):**
- Node size = cell type abundance
- Edge thickness = interaction strength  
- Colors = communication pathways
- Layout optimized for clarity

**2. Ligand-Receptor Pairs (`lr_pairs`):**
- Heatmap of LR interaction strengths
- Sender cell types (rows) × Receiver cell types (columns)
- Color intensity = statistical significance
- Top interacting pairs highlighted

**3. Spatial Interactions (`spatial_interaction`):**
- Distance-dependent communication strength
- Spatial maps showing interaction hotspots
- Decay curves for interaction vs. distance
- Local vs. long-range communication patterns

**4. Pathway Enrichment (`pathway_enrichment`):**
- Active signaling pathways in each region
- GSEA-style enrichment plots
- Pathway activity spatial maps
- Functional annotation of communications

**Communication Insights:**
- **Strongest interactions**: Neuron ↔ Astrocyte (neurotropic signaling)
- **Spatial hotspots**: Communication peaks at domain boundaries  
- **Active pathways**: Wnt signaling, Neurotrophin signaling, ECM interactions
- **Long-range communication**: Minimal beyond 100μm distance
- **Functional relevance**: Networks match known brain connectivity

Your tissue shows biologically relevant communication patterns that align with known neuroscience!"

---

## Step 7: Interactive and Dynamic Visualizations

### 💬 Exploring Interactively

---

**👤 You:** "These static figures are beautiful! But sometimes I want to explore the data interactively. Can you create something I can interact with?"

**🤖 ChatSpatial:** "Great idea! Let me create some interactive visualizations for exploration."

*[Interactive plot interface appears]*

**🤖 ChatSpatial:** "Here's your interactive exploration dashboard! 🖱️

**Interactive Features:**
- **Zoom and pan**: Explore tissue regions in detail
- **Gene selector**: Switch between any gene instantly
- **Hover information**: Get spot details on mouseover
- **Color scale adjustment**: Modify visualization on the fly
- **Region selection**: Click to highlight areas of interest

**How to Use:**
1. **Select a gene** from the dropdown menu
2. **Zoom in** on interesting regions
3. **Hover over spots** to see expression values
4. **Click spots** to get detailed information
5. **Adjust colors** to highlight patterns

This lets you explore your data dynamically and discover patterns that might not be obvious in static plots!"

**👤 You:** "This is amazing for exploration! Can you show me some of the specialized analysis plots? I'm particularly interested in spatial statistics and spatially variable genes."

**🤖 ChatSpatial:** "Excellent choice! Let me demonstrate some of our advanced spatial analysis capabilities."

*[Multiple specialized plots appear]*

**🤖 ChatSpatial:** "Here are several advanced spatial analysis visualizations! 📊🔬

**1. Spatially Variable Genes (SPARK-X/SpatialDE):**
- **Gene identification**: Statistical testing for spatial patterns
- **Spatial correlation**: Quantifies spatial expression relationships
- **Significance testing**: FDR-corrected p-values for reliability

**2. Spatial Statistics (12 sub-types available):**
- **Neighborhood enrichment**: Cell type co-occurrence patterns
- **Moran's I autocorrelation**: Spatial clustering of gene expression
- **Ripley's K function**: Point pattern analysis for cell distributions
- **Getis-Ord statistics**: Local hotspot detection
- **Centrality measures**: Network importance in spatial context
- **Co-occurrence analysis**: Spatial interaction patterns

**3. Communication Networks:**
- **Ligand-receptor pairs**: Interaction strength between cell types
- **Spatial interactions**: Distance-dependent communication
- **Network visualization**: Cell type interaction graphs

**4. Quality Control Plots:**
- **Integration assessment**: Batch effect visualization
- **Deconvolution results**: Cell type proportion maps

**Key Insights from These Plots:**
- Your tissue shows strong spatial organization (high Moran's I)
- 234 spatially variable genes identified (SPARK-X)
- Cell communication is highest at domain boundaries
- Integration successfully removed batch effects
- Neuronal markers show expected spatial patterns

These specialized plots reveal patterns that basic visualizations might miss!"

---

## Common Visualization Requests

### 💬 Popular Plot Types (All 20 Available)

**Basic Visualization:**
- **"Can you create a heatmap of the top genes in each region?"** → Clustered heatmaps with hierarchical organization
- **"Show me a violin plot comparing gene expression between domains"** → Beautiful violin plots with statistical comparisons
- **"Create a spatial map of my tissue"** → Spatial gene expression plots with tissue architecture
- **"Make a UMAP plot colored by cell types"** → UMAP embeddings with annotations

**Advanced Gene Analysis:**
- **"I want to see correlations between different genes"** → Gene correlation matrices and scatter plots
- **"Create a multi-gene panel showing key markers"** → Multi-panel spatial plots for gene comparison
- **"Find spatially variable genes"** → SPARK-X or SpatialDE analysis with statistical significance

**Cell Communication & Interactions:**
- **"Visualize ligand-receptor interactions"** → LR pair analysis and interaction networks
- **"Show me spatial cell-cell interactions"** → Spatial interaction maps and communication networks
- **"Create a cell communication network"** → Network diagrams showing intercellular signaling

**Trajectory & Development:**
- **"Show me RNA velocity patterns"** → Velocity field overlays and developmental arrows
- **"Create trajectory plots for pseudotime"** → Developmental pathways and differentiation routes

**Spatial Analysis (12 Sub-types):**
- **"Perform neighborhood enrichment analysis"** → Co-occurrence and spatial clustering patterns
- **"Show me Moran's I spatial autocorrelation"** → Spatial statistics and pattern detection
- **"Create Ripley's K function plots"** → Point pattern analysis
- **"Analyze centrality measures"** → Network centrality in spatial context
- **"Show Getis-Ord local statistics"** → Hotspot detection and clustering

**Quality Control & Integration:**
- **"Check my integration quality"** → Batch effect assessment and integration QC
- **"Show deconvolution results"** → Cell type proportion maps

**Enrichment Analysis:**
- **"Visualize pathway enrichment"** → GSEA results and pathway activity
- **"Create spatial enrichment maps"** → Functional enrichment in spatial context

### 💬 Customization Options

**Colors:** "Use a blue-to-red color scheme" / "Make it colorblind-friendly"
**Size:** "Make the plot larger" / "Create a small version for supplementary"
**Style:** "Use a minimal theme" / "Make it more colorful and engaging"
**Labels:** "Add more detailed labels" / "Simplify the text"
**Layout:** "Arrange in a 2x2 grid" / "Make one panel dominant"

## Tips for Great Visualizations

### 🎨 Design Principles

**1. Tell a Story**
- Start with overview, then details
- Use consistent color schemes
- Guide the viewer's eye logically

**2. Consider Your Audience**
- **Scientists**: Include technical details and statistics
- **General audience**: Focus on clear, simple messages
- **Journals**: Follow publication guidelines

**3. Accessibility Matters**
- Use colorblind-friendly palettes
- Ensure sufficient contrast
- Provide alternative visual cues beyond color

**4. Quality Standards**
- **Presentations**: 150-300 DPI
- **Publications**: 300+ DPI
- **Web**: Optimized file sizes

### 🚀 Advanced Techniques

**"Create a figure showing before/after treatment comparison"**
**"Make a time-lapse animation of spatial changes"**
**"Generate a poster-sized summary of all analyses"**
**"Create interactive plots for a web dashboard"**

## Next Steps

### 📊 Continue Visualizing

**Ready to create more visualizations?** Try these:
- **"Create publication figures for all my analyses"**
- **"Make an interactive dashboard for data exploration"**
- **"Generate presentation slides with key visualizations"**
- **"Create animated explanations of my results"**

### 📚 Learn More

- [Basic Spatial Analysis](basic_spatial_statistics.md) - Generate data to visualize
- [Cell Type Annotation](../analysis/cell_type_annotation.md) - Create cell type plots
- [Cell Communication Analysis](../analysis/cell_communication_analysis.md) - Visualize interactions
- [Spatial Statistics](spatial_statistics.md) - Explore spatial statistics and pattern analysis
- [API Reference](../../reference/api/README.md) - Complete visualization parameters

## The Complete Visualization Arsenal

### 🎯 All 17 Plot Types at Your Fingertips

ChatSpatial provides comprehensive visualization capabilities:

**Core Spatial Plots (4 types):**
1. `spatial` - Gene expression in tissue context
2. `spatial_domains` - Tissue architecture and regions
3. `umap` - Dimensionality reduction with annotations
4. `heatmap` - Multi-gene expression matrices

**Statistical Analysis (4 types):**
5. `violin` - Distribution comparisons
6. `gene_correlation` - Co-expression analysis
7. `spatial_statistics` - 12 statistical sub-types (neighborhood, Moran's I, Ripley's K, etc.)
8. `spatial_enrichment` - Functional enrichment mapping

**Cell Communication (4 types):**
9. `cell_communication` - Interaction networks
10. `lr_pairs` - Ligand-receptor analysis
11. `spatial_interaction` - Distance-dependent communication
12. `multi_gene` - Multi-panel gene comparisons

**Trajectory & Dynamics (2 types):**
13. `rna_velocity` - Velocity field visualization
14. `trajectory` - Developmental pathways

**Quality Control (3 types):**
15. `deconvolution` - Cell type proportions
16. `pathway_enrichment` - GSEA and pathway activity
17. `batch_integration` - Batch correction assessment

### 🚀 Pro Tips for Maximum Impact

**For Presentations:**
- Use `multi_gene` plots to compare key markers
- Combine `spatial_domains` with `cell_communication` for complete story
- Include `batch_integration` to show data quality

**For Publications:**
- Start with `spatial` overview, then zoom into specific analyses
- Use `find_spatial_genes` (SPARK-X/SpatialDE) for novel spatial gene discovery
- Include statistical plots (`spatial_statistics`) for rigor

**For Exploration:**
- Try `gene_correlation` to find co-expressed modules
- Use `spatial_interaction` to discover communication hotspots
- Experiment with different `spatial_statistics` sub-types

**Remember**: Great visualizations don't just show data - they reveal insights and tell compelling stories about your tissue's biology!
