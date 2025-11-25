# Examples

All examples are **natural language commands** you type in Claude chat. No coding required!

---

## Basic Spatial Analysis

### Load and Explore Data
```text
"Load my Visium dataset from /path/to/data.h5ad"
"Show me basic statistics and quality metrics"
"What genes are most expressed?"
"Create a spatial plot of the tissue"
```

### Quality Control
```text
"Filter spots with fewer than 200 genes"
"Remove mitochondrial genes"
"Show me quality control metrics"
```

### Preprocessing
```text
"Normalize and preprocess the data"
"Find the top 2000 highly variable genes"
"Perform PCA and UMAP dimensionality reduction"
```

### Clustering
```text
"Cluster spots using Leiden algorithm"
"Try different clustering resolutions"
"Show me the clusters on spatial coordinates"
"Find marker genes for each cluster"
```

---

## Spatial Domain Analysis

### Identify Tissue Regions
```text
"Identify spatial domains using SpaGCN"
"Find distinct tissue regions with STAGATE"
"Show me the spatial domain boundaries"
```

### Compare Domains
```text
"What are the differences between domain 1 and domain 3?"
"Find genes that distinguish different spatial regions"
"Create a heatmap comparing all domains"
```

### Spatial Patterns
```text
"Find spatially variable genes"
"Detect spatial hotspots using Getis-Ord"
"Analyze spatial autocorrelation with Moran's I"
```

---

## Cell Type Analysis

### Cell Type Annotation
```text
"Annotate cell types using Tangram with reference data"
"Use scANVI to transfer labels from reference to spatial"
"Identify cell types with CellAssign using marker genes"
"Map cell type distributions across the tissue"
```

### Cell Type Markers
```text
"Find marker genes for each cell type"
"Show me the top markers for T cells"
"Create a dotplot of cell type markers"
```

---

## Deconvolution Analysis

### Spatial Deconvolution
```text
"Deconvolve spatial spots using Cell2location"
"Estimate cell type proportions with RCTD"
"Use Tangram to map single cells to spatial locations"
```

### Analyze Composition
```text
"Show me the dominant cell type in each spot"
"Calculate cell type diversity across tissue regions"
"Visualize cell type proportions as pie charts"
"Find spatial patterns in cell type composition"
```

---

## Cell Communication Analysis

### Identify Interactions
```text
"Analyze cell-cell communication using LIANA"
"Find ligand-receptor pairs between cell types"
"Detect spatial communication patterns"
```

### Visualize Networks
```text
"Show me communication networks between cell types"
"Create a chord diagram of interactions"
"Which cell types communicate most?"
```

### Spatial Context
```text
"Analyze communication in spatial neighborhoods"
"Find interactions enriched at tissue boundaries"
"Compare communication patterns across regions"
```

---

## Advanced Analysis

### RNA Velocity & Trajectories
```text
"Analyze RNA velocity with scVelo"
"Infer developmental trajectories using CellRank"
"Find trajectory endpoints and branch points"
"Overlay velocity vectors on spatial coordinates"
```

### Pathway Enrichment
```text
"Run pathway enrichment analysis on marker genes"
"Find enriched GO terms in spatial domains"
"Perform spatial enrichment mapping"
"Compare pathway activities across regions"
```

### Spatial Statistics
```text
"Calculate neighborhood enrichment between cell types"
"Analyze co-occurrence patterns"
"Compute Ripley's K function for point patterns"
"Find spatially autocorrelated genes"
```

### Copy Number Variation
```text
"Analyze copy number variations using infercnvpy"
"Detect CNV patterns in malignant cells"
"Compare CNV profiles across spatial regions"
```

---

## Multi-Sample Analysis

### Batch Integration
```text
"Integrate multiple Visium samples using Harmony"
"Correct batch effects with scVI"
"Align samples using PASTE spatial registration"
```

### Comparative Analysis
```text
"Compare spatial domains across samples"
"Find conserved marker genes across batches"
"Identify sample-specific spatial patterns"
```

---

## Visualization

### Spatial Plots
```text
"Create spatial plot colored by gene expression"
"Show multiple genes side by side"
"Overlay cell types on tissue image"
"Add cluster outlines to spatial plot"
```

### Expression Heatmaps
```text
"Generate heatmap of top variable genes"
"Show marker gene expression across clusters"
"Create annotated heatmap with cell type labels"
```

### Dimension Reduction
```text
"Plot UMAP colored by clusters"
"Show spatial coordinates in low-dimensional space"
"Create interactive UMAP with gene expression"
```

### Specialized Plots
```text
"Visualize deconvolution results as scatterpie"
"Create stacked bar chart of cell type proportions"
"Show communication network as chord diagram"
"Plot trajectory with pseudotime gradient"
```

---

## Complete Workflows

### Workflow 1: Basic Spatial Analysis (5 minutes)
```text
1. "Load /path/to/visium_data.h5ad"
2. "Preprocess and cluster the data"
3. "Identify spatial domains using SpaGCN"
4. "Find marker genes for each domain"
5. "Create visualizations showing spatial patterns"
```

### Workflow 2: Cell Type Analysis (10 minutes)
```text
1. "Load spatial data and reference scRNA-seq"
2. "Annotate cell types using Tangram"
3. "Analyze spatial distribution of each cell type"
4. "Find cell type-specific marker genes"
5. "Visualize results with spatial plots and heatmaps"
```

### Workflow 3: Deconvolution + Communication (15 minutes)
```text
1. "Load spatial and reference data"
2. "Deconvolve spatial data with Cell2location"
3. "Analyze cell-cell communication using LIANA"
4. "Find spatially variable genes"
5. "Create comprehensive visualization of results"
```

### Workflow 4: Advanced Multi-Method (20 minutes)
```text
1. "Load and preprocess spatial data"
2. "Identify spatial domains with SpaGCN"
3. "Annotate cell types with scANVI"
4. "Deconvolve spots using RCTD"
5. "Analyze cell communication patterns"
6. "Infer trajectories with CellRank"
7. "Run pathway enrichment on all results"
8. "Generate publication-ready figures"
```

---

## Tips for Best Results

### Data Loading
- Always use absolute paths: `/Users/name/data.h5ad`
- Verify spatial coordinates exist in data
- Check data quality before analysis

### Analysis Strategy
- Start with preprocessing and clustering
- Try multiple methods for robust results
- Use reference data when available
- Validate findings with marker genes

### Asking Questions
- Be specific about what you want to see
- Mention method names if you have preferences
- Ask for visualizations to understand results
- Request comparisons between conditions

### Troubleshooting
- If analysis fails, try preprocessing first
- Check data format and required fields
- See [Troubleshooting Guide](advanced/troubleshooting.md)
- Ask Claude to explain errors

---

## Related Resources

- **[Quick Start](quickstart.md)** - Get ChatSpatial running in 5 minutes
- **[Methods Reference](advanced/methods-reference.md)** - Complete tool documentation
- **[FAQ](advanced/faq.md)** - Common questions and answers
- **[GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)** - Report problems or request features

---

**Remember**: ChatSpatial understands natural language. Don't memorize commands - just describe what you want to do!
