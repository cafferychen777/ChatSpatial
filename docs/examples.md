# Examples

Natural language commands for spatial transcriptomics analysis.

---

## Standard Workflow

A typical analysis follows this flow:

```
Load → Preprocess → Analyze → Visualize
```

### 1. Load Your Data

```
"Load /path/to/spatial_data.h5ad"
"Load my Visium data from /path/to/visium_folder"
```

### 2. Preprocess

```
"Preprocess the data"
"Normalize with log transformation and find 2000 variable genes"
"Preprocess with SCTransform normalization"
```

### 3. Analyze

Choose your analysis type below.

### 4. Visualize

```
"Show the spatial plot"
"Visualize CD3D expression on tissue"
"Create a UMAP colored by clusters"
```

---

## Analysis Types

### Spatial Domains

Identify tissue regions and niches.

```
"Identify spatial domains"
"Find 7 spatial domains using SpaGCN"
"Cluster the tissue into regions with STAGATE"
"Use Leiden clustering with resolution 0.5"
```

**Methods**: SpaGCN (default), STAGATE, GraphST, Leiden, Louvain

---

### Cell Type Annotation

Assign cell types to spots or cells.

```
"Annotate cell types using the reference dataset"
"Transfer labels from reference with Tangram"
"Use scANVI for label transfer"
"Annotate with marker genes using CellAssign"
```

**Methods**: Tangram, scANVI, CellAssign, scType, SingleR, mLLMCelltype

**Requires**: Reference dataset with cell type labels (for transfer methods)

---

### Deconvolution

Estimate cell type proportions in each spot.

```
"Deconvolve the spatial data"
"Estimate cell type proportions with FlashDeconv"
"Use Cell2location for deconvolution"
"Run RCTD deconvolution"
```

**Methods**: FlashDeconv (default, fastest), Cell2location, RCTD, DestVI, Stereoscope, Tangram, SPOTlight, CARD

**Requires**: Reference single-cell dataset with cell type annotations

---

### Spatial Statistics

Analyze spatial patterns and autocorrelation.

```
"Analyze spatial autocorrelation"
"Calculate Moran's I for marker genes"
"Find spatial hotspots with Getis-Ord"
"Compute neighborhood enrichment"
"Analyze co-occurrence of cell types"
```

**Methods**: Moran's I, Local Moran's I, Geary's C, Getis-Ord Gi*, Ripley's K, neighborhood enrichment, co-occurrence

---

### Spatially Variable Genes

Find genes with spatial expression patterns.

```
"Find spatially variable genes"
"Identify spatial genes with SpatialDE"
"Use SPARK-X to find spatial patterns"
```

**Methods**: SPARK-X (default, fast), SpatialDE

---

### Differential Expression

Compare gene expression between groups.

```
"Find marker genes for cluster 0"
"Compare gene expression between tumor and normal"
"Find differentially expressed genes in domain 3"
```

---

### Condition Comparison

Compare experimental conditions with proper statistics.

```
"Compare treatment vs control across patients"
"Find genes differentially expressed between conditions"
"Analyze condition effects stratified by cell type"
```

**Requires**: Sample/patient identifiers for pseudobulk analysis

---

### Cell Communication

Analyze ligand-receptor interactions.

```
"Analyze cell-cell communication"
"Find ligand-receptor interactions with LIANA"
"Identify spatial communication patterns"
"Which cell types are communicating?"
```

**Methods**: FastCCC (default, fastest), LIANA, CellPhoneDB, CellChat

**Requires**: Cell type annotations

---

### RNA Velocity

Understand cellular dynamics.

```
"Analyze RNA velocity"
"Run scVelo velocity analysis"
"Use VeloVI for velocity estimation"
```

**Methods**: scVelo (deterministic/stochastic/dynamical), VeloVI

**Requires**: Spliced and unspliced count layers

---

### Trajectory Analysis

Infer developmental trajectories.

```
"Infer cellular trajectories"
"Calculate pseudotime with Palantir"
"Use CellRank for fate mapping"
"Compute diffusion pseudotime"
```

**Methods**: CellRank (requires velocity), Palantir, DPT

---

### Pathway Enrichment

Find enriched biological pathways.

```
"Perform pathway enrichment analysis"
"Find enriched GO terms"
"Analyze KEGG pathway enrichment"
"Run GSEA on marker genes"
```

**Methods**: ORA (default), GSEA, ssGSEA, Enrichr

---

### CNV Analysis

Detect copy number variations.

```
"Detect copy number variations"
"Analyze CNV using immune cells as reference"
"Find chromosomal alterations in tumor cells"
```

**Methods**: inferCNVpy (default), Numbat

**Requires**: Normal cell types as reference

---

### Multi-Sample Integration

Combine multiple datasets.

```
"Integrate these three samples"
"Remove batch effects with Harmony"
"Combine datasets using scVI"
```

**Methods**: Harmony (default), BBKNN, Scanorama, scVI

---

### Spatial Registration

Align tissue sections.

```
"Align these two tissue sections"
"Register spatial slices for 3D reconstruction"
```

**Methods**: PASTE, STalign

---

## Visualization Examples

### Basic Plots

```
"Show spatial expression of CD3D"
"Create UMAP plot"
"Plot violin of marker genes by cluster"
"Generate heatmap of top markers"
```

### Deconvolution Results

```
"Show cell type proportions on tissue"
"Create pie charts of cell composition"
"Visualize dominant cell type per spot"
```

### Communication Results

```
"Show ligand-receptor dotplot"
"Visualize communication network"
"Plot top interacting cell types"
```

### Spatial Statistics

```
"Show neighborhood enrichment heatmap"
"Visualize spatial hotspots"
"Plot Moran's I results"
```

---

## Complete Workflows

### Basic Spatial Analysis (5 min)

```
1. "Load /path/to/visium_data.h5ad"
2. "Preprocess the data"
3. "Identify spatial domains"
4. "Find marker genes for each domain"
5. "Visualize the domains on tissue"
```

### Deconvolution Workflow (10 min)

```
1. "Load spatial data from /path/to/spatial.h5ad"
2. "Load reference data from /path/to/reference.h5ad"
3. "Preprocess both datasets"
4. "Deconvolve using the reference"
5. "Show cell type proportions on tissue"
```

### Cell Communication Workflow (10 min)

```
1. "Load the spatial data"
2. "Preprocess with clustering"
3. "Annotate cell types" (or use existing annotations)
4. "Analyze cell-cell communication"
5. "Show the communication network"
```

### Trajectory Workflow (15 min)

```
1. "Load the data"
2. "Preprocess the data"
3. "Analyze RNA velocity"
4. "Infer trajectories with CellRank"
5. "Visualize velocity streams on tissue"
```

---

## Tips

**Be specific when needed**
- General: "Analyze the data" → ChatSpatial chooses defaults
- Specific: "Use SpaGCN with 7 domains" → ChatSpatial uses your settings

**Chain commands naturally**
- "Load the data, preprocess it, and identify spatial domains"

**Reference previous results**
- "Find markers for the domains we just identified"
- "Visualize the deconvolution results"

**Ask for help**
- "What methods are available for deconvolution?"
- "How should I preprocess my data?"
