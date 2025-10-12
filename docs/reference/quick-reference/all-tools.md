# 🔧 ChatSpatial MCP Tools Quick Reference

> **Quick Link**: Keep this page bookmarked for instant access to all 15 MCP tools.

## 📋 Tool Categories Overview

| Category | Tools | Use Cases |
|----------|-------|-----------|
| **Data Management** | `load_data`, `preprocess_data` | Loading and preparing datasets |
| **Analysis** | `analyze_spatial_data`, `find_markers`, `find_spatial_genes` | Core spatial analysis |
| **Annotation** | `annotate_cells` | Cell type identification |
| **Advanced Analysis** | `analyze_velocity_data`, `analyze_trajectory_data`, `analyze_cell_communication`, `analyze_enrichment` | RNA velocity, trajectories, communication |
| **Integration** | `integrate_samples`, `register_spatial_data` | Multi-sample and alignment |
| **Deconvolution** | `deconvolve_data`, `identify_spatial_domains` | Spot deconvolution and domain identification |
| **Visualization** | `visualize_data` | All plotting and visualization |

**Total: 14 Core MCP Tools** covering the complete spatial transcriptomics analysis pipeline.

---

## 🔵 Data Management Tools

### `load_data`
**Purpose**: Load spatial transcriptomics data from various formats  
**Difficulty**: 🟢 Beginner

**Key Parameters**:
- `data_path`: Path to data file/directory
- `data_type`: `"auto"`, `"10x_visium"`, `"slide_seq"`, `"merfish"`, `"seqfish"`, `"h5ad"`
- `name`: Optional dataset name

**Example Queries**:
```
"Load my Visium data from /path/to/data"
"Import the 10x Genomics dataset"
"Load multiple H5AD files"
```

**Returns**: Dataset information with ID for further analysis

---

### `preprocess_data`
**Purpose**: Normalize, filter, and prepare data for analysis  
**Difficulty**: 🟡 Intermediate

**Key Parameters**:
- `data_id`: Dataset ID from load_data
- `params.normalization_method`: `"log"`, `"sct"`, `"pearson_residuals"`, `"scvi"`
- `params.n_top_genes`: Number of highly variable genes (default: 3000)
- `params.clustering_resolution`: Leiden clustering resolution (auto-adaptive)

**Example Queries**:
```
"Preprocess the data with SCTransform normalization"
"Filter cells and normalize using log transformation"  
"Use scVI for advanced preprocessing"
```

**Returns**: Preprocessed data with QC metrics

---

## 🔴 Core Analysis Tools

### `analyze_spatial_data`
**Purpose**: Analyze spatial patterns and relationships  
**Difficulty**: 🟡 Intermediate

**Key Parameters**:
- `data_id`: Dataset ID
- `params.analysis_type`: 
  - `"moran"` - Global Moran's I
  - `"local_moran"` - Local Moran's I (LISA) 
  - `"geary"` - Geary's C
  - `"getis_ord"` - Getis-Ord Gi*
  - `"neighborhood"` - Neighborhood enrichment
  - `"co_occurrence"` - Co-occurrence patterns
  - `"ripley"` - Ripley's K/L
  - `"centrality"` - Graph centrality
  - `"bivariate_moran"` - Gene pair correlation
  - `"join_count"` - Categorical autocorrelation
  - `"network_properties"` - Network analysis
  - `"spatial_centrality"` - Spatial centrality
- `params.genes`: List of genes to analyze (unified parameter)
- `params.n_neighbors`: Number of spatial neighbors (default: 6)

**Example Queries**:
```
"Analyze spatial autocorrelation for CD3D"
"Find local spatial hotspots for these genes"
"Calculate Geary's C for my marker genes"
"Find spatial hotspots using Getis-Ord"
"Calculate neighborhood enrichment"
"Analyze bivariate spatial correlation between gene pairs"
```

**Returns**: Spatial statistics and significance values

---

### `find_markers`
**Purpose**: Identify differentially expressed genes between groups  
**Difficulty**: 🟢 Beginner

**Key Parameters**:
- `data_id`: Dataset ID
- `group_key`: Column defining groups (e.g., "leiden")
- `group1`/`group2`: Specific groups to compare
- `method`: `"wilcoxon"`, `"t-test"`, `"logreg"`
- `n_top_genes`: Top DEGs to return (default: 25)

**Example Queries**:
```
"Find marker genes for cluster 0"
"Compare tumor vs normal regions"
"Identify top 50 differential genes"
```

**Returns**: Ranked list of marker genes with statistics

---

### `find_spatial_genes`
**Purpose**: Identify spatially variable genes  
**Difficulty**: 🔴 Advanced

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"gaston"`, `"spatialde"`, `"spark"`
- `params.n_top_genes`: Number of top SVGs to return
- `params.epochs`: Training epochs for GASTON (default: 10000)

**Example Queries**:
```
"Find spatially variable genes using GASTON"
"Identify spatial patterns with SpatialDE"
"Use SPARK to find spatial genes"
```

**Returns**: Ranked spatially variable genes

---


## 🟣 Cell Annotation Tools

### `annotate_cells`
**Purpose**: Identify cell types using multiple methods  
**Difficulty**: 🟡 Intermediate

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"marker_genes"`, `"tangram"`, `"scanvi"`, `"cellassign"`, `"sctype"`, `"mllmcelltype"`
- `params.reference_data_id`: Reference dataset ID (for transfer methods)
- `params.marker_genes`: Custom marker gene lists

**Example Queries**:
```
"Annotate cells using marker genes"
"Use Tangram with single-cell reference"
"Apply scType for automatic annotation"
```

**Returns**: Cell type assignments and confidence scores

---

## 🟠 Advanced Analysis Tools

### `analyze_velocity_data`
**Purpose**: RNA velocity analysis for cellular dynamics  
**Difficulty**: 🔴 Advanced

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"scvelo"` (standard), `"velovi"` (deep learning), `"sirv"` (reference-based)
- `params.mode`: `"dynamical"`, `"stochastic"` (for scVelo)
- `params.n_top_genes`: Velocity genes to use
- `params.velovi_n_epochs`: Training epochs (for VELOVI)

**Example Queries**:
```
"Analyze RNA velocity with scVelo"
"Use VELOVI for deep learning velocity computation"
"Calculate velocity with stochastic mode"
```

**Returns**: Velocity vectors and velocity graph

**Note**: This computes RNA velocity. For trajectory inference, use `analyze_trajectory_data` afterwards.

---

### `analyze_trajectory_data`
**Purpose**: Infer cellular trajectories and pseudotime  
**Difficulty**: 🔴 Advanced

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"dpt"` (diffusion), `"palantir"` (probabilistic), `"cellrank"` (velocity-based)
- `params.root_cells`: Starting cells for trajectory
- `params.spatial_weight`: Weight for spatial information

**Example Queries**:
```
"Infer differentiation trajectories"
"Calculate pseudotime with Palantir"
"Use CellRank for RNA velocity-based trajectories"
```

**Returns**: Pseudotime values and trajectory paths

**Note**: CellRank requires velocity data from `analyze_velocity_data`. Palantir and DPT work without velocity.

---

### `analyze_cell_communication`
**Purpose**: Cell-cell communication analysis  
**Difficulty**: 🔴 Advanced

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"liana"`, `"cellphonedb"`, `"cellchat_liana"`
- `params.species`: `"human"`, `"mouse"`, `"zebrafish"`
- `params.perform_spatial_analysis`: Spatial vs cluster-based analysis

**Example Queries**:
```
"Analyze ligand-receptor interactions with LIANA"
"Find spatial cell communication patterns"
"Identify tumor-immune interactions"
```

**Returns**: Significant ligand-receptor pairs and spatial patterns

---

### `analyze_enrichment`
**Purpose**: Gene set enrichment analysis  
**Difficulty**: 🟡 Intermediate

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"spatial_enrichmap"`, `"pathway_gsea"`, `"pathway_ora"`, `"pathway_enrichr"`
- `params.gene_set_database`: `"GO_Biological_Process"`, `"KEGG"`, `"Reactome"`
- `params.pvalue_cutoff`: Significance threshold (default: 0.05)

**Example Queries**:
```
"Perform spatial GSEA analysis"
"Find enriched pathways in tumor regions"
"Analyze GO term enrichment"
```

**Returns**: Enriched gene sets with statistics

---

## 🟡 Integration & Registration Tools

### `integrate_samples`
**Purpose**: Integrate multiple spatial samples  
**Difficulty**: 🔴 Advanced

**Key Parameters**:
- `data_ids`: List of dataset IDs
- `params.method`: `"harmony"`, `"bbknn"`, `"scanorama"`, `"scvi"`
- `params.batch_key`: Batch identifier column
- `params.n_pcs`: Principal components to use

**Example Queries**:
```
"Integrate three Visium samples with Harmony"
"Remove batch effects across datasets"
"Combine samples from different patients"
```

**Returns**: Integrated dataset with corrected embeddings

---

### `register_spatial_data`
**Purpose**: Align spatial sections or time points
**Difficulty**: 🔴 Advanced
**Status**: ⚠️ **IN DEVELOPMENT** - Core registration methods not yet implemented

**Key Parameters**:
- `source_id`: Source dataset ID
- `target_id`: Target dataset ID
- `method`: `"paste"`, `"stalign"` _(planned, not yet available)_
- `landmarks`: Optional alignment landmarks

**Example Queries**:
```
"Align consecutive tissue sections"
"Register brain slices for 3D reconstruction"
"Align samples to reference coordinates"
```

**Note**: While the tool interface exists, PASTE and STalign implementations are pending. Use alternative alignment tools or manual registration for now.

**Returns**: Transformation matrix and aligned coordinates _(when implemented)_

---

## 🟤 Deconvolution & Domain Tools

### `deconvolve_data`
**Purpose**: Deconvolve spots into cell type proportions  
**Difficulty**: 🔴 Advanced

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"cell2location"`, `"destvi"`, `"stereoscope"`, `"tangram"`, `"rctd"`, `"spotlight"`
- `params.reference_data_id`: Single-cell reference ID
- `params.n_epochs`: Training iterations

**Example Queries**:
```
"Deconvolve Visium spots with Cell2location"
"Estimate cell type proportions using RCTD"
"Use scRNA reference for deconvolution"
```

**Returns**: Cell type proportion estimates per spot

---

### `identify_spatial_domains`
**Purpose**: Find tissue domains and spatial niches  
**Difficulty**: 🟡 Intermediate

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"spagcn"`, `"leiden"`, `"louvain"`, `"stagate"`
- `params.n_domains`: Number of expected domains (default: 7)
- `params.resolution`: Clustering resolution

**Example Queries**:
```
"Identify spatial domains with SpaGCN"
"Find tissue niches using STAGATE"
"Cluster spots into spatial regions"
```

**Returns**: Domain assignments and spatial boundaries

---

## 🎨 Visualization Tools

### `visualize_data`
**Purpose**: Create all types of spatial plots and visualizations  
**Difficulty**: 🟢 Beginner

**Key Parameters**:
- `data_id`: Dataset ID
- `params.plot_type`: **15+ plot types available** (see below)
- `params.feature`: Gene(s) or metadata to visualize
- `params.colormap`: Color scheme (`"viridis"`, `"plasma"`, `"RdYlBu"`)

#### Available Plot Types:

**Spatial Plots**:
- `"spatial"` - Basic spatial gene expression
- `"spatial_domains"` - Spatial domain boundaries
- `"spatial_pie"` - Cell type proportion pie charts
- `"spatial_contour"` - Expression contour maps

**Dimensionality Reduction**:
- `"umap"` - UMAP embedding
- `"tsne"` - t-SNE embedding  
- `"pca"` - Principal component analysis

**Statistical Plots**:
- `"violin"` - Violin plots by group
- `"heatmap"` - Expression heatmaps
- `"dotplot"` - Dot plots for multiple genes
- `"matrixplot"` - Matrix-style heatmaps

**Analysis-Specific**:
- `"trajectory"` - Trajectory/pseudotime plots
- `"cell_communication"` - Ligand-receptor networks
- `"deconvolution"` - Cell type proportions
- `"spatial_analysis"` - Spatial statistics results
- `"gaston_domains"` - GASTON spatial domains
- `"gaston_genes"` - GASTON gene classifications

**Example Queries**:
```
"Show spatial expression of CD3D"
"Create UMAP plot colored by cell types"
"Plot spatial domains with boundaries"
"Generate trajectory visualization"
"Show ligand-receptor interaction network"
```

**Returns**: High-resolution image for display or saving

---

## 🚀 Quick Start Combinations

### Basic Analysis Workflow
```
1. load_data → 2. preprocess_data → 3. identify_spatial_domains → 4. visualize_data
```

### Cell Type Analysis
```
1. load_data → 2. preprocess_data → 3. annotate_cells → 4. find_markers → 5. visualize_data
```

### Communication Analysis  
```
1. load_data → 2. preprocess_data → 3. annotate_cells → 4. analyze_cell_communication → 5. visualize_data
```

### Multi-Sample Integration
```
1. load_data (×N) → 2. integrate_samples → 3. preprocess_data → 4. analyze_* → 5. visualize_data
```

---

## 💡 Pro Tips

### Natural Language Examples
- **Specific**: "Find spatial domains in my mouse brain data using SpaGCN"
- **General**: "Analyze the spatial patterns in this dataset"  
- **Complex**: "Compare cell communication between tumor core and invasive front"

### Parameter Guidelines
- Start with default parameters for initial exploration
- Most tools auto-adapt to your dataset size
- Use `context.info()` messages to understand what's happening
- Check return values for quality metrics and suggestions

### Memory Management  
- Large datasets (>50K spots): Consider subsampling for exploration
- GPU acceleration available for: scVI methods, STAGATE, deep learning tools
- Use sparse matrices automatically handled by the tools

### Result Resources
All analysis results are automatically saved as MCP resources:
- `spatial://datasets/{data_id}` - Dataset information
- `spatial://results/{data_id}/{analysis_type}` - Analysis results  
- `spatial://visualizations/{viz_id}` - Generated plots

---

## 🔍 Parameter Quick Reference

### Universal Parameters
These parameters are available across multiple tools:

| Parameter | Type | Default | Description | Used In |
|-----------|------|---------|-------------|----------|
| `data_id` | str | Required | Dataset identifier | All analysis tools |
| `context` | Context | Auto | MCP context object | All tools |
| `method` | str | Tool-specific | Analysis method to use | Most analysis tools |
| `n_neighbors` | int | 6 | Spatial neighbors count | Spatial analysis |
| `resolution` | float | 0.5 | Clustering resolution | Domain/clustering |
| `n_top_genes` | int | 2000-3000 | Top variable genes | Multiple tools |
| `random_state` | int | 42 | Random seed | Reproducible analyses |

### Data Type Support Matrix

| Tool | 10X Visium | Slide-seq | MERFISH | seqFISH | H5AD | Notes |
|------|------------|-----------|---------|---------|------|-------|
| `load_data` | ✅ | ✅ | ✅ | ✅ | ✅ | Universal loader |
| `preprocess_data` | ✅ | ✅ | ✅ | ✅ | ✅ | All spatial formats |
| `identify_spatial_domains` | ✅ | ✅ | ⚠️ | ⚠️ | ✅ | Best with continuous coordinates |
| `deconvolve_data` | ✅ | ✅ | ❌ | ❌ | ✅ | Requires spot-level data |
| `analyze_cell_communication` | ✅ | ✅ | ✅ | ✅ | ✅ | Universal |

### Method Compatibility

**Spatial Domain Methods:**
- `spagcn`: Best with histology images ⭐
- `stagate`: High-resolution, GPU recommended  
- `leiden`: Fast, CPU-friendly clustering
- `louvain`: Alternative clustering method

**Cell Annotation Methods:**
- `marker_genes`: Fast, requires known markers
- `tangram`: Best with scRNA reference ⭐
- `sctype`: Automatic, good for exploration
- `cellassign`: Probabilistic, custom markers

**Deconvolution Methods:**
- `cell2location`: GPU recommended, high accuracy ⭐
- `rctd`: CPU-friendly, robust
- `destvi`: Deep learning, experimental
- `stereoscope`: Fast, good baseline

### Performance Guidelines

| Dataset Size | Recommended Tools | Memory | Time | Notes |
|-------------|------------------|---------|-------|-------|
| Small (<5K spots) | Any method | 4GB+ | 5-15 min | Explore all options |
| Medium (5-20K spots) | Avoid GASTON, use SpaGCN | 8GB+ | 15-45 min | Standard workflows |
| Large (20-50K spots) | Leiden/Louvain clustering | 16GB+ | 30-90 min | Consider subsampling |
| XL (50K+ spots) | Leiden, chunked analysis | 32GB+ | 1-3 hours | Definitely subsample first |

---

## 🚀 Natural Language Query Examples

### 🎯 Beginner Queries
```text
"Load my Visium data and show me what's in it"
"Find spatial domains in my brain tissue"
"Annotate cell types using marker genes"
"Create a spatial plot of CD3D expression"
"Show me quality control metrics"
```

### 🎯 Intermediate Queries
```text
"Compare spatial domains using SpaGCN vs STAGATE"
"Deconvolve spots with Cell2location using my scRNA reference"
"Find ligand-receptor interactions between tumor and immune cells"
"Integrate three samples with batch correction"
"Perform pathway enrichment on spatially variable genes"
```

### 🎯 Advanced Queries
```text
"Build complete analysis pipeline: domains → annotation → communication → pathways"
"Compare cell communication networks between treated vs control samples"
"Register spatial sections and analyze developmental trajectories"
"Validate spatial domain results using multiple methods with statistical comparison"
"Create publication-ready figure panels with spatial domains, cell types, and communication networks"
```

### 🎯 Troubleshooting Queries
```text
"My clustering produced 50+ domains, how do I optimize parameters?"
"Cell2location failed with CUDA error, suggest CPU alternatives"
"Show me what analysis steps have been completed so far"
"My dataset is 100K spots, suggest memory-efficient workflow"
"Visualizations are not showing up, diagnose the issue"
```

---

**⚡ Remember**: This is your complete toolkit for spatial transcriptomics analysis. Each tool is designed to work together seamlessly through natural language interaction. Start with simple queries and build complexity as needed!