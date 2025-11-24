# 🔧 ChatSpatial MCP Tools Quick Reference

> **Quick Link**: Keep this page bookmarked for instant access to all 15 MCP tools.

## 📋 Tool Categories Overview

| Category | Tools | Use Cases |
|----------|-------|-----------|
| **Data Management** | `load_data`, `preprocess_data` | Loading and preparing datasets |
| **Analysis** | `analyze_spatial_data`, `find_markers`, `find_spatial_genes` | Core spatial analysis |
| **Annotation** | `annotate_cell_types` | Cell type identification |
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
- `data_type`: `"auto"` (default), `"10x_visium"`, `"slide_seq"`, `"merfish"`, `"seqfish"`, `"h5ad"`, `"other"`
- `name`: Optional dataset name

**Supported Formats**:
- **10x Visium**: Spatial Gene Expression directories (with `spatial/` folder)
- **H5AD**: Pre-processed AnnData files with spatial coordinates in `obsm['spatial']`
- **Slide-seq/MERFISH/seqFISH**: Custom formats with coordinate files
- **Other**: Generic spatial data (requires manual spatial coordinate specification)

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
- `params.normalization`: `"log"`, `"pearson_residuals"`, `"none"`
- `params.n_hvgs`: Number of highly variable genes (default: 2000)
- `params.n_pcs`: Number of principal components (default: 30)
- `params.n_neighbors`: Number of neighbors for graph (default: 15)
- `params.clustering_resolution`: Leiden clustering resolution (default: 1.0)
- `params.use_scvi_preprocessing`: Use scVI for batch correction (default: False)

**Normalization Options**:
- **log** (default): Standard log(x+1) normalization - robust and widely used
- **pearson_residuals**: GLM-based variance stabilization for UMI data - best for single-cell resolution
- **none**: Skip normalization (use when data is already normalized)
- ⚠️ **sct/scvi**: NOT directly available - use `use_scvi_preprocessing=True` for scVI

**Important**:
- Raw count data recommended for pearson_residuals
- For batch correction: use `use_scvi_preprocessing=True`
- Automatically saves `.raw` for downstream analyses requiring full gene set

**Example Queries**:
```
"Preprocess the data with log normalization"
"Filter cells and normalize using Pearson residuals"
"Use scVI for batch correction preprocessing"
```

**Returns**: Preprocessed data with QC metrics

---

## 🔴 Core Analysis Tools

### `analyze_spatial_data`
**Purpose**: Analyze spatial patterns and relationships  
**Difficulty**: 🟡 Intermediate

**Key Parameters**:
- `data_id`: Dataset ID
- `params.analysis_type`: **13 spatial analysis types** (default: `"neighborhood"`)
  - `"moran"` - Global Moran's I (gene-based)
  - `"local_moran"` - Local Moran's I / LISA (gene-based)
  - `"geary"` - Geary's C autocorrelation (gene-based)
  - `"getis_ord"` - Getis-Ord Gi* hotspot detection (gene-based)
  - `"neighborhood"` ⭐ - Neighborhood enrichment (requires `cluster_key`)
  - `"co_occurrence"` - Co-occurrence patterns (requires `cluster_key`)
  - `"ripley"` - Ripley's K/L point patterns (requires `cluster_key`)
  - `"centrality"` - Graph centrality (optional `cluster_key`)
  - `"bivariate_moran"` - Gene pair spatial correlation (gene-based)
  - `"join_count"` - Binary categorical autocorrelation (requires `cluster_key` with 2 categories)
  - `"local_join_count"` - Multi-category local join count (requires `cluster_key` with >2 categories)
  - `"network_properties"` - Network structure analysis (optional `cluster_key`)
  - `"spatial_centrality"` - Spatial importance measures (optional `cluster_key`)
- `params.cluster_key`: **REQUIRED** for group-based analyses (neighborhood, co_occurrence, ripley, join_count, local_join_count)
- `params.genes`: List of specific genes to analyze (for gene-based analyses)
- `params.n_neighbors`: Number of spatial neighbors (default: 8)
- `params.n_top_genes`: Top HVGs to analyze if genes not specified (default: 20)

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
- `params.method`: `"spatialde"`, `"sparkx"`
- `params.n_top_genes`: Number of top SVGs to return

**Example Queries**:
```
"Identify spatial patterns with SpatialDE"
"Use SPARK-X to find spatial genes"
"Find spatially variable genes"
```

**Returns**: Ranked spatially variable genes

---


## 🟣 Cell Annotation Tools

### `annotate_cell_types`
**Purpose**: Identify cell types using multiple methods
**Difficulty**: 🟡 Intermediate

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"tangram"`, `"scanvi"`, `"cellassign"`, `"mllmcelltype"`, `"sctype"`, `"singler"`
- `params.reference_data_id`: Reference dataset ID (for transfer methods)
- `params.cell_type_key`: **REQUIRED** for tangram/scanvi/singler - cell type column in reference data
- `params.marker_genes`: Custom marker gene dict for cellassign method

**Method Requirements**:
- **tangram/scanvi/singler**: Require `reference_data_id` + `cell_type_key`
- **cellassign**: Requires `marker_genes` parameter
- **sctype/mllmcelltype**: No reference needed (use built-in databases/LLM)

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
- `params.method`: `"scvelo"` (standard, tested ✅), `"velovi"` (deep learning, tested ✅)
- `params.mode`: `"dynamical"`, `"stochastic"`, `"deterministic"` (for scVelo)
- `params.n_top_genes`: Velocity genes to use
- `params.velovi_n_epochs`: Training epochs (for VeloVI)

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
- `params.cell_type_column`: **REQUIRED** - cell type/cluster column in spatial data
- `params.species`: `"human"`, `"mouse"`, `"zebrafish"` (default: "human")
- `params.liana_resource`: LR database - use `"mouseconsensus"` for mouse
- `params.data_source`: `"current"` or `"raw"` (use "raw" for full gene coverage)
- `params.perform_spatial_analysis`: Spatial vs cluster-based (default: True)

**Common Issues**:
- **"Too few features"**: Use `data_source="raw"` + correct species/resource
- **Missing connectivity**: Set `spatial_connectivity_handling="compute_with_params"`
- **Species mismatch**: Mouse data must use `species="mouse"` + `liana_resource="mouseconsensus"`

**Example Queries**:
```
"Analyze ligand-receptor interactions with LIANA"
"Find spatial cell communication patterns"
"Identify tumor-immune interactions"
```

**Returns**: Significant ligand-receptor pairs and spatial patterns

**Note**: Requires cell type annotations. Run `annotate_cell_types` or ensure `cell_type_column` exists in your data before cell communication analysis.

---

### `analyze_enrichment`
**Purpose**: Gene set enrichment analysis
**Difficulty**: 🟡 Intermediate

**Key Parameters**:
- `data_id`: Dataset ID
- `params.species`: **REQUIRED** - `"human"`, `"mouse"`, or `"zebrafish"`
- `params.method`: `"spatial_enrichmap"`, `"pathway_gsea"`, `"pathway_ora"`, `"pathway_enrichr"`
- `params.gene_set_database`: `"GO_Biological_Process"`, `"KEGG_Pathways"`, `"Reactome_Pathways"`
- `params.pvalue_cutoff`: Significance threshold (default: 0.05)

**Important**:
- **Species must match your data** (human=uppercase genes, mouse=capitalized)
- Mouse data: Use `species="mouse"` for correct gene name matching
- KEGG database is species-specific (KEGG_2021_Human, KEGG_2019_Mouse)

**Note**: Best used after `find_markers` (differential expression) or `find_spatial_genes` (spatial patterns) to enrich biologically meaningful gene sets

**Example Queries**:
```
"Perform spatial GSEA analysis"
"Find enriched pathways in tumor regions"
"Analyze GO term enrichment"
```

**Returns**: Enriched gene sets with statistics

---

### `analyze_cnv`
**Purpose**: Detect copy number variations (CNVs) from spatial transcriptomics
**Difficulty**: 🔴 Advanced

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"infercnvpy"`, `"numbat"`
- `params.reference_key`: Cell type/cluster column for reference cells
- `params.reference_categories`: List of normal cell types (e.g., `["T cells", "B cells"]`)
- `params.window_size`: CNV averaging window size (default: 100 genes)
- `params.step`: Sliding window step size (default: 10)
- `params.exclude_chromosomes`: Chromosomes to exclude (e.g., `["chrX", "chrY"]`)
- `params.cluster_cells`: Cluster cells by CNV pattern (default: False)

**Method Comparison**:
- **infercnvpy**: Expression-based, fast, GPU-friendly ⭐ (default)
- **numbat**: Haplotype-aware, requires allele data, more accurate (R-based)

**Example Queries**:
```
"Detect CNVs using immune cells as reference"
"Analyze copy number variations with inferCNVpy"
"Find chromosomal alterations in tumor cells"
"Use Numbat for haplotype-aware CNV detection"
```

**Returns**: CNV scores per cell/gene with clustering and visualization

**Note**: Requires reference (normal) cells for baseline. Numbat requires allele count data in adata layers.

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

**Key Parameters**:
- `source_id`: Source dataset ID
- `target_id`: Target dataset ID
- `method`: `"paste"`, `"stalign"`
- `landmarks`: Optional alignment landmarks

**Example Queries**:
```
"Align consecutive tissue sections"
"Register brain slices for 3D reconstruction"
"Align samples to reference coordinates"
```

**Returns**: Transformation matrix and aligned coordinates in `obsm['spatial_registered']`

**Note**:
- **PASTE**: Optimal transport-based alignment, supports multi-slice registration
- **STalign**: Diffeomorphic LDDMM registration, pairwise only

---

## 🟤 Deconvolution & Domain Tools

### `deconvolve_data`
**Purpose**: Deconvolve spots into cell type proportions
**Difficulty**: 🔴 Advanced

**Key Parameters**:
- `data_id`: Dataset ID
- `params.method`: `"cell2location"`, `"destvi"`, `"stereoscope"`, `"tangram"`, `"rctd"`, `"spotlight"`, `"card"`
- `params.reference_data_id`: Single-cell reference ID
- `params.cell_type_key`: **REQUIRED** - cell type column in reference data
- `params.n_epochs`: Training iterations (default: 30000 for cell2location)
- `params.n_cells_per_spot`: Expected cells per spot (default: 30)

**Important Notes**:
- All methods require reference single-cell data with cell type annotations
- `cell_type_key` must match a column in reference `adata.obs` (e.g., 'cell_type', 'annotation')
- Cell2location: GPU recommended, uses 2-stage training (250 + 30000 epochs)

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
- `params.method`: `"spagcn"`, `"leiden"`, `"louvain"`, `"stagate"`, `"graphst"`
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
- `params.plot_type`: **20 plot types available** (see complete list below)
- `params.feature`: Gene(s) or metadata column to visualize (str or List[str])
- `params.colormap`: Color scheme (default: `"viridis"`)
- `params.dpi`: Image resolution (default: 300 for publication quality)
- `params.analysis_type`: **REQUIRED** for `plot_type="spatial_statistics"` - see spatial analysis types
- `params.cluster_key`: **REQUIRED** for `plot_type="heatmap"` - grouping column in adata.obs

#### Complete List of 20 Plot Types:

**Basic Spatial & Dimensionality Reduction** (4 types):
- `"spatial"` ⭐ - Basic spatial gene/metadata expression (default)
- `"umap"` - UMAP dimensionality reduction embedding
- `"violin"` - Violin plots by group
- `"heatmap"` - Expression heatmap (requires `cluster_key`)

**Deconvolution** (1 type):
- `"deconvolution"` - Cell type proportion visualization

> **Note**: For spatial domains visualization, use `plot_type="spatial"` with `feature` set to the domain_key returned by `identify_spatial_domains` (e.g., `"spatial_domains_spagcn"`, `"spatial_domains_leiden"`)

**Trajectory & Velocity** (2 types):
- `"trajectory"` - Trajectory/pseudotime plots (CellRank/Palantir)
- `"rna_velocity"` - RNA velocity stream plots (scVelo/VeloVI)

**Cell Communication** (3 types):
- `"cell_communication"` - Ligand-receptor interaction networks
- `"lr_pairs"` - Specific L-R pair spatial visualization
- `"spatial_interaction"` - Spatial proximity-based interactions

**Gene Analysis** (2 types):
- `"multi_gene"` - Multi-gene panel visualization
- `"gene_correlation"` - Gene-gene correlation plots

**Enrichment & Pathways** (2 types):
- `"pathway_enrichment"` - Pathway enrichment results
- `"spatial_enrichment"` - Spatially enriched pathways

**Batch & Integration** (1 type):
- `"batch_integration"` - Integration quality assessment (UMAP/metrics)

**CNV Analysis** (2 types):
- `"cnv_heatmap"` - Copy number variation heatmap
- `"spatial_cnv"` - CNV spatial projection

**Advanced Deconvolution** (1 type):
- `"card_imputation"` - CARD high-resolution imputation results

**Spatial Statistics** (1 type):
- `"spatial_statistics"` - Spatial statistics visualization (requires `analysis_type`)
  - Available `analysis_type` values: `"neighborhood"`, `"co_occurrence"`, `"ripley"`, `"moran"`, `"centrality"`, `"getis_ord"`

#### Parameter Requirements by Plot Type:

| Plot Type | REQUIRED Parameters | Optional Key Parameters |
|-----------|---------------------|------------------------|
| `heatmap` | `cluster_key` | `feature`, `colormap` |
| `spatial_analysis` | `analysis_type` | `cluster_key` |
| `lr_pairs` | - | `lr_pairs` (list of tuples), `feature` (L^R format) |
| All others | - | `feature`, `colormap` |

**Example Queries**:
```
"Show spatial expression of CD3D"
"Create UMAP plot colored by cell types"
"Plot heatmap of top marker genes grouped by leiden clusters"
"Visualize neighborhood enrichment results"  # spatial_analysis with analysis_type="neighborhood"
"Show ligand-receptor pair Fn1^Cd79a spatially"
"Generate trajectory visualization with pseudotime"
"Display CNV heatmap for tumor cells"
```

**Returns**: High-resolution image (300 DPI) for display or publication

---

## 🚀 Quick Start Combinations

### Basic Analysis Workflow
```
1. load_data → 2. preprocess_data → 3. identify_spatial_domains → 4. visualize_data
```

### Cell Type Analysis
```
1. load_data → 2. preprocess_data → 3. annotate_cell_types → 4. find_markers → 5. visualize_data
```

### Communication Analysis  
```
1. load_data → 2. preprocess_data → 3. annotate_cell_types → 4. analyze_cell_communication → 5. visualize_data
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
- **GPU acceleration available for:**
  - **Annotation**: Tangram (tangram_device='cuda:0')
  - **Spatial Domains**: STAGATE, GraphST
  - **Deconvolution**: Cell2location ⭐, CARD
  - **Integration**: All scVI-based methods (scVI, scANVI)
  - **Velocity**: VeloVI, scVelo (partial)
  - **CNV**: inferCNVpy (GPU-friendly)
- Sparse matrices automatically handled by the tools

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
- `graphst`: Graph-based spatial transcriptomics

**Cell Annotation Methods:**
- `tangram`: Best with scRNA reference ⭐
- `scanvi`: Deep learning label transfer
- `cellassign`: Probabilistic, custom markers
- `sctype`: Automatic, good for exploration (R-based)
- `singler`: Reference-based annotation (R-based)
- `mllmcelltype`: LLM-based multimodal annotation

**Deconvolution Methods:**
- `cell2location`: GPU recommended, high accuracy ⭐
- `rctd`: CPU-friendly, robust (R-based)
- `destvi`: Deep learning, experimental
- `stereoscope`: Fast, good baseline
- `spotlight`: SPOTlight (R-based)
- `tangram`: Spatial mapping mode
- `card`: CARD deconvolution (R-based)

### Performance Guidelines

| Dataset Size | Recommended Tools | Memory | Time | Notes |
|-------------|------------------|---------|-------|-------|
| Small (<5K spots) | Any method | 4GB+ | 5-15 min | Explore all options |
| Medium (5-20K spots) | SpaGCN, SPARK-X | 8GB+ | 15-45 min | Standard workflows |
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