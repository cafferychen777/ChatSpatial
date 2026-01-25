# Methods Reference

All 20 ChatSpatial tools with parameters and options.

---

## Quick Reference

| Category | Tools |
|----------|-------|
| Data | `load_data`, `preprocess_data`, `compute_embeddings`, `export_data`, `reload_data` |
| Spatial | `analyze_spatial_statistics`, `find_spatial_genes`, `identify_spatial_domains` |
| Cells | `annotate_cell_types`, `deconvolve_data`, `analyze_cell_communication` |
| Genes | `find_markers`, `compare_conditions`, `analyze_enrichment` |
| Dynamics | `analyze_velocity_data`, `analyze_trajectory_data`, `analyze_cnv` |
| Multi-sample | `integrate_samples`, `register_spatial_data` |
| Output | `visualize_data` |

---

## Data Management

### load_data

Load spatial transcriptomics data.

| Parameter | Type | Description |
|-----------|------|-------------|
| `data_path` | str | Path to file or folder |
| `data_type` | str | `visium`, `xenium`, `slide_seq`, `merfish`, `seqfish`, `generic` |
| `name` | str | Optional dataset name |

**Supported formats**: H5AD, 10X Visium folders, H5, MTX

---

### preprocess_data

Normalize, filter, and prepare data.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `normalization` | `log` | `log`, `sct`, `pearson_residuals`, `scvi`, `none` |
| `n_hvgs` | 2000 | Highly variable genes |
| `n_pcs` | 30 | Principal components |
| `n_neighbors` | 15 | Neighbor graph |
| `clustering_resolution` | 1.0 | Leiden clustering |
| `filter_genes_min_cells` | 3 | Min cells per gene |
| `filter_cells_min_genes` | 30 | Min genes per cell |
| `filter_mito_pct` | 20.0 | Max mitochondrial % |
| `scale` | False | Scale to unit variance before PCA |

**Advanced options**:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `scrublet_enable` | False | Enable doublet detection (for single-cell resolution data) |
| `normalize_target_sum` | None | Target counts per cell (None=median, 1e4=Visium, 1e6=MERFISH) |
| `remove_mito_genes` | True | Exclude mito genes from HVG |
| `batch_key` | `batch` | Batch column for batch-aware normalization |

---

### compute_embeddings

Compute dimensionality reduction and clustering.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `compute_pca` | True | Compute PCA |
| `compute_umap` | True | Compute UMAP |
| `compute_clustering` | True | Leiden clustering |
| `compute_spatial_neighbors` | True | Spatial graph |
| `n_pcs` | 30 | Principal components |
| `clustering_resolution` | 1.0 | Clustering resolution |
| `force` | False | Recompute if exists |

---

### export_data / reload_data

Export dataset for external scripts, reload after modifications.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `data_id` | required | Dataset ID |
| `path` | auto | Custom path (default: `~/.chatspatial/active/`) |

---

## Spatial Analysis

### analyze_spatial_statistics

Analyze spatial patterns and autocorrelation.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `analysis_type` | `neighborhood` | See types below |
| `cluster_key` | None | Required for group-based analyses |
| `genes` | None | Specific genes to analyze |
| `n_top_genes` | 20 | Top HVGs to analyze (if genes not specified) |
| `n_neighbors` | 8 | Spatial neighbors |

**Analysis types**:

| Type | Category | Requires cluster_key |
|------|----------|---------------------|
| `moran` | Gene | No |
| `local_moran` | Gene | No |
| `geary` | Gene | No |
| `getis_ord` | Gene | No |
| `bivariate_moran` | Gene | No |
| `neighborhood` | Group | Yes |
| `co_occurrence` | Group | Yes |
| `ripley` | Group | Yes |
| `join_count` | Group | Yes |
| `centrality` | Network | Optional |

---

### find_spatial_genes

Identify spatially variable genes.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `sparkx` | `sparkx`, `spatialde` |
| `n_top_genes` | None | Top genes to return (None = all significant) |

---

### identify_spatial_domains

Find tissue domains and spatial niches.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `spagcn` | `spagcn`, `stagate`, `graphst`, `leiden`, `louvain` |
| `n_domains` | 7 | Expected number of domains |
| `resolution` | 0.5 | Clustering resolution |

---

## Cell Analysis

### annotate_cell_types

Assign cell types.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `tangram` | See methods below |
| `reference_data_id` | None | Reference dataset (for transfer methods) |
| `cell_type_key` | None | Cell type column in reference |
| `marker_genes` | None | Marker dict (for CellAssign) |

**Methods**:

| Method | Requires Reference | Notes |
|--------|-------------------|-------|
| `tangram` | Yes | Spatial mapping |
| `scanvi` | Yes | Deep learning transfer |
| `cellassign` | No | Marker-based |
| `sctype` | No | Automatic (R) |
| `singler` | No | Reference-based (R) |
| `mllmcelltype` | No | LLM-based |

---

### deconvolve_data

Estimate cell type proportions per spot.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `flashdeconv` | See methods below |
| `reference_data_id` | required | Reference dataset |
| `cell_type_key` | required | Cell type column in reference |

**Methods**:

| Method | Speed | GPU | Notes |
|--------|-------|-----|-------|
| `flashdeconv` | Fast | No | Default, recommended |
| `cell2location` | Slow | Yes | High accuracy |
| `rctd` | Fast | No | R-based |
| `destvi` | Medium | Yes | scvi-tools |
| `stereoscope` | Slow | Yes | Alternative DL |
| `tangram` | Medium | Yes | Spatial mapping |
| `spotlight` | Fast | No | R-based |
| `card` | Fast | No | R-based, imputation |

---

### analyze_cell_communication

Analyze ligand-receptor interactions.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `fastccc` | `fastccc`, `liana`, `cellphonedb`, `cellchat_r` |
| `species` | required | `human`, `mouse`, `zebrafish` |
| `cell_type_key` | required | Cell type column |
| `liana_resource` | `consensus` | LR database (`mouseconsensus` for mouse) |

---

## Gene Analysis

### find_markers

Find differentially expressed genes.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `group_key` | required | Grouping column |
| `group1` | None | First group (None = each vs rest) |
| `group2` | None | Second group |
| `method` | `wilcoxon` | `wilcoxon`, `t-test`, `t-test_overestim_var`, `logreg`, `pydeseq2` |
| `n_top_genes` | 50 | Top genes per group |

---

### compare_conditions

Compare experimental conditions (pseudobulk DESeq2).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `condition_key` | required | Condition column |
| `condition1` | required | Treatment group |
| `condition2` | required | Control group |
| `sample_key` | required | Sample/patient column |
| `cell_type_key` | None | Stratify by cell type |
| `n_top_genes` | 50 | Top DEGs |

---

### analyze_enrichment

Gene set enrichment analysis.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `species` | required | `human`, `mouse`, `zebrafish` |
| `method` | `pathway_ora` | `pathway_ora`, `pathway_gsea`, `pathway_ssgsea`, `spatial_enrichmap` |
| `gene_set_database` | `GO_Biological_Process` | See databases below |

**Databases**: `GO_Biological_Process`, `GO_Molecular_Function`, `KEGG_Pathways`, `Reactome_Pathways`, `MSigDB_Hallmark`

---

## Dynamics

### analyze_velocity_data

RNA velocity analysis.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `scvelo` | `scvelo`, `velovi` |
| `mode` | `stochastic` | `deterministic`, `stochastic`, `dynamical` |

**Requires**: `spliced` and `unspliced` layers

---

### analyze_trajectory_data

Trajectory and pseudotime inference.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `cellrank` | `cellrank`, `palantir`, `dpt` |
| `root_cells` | None | Starting cells |

**Note**: CellRank requires velocity data

---

### analyze_cnv

Copy number variation detection.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `infercnvpy` | `infercnvpy`, `numbat` |
| `reference_key` | required | Cell type column |
| `reference_categories` | required | Normal cell types |

---

## Multi-Sample

### integrate_samples

Batch integration.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `data_ids` | required | List of dataset IDs |
| `method` | `harmony` | `harmony`, `bbknn`, `scanorama`, `scvi` |
| `batch_key` | `batch` | Batch column |

---

### register_spatial_data

Align spatial sections.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `source_id` | required | Source dataset |
| `target_id` | required | Target dataset |
| `method` | `paste` | `paste`, `stalign` |

---

## Visualization

### visualize_data

Create all plot types.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `plot_type` | `feature` | See types below |
| `subtype` | None | Visualization variant |
| `feature` | None | Gene(s) or column to show |
| `basis` | `spatial` | `spatial`, `umap` |
| `cluster_key` | None | Grouping column |
| `colormap` | `coolwarm` | Color scheme |
| `dpi` | 300 | Resolution |
| `output_format` | `png` | `png`, `pdf`, `svg` |

**Plot types and subtypes**:

| Type | Subtypes | Use |
|------|----------|-----|
| `feature` | — | Gene/metadata on spatial or UMAP |
| `expression` | `heatmap`, `violin`, `dotplot`, `correlation` | Aggregated expression |
| `deconvolution` | `spatial_multi`, `pie`, `dominant`, `diversity`, `umap` | Cell proportions |
| `communication` | `dotplot`, `tileplot`, `circle_plot` | LR interactions |
| `interaction` | — | Spatial LR pairs |
| `trajectory` | `pseudotime`, `fate_map`, `gene_trends` | Pseudotime |
| `velocity` | `stream`, `phase`, `paga` | RNA velocity |
| `statistics` | `neighborhood`, `co_occurrence`, `ripley`, `moran` | Spatial stats |
| `enrichment` | `barplot`, `dotplot` | Pathway results |
| `cnv` | `heatmap`, `spatial` | CNV results |
| `integration` | `batch`, `cluster` | Integration QC |

---

## GPU Acceleration

Set `use_gpu=True` for these methods:

| Category | Methods |
|----------|---------|
| Preprocessing | scVI normalization |
| Annotation | Tangram, scANVI |
| Deconvolution | Cell2location, DestVI, Stereoscope, Tangram |
| Domains | STAGATE, GraphST |
| Velocity | VeloVI |
| Integration | scVI |
| CNV | inferCNVpy |

---

## Next Steps

- [Examples](../examples.md) — See methods in action
- [Concepts](../concepts.md) — When to use which method
- [Troubleshooting](troubleshooting.md) — Common issues and fixes
