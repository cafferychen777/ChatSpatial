# Methods Reference

This page is the canonical human-readable reference for ChatSpatial tool names,
method names, common defaults, accepted values, and user-facing parameter
behavior. MCP clients expose the full schema for method-specific advanced
options.

ChatSpatial's public interface is a set of 21 schema-validated MCP tools. Those
tools orchestrate 67 spatial transcriptomics methods across 15 analytical
categories. In this page, **tool** means the MCP entry point you or an AI client
can call; **method** means an algorithm or analysis backend selected through a
parameter such as `method`, `analysis_type`, `plot_type`, or `subtype`.

---

## Quick Reference

| Category | Tools |
|----------|-------|
| Data | `load_data`, `predict_spatial_expression_from_histology`, `preprocess_data`, `compute_embeddings`, `export_data`, `reload_data` |
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

### predict_spatial_expression_from_histology

Generate virtual spatial transcriptomics from pre-cut H&E tiles with
[DeepSpot-M](https://github.com/ratschlab/DeepSpotM). This is a second data
entry point, for archival or clinical material that carries no spatial assay.
It registers a new dataset and returns its `data_id`; `export_data` remains the
way to write it to disk.

| Parameter | Type | Description |
|-----------|------|-------------|
| `manifest_path` | str | Path to the v1 CSV coordinate manifest |
| `tile_directory` | str | Directory that every manifest `tile_path` resolves against |

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gene_embedding_source` | `scgpt` | `evo2`, `orthrus`, `prott5`, `scgpt`, `apertus` |
| `genes` | None | HGNC symbols to predict. None predicts the model's full panel |
| `model_repository` | `ratschlab/DeepSpotM` | Hugging Face repo id, or a local checkpoint directory |
| `model_revision` | None | Checkpoint branch, tag, or commit. None resolves to `main` |
| `batch_size` | 16 | Tiles per forward pass |
| `use_gpu` | True | Use CUDA when available, otherwise CPU |
| `name` | None | Display name for the registered dataset |
| `timeout` | 600 | Seconds allowed for model loading and inference |

#### Manifest v1

A CSV with these required columns:

| Column | Description |
|--------|-------------|
| `tile_id` | Unique identifier. Becomes the observation name |
| `tile_path` | Path relative to `tile_directory`. Must not escape it |
| `slide_id` | Manifest-wide slide identifier, repeated on every row |
| `x_px` | Tile upper-left x in the native level-0 slide frame |
| `y_px` | Tile upper-left y in the native level-0 slide frame |
| `mpp_x` | Manifest-wide microns per pixel along x, repeated on every row |
| `mpp_y` | Manifest-wide microns per pixel along y, repeated on every row |

The coordinate frame has a top-left origin, x increasing rightward and y
increasing downward. `slide_id`, `mpp_x` and `mpp_y` repeat on every row so the
manifest stays a portable single file, and the values must be identical across
rows. The fixed v1 schema defines the convention, so no companion metadata file
is needed. One slide and one tile directory per invocation.

`grid_row` and `grid_col` are optional generic lattice indices. They must appear
together, and they are written to `adata.obs` only when one affine mapping
relates them to the pixel coordinates (`x_px = stride * grid_col + offset`, and
likewise for rows). Manifests whose tiles do not come from a regular
axis-aligned lattice simply omit them and stay pixel-only. These are not Visium
`array_row`/`array_col` and are never presented as such.

```csv
tile_id,tile_path,slide_id,x_px,y_px,mpp_x,mpp_y,grid_row,grid_col
TCGA-XX_0_0,tiles/0_0.png,TCGA-XX,0,0,0.5,0.5,0,0
TCGA-XX_224_0,tiles/224_0.png,TCGA-XX,224,0,0.5,0.5,0,1
TCGA-XX_0_224,tiles/0_224.png,TCGA-XX,0,224,0.5,0.5,1,0
```

Every tile must decode as an RGB 224x224 image cut at roughly 20x (about 0.5
microns per pixel). ChatSpatial neither resamples tiles nor infers
magnification, so tiles outside that geometry are rejected rather than
converted. Missing files, duplicate tile identifiers, duplicate coordinates,
non-finite or negative coordinates, mixed slide identifiers, an inconsistent
lattice, and an empty manifest are all rejected with the offending rows named.
Manifest row order is preserved.

#### Output

`adata.X` holds predicted log1p-CPM once; it is not duplicated in a layer.
`adata.obsm["spatial"]` holds tile centers, `(x_px + width / 2, y_px + height /
2)`, which is the convention scanpy and squidpy expect.

Two metadata blocks describe the result. The generic one is backend
independent:

```python
adata.uns["chatspatial"]["expression"]
# {"schema_version": 1, "provenance": "predicted",
#  "units": "log1p_cpm", "producer": "deepspotm"}
```

`adata.uns["deepspotm"]` carries the backend detail: resolved model repository
and checkpoint revision, gene-embedding source, tile geometry and coordinate
convention, microns per pixel, portable manifest filename and SHA-256 digest,
installed package version, license identifiers, and the research-use-only notice.

#### Which tools accept these datasets

Tools branch on the generic provenance field, not on the producer, through the
shared helpers in `chatspatial.utils.provenance`.

**Accepted.** `compute_embeddings`, `visualize_data`,
`analyze_spatial_statistics`, `identify_spatial_domains`,
`analyze_trajectory_data`, `export_data` and `reload_data`.

**Rejected by `preprocess_data`.** The matrix is already log1p-CPM, so quality
control, count-based normalization and highly variable gene selection do not
apply. Run `compute_embeddings` directly.

**Rejected through the shared raw-count access paths.**
`get_raw_data_source(..., require_integer_counts=True)` and
`ensure_counts_layer` refuse before value-based count detection runs. This
covers count-requiring branches such as PyDESeq2, condition comparison,
SpatialDE, SPARK-X, count-based deconvolution, scANVI, Numbat and scVI. The
guard follows the count requirement, not the tool name: branches that
explicitly accept normalized expression can use the predicted matrix and the
helper reports it as non-count data even if a small prediction happens to look
integer-valued.

Normalized-expression consumers such as Scanpy marker testing, FlashS,
enrichment background/ranking preparation and compatible communication or
annotation branches therefore remain available. Their biological
interpretation is still model-predicted expression rather than measured assay
evidence, as recorded in the generic provenance block.

**Audit of paths that bypass the helpers.** `analyze_velocity_data` needs
`spliced` and `unspliced` layers, which a tile-derived dataset never has, so it
already fails with an accurate message. Four paths reach expression directly:
the CellPhoneDB export, InferCNVPy, Tangram annotation, and the PASTE branch of
`register_spatial_data`. CellPhoneDB still resolves through
`get_raw_data_source` during validation and may proceed when its selected
method accepts normalized expression. InferCNVPy and Tangram consume normalized
expression. PASTE is rejected explicitly because its current preparation would
normalize the predicted log1p-CPM matrix a second time; STalign remains the
registration path that does not renormalize expression.

#### Installation and licensing

```bash
pip install 'chatspatial[deepspotm]'
```

The extra is deliberately separate and is excluded from `[full]`, because
DeepSpot-M is non-commercial: the code is
[PolyForm Noncommercial 1.0.0](https://github.com/ratschlab/DeepSpotM) and the
published weights are CC-BY-NC-SA-4.0. ChatSpatial itself stays MIT; the
package is imported only when the tool runs.

The weights are gated. Request access at
[ratschlab/DeepSpotM](https://huggingface.co/ratschlab/DeepSpotM) and run
`huggingface-cli login` before the first prediction. Output is for research use
only, not for clinical or diagnostic use.

Method reference: [DeepSpot-M: a multimodal foundation model for
transcriptome-wide virtual spatial transcriptomics from
histology](https://doi.org/10.64898/2026.06.19.26356060), *medRxiv* (2026).

---

### preprocess_data

Normalize, filter, and prepare data.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `normalization` | `pearson_residuals` | `log`, `sct`, `pearson_residuals`, `scvi`, `none` |
| `n_hvgs` | 2000 | Highly variable genes |
| `filter_genes_min_cells` | 3 | Min cells per gene |
| `filter_cells_min_genes` | 30 | Min genes per cell |
| `filter_mito_pct` | 20.0 | Max mitochondrial % |
| `scale` | False | Scale to unit variance before PCA |

**Advanced options**:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_scrublet` | False | Enable doublet detection (for single-cell resolution data) |
| `normalize_target_sum` | None | Target counts per cell (None=median, 1e4=Visium, 1e6=MERFISH) |
| `remove_mito_genes` | True | Exclude mito genes from HVG |
| `batch_key` | `batch` | Batch column for batch-aware normalization |

PCA, neighbor, UMAP, and clustering parameters belong to
`compute_embeddings`; `preprocess_data` does not accept or retain them.

`preprocess_data` does not compute PCA, UMAP, clustering, or neighbor graphs.
Call `compute_embeddings` explicitly before cluster- or embedding-dependent tools.

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
| `path` | auto | Custom path (default: `~/.chatspatial/active/{data_id}.h5ad`) |

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
| `local_join_count` | Group | Yes |
| `centrality` | Network | Yes |
| `network_properties` | Network | Yes |
| `spatial_centrality` | Network | Yes |

---

### find_spatial_genes

Identify spatially variable genes.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `flashs` | `sparkx`, `flashs`, `spatialde` |
| `n_top_genes` | None | Top genes to return (None = all significant) |

---

### identify_spatial_domains

Find tissue domains and spatial niches.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `spagcn` | `spagcn`, `stagate`, `graphst`, `banksy`, `aestetik`, `leiden`, `louvain` |
| `n_domains` | 7 | Expected number of domains |
| `resolution` | 0.5 | Clustering resolution |

**AESTETIK** fuses a precomputed expression embedding, a precomputed per-spot morphology embedding, and the spatial neighborhood grid. It reads both representations from `adata.obsm` and does not extract morphology features from tissue images. It also needs discrete lattice coordinates in `adata.obs`, either `x_array`/`y_array` or the Visium `array_row`/`array_col` columns.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `aestetik_transcriptomics_key` | `X_pca` | obsm key with the expression embedding |
| `aestetik_morphology_key` | `X_pca_morphology` | obsm key with the morphology embedding |
| `aestetik_morphology_weight` | 1.5 | Morphology weight in the joint loss |
| `aestetik_window_size` | 3 | Odd side length of the neighborhood grid |
| `aestetik_clustering_method` | `kmeans` | `bgm`, `kmeans`, `louvain`, `leiden`. `kmeans`/`bgm` use `n_domains`; `leiden`/`louvain` use `resolution` |
| `aestetik_latent_dim` | 16 | Latent embedding dimension |
| `aestetik_max_epochs` | 100 | Training epochs |
| `aestetik_random_seed` | 2023 | Random seed |

Install with `pip install 'chatspatial[aestetik]'` (Python < 3.14). Method reference: [Representation learning for multi-modal spatially resolved transcriptomics data](https://doi.org/10.1093/bioinformatics/btag316), *Bioinformatics* (2026).

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
| `sctype_tissue` | None | scType tissue name, required unless custom markers are provided |
| `sctype_db_` | None | Local scType database path, or remote URL when explicitly allowed |
| `sctype_custom_markers` | None | Custom scType marker sets |
| `sctype_scaled` | True | Whether scType should treat the expression matrix as scaled |
| `sctype_allow_remote` | False | One-off opt-in to load scType remote R scripts and default marker database |
| `sctype_allow_runtime_r_install` | False | One-off opt-in to install missing R packages at runtime |

**scType remote resources**: by default, scType does not load remote R scripts or the remote default marker database. For one-off exploratory runs, pass `sctype_allow_remote=true`. For production or offline workflows, prefer local R scripts via `CHATSPATIAL_SCTYPE_R_DIR` and a local `sctype_db_` path.

```json
{
  "method": "sctype",
  "sctype_tissue": "Immune system",
  "sctype_allow_remote": true
}
```

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
| `method` | `spatial_enrichmap` | `spatial_enrichmap`, `pathway_gsea`, `pathway_ora`, `pathway_enrichr`, `pathway_ssgsea` |
| `gene_set_database` | `GO_Biological_Process` | See databases below |

**Databases**: `GO_Biological_Process`, `GO_Molecular_Function`, `GO_Cellular_Component`, `KEGG_Pathways`, `Reactome_Pathways`, `MSigDB_Hallmark`, `Cell_Type_Markers`

---

## Dynamics

### analyze_velocity_data

RNA velocity analysis.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `scvelo` | `scvelo`, `velovi` |
| `scvelo_mode` | `stochastic` | `deterministic`, `stochastic`, `dynamical` |

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
| `data_ids` | required | Two or more distinct dataset IDs |
| `method` | `harmony` | `harmony`, `bbknn`, `scanorama`, `scvi` |
| `batch_key` | `batch` | Batch column |

---

### register_spatial_data

Align spatial sections.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `source_id` | required | Source dataset; must differ from `target_id` |
| `target_id` | required | Target dataset; must differ from `source_id` |
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
| `output_format` | `png` | `png`, `pdf`, `svg`, `eps`, `tiff`, `jpg` |

**Plot types and subtypes**:

| Type | Subtypes | Use |
|------|----------|-----|
| `feature` | — | Gene/metadata on spatial or UMAP |
| `expression` | `heatmap`, `violin`, `dotplot`, `correlation` | Aggregated expression |
| `deconvolution` | `spatial_multi`, `pie`, `dominant`, `diversity`, `umap`, `imputation` | Cell proportions |
| `communication` | `dotplot`, `tileplot`, `circle_plot` | LR interactions |
| `interaction` | — | Spatial LR pairs |
| `trajectory` | `pseudotime`, `circular`, `fate_map`, `gene_trends`, `fate_heatmap`, `palantir` | Pseudotime |
| `velocity` | `stream`, `phase`, `proportions`, `heatmap`, `paga` | RNA velocity |
| `statistics` | `neighborhood`, `co_occurrence`, `ripley`, `moran`, `centrality`, `getis_ord` | Spatial stats |
| `enrichment` | `barplot`, `dotplot` | Pathway results |
| `cnv` | `heatmap`, `spatial` | CNV results |
| `integration` | `batch`, `cluster`, `highlight` | Integration QC |

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
