# ChatSpatial Reproducibility

Scripts, result data, and supplementary tables for reproducing the experiments in:

> **ChatSpatial: Schema-Enforced Orchestration of Cross-Ecosystem Spatial Transcriptomics Workflows**
>
> Chen Yang, Xianyang Zhang, Jun Chen

This repository accompanies the main [ChatSpatial package](https://github.com/cafferychen777/ChatSpatial) and contains all analysis scripts used to generate the figures and results reported in the manuscript.

---

## Repository Structure

```
ChatSpatial-Reproducibility/
├── README.md
├── .gitignore
├── scripts/                           # All experiment and analysis scripts
│   ├── reproducibility_analysis.py    # Schema constraint coverage analysis
│   ├── determinism_experiment.py      # Single-model determinism experiment
│   ├── determinism_multimodel.py      # Cross-model determinism experiment
│   ├── codegen_hallucination.py       # Code generation baseline comparison
│   ├── ablation_invocation.py         # Schema-enforcement ablation: invocation level
│   ├── ablation_e2e.py               # Schema-enforcement ablation: end-to-end execution
│   ├── ablation_analysis.py          # Ablation output concordance analysis
│   ├── prompt_sensitivity.py         # Prompt sensitivity experiment
│   ├── prompt_sensitivity_analysis.py # Prompt sensitivity analysis
│   ├── cross_system_comparison.py    # Cross-system invocation comparison
│   ├── cross_system_analysis.py      # Cross-system analysis and metrics
│   ├── e2e_benchmark.py              # End-to-end execution benchmark (3 systems)
│   ├── e2e_benchmark_analysis.py     # End-to-end benchmark analysis
│   ├── e2e_driver_stagent.py         # STAgent subprocess driver
│   ├── e2e_driver_spatialagent.py    # SpatialAgent subprocess driver
│   ├── dlpfc_prepare_data.py         # DLPFC ground-truth data preparation
│   ├── dlpfc_benchmark.py            # DLPFC ground-truth benchmark (3 systems)
│   ├── dlpfc_benchmark_analysis.py   # DLPFC benchmark analysis
│   ├── casestudy_reproducibility.py  # Case-study workflow concordance
│   ├── compute_case_study_stats.py   # Effect sizes and confidence intervals
│   └── build_supplementary_tables.py # Generate Supplementary Tables S1-S4
├── data/                              # Experiment results (CSV/TXT summaries)
│   ├── reproducibility_analysis.csv
│   ├── determinism_experiment.csv
│   ├── determinism_multimodel.csv
│   ├── codegen_hallucination.csv
│   ├── ablation/
│   │   ├── invocation/               # Invocation-level ablation results
│   │   ├── e2e/                      # End-to-end ablation results
│   │   └── analysis/                 # Concordance metrics and statistical tests
│   ├── sensitivity/                   # Prompt sensitivity results
│   ├── cross_system/                  # Cross-system comparison results
│   ├── e2e_benchmark/                 # End-to-end benchmark results
│   ├── dlpfc_benchmark/               # DLPFC ground-truth benchmark results
│   └── casestudy_reproducibility/     # Case-study concordance results
└── supplementary_tables/              # Supplementary Tables S1-S4
    ├── Supplementary_Table_1_Effect_Sizes.csv
    ├── Supplementary_Table_2_Integrated_Methods.csv
    ├── Supplementary_Table_3_Test_Scenarios.csv
    └── Supplementary_Table_4_AI_Agent_Comparison.csv
```

---

## Requirements

### Software

- Python 3.11+
- [ChatSpatial](https://github.com/cafferychen777/ChatSpatial) (`pip install chatspatial` or `pip install chatspatial[full]` for R-based methods)
- Additional Python packages: `requests`, `numpy`, `scipy`, `scikit-learn`, `pandas`, `anndata`, `scanpy`
- For R-based analyses: `rpy2`, R 4.1+, and the R packages `CellChat`, `RCTD`/`spacexr`, `spatialLIBD`

### API Keys

The LLM-based experiments require API keys set as environment variables:

```bash
export GEMINI_API_KEY="..."       # Google Gemini
export ANTHROPIC_API_KEY="..."    # Anthropic Claude
export OPENAI_API_KEY="..."       # OpenAI GPT
```

### Installation

```bash
# Clone this repository
git clone https://github.com/cafferychen777/ChatSpatial-Reproducibility.git
cd ChatSpatial-Reproducibility

# Install ChatSpatial and dependencies
pip install chatspatial requests numpy scipy scikit-learn pandas anndata scanpy

# For full environment including R-based methods
pip install chatspatial[full]
```

---

## Datasets

All datasets used in this study are publicly available. No custom datasets were generated. Below we document each dataset's provenance, access method, and role in the manuscript.

### Case Study 1: Oral Squamous Cell Carcinoma (OSCC)

| Attribute | Value |
|-----------|-------|
| GEO Accession | [GSE208253](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208253) |
| Original Paper | Arora et al., "Spatial transcriptomics reveals distinct and conserved tumor core and edge architectures that predict survival and targeted therapy response," *Nature Communications* 14, 5029 (2023) |
| Technology | 10x Genomics Visium |
| Samples | 12 Visium spatial transcriptomics slides from OSCC patients |
| Tissue | Oral squamous cell carcinoma (primary tumor resections) |
| Spot Count | ~1,000-4,500 spots per sample (after quality filtering at 200 genes/spot) |
| Gene Count | ~15,000 genes per sample |
| Data Format | Space Ranger output (filtered feature-barcode matrices + spatial coordinates) |
| Download | GEO supplementary files or Zenodo mirror (DOI: [10.5281/zenodo.8079095](https://doi.org/10.5281/zenodo.8079095)) |
| Manuscript Role | Case study (Results Section 2.3.1, Figure 2): spatial domain identification (Tumor Core vs. Leading Edge), CARD deconvolution with scRNA-seq reference, CellChat ligand-receptor analysis, Moran's I spatial autocorrelation |

This dataset is used in:
- `casestudy_reproducibility.py` (CARD deconvolution concordance, 160 trials)
- `compute_case_study_stats.py` (effect sizes for enrichment analysis)
- `ablation_e2e.py` (end-to-end ablation on Sample 1)

### HNSCC Single-Cell RNA-seq Reference

| Attribute | Value |
|-----------|-------|
| GEO Accession | [GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322) |
| Original Paper | Puram et al., "Single-cell transcriptomic analysis of primary and metastatic tumor ecosystems in head and neck cancer," *Cell* 171(7), 1611-1624.e24 (2017) |
| Technology | Smart-seq2 scRNA-seq |
| Cells | 5,902 cells from 18 HNSCC patients |
| Cell Types | Tumor cells, fibroblasts, macrophages, T cells, B cells, dendritic cells, mast cells, endothelial cells, myocytes |
| Manuscript Role | Reference atlas for CARD deconvolution in the OSCC case study; cell type annotations used to estimate spot-level cell type proportions |

### Case Study 2: High-Grade Serous Ovarian Carcinoma (HGSOC)

| Attribute | Value |
|-----------|-------|
| GEO Accession | [GSE211956](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211956) |
| Original Paper | Denisenko et al., "Spatial transcriptomics reveals discrete tumour microenvironments and autocrine loops within ovarian cancer subclones," *Nature Communications* 15, 2860 (2024) |
| Technology | 10x Genomics Visium |
| Samples | 8 Visium slides from HGSOC patients (P1-P8), spanning different chemotherapy response categories |
| Companion scRNA-seq | Single-cell reference with 8 consolidated cell types (tumor, macrophage, fibroblast, T cell, B/plasma cell, endothelial, mesothelial, myofibroblast) |
| Additional Platform | NanoString CosMx Spatial Molecular Imager (validation; data on Zenodo DOI: [10.5281/zenodo.10048057](https://doi.org/10.5281/zenodo.10048057)) |
| Manuscript Role | Case study (Results Section 2.3.2, Figure 3): multi-sample RCTD deconvolution, infercnvpy copy-number inference, cross-patient microenvironment comparison, Moran's I spatial autocorrelation |

### Human Lymph Node (End-to-End Benchmark)

| Attribute | Value |
|-----------|-------|
| Source | [10x Genomics Public Datasets](https://www.10xgenomics.com/datasets) |
| Dataset | Human Lymph Node (Visium Spatial Gene Expression) |
| Technology | 10x Genomics Visium |
| Spots | 4,035 spots |
| Genes | 36,601 genes |
| Manuscript Role | End-to-end execution benchmark dataset (Results Section 2.4, Figure 4h-i); chosen because it was not used in any case study, ensuring unbiased cross-system evaluation |

This dataset is used in:
- `e2e_benchmark.py` (3 systems x 3 tasks x 5 replicates = 45 trials)

### DLPFC Visium (Ground-Truth Benchmark)

| Attribute | Value |
|-----------|-------|
| Source | [spatialLIBD](http://spatial.libd.org/spatialLIBD/) R/Bioconductor package |
| Original Paper | Maynard et al., "Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex," *Nature Neuroscience* 24, 425-436 (2021) |
| GEO Accession | [GSE158328](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158328) |
| Technology | 10x Genomics Visium |
| Total Samples | 12 samples from 3 donors (4 samples per donor) |
| Samples Used | 3 samples (151673, 151507, 151669; one per donor) |
| Ground Truth | Expert-annotated cortical layer labels (L1-L6 and white matter, 7 domains) |
| Manuscript Role | Ground-truth benchmark (Results Section 2.4, Figure 4j-k): ARI against expert annotations evaluates biological accuracy, not just consistency |

Ground-truth preparation:
1. `dlpfc_prepare_data.py` downloads samples via `spatialLIBD::fetch_data()`, converts to h5ad via `zellkonverter`, performs standard preprocessing (normalization, HVG selection, PCA, spatial neighbors), and **removes** ground-truth columns from the data provided to systems
2. Ground-truth labels are stored separately in `data/dlpfc_benchmark/ground_truth/` and used only by the analysis script for post-hoc ARI/NMI computation

### Functional Validation Datasets

The following datasets were used for functional validation (Results Section 2.2, Supplementary Table 3) to test ChatSpatial's implementation coverage across diverse spatial transcriptomics platforms:

| Dataset | Accession | Technology | Reference | Spots/Cells | Genes |
|---------|-----------|------------|-----------|-------------|-------|
| SPOTS Benchmark | [GSE198353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198353) | Visium + protein | Ben-Chetrit, Niu et al., *Nat. Biotechnol.* 41, 788-793 (2023) | ~4,000 spots | ~33,000 |
| Visium Multi-Sample (Breast) | [GSE254652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254652) | Visium | See Supplementary Table S2 | Multiple samples | ~36,000 |
| Visium Multi-Sample (Brain) | [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275) | Visium | See Supplementary Table S2 | Multiple samples | ~36,000 |
| MERFISH (Hypothalamus) | [GSE113576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576) | MERFISH | Moffitt et al., *Science* 362, eaau5324 (2018) | ~73,000 cells | 155 genes |
| seqFISH (Mouse Embryo) | [GSE133244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133244) | seqFISH | Lohoff et al., *Nat. Biotechnol.* 40, 74-85 (2022) | ~19,000 cells | 351 genes |
| STARmap (Mouse VISp) | [Wang Lab Data Portal](https://www.wangxiaolab.org/data-portal-1) | STARmap | Wang et al., *Science* 361, eaat5691 (2018) | ~300 spots | 150 genes |
| Slide-seq (Mouse Cerebellum) | [SCP354](https://singlecell.broadinstitute.org/single_cell/study/SCP354) | Slide-seq | Rodriques et al., *Science* 363, 1463-1467 (2019) | ~40,000 beads | ~21,000 |
| Xenium (Various tissues) | [10x Genomics](https://www.10xgenomics.com/datasets) | Xenium | 10x Genomics public dataset | ~150,000 cells | 280+ genes |

---

## Experiment Protocols

### 1. Schema Constraint Coverage Analysis

**Script:** `scripts/reproducibility_analysis.py`
**Paper Section:** Results Section 2.1
**Trial Count:** Static analysis (no API calls)

Performs static analysis of ChatSpatial's 20 MCP tool definitions to quantify the constraint surface: number of typed parameters, Literal enumerations, numeric bounds, defaults, and natural-language descriptions per tool. Validates that the schema design provides sufficient constraint coverage to guide LLM parameter selection.

**Output:** `data/reproducibility_analysis.csv`

### 2. Single-Model Determinism

**Script:** `scripts/determinism_experiment.py`
**Paper Section:** Results Section 2.4 (preliminary experiment)
**Design:** 20 prompts x 20 trials per prompt x 1 model (Gemini 2.0 Flash) at T=1.0
**Trial Count:** 400

Tests whether a single LLM produces deterministic tool invocations under schema enforcement across repeated trials of the same prompt. Measures tool selection consistency and parameter concordance.

**Output:** `data/determinism_experiment.csv`

### 3. Cross-Model Determinism

**Script:** `scripts/determinism_multimodel.py`
**Paper Section:** Results Section 2.4
**Design:** 8 prompts x 3 models x 10 replicates at T=1.0
**Models:** Gemini 2.5 Flash, Claude Haiku 4.5, GPT-5.4
**Trial Count:** 240

Extends the single-model experiment to test whether schema enforcement produces consistent invocations *across* different LLM providers. Measures cross-model tool selection accuracy and parameter agreement.

**Output:** `data/determinism_multimodel.csv`

### 4. Code Generation Baseline

**Script:** `scripts/codegen_hallucination.py`
**Paper Section:** Results Section 2.4
**Design:** Same 8 prompts x 3 models x 10 replicates at T=1.0
**Trial Count:** 240

Baseline comparison where the same LLMs are asked to generate free-form Python code for the same analytical tasks, without schema enforcement. Measures code correctness, import validity, and API usage accuracy to quantify the error rate of unconstrained code generation.

**Output:** `data/codegen_hallucination.csv`

### 5. Schema-Enforcement Ablation

**Scripts:**
- `scripts/ablation_invocation.py` (Part 1: invocation-level)
- `scripts/ablation_e2e.py` (Part 2: end-to-end execution)
- `scripts/ablation_analysis.py` (Part 3: concordance analysis)

**Paper Section:** Results Section 2.4, Figure 4a-c

#### Part 1: Invocation-Level (720 trials)

**Design:** 8 prompts x 3 models x 10 replicates x 3 conditions = 720 trials

Three schema conditions isolate the causal contribution of schema enforcement:

| Condition | Description |
|-----------|-------------|
| **Full schema** | Types + Literal enumerations + bounds + defaults + natural-language descriptions |
| **Bare schema** | Types + enumerations + bounds + defaults, *no* descriptions |
| **No schema** | Tool names + one-sentence summaries only, *no* typed parameters |

All conditions share the same JSON response format, so differences reflect schema content rather than output structure. Each response is validated against ChatSpatial's Pydantic parameter models.

**Metrics:** Parse success rate, tool selection accuracy, Pydantic validation rate

**Output:** `data/ablation/invocation/`

#### Part 2: End-to-End Execution (270 trials)

**Design:** 3 tasks x 3 models x 3 conditions x 10 replicates = 270 trials

Extends the ablation to executable tasks on an OSCC Visium dataset (GSE208253, Sample 1; 1,159 spots x 15,215 genes):

| Task | Method | Output Metric |
|------|--------|---------------|
| Spatial domain identification | Leiden clustering | ARI, NMI |
| Spatially variable gene detection | FlashSV | Jaccard@100 |
| Cell-type deconvolution | FlashDeconv | Pearson r |

Each trial operates on an independent data copy to prevent cross-contamination.

**Output:** `data/ablation/e2e/`

#### Part 3: Concordance Analysis

Computes pairwise cross-model concordance metrics within each condition x task cell, with 95% percentile bootstrap CIs (10,000 resamples, seed=42). Statistical comparisons use Kruskal-Wallis tests with pairwise Mann-Whitney U and Bonferroni correction.

**Output:** `data/ablation/analysis/`

### 6. Prompt Sensitivity

**Scripts:**
- `scripts/prompt_sensitivity.py`
- `scripts/prompt_sensitivity_analysis.py`

**Paper Section:** Results Section 2.4, Figure 4d-e
**Design:** 5 prompt groups x 5 paraphrases x 3 models x 5 replicates = 375 trials

Tests robustness to natural-language variation: five semantically equivalent rephrasings of each analytical task are presented under the full-schema condition. Measures cross-paraphrase constrained-parameter concordance and compares with within-paraphrase consistency from the ablation.

**Prompt Groups:**
1. Spatial domain identification
2. Cell-type deconvolution
3. Spatially variable gene detection
4. Cell-cell communication
5. Moran's I spatial autocorrelation

**Output:** `data/sensitivity/`

### 7. Cross-System Invocation Comparison

**Scripts:**
- `scripts/cross_system_comparison.py`
- `scripts/cross_system_analysis.py`

**Paper Section:** Results Section 2.4, Figure 4f-g
**Design:** 8 prompts x 3 models x 10 replicates x 2 conditions = 480 new trials (1,200 total across 5 conditions including ChatSpatial ablation data)

Compares tool selection quality across three spatial transcriptomics AI frameworks:

| System | Tools | Architecture |
|--------|-------|-------------|
| ChatSpatial (3 conditions) | 20 MCP tools | Schema-enforced typed parameters |
| STAgent | 8 tools (incl. general code execution) | LangGraph agent with tool retrieval |
| SpatialAgent | 72 specialized tools | LangChain agent with category-based routing |

System contexts (tool catalogs, system prompts, routing instructions) were extracted from each framework's source code and presented to the same LLMs with the same prompts. Tool selection accuracy was evaluated by domain-expert annotation.

**Output:** `data/cross_system/`

### 8. End-to-End Execution Benchmark

**Scripts:**
- `scripts/e2e_benchmark.py`
- `scripts/e2e_benchmark_analysis.py`
- `scripts/e2e_driver_stagent.py` (subprocess driver)
- `scripts/e2e_driver_spatialagent.py` (subprocess driver)

**Paper Section:** Results Section 2.4, Figure 4h-i
**Design:** 3 tasks x 3 systems x 5 replicates = 45 trials
**LLM:** Claude Sonnet 4, T=1.0
**Dataset:** Human lymph node Visium (4,035 spots x 36,601 genes)

Runs all three systems end-to-end on the same dataset with identical prompts:

| Task | ChatSpatial Metric | Concordance Measure |
|------|-------------------|---------------------|
| Spatial domain identification | Leiden clustering | Pairwise ARI |
| SVG detection | Method-dependent | Jaccard@100 |
| Cell-cell communication | CellPhoneDB/CellChat | Jaccard@50 (LR pairs) |

Each system runs through its native API with a 10-minute timeout per trial. Systems are isolated in separate virtual environments.

**Output:** `data/e2e_benchmark/`

### 9. DLPFC Ground-Truth Benchmark

**Scripts:**
- `scripts/dlpfc_prepare_data.py` (data download and preparation)
- `scripts/dlpfc_benchmark.py` (benchmark execution)
- `scripts/dlpfc_benchmark_analysis.py` (analysis)

**Paper Section:** Results Section 2.4, Figure 4j-k
**Design:** 3 samples x 3 systems x 10 replicates = 90 trials
**LLM:** Claude Sonnet 4, T=1.0
**Dataset:** DLPFC Visium (Maynard et al. 2021, spatialLIBD)
**Prompt:** Open-ended: "This is human dorsolateral prefrontal cortex (DLPFC) Visium data. Identify spatial domains in this dataset."

Ground-truth concordance is measured by ARI and NMI against expert-annotated cortical layers (L1-L6 and white matter). Cross-replicate concordance is measured by pairwise ARI. All 95% CIs use hierarchical bootstrap (resampling samples, then replicates within; 10,000 iterations) to account for nested trial structure.

**Output:** `data/dlpfc_benchmark/`

### 10. Case-Study Workflow Concordance

**Script:** `scripts/casestudy_reproducibility.py`
**Paper Section:** Results Section 2.4
**Design:** 2 samples x 2 conditions x 4 models x 10 replicates = 160 trials

Bridges the ablation to the specific case-study workflow by repeating the OSCC CARD deconvolution (Step 2 from Section 2.3.1) under full-schema and no-schema conditions. CARD deconvolution is deterministic given fixed parameters, so any output variation is attributable entirely to schema-condition-induced parameter variation.

**Models:** Gemini 2.5 Flash, Claude Haiku 4.5, GPT-5.4, Claude Sonnet 4.5
**Samples:** OSCC Sample 1 and Sample 9
**Metric:** Pairwise Pearson r of flattened proportion matrices with bootstrap 95% CIs

**Output:** `data/casestudy_reproducibility/`

---

## Running the Scripts

> **Note on script paths:** These scripts were developed and executed in the ChatSpatial development workspace, where they sit alongside the main `code/` package and `benchmarks/` directories. Some scripts (particularly Group 3: end-to-end benchmarks) reference sibling directories for dataset paths and competitor framework installations. When running outside the development workspace, ensure ChatSpatial is installed via pip (`pip install chatspatial`) and adjust data paths as needed. Several scripts share utility functions via sibling imports (e.g., `prompt_sensitivity.py` imports from `ablation_invocation.py`); run them from the repository root (`python scripts/<name>.py`) so that Python resolves these imports correctly.

### Group 1: Static Analysis (no API keys needed)

```bash
python scripts/reproducibility_analysis.py
```

Requires only ChatSpatial installed. Outputs `data/reproducibility_analysis.csv`.

### Group 2: Invocation-Level Experiments (API keys required)

```bash
# Single-model determinism (Gemini only)
python scripts/determinism_experiment.py

# Cross-model determinism (all 3 API keys)
python scripts/determinism_multimodel.py

# Code generation baseline (all 3 API keys)
python scripts/codegen_hallucination.py

# Schema-enforcement ablation: invocation level (all 3 API keys)
python scripts/ablation_invocation.py

# Prompt sensitivity (all 3 API keys)
python scripts/prompt_sensitivity.py

# Cross-system invocation comparison (all 3 API keys)
python scripts/cross_system_comparison.py
```

These scripts make LLM API calls and include incremental JSONL checkpointing for resume-on-failure. Rate limiting and retry logic are built in.

### Group 3: End-to-End Execution (datasets + API keys required)

```bash
# Prepare DLPFC benchmark data (requires R + spatialLIBD)
python scripts/dlpfc_prepare_data.py

# End-to-end execution benchmark
python scripts/e2e_benchmark.py

# DLPFC ground-truth benchmark
python scripts/dlpfc_benchmark.py

# Case-study workflow concordance
python scripts/casestudy_reproducibility.py
```

These scripts require both API keys and local dataset files. They also require ChatSpatial and (for multi-system benchmarks) STAgent and SpatialAgent installed in separate virtual environments.

### Group 4: Analysis Scripts (no API keys needed)

```bash
# Ablation concordance analysis
python scripts/ablation_analysis.py

# Prompt sensitivity analysis
python scripts/prompt_sensitivity_analysis.py

# Cross-system analysis
python scripts/cross_system_analysis.py

# End-to-end benchmark analysis
python scripts/e2e_benchmark_analysis.py

# DLPFC benchmark analysis
python scripts/dlpfc_benchmark_analysis.py

# Case-study effect sizes
python scripts/compute_case_study_stats.py

# Generate supplementary tables
python scripts/build_supplementary_tables.py
```

Analysis scripts read from `data/` and write summary CSVs and text files. They do not make API calls.

---

## Output Data Description

### Top-Level Result Files

| File | Description | Key Columns |
|------|-------------|-------------|
| `reproducibility_analysis.csv` | Per-tool schema constraint statistics | tool, n_params, n_literal, n_bounded, n_described |
| `determinism_experiment.csv` | Per-trial determinism results (single model) | prompt, rep, tool_name, params_json, match |
| `determinism_multimodel.csv` | Per-trial cross-model determinism | model, prompt, rep, tool_name, valid |
| `codegen_hallucination.csv` | Per-trial code generation quality | model, prompt, rep, imports_valid, api_correct |

### Ablation Results (`data/ablation/`)

| File | Description |
|------|-------------|
| `invocation/ablation_invocation_results.csv` | Per-trial validation results across 3 conditions |
| `invocation/ablation_invocation_summary.txt` | Human-readable summary statistics |
| `e2e/ablation_e2e_results.csv` | Per-trial execution results for 3 tasks |
| `analysis/ablation_metrics_bootstrap.csv` | Concordance metrics with 95% bootstrap CIs |
| `analysis/ablation_statistical_tests.csv` | Kruskal-Wallis and pairwise Mann-Whitney tests |
| `analysis/ablation_summary.txt` | Human-readable analysis summary |

### Sensitivity Results (`data/sensitivity/`)

| File | Description |
|------|-------------|
| `prompt_sensitivity_results.csv` | Per-trial results across 5 prompt groups |
| `prompt_sensitivity_summary.txt` | Summary statistics |
| `analysis/sensitivity_cross_paraphrase.csv` | Cross-paraphrase concordance metrics |
| `analysis/sensitivity_summary.txt` | Analysis summary |

### Cross-System Results (`data/cross_system/`)

| File | Description |
|------|-------------|
| `cross_system_results.csv` | Per-trial tool selection results, 5 conditions |
| `cross_system_summary.txt` | Summary statistics |
| `analysis/cross_system_aggregate.csv` | Per-condition aggregates with bootstrap CIs |
| `analysis/cross_system_metrics.csv` | Per condition x model x prompt metrics |
| `analysis/cross_system_summary.txt` | Analysis summary |

### E2E Benchmark Results (`data/e2e_benchmark/`)

| File | Description |
|------|-------------|
| `e2e_benchmark_results.csv` | Per-trial execution results (3 systems x 3 tasks) |
| `e2e_benchmark_summary.txt` | Summary with concordance metrics |

### DLPFC Benchmark Results (`data/dlpfc_benchmark/`)

| File | Description |
|------|-------------|
| `dlpfc_benchmark_results.csv` | Aggregate results per system |
| `dlpfc_benchmark_per_trial.csv` | Per-trial ARI and NMI against ground truth |
| `dlpfc_benchmark_summary.txt` | Summary with hierarchical bootstrap CIs |
| `ground_truth/151507_labels.csv` | Expert annotations for sample 151507 |
| `ground_truth/151669_labels.csv` | Expert annotations for sample 151669 |
| `ground_truth/151673_labels.csv` | Expert annotations for sample 151673 |

### Case-Study Concordance (`data/casestudy_reproducibility/`)

| File | Description |
|------|-------------|
| `casestudy_repro_metrics.csv` | Pairwise Pearson r across conditions and models |
| `casestudy_repro_summary.txt` | Summary statistics |

---

## Supplementary Tables

| Table | File | Description |
|-------|------|-------------|
| S1 | `Supplementary_Table_1_Effect_Sizes.csv` | Sample-level effect sizes for OSCC TC vs. LE enrichment (Wilcoxon signed-rank, log2 fold change, 95% CIs) |
| S2 | `Supplementary_Table_2_Integrated_Methods.csv` | Complete catalog of 65 integrated analytical methods across 15 categories, with versions and citations |
| S3 | `Supplementary_Table_3_Test_Scenarios.csv` | 31 predefined test scenarios covering data handling (5), core analysis (11), conversational workflows (5), scalability (7), and known limitations (3) |
| S4 | `Supplementary_Table_4_AI_Agent_Comparison.csv` | Feature comparison of ChatSpatial, STAgent, and SpatialAgent across architecture, tool count, validation, and error handling |

---

## Reproducibility Notes

- All random seeds are fixed (seed=42 for bootstrap resampling) where applicable.
- LLM experiments use T=1.0 (temperature) to maximize stochastic variation and stress-test consistency.
- All scripts include incremental JSONL checkpointing: if a run is interrupted, re-running the script resumes from the last completed trial.
- Bootstrap confidence intervals use 10,000 resamples throughout.
- The DLPFC benchmark uses hierarchical bootstrap (samples, then replicates within) to account for nested structure.
- Cross-system benchmarks isolate each system in separate virtual environments with independent timeouts (10-15 minutes per trial).

---

## Related Repositories

- **ChatSpatial Package:** [github.com/cafferychen777/ChatSpatial](https://github.com/cafferychen777/ChatSpatial) — main package with MCP server, tool implementations, and documentation
- **Documentation:** [docs.cafferyang.com](https://docs.cafferyang.com/) — comprehensive user guide and API reference
- **Docker Image:** `ghcr.io/cafferychen777/chatspatial:v1.2.6` — pre-built environment with all dependencies

## License

MIT
