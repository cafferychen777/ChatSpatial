#!/usr/bin/env python3
"""
Cross-System Invocation Comparison

Compares tool selection quality across three spatial transcriptomics frameworks:
  - ChatSpatial (full_schema, bare_schema, no_schema — reuse existing ablation data)
  - STAgent (system prompt + 8 tools — new API calls)
  - SpatialAgent (system prompt + 72 tools — new API calls)

Trial matrix (new): 8 prompts × 3 models × 10 reps × 2 conditions = 480 API calls.
ChatSpatial data (720 existing trials) loaded from ablation results for comparison.

Same prompts, LLMs, T=1.0, and JSON response format across all conditions.
The only variable is the system context (tool catalog + guidance).
"""

import csv
import json
import os
import time
from collections import Counter, defaultdict
from pathlib import Path

import requests

# Reuse infrastructure from ablation experiment
from ablation_invocation import (
    CALLERS as _BASE_CALLERS,
    DELAY,
    MODELS,
    N_TRIALS,
    RESPONSE_INSTRUCTION,
    TEMPERATURE,
    TEST_PROMPTS,
    append_raw,
    compute_consistency,
    load_completed_trials,
    parse_json,
    GEMINI_KEY,
    ANTHROPIC_KEY,
    OPENAI_API_URL,
    OPENAI_MODEL,
    MAX_RETRIES,
    RETRY_BACKOFF,
)


# ---------------------------------------------------------------------------
# Override API callers with higher max output tokens
# STAgent/SpatialAgent contexts cause LLMs to generate long code strings
# in python_repl_tool parameters, exceeding the ablation's 512-token limit.
# ---------------------------------------------------------------------------

_MAX_OUTPUT = 4096  # 8x ablation default — accommodates code-gen responses


def _call_gemini_extended(system: str, prompt: str) -> str | None:
    """Gemini caller with extended maxOutputTokens."""
    for attempt in range(MAX_RETRIES):
        try:
            payload = {
                "contents": [
                    {"role": "user", "parts": [{"text": system + RESPONSE_INSTRUCTION + "\nUser: " + prompt}]}
                ],
                "generationConfig": {"temperature": TEMPERATURE, "maxOutputTokens": 8192},
            }
            r = requests.post(
                f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={GEMINI_KEY}",
                json=payload, timeout=90,
            )
            r.raise_for_status()
            return r.json()["candidates"][0]["content"]["parts"][0]["text"]
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF * (attempt + 1))
            else:
                print(f"  Gemini error after {MAX_RETRIES} retries: {e}")
                return None


def _call_anthropic_extended(system: str, prompt: str) -> str | None:
    """Anthropic caller with extended max_tokens."""
    for attempt in range(MAX_RETRIES):
        try:
            payload = {
                "model": "claude-haiku-4-5-20251001",
                "max_tokens": _MAX_OUTPUT,
                "temperature": TEMPERATURE,
                "system": system + RESPONSE_INSTRUCTION,
                "messages": [{"role": "user", "content": prompt}],
            }
            r = requests.post(
                "https://api.anthropic.com/v1/messages",
                headers={
                    "x-api-key": ANTHROPIC_KEY,
                    "anthropic-version": "2023-06-01",
                    "content-type": "application/json",
                },
                json=payload, timeout=90,
            )
            r.raise_for_status()
            return r.json()["content"][0]["text"]
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF * (attempt + 1))
            else:
                print(f"  Anthropic error after {MAX_RETRIES} retries: {e}")
                return None


def _call_openai_extended(system: str, prompt: str) -> str | None:
    """OpenAI GPT caller with extended max_tokens."""
    openai_key = os.environ.get("OPENAI_API_KEY", "")
    if not openai_key or not OPENAI_API_URL:
        print(
            "  GPT-5.4 gateway not configured; set OPENAI_API_KEY and "
            "OPENAI_API_URL to run this model."
        )
        return None
    for attempt in range(MAX_RETRIES):
        try:
            payload = {
                "model": OPENAI_MODEL,
                "max_tokens": _MAX_OUTPUT,
                "temperature": TEMPERATURE,
                "system": system + RESPONSE_INSTRUCTION,
                "messages": [{"role": "user", "content": prompt}],
            }
            r = requests.post(
                OPENAI_API_URL,
                headers={
                    "x-api-key": openai_key,
                    "anthropic-version": "2023-06-01",
                    "content-type": "application/json",
                },
                json=payload, timeout=120,
            )
            r.raise_for_status()
            content_blocks = r.json()["content"]
            for block in content_blocks:
                if block["type"] == "text":
                    return block["text"]
            return None
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF * (attempt + 1))
            else:
                print(f"  OpenAI/GPT error after {MAX_RETRIES} retries: {e}")
                return None


CALLERS = {
    "gemini": _call_gemini_extended,
    "anthropic": _call_anthropic_extended,
    "openai": _call_openai_extended,
}

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

N_REPS = N_TRIALS  # 10, same as ablation
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR.parent / "data" / "cross_system"
DATA_DIR.mkdir(parents=True, exist_ok=True)
RAW_PATH = DATA_DIR / "cross_system_raw.jsonl"
CSV_PATH = DATA_DIR / "cross_system_results.csv"
SUMMARY_PATH = DATA_DIR / "cross_system_summary.txt"
ABLATION_RAW = SCRIPT_DIR.parent / "data" / "ablation" / "invocation" / "ablation_invocation_raw.jsonl"

# ---------------------------------------------------------------------------
# STAgent System Context
# Extracted from STAgent/src/prompt.py (tool descriptions + pipeline).
# Dataset-specific context and repetitive execution instructions removed
# to isolate tool selection guidance only.
# ---------------------------------------------------------------------------

STAGENT_CONTEXT = """\
Spatial Transcriptomics AI Agent

This AI agent specializes in analyzing spatial transcriptomics data through a \
systematic pipeline. It utilizes a set of tools for data exploration, \
visualization, and biological interpretation.

---

Available Tools:
1. python_repl_tool:
   - Executes Python code in a live Python shell
   - Returns printed outputs and generated visualizations
   - Input: Valid Python commands
   - Output: Execution results and plot file paths

2. google_scholar_search:
   - Retrieves academic articles and summaries
   - Input: Research topic or biological query
   - Output: Article titles, authors, and summaries
   - Usage: For literature-backed biological interpretation

3. squidpy_rag_agent:
   - Provides guidance on Squidpy usage and spatial analysis workflows
   - Input: Questions about Squidpy functions
   - Output: Code examples and explanations
   - Usage: For spatial analysis using Squidpy (e.g., spatial_autocorr, nhood_enrichment)

4. visualize_umap:
   - Creates UMAP plots colored by cell type
   - Shows clustering patterns of different cell populations

5. visualize_cell_type_composition:
   - Shows cell type proportions across samples
   - Displays stacked bar plots and heatmaps of composition changes

6. visualize_spatial_cell_type_map:
   - Creates spatial scatter plots of cell types in tissue context
   - Shows cell locations and spatial organization

7. visualize_cell_cell_interaction:
   - Analyzes cell type interaction patterns
   - Output: Neighborhood enrichment heatmaps
   - Reveals spatial relationships between cell types

8. report_tool:
   - Generates a summary report of the analysis

---

Pipeline Instructions:
1. Dimensionality Reduction: Use visualize_umap for cell type clustering
2. Cell Type Composition: Use visualize_cell_type_composition for proportions
3. Spatial Distribution: Use visualize_spatial_cell_type_map for tissue context
4. Cell-Cell Interaction: Use visualize_cell_cell_interaction for neighborhood patterns
5. Report: Use report_tool to summarize the analysis
"""

# ---------------------------------------------------------------------------
# SpatialAgent Tool Catalog (72 tools from 7 modules)
# Format mirrors Tool.to_text() from spatialagent/agent/tool_system.py.
# Condensed: name + 1-sentence description + parameter names/types.
# ---------------------------------------------------------------------------

_SPATIALAGENT_TOOLS = """\
## Database & Reference Tools

Tool: search_panglao
Description: Search PanglaoDB for marker genes of given cell types.
Parameters:
  - cell_types (str): Cell types to search — comma-separated or CSV path
  - organism (str): "Hs" (human) or "Mm" (mouse)
  - tissue (str): Target tissue context

Tool: search_czi_datasets
Description: Search CZI CELLxGENE Census for reference single-cell datasets \
using embedding-based ranking.
Parameters:
  - query (str): Condition/context query (e.g., "breast cancer")
  - n_datasets (int): Number of top datasets to return
  - organism (str): Filter by organism
  - tissue (str): Filter by tissue

Tool: search_cellmarker2
Description: Search CellMarker2 database for cell type marker genes.
Parameters:
  - cell_types (str): Cell types to search
  - organism (str): "Human" or "Mouse"
  - tissue (str): Target tissue

Tool: extract_czi_markers
Description: Download and process CZI reference datasets to extract cell types \
and marker genes.
Parameters:
  - save_path (str): Experiment directory
  - dataset_id (str): CZI dataset ID
  - iter_round (int): Iteration round (1-3)
  - organism (str): "Homo sapiens" or "Mus musculus"

Tool: download_czi_reference
Description: Download CZI reference scRNA-seq data for Harmony label transfer.
Parameters:
  - dataset_id (str): CZI dataset ID
  - organism (str): Species
  - save_path (str): Experiment directory

Tool: query_tissue_expression
Description: Query ARCHS4 for tissue-specific gene expression (median TPM).
Parameters:
  - gene (str): Gene symbol (e.g., "GFAP")
  - top_k (int): Number of top tissues to return

Tool: query_celltype_genesets
Description: Query Enrichr/PanglaoDB for cell type-specific gene sets in a tissue.
Parameters:
  - tissue (str): Tissue type (e.g., "brain")
  - top_k (int): Number of gene sets

Tool: validate_genes_expression
Description: Validate that genes are expressed in target tissue using ARCHS4.
Parameters:
  - genes (str): Comma-separated gene symbols
  - target_tissue (str): Target tissue

Tool: query_disease_genes
Description: Query disease-associated genes from GWAS Catalog and OpenTargets.
Parameters:
  - disease (str): Disease or trait to search
  - source (str): "opentargets", "gwas", or "all"
  - max_genes (int): Maximum genes to return

## Literature & Web Tools

Tool: query_pubmed
Description: Query PubMed for papers matching a search query.
Parameters:
  - query (str): PubMed search query
  - max_papers (int): Maximum papers to retrieve

Tool: query_arxiv
Description: Query arXiv for preprints matching a search query.
Parameters:
  - query (str): arXiv search query
  - max_papers (int): Maximum papers

Tool: search_semantic_scholar
Description: Search Semantic Scholar for academic papers with citation counts.
Parameters:
  - query (str): Search query
  - max_papers (int): Maximum papers

Tool: extract_url_content
Description: Extract text content from a webpage URL.
Parameters:
  - url (str): URL to extract
  - max_chars (int): Maximum characters

Tool: extract_pdf_content
Description: Extract text from a PDF file URL.
Parameters:
  - url (str): PDF URL

Tool: fetch_supplementary_from_doi
Description: Fetch paper content and supplementary materials from DOI.
Parameters:
  - doi (str): Paper DOI

Tool: web_search
Description: Search the web using server-side search.
Parameters:
  - query (str): Search query

## Preprocessing & Data Integration

Tool: preprocess_spatial_data
Description: Preprocess spatial transcriptomics data: QC filtering, \
normalization, PCA, UMAP, neighbors graph.
Parameters:
  - adata_path (str): Path to raw h5ad file
  - save_path (str): Output directory

Tool: harmony_transfer_labels
Description: Transfer cell type labels from scRNA-seq reference to spatial \
data using Harmony integration and MLP classifier.
Parameters:
  - adata_path (str): Preprocessed spatial data
  - ref_path (str): CZI reference scRNA data
  - save_path (str): Output directory

Tool: scanpy_ingest
Description: Transfer annotations from reference to query dataset using \
Scanpy ingest.
Parameters:
  - query_path (str): Query dataset path
  - ref_path (str): Reference dataset path
  - cell_type_key (str): Cell type column in reference

Tool: scanpy_bbknn
Description: Batch-balanced KNN integration for multi-batch spatial data.
Parameters:
  - adata_path (str): Multi-batch data path
  - batch_key (str): Batch column name

Tool: totalvi_integration
Description: Integrate RNA and protein data using totalVI (CITE-seq).
Parameters:
  - adata_path (str): Data with RNA and protein
  - protein_key (str): Protein data key in obsm

Tool: multivi_integration
Description: Integrate RNA and ATAC data using MultiVI (multiome).
Parameters:
  - adata_path (str): Data with RNA and ATAC
  - atac_key (str): ATAC data key in obsm

Tool: mofa_integration
Description: Multi-omics factor analysis using MOFA+.
Parameters:
  - adata_path (str): Multi-omics data
  - modalities (list): Modalities to analyze

## Spatial Clustering

Tool: run_utag_clustering
Description: Run UTAG spatial clustering to identify tissue niches by \
combining gene expression with spatial information via message passing.
Parameters:
  - adata_path (str): Annotated spatial data
  - save_path (str): Output directory
  - slide_key (str): Sample/slide column for per-sample UTAG
  - resolutions (list): Clustering resolutions (default: [0.05, 0.1, 0.3])

Tool: spagcn_clustering
Description: Run SpaGCN graph convolutional network for spatial domain \
identification.
Parameters:
  - adata_path (str): Spatial data
  - n_clusters (int): Number of spatial domains
  - save_path (str): Output directory

Tool: graphst_clustering
Description: Run GraphST graph-based method for spatial domain detection.
Parameters:
  - adata_path (str): Spatial data
  - n_clusters (int): Number of spatial domains

Tool: scanpy_score_genes
Description: Score cells/spots based on expression of a gene signature \
using Scanpy.
Parameters:
  - adata_path (str): Spatial data
  - genes (str): Comma-separated gene signature
  - score_name (str): Score column name

## Spatial Deconvolution

Tool: destvi_deconvolution
Description: Run DestVI for multi-resolution spatial deconvolution to \
estimate cell type composition.
Parameters:
  - adata_path (str): Spatial data
  - sc_adata_path (str): scRNA-seq reference

Tool: cell2location_mapping
Description: Run Cell2location to map cell type signatures from reference \
to spatial coordinates.
Parameters:
  - adata_path (str): Spatial data
  - sc_adata_path (str): scRNA-seq reference
  - cell_type_key (str): Cell type column in reference

Tool: stereoscope_deconvolution
Description: Run Stereoscope to estimate cell type abundances at each \
spatial spot.
Parameters:
  - adata_path (str): Spatial data
  - sc_adata_path (str): scRNA-seq reference

Tool: gimvi_imputation
Description: Run gimVI to impute gene expression for genes not measured \
in spatial data.
Parameters:
  - adata_path (str): Spatial data
  - sc_adata_path (str): scRNA-seq reference

## Tangram Spatial Mapping

Tool: tangram_preprocess
Description: Preprocess scRNA-seq and spatial data for Tangram cell mapping.
Parameters:
  - sc_adata_path (str): scRNA-seq data
  - sp_adata_path (str): Spatial data

Tool: tangram_map_cells
Description: Map single cells onto spatial locations using Tangram.
Parameters:
  - sc_adata_path (str): scRNA-seq data
  - sp_adata_path (str): Spatial data

Tool: tangram_project_annotations
Description: Project cell type annotations from scRNA-seq onto spatial data.
Parameters:
  - mapping_path (str): Tangram mapping results

Tool: tangram_project_genes
Description: Project gene expression from scRNA-seq onto spatial locations.
Parameters:
  - mapping_path (str): Tangram mapping results
  - genes (str): Comma-separated gene names

Tool: tangram_evaluate
Description: Evaluate Tangram mapping quality metrics.
Parameters:
  - mapping_path (str): Tangram mapping results

## Squidpy Spatial Analysis

Tool: squidpy_spatial_neighbors
Description: Build spatial neighbors graph for downstream spatial statistics.
Parameters:
  - adata_path (str): Spatial data
  - max_distance (float): Maximum neighbor distance

Tool: squidpy_nhood_enrichment
Description: Compute neighborhood enrichment to test non-random \
co-localization between cell types.
Parameters:
  - adata_path (str): Spatial data
  - cell_type_key (str): Cell type column

Tool: squidpy_co_occurrence
Description: Compute co-occurrence probability between cell types as \
function of distance.
Parameters:
  - adata_path (str): Spatial data
  - cell_type_key (str): Cell type column

Tool: squidpy_spatial_autocorr
Description: Compute spatial autocorrelation (Moran's I or Geary's C) \
to identify spatially variable genes.
Parameters:
  - adata_path (str): Spatial data
  - genes (str): Comma-separated gene names (default: all genes)

Tool: squidpy_ripley
Description: Compute Ripley's statistics for spatial point pattern analysis \
of cell types.
Parameters:
  - adata_path (str): Spatial data
  - cell_type_key (str): Cell type column

Tool: squidpy_centrality
Description: Compute centrality scores for cell type clusters in spatial graph.
Parameters:
  - adata_path (str): Spatial data
  - cell_type_key (str): Cell type column

Tool: squidpy_interaction_matrix
Description: Compute pairwise interaction matrix between cell type clusters.
Parameters:
  - adata_path (str): Spatial data
  - cell_type_key (str): Cell type column

Tool: squidpy_ligrec
Description: Run receptor-ligand analysis using Squidpy and OmniPath database.
Parameters:
  - adata_path (str): Spatial data
  - cell_type_key (str): Cell type column
  - organism (str): "human" or "mouse"

## Cell-Cell Communication: CellPhoneDB

Tool: cellphonedb_prepare
Description: Prepare AnnData for CellPhoneDB analysis format.
Parameters:
  - adata_path (str): Spatial data with cell types
  - cell_type_key (str): Cell type column

Tool: cellphonedb_analysis
Description: Run CellPhoneDB ligand-receptor interaction analysis \
(statistical or simple mode).
Parameters:
  - adata_path (str): Spatial data
  - cell_type_key (str): Cell type column
  - method (str): "statistical" or "simple"

Tool: cellphonedb_degs_analysis
Description: Run CellPhoneDB DEG-based analysis for condition-specific \
interactions.
Parameters:
  - adata_path (str): Spatial data with conditions
  - condition_column (str): Condition column
  - cell_type_key (str): Cell type column

Tool: cellphonedb_filter
Description: Search and filter CellPhoneDB results by gene, cell type, \
or pathway.
Parameters:
  - cpdb_results_path (str): CellPhoneDB results path
  - query (str): Search query

Tool: cellphonedb_plot
Description: Visualize CellPhoneDB results as heatmaps and network plots.
Parameters:
  - cpdb_results_path (str): CellPhoneDB results path
  - top_k (int): Top interactions to plot

## Cell-Cell Communication: LIANA

Tool: liana_tensor
Description: Run LIANA + Tensor-Cell2Cell for multi-sample tensor \
decomposition of cell-cell interactions.
Parameters:
  - adata_path (str): AnnData file
  - sample_key (str): Sample column
  - condition_key (str): Condition column
  - cell_type_key (str): Cell type column
  - organism (str): "human" or "mouse"

Tool: liana_inference
Description: Run LIANA rank aggregate for consensus ligand-receptor \
inference across multiple algorithms.
Parameters:
  - adata_path (str): Spatial data
  - cell_type_key (str): Cell type column
  - organism (str): "human" or "mouse"

Tool: liana_spatial
Description: Run LIANA spatial bivariate analysis of ligand-receptor \
expression across coordinates.
Parameters:
  - adata_path (str): Spatial data
  - cell_type_key (str): Cell type column
  - organism (str): "human" or "mouse"

Tool: liana_misty
Description: Run LIANA MISTy to learn how gene expression depends on \
spatial context.
Parameters:
  - adata_path (str): Spatial data
  - cell_type_key (str): Cell type column

Tool: liana_plot
Description: Visualize LIANA ligand-receptor interaction results.
Parameters:
  - liana_results_path (str): LIANA results path
  - top_k (int): Top interactions to plot

## Trajectory Analysis

Tool: scvelo_velocity
Description: Compute RNA velocity using scVelo from spliced/unspliced counts.
Parameters:
  - adata_path (str): Data with spliced/unspliced counts

Tool: scvelo_velocity_embedding
Description: Project RNA velocity onto embedding and create velocity \
stream plot.
Parameters:
  - adata_path (str): Data with velocity computed
  - embedding (str): "umap" or "pca"

Tool: cellrank_terminal_states
Description: Identify terminal/macrostates using CellRank from velocity \
or transition information.
Parameters:
  - adata_path (str): Data with velocity computed

Tool: cellrank_fate_probabilities
Description: Compute cell fate probabilities towards each terminal state \
using CellRank.
Parameters:
  - adata_path (str): Data with terminal states identified

Tool: paga_trajectory
Description: Compute PAGA partition-based graph abstraction for trajectory \
analysis.
Parameters:
  - adata_path (str): Annotated data
  - cluster_key (str): Cell type/cluster column

## Analytics & Summarization

Tool: aggregate_gene_voting
Description: Aggregate marker genes across cells/niches using LLM-based \
voting.
Parameters:
  - adata_path (str): Annotated data
  - group_by (str): Grouping column

Tool: infer_dynamics
Description: Differential expression analysis (Wilcoxon test) between \
two conditions.
Parameters:
  - adata_path (str): Spatial data
  - condition_column (str): Condition column
  - condition1 (str): First condition
  - condition2 (str): Second condition

Tool: summarize_conditions
Description: Summarize cell type distributions across experimental conditions.
Parameters:
  - adata_path (str): Spatial data
  - condition_key (str): Condition column

Tool: summarize_celltypes
Description: Summarize cell type distributions and marker genes in the dataset.
Parameters:
  - adata_path (str): Annotated data
  - cell_type_key (str): Cell type column

Tool: summarize_tissue_regions
Description: Summarize tissue regions and their cell type compositions.
Parameters:
  - adata_path (str): Spatial data
  - region_key (str): Region column

## Interpretation

Tool: annotate_cell_types
Description: Annotate cell type clusters using hierarchical two-level \
LLM-based approach.
Parameters:
  - adata_path (str): Preprocessed spatial data
  - transferred_celltype (str): Path to transferred cell type labels
  - data_info (str): Dataset description (e.g., "human heart MERFISH")

Tool: annotate_tissue_niches
Description: Annotate tissue niches using spatial domain and cell type \
composition with LLM.
Parameters:
  - adata_path (str): Spatial data with cell types
  - utag_csv (str): UTAG clustering results
  - data_info (str): Dataset description

Tool: interpret_figure
Description: Interpret a plotted figure using vision LLM for biological \
insights.
Parameters:
  - image_path (str): Figure image path
  - context (str): What was plotted

## Code Execution

Tool: execute_python
Description: Execute Python code in a stateful environment with \
pre-imported scientific libraries (numpy, pandas, scanpy, squidpy, \
matplotlib, seaborn).
Parameters:
  - code (str): Python code to execute

Tool: execute_bash
Description: Execute a Bash shell command.
Parameters:
  - command (str): Bash command

## Utility

Tool: inspect_tool_code
Description: Retrieve source code of a predefined tool for understanding \
or debugging.
Parameters:
  - tool_name (str): Tool name to inspect

Tool: report_subagent
Description: Generate publication-quality research report from all analysis \
artifacts.
Parameters:
  - user_query (str): Research question
  - data_info (str): Dataset description

Tool: verification_subagent
Description: Verify analysis conclusions against saved evidence (figures, \
data, observations).
Parameters:
  - user_query (str): Research question
  - conclusions (str): Conclusions to verify
  - data_info (str): Dataset description
"""

# SpatialAgent system context = role + tool catalog + analysis guidance.
# Execution-specific instructions (<act>/<conclude> tags) stripped;
# replaced by shared RESPONSE_INSTRUCTION for JSON output.
SPATIALAGENT_CONTEXT = f"""\
You are a computational biologist specialized in spatial transcriptomics \
analysis. Your task is to select the most appropriate tool and parameters \
for the given analysis request.

# Available Tools
{_SPATIALAGENT_TOOLS}
# Analysis Capabilities
- Reference databases: PanglaoDB, CellMarker2, CZI CELLxGENE Census
- Data processing: Scanpy, Harmony, quality control, normalization
- Statistical analysis: Leiden clustering, UTAG spatial clustering
- Cell-cell interactions: LIANA, CellPhoneDB, Tensor-Cell2cell, Squidpy ligrec
- Spatial analysis: Squidpy (Moran's I, Geary's C, Ripley's, neighborhood \
enrichment, co-occurrence)
- Deconvolution: DestVI, Cell2location, Stereoscope, Tangram, gimVI
- Spatial clustering: UTAG, SpaGCN, GraphST
- Trajectory: scVelo RNA velocity, CellRank, PAGA
- LLM interpretation: Hierarchical cell type annotation, tissue niche \
annotation, figure interpretation
- Code execution: Custom Python and Bash with persistent state

# Platform-Aware Analysis
- Spot-based platforms (Visium, Slide-seq, ST): Each spot contains multiple \
cells. Use cell type DECONVOLUTION (Cell2location, DestVI, Stereoscope) to \
estimate proportions.
- Single-cell platforms (MERFISH, Xenium, CosMx, SeqFISH): Each observation \
is one cell. Use cell type ANNOTATION via Harmony transfer or clustering.
"""

# ---------------------------------------------------------------------------
# Expected Tool Mappings (for scoring)
# Each prompt_id -> list of acceptable tool names for each system.
# Scoring: did the LLM select a tool the system was designed to use?
# ---------------------------------------------------------------------------

STAGENT_EXPECTED: dict[str, list[str]] = {
    # STAgent has no dedicated spatial domain tool -> code execution
    "spatial_domains": ["python_repl_tool", "squidpy_rag_agent"],
    # STAgent has no deconvolution tool -> code execution
    "deconvolution_rctd": ["python_repl_tool"],
    # SVGs via Squidpy spatial_autocorr -> RAG or code
    "svg_detection": ["python_repl_tool", "squidpy_rag_agent"],
    # STAgent has neighborhood enrichment viz but not L-R analysis
    "cellchat": ["python_repl_tool", "visualize_cell_cell_interaction"],
    # Moran's I via Squidpy -> RAG or code
    "moran_i": ["python_repl_tool", "squidpy_rag_agent"],
    # Direct match to spatial cell type visualization
    "spatial_plot": ["visualize_spatial_cell_type_map", "python_repl_tool"],
    # No trajectory tool -> code execution
    "trajectory": ["python_repl_tool"],
    # No CNV tool -> code execution
    "cnv": ["python_repl_tool"],
}

SPATIALAGENT_EXPECTED: dict[str, list[str]] = {
    # Spatial clustering tools
    "spatial_domains": [
        "run_utag_clustering", "spagcn_clustering", "graphst_clustering",
    ],
    # Deconvolution tools (prompt asks RCTD but SpatialAgent has DestVI etc.)
    "deconvolution_rctd": [
        "destvi_deconvolution", "cell2location_mapping",
        "stereoscope_deconvolution",
    ],
    # SVG detection via Squidpy
    "svg_detection": ["squidpy_spatial_autocorr"],
    # Cell communication tools (prompt asks CellChat; SpatialAgent has LIANA, CellPhoneDB)
    "cellchat": [
        "liana_inference", "liana_spatial", "cellphonedb_analysis",
        "cellphonedb_degs_analysis", "squidpy_ligrec",
    ],
    # Moran's I via Squidpy
    "moran_i": ["squidpy_spatial_autocorr"],
    # No dedicated visualization tool -> code execution
    "spatial_plot": ["execute_python"],
    # Trajectory tools
    "trajectory": [
        "scvelo_velocity", "cellrank_terminal_states",
        "paga_trajectory", "scvelo_velocity_embedding",
        "cellrank_fate_probabilities",
    ],
    # No CNV tool -> code execution
    "cnv": ["execute_python"],
}

# ---------------------------------------------------------------------------
# Tool Name Resolution
# ---------------------------------------------------------------------------

STAGENT_TOOL_NAMES = {
    "python_repl_tool", "google_scholar_search", "squidpy_rag_agent",
    "visualize_umap", "visualize_cell_type_composition",
    "visualize_spatial_cell_type_map", "visualize_cell_cell_interaction",
    "report_tool",
}

_STAGENT_ALIASES: dict[str, str] = {
    "python_repl": "python_repl_tool",
    "python": "python_repl_tool",
    "repl": "python_repl_tool",
    "code": "python_repl_tool",
    "execute_code": "python_repl_tool",
    "execute_python": "python_repl_tool",
    "run_python": "python_repl_tool",
    "code_execution": "python_repl_tool",
    "google_scholar": "google_scholar_search",
    "scholar_search": "google_scholar_search",
    "squidpy_rag": "squidpy_rag_agent",
    "squidpy": "squidpy_rag_agent",
    "rag_agent": "squidpy_rag_agent",
    "umap": "visualize_umap",
    "cell_type_composition": "visualize_cell_type_composition",
    "composition": "visualize_cell_type_composition",
    "spatial_cell_type_map": "visualize_spatial_cell_type_map",
    "spatial_map": "visualize_spatial_cell_type_map",
    "cell_cell_interaction": "visualize_cell_cell_interaction",
    "interaction": "visualize_cell_cell_interaction",
    "nhood_enrichment": "visualize_cell_cell_interaction",
    "report": "report_tool",
}

SPATIALAGENT_TOOL_NAMES = {
    # Database (9)
    "search_panglao", "search_czi_datasets", "search_cellmarker2",
    "extract_czi_markers", "download_czi_reference",
    "query_tissue_expression", "query_celltype_genesets",
    "validate_genes_expression", "query_disease_genes",
    # Literature (7)
    "query_pubmed", "query_arxiv", "search_semantic_scholar",
    "extract_url_content", "extract_pdf_content",
    "fetch_supplementary_from_doi", "web_search",
    # Preprocessing & Integration (7)
    "preprocess_spatial_data", "harmony_transfer_labels",
    "scanpy_ingest", "scanpy_bbknn",
    "totalvi_integration", "multivi_integration", "mofa_integration",
    # Spatial Clustering (4)
    "run_utag_clustering", "spagcn_clustering", "graphst_clustering",
    "scanpy_score_genes",
    # Deconvolution (4)
    "destvi_deconvolution", "cell2location_mapping",
    "stereoscope_deconvolution", "gimvi_imputation",
    # Tangram (5)
    "tangram_preprocess", "tangram_map_cells",
    "tangram_project_annotations", "tangram_project_genes",
    "tangram_evaluate",
    # Squidpy (8)
    "squidpy_spatial_neighbors", "squidpy_nhood_enrichment",
    "squidpy_co_occurrence", "squidpy_spatial_autocorr",
    "squidpy_ripley", "squidpy_centrality",
    "squidpy_interaction_matrix", "squidpy_ligrec",
    # CellPhoneDB (5)
    "cellphonedb_prepare", "cellphonedb_analysis",
    "cellphonedb_degs_analysis", "cellphonedb_filter",
    "cellphonedb_plot",
    # LIANA (5)
    "liana_tensor", "liana_inference", "liana_spatial",
    "liana_misty", "liana_plot",
    # Trajectory (5)
    "scvelo_velocity", "scvelo_velocity_embedding",
    "cellrank_terminal_states", "cellrank_fate_probabilities",
    "paga_trajectory",
    # Analytics (5)
    "aggregate_gene_voting", "infer_dynamics",
    "summarize_conditions", "summarize_celltypes",
    "summarize_tissue_regions",
    # Interpretation (3)
    "annotate_cell_types", "annotate_tissue_niches", "interpret_figure",
    # Coding (2)
    "execute_python", "execute_bash",
    # Utility (3)
    "inspect_tool_code", "report_subagent", "verification_subagent",
}

_SPATIALAGENT_ALIASES: dict[str, str] = {
    "python_repl": "execute_python",
    "python_repl_tool": "execute_python",
    "python": "execute_python",
    "run_python": "execute_python",
    "code_execution": "execute_python",
    "bash": "execute_bash",
    "run_bash": "execute_bash",
    "utag": "run_utag_clustering",
    "utag_clustering": "run_utag_clustering",
    "spagcn": "spagcn_clustering",
    "graphst": "graphst_clustering",
    "destvi": "destvi_deconvolution",
    "cell2location": "cell2location_mapping",
    "stereoscope": "stereoscope_deconvolution",
    "gimvi": "gimvi_imputation",
    "cellphonedb": "cellphonedb_analysis",
    "liana": "liana_inference",
    "spatial_autocorr": "squidpy_spatial_autocorr",
    "morans_i": "squidpy_spatial_autocorr",
    "moran_i": "squidpy_spatial_autocorr",
    "moran": "squidpy_spatial_autocorr",
    "nhood_enrichment": "squidpy_nhood_enrichment",
    "neighborhood_enrichment": "squidpy_nhood_enrichment",
    "co_occurrence": "squidpy_co_occurrence",
    "ripley": "squidpy_ripley",
    "centrality": "squidpy_centrality",
    "interaction_matrix": "squidpy_interaction_matrix",
    "ligrec": "squidpy_ligrec",
    "scvelo": "scvelo_velocity",
    "cellrank": "cellrank_terminal_states",
    "paga": "paga_trajectory",
    "preprocess": "preprocess_spatial_data",
    "harmony": "harmony_transfer_labels",
    "tangram": "tangram_map_cells",
    "pubmed": "query_pubmed",
    "arxiv": "query_arxiv",
    "annotate_cells": "annotate_cell_types",
    "annotate_niches": "annotate_tissue_niches",
}


def resolve_stagent_tool(raw_name: str | None) -> str | None:
    """Map raw tool name to STAgent canonical name."""
    if raw_name is None:
        return None
    name = raw_name.lower().strip().replace("-", "_")
    if name in STAGENT_TOOL_NAMES:
        return name
    return _STAGENT_ALIASES.get(name)


def resolve_spatialagent_tool(raw_name: str | None) -> str | None:
    """Map raw tool name to SpatialAgent canonical name."""
    if raw_name is None:
        return None
    name = raw_name.lower().strip().replace("-", "_")
    if name in SPATIALAGENT_TOOL_NAMES:
        return name
    return _SPATIALAGENT_ALIASES.get(name)


# ---------------------------------------------------------------------------
# Conditions and dispatch tables
# ---------------------------------------------------------------------------

CONDITIONS = {
    "stagent_context": STAGENT_CONTEXT,
    "spatialagent_context": SPATIALAGENT_CONTEXT,
}

RESOLVERS = {
    "stagent_context": resolve_stagent_tool,
    "spatialagent_context": resolve_spatialagent_tool,
}

EXPECTED_TOOLS = {
    "stagent_context": STAGENT_EXPECTED,
    "spatialagent_context": SPATIALAGENT_EXPECTED,
}

# ---------------------------------------------------------------------------
# Main Experiment
# ---------------------------------------------------------------------------

def run_experiment():
    print("Cross-System Invocation Comparison")
    print(f"Models: {list(MODELS.keys())}")
    print(f"New conditions: {list(CONDITIONS.keys())}")
    print(f"Prompts: {len(TEST_PROMPTS)}, Reps: {N_REPS}")
    total_new = len(CONDITIONS) * len(MODELS) * len(TEST_PROMPTS) * N_REPS
    print(f"Total new trials: {total_new}")
    print(f"SpatialAgent tools: {len(SPATIALAGENT_TOOL_NAMES)}")
    print(f"STAgent tools: {len(STAGENT_TOOL_NAMES)}")
    print()

    completed = load_completed_trials(RAW_PATH)
    print(f"Resuming: {len(completed)} trials already completed\n")

    count = 0
    skipped = 0

    for condition, context in CONDITIONS.items():
        resolver = RESOLVERS[condition]
        expected = EXPECTED_TOOLS[condition]

        for model_name, provider in MODELS.items():
            caller = CALLERS[provider]
            print(f"\n{'='*60}")
            print(f"Condition: {condition} | Model: {model_name}")
            print(f"{'='*60}")

            for test in TEST_PROMPTS:
                prompt_id = test["id"]
                prompt_text = test["prompt"]
                expected_tools = expected[prompt_id]

                print(f"  {prompt_id}: ", end="", flush=True)

                for rep in range(N_REPS):
                    trial_key = f"{condition}|{model_name}|{prompt_id}|{rep}"

                    if trial_key in completed:
                        skipped += 1
                        print("s", end="", flush=True)
                        continue

                    count += 1
                    raw_text = caller(context, prompt_text)
                    parsed = parse_json(raw_text) if raw_text else None

                    raw_tool = parsed.get("tool_name") if parsed else None
                    canonical_tool = resolver(raw_tool)
                    params = parsed.get("parameters", {}) if parsed else {}
                    tool_correct = (
                        canonical_tool in expected_tools
                        if canonical_tool
                        else False
                    )

                    # Also check if tool is valid (in the system's catalog)
                    if condition == "stagent_context":
                        tool_valid = canonical_tool in STAGENT_TOOL_NAMES
                    else:
                        tool_valid = canonical_tool in SPATIALAGENT_TOOL_NAMES
                    tool_valid = tool_valid if canonical_tool else False

                    record = {
                        "trial_key": trial_key,
                        "condition": condition,
                        "model": model_name,
                        "prompt_id": prompt_id,
                        "rep": rep,
                        "raw_response": (raw_text or "")[:2000],
                        "parse_success": parsed is not None,
                        "raw_tool_name": raw_tool,
                        "canonical_tool": canonical_tool,
                        "tool_valid": tool_valid,
                        "tool_correct": tool_correct,
                        "parameters": json.dumps(params) if params else "{}",
                    }
                    append_raw(RAW_PATH, record)
                    completed.add(trial_key)

                    if not parsed:
                        status = "x"
                    elif tool_correct:
                        status = "."
                    elif tool_valid:
                        status = "w"  # wrong but valid tool
                    else:
                        status = "?"  # unrecognized tool
                    print(status, end="", flush=True)
                    time.sleep(DELAY)

                print()

    print(f"\nDone: {count} new trials, {skipped} skipped (checkpoint)")
    print("Generating summary ...")
    generate_summary()


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def load_jsonl(path: Path) -> list[dict]:
    records = []
    if not path.exists():
        return records
    with open(path) as f:
        for line in f:
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return records


def generate_summary():
    """Aggregate results across all 5 conditions for comparison."""
    # Load new cross-system data
    new_records = load_jsonl(RAW_PATH)

    # Load existing ablation data
    abl_records = load_jsonl(ABLATION_RAW)

    # Combine (ablation uses slightly different field names)
    all_records = []
    for rec in abl_records:
        all_records.append({
            "condition": rec["condition"],
            "model": rec["model"],
            "prompt_id": rec["prompt_id"],
            "rep": rec["rep"],
            "parse_success": rec.get("parse_success", False),
            "canonical_tool": rec.get("canonical_tool"),
            "tool_correct": rec.get("tool_correct", False),
        })
    for rec in new_records:
        all_records.append({
            "condition": rec["condition"],
            "model": rec["model"],
            "prompt_id": rec["prompt_id"],
            "rep": rec["rep"],
            "parse_success": rec.get("parse_success", False),
            "canonical_tool": rec.get("canonical_tool"),
            "tool_correct": rec.get("tool_correct", False),
        })

    # Group by (condition, model, prompt_id)
    groups = defaultdict(list)
    for rec in all_records:
        key = (rec["condition"], rec["model"], rec["prompt_id"])
        groups[key].append(rec)

    # Compute per-group metrics
    rows = []
    for (condition, model, prompt_id), recs in sorted(groups.items()):
        n = len(recs)
        parse_rate = sum(1 for r in recs if r["parse_success"]) / n if n else 0
        correct_rate = sum(1 for r in recs if r["tool_correct"]) / n if n else 0

        # Tool consistency: mode agreement across reps
        tools = [r["canonical_tool"] for r in recs if r["parse_success"] and r["canonical_tool"]]
        tool_consistency = compute_consistency(tools) if tools else 0

        rows.append({
            "condition": condition,
            "model": model,
            "prompt_id": prompt_id,
            "n_trials": n,
            "parse_rate": round(parse_rate, 3),
            "tool_correct_rate": round(correct_rate, 3),
            "tool_consistency": round(tool_consistency, 3),
        })

    # Write CSV
    if rows:
        with open(CSV_PATH, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"CSV saved: {CSV_PATH}")

    # Write summary
    all_conditions = [
        "full_schema", "bare_schema", "no_schema",
        "stagent_context", "spatialagent_context",
    ]
    condition_labels = {
        "full_schema": "ChatSpatial (full)",
        "bare_schema": "ChatSpatial (bare)",
        "no_schema": "ChatSpatial (none)",
        "stagent_context": "STAgent",
        "spatialagent_context": "SpatialAgent",
    }

    lines = []
    lines.append("Cross-System Invocation Comparison — Summary")
    lines.append("=" * 70)
    lines.append("")
    lines.append(
        f"{'Condition':<25} {'Parse%':>8} {'ToolCorr%':>10} "
        f"{'Consist%':>10} {'N':>6}"
    )
    lines.append("-" * 70)

    for cond in all_conditions:
        cond_rows = [r for r in rows if r["condition"] == cond]
        if not cond_rows:
            continue
        avg_parse = sum(r["parse_rate"] for r in cond_rows) / len(cond_rows)
        avg_correct = sum(r["tool_correct_rate"] for r in cond_rows) / len(cond_rows)
        avg_consist = sum(r["tool_consistency"] for r in cond_rows) / len(cond_rows)
        total_n = sum(r["n_trials"] for r in cond_rows)
        label = condition_labels.get(cond, cond)
        lines.append(
            f"{label:<25} {avg_parse*100:7.1f}% {avg_correct*100:9.1f}% "
            f"{avg_consist*100:9.1f}% {total_n:>6}"
        )

    lines.append("")
    lines.append("Per-model breakdown:")
    lines.append("-" * 70)
    for cond in all_conditions:
        label = condition_labels.get(cond, cond)
        for model in MODELS:
            model_rows = [
                r for r in rows
                if r["condition"] == cond and r["model"] == model
            ]
            if not model_rows:
                continue
            avg_correct = (
                sum(r["tool_correct_rate"] for r in model_rows) / len(model_rows)
            )
            avg_consist = (
                sum(r["tool_consistency"] for r in model_rows) / len(model_rows)
            )
            lines.append(
                f"  {label:<23} {model:<30} "
                f"correct={avg_correct*100:.1f}% consist={avg_consist*100:.1f}%"
            )

    lines.append("")
    lines.append("Per-prompt breakdown (new conditions only):")
    lines.append("-" * 70)
    for cond in ["stagent_context", "spatialagent_context"]:
        label = condition_labels.get(cond, cond)
        lines.append(f"\n  {label}:")
        for test in TEST_PROMPTS:
            pid = test["id"]
            prompt_rows = [
                r for r in rows
                if r["condition"] == cond and r["prompt_id"] == pid
            ]
            if not prompt_rows:
                continue
            avg_correct = (
                sum(r["tool_correct_rate"] for r in prompt_rows) / len(prompt_rows)
            )
            lines.append(f"    {pid:<20} correct={avg_correct*100:.1f}%")

    lines.append("")

    summary_text = "\n".join(lines)
    with open(SUMMARY_PATH, "w") as f:
        f.write(summary_text)
    print(f"Summary saved: {SUMMARY_PATH}")
    print()
    print(summary_text)


if __name__ == "__main__":
    run_experiment()
