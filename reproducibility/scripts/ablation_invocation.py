#!/usr/bin/env python3
"""
Schema-Enforcement Ablation Study — Part 1: Invocation-Level Experiment

Tests whether schema enforcement causally improves LLM invocation consistency
by comparing three conditions:
  A) Full schema (types + enums + descriptions + method guidance)
  B) Bare schema (types + enums + bounds, NO descriptions)
  C) Minimal   (tool names + 1-sentence descriptions only, NO typed parameters)

Trial matrix: 8 prompts x 3 models x 10 reps x 3 conditions = 720 API calls.

Extends the existing determinism_multimodel.py experiment with:
  - Three schema conditions (the ablation variable)
  - Pydantic model validation (executability measurement)
  - Incremental JSONL checkpointing
"""

import csv
import json
import os
import re
import time
from collections import Counter
from pathlib import Path

import requests
from paths import load_env_file

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

load_env_file()

N_TRIALS = 10
TEMPERATURE = 1.0
DELAY = 1.0  # seconds between calls
MAX_RETRIES = 3
RETRY_BACKOFF = 2.0

GEMINI_KEY = os.environ.get("GEMINI_API_KEY", "")
ANTHROPIC_KEY = os.environ.get("ANTHROPIC_API_KEY", "")
OPENAI_KEY = os.environ.get("OPENAI_API_KEY", "")

# GPT-5.4 was run through an Anthropic Messages-compatible gateway.
OPENAI_API_URL = os.environ.get("OPENAI_API_URL", "")
OPENAI_MODEL = "gpt-5.4"

MODELS: dict[str, str] = {}
if GEMINI_KEY:
    MODELS["gemini-2.5-flash"] = "gemini"
if ANTHROPIC_KEY:
    MODELS["claude-haiku-4-5-20251001"] = "anthropic"
if OPENAI_KEY and OPENAI_API_URL:
    MODELS[OPENAI_MODEL] = "openai"

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR.parent / "data" / "ablation" / "invocation"
RAW_PATH = DATA_DIR / "ablation_invocation_raw.jsonl"
CSV_PATH = DATA_DIR / "ablation_invocation_results.csv"
SUMMARY_PATH = DATA_DIR / "ablation_invocation_summary.txt"

# ---------------------------------------------------------------------------
# Schema Variants (the ablation variable)
# ---------------------------------------------------------------------------

# Condition A: Full Schema — types + enums + descriptions + method guidance
# (same as the proven schema from determinism_multimodel.py)
SCHEMA_FULL = """
You are an AI assistant that selects the appropriate analysis tool and parameters for spatial transcriptomics workflows. You have access to the following MCP tools:

## Available Tools

### 1. identify_spatial_domains
Identify spatial domains and tissue architecture.
Parameters:
- method: one of ["spagcn", "leiden", "louvain", "stagate", "graphst", "banksy"] (default: "spagcn"). 'spagcn' uses histology images. 'banksy' uses spatial feature augmentation. 'stagate'/'graphst' use deep learning for spatial embedding.
- n_domains: integer 1-50 (default: 7). Number of spatial domains to identify.
- resolution: float (default: 0.5). Leiden/Louvain clustering resolution. Higher values produce more clusters.

### 2. deconvolve_data
Deconvolve spatial spots to estimate cell type proportions.
Parameters:
- method: one of ["flashdeconv", "cell2location", "rctd", "destvi", "stereoscope", "spotlight", "tangram", "card"] (default: "flashdeconv"). All methods require reference_data_id and cell_type_key.
- reference_data_id: string (REQUIRED). ID of loaded single-cell reference dataset.
- cell_type_key: string (REQUIRED). Column in reference data with cell type labels.
- use_gpu: boolean (default: false). GPU acceleration. Supported by cell2location, destvi, stereoscope, tangram.

### 3. find_spatial_genes
Identify spatially variable genes.
Parameters:
- method: one of ["spatialde", "sparkx", "flashs"] (default: "sparkx"). 'sparkx' requires R. 'flashs' is Python-native and ultra-fast. 'spatialde' uses Gaussian process kernels.
- n_top_genes: integer or null (default: null = all significant). Number of top genes to return.
- filter_mt_genes: boolean (default: true). Filter mitochondrial genes from results.

### 4. analyze_cell_communication
Analyze cell-cell communication patterns.
Parameters:
- method: one of ["liana", "cellphonedb", "cellchat_r", "fastccc"] (default: "fastccc"). 'cellchat_r' calls R-native CellChat. 'liana' supports multiple backends. 'fastccc' is Python-native and fast.
- species: one of ["human", "mouse", "zebrafish"] (REQUIRED). Species for ligand-receptor database.
- cell_type_key: string (REQUIRED). Column with cell type annotations.

### 5. analyze_spatial_statistics
Analyze spatial statistics and autocorrelation.
Parameters:
- analysis_type: one of ["neighborhood", "co_occurrence", "ripley", "moran", "local_moran", "geary", "centrality", "getis_ord"] (default: "neighborhood"). Type of spatial statistic to compute.
- genes: list of strings or null. Specific genes to analyze.
- n_top_genes: integer 1-500 (default: 20). Number of top genes for gene-level statistics.

### 6. visualize_data
Visualize spatial transcriptomics data.
Parameters:
- plot_type: one of ["feature", "expression", "deconvolution", "communication", "trajectory", "velocity", "statistics", "enrichment", "cnv"] (default: "feature"). Type of visualization.
- feature: string or list of strings or null. Feature(s) to visualize.
- basis: one of ["spatial", "umap", "pca"] (default: "spatial"). Coordinate space for plotting.

### 7. analyze_trajectory_data
Infer cellular trajectories and pseudotime.
Parameters:
- method: one of ["cellrank", "palantir", "dpt"] (default: "cellrank"). 'cellrank' uses velocity + expression. 'palantir' uses diffusion maps. 'dpt' is simpler diffusion pseudotime.
- spatial_weight: float 0-1 (default: 0.5). Weight of spatial coordinates in trajectory computation.

### 8. analyze_cnv
Analyze copy number variations.
Parameters:
- method: one of ["infercnvpy", "numbat"] (default: "infercnvpy"). 'infercnvpy' is Python-native. 'numbat' requires R.
- reference_key: string (REQUIRED). Column in adata.obs for identifying reference cells.
- reference_categories: list of strings (REQUIRED). Cell types to use as normal reference.
"""

# Condition B: Bare Schema — types + enums + bounds, NO descriptions
SCHEMA_BARE = """
You are an AI assistant that selects the appropriate analysis tool and parameters for spatial transcriptomics workflows. You have access to the following MCP tools:

## Available Tools

### 1. identify_spatial_domains
Parameters:
- method: one of ["spagcn", "leiden", "louvain", "stagate", "graphst", "banksy"] (default: "spagcn")
- n_domains: integer 1-50 (default: 7)
- resolution: float (default: 0.5)

### 2. deconvolve_data
Parameters:
- method: one of ["flashdeconv", "cell2location", "rctd", "destvi", "stereoscope", "spotlight", "tangram", "card"] (default: "flashdeconv")
- reference_data_id: string (REQUIRED)
- cell_type_key: string (REQUIRED)
- use_gpu: boolean (default: false)

### 3. find_spatial_genes
Parameters:
- method: one of ["spatialde", "sparkx", "flashs"] (default: "sparkx")
- n_top_genes: integer or null (default: null)
- filter_mt_genes: boolean (default: true)

### 4. analyze_cell_communication
Parameters:
- method: one of ["liana", "cellphonedb", "cellchat_r", "fastccc"] (default: "fastccc")
- species: one of ["human", "mouse", "zebrafish"] (REQUIRED)
- cell_type_key: string (REQUIRED)

### 5. analyze_spatial_statistics
Parameters:
- analysis_type: one of ["neighborhood", "co_occurrence", "ripley", "moran", "local_moran", "geary", "centrality", "getis_ord"] (default: "neighborhood")
- genes: list of strings or null
- n_top_genes: integer 1-500 (default: 20)

### 6. visualize_data
Parameters:
- plot_type: one of ["feature", "expression", "deconvolution", "communication", "trajectory", "velocity", "statistics", "enrichment", "cnv"] (default: "feature")
- feature: string or list of strings or null
- basis: one of ["spatial", "umap", "pca"] (default: "spatial")

### 7. analyze_trajectory_data
Parameters:
- method: one of ["cellrank", "palantir", "dpt"] (default: "cellrank")
- spatial_weight: float 0-1 (default: 0.5)

### 8. analyze_cnv
Parameters:
- method: one of ["infercnvpy", "numbat"] (default: "infercnvpy")
- reference_key: string (REQUIRED)
- reference_categories: list of strings (REQUIRED)
"""

# Condition C: Minimal — tool names + 1-sentence descriptions, NO typed parameters
SCHEMA_MINIMAL = """
You are an AI assistant that selects the appropriate analysis tool and parameters for spatial transcriptomics workflows. You have access to the following tools:

1. identify_spatial_domains — Identify spatial domains and tissue architecture in spatial transcriptomics data.
2. deconvolve_data — Estimate cell type proportions at each spatial location using a reference dataset.
3. find_spatial_genes — Identify genes with statistically significant spatial expression patterns.
4. analyze_cell_communication — Analyze cell-cell communication and ligand-receptor interactions.
5. analyze_spatial_statistics — Compute spatial statistics and autocorrelation measures.
6. visualize_data — Visualize spatial transcriptomics data and analysis results.
7. analyze_trajectory_data — Infer cellular trajectories and pseudotime ordering.
8. analyze_cnv — Analyze copy number variations from spatial transcriptomics data.
"""

SCHEMA_VARIANTS = {
    "full_schema": SCHEMA_FULL,
    "bare_schema": SCHEMA_BARE,
    "no_schema": SCHEMA_MINIMAL,
}

RESPONSE_INSTRUCTION = """
Based on the user's request, select the most appropriate tool and parameters.
Respond with ONLY a JSON object (no markdown, no explanation):
{
  "tool_name": "<tool_name>",
  "parameters": {
    "<param1>": <value1>,
    ...
  }
}
Include only parameters that differ from defaults or are required.
"""

# ---------------------------------------------------------------------------
# Test Prompts (same 8 as determinism_multimodel.py)
# ---------------------------------------------------------------------------

TEST_PROMPTS = [
    {
        "id": "spatial_domains",
        "prompt": "Identify spatial domains in this Visium dataset using deep learning",
        "expected_tool": "identify_spatial_domains",
        "constrained_params": ["method"],
        "free_params": ["n_domains", "resolution"],
    },
    {
        "id": "deconvolution_rctd",
        "prompt": "Deconvolve cell types using RCTD with the reference I loaded (ref_id='sc_ref', cell_type_key='celltype')",
        "expected_tool": "deconvolve_data",
        "constrained_params": ["method"],
        "free_params": ["reference_data_id", "cell_type_key"],
    },
    {
        "id": "svg_detection",
        "prompt": "Find spatially variable genes in this MERFISH dataset",
        "expected_tool": "find_spatial_genes",
        "constrained_params": ["method"],
        "free_params": ["n_top_genes"],
    },
    {
        "id": "cellchat",
        "prompt": "Run CellChat cell communication analysis on this human cancer data (cell_type_key='cell_type')",
        "expected_tool": "analyze_cell_communication",
        "constrained_params": ["method", "species"],
        "free_params": ["cell_type_key"],
    },
    {
        "id": "moran_i",
        "prompt": "Compute Moran's I spatial autocorrelation for the top marker genes",
        "expected_tool": "analyze_spatial_statistics",
        "constrained_params": ["analysis_type"],
        "free_params": ["n_top_genes"],
    },
    {
        "id": "spatial_plot",
        "prompt": "Create a spatial plot colored by cell type annotations in 'cell_type' column",
        "expected_tool": "visualize_data",
        "constrained_params": ["plot_type", "basis"],
        "free_params": ["feature"],
    },
    {
        "id": "trajectory",
        "prompt": "Perform trajectory analysis using CellRank",
        "expected_tool": "analyze_trajectory_data",
        "constrained_params": ["method"],
        "free_params": ["spatial_weight"],
    },
    {
        "id": "cnv",
        "prompt": "Infer CNV in this tumor, using immune cells (T cells, B cells) as reference (cell_type column='cell_type')",
        "expected_tool": "analyze_cnv",
        "constrained_params": ["method"],
        "free_params": ["reference_key", "reference_categories"],
    },
]

# Tool name aliases for Condition C (generous mapping)
TOOL_ALIASES: dict[str, list[str]] = {
    "identify_spatial_domains": [
        "identify_spatial_domains",
        "spatial_domains",
        "spatial_clustering",
        "cluster_spatial",
        "find_domains",
        "find_spatial_domains",
    ],
    "deconvolve_data": [
        "deconvolve_data",
        "deconvolution",
        "deconvolve",
        "cell_type_deconvolution",
        "estimate_proportions",
        "estimate_cell_type_proportions",
    ],
    "find_spatial_genes": [
        "find_spatial_genes",
        "spatial_genes",
        "svg_detection",
        "find_svg",
        "spatially_variable_genes",
        "identify_spatial_genes",
        "spatial_variable_genes",
        "find_spatially_variable_genes",
    ],
    "analyze_cell_communication": [
        "analyze_cell_communication",
        "cell_communication",
        "cellchat",
        "cell_chat",
        "ligand_receptor",
        "analyze_communication",
        "cell_cell_communication",
    ],
    "analyze_spatial_statistics": [
        "analyze_spatial_statistics",
        "spatial_statistics",
        "spatial_autocorrelation",
        "morans_i",
        "moran_i",
        "compute_spatial_statistics",
    ],
    "visualize_data": [
        "visualize_data",
        "visualize",
        "plot",
        "spatial_plot",
        "create_plot",
        "visualization",
    ],
    "analyze_trajectory_data": [
        "analyze_trajectory_data",
        "trajectory",
        "trajectory_analysis",
        "analyze_trajectory",
        "pseudotime",
    ],
    "analyze_cnv": [
        "analyze_cnv",
        "cnv",
        "cnv_analysis",
        "infer_cnv",
        "copy_number",
        "cnv_inference",
    ],
}

# Build reverse lookup
_ALIAS_TO_CANONICAL: dict[str, str] = {}
for canonical, aliases in TOOL_ALIASES.items():
    for alias in aliases:
        _ALIAS_TO_CANONICAL[alias.lower()] = canonical


def resolve_tool_name(raw_name: str | None) -> str | None:
    """Map a raw tool name to canonical name via aliases. Returns None if unresolvable."""
    if raw_name is None:
        return None
    return _ALIAS_TO_CANONICAL.get(raw_name.lower().strip())


# ---------------------------------------------------------------------------
# Pydantic Validation
# ---------------------------------------------------------------------------

try:
    from chatspatial.models.data import (
        CellCommunicationParameters,
        CNVParameters,
        DeconvolutionParameters,
        SpatialDomainParameters,
        SpatialStatisticsParameters,
        SpatialVariableGenesParameters,
        TrajectoryParameters,
        VisualizationParameters,
    )

    TOOL_PARAM_MODELS: dict[str, type] = {
        "identify_spatial_domains": SpatialDomainParameters,
        "deconvolve_data": DeconvolutionParameters,
        "find_spatial_genes": SpatialVariableGenesParameters,
        "analyze_cell_communication": CellCommunicationParameters,
        "analyze_spatial_statistics": SpatialStatisticsParameters,
        "visualize_data": VisualizationParameters,
        "analyze_trajectory_data": TrajectoryParameters,
        "analyze_cnv": CNVParameters,
    }
    PYDANTIC_AVAILABLE = True
except ImportError:
    TOOL_PARAM_MODELS = {}
    PYDANTIC_AVAILABLE = False
    print("WARNING: chatspatial not importable, Pydantic validation disabled")


def validate_with_pydantic(tool_name: str, params: dict) -> tuple[bool, str]:
    """Attempt Pydantic validation. Returns (success, error_message)."""
    if not PYDANTIC_AVAILABLE:
        return False, "pydantic_unavailable"
    model_cls = TOOL_PARAM_MODELS.get(tool_name)
    if model_cls is None:
        return False, f"unknown_tool:{tool_name}"
    try:
        model_cls(**params)
        return True, ""
    except Exception as e:
        return False, str(e)[:500]


# ---------------------------------------------------------------------------
# API Callers (adapted from determinism_multimodel.py)
# ---------------------------------------------------------------------------


def call_gemini(system: str, prompt: str) -> str | None:
    for attempt in range(MAX_RETRIES):
        try:
            payload = {
                "contents": [
                    {
                        "role": "user",
                        "parts": [
                            {
                                "text": system
                                + RESPONSE_INSTRUCTION
                                + "\nUser: "
                                + prompt
                            }
                        ],
                    }
                ],
                "generationConfig": {
                    "temperature": TEMPERATURE,
                    "maxOutputTokens": 512,
                },
            }
            r = requests.post(
                f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={GEMINI_KEY}",
                json=payload,
                timeout=60,
            )
            r.raise_for_status()
            return r.json()["candidates"][0]["content"]["parts"][0]["text"]
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF * (attempt + 1))
            else:
                print(f"  Gemini error after {MAX_RETRIES} retries: {e}")
                return None


def call_anthropic(system: str, prompt: str) -> str | None:
    for attempt in range(MAX_RETRIES):
        try:
            payload = {
                "model": "claude-haiku-4-5-20251001",
                "max_tokens": 512,
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
                json=payload,
                timeout=60,
            )
            r.raise_for_status()
            return r.json()["content"][0]["text"]
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF * (attempt + 1))
            else:
                print(f"  Anthropic error after {MAX_RETRIES} retries: {e}")
                return None


def call_openai(system: str, prompt: str) -> str | None:
    """Call OpenAI GPT model via Anthropic Messages-compatible gateway."""
    if not OPENAI_KEY or not OPENAI_API_URL:
        print(
            "  GPT-5.4 gateway not configured; set OPENAI_API_KEY and "
            "OPENAI_API_URL to run this model."
        )
        return None
    for attempt in range(MAX_RETRIES):
        try:
            payload = {
                "model": OPENAI_MODEL,
                "max_tokens": 512,
                "temperature": TEMPERATURE,
                "system": system + RESPONSE_INSTRUCTION,
                "messages": [{"role": "user", "content": prompt}],
            }
            r = requests.post(
                OPENAI_API_URL,
                headers={
                    "x-api-key": OPENAI_KEY,
                    "anthropic-version": "2023-06-01",
                    "content-type": "application/json",
                },
                json=payload,
                timeout=120,
            )
            r.raise_for_status()
            # Extract text from response (skip thinking blocks)
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
    "gemini": call_gemini,
    "anthropic": call_anthropic,
    "openai": call_openai,
}


def parse_json(text: str) -> dict | None:
    """Parse JSON from LLM response, handling markdown blocks and thinking tags."""
    text = text.strip()
    # Strip thinking tags (Gemini 2.5)
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL)
    # Strip markdown code blocks
    text = re.sub(r"^```(?:json)?\s*", "", text)
    text = re.sub(r"\s*```$", "", text)
    text = text.strip()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        match = re.search(r"\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}", text)
        if match:
            try:
                return json.loads(match.group())
            except json.JSONDecodeError:
                pass
    return None


# ---------------------------------------------------------------------------
# Checkpoint / Resume
# ---------------------------------------------------------------------------


def load_completed_trials(path: Path) -> set[str]:
    """Load trial IDs already completed from JSONL checkpoint."""
    completed = set()
    if path.exists():
        with open(path) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    completed.add(rec["trial_key"])
                except (json.JSONDecodeError, KeyError):
                    continue
    return completed


def append_raw(path: Path, record: dict) -> None:
    """Append a single record to JSONL checkpoint."""
    with open(path, "a") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")


# ---------------------------------------------------------------------------
# Consistency Metrics
# ---------------------------------------------------------------------------


def compute_consistency(values: list) -> float:
    """Fraction of values agreeing with the mode."""
    if not values:
        return 0.0
    normed = [
        json.dumps(v, sort_keys=True) if isinstance(v, (list, dict)) else str(v)
        for v in values
    ]
    return Counter(normed).most_common(1)[0][1] / len(normed)


# ---------------------------------------------------------------------------
# Main Experiment
# ---------------------------------------------------------------------------


def run_experiment():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    print("Schema-Enforcement Ablation — Part 1: Invocation-Level")
    print(f"Models: {list(MODELS.keys())}")
    print(f"Conditions: {list(SCHEMA_VARIANTS.keys())}")
    print(f"Prompts: {len(TEST_PROMPTS)}, Reps: {N_TRIALS}, Temp: {TEMPERATURE}")
    total = len(MODELS) * len(SCHEMA_VARIANTS) * len(TEST_PROMPTS) * N_TRIALS
    print(f"Total trials: {total}")
    print(f"Pydantic validation: {'enabled' if PYDANTIC_AVAILABLE else 'DISABLED'}")
    print()

    completed = load_completed_trials(RAW_PATH)
    print(f"Resuming: {len(completed)} trials already completed\n")

    for condition, schema_text in SCHEMA_VARIANTS.items():
        for model_name, provider in MODELS.items():
            caller = CALLERS[provider]
            print(f"\n{'='*60}")
            print(f"Condition: {condition} | Model: {model_name}")
            print(f"{'='*60}")

            for test in TEST_PROMPTS:
                prompt_id = test["id"]
                prompt = test["prompt"]
                print(f"  {prompt_id}: ", end="", flush=True)

                responses = []
                for rep in range(N_TRIALS):
                    trial_key = f"{condition}|{model_name}|{prompt_id}|{rep}"

                    if trial_key in completed:
                        print("s", end="", flush=True)
                        continue

                    raw_text = caller(schema_text, prompt)
                    parsed = parse_json(raw_text) if raw_text else None

                    # Extract tool name
                    raw_tool = parsed.get("tool_name") if parsed else None
                    canonical_tool = resolve_tool_name(raw_tool)
                    params = parsed.get("parameters", {}) if parsed else {}

                    # Pydantic validation
                    pydantic_valid = False
                    pydantic_errors = ""
                    if canonical_tool and params is not None:
                        pydantic_valid, pydantic_errors = validate_with_pydantic(
                            canonical_tool, params
                        )

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
                        "tool_correct": canonical_tool == test["expected_tool"],
                        "parameters": json.dumps(params) if params else "{}",
                        "pydantic_valid": pydantic_valid,
                        "pydantic_errors": pydantic_errors[:300],
                    }
                    append_raw(RAW_PATH, record)
                    completed.add(trial_key)
                    responses.append(record)

                    print("." if parsed else "x", end="", flush=True)
                    time.sleep(DELAY)

                print(f" ({len(responses)} new)")

    # ---- Analyze all results from checkpoint file ----
    print("\n\nAnalyzing results from checkpoint...")
    all_records = []
    with open(RAW_PATH) as f:
        for line in f:
            try:
                all_records.append(json.loads(line))
            except json.JSONDecodeError:
                continue

    # Group by (condition, model, prompt_id)
    from collections import defaultdict

    groups = defaultdict(list)
    for rec in all_records:
        key = (rec["condition"], rec["model"], rec["prompt_id"])
        groups[key].append(rec)

    # Compute per-group metrics
    summary_rows = []
    for (condition, model, prompt_id), recs in sorted(groups.items()):
        test = next((t for t in TEST_PROMPTS if t["id"] == prompt_id), None)
        if test is None:
            continue

        n = len(recs)
        valid = [r for r in recs if r["parse_success"]]
        parse_rate = len(valid) / n if n else 0

        # Tool selection consistency (on parsed responses)
        tools = [r["canonical_tool"] for r in valid if r["canonical_tool"]]
        tool_consistency = compute_consistency(tools) if tools else 0

        # Tool correctness
        correct = (
            sum(1 for r in valid if r["tool_correct"]) / len(valid) if valid else 0
        )

        # Pydantic validation rate
        pydantic_rate = sum(1 for r in recs if r["pydantic_valid"]) / n if n else 0

        # Constrained parameter consistency
        constrained_scores = []
        for p in test["constrained_params"]:
            vals = []
            for r in valid:
                params = json.loads(r["parameters"])
                if p in params:
                    vals.append(params[p])
            if vals:
                constrained_scores.append(compute_consistency(vals))
        avg_constrained = (
            sum(constrained_scores) / len(constrained_scores)
            if constrained_scores
            else 0
        )

        # Free parameter consistency
        free_scores = []
        for p in test["free_params"]:
            vals = []
            for r in valid:
                params = json.loads(r["parameters"])
                if p in params:
                    vals.append(params[p])
            if vals:
                free_scores.append(compute_consistency(vals))
        avg_free = sum(free_scores) / len(free_scores) if free_scores else 0

        summary_rows.append(
            {
                "condition": condition,
                "model": model,
                "prompt_id": prompt_id,
                "n_trials": n,
                "parse_rate": round(parse_rate, 3),
                "tool_consistency": round(tool_consistency, 3),
                "tool_correctness": round(correct, 3),
                "pydantic_valid_rate": round(pydantic_rate, 3),
                "constrained_param_consistency": round(avg_constrained, 3),
                "free_param_consistency": round(avg_free, 3),
            }
        )

    # Write CSV
    if summary_rows:
        with open(CSV_PATH, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=summary_rows[0].keys(), lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"CSV saved: {CSV_PATH}")

    # Write summary
    with open(SUMMARY_PATH, "w") as f:
        f.write("Schema-Enforcement Ablation — Part 1: Invocation-Level Summary\n")
        f.write(f"{'='*70}\n")
        f.write(f"Trials per cell: {N_TRIALS} | Temperature: {TEMPERATURE}\n")
        f.write(f"Models: {list(MODELS.keys())}\n")
        f.write(f"Conditions: {list(SCHEMA_VARIANTS.keys())}\n\n")

        for condition in SCHEMA_VARIANTS:
            cond_rows = [r for r in summary_rows if r["condition"] == condition]
            if not cond_rows:
                continue
            avg_parse = sum(r["parse_rate"] for r in cond_rows) / len(cond_rows)
            avg_tool = sum(r["tool_consistency"] for r in cond_rows) / len(cond_rows)
            avg_correct = sum(r["tool_correctness"] for r in cond_rows) / len(cond_rows)
            avg_pydantic = sum(r["pydantic_valid_rate"] for r in cond_rows) / len(
                cond_rows
            )
            avg_const = sum(
                r["constrained_param_consistency"] for r in cond_rows
            ) / len(cond_rows)
            avg_free = sum(r["free_param_consistency"] for r in cond_rows) / len(
                cond_rows
            )

            f.write(f"\n{condition.upper()}\n")
            f.write(f"{'-'*50}\n")
            f.write(f"  Parse rate:                    {avg_parse:.1%}\n")
            f.write(f"  Tool selection consistency:     {avg_tool:.1%}\n")
            f.write(f"  Tool correctness:              {avg_correct:.1%}\n")
            f.write(f"  Pydantic validation rate:      {avg_pydantic:.1%}\n")
            f.write(f"  Constrained param consistency: {avg_const:.1%}\n")
            f.write(f"  Free-text param consistency:   {avg_free:.1%}\n")

            # Per-model breakdown
            for model in MODELS:
                model_rows = [r for r in cond_rows if r["model"] == model]
                if not model_rows:
                    continue
                mp = sum(r["pydantic_valid_rate"] for r in model_rows) / len(model_rows)
                mc = sum(r["constrained_param_consistency"] for r in model_rows) / len(
                    model_rows
                )
                f.write(f"    {model}: pydantic={mp:.1%} constrained={mc:.1%}\n")

        # Cross-condition comparison (the key result)
        f.write(f"\n\n{'='*70}\n")
        f.write("KEY ABLATION RESULTS\n")
        f.write(f"{'='*70}\n\n")
        f.write(f"{'Metric':<35} {'Full':>10} {'Bare':>10} {'Minimal':>10}\n")
        f.write(f"{'-'*65}\n")
        for metric in [
            "parse_rate",
            "tool_consistency",
            "tool_correctness",
            "pydantic_valid_rate",
            "constrained_param_consistency",
        ]:
            vals = {}
            for cond in SCHEMA_VARIANTS:
                rows = [r for r in summary_rows if r["condition"] == cond]
                vals[cond] = sum(r[metric] for r in rows) / len(rows) if rows else 0
            f.write(
                f"{metric:<35} {vals.get('full_schema',0):>9.1%} "
                f"{vals.get('bare_schema',0):>9.1%} "
                f"{vals.get('no_schema',0):>9.1%}\n"
            )

    print(f"Summary saved: {SUMMARY_PATH}")
    print("\nDone!")


if __name__ == "__main__":
    run_experiment()
