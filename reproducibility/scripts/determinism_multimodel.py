#!/usr/bin/env python3
"""
Multi-Model Determinism Experiment for ChatSpatial Paper

Quantifies how consistently different LLMs select the same MCP tool and parameters
when given the same user prompt with schema constraints.

Models tested: Gemini 2.5 Flash, Claude Haiku 4.5, GPT-5 Mini
All at temperature=1.0 (worst case for determinism).
"""

import csv
import json
import os
import re
import time
from collections import Counter
from pathlib import Path

import requests

# ============================================================
# Configuration
# ============================================================

N_TRIALS = 10  # per model per prompt
TEMPERATURE = 1.0
DELAY = 1.0  # seconds between calls

# API Keys
GEMINI_KEY = os.environ.get("GEMINI_API_KEY", "")
ANTHROPIC_KEY = os.environ.get("ANTHROPIC_API_KEY", "")
OPENAI_KEY = os.environ.get("OPENAI_API_KEY", "")

MODELS = {}
if GEMINI_KEY:
    MODELS["gemini-2.5-flash"] = "gemini"
if ANTHROPIC_KEY:
    MODELS["claude-haiku-4-5-20251001"] = "anthropic"
if OPENAI_KEY:
    MODELS["gpt-5-mini"] = "openai"

# Output
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR.parent / "data"
CSV_PATH = DATA_DIR / "determinism_multimodel.csv"
SUMMARY_PATH = DATA_DIR / "determinism_multimodel_summary.txt"

# ============================================================
# MCP Tool Schema
# ============================================================

TOOL_SCHEMAS = """
You are an AI assistant that selects the appropriate analysis tool and parameters for spatial transcriptomics workflows. You have access to the following MCP tools:

## Available Tools

### 1. identify_spatial_domains
Identify spatial domains and tissue architecture.
Parameters:
- method: one of ["spagcn", "leiden", "louvain", "stagate", "graphst", "banksy"] (default: "spagcn")
- n_domains: integer 1-50 (default: 7)
- resolution: float (default: 0.5)

### 2. deconvolve_data
Deconvolve spatial spots to estimate cell type proportions.
Parameters:
- method: one of ["flashdeconv", "cell2location", "rctd", "destvi", "stereoscope", "spotlight", "tangram", "card"] (default: "flashdeconv")
- reference_data_id: string (REQUIRED)
- cell_type_key: string (REQUIRED)
- use_gpu: boolean (default: false)

### 3. find_spatial_genes
Identify spatially variable genes.
Parameters:
- method: one of ["spatialde", "sparkx", "flashs"] (default: "sparkx")
- n_top_genes: integer or null (default: null = all significant)
- filter_mt_genes: boolean (default: true)

### 4. analyze_cell_communication
Analyze cell-cell communication patterns.
Parameters:
- method: one of ["liana", "cellphonedb", "cellchat_r", "fastccc"] (default: "fastccc")
- species: one of ["human", "mouse", "zebrafish"] (REQUIRED)
- cell_type_key: string (REQUIRED)

### 5. analyze_spatial_statistics
Analyze spatial statistics and autocorrelation.
Parameters:
- analysis_type: one of ["neighborhood", "co_occurrence", "ripley", "moran", "local_moran", "geary", "centrality", "getis_ord"] (default: "neighborhood")
- genes: list of strings or null
- n_top_genes: integer 1-500 (default: 20)

### 6. visualize_data
Visualize spatial transcriptomics data.
Parameters:
- plot_type: one of ["feature", "expression", "deconvolution", "communication", "trajectory", "velocity", "statistics", "enrichment", "cnv"] (default: "feature")
- feature: string or list of strings or null
- basis: one of ["spatial", "umap", "pca"] (default: "spatial")

### 7. analyze_trajectory_data
Infer cellular trajectories and pseudotime.
Parameters:
- method: one of ["cellrank", "palantir", "dpt"] (default: "cellrank")
- spatial_weight: float 0-1 (default: 0.5)

### 8. analyze_cnv
Analyze copy number variations.
Parameters:
- method: one of ["infercnvpy", "numbat"] (default: "infercnvpy")
- reference_key: string (REQUIRED)
- reference_categories: list of strings (REQUIRED)
"""

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

# ============================================================
# Test Prompts
# ============================================================

TEST_PROMPTS = [
    {
        "id": "spatial_domains",
        "prompt": "Identify spatial domains in this Visium dataset using deep learning",
        "constrained_params": ["method"],
        "free_params": ["n_domains", "resolution"],
    },
    {
        "id": "deconvolution_rctd",
        "prompt": "Deconvolve cell types using RCTD with the reference I loaded (ref_id='sc_ref', cell_type_key='celltype')",
        "constrained_params": ["method"],
        "free_params": ["reference_data_id", "cell_type_key"],
    },
    {
        "id": "svg_detection",
        "prompt": "Find spatially variable genes in this MERFISH dataset",
        "constrained_params": ["method"],
        "free_params": ["n_top_genes"],
    },
    {
        "id": "cellchat",
        "prompt": "Run CellChat cell communication analysis on this human cancer data (cell_type_key='cell_type')",
        "constrained_params": ["method", "species"],
        "free_params": ["cell_type_key"],
    },
    {
        "id": "moran_i",
        "prompt": "Compute Moran's I spatial autocorrelation for the top marker genes",
        "constrained_params": ["analysis_type"],
        "free_params": ["n_top_genes"],
    },
    {
        "id": "spatial_plot",
        "prompt": "Create a spatial plot colored by cell type annotations in 'cell_type' column",
        "constrained_params": ["plot_type", "basis"],
        "free_params": ["feature"],
    },
    {
        "id": "trajectory",
        "prompt": "Perform trajectory analysis using CellRank",
        "constrained_params": ["method"],
        "free_params": ["spatial_weight"],
    },
    {
        "id": "cnv",
        "prompt": "Infer CNV in this tumor, using immune cells (T cells, B cells) as reference (cell_type column='cell_type')",
        "constrained_params": ["method"],
        "free_params": ["reference_key", "reference_categories"],
    },
]


# ============================================================
# API Callers
# ============================================================

def call_gemini(prompt: str) -> dict | None:
    payload = {
        "contents": [{"role": "user", "parts": [{"text": TOOL_SCHEMAS + RESPONSE_INSTRUCTION + "\nUser: " + prompt}]}],
        "generationConfig": {"temperature": TEMPERATURE, "maxOutputTokens": 512},
    }
    try:
        r = requests.post(
            f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={GEMINI_KEY}",
            # Gemini 2.5 Flash
            json=payload, timeout=30
        )
        r.raise_for_status()
        text = r.json()["candidates"][0]["content"]["parts"][0]["text"]
        return parse_json(text)
    except Exception as e:
        print(f"  Gemini error: {e}")
        return None


def call_anthropic(prompt: str) -> dict | None:
    payload = {
        "model": "claude-haiku-4-5-20251001",
        "max_tokens": 512,
        "temperature": TEMPERATURE,
        "system": TOOL_SCHEMAS + RESPONSE_INSTRUCTION,
        "messages": [{"role": "user", "content": prompt}],
    }
    try:
        r = requests.post(
            "https://api.anthropic.com/v1/messages",
            headers={
                "x-api-key": ANTHROPIC_KEY,
                "anthropic-version": "2023-06-01",
                "content-type": "application/json",
            },
            json=payload, timeout=30
        )
        r.raise_for_status()
        text = r.json()["content"][0]["text"]
        return parse_json(text)
    except Exception as e:
        print(f"  Anthropic error: {e}")
        return None


def call_openai(prompt: str) -> dict | None:
    payload = {
        "model": "gpt-5-mini",
        "max_completion_tokens": 512,
        "temperature": TEMPERATURE,
        "messages": [
            {"role": "system", "content": TOOL_SCHEMAS + RESPONSE_INSTRUCTION},
            {"role": "user", "content": prompt},
        ],
    }
    try:
        r = requests.post(
            "https://api.openai.com/v1/chat/completions",
            headers={"Authorization": f"Bearer {OPENAI_KEY}", "Content-Type": "application/json"},
            json=payload, timeout=30
        )
        r.raise_for_status()
        text = r.json()["choices"][0]["message"]["content"]
        return parse_json(text)
    except Exception as e:
        print(f"  OpenAI error: {e}")
        return None


def parse_json(text: str) -> dict | None:
    text = text.strip()
    # Remove markdown code blocks
    text = re.sub(r"^```(?:json)?\s*", "", text)
    text = re.sub(r"\s*```$", "", text)
    text = text.strip()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        # Try to find JSON in the text
        match = re.search(r'\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}', text)
        if match:
            try:
                return json.loads(match.group())
            except json.JSONDecodeError:
                pass
    return None


CALLERS = {
    "gemini": call_gemini,
    "anthropic": call_anthropic,
    "openai": call_openai,
}


# ============================================================
# Analysis
# ============================================================

def compute_consistency(responses: list[dict | None], key_path: str) -> float:
    """Compute % of valid responses that agree on a value."""
    values = []
    for r in responses:
        if r is None:
            continue
        if key_path == "tool_name":
            v = r.get("tool_name")
        else:
            v = r.get("parameters", {}).get(key_path)
        if v is not None:
            # Normalize: convert to string for comparison
            values.append(json.dumps(v, sort_keys=True) if isinstance(v, (list, dict)) else str(v))

    if not values:
        return 0.0

    most_common_count = Counter(values).most_common(1)[0][1]
    return most_common_count / len(values)


def run_experiment():
    print(f"Running determinism experiment with {len(MODELS)} models, {N_TRIALS} trials each")
    print(f"Models: {list(MODELS.keys())}")
    print(f"Temperature: {TEMPERATURE}")
    print()

    all_results = []

    for model_name, provider in MODELS.items():
        caller = CALLERS[provider]
        print(f"\n{'='*60}")
        print(f"Model: {model_name}")
        print(f"{'='*60}")

        for test in TEST_PROMPTS:
            prompt_id = test["id"]
            prompt = test["prompt"]
            print(f"\n  Prompt: {prompt_id} ({N_TRIALS} trials)...", end="", flush=True)

            responses = []
            for i in range(N_TRIALS):
                r = caller(prompt)
                responses.append(r)
                print(".", end="", flush=True)
                time.sleep(DELAY)

            valid = [r for r in responses if r is not None]
            parse_rate = len(valid) / N_TRIALS

            # Tool selection consistency
            tool_consistency = compute_consistency(responses, "tool_name")

            # Constrained parameter consistency
            constrained_scores = []
            for p in test["constrained_params"]:
                constrained_scores.append(compute_consistency(responses, p))
            avg_constrained = sum(constrained_scores) / len(constrained_scores) if constrained_scores else 0

            # Free parameter consistency
            free_scores = []
            for p in test["free_params"]:
                free_scores.append(compute_consistency(responses, p))
            avg_free = sum(free_scores) / len(free_scores) if free_scores else 0

            # Overall determinism
            overall = (tool_consistency + avg_constrained) / 2

            result = {
                "model": model_name,
                "prompt_id": prompt_id,
                "n_trials": N_TRIALS,
                "valid_responses": len(valid),
                "parse_rate": round(parse_rate, 3),
                "tool_consistency": round(tool_consistency, 3),
                "constrained_param_consistency": round(avg_constrained, 3),
                "free_param_consistency": round(avg_free, 3),
                "overall_determinism": round(overall, 3),
            }
            all_results.append(result)

            print(f" tool={tool_consistency:.0%} constrained={avg_constrained:.0%} free={avg_free:.0%}")

    # Write CSV
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_results[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nCSV saved: {CSV_PATH}")

    # Write summary
    with open(SUMMARY_PATH, "w") as f:
        f.write("Determinism Experiment Summary\n")
        f.write(f"{'='*60}\n")
        f.write(f"Trials per prompt per model: {N_TRIALS}\n")
        f.write(f"Temperature: {TEMPERATURE}\n")
        f.write(f"Models: {list(MODELS.keys())}\n\n")

        for model_name in MODELS:
            model_results = [r for r in all_results if r["model"] == model_name]
            avg_tool = sum(r["tool_consistency"] for r in model_results) / len(model_results)
            avg_constrained = sum(r["constrained_param_consistency"] for r in model_results) / len(model_results)
            avg_free = sum(r["free_param_consistency"] for r in model_results) / len(model_results)
            avg_overall = sum(r["overall_determinism"] for r in model_results) / len(model_results)
            avg_parse = sum(r["parse_rate"] for r in model_results) / len(model_results)

            f.write(f"\n{model_name}\n")
            f.write(f"{'-'*40}\n")
            f.write(f"  Parse rate:                    {avg_parse:.1%}\n")
            f.write(f"  Tool selection consistency:     {avg_tool:.1%}\n")
            f.write(f"  Constrained param consistency:  {avg_constrained:.1%}\n")
            f.write(f"  Free-text param consistency:    {avg_free:.1%}\n")
            f.write(f"  Overall determinism score:      {avg_overall:.1%}\n")

        # Cross-model averages
        f.write(f"\n\nCross-Model Averages\n")
        f.write(f"{'='*40}\n")
        all_tool = sum(r["tool_consistency"] for r in all_results) / len(all_results)
        all_constrained = sum(r["constrained_param_consistency"] for r in all_results) / len(all_results)
        all_free = sum(r["free_param_consistency"] for r in all_results) / len(all_results)
        all_overall = sum(r["overall_determinism"] for r in all_results) / len(all_results)
        f.write(f"  Tool selection:     {all_tool:.1%}\n")
        f.write(f"  Constrained params: {all_constrained:.1%}\n")
        f.write(f"  Free-text params:   {all_free:.1%}\n")
        f.write(f"  Overall:            {all_overall:.1%}\n")

        # Key finding for paper
        f.write(f"\n\nKey Finding for Paper\n")
        f.write(f"{'='*40}\n")
        f.write(f"Schema-constrained parameters showed {all_constrained:.1%} consistency\n")
        f.write(f"across {len(MODELS)} models at temperature={TEMPERATURE} (worst case),\n")
        f.write(f"compared to {all_free:.1%} for free-text parameters.\n")
        f.write(f"Tool selection was {all_tool:.1%} consistent.\n")

    print(f"Summary saved: {SUMMARY_PATH}")


if __name__ == "__main__":
    run_experiment()
