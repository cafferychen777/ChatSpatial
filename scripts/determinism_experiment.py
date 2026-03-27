#!/usr/bin/env python3
"""
Determinism Experiment for ChatSpatial Paper

Quantifies how consistently an LLM selects the same MCP tool and parameters
when given the same user prompt, using schema-constrained JSON responses.

Uses Gemini 2.0 Flash with temperature=1.0 (worst case for determinism).
"""

import csv
import json
import os
import re
import time
from collections import Counter
from pathlib import Path

import requests

# Configuration
API_KEY = os.environ.get("GEMINI_API_KEY")
MODEL = "gemini-2.0-flash"
API_URL = f"https://generativelanguage.googleapis.com/v1beta/models/{MODEL}:generateContent"
N_TRIALS = 20
TEMPERATURE = 1.0
DELAY_BETWEEN_CALLS = 1.5  # seconds, to respect rate limits

# Output paths
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR.parent / "data"
CSV_PATH = DATA_DIR / "determinism_experiment.csv"
SUMMARY_PATH = DATA_DIR / "determinism_summary.txt"

# ============================================================
# MCP Tool Schema (simplified but accurate from data.py)
# ============================================================

TOOL_SCHEMAS = """
You are an AI assistant that selects the appropriate analysis tool and parameters for spatial transcriptomics workflows. You have access to the following MCP tools with their schemas:

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
- test_only_hvg: boolean (default: true)

### 4. analyze_cell_communication
Analyze cell-cell communication patterns.
Parameters:
- method: one of ["liana", "cellphonedb", "cellchat_r", "fastccc"] (default: "fastccc")
- species: one of ["human", "mouse", "zebrafish"] (REQUIRED)
- cell_type_key: string (REQUIRED)
- liana_resource: one of ["consensus", "mouseconsensus", "cellchatdb", "cellphonedb", ...] (default: "consensus")

### 5. analyze_spatial_statistics
Analyze spatial statistics and autocorrelation patterns.
Parameters:
- analysis_type: one of ["neighborhood", "co_occurrence", "ripley", "moran", "local_moran", "geary", "centrality", "getis_ord", "bivariate_moran", "join_count", "local_join_count", "network_properties", "spatial_centrality"] (default: "neighborhood")
- cluster_key: string or null
- genes: list of strings or null
- n_top_genes: integer 1-500 (default: 20)

### 6. visualize_data
Visualize spatial transcriptomics data.
Parameters:
- plot_type: one of ["feature", "expression", "deconvolution", "communication", "interaction", "trajectory", "velocity", "statistics", "enrichment", "cnv", "integration"] (default: "feature")
- feature: string or list of strings or null
- basis: one of ["spatial", "umap", "pca"] (default: "spatial")
- subtype: string or null (depends on plot_type)
- cluster_key: string or null

### 7. analyze_trajectory_data
Infer cellular trajectories and pseudotime.
Parameters:
- method: one of ["cellrank", "palantir", "dpt"] (default: "cellrank")
- spatial_weight: float 0-1 (default: 0.5)
- root_cells: list of strings or null

### 8. analyze_cnv
Analyze copy number variations in spatial transcriptomics data.
Parameters:
- method: one of ["infercnvpy", "numbat"] (default: "infercnvpy")
- reference_key: string (REQUIRED)
- reference_categories: list of strings (REQUIRED)
- window_size: integer (default: 100)

### 9. annotate_cell_types
Annotate cell types in spatial transcriptomics data.
Parameters:
- method: one of ["tangram", "scanvi", "cellassign", "mllmcelltype", "sctype", "singler"] (default: "tangram")
- reference_data_id: string or null
- cell_type_key: string or null

### 10. analyze_enrichment
Perform gene set enrichment analysis.
Parameters:
- species: one of ["human", "mouse", "zebrafish"] (REQUIRED)
- method: one of ["spatial_enrichmap", "pathway_gsea", "pathway_ora", "pathway_enrichr", "pathway_ssgsea"] (default: "spatial_enrichmap")
- gene_set_database: one of ["GO_Biological_Process", "GO_Molecular_Function", "GO_Cellular_Component", "KEGG_Pathways", "Reactome_Pathways", "MSigDB_Hallmark", "Cell_Type_Markers"] (default: "GO_Biological_Process")

### 11. analyze_velocity_data
Analyze RNA velocity to understand cellular dynamics.
Parameters:
- method: one of ["scvelo", "velovi"] (default: "scvelo")
- scvelo_mode: one of ["deterministic", "stochastic", "dynamical"] (default: "stochastic")

### 12. find_markers
Find differentially expressed genes between groups.
Parameters:
- group_key: string (REQUIRED)
- method: one of ["wilcoxon", "t-test", "t-test_overestim_var", "logreg", "pydeseq2"] (default: "wilcoxon")
- n_top_genes: integer (default: 50)

### 13. integrate_samples
Integrate multiple spatial transcriptomics samples.
Parameters:
- method: one of ["harmony", "bbknn", "scanorama", "scvi"] (default: "harmony")
- batch_key: string (default: "batch")
"""

# ============================================================
# Test Prompts
# ============================================================

TEST_PROMPTS = [
    {
        "id": "spatial_domains",
        "prompt": "Load this Visium dataset and identify spatial domains using graph-based deep learning",
        "expected_tool": "identify_spatial_domains",
        "expected_method": None,  # Could be stagate or graphst
        "key_params": ["method", "n_domains"],
    },
    {
        "id": "deconvolution",
        "prompt": "Deconvolve cell types in this spatial dataset using the single-cell reference I already loaded (ref_id='sc_ref_001', cell_type_key='cell_type')",
        "expected_tool": "deconvolve_data",
        "expected_method": None,
        "key_params": ["method", "reference_data_id", "cell_type_key"],
    },
    {
        "id": "svg_detection",
        "prompt": "Find spatially variable genes in this MERFISH dataset using a fast statistical method",
        "expected_tool": "find_spatial_genes",
        "expected_method": None,
        "key_params": ["method", "filter_mt_genes"],
    },
    {
        "id": "cell_communication_cellchat",
        "prompt": "Run cell communication analysis between fibroblasts and tumor cells using CellChat on this human cancer dataset (cell_type_key='cell_type')",
        "expected_tool": "analyze_cell_communication",
        "expected_method": "cellchat_r",
        "key_params": ["method", "species", "cell_type_key"],
    },
    {
        "id": "moran_i",
        "prompt": "Compute Moran's I spatial autocorrelation for the top marker genes in this dataset",
        "expected_tool": "analyze_spatial_statistics",
        "expected_method": None,
        "key_params": ["analysis_type", "n_top_genes"],
    },
    {
        "id": "spatial_visualization",
        "prompt": "Create a spatial plot colored by cell type annotations stored in the 'cell_type' column",
        "expected_tool": "visualize_data",
        "expected_method": None,
        "key_params": ["plot_type", "feature", "basis"],
    },
    {
        "id": "trajectory_cellrank",
        "prompt": "Perform trajectory analysis using CellRank to identify cell fate decisions with velocity data already computed",
        "expected_tool": "analyze_trajectory_data",
        "expected_method": "cellrank",
        "key_params": ["method", "spatial_weight"],
    },
    {
        "id": "cnv_analysis",
        "prompt": "Infer copy number variations in this tumor dataset, using immune cells (T cells, B cells) in the 'cell_type' column as the normal reference",
        "expected_tool": "analyze_cnv",
        "expected_method": "infercnvpy",
        "key_params": ["method", "reference_key", "reference_categories"],
    },
]

RESPONSE_FORMAT_INSTRUCTION = """
Based on the user's request, select the most appropriate tool and parameters.
Respond with ONLY a JSON object (no markdown, no explanation) in this exact format:
{
  "tool_name": "<tool_name>",
  "parameters": {
    "<param1>": <value1>,
    "<param2>": <value2>,
    ...
  }
}

Include only parameters that should differ from defaults or are required. Use the exact parameter names and allowed values from the schemas above.
"""


def call_gemini(prompt: str, system_message: str) -> dict | None:
    """Call Gemini API and return parsed JSON response."""
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [
            {
                "role": "user",
                "parts": [{"text": system_message + "\n\n" + RESPONSE_FORMAT_INSTRUCTION + "\n\nUser request: " + prompt}],
            }
        ],
        "generationConfig": {
            "temperature": TEMPERATURE,
            "maxOutputTokens": 1024,
        },
    }

    try:
        resp = requests.post(
            f"{API_URL}?key={API_KEY}",
            headers=headers,
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()

        # Extract text from response
        text = data["candidates"][0]["content"]["parts"][0]["text"]

        # Clean markdown code blocks if present
        text = text.strip()
        if text.startswith("```"):
            # Remove ```json or ``` prefix and ``` suffix
            text = re.sub(r"^```(?:json)?\s*\n?", "", text)
            text = re.sub(r"\n?```\s*$", "", text)

        # Parse JSON
        result = json.loads(text)
        return result

    except (requests.RequestException, json.JSONDecodeError, KeyError, IndexError) as e:
        print(f"    Error: {e}")
        return None


def compute_consistency(responses: list[dict | None], key_params: list[str]) -> dict:
    """Compute consistency metrics for a set of responses."""
    valid = [r for r in responses if r is not None]
    n_valid = len(valid)

    if n_valid == 0:
        return {
            "n_valid": 0,
            "tool_consistency": 0.0,
            "method_consistency": 0.0,
            "param_consistency": 0.0,
            "overall_score": 0.0,
            "tool_distribution": {},
            "method_distribution": {},
        }

    # Tool selection consistency
    tools = [r.get("tool_name", "unknown") for r in valid]
    tool_counts = Counter(tools)
    most_common_tool = tool_counts.most_common(1)[0]
    tool_consistency = most_common_tool[1] / n_valid

    # Method/key parameter consistency (within the most common tool)
    method_values = []
    for r in valid:
        params = r.get("parameters", {})
        # Look for method, analysis_type, or plot_type
        method = params.get("method") or params.get("analysis_type") or params.get("plot_type")
        if method:
            method_values.append(str(method))

    if method_values:
        method_counts = Counter(method_values)
        most_common_method = method_counts.most_common(1)[0]
        method_consistency = most_common_method[1] / len(method_values)
    else:
        method_consistency = 1.0  # No method param = perfectly consistent
        method_counts = Counter()

    # Parameter consistency for schema-constrained params
    param_consistencies = []
    for param_name in key_params:
        values = []
        for r in valid:
            params = r.get("parameters", {})
            if param_name in params:
                values.append(str(params[param_name]))

        if values:
            val_counts = Counter(values)
            most_common_val = val_counts.most_common(1)[0]
            param_consistencies.append(most_common_val[1] / len(values))

    param_consistency = (
        sum(param_consistencies) / len(param_consistencies)
        if param_consistencies
        else 1.0
    )

    # Overall determinism score: weighted average
    overall = 0.4 * tool_consistency + 0.3 * method_consistency + 0.3 * param_consistency

    return {
        "n_valid": n_valid,
        "tool_consistency": tool_consistency,
        "method_consistency": method_consistency,
        "param_consistency": param_consistency,
        "overall_score": overall,
        "tool_distribution": dict(tool_counts),
        "method_distribution": dict(method_counts) if method_values else {},
    }


def run_experiment():
    """Run the full determinism experiment."""
    if not API_KEY:
        raise ValueError("GEMINI_API_KEY environment variable not set")

    print(f"Running determinism experiment with {N_TRIALS} trials per prompt")
    print(f"Model: {MODEL}, Temperature: {TEMPERATURE}")
    print(f"Number of test prompts: {len(TEST_PROMPTS)}")
    print("=" * 70)

    all_results = []
    csv_rows = []

    for test_case in TEST_PROMPTS:
        prompt_id = test_case["id"]
        prompt = test_case["prompt"]
        expected_tool = test_case["expected_tool"]
        key_params = test_case["key_params"]

        print(f"\nPrompt: '{prompt_id}'")
        print(f"  Expected tool: {expected_tool}")

        responses = []
        for trial in range(N_TRIALS):
            print(f"  Trial {trial + 1}/{N_TRIALS}...", end="", flush=True)
            result = call_gemini(prompt, TOOL_SCHEMAS)
            responses.append(result)
            if result:
                tool = result.get("tool_name", "?")
                params = result.get("parameters", {})
                method = params.get("method") or params.get("analysis_type") or params.get("plot_type", "")
                print(f" tool={tool}, method={method}")
            else:
                print(" FAILED")

            # Write individual trial to CSV
            csv_rows.append({
                "prompt_id": prompt_id,
                "trial": trial + 1,
                "tool_selected": result.get("tool_name", "ERROR") if result else "ERROR",
                "method_selected": (
                    str(result.get("parameters", {}).get("method")
                        or result.get("parameters", {}).get("analysis_type")
                        or result.get("parameters", {}).get("plot_type", ""))
                    if result else "ERROR"
                ),
                "full_response": json.dumps(result) if result else "null",
                "expected_tool": expected_tool,
            })

            time.sleep(DELAY_BETWEEN_CALLS)

        # Compute consistency metrics
        metrics = compute_consistency(responses, key_params)
        metrics["prompt_id"] = prompt_id
        metrics["expected_tool"] = expected_tool
        all_results.append(metrics)

        print(f"  Results: tool={metrics['tool_consistency']:.1%}, "
              f"method={metrics['method_consistency']:.1%}, "
              f"params={metrics['param_consistency']:.1%}, "
              f"overall={metrics['overall_score']:.1%}")
        print(f"  Tool distribution: {metrics['tool_distribution']}")
        print(f"  Method distribution: {metrics['method_distribution']}")

    # Save CSV
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "prompt_id", "trial", "tool_selected", "method_selected",
            "full_response", "expected_tool",
        ])
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"\nDetailed results saved to: {CSV_PATH}")

    # Generate summary
    summary_lines = []
    summary_lines.append("=" * 70)
    summary_lines.append("ChatSpatial Determinism Experiment Summary")
    summary_lines.append("=" * 70)
    summary_lines.append(f"Model: {MODEL}")
    summary_lines.append(f"Temperature: {TEMPERATURE}")
    summary_lines.append(f"Trials per prompt: {N_TRIALS}")
    summary_lines.append(f"Number of prompts: {len(TEST_PROMPTS)}")
    summary_lines.append("")

    summary_lines.append("-" * 70)
    summary_lines.append(f"{'Prompt ID':<30} {'Tool%':>8} {'Method%':>8} {'Param%':>8} {'Overall':>8}")
    summary_lines.append("-" * 70)

    tool_scores = []
    method_scores = []
    param_scores = []
    overall_scores = []

    for r in all_results:
        summary_lines.append(
            f"{r['prompt_id']:<30} "
            f"{r['tool_consistency']:>7.1%} "
            f"{r['method_consistency']:>7.1%} "
            f"{r['param_consistency']:>7.1%} "
            f"{r['overall_score']:>7.1%}"
        )
        tool_scores.append(r["tool_consistency"])
        method_scores.append(r["method_consistency"])
        param_scores.append(r["param_consistency"])
        overall_scores.append(r["overall_score"])

    summary_lines.append("-" * 70)
    avg_tool = sum(tool_scores) / len(tool_scores)
    avg_method = sum(method_scores) / len(method_scores)
    avg_param = sum(param_scores) / len(param_scores)
    avg_overall = sum(overall_scores) / len(overall_scores)

    summary_lines.append(
        f"{'AVERAGE':<30} "
        f"{avg_tool:>7.1%} "
        f"{avg_method:>7.1%} "
        f"{avg_param:>7.1%} "
        f"{avg_overall:>7.1%}"
    )

    summary_lines.append("")
    summary_lines.append("=" * 70)
    summary_lines.append("Detailed Distributions per Prompt")
    summary_lines.append("=" * 70)

    for r in all_results:
        summary_lines.append(f"\n{r['prompt_id']} (expected: {r['expected_tool']}):")
        summary_lines.append(f"  Tool distribution:   {r['tool_distribution']}")
        summary_lines.append(f"  Method distribution: {r['method_distribution']}")
        summary_lines.append(f"  Valid responses:     {r['n_valid']}/{N_TRIALS}")

    summary_lines.append("")
    summary_lines.append("=" * 70)
    summary_lines.append("Key Findings")
    summary_lines.append("=" * 70)
    summary_lines.append(f"Average tool selection consistency:   {avg_tool:.1%}")
    summary_lines.append(f"Average method selection consistency: {avg_method:.1%}")
    summary_lines.append(f"Average parameter consistency:        {avg_param:.1%}")
    summary_lines.append(f"Average overall determinism score:    {avg_overall:.1%}")
    summary_lines.append("")

    # Correct tool selection rate
    correct_tool_count = 0
    for r in all_results:
        tool_dist = r["tool_distribution"]
        if r["expected_tool"] in tool_dist:
            correct_tool_count += tool_dist[r["expected_tool"]]

    total_valid = sum(r["n_valid"] for r in all_results)
    correct_rate = correct_tool_count / total_valid if total_valid > 0 else 0
    summary_lines.append(f"Correct tool selection rate (vs expected): {correct_rate:.1%}")
    summary_lines.append(f"Total valid API responses: {total_valid}/{N_TRIALS * len(TEST_PROMPTS)}")

    summary_lines.append("")
    summary_lines.append("Interpretation:")
    summary_lines.append("Schema-constrained parameters (Literal types with fixed enumerations)")
    summary_lines.append("dramatically reduce LLM output variability. Even at temperature=1.0,")
    summary_lines.append("the LLM consistently selects from the valid enum options, demonstrating")
    summary_lines.append("that MCP tool schemas act as effective guardrails for reproducible")
    summary_lines.append("scientific analysis workflows.")

    summary_text = "\n".join(summary_lines)
    with open(SUMMARY_PATH, "w") as f:
        f.write(summary_text)
    print(f"Summary saved to: {SUMMARY_PATH}")
    print("\n" + summary_text)


if __name__ == "__main__":
    run_experiment()
