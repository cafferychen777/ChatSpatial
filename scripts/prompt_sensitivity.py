#!/usr/bin/env python3
"""
Prompt Sensitivity Experiment

Tests whether semantically equivalent paraphrases of user requests produce
consistent tool selection and parameters under the full schema condition.

Trial matrix: 5 prompt groups x 5 paraphrases x 3 models x 5 reps = 375 API calls.

Complements the schema-enforcement ablation: ablation varies the schema,
this experiment varies the user input. Together they bound two axes of variation.
"""

import csv
import json
import time
from collections import Counter
from pathlib import Path

# Reuse infrastructure from ablation experiment
from ablation_invocation import (
    CALLERS,
    DELAY,
    MODELS,
    RESPONSE_INSTRUCTION,
    SCHEMA_FULL,
    append_raw,
    compute_consistency,
    load_completed_trials,
    parse_json,
    resolve_tool_name,
    validate_with_pydantic,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

N_REPS = 5
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR.parent / "data" / "sensitivity"
DATA_DIR.mkdir(parents=True, exist_ok=True)
RAW_PATH = DATA_DIR / "prompt_sensitivity_raw.jsonl"
CSV_PATH = DATA_DIR / "prompt_sensitivity_results.csv"
SUMMARY_PATH = DATA_DIR / "prompt_sensitivity_summary.txt"

# ---------------------------------------------------------------------------
# Prompt Groups: 5 groups x 5 paraphrases each
# ---------------------------------------------------------------------------

PROMPT_GROUPS = [
    {
        "id": "spatial_domains",
        "expected_tool": "identify_spatial_domains",
        "constrained_params": ["method"],
        "free_params": ["n_domains", "resolution"],
        "paraphrases": [
            "Identify spatial domains in this Visium dataset using deep learning",
            "Use a deep-learning approach to find tissue regions in this spatial transcriptomics slide",
            "Segment this Visium section into spatial domains with a neural network method",
            "Apply deep learning to discover spatially coherent tissue domains in my Visium data",
            "Run spatial domain identification on this Visium sample -- I'd like a deep learning method",
        ],
    },
    {
        "id": "deconvolution",
        "expected_tool": "deconvolve_data",
        "constrained_params": ["method"],
        "free_params": ["reference_data_id", "cell_type_key"],
        "paraphrases": [
            "Deconvolve cell types using RCTD with the reference I loaded (ref_id='sc_ref', cell_type_key='celltype')",
            "Estimate cell type proportions via RCTD -- reference dataset is 'sc_ref', cell types in 'celltype' column",
            "Run RCTD deconvolution on my spatial data using reference_data_id='sc_ref' and cell_type_key='celltype'",
            "I want to deconvolve this tissue section with RCTD. The reference is loaded as 'sc_ref', cell type labels are in 'celltype'",
            "Use the RCTD method to infer cell type fractions. Reference: ref_id='sc_ref', cell_type_key='celltype'",
        ],
    },
    {
        "id": "svg_detection",
        "expected_tool": "find_spatial_genes",
        "constrained_params": ["method"],
        "free_params": ["n_top_genes"],
        "paraphrases": [
            "Find spatially variable genes in this MERFISH dataset",
            "Detect genes with significant spatial expression patterns in my MERFISH data",
            "Identify spatially variable genes from this MERFISH experiment",
            "Which genes show spatial variation in this MERFISH dataset? Run the analysis",
            "Perform spatially variable gene detection on this MERFISH slide",
        ],
    },
    {
        "id": "cell_communication",
        "expected_tool": "analyze_cell_communication",
        "constrained_params": ["method", "species"],
        "free_params": ["cell_type_key"],
        "paraphrases": [
            "Run CellChat cell communication analysis on this human cancer data (cell_type_key='cell_type')",
            "Analyze ligand-receptor interactions in this human tumor dataset using CellChat -- cell types are in 'cell_type'",
            "I need cell-cell communication analysis with CellChat on human cancer tissue, cell_type_key='cell_type'",
            "Use CellChat to find signaling pathways in this human cancer spatial dataset (cell type column: 'cell_type')",
            "Run a CellChat-based cell communication analysis -- species is human, cell_type_key='cell_type'",
        ],
    },
    {
        "id": "moran_i",
        "expected_tool": "analyze_spatial_statistics",
        "constrained_params": ["analysis_type"],
        "free_params": ["n_top_genes"],
        "paraphrases": [
            "Compute Moran's I spatial autocorrelation for the top marker genes",
            "Calculate Moran's I statistic for the leading marker genes in this dataset",
            "Run spatial autocorrelation analysis (Moran's I) on the top marker genes",
            "I want Moran's I spatial autocorrelation computed for the top genes",
            "Measure spatial autocorrelation using Moran's I for the most significant marker genes",
        ],
    },
]

# ---------------------------------------------------------------------------
# Main Experiment
# ---------------------------------------------------------------------------

def run_experiment():
    print("Prompt Sensitivity Experiment")
    print(f"Models: {list(MODELS.keys())}")
    print(f"Groups: {[g['id'] for g in PROMPT_GROUPS]}")
    print(f"Paraphrases per group: 5, Reps per paraphrase: {N_REPS}")
    print(f"Total trials: {5 * 5 * len(MODELS) * N_REPS}")
    print()

    completed = load_completed_trials(RAW_PATH)
    print(f"Resuming: {len(completed)} trials already completed")

    system_prompt = SCHEMA_FULL  # Full schema condition only

    total = 0
    skipped = 0

    for group in PROMPT_GROUPS:
        group_id = group["id"]
        expected_tool = group["expected_tool"]

        for para_idx, paraphrase in enumerate(group["paraphrases"]):
            for model_name, provider in MODELS.items():
                caller = CALLERS[provider]

                for rep in range(N_REPS):
                    trial_key = f"{group_id}|p{para_idx}|{model_name}|r{rep}"

                    if trial_key in completed:
                        skipped += 1
                        continue

                    total += 1
                    print(f"  [{total}] {group_id} p{para_idx} {model_name} rep{rep} ... ", end="", flush=True)

                    raw_text = caller(system_prompt, paraphrase)

                    if raw_text is None:
                        record = {
                            "trial_key": trial_key,
                            "group_id": group_id,
                            "paraphrase_idx": para_idx,
                            "model": model_name,
                            "rep": rep,
                            "raw_text": None,
                            "parse_ok": False,
                            "tool_name_raw": None,
                            "tool_name_resolved": None,
                            "tool_correct": False,
                            "params": {},
                            "pydantic_valid": False,
                            "pydantic_error": "no_response",
                        }
                        append_raw(RAW_PATH, record)
                        print("FAIL (no response)")
                        time.sleep(DELAY)
                        continue

                    parsed = parse_json(raw_text)
                    parse_ok = parsed is not None

                    tool_name_raw = None
                    tool_name_resolved = None
                    tool_correct = False
                    params = {}
                    pydantic_valid = False
                    pydantic_error = ""

                    if parse_ok:
                        tool_name_raw = parsed.get("tool_name")
                        tool_name_resolved = resolve_tool_name(tool_name_raw)
                        tool_correct = tool_name_resolved == expected_tool
                        params = parsed.get("parameters", {})

                        if tool_name_resolved:
                            pydantic_valid, pydantic_error = validate_with_pydantic(
                                tool_name_resolved, params
                            )

                    record = {
                        "trial_key": trial_key,
                        "group_id": group_id,
                        "paraphrase_idx": para_idx,
                        "model": model_name,
                        "rep": rep,
                        "raw_text": raw_text[:2000],
                        "parse_ok": parse_ok,
                        "tool_name_raw": tool_name_raw,
                        "tool_name_resolved": tool_name_resolved,
                        "tool_correct": tool_correct,
                        "params": params,
                        "pydantic_valid": pydantic_valid,
                        "pydantic_error": pydantic_error,
                    }
                    append_raw(RAW_PATH, record)

                    status = "OK" if tool_correct else f"WRONG({tool_name_raw})"
                    if not parse_ok:
                        status = "PARSE_FAIL"
                    print(status)

                    time.sleep(DELAY)

    print(f"\nDone: {total} new trials, {skipped} skipped (checkpoint)")
    print("Generating summary CSV ...")
    generate_summary()


# ---------------------------------------------------------------------------
# Summary Generation
# ---------------------------------------------------------------------------

def generate_summary():
    """Aggregate raw JSONL into per-(group, paraphrase_idx, model) CSV rows."""
    # Load all raw records
    records = []
    with open(RAW_PATH) as f:
        for line in f:
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue

    # Group by (group_id, paraphrase_idx, model)
    from collections import defaultdict
    cells = defaultdict(list)
    for rec in records:
        key = (rec["group_id"], rec["paraphrase_idx"], rec["model"])
        cells[key].append(rec)

    rows = []
    for (group_id, para_idx, model), trials in sorted(cells.items()):
        n = len(trials)
        parse_rate = sum(1 for t in trials if t["parse_ok"]) / n if n > 0 else 0
        tool_correct_rate = sum(1 for t in trials if t["tool_correct"]) / n if n > 0 else 0
        pydantic_rate = sum(1 for t in trials if t["pydantic_valid"]) / n if n > 0 else 0

        # Constrained param values (for consistency computation)
        group = next(g for g in PROMPT_GROUPS if g["id"] == group_id)
        constrained_vals = {}
        for param_name in group["constrained_params"]:
            vals = [t["params"].get(param_name) for t in trials if t["tool_correct"] and t["parse_ok"]]
            constrained_vals[param_name] = compute_consistency(vals) if vals else 0.0

        constrained_consistency = (
            sum(constrained_vals.values()) / len(constrained_vals)
            if constrained_vals else 0.0
        )

        rows.append({
            "group_id": group_id,
            "paraphrase_idx": para_idx,
            "model": model,
            "n_trials": n,
            "parse_rate": round(parse_rate, 4),
            "tool_correct_rate": round(tool_correct_rate, 4),
            "pydantic_valid_rate": round(pydantic_rate, 4),
            "constrained_param_consistency": round(constrained_consistency, 4),
        })

    # Write CSV
    fieldnames = [
        "group_id", "paraphrase_idx", "model", "n_trials",
        "parse_rate", "tool_correct_rate", "pydantic_valid_rate",
        "constrained_param_consistency",
    ]
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {CSV_PATH}")

    # Write summary text
    lines = []
    lines.append("Prompt Sensitivity Experiment -- Summary")
    lines.append("=" * 60)
    lines.append("")

    # Overall stats
    total_trials = len(records)
    total_parse = sum(1 for r in records if r["parse_ok"])
    total_correct = sum(1 for r in records if r["tool_correct"])
    total_pydantic = sum(1 for r in records if r["pydantic_valid"])

    lines.append(f"Total trials: {total_trials}")
    lines.append(f"Parse rate: {total_parse / total_trials * 100:.1f}%")
    lines.append(f"Tool selection accuracy: {total_correct / total_trials * 100:.1f}%")
    lines.append(f"Pydantic validation rate: {total_pydantic / total_trials * 100:.1f}%")
    lines.append("")

    # Per-model summary
    lines.append("PER-MODEL SUMMARY")
    lines.append("-" * 40)
    for model in MODELS:
        model_trials = [r for r in records if r["model"] == model]
        n = len(model_trials)
        if n == 0:
            continue
        correct = sum(1 for t in model_trials if t["tool_correct"])
        pydantic = sum(1 for t in model_trials if t["pydantic_valid"])
        lines.append(f"  {model}: tool_correct={correct}/{n} ({correct/n*100:.1f}%), "
                      f"pydantic={pydantic}/{n} ({pydantic/n*100:.1f}%)")
    lines.append("")

    # Per-group summary
    lines.append("PER-GROUP SUMMARY")
    lines.append("-" * 40)
    for group in PROMPT_GROUPS:
        gid = group["id"]
        group_trials = [r for r in records if r["group_id"] == gid]
        n = len(group_trials)
        if n == 0:
            continue
        correct = sum(1 for t in group_trials if t["tool_correct"])
        pydantic = sum(1 for t in group_trials if t["pydantic_valid"])
        lines.append(f"  {gid}: tool_correct={correct}/{n} ({correct/n*100:.1f}%), "
                      f"pydantic={pydantic}/{n} ({pydantic/n*100:.1f}%)")
    lines.append("")

    # Cross-paraphrase consistency per group per model
    lines.append("CROSS-PARAPHRASE TOOL SELECTION (per model x group)")
    lines.append("-" * 60)
    for group in PROMPT_GROUPS:
        gid = group["id"]
        for model in MODELS:
            group_model = [r for r in records if r["group_id"] == gid and r["model"] == model]
            n = len(group_model)
            correct = sum(1 for t in group_model if t["tool_correct"])
            lines.append(f"  {gid} x {model}: {correct}/{n}")
    lines.append("")

    summary_text = "\n".join(lines)
    with open(SUMMARY_PATH, "w") as f:
        f.write(summary_text)
    print(f"Saved: {SUMMARY_PATH}")
    print(summary_text)


if __name__ == "__main__":
    run_experiment()
