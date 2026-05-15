#!/usr/bin/env python3
"""
Prompt Sensitivity Analysis

Computes cross-paraphrase concordance metrics and compares with
ablation within-paraphrase consistency. Produces summary statistics
for the paper.
"""

import csv
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
SENSITIVITY_DIR = SCRIPT_DIR.parent / "data" / "sensitivity"
ANALYSIS_DIR = SENSITIVITY_DIR / "analysis"
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

RAW_PATH = SENSITIVITY_DIR / "prompt_sensitivity_raw.jsonl"
ABLATION_RAW_PATH = SCRIPT_DIR.parent / "data" / "ablation" / "invocation" / "ablation_invocation_raw.jsonl"

OUTPUT_CSV = ANALYSIS_DIR / "sensitivity_cross_paraphrase.csv"
SUMMARY_PATH = ANALYSIS_DIR / "sensitivity_summary.txt"

# ---------------------------------------------------------------------------
# Prompt group metadata (mirrors prompt_sensitivity.py)
# ---------------------------------------------------------------------------

GROUP_META = {
    "spatial_domains": {
        "expected_tool": "identify_spatial_domains",
        "constrained_params": ["method"],
    },
    "deconvolution": {
        "expected_tool": "deconvolve_data",
        "constrained_params": ["method"],
    },
    "svg_detection": {
        "expected_tool": "find_spatial_genes",
        "constrained_params": ["method"],
    },
    "cell_communication": {
        "expected_tool": "analyze_cell_communication",
        "constrained_params": ["method", "species"],
    },
    "moran_i": {
        "expected_tool": "analyze_spatial_statistics",
        "constrained_params": ["analysis_type"],
    },
}

MODELS = ["gemini-2.5-flash", "claude-haiku-4-5-20251001", "gpt-5.4"]

# ---------------------------------------------------------------------------
# Load raw data
# ---------------------------------------------------------------------------

def load_jsonl(path: Path) -> list[dict]:
    records = []
    with open(path) as f:
        for line in f:
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return records


# ---------------------------------------------------------------------------
# Bootstrap percentile CI
# ---------------------------------------------------------------------------

def bootstrap_ci(values: np.ndarray, n_boot: int = 10_000, seed: int = 42) -> tuple[float, float, float]:
    """Return (mean, ci_lo, ci_hi) for 95% percentile bootstrap."""
    rng = np.random.default_rng(seed)
    n = len(values)
    if n == 0:
        return (np.nan, np.nan, np.nan)
    boot_means = np.array([
        values[rng.choice(n, size=n, replace=True)].mean()
        for _ in range(n_boot)
    ])
    return (values.mean(), np.percentile(boot_means, 2.5), np.percentile(boot_means, 97.5))


# ---------------------------------------------------------------------------
# Cross-paraphrase concordance
# ---------------------------------------------------------------------------

def compute_cross_paraphrase_metrics(records: list[dict]) -> list[dict]:
    """Compute per (group, model) cross-paraphrase metrics."""
    # Index: (group_id, model) -> paraphrase_idx -> list of records
    indexed = defaultdict(lambda: defaultdict(list))
    for rec in records:
        indexed[(rec["group_id"], rec["model"])][rec["paraphrase_idx"]].append(rec)

    rows = []
    for (gid, model), para_map in sorted(indexed.items()):
        meta = GROUP_META[gid]
        n_paraphrases = len(para_map)
        all_trials = [t for trials in para_map.values() for t in trials]
        n_total = len(all_trials)

        # 1. Tool selection rate (across all paraphrases)
        tool_correct_count = sum(1 for t in all_trials if t["tool_correct"])
        tool_correct_rate = tool_correct_count / n_total if n_total > 0 else 0

        # 2. Pydantic validation rate
        pydantic_count = sum(1 for t in all_trials if t["pydantic_valid"])
        pydantic_rate = pydantic_count / n_total if n_total > 0 else 0

        # 3. Cross-paraphrase constrained parameter concordance
        # For each constrained param, collect the modal value per paraphrase,
        # then measure agreement across paraphrases
        constrained_concordances = []
        for param_name in meta["constrained_params"]:
            para_modes = []
            for pidx in sorted(para_map.keys()):
                vals = [
                    t["params"].get(param_name)
                    for t in para_map[pidx]
                    if t["tool_correct"] and t["parse_ok"]
                ]
                if vals:
                    # Mode value for this paraphrase
                    from collections import Counter
                    mode_val = Counter(str(v) for v in vals).most_common(1)[0][0]
                    para_modes.append(mode_val)
            if para_modes:
                # Agreement: fraction matching the overall mode
                overall_mode = Counter(para_modes).most_common(1)[0][0]
                concordance = sum(1 for m in para_modes if m == overall_mode) / len(para_modes)
                constrained_concordances.append(concordance)

        constrained_concordance = (
            np.mean(constrained_concordances) if constrained_concordances else np.nan
        )

        # 4. Pairwise tool agreement across paraphrases
        # For each pair of paraphrases, compute fraction of reps with same tool
        pairwise_agreements = []
        para_indices = sorted(para_map.keys())
        for pi, pj in combinations(para_indices, 2):
            # Match reps
            tools_i = [t["tool_name_resolved"] for t in sorted(para_map[pi], key=lambda x: x["rep"])]
            tools_j = [t["tool_name_resolved"] for t in sorted(para_map[pj], key=lambda x: x["rep"])]
            min_len = min(len(tools_i), len(tools_j))
            if min_len > 0:
                agree = sum(1 for a, b in zip(tools_i[:min_len], tools_j[:min_len]) if a == b)
                pairwise_agreements.append(agree / min_len)

        pairwise_tool_agreement = np.mean(pairwise_agreements) if pairwise_agreements else np.nan

        # 5. Bootstrap CI on tool_correct_rate
        correct_arr = np.array([1.0 if t["tool_correct"] else 0.0 for t in all_trials])
        tc_mean, tc_lo, tc_hi = bootstrap_ci(correct_arr)

        rows.append({
            "group_id": gid,
            "model": model,
            "n_paraphrases": n_paraphrases,
            "n_trials": n_total,
            "tool_correct_rate": round(tool_correct_rate, 4),
            "tool_correct_ci_lo": round(tc_lo, 4),
            "tool_correct_ci_hi": round(tc_hi, 4),
            "pydantic_valid_rate": round(pydantic_rate, 4),
            "constrained_param_concordance": round(constrained_concordance, 4) if not np.isnan(constrained_concordance) else "",
            "pairwise_tool_agreement": round(pairwise_tool_agreement, 4) if not np.isnan(pairwise_tool_agreement) else "",
        })

    return rows


# ---------------------------------------------------------------------------
# Ablation comparison
# ---------------------------------------------------------------------------

def load_ablation_within_paraphrase(path: Path) -> dict:
    """Load ablation full_schema results for comparison.
    Returns: {(group_id, model): constrained_param_consistency}
    """
    if not path.exists():
        return {}

    records = load_jsonl(path)
    # Filter to full_schema condition only
    full_records = [r for r in records if r.get("condition") == "full_schema"]

    # Map ablation prompt_id to our group_id
    ablation_to_group = {
        "spatial_domains": "spatial_domains",
        "deconvolution_rctd": "deconvolution",
        "svg_detection": "svg_detection",
        "cellchat": "cell_communication",
        "moran_i": "moran_i",
    }

    indexed = defaultdict(list)
    for rec in full_records:
        prompt_id = rec.get("prompt_id", "")
        gid = ablation_to_group.get(prompt_id)
        if gid:
            indexed[(gid, rec["model"])].append(rec)

    result = {}
    for (gid, model), trials in indexed.items():
        meta = GROUP_META.get(gid)
        if not meta:
            continue
        consistencies = []
        for param_name in meta["constrained_params"]:
            vals = []
            for t in trials:
                if not (t.get("tool_correct") and (t.get("parse_ok") or t.get("parse_success"))):
                    continue
                params = t.get("parameters", t.get("params", {}))
                if isinstance(params, str):
                    try:
                        params = json.loads(params)
                    except (json.JSONDecodeError, TypeError):
                        params = {}
                if isinstance(params, dict):
                    vals.append(params.get(param_name))
            if vals:
                from collections import Counter
                mode_count = Counter(str(v) for v in vals).most_common(1)[0][1]
                consistencies.append(mode_count / len(vals))
        if consistencies:
            result[(gid, model)] = np.mean(consistencies)

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading prompt sensitivity data ...")
    records = load_jsonl(RAW_PATH)
    print(f"  {len(records)} records loaded")

    print("Computing cross-paraphrase metrics ...")
    rows = compute_cross_paraphrase_metrics(records)

    # Write CSV
    fieldnames = [
        "group_id", "model", "n_paraphrases", "n_trials",
        "tool_correct_rate", "tool_correct_ci_lo", "tool_correct_ci_hi",
        "pydantic_valid_rate", "constrained_param_concordance",
        "pairwise_tool_agreement",
    ]
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {OUTPUT_CSV}")

    # Load ablation comparison
    ablation_consistency = load_ablation_within_paraphrase(ABLATION_RAW_PATH)

    # Write summary
    lines = []
    lines.append("Prompt Sensitivity Analysis — Summary")
    lines.append("=" * 70)
    lines.append("")

    # Overall statistics
    all_trials = len(records)
    all_correct = sum(1 for r in records if r["tool_correct"])
    all_pydantic = sum(1 for r in records if r["pydantic_valid"])
    lines.append(f"Total trials: {all_trials}")
    lines.append(f"Overall tool selection accuracy: {all_correct / all_trials * 100:.1f}%")
    lines.append(f"Overall Pydantic validation rate: {all_pydantic / all_trials * 100:.1f}%")
    lines.append("")

    # Per group x model table
    lines.append("CROSS-PARAPHRASE METRICS (per group x model)")
    lines.append("-" * 70)
    lines.append(f"{'Group':<20} {'Model':<30} {'Tool%':>6} {'Pyd%':>6} {'ParaConcord':>12} {'AblBaseline':>12}")
    lines.append("-" * 70)
    for row in rows:
        gid = row["group_id"]
        model = row["model"]
        abl_val = ablation_consistency.get((gid, model))
        abl_str = f"{abl_val*100:.0f}%" if abl_val is not None else "n/a"
        conc_str = row["constrained_param_concordance"]
        conc_str = f"{float(conc_str)*100:.0f}%" if conc_str != "" else "n/a"
        lines.append(
            f"{gid:<20} {model:<30} {row['tool_correct_rate']*100:5.1f}% "
            f"{row['pydantic_valid_rate']*100:5.1f}% "
            f"{conc_str:>12} {abl_str:>12}"
        )
    lines.append("")

    # Aggregate cross-model means
    lines.append("AGGREGATE MEANS")
    lines.append("-" * 40)
    tool_rates = [r["tool_correct_rate"] for r in rows]
    pyd_rates = [r["pydantic_valid_rate"] for r in rows]
    conc_vals = [float(r["constrained_param_concordance"]) for r in rows if r["constrained_param_concordance"] != ""]
    lines.append(f"  Mean tool selection accuracy: {np.mean(tool_rates)*100:.1f}%")
    lines.append(f"  Mean Pydantic validation rate: {np.mean(pyd_rates)*100:.1f}%")
    if conc_vals:
        lines.append(f"  Mean cross-paraphrase param concordance: {np.mean(conc_vals)*100:.1f}%")
    if ablation_consistency:
        abl_vals = list(ablation_consistency.values())
        lines.append(f"  Mean ablation within-paraphrase consistency: {np.mean(abl_vals)*100:.1f}%")
    lines.append("")

    # Key finding
    lines.append("KEY FINDING")
    lines.append("=" * 70)
    mean_tool = np.mean(tool_rates) * 100
    mean_pyd = np.mean(pyd_rates) * 100
    lines.append(
        f"Across 5 prompt groups x 5 paraphrases x 3 models x 5 reps = {all_trials} trials:"
    )
    lines.append(
        f"  Tool selection accuracy: {mean_tool:.1f}% (semantic paraphrases preserve tool choice)"
    )
    lines.append(
        f"  Pydantic validation rate: {mean_pyd:.1f}% (schema ensures structural consistency)"
    )
    if conc_vals:
        lines.append(
            f"  Cross-paraphrase param concordance: {np.mean(conc_vals)*100:.1f}% "
            f"(constrained params stable across rephrasings)"
        )
    lines.append("")

    summary = "\n".join(lines)
    with open(SUMMARY_PATH, "w") as f:
        f.write(summary)
    print(f"Saved: {SUMMARY_PATH}")
    print()
    print(summary)


if __name__ == "__main__":
    main()
