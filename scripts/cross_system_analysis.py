#!/usr/bin/env python3
"""
Cross-System Invocation Comparison — Analysis

Loads cross-system raw data + existing ablation data, computes per-condition
metrics with bootstrap CIs, and generates summary comparing all 5 conditions.

Metrics:
  - Tool selection accuracy (fraction selecting an expected tool)
  - Parse success rate
  - Cross-rep tool consistency (modal agreement)
  - Tool diversity (unique tools selected per prompt × model cell)

Outputs:
  - cross_system_metrics.csv      (per condition × model × prompt)
  - cross_system_aggregate.csv    (per condition, with bootstrap CIs)
  - cross_system_summary.txt      (human-readable summary)
"""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR.parent / "data" / "cross_system"
ANALYSIS_DIR = DATA_DIR / "analysis"
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

CROSS_RAW = DATA_DIR / "cross_system_raw.jsonl"
ABLATION_RAW = SCRIPT_DIR.parent / "data" / "ablation" / "invocation" / "ablation_invocation_raw.jsonl"

METRICS_CSV = ANALYSIS_DIR / "cross_system_metrics.csv"
AGGREGATE_CSV = ANALYSIS_DIR / "cross_system_aggregate.csv"
SUMMARY_PATH = ANALYSIS_DIR / "cross_system_summary.txt"

ALL_CONDITIONS = [
    "full_schema", "bare_schema", "no_schema",
    "stagent_context", "spatialagent_context",
]
CONDITION_LABELS = {
    "full_schema": "ChatSpatial (full)",
    "bare_schema": "ChatSpatial (bare)",
    "no_schema": "ChatSpatial (none)",
    "stagent_context": "STAgent",
    "spatialagent_context": "SpatialAgent",
}

PROMPT_IDS = [
    "spatial_domains", "deconvolution_rctd", "svg_detection", "cellchat",
    "moran_i", "spatial_plot", "trajectory", "cnv",
]


def committed_outputs_available() -> bool:
    return all(path.exists() and path.stat().st_size > 0 for path in (
        METRICS_CSV,
        AGGREGATE_CSV,
        SUMMARY_PATH,
    ))

# ---------------------------------------------------------------------------
# Bootstrap CI
# ---------------------------------------------------------------------------

def bootstrap_ci(
    values: list[float],
    n_bootstrap: int = 10_000,
    ci: float = 0.95,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Percentile bootstrap for mean. Returns (mean, ci_lower, ci_upper)."""
    if not values:
        return (np.nan, np.nan, np.nan)
    arr = np.array(values, dtype=float)
    if len(arr) == 1:
        return (arr[0], arr[0], arr[0])
    rng = np.random.default_rng(seed)
    boot_means = [
        np.mean(rng.choice(arr, size=len(arr), replace=True))
        for _ in range(n_bootstrap)
    ]
    alpha = (1 - ci) / 2
    return (
        float(np.mean(arr)),
        float(np.percentile(boot_means, alpha * 100)),
        float(np.percentile(boot_means, (1 - alpha) * 100)),
    )


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_jsonl(path: Path) -> list[dict]:
    records = []
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return records
    with open(path) as f:
        for line in f:
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return records


def normalize_records(raw: list[dict]) -> list[dict]:
    """Normalize field names across ablation and cross-system data."""
    out = []
    for rec in raw:
        out.append({
            "condition": rec["condition"],
            "model": rec["model"],
            "prompt_id": rec["prompt_id"],
            "rep": rec.get("rep", 0),
            "parse_success": rec.get("parse_success", rec.get("parse_ok", False)),
            "canonical_tool": rec.get("canonical_tool", rec.get("tool_name_resolved")),
            "tool_correct": rec.get("tool_correct", False),
        })
    return out


# ---------------------------------------------------------------------------
# Compute metrics
# ---------------------------------------------------------------------------

def compute_metrics(records: list[dict]) -> list[dict]:
    """Per (condition, model, prompt_id) metrics."""
    groups = defaultdict(list)
    for rec in records:
        key = (rec["condition"], rec["model"], rec["prompt_id"])
        groups[key].append(rec)

    rows = []
    for (condition, model, prompt_id), recs in sorted(groups.items()):
        n = len(recs)
        parsed = [r for r in recs if r["parse_success"]]
        correct = [r for r in recs if r["tool_correct"]]

        parse_rate = len(parsed) / n if n else 0
        correct_rate = len(correct) / n if n else 0

        # Tool consistency: fraction agreeing with mode
        tools = [r["canonical_tool"] for r in parsed if r["canonical_tool"]]
        if tools:
            mode_count = Counter(tools).most_common(1)[0][1]
            consistency = mode_count / len(tools)
        else:
            consistency = 0

        # Tool diversity: number of unique tools selected
        unique_tools = len(set(tools)) if tools else 0

        rows.append({
            "condition": condition,
            "model": model,
            "prompt_id": prompt_id,
            "n_trials": n,
            "parse_rate": round(parse_rate, 4),
            "tool_correct_rate": round(correct_rate, 4),
            "tool_consistency": round(consistency, 4),
            "tool_diversity": unique_tools,
        })

    return rows


def compute_aggregate(metric_rows: list[dict]) -> list[dict]:
    """Per-condition aggregate with bootstrap CIs on tool_correct_rate."""
    agg = []
    for cond in ALL_CONDITIONS:
        cond_rows = [r for r in metric_rows if r["condition"] == cond]
        if not cond_rows:
            continue

        # Collect all per-cell correct rates for bootstrap
        correct_rates = [r["tool_correct_rate"] for r in cond_rows]
        consistency_rates = [r["tool_consistency"] for r in cond_rows]
        parse_rates = [r["parse_rate"] for r in cond_rows]

        corr_mean, corr_lo, corr_hi = bootstrap_ci(correct_rates)
        cons_mean, cons_lo, cons_hi = bootstrap_ci(consistency_rates)
        parse_mean, parse_lo, parse_hi = bootstrap_ci(parse_rates)

        total_n = sum(r["n_trials"] for r in cond_rows)
        avg_diversity = np.mean([r["tool_diversity"] for r in cond_rows])

        agg.append({
            "condition": cond,
            "label": CONDITION_LABELS.get(cond, cond),
            "n_trials": total_n,
            "n_cells": len(cond_rows),
            "parse_rate_mean": round(parse_mean, 4),
            "parse_rate_ci_lo": round(parse_lo, 4),
            "parse_rate_ci_hi": round(parse_hi, 4),
            "tool_correct_mean": round(corr_mean, 4),
            "tool_correct_ci_lo": round(corr_lo, 4),
            "tool_correct_ci_hi": round(corr_hi, 4),
            "consistency_mean": round(cons_mean, 4),
            "consistency_ci_lo": round(cons_lo, 4),
            "consistency_ci_hi": round(cons_hi, 4),
            "avg_tool_diversity": round(avg_diversity, 2),
        })

    return agg


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading data ...")
    missing_raw = [path for path in (ABLATION_RAW, CROSS_RAW) if not path.exists()]
    if missing_raw:
        if committed_outputs_available():
            print("Raw checkpoints not found:")
            for path in missing_raw:
                print(f"  {path}")
            print("Using committed aggregate outputs:")
            print(f"  {METRICS_CSV}")
            print(f"  {AGGREGATE_CSV}")
            print(f"  {SUMMARY_PATH}")
            print(
                "To recompute raw-level metrics, run scripts/ablation_invocation.py "
                "and scripts/cross_system_comparison.py first."
            )
            return
        raise FileNotFoundError(
            "Raw JSONL checkpoints are required to recompute cross-system metrics. "
            "Run scripts/ablation_invocation.py and scripts/cross_system_comparison.py first, "
            "or use the committed aggregate CSV/TXT outputs."
        )

    abl_raw = load_jsonl(ABLATION_RAW)
    cross_raw = load_jsonl(CROSS_RAW)
    print(f"  Ablation records: {len(abl_raw)}")
    print(f"  Cross-system records: {len(cross_raw)}")

    all_records = normalize_records(abl_raw) + normalize_records(cross_raw)
    print(f"  Total records: {len(all_records)}")

    # Per-cell metrics
    print("\nComputing per-cell metrics ...")
    metric_rows = compute_metrics(all_records)
    if not metric_rows:
        raise ValueError("No metric rows were produced from the raw checkpoints.")

    with open(METRICS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=metric_rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(metric_rows)
    print(f"  Saved: {METRICS_CSV}")

    # Aggregate with CIs
    print("Computing aggregate metrics with bootstrap CIs ...")
    agg_rows = compute_aggregate(metric_rows)

    with open(AGGREGATE_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=agg_rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(agg_rows)
    print(f"  Saved: {AGGREGATE_CSV}")

    # Summary
    lines = []
    lines.append("Cross-System Invocation Comparison — Analysis Summary")
    lines.append("=" * 75)
    lines.append("")

    # Aggregate table
    lines.append(
        f"{'Condition':<25} {'N':>5} {'Parse%':>8} "
        f"{'Correct%':>9} {'95% CI':>14} {'Consist%':>10} {'Diversity':>10}"
    )
    lines.append("-" * 85)
    for row in agg_rows:
        ci_str = f"[{row['tool_correct_ci_lo']*100:.1f}, {row['tool_correct_ci_hi']*100:.1f}]"
        lines.append(
            f"{row['label']:<25} {row['n_trials']:>5} "
            f"{row['parse_rate_mean']*100:7.1f}% "
            f"{row['tool_correct_mean']*100:8.1f}% "
            f"{ci_str:>14} "
            f"{row['consistency_mean']*100:9.1f}% "
            f"{row['avg_tool_diversity']:>9.1f}"
        )

    lines.append("")

    # Per-prompt breakdown for new conditions
    for cond in ["stagent_context", "spatialagent_context"]:
        label = CONDITION_LABELS.get(cond, cond)
        lines.append(f"\n{label} — Per-prompt tool selection accuracy:")
        lines.append("-" * 50)
        for pid in PROMPT_IDS:
            prompt_rows = [
                r for r in metric_rows
                if r["condition"] == cond and r["prompt_id"] == pid
            ]
            if not prompt_rows:
                continue
            avg = np.mean([r["tool_correct_rate"] for r in prompt_rows])
            # Show which tools were actually selected
            tools_selected = []
            for rec in all_records:
                if (rec["condition"] == cond and rec["prompt_id"] == pid
                        and rec["parse_success"] and rec["canonical_tool"]):
                    tools_selected.append(rec["canonical_tool"])
            top_tools = Counter(tools_selected).most_common(3)
            top_str = ", ".join(f"{t}({c})" for t, c in top_tools)
            lines.append(f"  {pid:<20} {avg*100:5.1f}%  top: {top_str}")

    lines.append("")

    # Key findings
    lines.append("\nKEY FINDINGS")
    lines.append("=" * 75)

    cs_full = next((r for r in agg_rows if r["condition"] == "full_schema"), None)
    stagent = next((r for r in agg_rows if r["condition"] == "stagent_context"), None)
    spatial = next((r for r in agg_rows if r["condition"] == "spatialagent_context"), None)

    if cs_full and stagent:
        lines.append(
            f"1. ChatSpatial full schema: {cs_full['tool_correct_mean']*100:.1f}% tool accuracy "
            f"vs STAgent: {stagent['tool_correct_mean']*100:.1f}%"
        )
    if cs_full and spatial:
        lines.append(
            f"2. ChatSpatial full schema: {cs_full['tool_correct_mean']*100:.1f}% tool accuracy "
            f"vs SpatialAgent: {spatial['tool_correct_mean']*100:.1f}%"
        )
    if stagent:
        lines.append(
            f"3. STAgent consistency: {stagent['consistency_mean']*100:.1f}% "
            f"(trivially high — most tasks route to python_repl_tool)"
        )
    if spatial:
        lines.append(
            f"4. SpatialAgent diversity: {spatial['avg_tool_diversity']:.1f} tools/cell "
            f"(72 tools makes selection harder)"
        )
    lines.append("")

    summary_text = "\n".join(lines)
    with open(SUMMARY_PATH, "w") as f:
        f.write(summary_text)
    print(f"\n  Saved: {SUMMARY_PATH}")
    print()
    print(summary_text)


if __name__ == "__main__":
    main()
