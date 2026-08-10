#!/usr/bin/env python3
"""
Schema-Enforcement Ablation Study — Part 3: Analysis

Computes concordance metrics from Part 2 E2E outputs, bootstrap CIs,
and statistical tests comparing the three ablation conditions.

Metrics:
  - Spatial domains: pairwise ARI and NMI
  - SVG detection:   Jaccard similarity of top-K gene lists
  - Deconvolution:   Pearson correlation of proportion matrices

Outputs:
  - ablation_metrics_bootstrap.csv    (all metrics with 95% CIs)
  - ablation_statistical_tests.csv    (Kruskal-Wallis + pairwise Mann-Whitney)
  - ablation_summary.txt              (human-readable summary for paper)
"""

import csv
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import kruskal, mannwhitneyu, pearsonr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
E2E_DIR = SCRIPT_DIR.parent / "data" / "ablation" / "e2e"
INV_DIR = SCRIPT_DIR.parent / "data" / "ablation" / "invocation"
OUT_DIR = SCRIPT_DIR.parent / "data" / "ablation" / "analysis"

E2E_RAW = E2E_DIR / "ablation_e2e_raw.jsonl"
E2E_RESULTS = E2E_DIR / "ablation_e2e_results.csv"
INV_RESULTS = INV_DIR / "ablation_invocation_results.csv"

METRICS_CSV = OUT_DIR / "ablation_metrics_bootstrap.csv"
TESTS_CSV = OUT_DIR / "ablation_statistical_tests.csv"
SUMMARY_PATH = OUT_DIR / "ablation_summary.txt"

CONDITIONS = ["full_schema", "bare_schema", "no_schema"]
TASKS = ["spatial_domains", "svg_detection", "deconvolution"]
TOP_K_SVG = 100  # Jaccard on top-K genes


def committed_outputs_available() -> bool:
    return all(
        path.exists() and path.stat().st_size > 0
        for path in (
            METRICS_CSV,
            TESTS_CSV,
            SUMMARY_PATH,
        )
    )


# ---------------------------------------------------------------------------
# Bootstrap CI
# ---------------------------------------------------------------------------


def bootstrap_ci(
    values: list[float],
    n_bootstrap: int = 10000,
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
# Concordance Metrics
# ---------------------------------------------------------------------------


def pairwise_ari_nmi(label_lists: list[list]) -> dict:
    """Compute pairwise ARI and NMI for lists of cluster labels."""
    aris, nmis = [], []
    for i, j in combinations(range(len(label_lists)), 2):
        aris.append(adjusted_rand_score(label_lists[i], label_lists[j]))
        nmis.append(
            normalized_mutual_info_score(
                label_lists[i], label_lists[j], average_method="arithmetic"
            )
        )
    return {"ari_values": aris, "nmi_values": nmis}


def jaccard_top_k(gene_lists: list[list[str]], k: int = TOP_K_SVG) -> list[float]:
    """Compute pairwise Jaccard similarity of top-K gene lists."""
    jaccards = []
    for i, j in combinations(range(len(gene_lists)), 2):
        set_i = set(gene_lists[i][:k])
        set_j = set(gene_lists[j][:k])
        union = len(set_i | set_j)
        jaccards.append(len(set_i & set_j) / union if union > 0 else 0.0)
    return jaccards


def pairwise_pearson(prop_matrices: list[np.ndarray]) -> list[float]:
    """Compute pairwise Pearson r of flattened proportion matrices."""
    correlations = []
    for i, j in combinations(range(len(prop_matrices)), 2):
        flat_i = prop_matrices[i].flatten()
        flat_j = prop_matrices[j].flatten()
        if len(flat_i) < 2 or np.std(flat_i) == 0 or np.std(flat_j) == 0:
            correlations.append(1.0 if np.array_equal(flat_i, flat_j) else 0.0)
        else:
            r, _ = pearsonr(flat_i, flat_j)
            correlations.append(r)
    return correlations


# ---------------------------------------------------------------------------
# Load E2E Outputs
# ---------------------------------------------------------------------------


def load_e2e_outputs() -> dict:
    """Load E2E experiment raw records grouped by (condition, model, task)."""
    groups = defaultdict(list)
    if not E2E_RAW.exists():
        print(f"WARNING: {E2E_RAW} not found. Run ablation_e2e.py first.")
        return groups

    with open(E2E_RAW) as f:
        for line in f:
            try:
                rec = json.loads(line)
                key = (rec["condition"], rec.get("model", ""), rec["task_id"])
                groups[key].append(rec)
            except (json.JSONDecodeError, KeyError):
                continue
    return groups


def load_artifacts(records: list[dict], task_id: str):
    """Load output artifacts from saved JSON files for successful executions."""
    artifacts = []
    for rec in records:
        if not rec.get("exec_success"):
            continue
        path = rec.get("artifact_path", "")
        if not path or not Path(path).exists():
            continue
        with open(path) as f:
            data = json.load(f)
        artifacts.append(data)
    return artifacts


# ---------------------------------------------------------------------------
# Compute Per-Cell Concordance
# ---------------------------------------------------------------------------


def compute_concordance(artifacts: list[dict], task_id: str) -> dict:
    """Compute concordance metrics for a list of artifacts from the same cell."""
    if len(artifacts) < 2:
        return {"metric": task_id, "n_runs": len(artifacts), "n_pairs": 0}

    if task_id == "spatial_domains":
        labels = [a["domain_labels"] for a in artifacts]
        pw = pairwise_ari_nmi(labels)
        ari_mean, ari_lo, ari_hi = bootstrap_ci(pw["ari_values"])
        nmi_mean, nmi_lo, nmi_hi = bootstrap_ci(pw["nmi_values"])
        return {
            "metric": "spatial_domains",
            "n_runs": len(labels),
            "n_pairs": len(pw["ari_values"]),
            "ari_mean": round(ari_mean, 4),
            "ari_ci_lo": round(ari_lo, 4),
            "ari_ci_hi": round(ari_hi, 4),
            "nmi_mean": round(nmi_mean, 4),
            "nmi_ci_lo": round(nmi_lo, 4),
            "nmi_ci_hi": round(nmi_hi, 4),
            "ari_values": pw["ari_values"],
            "nmi_values": pw["nmi_values"],
        }

    elif task_id == "svg_detection":
        gene_lists = [a["spatial_genes"] for a in artifacts]
        jaccards = jaccard_top_k(gene_lists, TOP_K_SVG)
        j_mean, j_lo, j_hi = bootstrap_ci(jaccards)
        return {
            "metric": "svg_detection",
            "n_runs": len(gene_lists),
            "n_pairs": len(jaccards),
            "jaccard_mean": round(j_mean, 4),
            "jaccard_ci_lo": round(j_lo, 4),
            "jaccard_ci_hi": round(j_hi, 4),
            "jaccard_values": jaccards,
        }

    elif task_id == "deconvolution":
        props = [np.array(a["proportions"]) for a in artifacts]
        pearsons = pairwise_pearson(props)
        p_mean, p_lo, p_hi = bootstrap_ci(pearsons)
        return {
            "metric": "deconvolution",
            "n_runs": len(props),
            "n_pairs": len(pearsons),
            "pearson_mean": round(p_mean, 4),
            "pearson_ci_lo": round(p_lo, 4),
            "pearson_ci_hi": round(p_hi, 4),
            "pearson_values": pearsons,
        }

    return {"metric": task_id, "n_runs": len(artifacts), "n_pairs": 0}


# ---------------------------------------------------------------------------
# Statistical Tests
# ---------------------------------------------------------------------------


def run_statistical_tests(concordance_by_condition: dict) -> list[dict]:
    """Kruskal-Wallis across 3 conditions + pairwise Mann-Whitney U."""
    test_results = []

    # Determine which metric keys to use
    for metric_key in ["ari_values", "nmi_values", "jaccard_values", "pearson_values"]:
        data_by_cond = {}
        for cond in CONDITIONS:
            vals = concordance_by_condition.get(cond, {}).get(metric_key, [])
            if vals:
                data_by_cond[cond] = vals

        if len(data_by_cond) < 2:
            continue

        # Kruskal-Wallis
        groups = list(data_by_cond.values())
        if len(groups) >= 2 and all(len(g) >= 2 for g in groups):
            try:
                kw_stat, kw_p = kruskal(*groups)
            except ValueError:
                kw_stat, kw_p = np.nan, np.nan
        else:
            kw_stat, kw_p = np.nan, np.nan

        test_results.append(
            {
                "metric": metric_key.replace("_values", ""),
                "test": "kruskal_wallis",
                "groups": list(data_by_cond.keys()),
                "statistic": round(kw_stat, 4) if not np.isnan(kw_stat) else "NA",
                "p_value": round(kw_p, 6) if not np.isnan(kw_p) else "NA",
                "significant": kw_p < 0.05 if not np.isnan(kw_p) else "NA",
            }
        )

        # Pairwise Mann-Whitney U with Bonferroni correction
        cond_list = list(data_by_cond.keys())
        n_tests = len(list(combinations(cond_list, 2)))
        for c1, c2 in combinations(cond_list, 2):
            g1 = data_by_cond[c1]
            g2 = data_by_cond[c2]
            if len(g1) >= 2 and len(g2) >= 2:
                try:
                    u_stat, u_p = mannwhitneyu(g1, g2, alternative="two-sided")
                except ValueError:
                    u_stat, u_p = np.nan, np.nan
            else:
                u_stat, u_p = np.nan, np.nan

            bonferroni_alpha = 0.05 / n_tests if n_tests > 0 else 0.05
            test_results.append(
                {
                    "metric": metric_key.replace("_values", ""),
                    "test": f"mann_whitney_{c1}_vs_{c2}",
                    "groups": [c1, c2],
                    "statistic": round(u_stat, 4) if not np.isnan(u_stat) else "NA",
                    "p_value": round(u_p, 6) if not np.isnan(u_p) else "NA",
                    "significant": (
                        u_p < bonferroni_alpha if not np.isnan(u_p) else "NA"
                    ),
                    "bonferroni_alpha": round(bonferroni_alpha, 4),
                }
            )

    return test_results


# ---------------------------------------------------------------------------
# Main Analysis
# ---------------------------------------------------------------------------


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Schema-Enforcement Ablation — Part 3: Analysis")
    print(f"{'='*60}\n")

    if not E2E_RAW.exists():
        if committed_outputs_available():
            print(f"Raw checkpoint not found: {E2E_RAW}")
            print("Using committed aggregate outputs:")
            print(f"  {METRICS_CSV}")
            print(f"  {TESTS_CSV}")
            print(f"  {SUMMARY_PATH}")
            print(
                "To recompute raw-level concordance, run scripts/ablation_e2e.py first."
            )
            return
        raise FileNotFoundError(
            f"{E2E_RAW} is required to recompute ablation concordance. "
            "Run scripts/ablation_e2e.py first, or use the committed aggregate CSV/TXT outputs."
        )

    # ---- Load invocation-level results ----
    inv_data = {}
    if INV_RESULTS.exists():
        with open(INV_RESULTS) as f:
            reader = csv.DictReader(f)
            for row in reader:
                cond = row["condition"]
                if cond not in inv_data:
                    inv_data[cond] = []
                inv_data[cond].append(row)
        print(
            f"Loaded invocation results: {sum(len(v) for v in inv_data.values())} rows"
        )
    else:
        print("No invocation results found.")

    # ---- Load E2E results ----
    groups = load_e2e_outputs()
    print(
        f"Loaded E2E results: {sum(len(v) for v in groups.values())} records "
        f"across {len(groups)} cells\n"
    )

    # ---- Compute concordance per (condition, task) aggregated across models ----
    concordance_results = []  # For CSV
    concordance_by_cond_task = defaultdict(lambda: defaultdict(list))  # For stats

    for task_id in TASKS:
        print(f"\n--- {task_id.upper()} ---")
        for cond in CONDITIONS:
            # Aggregate across all models for this condition+task
            all_artifacts = []
            for model in ["gemini-2.5-flash", "claude-haiku-4-5-20251001", "gpt-5.4"]:
                key = (cond, model, task_id)
                if key in groups:
                    arts = load_artifacts(groups[key], task_id)
                    all_artifacts.extend(arts)

            # Also compute per-model concordance for detail
            per_model_concordances = []
            for model in ["gemini-2.5-flash", "claude-haiku-4-5-20251001", "gpt-5.4"]:
                key = (cond, model, task_id)
                if key not in groups:
                    continue
                arts = load_artifacts(groups[key], task_id)
                if len(arts) >= 2:
                    conc = compute_concordance(arts, task_id)
                    conc["condition"] = cond
                    conc["model"] = model
                    conc["task"] = task_id
                    concordance_results.append(conc)
                    per_model_concordances.append(conc)

            # Aggregate concordance across all models
            if len(all_artifacts) >= 2:
                agg = compute_concordance(all_artifacts, task_id)
                agg["condition"] = cond
                agg["model"] = "all_models"
                agg["task"] = task_id
                concordance_results.append(agg)

                # Store raw values for statistical tests
                for k in [
                    "ari_values",
                    "nmi_values",
                    "jaccard_values",
                    "pearson_values",
                ]:
                    if k in agg:
                        concordance_by_cond_task[(task_id, cond)][k] = agg[k]

            exec_total = sum(
                sum(
                    1
                    for r in groups.get((cond, m, task_id), [])
                    if r.get("exec_success")
                )
                for m in ["gemini-2.5-flash", "claude-haiku-4-5-20251001", "gpt-5.4"]
            )
            total = sum(
                len(groups.get((cond, m, task_id), []))
                for m in ["gemini-2.5-flash", "claude-haiku-4-5-20251001", "gpt-5.4"]
            )
            exec_rate = exec_total / total if total else 0
            print(
                f"  {cond}: {exec_total}/{total} exec ({exec_rate:.0%}), {len(all_artifacts)} artifacts"
            )

    # ---- Write concordance CSV ----
    if concordance_results:
        # Flatten: remove raw value lists for CSV
        csv_rows = []
        for r in concordance_results:
            row = {k: v for k, v in r.items() if not k.endswith("_values")}
            csv_rows.append(row)

        # Collect all fieldnames
        all_fields = set()
        for row in csv_rows:
            all_fields.update(row.keys())
        ordered_fields = ["condition", "model", "task", "metric", "n_runs", "n_pairs"]
        for f in sorted(all_fields):
            if f not in ordered_fields:
                ordered_fields.append(f)

        with open(METRICS_CSV, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=ordered_fields,
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"\nMetrics CSV: {METRICS_CSV}")

    # ---- Statistical tests per task ----
    print(f"\n{'='*60}")
    print("STATISTICAL TESTS")
    print(f"{'='*60}")

    all_tests = []
    for task_id in TASKS:
        cond_data = {}
        for cond in CONDITIONS:
            cond_data[cond] = concordance_by_cond_task.get((task_id, cond), {})
        if cond_data:
            tests = run_statistical_tests(cond_data)
            for t in tests:
                t["task"] = task_id
            all_tests.extend(tests)
            for t in tests:
                print(
                    f"  [{task_id}] {t['test']}: stat={t['statistic']} p={t['p_value']} sig={t['significant']}"
                )

    if all_tests:
        test_fields = [
            "task",
            "metric",
            "test",
            "groups",
            "statistic",
            "p_value",
            "significant",
        ]
        extra = set()
        for t in all_tests:
            extra.update(t.keys())
        for f in sorted(extra):
            if f not in test_fields:
                test_fields.append(f)

        with open(TESTS_CSV, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=test_fields,
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()
            for t in all_tests:
                row = {
                    k: (json.dumps(v) if isinstance(v, list) else v)
                    for k, v in t.items()
                }
                writer.writerows([row])
        print(f"\nTests CSV: {TESTS_CSV}")

    # ---- Human-readable summary ----
    with open(SUMMARY_PATH, "w") as f:
        f.write("Schema-Enforcement Ablation Study — Summary\n")
        f.write(f"{'='*70}\n\n")

        # Invocation-level summary
        if inv_data:
            f.write("PART 1: INVOCATION-LEVEL ABLATION\n")
            f.write(f"{'-'*50}\n")
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
                for cond in CONDITIONS:
                    rows = inv_data.get(cond, [])
                    if rows:
                        vals[cond] = sum(float(r[metric]) for r in rows) / len(rows)
                    else:
                        vals[cond] = 0
                f.write(
                    f"{metric:<35} {vals.get('full_schema',0):>9.1%} "
                    f"{vals.get('bare_schema',0):>9.1%} "
                    f"{vals.get('no_schema',0):>9.1%}\n"
                )
            f.write("\n")

        # E2E summary
        f.write("\nPART 2: END-TO-END OUTPUT ABLATION\n")
        f.write(f"{'-'*50}\n\n")

        for task_id in TASKS:
            f.write(f"\n  {task_id.upper()}\n")
            f.write(f"  {'~'*40}\n")
            for cond in CONDITIONS:
                agg = next(
                    (
                        r
                        for r in concordance_results
                        if r.get("condition") == cond
                        and r.get("model") == "all_models"
                        and r.get("task") == task_id
                    ),
                    None,
                )
                if agg is None:
                    f.write(f"  {cond}: no data\n")
                    continue

                f.write(
                    f"  {cond} (n_runs={agg.get('n_runs', 0)}, n_pairs={agg.get('n_pairs', 0)}):\n"
                )
                if task_id == "spatial_domains":
                    f.write(
                        f"    ARI: {agg.get('ari_mean','NA')} [{agg.get('ari_ci_lo','NA')}, {agg.get('ari_ci_hi','NA')}]\n"
                    )
                    f.write(
                        f"    NMI: {agg.get('nmi_mean','NA')} [{agg.get('nmi_ci_lo','NA')}, {agg.get('nmi_ci_hi','NA')}]\n"
                    )
                elif task_id == "svg_detection":
                    f.write(
                        f"    Jaccard@{TOP_K_SVG}: {agg.get('jaccard_mean','NA')} [{agg.get('jaccard_ci_lo','NA')}, {agg.get('jaccard_ci_hi','NA')}]\n"
                    )
                elif task_id == "deconvolution":
                    f.write(
                        f"    Pearson r: {agg.get('pearson_mean','NA')} [{agg.get('pearson_ci_lo','NA')}, {agg.get('pearson_ci_hi','NA')}]\n"
                    )

        # Key finding
        f.write(f"\n\n{'='*70}\n")
        f.write("KEY FINDING\n")
        f.write(f"{'='*70}\n\n")
        f.write("If Full Schema >> Bare Schema >> Minimal, this demonstrates:\n")
        f.write(
            "  1. Schema constraints (types + enums) improve invocation consistency\n"
        )
        f.write("  2. Method guidance (descriptions) provide additional improvement\n")
        f.write(
            "  3. The combination produces the most reproducible end-to-end outputs\n"
        )
        f.write(
            "\nThis establishes a causal chain from schema design to output concordance.\n"
        )

    print(f"\nSummary: {SUMMARY_PATH}")
    print("\nDone!")


if __name__ == "__main__":
    main()
