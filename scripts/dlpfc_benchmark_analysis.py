#!/usr/bin/env python3
"""
DLPFC Ground-Truth Benchmark Analysis

Reads dlpfc_benchmark_raw.jsonl and computes:
  1. Execution success rates per system
  2. ARI and NMI against expert-annotated ground truth
  3. Cross-replicate pairwise ARI (consistency)
  4. Bootstrap 95% CIs for all metrics
"""

import csv
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR.parent / "data" / "dlpfc_benchmark"
SAMPLES_DIR = DATA_DIR / "samples"
GT_DIR = DATA_DIR / "ground_truth"
OUTPUT_DIR = DATA_DIR / "outputs"
RAW_PATH = DATA_DIR / "dlpfc_benchmark_raw.jsonl"
RESULTS_CSV = DATA_DIR / "dlpfc_benchmark_results.csv"
PER_TRIAL_CSV = DATA_DIR / "dlpfc_benchmark_per_trial.csv"
SUMMARY_PATH = DATA_DIR / "dlpfc_benchmark_summary.txt"

SAMPLE_IDS = ["151673", "151507", "151669"]


def committed_outputs_available() -> bool:
    return all(path.exists() and path.stat().st_size > 0 for path in (
        RESULTS_CSV,
        PER_TRIAL_CSV,
        SUMMARY_PATH,
    ))

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_raw():
    records = []
    with open(RAW_PATH) as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def load_ground_truth(sample_id: str) -> pd.Series:
    """Load expert layer annotations for a sample."""
    gt_path = GT_DIR / f"{sample_id}_labels.csv"
    gt = pd.read_csv(gt_path, index_col=0)
    return gt["ground_truth"]


# ---------------------------------------------------------------------------
# Label extraction
# ---------------------------------------------------------------------------

def extract_labels(system: str, sample_id: str, rep: int) -> dict | None:
    """Extract predicted cluster labels from a trial's output.

    Returns dict mapping spot barcode → label string.
    """
    trial_dir = OUTPUT_DIR / system / sample_id / f"rep_{rep}"

    def _from_adata(adata):
        domain_cols = [c for c in adata.obs.columns
                       if any(kw in c.lower()
                              for kw in ["domain", "leiden", "cluster", "louvain"])]
        if domain_cols:
            return dict(zip(adata.obs_names, adata.obs[domain_cols[0]].astype(str)))
        return None

    # Check h5ad files
    for f in sorted(trial_dir.glob("*.h5ad")):
        try:
            adata = ad.read_h5ad(f)
            result = _from_adata(adata)
            if result:
                return result
        except Exception:
            pass

    return None


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def compute_gt_metrics(pred_labels: dict, gt_series: pd.Series) -> dict:
    """Compute ARI and NMI between predicted labels and ground truth.

    Only computes on spots present in both and with non-NA ground truth.
    """
    common = sorted(set(pred_labels.keys()) & set(gt_series.index))
    # Filter out NA ground truth
    common = [s for s in common if pd.notna(gt_series[s])]

    if len(common) < 10:
        return {"ari": np.nan, "nmi": np.nan, "n_spots": len(common)}

    pred = [pred_labels[s] for s in common]
    true = [str(gt_series[s]) for s in common]

    return {
        "ari": adjusted_rand_score(true, pred),
        "nmi": normalized_mutual_info_score(true, pred),
        "n_spots": len(common),
    }


def pairwise_ari(labels_list: list[dict]) -> list[float]:
    """Compute pairwise ARI between all pairs of label dicts."""
    aris = []
    for a, b in combinations(range(len(labels_list)), 2):
        da, db = labels_list[a], labels_list[b]
        common = sorted(set(da.keys()) & set(db.keys()))
        if len(common) < 10:
            continue
        la = [da[k] for k in common]
        lb = [db[k] for k in common]
        aris.append(adjusted_rand_score(la, lb))
    return aris


def bootstrap_ci(values, n_boot=10_000, seed=42):
    rng = np.random.default_rng(seed)
    arr = np.array(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    n = len(arr)
    if n == 0:
        return (np.nan, np.nan, np.nan)
    boots = np.array([arr[rng.choice(n, n, replace=True)].mean() for _ in range(n_boot)])
    return (arr.mean(), np.percentile(boots, 2.5), np.percentile(boots, 97.5))


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyze():
    if not RAW_PATH.exists():
        if committed_outputs_available():
            print(f"Raw checkpoint not found: {RAW_PATH}")
            print("Using committed aggregate outputs:")
            print(f"  {RESULTS_CSV}")
            print(f"  {PER_TRIAL_CSV}")
            print(f"  {SUMMARY_PATH}")
            print("To recompute raw-level metrics, run scripts/dlpfc_benchmark.py first.")
            return
        raise FileNotFoundError(
            f"{RAW_PATH} is required to recompute DLPFC benchmark metrics. "
            "Run scripts/dlpfc_benchmark.py first, or use the committed aggregate CSV/TXT outputs."
        )

    records = load_raw()
    print(f"Loaded {len(records)} trial records")

    # Load ground truth for all samples
    gt_data = {sid: load_ground_truth(sid) for sid in SAMPLE_IDS}

    # Group by (system, sample_id)
    grouped = defaultdict(list)
    for r in records:
        grouped[(r["system"], r["sample_id"])].append(r)

    results = []
    summary_lines = []
    summary_lines.append("=" * 70)
    summary_lines.append("DLPFC GROUND-TRUTH BENCHMARK ANALYSIS")
    summary_lines.append("=" * 70)

    # Per-system aggregation
    system_aris = defaultdict(list)
    system_nmis = defaultdict(list)
    system_cross_rep_aris = defaultdict(list)

    for (system, sample_id) in sorted(grouped.keys()):
        trials = grouped[(system, sample_id)]
        n_total = len(trials)
        n_success = sum(1 for t in trials if t["success"])
        success_rate = n_success / n_total if n_total > 0 else 0

        summary_lines.append(f"\n{system}/{sample_id}:")
        summary_lines.append(f"  Success: {n_success}/{n_total} ({success_rate*100:.0f}%)")

        # Extract labels and compute ground-truth metrics
        trial_aris = []
        trial_nmis = []
        trial_labels = []
        methods = []

        for t in trials:
            if not t["success"]:
                continue

            labels = extract_labels(system, sample_id, t["rep"])
            if labels is None:
                summary_lines.append(f"    rep {t['rep']}: labels not extractable")
                continue

            trial_labels.append(labels)

            # Ground truth comparison
            gt_metrics = compute_gt_metrics(labels, gt_data[sample_id])
            trial_aris.append(gt_metrics["ari"])
            trial_nmis.append(gt_metrics["nmi"])

            # Track methods
            m = t.get("method", "")
            if m:
                methods.append(m)

        # Ground-truth ARI
        if trial_aris:
            mean_ari, lo_ari, hi_ari = bootstrap_ci(trial_aris)
            summary_lines.append(f"  GT ARI: {mean_ari:.3f} [{lo_ari:.3f}, {hi_ari:.3f}] (n={len(trial_aris)})")
            system_aris[system].extend(trial_aris)
        else:
            mean_ari = lo_ari = hi_ari = np.nan
            summary_lines.append(f"  GT ARI: no data")

        # Ground-truth NMI
        if trial_nmis:
            mean_nmi, lo_nmi, hi_nmi = bootstrap_ci(trial_nmis)
            summary_lines.append(f"  GT NMI: {mean_nmi:.3f} [{lo_nmi:.3f}, {hi_nmi:.3f}]")
            system_nmis[system].extend(trial_nmis)
        else:
            mean_nmi = lo_nmi = hi_nmi = np.nan

        # Cross-replicate ARI
        cross_rep_pairs = pairwise_ari(trial_labels) if len(trial_labels) >= 2 else []
        if cross_rep_pairs:
            cr_mean, cr_lo, cr_hi = bootstrap_ci(cross_rep_pairs)
            summary_lines.append(f"  Cross-rep ARI: {cr_mean:.3f} [{cr_lo:.3f}, {cr_hi:.3f}] ({len(cross_rep_pairs)} pairs)")
            system_cross_rep_aris[system].extend(cross_rep_pairs)
        else:
            cr_mean = cr_lo = cr_hi = np.nan

        # Methods used
        if methods:
            from collections import Counter
            method_counts = Counter(methods)
            summary_lines.append(f"  Methods: {dict(method_counts)}")

        # Failures
        failures = [t.get("error", "unknown")[:60] for t in trials if not t["success"]]
        if failures:
            summary_lines.append(f"  Failures: {failures[:5]}")

        results.append({
            "system": system,
            "sample_id": sample_id,
            "n_total": n_total,
            "n_success": n_success,
            "success_rate": success_rate,
            "gt_ari_mean": mean_ari,
            "gt_ari_ci_lo": lo_ari,
            "gt_ari_ci_hi": hi_ari,
            "gt_nmi_mean": mean_nmi,
            "gt_nmi_ci_lo": lo_nmi,
            "gt_nmi_ci_hi": hi_nmi,
            "cross_rep_ari_mean": cr_mean,
            "cross_rep_ari_ci_lo": cr_lo,
            "cross_rep_ari_ci_hi": cr_hi,
            "n_cross_rep_pairs": len(cross_rep_pairs),
        })

    # Pooled per-system summary
    summary_lines.append("\n" + "=" * 70)
    summary_lines.append("POOLED SUMMARY (across all samples)")
    summary_lines.append("=" * 70)

    for system in sorted(system_aris.keys()):
        aris = system_aris[system]
        nmis = system_nmis[system]
        cr_aris = system_cross_rep_aris[system]

        mean_ari, lo_ari, hi_ari = bootstrap_ci(aris)
        mean_nmi, lo_nmi, hi_nmi = bootstrap_ci(nmis)
        cr_mean, cr_lo, cr_hi = bootstrap_ci(cr_aris) if cr_aris else (np.nan, np.nan, np.nan)

        n_total = sum(len(grouped[(system, sid)]) for sid in SAMPLE_IDS)
        n_success = sum(1 for sid in SAMPLE_IDS for t in grouped.get((system, sid), []) if t["success"])

        summary_lines.append(f"\n{system}:")
        summary_lines.append(f"  Success: {n_success}/{n_total}")
        summary_lines.append(f"  GT ARI: {mean_ari:.3f} [{lo_ari:.3f}, {hi_ari:.3f}] (n={len(aris)})")
        summary_lines.append(f"  GT NMI: {mean_nmi:.3f} [{lo_nmi:.3f}, {hi_nmi:.3f}]")
        summary_lines.append(f"  Cross-rep ARI: {cr_mean:.3f} [{cr_lo:.3f}, {cr_hi:.3f}] ({len(cr_aris)} pairs)")

    # Also save per-trial detail for figure scripts
    per_trial_rows = []
    for (system, sample_id) in sorted(grouped.keys()):
        trials = grouped[(system, sample_id)]
        for t in trials:
            if not t["success"]:
                per_trial_rows.append({
                    "system": system, "sample_id": sample_id, "rep": t["rep"],
                    "success": False, "ari": np.nan, "nmi": np.nan,
                    "method": t.get("method", ""),
                })
                continue
            labels = extract_labels(system, sample_id, t["rep"])
            if labels is None:
                per_trial_rows.append({
                    "system": system, "sample_id": sample_id, "rep": t["rep"],
                    "success": True, "ari": np.nan, "nmi": np.nan,
                    "method": t.get("method", ""),
                })
                continue
            gt_metrics = compute_gt_metrics(labels, gt_data[sample_id])
            per_trial_rows.append({
                "system": system, "sample_id": sample_id, "rep": t["rep"],
                "success": True,
                "ari": gt_metrics["ari"], "nmi": gt_metrics["nmi"],
                "n_spots": gt_metrics["n_spots"],
                "method": t.get("method", ""),
            })

    pd.DataFrame(per_trial_rows).to_csv(PER_TRIAL_CSV, index=False)
    print(f"Wrote {PER_TRIAL_CSV}")

    # Write results CSV
    if results:
        with open(RESULTS_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(results)
        print(f"Wrote {RESULTS_CSV}")

    # Write summary
    summary_text = "\n".join(summary_lines)
    with open(SUMMARY_PATH, "w") as f:
        f.write(summary_text + "\n")
    print(f"Wrote {SUMMARY_PATH}")
    print()
    print(summary_text)


if __name__ == "__main__":
    analyze()
