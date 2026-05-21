#!/usr/bin/env python3
"""
End-to-End Benchmark Analysis

Reads the raw JSONL from e2e_benchmark.py and computes:
  1. Execution success rates
  2. Method consistency across replicates
  3. Output concordance:
     - Spatial domains: pairwise ARI between replicates
     - SVG detection: pairwise Jaccard@100 between replicates
     - Cell communication: pairwise Jaccard@50 between ligand-receptor lists
"""

import csv
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR.parent / "data" / "e2e_benchmark"
OUTPUT_DIR = DATA_DIR / "outputs"
RAW_PATH = DATA_DIR / "e2e_benchmark_raw.jsonl"
RESULTS_CSV = DATA_DIR / "e2e_benchmark_results.csv"
SUMMARY_PATH = DATA_DIR / "e2e_benchmark_summary.txt"


def committed_outputs_available() -> bool:
    return all(path.exists() and path.stat().st_size > 0 for path in (
        RESULTS_CSV,
        SUMMARY_PATH,
    ))

# ---------------------------------------------------------------------------
# Load raw data
# ---------------------------------------------------------------------------

def load_raw():
    records = []
    with open(RAW_PATH) as f:
        for line in f:
            records.append(json.loads(line))
    return records


# ---------------------------------------------------------------------------
# Output extraction
# ---------------------------------------------------------------------------

def extract_spatial_domain_labels(system: str, task_id: str, rep: int) -> dict | None:
    """Extract cluster labels from a spatial domains trial.

    Returns dict mapping spot barcode/index → label string, to handle
    different spot counts across runs.
    """
    trial_dir = OUTPUT_DIR / system / task_id / f"rep_{rep}"

    def _extract_from_adata(adata):
        # Prefer columns with domain/leiden/cluster in the name
        domain_cols = [c for c in adata.obs.columns if any(kw in c.lower() for kw in ["domain", "leiden", "cluster", "louvain"])]
        if domain_cols:
            return dict(zip(adata.obs_names, adata.obs[domain_cols[0]].astype(str)))
        return None

    # ChatSpatial saves h5ad with domain labels
    h5ad = trial_dir / "result.h5ad"
    if h5ad.exists():
        adata = ad.read_h5ad(h5ad)
        result = _extract_from_adata(adata)
        if result:
            return result

    # STAgent/SpatialAgent may have saved h5ad
    for f in sorted(trial_dir.glob("*.h5ad")):
        adata = ad.read_h5ad(f)
        result = _extract_from_adata(adata)
        if result:
            return result

    return None


def extract_svg_list(system: str, task_id: str, rep: int) -> list[str] | None:
    """Extract top spatially variable genes from an SVG trial."""
    trial_dir = OUTPUT_DIR / system / task_id / f"rep_{rep}"

    # Check CSV first
    csv_path = trial_dir / "svg_results.csv"
    if csv_path.exists():
        df = pd.read_csv(csv_path)
        gene_col = [c for c in df.columns if c.lower() in ("gene", "genes", "gene_name")]
        if gene_col:
            return df[gene_col[0]].tolist()[:100]
        # First column
        return df.iloc[:, 0].tolist()[:100]

    # Check h5ad
    for f in sorted(trial_dir.glob("*.h5ad")):
        adata = ad.read_h5ad(f)
        if "morans_I" in adata.var.columns:
            top = adata.var.sort_values("morans_I", ascending=False).head(100).index.tolist()
            return top
        if "highly_variable" in adata.var.columns:
            top = adata.var[adata.var["highly_variable"]].index.tolist()[:100]
            if top:
                return top

    # Check any CSV in the directory
    for f in sorted(trial_dir.glob("*.csv")):
        try:
            df = pd.read_csv(f)
            if len(df) >= 10 and len(df.columns) >= 1:
                return df.iloc[:, 0].tolist()[:100]
        except Exception:
            pass

    return None


def extract_lr_pairs(system: str, task_id: str, rep: int) -> list[str] | None:
    """Extract top ligand-receptor pairs from a cell communication trial.

    Returns list of normalized "LIGAND^RECEPTOR" strings for Jaccard comparison.
    """
    trial_dir = OUTPUT_DIR / system / task_id / f"rep_{rep}"

    for f in sorted(trial_dir.glob("*.csv")):
        try:
            df = pd.read_csv(f)
            cols_lower = {c.lower(): c for c in df.columns}
            if "ligand" in cols_lower and "receptor" in cols_lower:
                pairs = [
                    f"{row[cols_lower['ligand']]}^{row[cols_lower['receptor']]}".upper()
                    for _, row in df.head(50).iterrows()
                ]
                return pairs
            # CellPhoneDB format
            if "interacting_pair" in cols_lower:
                return [
                    p.upper()
                    for p in df[cols_lower["interacting_pair"]].head(50).tolist()
                ]
        except Exception:
            pass

    return None


# ---------------------------------------------------------------------------
# Concordance metrics
# ---------------------------------------------------------------------------

def pairwise_ari(labels_dict_list: list[dict]) -> list[float]:
    """Compute pairwise ARI between all pairs of label dictionaries.

    Each element is {spot_index: label} to handle different spot counts.
    ARI is computed only on the intersection of spots.
    """
    aris = []
    for a, b in combinations(range(len(labels_dict_list)), 2):
        da, db = labels_dict_list[a], labels_dict_list[b]
        common = sorted(set(da.keys()) & set(db.keys()))
        if len(common) < 10:
            continue
        la = [da[k] for k in common]
        lb = [db[k] for k in common]
        ari = adjusted_rand_score(la, lb)
        aris.append(ari)
    return aris


def pairwise_jaccard(gene_lists: list[list[str]], k: int = 100) -> list[float]:
    """Compute pairwise Jaccard similarity between top-k gene lists."""
    jaccards = []
    for a, b in combinations(range(len(gene_lists)), 2):
        set_a = set(gene_lists[a][:k])
        set_b = set(gene_lists[b][:k])
        if not set_a and not set_b:
            jaccards.append(1.0)
        elif not set_a or not set_b:
            jaccards.append(0.0)
        else:
            jaccards.append(len(set_a & set_b) / len(set_a | set_b))
    return jaccards


def bootstrap_ci(values, n_boot=10_000, seed=42):
    rng = np.random.default_rng(seed)
    arr = np.array(values, dtype=float)
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
            print(f"  {SUMMARY_PATH}")
            print("To recompute raw-level metrics, run scripts/e2e_benchmark.py first.")
            return
        raise FileNotFoundError(
            f"{RAW_PATH} is required to recompute end-to-end benchmark metrics. "
            "Run scripts/e2e_benchmark.py first, or use the committed aggregate CSV/TXT outputs."
        )

    records = load_raw()

    # Group by (system, task)
    grouped = defaultdict(list)
    for r in records:
        grouped[(r["system"], r["task_id"])].append(r)

    results = []
    summary_lines = []

    summary_lines.append("=" * 70)
    summary_lines.append("END-TO-END BENCHMARK ANALYSIS")
    summary_lines.append("=" * 70)

    for (system, task_id) in sorted(grouped.keys()):
        trials = grouped[(system, task_id)]
        n_total = len(trials)
        n_success = sum(1 for t in trials if t["success"])
        success_rate = n_success / n_total if n_total > 0 else 0

        # Wall time for successful trials
        times = [t["wall_time"] for t in trials if t["success"]]
        avg_time = np.mean(times) if times else 0

        # Failure taxonomy
        failures = [t.get("error", "unknown")[:60] for t in trials if not t["success"]]

        summary_lines.append(f"\n{system}/{task_id}:")
        summary_lines.append(f"  Success: {n_success}/{n_total} ({success_rate*100:.0f}%)")
        summary_lines.append(f"  Avg time (success): {avg_time:.1f}s")
        if failures:
            summary_lines.append(f"  Failures: {failures}")

        # Output concordance
        concordance = np.nan
        concordance_pairs = []

        if task_id == "spatial_domains":
            labels = []
            for t in trials:
                if t["success"]:
                    l = extract_spatial_domain_labels(system, task_id, t["rep"])
                    if l is not None:
                        labels.append(l)
            if len(labels) >= 2:
                concordance_pairs = pairwise_ari(labels)
                mean, lo, hi = bootstrap_ci(concordance_pairs)
                concordance = mean
                summary_lines.append(f"  ARI concordance: {mean:.3f} [{lo:.3f}, {hi:.3f}] ({len(concordance_pairs)} pairs)")
            else:
                summary_lines.append(f"  ARI concordance: insufficient data ({len(labels)} outputs)")

        elif task_id == "svg_detection":
            gene_lists = []
            for t in trials:
                if t["success"]:
                    gl = extract_svg_list(system, task_id, t["rep"])
                    if gl is not None:
                        gene_lists.append(gl)
            if len(gene_lists) >= 2:
                concordance_pairs = pairwise_jaccard(gene_lists)
                mean, lo, hi = bootstrap_ci(concordance_pairs)
                concordance = mean
                summary_lines.append(f"  Jaccard@100 concordance: {mean:.3f} [{lo:.3f}, {hi:.3f}] ({len(concordance_pairs)} pairs)")
            else:
                summary_lines.append(f"  Jaccard@100 concordance: insufficient data ({len(gene_lists)} outputs)")

        elif task_id == "cell_communication":
            lr_lists = []
            for t in trials:
                if t["success"]:
                    lp = extract_lr_pairs(system, task_id, t["rep"])
                    if lp is not None:
                        lr_lists.append(lp)
            if len(lr_lists) >= 2:
                concordance_pairs = pairwise_jaccard(lr_lists, k=50)
                mean, lo, hi = bootstrap_ci(concordance_pairs)
                concordance = mean
                summary_lines.append(f"  Jaccard@50 concordance: {mean:.3f} [{lo:.3f}, {hi:.3f}] ({len(concordance_pairs)} pairs)")
            else:
                summary_lines.append(f"  Jaccard@50 concordance: insufficient data ({len(lr_lists)} outputs)")

        # Method consistency (for ChatSpatial, always the same; for others, check logs)
        methods = []
        for t in trials:
            if t["success"]:
                m = t.get("method", "")
                if m:
                    methods.append(m)
        method_consistency = 1.0 if len(set(methods)) <= 1 else len(set(methods)) / len(methods)

        results.append({
            "system": system,
            "task_id": task_id,
            "n_total": n_total,
            "n_success": n_success,
            "success_rate": success_rate,
            "avg_time": avg_time,
            "concordance": concordance,
            "concordance_n_pairs": len(concordance_pairs),
            "method_consistency": method_consistency,
        })

    # Write CSV
    if not results:
        raise ValueError("No benchmark result rows were produced from the raw checkpoint.")

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
