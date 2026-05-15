#!/usr/bin/env python3
"""
Compute effect sizes and confidence intervals for case study figures.

Outputs:
  - paper/data/Case2_OSCC/effect_sizes_enrichment.csv
  - paper/data/Case2_OSCC/consensus_binomial_ci.csv
  - paper/data/Case2_OSCC/moranI_cross_sample_stats.csv
  - paper/data/Case1_HGSOC/cliffs_delta_cnv_vs_background.csv
  - paper/data/Case1_HGSOC/cnv_score_per_patient_stats.csv
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

# ── Paths ────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
PAPER_DIR = SCRIPT_DIR.parent
CODE_DIR = PAPER_DIR.parent / "code"

OSCC_H5AD_DIR = (
    CODE_DIR
    / "data"
    / "geo_data"
    / "analyzed_projects"
    / "GSE208253"
    / "processed_h5ad"
)
OSCC_RESULTS_DIR = (
    CODE_DIR
    / "data"
    / "geo_data"
    / "analyzed_projects"
    / "GSE208253"
    / "server_results"
)
OSCC_OUT_DIR = PAPER_DIR / "data" / "Case2_OSCC"

HGSOC_H5AD_DIR = PAPER_DIR / "data" / "Case1_HGSOC" / "Visium_replication"
HGSOC_OUT_DIR = PAPER_DIR / "data" / "Case1_HGSOC"

SEED = 42
N_BOOT = 10_000


# ── Utilities ────────────────────────────────────────────────────

def rank_biserial_r(x: np.ndarray, y: np.ndarray) -> float:
    """Rank-biserial correlation from Mann-Whitney U."""
    u_stat, _ = stats.mannwhitneyu(x, y, alternative="two-sided")
    n1, n2 = len(x), len(y)
    return 1 - 2 * u_stat / (n1 * n2)


def bootstrap_ci(
    x: np.ndarray,
    y: np.ndarray,
    stat_func,
    n_boot: int = N_BOOT,
    seed: int = SEED,
) -> tuple[float, float]:
    """Bootstrap 95% CI (percentile method) for a two-sample statistic."""
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        bx = rng.choice(x, len(x), replace=True)
        by = rng.choice(y, len(y), replace=True)
        boots[i] = stat_func(bx, by)
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def clopper_pearson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Exact binomial (Clopper-Pearson) confidence interval."""
    if k == 0:
        lo = 0.0
    else:
        lo = stats.beta.ppf(alpha / 2, k, n - k + 1)
    if k == n:
        hi = 1.0
    else:
        hi = stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    return float(lo), float(hi)


def bootstrap_median_ci(
    x: np.ndarray,
    n_boot: int = N_BOOT,
    seed: int = SEED,
) -> tuple[float, float]:
    """Bootstrap 95% CI for the median (one-sample, percentile method)."""
    rng = np.random.default_rng(seed)
    medians = np.empty(n_boot)
    for i in range(n_boot):
        medians[i] = np.median(rng.choice(x, len(x), replace=True))
    return float(np.percentile(medians, 2.5)), float(np.percentile(medians, 97.5))


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Cliff's delta effect size for two independent groups."""
    n1, n2 = len(x), len(y)
    # Vectorized: count concordant and discordant pairs
    more = np.sum(x[:, None] > y[None, :])
    less = np.sum(x[:, None] < y[None, :])
    return float((more - less) / (n1 * n2))


# ═══════════════════════════════════════════════════════════════
# 1. OSCC: Rank-biserial r + CIs for TC vs LE enrichment
# ═══════════════════════════════════════════════════════════════

def compute_oscc_effect_sizes():
    print("Computing OSCC enrichment effect sizes...")

    # Load existing enrichment summary for p-values and log2fc
    enrich_df = pd.read_csv(OSCC_RESULTS_DIR / "per_sample_tc_le_enrichment.csv")

    rows = []
    for sample_id in sorted(enrich_df["sample"].unique()):
        h5ad_path = OSCC_H5AD_DIR / f"{sample_id}_processed.h5ad"
        if not h5ad_path.exists():
            print(f"  WARNING: {h5ad_path} not found, skipping")
            continue

        adata = ad.read_h5ad(h5ad_path)
        tc_mask = adata.obs["tc_le_region"] == "TC"
        le_mask = adata.obs["tc_le_region"] == "LE"
        n_tc = int(tc_mask.sum())
        n_le = int(le_mask.sum())

        if n_tc < 5 or n_le < 5:
            continue

        deconv_cols = [c for c in adata.obs.columns if c.startswith("deconvolution_card_")]

        for col in deconv_cols:
            cell_type = col.replace("deconvolution_card_", "")
            tc_vals = adata.obs.loc[tc_mask, col].values.astype(float)
            le_vals = adata.obs.loc[le_mask, col].values.astype(float)

            # Skip if all zeros
            if tc_vals.sum() == 0 and le_vals.sum() == 0:
                rows.append({
                    "sample": sample_id,
                    "cell_type": cell_type,
                    "n_TC": n_tc,
                    "n_LE": n_le,
                    "rank_biserial_r": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                })
                continue

            r = rank_biserial_r(tc_vals, le_vals)
            ci_lo, ci_hi = bootstrap_ci(tc_vals, le_vals, rank_biserial_r)

            rows.append({
                "sample": sample_id,
                "cell_type": cell_type,
                "n_TC": n_tc,
                "n_LE": n_le,
                "rank_biserial_r": round(r, 4),
                "ci_lower": round(ci_lo, 4),
                "ci_upper": round(ci_hi, 4),
            })

        print(f"  {sample_id}: {n_tc} TC, {n_le} LE spots")

    effect_df = pd.DataFrame(rows)

    # Merge with existing enrichment data for log2fc and FDR
    merged = effect_df.merge(
        enrich_df[["sample", "cell_type", "log2fc", "pvalue_fdr"]],
        on=["sample", "cell_type"],
        how="left",
    )

    out_path = OSCC_OUT_DIR / "effect_sizes_enrichment.csv"
    merged.to_csv(out_path, index=False)
    print(f"  Wrote {out_path} ({len(merged)} rows)")
    return merged


# ═══════════════════════════════════════════════════════════════
# 2. OSCC: Clopper-Pearson CIs for consensus proportions
# ═══════════════════════════════════════════════════════════════

def compute_oscc_consensus_ci():
    print("Computing OSCC consensus binomial CIs...")
    consensus = pd.read_csv(OSCC_RESULTS_DIR / "consensus_tc_le_patterns.csv", index_col=0)

    rows = []
    for cell_type in consensus.index:
        n = int(consensus.loc[cell_type, "total_samples"])

        for direction in ["le", "tc"]:
            k = int(consensus.loc[cell_type, f"{direction}_enriched_count"])
            lo, hi = clopper_pearson_ci(k, n)
            rows.append({
                "cell_type": cell_type,
                "direction": direction.upper(),
                "k": k,
                "n": n,
                "proportion": round(k / n, 4),
                "ci_lower": round(lo, 4),
                "ci_upper": round(hi, 4),
            })

    ci_df = pd.DataFrame(rows)
    out_path = OSCC_OUT_DIR / "consensus_binomial_ci.csv"
    ci_df.to_csv(out_path, index=False)
    print(f"  Wrote {out_path} ({len(ci_df)} rows)")
    return ci_df


# ═══════════════════════════════════════════════════════════════
# 3. OSCC: Sample-level Wilcoxon signed-rank tests
# ═══════════════════════════════════════════════════════════════

def compute_oscc_sample_level_tests():
    """Sample-level inference: Wilcoxon signed-rank on per-patient log2FC.

    For each cell type, collects sample-level log2FCs (s6 excluded: no LE),
    runs Wilcoxon signed-rank (H0: median log2FC = 0) and sign test,
    computes median log2FC with bootstrap 95% CI.
    """
    from statsmodels.stats.multitest import multipletests

    print("Computing OSCC sample-level tests (Wilcoxon signed-rank)...")

    enrich_df = pd.read_csv(OSCC_RESULTS_DIR / "per_sample_tc_le_enrichment.csv")

    # Exclude s6 (no LE spots — le_mean is NaN for all cell types)
    enrich_df = enrich_df[enrich_df["sample"] != "s6"].copy()

    cell_types = sorted(enrich_df["cell_type"].unique())
    rows = []

    for ct in cell_types:
        ct_data = enrich_df[enrich_df["cell_type"] == ct]
        log2fcs = ct_data["log2fc"].dropna().values.astype(float)
        n_samples = len(log2fcs)

        # Nonzero values for Wilcoxon (it excludes exact zeros)
        nonzero = log2fcs[log2fcs != 0.0]
        n_nonzero = len(nonzero)

        # Wilcoxon signed-rank test
        if n_nonzero >= 3:
            w_stat, w_p = stats.wilcoxon(nonzero, alternative="two-sided")
        else:
            w_stat, w_p = np.nan, np.nan

        # Sign test
        n_pos = int(np.sum(nonzero > 0))
        n_neg = int(np.sum(nonzero < 0))
        if n_nonzero >= 1:
            sign_result = stats.binomtest(n_pos, n_nonzero, 0.5)
            sign_p = sign_result.pvalue
        else:
            sign_p = np.nan

        # Median and bootstrap CI
        if n_samples >= 2:
            med = float(np.median(log2fcs))
            ci_lo, ci_hi = bootstrap_median_ci(log2fcs)
        elif n_samples == 1:
            med = float(log2fcs[0])
            ci_lo, ci_hi = np.nan, np.nan
        else:
            med, ci_lo, ci_hi = np.nan, np.nan, np.nan

        direction = "LE" if med > 0 else ("TC" if med < 0 else "--")
        if np.isnan(med):
            direction = "--"

        rows.append({
            "cell_type": ct,
            "n_samples": n_samples,
            "n_nonzero": n_nonzero,
            "median_log2fc": round(med, 4) if not np.isnan(med) else np.nan,
            "median_ci_lower": round(ci_lo, 4) if not np.isnan(ci_lo) else np.nan,
            "median_ci_upper": round(ci_hi, 4) if not np.isnan(ci_hi) else np.nan,
            "wilcoxon_stat": round(w_stat, 1) if not np.isnan(w_stat) else np.nan,
            "wilcoxon_p": w_p,
            "sign_n_pos": n_pos,
            "sign_n_neg": n_neg,
            "sign_p": sign_p,
            "direction": direction,
        })

    result_df = pd.DataFrame(rows)

    # BH-FDR across Wilcoxon p-values
    valid = result_df["wilcoxon_p"].notna()
    result_df["wilcoxon_p_fdr"] = np.nan
    if valid.sum() > 0:
        _, fdr, _, _ = multipletests(
            result_df.loc[valid, "wilcoxon_p"].values, method="fdr_bh"
        )
        result_df.loc[valid, "wilcoxon_p_fdr"] = fdr

    out_path = OSCC_OUT_DIR / "sample_level_tests.csv"
    result_df.to_csv(out_path, index=False)
    print(f"  Wrote {out_path} ({len(result_df)} rows)")
    print(result_df.to_string(index=False))
    return result_df


# ═══════════════════════════════════════════════════════════════
# 4. OSCC: Effect size summary (all samples, not just significant)
# ═══════════════════════════════════════════════════════════════

def compute_oscc_effect_size_summary(effect_df: pd.DataFrame):
    """Summarize rank-biserial |r| across ALL samples (primary)
    and significant-sample subset (secondary) per cell type."""
    print("Computing OSCC effect size summary...")

    consensus_path = OSCC_RESULTS_DIR / "consensus_tc_le_patterns.csv"
    consensus = pd.read_csv(consensus_path, index_col=0)

    rows = []
    for ct in sorted(effect_df["cell_type"].unique()):
        ct_data = effect_df[effect_df["cell_type"] == ct]

        # All samples with non-NaN r
        valid = ct_data["rank_biserial_r"].dropna()
        n_all = len(valid)
        mean_abs_r_all = float(valid.abs().mean()) if n_all > 0 else np.nan

        # Significant samples only
        sig = ct_data[ct_data["pvalue_fdr"] < 0.05]["rank_biserial_r"].dropna()
        n_sig = len(sig)
        mean_abs_r_sig = float(sig.abs().mean()) if n_sig > 0 else np.nan

        # Direction and consensus from the consensus file
        if ct in consensus.index:
            le_pct = consensus.loc[ct, "le_percentage"]
            tc_pct = consensus.loc[ct, "tc_percentage"]
            if le_pct > tc_pct:
                direction = "LE"
                consensus_pct = le_pct
            elif tc_pct > le_pct:
                direction = "TC"
                consensus_pct = tc_pct
            else:
                direction = "--"
                consensus_pct = 0
        else:
            direction, consensus_pct = "--", 0

        rows.append({
            "cell_type": ct,
            "direction": direction,
            "n_all": n_all,
            "mean_abs_r_all": round(mean_abs_r_all, 4) if not np.isnan(mean_abs_r_all) else np.nan,
            "n_significant": n_sig,
            "mean_abs_r_sig": round(mean_abs_r_sig, 4) if not np.isnan(mean_abs_r_sig) else np.nan,
            "consensus_pct": round(consensus_pct, 1),
        })

    summary_df = pd.DataFrame(rows)
    out_path = OSCC_OUT_DIR / "effect_size_summary.csv"
    summary_df.to_csv(out_path, index=False)
    print(f"  Wrote {out_path}")
    print(summary_df.to_string(index=False))
    return summary_df


# ═══════════════════════════════════════════════════════════════
# 5. OSCC: Cross-sample Moran's I statistics
# ═══════════════════════════════════════════════════════════════

def compute_oscc_moranI_stats():
    print("Computing OSCC cross-sample Moran's I stats...")
    moranI = pd.read_csv(OSCC_RESULTS_DIR / "moranI_results_all_samples.csv")

    targets = ["FN1", "COL1A1", "SPRR1B", "KRT6C", "SPRR2E", "LAMC2", "ITGA5", "CLDN4"]

    rows = []
    for gene in targets:
        sub = moranI[moranI["gene"] == gene]
        n = len(sub)
        vals = sub["I"].values

        if n < 2:
            rows.append({
                "gene": gene,
                "n_samples": n,
                "mean_I": round(vals[0], 4) if n == 1 else np.nan,
                "sd_I": np.nan,
                "ci_lower": np.nan,
                "ci_upper": np.nan,
            })
            continue

        mean_I = vals.mean()
        sd_I = vals.std(ddof=1)
        # t-distribution CI for the mean
        se = sd_I / np.sqrt(n)
        t_crit = stats.t.ppf(0.975, df=n - 1)
        ci_lo = mean_I - t_crit * se
        ci_hi = mean_I + t_crit * se

        rows.append({
            "gene": gene,
            "n_samples": n,
            "mean_I": round(mean_I, 4),
            "sd_I": round(sd_I, 4),
            "ci_lower": round(ci_lo, 4),
            "ci_upper": round(ci_hi, 4),
        })

    stats_df = pd.DataFrame(rows)
    out_path = OSCC_OUT_DIR / "moranI_cross_sample_stats.csv"
    stats_df.to_csv(out_path, index=False)
    print(f"  Wrote {out_path}")
    print(stats_df.to_string(index=False))
    return stats_df


# ═══════════════════════════════════════════════════════════════
# 4. HGSOC: Cliff's delta for CNV markers vs background
# ═══════════════════════════════════════════════════════════════

def compute_hgsoc_cliffs_delta():
    print("Computing HGSOC Cliff's delta (CNV markers vs background)...")

    svg_path = HGSOC_OUT_DIR / "svg_results_P3.csv"
    svg_df = pd.read_csv(svg_path, index_col=0)

    # Load P3 h5ad to get CNV cluster markers
    h5ad_path = HGSOC_H5AD_DIR / "P3_cnv_infercnv.h5ad"
    adata = ad.read_h5ad(h5ad_path)

    # Identify CNV cluster marker genes (from cnv_leiden clusters)
    # Use differential expression between CNV clusters if available
    cnv_marker_genes = []
    if "rank_genes_groups" in adata.uns:
        rgg = adata.uns["rank_genes_groups"]
        names = rgg["names"]
        for group in names.dtype.names:
            cnv_marker_genes.extend(names[group][:20].tolist())
        cnv_marker_genes = list(set(cnv_marker_genes))

    if not cnv_marker_genes:
        # Fallback: use cnv_leiden to run quick DE
        print("  No rank_genes_groups found; computing DE for CNV clusters...")
        import scanpy as sc
        adata_sub = adata[adata.obs["cnv_leiden"].notna()].copy()
        if "cnv_leiden" in adata_sub.obs.columns:
            sc.tl.rank_genes_groups(adata_sub, groupby="cnv_leiden", method="wilcoxon", n_genes=20)
            rgg = adata_sub.uns["rank_genes_groups"]
            names = rgg["names"]
            for group in names.dtype.names:
                cnv_marker_genes.extend(names[group][:20].tolist())
            cnv_marker_genes = list(set(cnv_marker_genes))

    # Filter to genes present in SVG results
    cnv_in_svg = [g for g in cnv_marker_genes if g in svg_df.index]
    background = [g for g in svg_df.index if g not in cnv_in_svg]

    print(f"  CNV markers in SVG results: {len(cnv_in_svg)}")
    print(f"  Background genes: {len(background)}")

    if len(cnv_in_svg) < 5:
        print("  WARNING: Too few CNV markers for meaningful comparison")
        return None

    cnv_vals = svg_df.loc[cnv_in_svg, "I"].values.astype(float)
    bg_vals = svg_df.loc[background, "I"].values.astype(float)

    delta = cliffs_delta(cnv_vals, bg_vals)
    ci_lo, ci_hi = bootstrap_ci(cnv_vals, bg_vals, cliffs_delta)

    # Also compute Mann-Whitney for verification
    u_stat, p_val = stats.mannwhitneyu(cnv_vals, bg_vals, alternative="two-sided")

    result = {
        "n_cnv_markers": len(cnv_in_svg),
        "n_background": len(background),
        "cnv_mean_I": round(cnv_vals.mean(), 4),
        "bg_mean_I": round(bg_vals.mean(), 4),
        "cliffs_delta": round(delta, 4),
        "ci_lower": round(ci_lo, 4),
        "ci_upper": round(ci_hi, 4),
        "mannwhitney_U": round(u_stat, 1),
        "mannwhitney_p": p_val,
    }

    result_df = pd.DataFrame([result])
    out_path = HGSOC_OUT_DIR / "cliffs_delta_cnv_vs_background.csv"
    result_df.to_csv(out_path, index=False)
    print(f"  Wrote {out_path}")
    print(f"  Cliff's delta = {delta:.4f} [{ci_lo:.4f}, {ci_hi:.4f}], p = {p_val:.4g}")
    return result


# ═══════════════════════════════════════════════════════════════
# 5. HGSOC: Per-patient CNV score mean + SD
# ═══════════════════════════════════════════════════════════════

def compute_hgsoc_cnv_stats():
    print("Computing HGSOC per-patient CNV score statistics...")

    rows = []
    for p in range(1, 9):
        h5ad_path = HGSOC_H5AD_DIR / f"P{p}_cnv_infercnv.h5ad"
        if not h5ad_path.exists():
            continue

        adata = ad.read_h5ad(h5ad_path)
        if "cnv_score" not in adata.obs.columns:
            print(f"  WARNING: P{p} missing cnv_score column")
            continue

        scores = adata.obs["cnv_score"].dropna().values.astype(float)
        rows.append({
            "patient": f"P{p}",
            "n_spots": len(scores),
            "mean_cnv_score": round(scores.mean(), 6),
            "sd_cnv_score": round(scores.std(ddof=1), 6),
            "median_cnv_score": round(np.median(scores), 6),
            "iqr_lower": round(np.percentile(scores, 25), 6),
            "iqr_upper": round(np.percentile(scores, 75), 6),
        })
        print(f"  P{p}: {len(scores)} spots, mean={scores.mean():.4f}, SD={scores.std(ddof=1):.4f}")

    cnv_df = pd.DataFrame(rows)
    out_path = HGSOC_OUT_DIR / "cnv_score_per_patient_stats.csv"
    cnv_df.to_csv(out_path, index=False)
    print(f"  Wrote {out_path}")
    return cnv_df


# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    merged = compute_oscc_effect_sizes()
    print()
    compute_oscc_consensus_ci()
    print()
    compute_oscc_sample_level_tests()
    print()
    compute_oscc_effect_size_summary(merged)
    print()
    compute_oscc_moranI_stats()
    print()
    compute_hgsoc_cliffs_delta()
    print()
    compute_hgsoc_cnv_stats()
    print("\nDone.")
