#!/usr/bin/env python3
"""
DLPFC Ground-Truth Benchmark — Data Preparation

Downloads the Maynard et al. 2021 DLPFC dataset via spatialLIBD::fetch_data,
extracts three samples (one per donor), saves as h5ad via zellkonverter,
then post-processes in Python to remove ground-truth columns and preprocess.

Output structure:
  paper/data/dlpfc_benchmark/
    samples/{sample_id}.h5ad          — benchmark input (NO ground truth)
    ground_truth/{sample_id}_labels.csv — expert annotations
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR.parent / "data" / "dlpfc_benchmark"
SAMPLES_DIR = DATA_DIR / "samples"
GT_DIR = DATA_DIR / "ground_truth"

SAMPLES_DIR.mkdir(parents=True, exist_ok=True)
GT_DIR.mkdir(parents=True, exist_ok=True)

# Target samples: one per donor
TARGET_SAMPLES = ["151673", "151507", "151669"]


def build_r_script(export_dir: Path) -> str:
    """Generate R script that downloads DLPFC data and exports per-sample h5ad."""
    return f"""
suppressPackageStartupMessages(library(spatialLIBD))
suppressPackageStartupMessages(library(SpatialExperiment))
suppressPackageStartupMessages(library(zellkonverter))

export_dir <- "{export_dir}"
target_samples <- c("{TARGET_SAMPLES[0]}", "{TARGET_SAMPLES[1]}", "{TARGET_SAMPLES[2]}")

cat("Downloading full DLPFC dataset via spatialLIBD::fetch_data...\\n")
spe <- fetch_data("spe")

cat("Available samples:", paste(unique(spe$sample_id), collapse=", "), "\\n")
cat("Dimensions:", nrow(spe), "genes x", ncol(spe), "spots\\n")

for (sid in target_samples) {{
    cat("\\nProcessing sample", sid, "...\\n")
    sub <- spe[, spe$sample_id == sid]
    n_spots <- ncol(sub)
    cat("  Spots:", n_spots, "\\n")
    if (n_spots == 0) {{
        cat("  ERROR: No spots found\\n")
        next
    }}

    # Save ground truth labels as CSV
    gt <- data.frame(
        spot = colnames(sub),
        ground_truth = sub$layer_guess_reordered,
        stringsAsFactors = FALSE
    )
    write.csv(gt, file.path(export_dir, paste0(sid, "_gt.csv")),
              row.names = FALSE, quote = FALSE)

    # Convert to SingleCellExperiment for zellkonverter
    sce <- as(sub, "SingleCellExperiment")
    h5ad_path <- file.path(export_dir, paste0(sid, ".h5ad"))
    writeH5AD(sce, h5ad_path)
    cat("  Saved:", h5ad_path, "\\n")
}}
cat("\\nDone.\\n")
"""


def run_r_export(export_dir: Path):
    """Run R script to download and export DLPFC data."""
    r_script_path = export_dir / "export_dlpfc.R"
    r_script_path.write_text(build_r_script(export_dir))

    cmd = ["R", "--no-save", "--no-restore", "-f", str(r_script_path)]
    print(f"Running R export...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)

    if result.returncode != 0:
        print("R STDOUT:", result.stdout[-3000:])
        print("R STDERR:", result.stderr[-2000:])
        raise RuntimeError(f"R export failed with code {result.returncode}")

    print(result.stdout[-2000:])


def postprocess_sample(export_dir: Path, sample_id: str):
    """Load R-exported h5ad, remove ground truth, preprocess, save."""
    print(f"\nPost-processing {sample_id}...")

    # Load h5ad from R export
    raw_path = export_dir / f"{sample_id}.h5ad"
    adata = ad.read_h5ad(raw_path)
    print(f"  Raw shape: {adata.shape}")
    print(f"  obs columns: {list(adata.obs.columns)}")

    # Save ground truth from CSV (more reliable than extracting from h5ad)
    gt_csv = export_dir / f"{sample_id}_gt.csv"
    gt_df = pd.read_csv(gt_csv)
    gt_df = gt_df.set_index("spot")
    gt_path = GT_DIR / f"{sample_id}_labels.csv"
    gt_df.to_csv(gt_path)

    n_labeled = gt_df["ground_truth"].notna().sum()
    print(f"  Ground truth: {n_labeled}/{len(gt_df)} labeled spots")
    print(f"  Label distribution:")
    for label, count in gt_df["ground_truth"].value_counts().items():
        print(f"    {label}: {count}")
    n_na = gt_df["ground_truth"].isna().sum()
    if n_na > 0:
        print(f"    NA: {n_na}")

    # Remove all annotation/leakage columns
    keep_cols = [c for c in adata.obs.columns
                 if c in ("barcode_id", "sample_id", "in_tissue",
                          "array_row", "array_col", "cell_count",
                          "expr_chrM", "expr_chrM_ratio")]
    adata.obs = adata.obs[keep_cols] if keep_cols else pd.DataFrame(index=adata.obs_names)
    print(f"  Cleaned obs columns: {list(adata.obs.columns)}")

    # Ensure we have raw counts in X
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    # Store counts layer
    adata.layers["counts"] = adata.X.copy()

    # Standard preprocessing
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                                layer="counts")
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

    # Ensure spatial coordinates exist
    if "spatial" not in adata.obsm:
        # Try to reconstruct from obs
        if "array_row" in adata.obs.columns and "array_col" in adata.obs.columns:
            adata.obsm["spatial"] = adata.obs[["array_col", "array_row"]].values.astype(float)
            print("  WARNING: spatial coords reconstructed from array_row/col")

    # Save
    h5ad_path = SAMPLES_DIR / f"{sample_id}.h5ad"
    adata.write_h5ad(h5ad_path)
    print(f"  Saved: {h5ad_path}")
    print(f"  Final shape: {adata.shape[0]} spots x {adata.shape[1]} genes")
    print(f"  Has spatial: {'spatial' in adata.obsm}")


def main():
    print("=" * 60)
    print("DLPFC Ground-Truth Benchmark — Data Preparation")
    print("=" * 60)

    # Check if data already exists
    existing = [s for s in TARGET_SAMPLES
                if (SAMPLES_DIR / f"{s}.h5ad").exists()
                and (GT_DIR / f"{s}_labels.csv").exists()]
    if len(existing) == len(TARGET_SAMPLES):
        print("All samples already prepared. Verifying...")
        for sid in TARGET_SAMPLES:
            adata = ad.read_h5ad(SAMPLES_DIR / f"{sid}.h5ad")
            gt = pd.read_csv(GT_DIR / f"{sid}_labels.csv", index_col=0)
            print(f"  {sid}: {adata.shape[0]} spots x {adata.shape[1]} genes, "
                  f"{gt['ground_truth'].notna().sum()} labeled, "
                  f"spatial={'spatial' in adata.obsm}")
            leak = [c for c in adata.obs.columns
                    if any(kw in c.lower() for kw in ["layer", "ground", "manual"])]
            if leak:
                print(f"    WARNING: leakage columns: {leak}")
        print("All OK. Delete files to re-download.")
        return

    # Step 1: R download and export
    with tempfile.TemporaryDirectory(prefix="dlpfc_") as tmpdir:
        export_dir = Path(tmpdir)
        print(f"\nExport directory: {export_dir}")

        run_r_export(export_dir)

        # Step 2: Post-process each sample
        for sid in TARGET_SAMPLES:
            h5ad_file = export_dir / f"{sid}.h5ad"
            if not h5ad_file.exists():
                print(f"WARNING: {sid}.h5ad not found, skipping")
                continue
            postprocess_sample(export_dir, sid)

    # Final verification
    print("\n" + "=" * 60)
    print("Verification")
    print("=" * 60)
    for sid in TARGET_SAMPLES:
        h5ad_path = SAMPLES_DIR / f"{sid}.h5ad"
        gt_path = GT_DIR / f"{sid}_labels.csv"
        if h5ad_path.exists() and gt_path.exists():
            adata = ad.read_h5ad(h5ad_path)
            gt = pd.read_csv(gt_path, index_col=0)
            print(f"  {sid}: {adata.shape[0]} spots x {adata.shape[1]} genes, "
                  f"{gt['ground_truth'].notna().sum()} labeled spots, "
                  f"spatial={'spatial' in adata.obsm}")
        else:
            print(f"  {sid}: MISSING")


if __name__ == "__main__":
    main()
