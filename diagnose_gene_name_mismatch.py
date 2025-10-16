"""
Diagnose gene name mismatch between rank_genes_groups and adata.var_names
"""
import scanpy as sc

print("Loading data...")
adata = sc.read_h5ad("outputs/card_spatial_with_analysis.h5ad")

# Extract DEGs
result = adata.uns["rank_genes_groups"]
names = result["names"]
pvals = result.get("pvals_adj")

degs = []
for group_name in names.dtype.names:
    for i in range(len(names)):
        if pvals[group_name][i] < 0.05:
            gene = names[group_name][i]
            if gene not in degs:
                degs.append(gene)

print(f"\nExtracted {len(degs)} DEGs")
print(f"First 20 DEGs: {degs[:20]}")

# Check var_names
print(f"\nadata.var_names (first 20): {adata.var_names[:20].tolist()}")

# Check intersection
background_genes = set(adata.var_names)
query_genes_set = set(degs)

print(f"\n" + "=" * 80)
print("INTERSECTION CHECK")
print("=" * 80)
print(f"DEGs: {len(degs)}")
print(f"Background (var_names): {len(background_genes)}")
print(f"Intersection: {len(query_genes_set & background_genes)}")

# Check for specific examples
print(f"\n" + "=" * 80)
print("EXAMPLE GENE CHECKS")
print("=" * 80)

test_degs = degs[:10]
for gene in test_degs:
    in_var = gene in background_genes
    print(f"  {gene}: {'✓' if in_var else '✗'} in adata.var_names")

# Check if there's a case sensitivity issue
print(f"\n" + "=" * 80)
print("CASE SENSITIVITY CHECK")
print("=" * 80)

background_upper = {g.upper() for g in background_genes}
degs_upper = {g.upper() for g in degs}
overlap_upper = degs_upper & background_upper

print(f"Case-insensitive overlap: {len(overlap_upper)}")

if len(overlap_upper) > 0:
    print("\n✓ Found overlap with case-insensitive matching!")
    print("  This suggests case mismatch between DEG names and var_names")
else:
    print("\n✗ No overlap even with case-insensitive matching")

# Check if there's a filtering issue
print(f"\n" + "=" * 80)
print("GENE FILTERING CHECK")
print("=" * 80)

print(f"\nadata.var keys: {list(adata.var.columns)}")

if "highly_variable" in adata.var:
    hvg = adata.var_names[adata.var["highly_variable"]].tolist()
    print(f"Highly variable genes: {len(hvg)}")
    print(f"First 20 HVGs: {hvg[:20]}")

    # Check if DEGs are in HVG
    degs_in_hvg = set(degs) & set(hvg)
    print(f"\nDEGs in HVG: {len(degs_in_hvg)}")

# Check raw data
print(f"\n" + "=" * 80)
print("RAW DATA CHECK")
print("=" * 80)

if adata.raw is not None:
    print(f"✓ adata.raw exists")
    print(f"  Raw shape: {adata.raw.n_obs} × {adata.raw.n_vars}")
    print(f"  Raw var_names (first 20): {adata.raw.var_names[:20].tolist()}")

    raw_background = set(adata.raw.var_names)
    overlap_with_raw = query_genes_set & raw_background
    print(f"\n  DEGs ∩ raw.var_names: {len(overlap_with_raw)}")

    if len(overlap_with_raw) > 0:
        print("\n  ✓ DEGs match adata.raw.var_names!")
        print("  ❌ BUG: rank_genes_groups was run on raw data,")
        print("         but enrichment analysis uses filtered data!")
else:
    print("✗ adata.raw is None")

print(f"\n" + "=" * 80)
print("DIAGNOSIS")
print("=" * 80)

if len(overlap_with_raw) > 0:
    print("✓ ROOT CAUSE IDENTIFIED:")
    print("  rank_genes_groups references genes from adata.raw")
    print("  but enrichment analysis uses adata (filtered)")
    print("\n  SOLUTION: Use adata.raw.var_names as background,")
    print("           or re-run rank_genes_groups after filtering")
elif len(overlap_upper) > 0:
    print("✓ ROOT CAUSE IDENTIFIED:")
    print("  Gene names have case mismatch")
    print("\n  SOLUTION: Normalize gene names to same case")
else:
    print("❌ UNKNOWN ISSUE - genes don't match even with:")
    print("  - Raw data check")
    print("  - Case-insensitive check")
