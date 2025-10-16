"""
Deep dive investigation into enrichment analysis zero results problem.

This script will trace the entire ORA analysis flow to identify the exact failure point.
"""
import scanpy as sc
import gseapy as gp
import numpy as np
from scipy import stats

# Load data
print("=" * 80)
print("STEP 1: Load spatial data and run preprocessing")
print("=" * 80)

adata = sc.read_h5ad("chatspatial/data/card_example/card_spatial.h5ad")
print(f"\nData shape: {adata.n_obs} × {adata.n_vars}")
print(f"Gene names (first 10): {adata.var_names[:10].tolist()}")

# Preprocess
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Clustering
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.leiden(adata, resolution=1.0)

print(f"\nAfter preprocessing: {adata.n_obs} × {adata.n_vars}")
print(f"N clusters: {len(adata.obs['leiden'].unique())}")

# Run rank_genes_groups
print("\n" + "=" * 80)
print("STEP 2: Run rank_genes_groups")
print("=" * 80)

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", n_genes=100)

print("\nrank_genes_groups keys:", list(adata.uns['rank_genes_groups'].keys()))

# Extract DEGs using the same logic as enrichment.py
print("\n" + "=" * 80)
print("STEP 3: Extract DEGs (using same logic as enrichment.py)")
print("=" * 80)

result = adata.uns["rank_genes_groups"]
names = result["names"]
print(f"\nnames type: {type(names)}")
print(f"names dtype: {names.dtype}")
print(f"names shape: {names.shape}")
print(f"names.dtype.names (groups): {names.dtype.names}")

# Extract using corrected logic
pvals = result.get("pvals_adj", result.get("pvals"))
pvalue_threshold = 0.05

degs = []
for group_name in names.dtype.names:
    for i in range(len(names)):
        if pvals is not None:
            if pvals[group_name][i] < pvalue_threshold:
                gene = names[group_name][i]
                if gene not in degs:
                    degs.append(gene)

print(f"\n✓ Extracted {len(degs)} DEGs (p < {pvalue_threshold})")
print(f"  First 20 DEGs: {degs[:20]}")

# Check query genes
print("\n" + "=" * 80)
print("STEP 4: Check query genes vs background")
print("=" * 80)

background_genes = set(adata.var_names)
query_genes = set(degs) & background_genes

print(f"\nBackground genes: {len(background_genes)}")
print(f"DEG list: {len(degs)}")
print(f"Query genes (intersection): {len(query_genes)}")
print(f"  First 20 query genes: {list(query_genes)[:20]}")

if len(query_genes) == 0:
    print("\n❌ CRITICAL: Query genes is EMPTY!")
    print("   This means no DEGs overlap with adata.var_names")
else:
    print(f"\n✓ Query genes: {len(query_genes)} genes")

# Load KEGG gene sets
print("\n" + "=" * 80)
print("STEP 5: Load KEGG gene sets")
print("=" * 80)

organism = "human"
gene_sets_dict = gp.get_library("KEGG_2021_Human", organism=organism)

print(f"\nLoaded {len(gene_sets_dict)} KEGG pathways")
print(f"First 5 pathway names: {list(gene_sets_dict.keys())[:5]}")

# Check gene format in gene sets
first_pathway = list(gene_sets_dict.keys())[0]
first_genes = gene_sets_dict[first_pathway]
print(f"\nFirst pathway: {first_pathway}")
print(f"  N genes: {len(first_genes)}")
print(f"  First 10 genes: {first_genes[:10]}")

# Check gene format matching
print("\n" + "=" * 80)
print("STEP 6: Check gene format matching")
print("=" * 80)

sample_query_genes = list(query_genes)[:10] if len(query_genes) > 0 else []
sample_pathway_genes = first_genes[:10]

print(f"\nSample query genes: {sample_query_genes}")
print(f"Sample pathway genes: {sample_pathway_genes}")

# Check case sensitivity
if len(query_genes) > 0:
    query_upper = {g.upper() for g in query_genes}
    pathway_upper = {g.upper() for g in first_genes}
    overlap_upper = query_upper & pathway_upper
    print(f"\nCase-insensitive overlap (first pathway): {len(overlap_upper)} genes")

# Check for dot vs hyphen (e.g., MT.CO1 vs MT-CO1)
if len(query_genes) > 0:
    query_normalized = {g.replace('.', '-') for g in query_genes}
    pathway_normalized = {g.replace('.', '-') for g in first_genes}
    overlap_normalized = query_normalized & pathway_normalized
    print(f"Dot→Hyphen normalized overlap (first pathway): {len(overlap_normalized)} genes")

# Perform ORA for one pathway
print("\n" + "=" * 80)
print("STEP 7: Perform ORA for one pathway (manual)")
print("=" * 80)

if len(query_genes) > 0:
    min_size = 10
    max_size = 500

    # Filter gene sets by size
    filtered_sets = {}
    for name, genes in gene_sets_dict.items():
        if min_size <= len(genes) <= max_size:
            filtered_sets[name] = genes

    print(f"\nFiltered to {len(filtered_sets)} pathways (size {min_size}-{max_size})")

    # Test first pathway
    test_pathway = list(filtered_sets.keys())[0]
    test_genes = set(filtered_sets[test_pathway]) & background_genes

    print(f"\nTesting pathway: {test_pathway}")
    print(f"  Pathway genes: {len(filtered_sets[test_pathway])}")
    print(f"  Pathway genes in background: {len(test_genes)}")

    # Hypergeometric test
    a = len(query_genes & test_genes)
    b = len(query_genes - test_genes)
    c = len(test_genes - query_genes)
    d = len(background_genes - query_genes - test_genes)

    print(f"\nContingency table:")
    print(f"  a (query ∩ pathway): {a}")
    print(f"  b (query - pathway): {b}")
    print(f"  c (pathway - query): {c}")
    print(f"  d (neither): {d}")

    odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

    print(f"\nFisher's exact test:")
    print(f"  Odds ratio: {odds_ratio:.4f}")
    print(f"  P-value: {p_value:.4e}")

    if a == 0:
        print("\n❌ CRITICAL: No overlap between query genes and pathway genes!")
    else:
        print(f"\n✓ Found {a} overlapping genes")
        overlapping = list(query_genes & test_genes)[:10]
        print(f"  Overlapping genes (first 10): {overlapping}")

    # Test all pathways
    print("\n" + "=" * 80)
    print("STEP 8: Test all pathways")
    print("=" * 80)

    significant_pathways = []
    total_tested = 0
    max_overlap = 0

    for pathway_name, pathway_genes in filtered_sets.items():
        gs_genes_set = set(pathway_genes) & background_genes

        a = len(query_genes & gs_genes_set)
        max_overlap = max(max_overlap, a)

        if a > 0:
            b = len(query_genes - gs_genes_set)
            c = len(gs_genes_set - query_genes)
            d = len(background_genes - query_genes - gs_genes_set)

            odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

            if p_value < 0.05:
                significant_pathways.append((pathway_name, a, p_value))

            total_tested += 1

    print(f"\nTested {total_tested} pathways")
    print(f"Max overlap found: {max_overlap} genes")
    print(f"Significant pathways (p < 0.05, before FDR): {len(significant_pathways)}")

    if len(significant_pathways) > 0:
        print("\nTop 10 significant pathways:")
        for name, overlap, pval in sorted(significant_pathways, key=lambda x: x[2])[:10]:
            print(f"  {name}: {overlap} genes, p={pval:.4e}")
    else:
        print("\n❌ NO significant pathways found!")
        print("\nDiagnosis:")
        if max_overlap == 0:
            print("  - NO overlap between any pathway and query genes")
            print("  - Likely gene ID format mismatch")
        elif max_overlap < 3:
            print(f"  - Max overlap only {max_overlap} genes (too few for significance)")
            print("  - Query gene list may be too small or unrelated to pathways")
else:
    print("\n❌ Cannot perform ORA - query_genes is empty!")

print("\n" + "=" * 80)
print("DIAGNOSIS COMPLETE")
print("=" * 80)
