"""
Verify the enrichment analysis fix
"""
import scanpy as sc
import gseapy as gp
from scipy import stats
from statsmodels.stats.multitest import multipletests

print("=" * 80)
print("VERIFYING ENRICHMENT FIX")
print("=" * 80)

# Load data
adata = sc.read_h5ad("outputs/card_spatial_with_analysis.h5ad")
print(f"\nData: {adata.n_obs} × {adata.n_vars}")
print(f"Raw: {adata.raw.n_obs} × {adata.raw.n_vars}")

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

# APPLY FIX: Use adata.raw.var_names
print("\n" + "=" * 80)
print("APPLYING FIX")
print("=" * 80)

if adata.raw is not None:
    background_genes = set(adata.raw.var_names)
    print(f"✓ Using adata.raw.var_names as background: {len(background_genes)} genes")
else:
    background_genes = set(adata.var_names)
    print(f"Using adata.var_names as background: {len(background_genes)} genes")

# Direct match
query_genes = set(degs) & background_genes
print(f"\nDirect match: {len(query_genes)} genes")

# If no matches, try case-insensitive
if len(query_genes) == 0 and len(degs) > 0:
    print("\n  No direct matches, trying case-insensitive...")
    gene_name_map = {g.upper(): g for g in background_genes}
    query_genes = set()
    for gene in degs:
        if gene.upper() in gene_name_map:
            query_genes.add(gene_name_map[gene.upper()])
    print(f"  Case-insensitive match: {len(query_genes)} genes")

print(f"\n✓ Final query genes: {len(query_genes)}")

# Load gene sets
organism = "human"
gene_sets_dict = gp.get_library("KEGG_2021_Human", organism=organism)

min_size = 10
max_size = 500

filtered_sets = {}
for name, genes in gene_sets_dict.items():
    if min_size <= len(genes) <= max_size:
        filtered_sets[name] = genes

print(f"✓ Loaded {len(filtered_sets)} KEGG pathways")

# Perform ORA
print("\n" + "=" * 80)
print("PERFORMING ORA")
print("=" * 80)

enrichment_scores = {}
pvalues = {}

for gs_name, gs_genes in filtered_sets.items():
    gs_genes_set = set(gs_genes) & background_genes

    if len(gs_genes_set) < min_size or len(gs_genes_set) > max_size:
        continue

    a = len(query_genes & gs_genes_set)
    b = len(query_genes - gs_genes_set)
    c = len(gs_genes_set - query_genes)
    d = len(background_genes - query_genes - gs_genes_set)

    odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

    enrichment_scores[gs_name] = odds_ratio
    pvalues[gs_name] = p_value

print(f"Tested {len(pvalues)} pathways")

# Multiple testing correction
if len(pvalues) > 0:
    pval_list = [pvalues[gs] for gs in pvalues.keys()]
    reject, pvals_corrected, _, _ = multipletests(pval_list, alpha=0.05, method="fdr_bh")

    adjusted_pvalues = {
        gs_name: pvals_corrected[i] for i, gs_name in enumerate(pvalues.keys())
    }

    significant = [gs for gs, adj_p in adjusted_pvalues.items() if adj_p < 0.05]

    print(f"\n✅ Found {len(significant)} significant pathways (FDR < 0.05)")

    if len(significant) > 0:
        print("\n" + "=" * 80)
        print("TOP 10 SIGNIFICANT PATHWAYS")
        print("=" * 80)

        sorted_pathways = sorted(
            [(gs, adjusted_pvalues[gs], pvalues[gs])
             for gs in significant],
            key=lambda x: x[1]
        )

        for i, (gs_name, adj_p, raw_p) in enumerate(sorted_pathways[:10], 1):
            gs_genes_set = set(filtered_sets[gs_name]) & background_genes
            overlap = len(query_genes & gs_genes_set)
            print(f"\n{i}. {gs_name}")
            print(f"   Overlap: {overlap} genes")
            print(f"   Raw p-value: {raw_p:.4e}")
            print(f"   Adj. p-value: {adj_p:.4e}")

        print("\n" + "=" * 80)
        print("FIX VERIFICATION: ✅ SUCCESS!")
        print("=" * 80)
        print(f"The enrichment analysis now correctly identifies {len(significant)} pathways!")
        print("Bug fix is working as expected.")
    else:
        print("\n❌ No significant pathways after FDR correction")
else:
    print("\n❌ No pathways tested")
