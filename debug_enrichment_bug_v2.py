"""Debug script v2 - Fix numpy recarray indexing"""
import scanpy as sc
import numpy as np

# Load the data
print("=" * 80)
print("Understanding numpy recarray structure")
print("=" * 80)

adata = sc.read_h5ad("chatspatial/data/card_example/card_spatial.h5ad")

# Create random clustering
np.random.seed(42)
adata.obs['leiden'] = np.random.choice(['0', '1', '2', '3', '4', '5', '6'], size=adata.n_obs)

# Run rank_genes_groups
sc.tl.rank_genes_groups(
    adata,
    groupby='leiden',
    method='wilcoxon',
    n_genes=100
)

result = adata.uns['rank_genes_groups']
names = result['names']
pvals_adj = result['pvals_adj']

print(f"\nRecarray structure:")
print(f"  Type: {type(names)}")
print(f"  Shape: {names.shape}")
print(f"  Dtype: {names.dtype}")
print(f"  Field names (groups): {names.dtype.names}")
print(f"  Number of groups: {len(names.dtype.names)}")

print(f"\nWhat does names[0] mean?")
print(f"  names[0] = {names[0]}")
print(f"  Type: {type(names[0])}")
print(f"  Length: {len(names[0])}")
print(f"  Interpretation: First row, containing top gene for each group")

print(f"\nWhat does names['0'] mean?")
group_0_genes = names['0']  # Get all genes for group '0'
print(f"  names['0'] shape: {group_0_genes.shape}")
print(f"  names['0'] length: {len(group_0_genes)}")
print(f"  First 10 genes for group '0': {group_0_genes[:10]}")

print(f"\n" + "=" * 80)
print("CURRENT BUGGY LOGIC:")
print("=" * 80)

print("\nCode: for i in range(len(names[0])):")
print(f"  This loops {len(names[0])} times (number of groups)")
print(f"  names[0][i] gives the i-th group's top gene")

degs_buggy = []
pvalue_threshold = 0.05

for i in range(len(names[0])):
    # Access p-value for position 0, field i
    # But pvals_adj[0] is a record with fields for each group!
    try:
        pval = pvals_adj[0][i]  # This gets the i-th field of the first record
        gene = names[0][i]
        degs_buggy.append(gene)
        print(f"  i={i}: gene={gene}, pval={pval}")
    except Exception as e:
        print(f"  i={i}: Error - {e}")

print(f"\nBuggy result: {len(degs_buggy)} genes selected")

print(f"\n" + "=" * 80)
print("CORRECT LOGIC - Get all genes for first group:")
print("=" * 80)

degs_correct = []
group_name = names.dtype.names[0]  # Get first group name

print(f"\nGetting genes for group '{group_name}':")
for i in range(len(names)):
    pval = pvals_adj[group_name][i]
    gene = names[group_name][i]
    if pval < pvalue_threshold:
        degs_correct.append(gene)

print(f"  Total genes with p < {pvalue_threshold}: {len(degs_correct)}")
print(f"  First 10: {degs_correct[:10]}")

print(f"\n" + "=" * 80)
print("ALTERNATIVE CORRECT LOGIC - Get genes from all groups:")
print("=" * 80)

degs_all_groups = []
for group_name in names.dtype.names:
    print(f"\nGroup '{group_name}':")
    group_genes = []
    for i in range(len(names)):
        pval = pvals_adj[group_name][i]
        gene = names[group_name][i]
        if pval < pvalue_threshold:
            group_genes.append(gene)
    print(f"  Genes with p < {pvalue_threshold}: {len(group_genes)}")
    degs_all_groups.extend(group_genes)

# Remove duplicates
degs_all_groups = list(set(degs_all_groups))
print(f"\nTotal unique genes from all groups: {len(degs_all_groups)}")

print(f"\n" + "=" * 80)
print("SUMMARY:")
print("=" * 80)
print(f"Buggy logic: {len(degs_buggy)} genes (equal to number of groups)")
print(f"Correct logic (first group): {len(degs_correct)} genes")
print(f"Correct logic (all groups): {len(degs_all_groups)} genes")
