"""Debug script to trace enrichment bug - why only 6 genes?"""
import scanpy as sc
import sys

# Load the data
print("=" * 80)
print("STEP 1: Load data")
print("=" * 80)
adata = sc.read_h5ad("chatspatial/data/card_example/card_spatial.h5ad")
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"Data type: {adata.X.dtype}")
print(f"Data contains negative values: {(adata.X < 0).sum() > 0 if hasattr(adata.X, 'sum') else 'sparse matrix'}")

# Check if data is log-normalized
if hasattr(adata.X, 'todense'):
    x_dense = adata.X.todense()
    has_negative = (x_dense < 0).any()
else:
    has_negative = (adata.X < 0).any()
print(f"Has negative values (log-normalized): {has_negative}")

print("\n" + "=" * 80)
print("STEP 2: Run find_markers simulation")
print("=" * 80)

# Simulate preprocessing (assign leiden clusters)
import numpy as np
if 'leiden' not in adata.obs:
    # Random clustering for testing
    np.random.seed(42)
    adata.obs['leiden'] = np.random.choice(['0', '1', '2', '3', '4', '5', '6'], size=adata.n_obs)
    print(f"Created random leiden clusters: {adata.obs['leiden'].value_counts().to_dict()}")
else:
    print(f"Existing leiden clusters: {adata.obs['leiden'].value_counts().to_dict()}")

# Run rank_genes_groups (simulating find_markers)
print("\nRunning sc.tl.rank_genes_groups...")
try:
    sc.tl.rank_genes_groups(
        adata,
        groupby='leiden',
        method='wilcoxon',
        n_genes=100
    )
    print("✓ rank_genes_groups completed successfully")
except Exception as e:
    print(f"❌ rank_genes_groups failed: {e}")
    sys.exit(1)

# Check results
print("\n" + "=" * 80)
print("STEP 3: Check rank_genes_groups results")
print("=" * 80)

if 'rank_genes_groups' in adata.uns:
    result = adata.uns['rank_genes_groups']
    print(f"✓ rank_genes_groups exists in adata.uns")
    print(f"  Keys: {list(result.keys())}")

    names = result['names']
    print(f"  Gene names type: {type(names)}")
    print(f"  Gene names shape: {names.shape}")
    print(f"  Number of genes per group: {len(names[0])}")

    pvals = result.get('pvals_adj', result.get('pvals'))
    if pvals is not None:
        print(f"  P-values available: YES")
        print(f"  P-values type: {type(pvals)}")
        print(f"  P-values shape: {pvals.shape}")

        # Count genes passing different thresholds
        first_group_pvals = pvals[0]
        n_sig_001 = (first_group_pvals < 0.01).sum()
        n_sig_005 = (first_group_pvals < 0.05).sum()
        n_sig_010 = (first_group_pvals < 0.10).sum()

        print(f"\n  Genes passing p-value thresholds (first group):")
        print(f"    p < 0.01: {n_sig_001}")
        print(f"    p < 0.05: {n_sig_005}")
        print(f"    p < 0.10: {n_sig_010}")
    else:
        print(f"  P-values available: NO")
else:
    print("❌ rank_genes_groups NOT found in adata.uns")

print("\n" + "=" * 80)
print("STEP 4: Simulate enrichment gene selection logic")
print("=" * 80)

# Simulate the exact logic from perform_ora
if "rank_genes_groups" in adata.uns:
    result = adata.uns["rank_genes_groups"]
    names = result["names"]

    pvals = None
    if "pvals_adj" in result:
        pvals = result["pvals_adj"]
    elif "pvals" in result:
        pvals = result["pvals"]

    degs = []
    pvalue_threshold = 0.05  # Default from perform_ora

    for i in range(len(names[0])):
        if pvals is not None:
            if pvals[0][i] < pvalue_threshold:
                degs.append(names[0][i])
        else:
            if i < 100:
                degs.append(names[0][i])

    print(f"Gene selection logic results:")
    print(f"  P-values available: {pvals is not None}")
    print(f"  P-value threshold: {pvalue_threshold}")
    print(f"  Number of genes selected: {len(degs)}")
    print(f"  First 10 genes: {degs[:10]}")

    # Now check background intersection
    background_genes = set(adata.var_names)
    query_genes = set(degs) & background_genes

    print(f"\n  Background genes: {len(background_genes)}")
    print(f"  Query genes (before intersection): {len(degs)}")
    print(f"  Query genes (after intersection): {len(query_genes)}")

    if len(query_genes) < len(degs):
        print(f"\n  ⚠️ Lost genes during intersection!")
        missing = set(degs) - background_genes
        print(f"  Missing genes (in degs but not in adata.var_names):")
        for gene in list(missing)[:10]:
            print(f"    - {gene}")

else:
    print("Falling back to variance-based selection...")
    import numpy as np

    if hasattr(adata.X, "todense"):
        var_scores = np.array(adata.X.todense().var(axis=0)).flatten()
    else:
        var_scores = np.array(adata.X.var(axis=0)).flatten()

    print(f"  Variance scores shape: {var_scores.shape}")
    print(f"  Variance scores range: {var_scores.min():.4f} to {var_scores.max():.4f}")

    top_indices = np.argsort(var_scores)[-500:]
    gene_list = adata.var_names[top_indices].tolist()

    print(f"  Top 500 variable genes selected: {len(gene_list)}")
    print(f"  First 10: {gene_list[:10]}")

print("\n" + "=" * 80)
print("DIAGNOSIS COMPLETE")
print("=" * 80)
