"""Debug script to check if rank_genes_groups is stored"""
import scanpy as sc
import numpy as np

# Load the data
adata = sc.read_h5ad("chatspatial/data/card_example/card_spatial.h5ad")

print("="*60)
print("CHECKING RANK_GENES_GROUPS STORAGE")
print("="*60)

# Check if rank_genes_groups exists
print(f"\n1. rank_genes_groups exists: {'rank_genes_groups' in adata.uns}")

if 'rank_genes_groups' in adata.uns:
    result = adata.uns['rank_genes_groups']
    print(f"\n2. rank_genes_groups structure:")
    print(f"   Keys: {list(result.keys())}")

    if 'names' in result:
        names = result['names']
        print(f"\n3. Gene names:")
        print(f"   Type: {type(names)}")
        print(f"   Shape: {names.shape if hasattr(names, 'shape') else 'N/A'}")
        print(f"   Groups: {names.dtype.names if hasattr(names.dtype, 'names') else 'N/A'}")

    if 'pvals_adj' in result or 'pvals' in result:
        pvals = result.get('pvals_adj', result.get('pvals'))
        print(f"\n4. P-values:")
        print(f"   Key used: {'pvals_adj' if 'pvals_adj' in result else 'pvals'}")
        print(f"   Type: {type(pvals)}")
        print(f"   Shape: {pvals.shape if hasattr(pvals, 'shape') else 'N/A'}")

        # Count significant genes
        if hasattr(pvals, 'shape') and len(pvals.shape) > 0:
            first_group_pvals = pvals[0] if pvals.ndim > 1 else pvals
            n_sig_001 = np.sum(first_group_pvals < 0.01)
            n_sig_005 = np.sum(first_group_pvals < 0.05)
            n_sig_010 = np.sum(first_group_pvals < 0.10)

            print(f"\n5. Significant genes (first group):")
            print(f"   p < 0.01: {n_sig_001}")
            print(f"   p < 0.05: {n_sig_005}")
            print(f"   p < 0.10: {n_sig_010}")

            # Show top 10 p-values
            print(f"\n6. Top 10 smallest p-values:")
            sorted_idx = np.argsort(first_group_pvals)[:10]
            for idx in sorted_idx:
                gene = names[0][idx] if hasattr(names, '__getitem__') else f"gene_{idx}"
                pval = first_group_pvals[idx]
                print(f"   {gene}: {pval:.2e}")
else:
    print("\n❌ rank_genes_groups NOT FOUND in adata.uns")
    print(f"   Available keys in adata.uns: {list(adata.uns.keys())}")
    print(f"\n   This means find_markers results were NOT persisted!")
