"""
Diagnose Tangram clusters mode issue
Check what's in adata_sp after Tangram annotation
"""

import scanpy as sc

# Load the annotated spatial data
adata_sp = sc.read_h5ad("/Users/apple/Research/SpatialTrans_MCP/chatspatial/outputs/tangram_test_data_2.h5ad")

print("=" * 80)
print("Tangram Clusters Mode Diagnosis")
print("=" * 80)

print("\n1. adata_sp.obsm keys (should contain tangram_ct_pred):")
print(f"   Available keys: {list(adata_sp.obsm.keys())}")
print(f"   Has 'tangram_ct_pred': {'tangram_ct_pred' in adata_sp.obsm}")

if "tangram_ct_pred" in adata_sp.obsm:
    print("\n2. tangram_ct_pred content:")
    ct_pred = adata_sp.obsm["tangram_ct_pred"]
    print(f"   Shape: {ct_pred.shape}")
    print(f"   Columns (cell types): {list(ct_pred.columns)[:5]}... (total: {len(ct_pred.columns)})")
    print(f"   Value range: [{ct_pred.values.min():.3f}, {ct_pred.values.max():.3f}]")
    print(f"   First row sum: {ct_pred.iloc[0].sum():.3f}")
    print(f"   All row sums: [{ct_pred.sum(axis=1).min():.3f}, {ct_pred.sum(axis=1).max():.3f}]")

    print("\n3. First 3 spots' predictions:")
    print(ct_pred.head(3))

    print("\n4. Max probability per spot (should be used for cell type assignment):")
    max_probs = ct_pred.max(axis=1)
    print(f"   Range: [{max_probs.min():.6f}, {max_probs.max():.6f}]")
    print(f"   Mean: {max_probs.mean():.6f}")

    print("\n5. Assigned cell types (argmax):")
    assigned_types = ct_pred.idxmax(axis=1)
    print(f"   First 10 assignments: {list(assigned_types[:10])}")
    print(f"   Unique types assigned: {assigned_types.nunique()}")
    print(f"   Value counts:\n{assigned_types.value_counts()}")

else:
    print("\n❌ ERROR: tangram_ct_pred NOT FOUND in adata_sp.obsm")
    print("   This explains why counts and confidence_scores are empty!")

print("\n6. adata_sp.obs columns (should contain cell_type_tangram):")
print(f"   Available columns: {list(adata_sp.obs.columns)}")
print(f"   Has 'cell_type_tangram': {'cell_type_tangram' in adata_sp.obs}")

if "cell_type_tangram" in adata_sp.obs:
    print("\n7. cell_type_tangram content:")
    cell_types = adata_sp.obs["cell_type_tangram"]
    print(f"   Unique values: {cell_types.nunique()}")
    print(f"   Value counts:\n{cell_types.value_counts()}")
else:
    print("\n❌ ERROR: cell_type_tangram NOT FOUND in adata_sp.obs")

print("\n8. adata_sp.uns keys (may contain Tangram metadata):")
print(f"   Available keys: {list(adata_sp.uns.keys())}")

print("\n" + "=" * 80)
print("Diagnosis Complete")
print("=" * 80)
