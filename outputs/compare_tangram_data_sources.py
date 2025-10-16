"""
Compare Tangram results with different data sources
Test 1: Using preprocessed data (scaled/normalized)
Test 2: Using raw data (counts)
"""

import scanpy as sc
import tangram as tg
import numpy as np

print("=" * 80)
print("Tangram Data Source Comparison Test")
print("=" * 80)

# Load preprocessed data
adata_sp_processed = sc.read_h5ad("/Users/apple/Research/SpatialTrans_MCP/chatspatial/outputs/tangram_test_data_2.h5ad")
adata_sc_original = sc.read_h5ad("/Users/apple/Research/SpatialTrans_MCP/chatspatial/chatspatial/data/card_example/card_reference.h5ad")

print("\n1. Data loaded:")
print(f"   Spatial (processed): {adata_sp_processed.shape}")
print(f"   Spatial has raw: {adata_sp_processed.raw is not None}")
print(f"   Reference: {adata_sc_original.shape}")

# Check genes BEFORE filtering
print(f"\n   Checking gene names before filtering:")
print(f"   Spatial genes sample: {list(adata_sp_processed.var_names[:5])}")
print(f"   Reference genes sample: {list(adata_sc_original.var_names[:5])}")

overlap_before = set(adata_sp_processed.var_names) & set(adata_sc_original.var_names)
print(f"   Overlapping genes before filtering: {len(overlap_before)}")

if len(overlap_before) == 0:
    print("\n   ❌ ERROR: No overlapping genes even before filtering!")
    print("   This suggests gene name format mismatch (e.g., case sensitivity)")
    # Check if case-insensitive matching helps
    sp_genes_upper = {g.upper() for g in adata_sp_processed.var_names}
    sc_genes_upper = {g.upper() for g in adata_sc_original.var_names}
    overlap_upper = sp_genes_upper & sc_genes_upper
    print(f"   Case-insensitive overlapping genes: {len(overlap_upper)}")
    if len(overlap_upper) > 0:
        print(f"   ✓ Gene name case mismatch detected! Fixing...")
        # Convert spatial genes to match reference case
        # This is a temporary fix for testing
        exit(1)

# Preprocess reference data
sc.pp.filter_genes(adata_sc_original, min_cells=10)
if 'cellType' not in adata_sc_original.obs:
    print("   ❌ ERROR: cellType not found in reference")
    exit(1)

# ========================================
# TEST 1: Using PREPROCESSED data (current ChatSpatial approach)
# ========================================
print("\n" + "=" * 80)
print("TEST 1: Using PREPROCESSED (scaled) data")
print("=" * 80)

adata_sp_test1 = adata_sp_processed.copy()
adata_sc_test1 = adata_sc_original.copy()

# Find overlapping genes
overlap_genes_test1 = list(set(adata_sc_test1.var_names) & set(adata_sp_test1.var_names))
print(f"\n2. Overlapping genes: {len(overlap_genes_test1)}")

if len(overlap_genes_test1) < 100:
    print(f"   ❌ Too few overlapping genes")
    exit(1)

# Subset to overlapping genes
adata_sc_test1 = adata_sc_test1[:, overlap_genes_test1]
adata_sp_test1 = adata_sp_test1[:, overlap_genes_test1]

print(f"   Spatial X type: {type(adata_sp_test1.X)}")
print(f"   Spatial X range: [{np.min(adata_sp_test1.X):.3f}, {np.max(adata_sp_test1.X):.3f}]")

# Run Tangram preprocessing
print("\n3. Running tg.pp_adatas...")
tg.pp_adatas(adata_sc_test1, adata_sp_test1, genes=None)

# Run Tangram mapping
print("\n4. Running Tangram mapping (clusters mode, 100 epochs)...")
ad_map_test1 = tg.map_cells_to_space(
    adata_sc_test1,
    adata_sp_test1,
    mode="clusters",
    cluster_label="cellType",
    device="cpu",
    num_epochs=100,  # Reduced for speed
)

# Project annotations
print("\n5. Projecting cell annotations...")
tg.project_cell_annotations(ad_map_test1, adata_sp_test1, annotation="cellType")

# Check results
if "tangram_ct_pred" in adata_sp_test1.obsm:
    ct_pred_test1 = adata_sp_test1.obsm["tangram_ct_pred"]
    print(f"\n6. Results (PREPROCESSED data):")
    print(f"   Shape: {ct_pred_test1.shape}")
    print(f"   Value range: [{ct_pred_test1.values.min():.6f}, {ct_pred_test1.values.max():.6f}]")
    print(f"   Contains NaN: {ct_pred_test1.isna().any().any()}")
    print(f"   All NaN: {ct_pred_test1.isna().all().all()}")
    if not ct_pred_test1.isna().all().all():
        print(f"   ✓ Valid predictions!")
        print(f"   First row sum: {ct_pred_test1.iloc[0].sum():.6f}")
    else:
        print(f"   ❌ ALL NaN - This is the bug!")
else:
    print(f"\n6. ❌ tangram_ct_pred not found")

# ========================================
# TEST 2: Using RAW (counts) data
# ========================================
print("\n" + "=" * 80)
print("TEST 2: Using RAW (counts) data")
print("=" * 80)

if adata_sp_processed.raw is None:
    print("\n❌ No raw data available in spatial dataset")
    print("   This test cannot be performed")
else:
    adata_sp_test2 = adata_sp_processed.raw.to_adata()
    adata_sp_test2.obsm['spatial'] = adata_sp_processed.obsm['spatial'].copy()
    adata_sc_test2 = adata_sc_original.copy()

    # Find overlapping genes
    overlap_genes_test2 = list(set(adata_sc_test2.var_names) & set(adata_sp_test2.var_names))
    print(f"\n2. Overlapping genes: {len(overlap_genes_test2)}")

    # Subset to overlapping genes
    adata_sc_test2 = adata_sc_test2[:, overlap_genes_test2]
    adata_sp_test2 = adata_sp_test2[:, overlap_genes_test2]

    print(f"   Spatial X type: {type(adata_sp_test2.X)}")
    if hasattr(adata_sp_test2.X, 'toarray'):
        X_dense = adata_sp_test2.X.toarray()
    else:
        X_dense = adata_sp_test2.X
    print(f"   Spatial X range: [{np.min(X_dense):.3f}, {np.max(X_dense):.3f}]")

    # Run Tangram preprocessing
    print("\n3. Running tg.pp_adatas...")
    tg.pp_adatas(adata_sc_test2, adata_sp_test2, genes=None)

    # Run Tangram mapping
    print("\n4. Running Tangram mapping (clusters mode, 100 epochs)...")
    ad_map_test2 = tg.map_cells_to_space(
        adata_sc_test2,
        adata_sp_test2,
        mode="clusters",
        cluster_label="cellType",
        device="cpu",
        num_epochs=100,
    )

    # Project annotations
    print("\n5. Projecting cell annotations...")
    tg.project_cell_annotations(ad_map_test2, adata_sp_test2, annotation="cellType")

    # Check results
    if "tangram_ct_pred" in adata_sp_test2.obsm:
        ct_pred_test2 = adata_sp_test2.obsm["tangram_ct_pred"]
        print(f"\n6. Results (RAW data):")
        print(f"   Shape: {ct_pred_test2.shape}")
        print(f"   Value range: [{ct_pred_test2.values.min():.6f}, {ct_pred_test2.values.max():.6f}]")
        print(f"   Contains NaN: {ct_pred_test2.isna().any().any()}")
        print(f"   All NaN: {ct_pred_test2.isna().all().all()}")
        if not ct_pred_test2.isna().all().all():
            print(f"   ✓ Valid predictions!")
            print(f"   First row sum: {ct_pred_test2.iloc[0].sum():.6f}")
        else:
            print(f"   ❌ ALL NaN - Raw data also fails!")
    else:
        print(f"\n6. ❌ tangram_ct_pred not found")

print("\n" + "=" * 80)
print("Comparison Complete")
print("=" * 80)
