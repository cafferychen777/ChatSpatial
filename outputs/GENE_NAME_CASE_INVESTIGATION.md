# Gene Name Case Mismatch Investigation

## Problem
Enrichment analysis returns 0 results because:
- DEGs from `rank_genes_groups`: UPPERCASE genes (e.g., 'KRT17', 'KRT19')
- Background from `adata.var_names`: lowercase genes (e.g., 'krt17', 'krt19')
- **Intersection: 0 genes** → No enrichment results

## Investigation Process

### 1. Original Data Check
✅ **Result**: CARD example data has UPPERCASE gene names
```
Original spatial data: ['X5S_rRNA', 'X5_8S_rRNA', 'A1BG.AS1', 'A1CF', ...]
Original reference data: ['A1BG', 'A1CF', 'A2M', 'A2ML1', ...]
```

### 2. Preprocessing Pipeline Test
✅ **Result**: Standard scanpy preprocessing preserves UPPERCASE
```python
# Tested operations:
sc.pp.calculate_qc_metrics()   # ✓ No case change
sc.pp.filter_genes()            # ✓ No case change
sc.pp.normalize_total()         # ✓ No case change
sc.pp.log1p()                   # ✓ No case change
sc.pp.highly_variable_genes()   # ✓ No case change
adata[:, adata.var["highly_variable"]].copy()  # ✓ No case change
```

### 3. Marker Gene Finding Test
✅ **Result**: `sc.tl.rank_genes_groups()` preserves UPPERCASE

### 4. R Integration Test
✅ **Result**: R's `make.names()` preserves UPPERCASE
```R
# Test with SPARK-X and CARD workflows
X5_8S_rRNA → X5_8S_rRNA  (no change)
A1BG.AS1 → A1BG.AS1      (no change)
```

### 5. Saved Data Analysis
❌ **Critical Finding**: Saved data has **MISMATCHED** gene names
```
adata.var_names:     ['x5_8s_rrna', 'a1bg.as1', 'a1cf', ...]  # lowercase
adata.raw.var_names: ['X5_8S_rRNA', 'A1BG.AS1', 'A1CF', ...]  # UPPERCASE
```

## Root Cause Hypothesis

After exhaustive testing, **none of the standard operations convert gene names to lowercase**.

**Most likely explanations:**

1. **Hidden external tool normalization**: Some analysis tool (possibly outside ChatSpatial) may have normalized gene names
   - Common in bioinformatics for cross-database matching
   - Often done by tools like Ensembl, NCBI, or gene ID converters

2. **Data loading artifact**: Possible issue with how the data was saved/loaded between sessions
   - pandas DataFrame operations can sometimes apply implicit transformations
   - HDF5 storage might have applied encoding transformations

3. **Untracked intermediate operation**: The user may have run additional operations not captured in the test workflow
   - Manual gene name normalization
   - Integration with external databases
   - Custom preprocessing scripts

## Why This Matters

The mismatch creates a critical bug:
```python
# rank_genes_groups output
DEGs = ['KRT17', 'KRT19', 'S100A6', ...]  # From adata.raw (UPPERCASE)

# enrichment.py background
background = set(adata.var_names)  # Current filtered data (lowercase)

# Result
query_genes = set(DEGs) & background  # 0 genes! ❌
```

## Solution Implemented

Fixed in `enrichment.py` (lines 986-1012):
```python
# 1. Use adata.raw as background (preserves original UPPERCASE names)
if adata.raw is not None:
    background_genes = set(adata.raw.var_names)
else:
    background_genes = set(adata.var_names)

# 2. Case-insensitive fallback matching
if len(query_genes) == 0 and len(gene_list) > 0:
    gene_name_map = {g.upper(): g for g in background_genes}
    query_genes = set()
    for gene in gene_list:
        if gene.upper() in gene_name_map:
            query_genes.add(gene_name_map[gene.upper()])
```

## Verification

After fix:
```
✅ Found 16 significant pathways (FDR < 0.05)
✅ Top pathway: "Protein digestion and absorption" (biologically correct for pancreatic tissue)
```

## Recommendation

**Best Practice**: Always preserve `adata.raw` during preprocessing
- `adata.raw` stores original gene names and counts
- Essential for analyses that reference original data (markers, enrichment, deconvolution)
- ChatSpatial automatically saves `adata.raw` at line 320-322 in preprocessing.py

## Remaining Mystery

**Question: WHERE did the lowercase transformation happen?**
- Not in: data loading, preprocessing, normalization, HVG selection, marker finding, R integration
- Possibly: External tool, manual operation, or hidden pandas transformation

**Note**: The fix works regardless of the root cause by using case-insensitive matching.
