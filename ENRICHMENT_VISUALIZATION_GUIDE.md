# Enrichment Analysis Visualization Guide

## Overview

ChatSpatial separates analysis and visualization into distinct steps. After running enrichment analysis, you can visualize the results using the visualization tools.

## Workflow

### Step 1: Run Enrichment Analysis

```python
# Run GSEA analysis
analyze_enrichment(
    data_id="data_1",
    params={
        "method": "gsea",
        "pvalue_cutoff": 0.05,
        "plot_top_terms": 15,
        "gene_set_database": "GO_Biological_Process"
    }
)
```

This saves results to `adata.uns['gsea_results']` (or `ora_results` for ORA).

### Step 2: Visualize Results

```python
# Create visualization
create_visualization(
    data_id="data_1",
    params={
        "plot_type": "gsea",
        "gsea_plot_type": "barplot",  # Options: barplot, enrichment_plot, dotplot
        "n_top_pathways": 10,
        "title": "GSEA Results - Top Enriched Pathways"
    }
)
```

## Visualization Types

### 1. Barplot (Default)
- Shows top enriched pathways as horizontal bars
- Color-coded by enrichment score
- Best for presenting top results

### 2. Enrichment Plot
- Classic GSEA running enrichment score plot
- Shows enrichment profile for a specific pathway
- Use `feature` parameter to specify pathway

### 3. Dotplot
- Shows enrichment across multiple conditions
- Size represents significance
- Color represents enrichment score

## Storage Format

Results are stored in `adata.uns` with the following keys:
- `gsea_results`: DataFrame with columns: Term, ES, NES, NOM p-val, FDR q-val
- `ora_results`: DataFrame with columns: pathway, odds_ratio, pvalue, adjusted_pvalue

The visualization tools automatically detect and use these results.

## Example Usage in Claude

1. Run enrichment analysis:
   ```
   User: Run GSEA enrichment analysis on my data
   ```

2. Visualize results:
   ```
   User: Show me a barplot of the top 10 enriched pathways
   ```

The visualization will automatically use the stored results from step 1.