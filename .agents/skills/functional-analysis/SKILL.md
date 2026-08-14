---
name: functional-analysis
description: |
  Understand the biological meaning of gene sets through pathway and functional enrichment analysis.
  Use when user has gene lists (DEGs, markers, SVGs) and wants to know what pathways or functions they represent.
  Triggers: "pathway analysis", "enrichment", "GSEA", "GO terms", "what pathways",
  "biological function", "gene set", "KEGG", "Reactome", "functional annotation".
context: fork
agent: general-purpose
allowed-tools: Read, Grep, Glob
---

# Functional Analysis

## Overview

This skill answers: **What biological processes/pathways are represented by this gene set?**

Two main approaches:
1. **Over-representation Analysis (ORA)**: Test if genes are enriched in pathways
2. **Gene Set Enrichment Analysis (GSEA)**: Use full ranked gene list

Plus specialized analysis:
3. **CNV Analysis**: Copy number variation in spatial context (tumor analysis)

## Decision Tree: Which Method?

```
Q: What type of gene list do you have?
│
├─ Discrete gene list (DEGs, markers)
│   └─ Over-Representation Analysis (ORA)
│       ├─ Input: Gene list (up/down regulated)
│       ├─ Tests: Hypergeometric / Fisher's exact
│       └─ Output: Enriched pathways with p-values
│
├─ Ranked gene list (all genes with scores)
│   └─ Gene Set Enrichment Analysis (GSEA)
│       ├─ Input: All genes ranked by fold change/statistic
│       ├─ Tests: Running sum statistic
│       └─ Output: Enrichment scores, leading edge genes
│
└─ Spatial tumor data
    └─ CNV Analysis
        ├─ Infer copy number from expression
        ├─ Map CNV spatially
        └─ Identify subclones
```

## Over-Representation Analysis (ORA)

### When to Use

- You have a discrete gene list (e.g., significant DEGs)
- Want to know what pathways are enriched
- Quick, interpretable results needed

### Databases

| Database | Content | Best For |
|----------|---------|----------|
| GO (Gene Ontology) | BP, MF, CC terms | Broad functional annotation |
| KEGG | Metabolic/signaling pathways | Pathway diagrams |
| Reactome | Curated pathways | Detailed mechanisms |
| MSigDB | Hallmarks, curated sets | Cancer, immunology |
| WikiPathways | Community pathways | Specific processes |

### Workflow

#### Step 1: Prepare Gene List

Requirements:
- Gene symbols or IDs (consistent naming)
- Typically 50-500 genes works best
- Separate up/down regulated if directional

#### Step 2: Run ORA

Use `analyze_enrichment` tool with:
- `gene_list`: Your genes of interest
- `database`: "GO_BP", "KEGG", "Reactome", etc.
- `background`: All detected genes (important!)

#### Step 3: Interpret Results

**Key outputs:**
- `term`: Pathway/GO term name
- `pvalue` / `padj`: Significance
- `odds_ratio`: Enrichment strength
- `gene_count`: Genes in overlap
- `genes`: Specific overlapping genes

**Significance criteria:**
```
Typically significant:
- padj < 0.05
- gene_count >= 3 (robust overlap)
- odds_ratio > 2 (meaningful enrichment)
```

## Gene Set Enrichment Analysis (GSEA)

### When to Use

- You have scores for all genes (not just significant ones)
- Want to detect subtle but coordinated changes
- Directional enrichment matters (up vs down)

### Advantages over ORA

- Uses information from all genes
- Detects subtle coordinated changes
- No arbitrary threshold for "significant" genes
- Better statistical properties

### Workflow

#### Step 1: Rank Genes

Rank all genes by:
- Log fold change (signed)
- -log10(p) × sign(logFC) (recommended)
- Test statistic

#### Step 2: Run GSEA

Use `analyze_enrichment` tool with:
- `method`: "gsea"
- `ranking`: Pre-computed gene ranking
- `database`: Gene set database

#### Step 3: Interpret Results

**Key outputs:**
- `NES`: Normalized Enrichment Score
  - Positive: Enriched in highly ranked genes
  - Negative: Enriched in lowly ranked genes
- `pvalue` / `padj`: Significance
- `leading_edge`: Core genes driving enrichment

**Enrichment plot interpretation:**
- Running enrichment score curve
- Peaks show where gene set is concentrated
- Leading edge = genes before peak

## CNV Analysis (Tumor Spatial Data)

### Purpose

Infer copy number variations from spatial transcriptomics for:
- Tumor heterogeneity mapping
- Clone identification
- Spatial organization of subclones

### Methods

| Method | Approach | Best For |
|--------|----------|----------|
| inferCNV | Reference-based | scRNA-seq, validated |
| CopyKAT | Reference-free | No normal reference |
| SCEVAN | Fast, integrated | Large datasets |

### Workflow

#### Step 1: Define Reference

For CNV inference, need:
- Normal cells as reference (from annotation)
- Or use reference-free method

#### Step 2: Run CNV Analysis

Use `analyze_cnv` tool with:
- `method`: "infercnv", "copykat", or "scevan"
- `reference_group`: Normal cell label (if applicable)

#### Step 3: Interpret Results

**CNV outputs:**
- `cnv_score`: Per-cell CNV intensity
- `clone_id`: Inferred subclone assignment
- Chromosome-level gains/losses

**Spatial visualization:**
- Map CNV scores spatially
- Identify clone spatial organization
- Correlate with histology

## Visualization

### Enrichment Visualizations

| Plot Type | Use Case |
|-----------|----------|
| Bar plot | Top enriched terms |
| Dot plot | Significance + gene count |
| Network | Term relationships |
| GSEA plot | Enrichment curve |
| Heatmap | Gene-pathway matrix |

Use `visualize_data` with `plot_type="enrichment"`.

### CNV Visualizations

| Plot Type | Use Case |
|-----------|----------|
| Heatmap | CNV across chromosomes |
| Spatial | CNV score on tissue |
| Clone map | Subclone distribution |

## Common Analysis Scenarios

### Scenario 1: Understand DEG Function
```
Goal: What pathways are activated in disease?
1. Get DEG list from differential analysis
2. Split into up/down regulated
3. Run ORA on each separately
4. Compare enriched pathways
```

### Scenario 2: Pathway-Level Changes
```
Goal: Are cancer hallmarks affected?
1. Rank all genes by fold change
2. Run GSEA with MSigDB Hallmarks
3. Identify significantly enriched hallmarks
4. Examine leading edge genes
```

### Scenario 3: Tumor Clone Mapping
```
Goal: Map tumor heterogeneity spatially
1. Annotate normal cells (if present)
2. Run CNV inference
3. Cluster by CNV profile to identify clones
4. Map clones spatially
5. Correlate with histology
```

### Scenario 4: Immune Signature Analysis
```
Goal: What immune processes are active where?
1. Define immune gene signatures
2. Score each spot for signatures
3. Map signature scores spatially
4. Correlate with cell composition
```

## Best Practices

### Gene List Quality

- **Remove duplicates**: Each gene counted once
- **Consistent IDs**: Use symbols or Ensembl, not mixed
- **Appropriate background**: All expressed genes, not all genes
- **Reasonable size**: 50-500 genes typically best for ORA

### Multiple Testing

- Always use adjusted p-values
- Report number of tests performed
- Consider redundancy in GO terms

### Interpretation Caution

- Enrichment ≠ causation
- Multiple pathways may overlap
- Validate key findings experimentally

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| No significant terms | Gene list too small/diverse | Use GSEA, check gene IDs |
| Too many terms | Loose threshold | Tighten padj, filter redundant |
| Unexpected results | Wrong background | Use detected genes as background |
| Missing expected pathway | Gene naming mismatch | Check ID format (symbol vs Ensembl) |

## Output Summary

Successful functional analysis provides:
- **Enriched pathways**: Ranked by significance
- **Gene-pathway mapping**: Which genes in which pathways
- **Biological interpretation**: What processes are affected
- **CNV landscape**: Tumor heterogeneity map (if applicable)
- **Downstream hypotheses**: Mechanisms to investigate
