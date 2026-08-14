---
name: cell-interaction
description: |
  Analyze cell-cell communication patterns and ligand-receptor interactions in spatial context.
  Use when user wants to understand how cells communicate, signaling pathways, or intercellular interactions.
  Triggers: "cell communication", "ligand receptor", "cell-cell interaction", "signaling",
  "how do cells talk", "communication analysis", "CellChat", "LIANA".
context: fork
agent: general-purpose
allowed-tools: Read, Grep, Glob
---

# Cell-Cell Interaction Analysis

## Overview

This skill answers: **How do cells communicate with each other in this tissue?**

Cell communication fundamentally requires:
1. A **sender cell** expressing a ligand
2. A **receiver cell** expressing a receptor
3. **Spatial proximity** enabling the interaction

## Prerequisites

Before running cell communication analysis:
- [ ] Cell type annotations available (from deconvolution or direct annotation)
- [ ] Species identified (human or mouse)
- [ ] Biological question defined (specific pathways? global patterns?)

## Method Selection

### Decision Framework

```
Q: What level of analysis do you need?
│
├─ Ligand-receptor pairs only
│   └─ LIANA+ (Python) - Multi-method consensus, recommended default
│
├─ Signaling pathway analysis
│   └─ CellChat (R) - Pathway-level insights, communication networks
│
├─ Fast computation needed
│   └─ FastCCC (C++) - Human only, very fast
│
└─ Classic/benchmark comparison
    └─ CellPhoneDB - Original method, good for comparisons
```

### Method Comparison

| Aspect | LIANA+ | CellChat | CellPhoneDB | FastCCC |
|--------|--------|----------|-------------|---------|
| Language | Python | R | Python | C++ |
| Speed | Fast | Moderate | Moderate | Very Fast |
| Output | LR pairs | Pathways + Networks | LR pairs | LR pairs |
| Species | Human/Mouse | Human/Mouse | Human | Human |
| Strength | Multi-method consensus | Pathway interpretation | Literature standard | Scale |

## Workflow

### Step 1: Verify Cell Type Annotations

```
Required in adata.obs:
- Cell type column (e.g., 'cell_type', 'annotation')
- Clean labels (no "unknown", "unassigned")
- Biologically meaningful categories
```

### Step 2: Configure Species

```
Species-specific databases:
- Human: liana_resource="consensus" (default)
- Mouse: liana_resource="mouseconsensus"
```

### Step 3: Run Analysis

Use `analyze_cell_communication` tool with:
- `cell_type_key`: Column name for cell type annotations
- `species`: "human" or "mouse"
- `method`: Selected from decision framework

### Step 4: Interpret Results

**Key outputs:**
- **Ligand-receptor pairs**: Specific molecular interactions
- **Interaction scores**: Strength/significance of communication
- **Cell type pairs**: Which cells are talking to which
- **Pathways** (CellChat): Grouped signaling pathways

### Step 5: Visualize

Use `visualize_data` with `plot_type="communication"`:
- `subtype="dotplot"`: Overview of all interactions
- `subtype="circle_plot"`: Network of cell type communication
- `subtype="tileplot"`: Heatmap of interaction strength

## Spatial Considerations

### Adding Spatial Context

Standard cell communication ignores distance. To incorporate spatial proximity:

1. **Pre-filter by spatial neighbors**: Only consider cells that are spatially adjacent
2. **Use spatial statistics**: Correlate communication scores with co-localization

### Interpreting Spatial Patterns

Questions to consider:
- Are communicating cell types spatially co-localized?
- Do interaction hotspots correspond to specific tissue regions?
- How does communication differ across spatial domains?

## Common Analysis Scenarios

### Scenario 1: Tumor Microenvironment
```
Focus on:
- Tumor-immune interactions
- CAF-tumor communication
- Immune checkpoint pathways (PD-1/PD-L1)
```

### Scenario 2: Development/Differentiation
```
Focus on:
- Signaling gradients (Wnt, BMP, Notch)
- Niche-stem cell interactions
- Temporal progression of communication
```

### Scenario 3: Tissue Architecture
```
Focus on:
- ECM-cell interactions
- Structural organization signals
- Regional communication patterns
```

## Output Interpretation Guide

### Significant Interactions

Look for:
- High specificity scores (interaction specific to certain cell pairs)
- Biological relevance (known pathways in tissue type)
- Spatial coherence (co-localization of interacting cells)

### Common Patterns

| Pattern | Interpretation |
|---------|---------------|
| Hub cell type | Central coordinator of tissue communication |
| Reciprocal pairs | Bidirectional signaling between cell types |
| Pathway enrichment | Dominant signaling mechanism |
| Spatial clustering | Localized communication zones |

## Troubleshooting

| Issue | Solution |
|-------|----------|
| No significant interactions | Lower threshold, check cell type balance |
| Too many interactions | Increase stringency, focus on specific pathways |
| CellChat R errors | Verify rpy2 installation, check R packages |
| Species mismatch | Ensure correct database for species |
