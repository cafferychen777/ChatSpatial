---
name: cell-dynamics
description: |
  Understand cellular differentiation, state transitions, and dynamic processes through trajectory and velocity analysis.
  Use when user wants to study development, differentiation, cell fate decisions, or temporal dynamics.
  Triggers: "trajectory", "pseudotime", "RNA velocity", "differentiation", "cell fate",
  "how do cells develop", "lineage", "dynamics", "scVelo", "CellRank".
context: fork
agent: general-purpose
allowed-tools: Read, Grep, Glob
---

# Cell Dynamics Analysis

## Overview

This skill answers: **How do cells change over time, and what are their fate decisions?**

Two complementary approaches:
1. **RNA Velocity**: Infers directionality from spliced/unspliced ratios
2. **Pseudotime/Trajectory**: Orders cells along developmental paths

## Decision Tree: Which Approach?

```
Q: Does your data have spliced/unspliced information?
│
├─ YES (velocyto/kallisto processed)
│   │
│   └─ RNA Velocity Analysis
│       ├─ scVelo (deterministic) - Fast, good default
│       ├─ scVelo (dynamical) - More accurate, slower
│       └─ VeloVI - Deep learning, handles noise
│       │
│       └─ Then: CellRank for fate probability
│
└─ NO (standard scRNA-seq/spatial)
    │
    └─ Pseudotime Analysis
        ├─ Palantir - Multi-lineage, fate probabilities
        ├─ DPT (Diffusion Pseudotime) - Classic, fast
        └─ PAGA - Trajectory + clustering
```

## RNA Velocity Workflow

### Prerequisites

Data must contain:
- `adata.layers['spliced']`: Spliced counts
- `adata.layers['unspliced']`: Unspliced counts

Check with: `'spliced' in adata.layers and 'unspliced' in adata.layers`

### Step 1: Preprocessing for Velocity

```
Velocity-specific preprocessing:
1. Filter genes by spliced/unspliced detection
2. Normalize spliced and unspliced separately
3. Compute moments (first/second order)
```

### Step 2: Compute Velocity

Use `analyze_velocity_data` tool with method selection:

| Mode | Speed | Accuracy | When to Use |
|------|-------|----------|-------------|
| deterministic | Fast | Good | Initial exploration |
| stochastic | Moderate | Better | Publication, noisy data |
| dynamical | Slow | Best | Complex dynamics, time recovery |

### Step 3: Visualize Velocity

Use `visualize_data` with `plot_type="velocity"`:
- `subtype="stream"`: Velocity streamlines on embedding
- `subtype="phase"`: Phase portraits for specific genes
- `subtype="proportions"`: Spliced/unspliced ratios
- `subtype="heatmap"`: Velocity across pseudotime

### Step 4: Fate Analysis with CellRank

After velocity computation:
```
1. Use analyze_trajectory_data with method="cellrank"
2. Identify initial/terminal states
3. Compute fate probabilities
4. Identify lineage driver genes
```

## Pseudotime Workflow (No Velocity Data)

### Method Selection

| Method | Strengths | When to Use |
|--------|-----------|-------------|
| Palantir | Multi-lineage, fate probs | Complex differentiation |
| DPT | Fast, robust | Simple trajectories |
| PAGA | Topology discovery | Unknown structure |

### Step 1: Define Root Cell

Critical for pseudotime:
```
Options for defining root:
1. Known marker gene expression (e.g., stem cell markers)
2. Specific cluster known to be progenitor
3. User-specified cell barcode
```

### Step 2: Compute Trajectory

Use `analyze_trajectory_data` tool with:
- `method`: "palantir", "dpt", or "cellrank"
- `root_cell` or `root_cluster`: Starting point

### Step 3: Visualize Trajectory

Use `visualize_data` with `plot_type="trajectory"`:
- `subtype="pseudotime"`: Cells colored by pseudotime
- `subtype="fate_map"`: Fate probability heatmap
- `subtype="gene_trends"`: Expression along trajectory
- `subtype="circular"`: Circular trajectory plot

## Biological Interpretation

### What Velocity Tells You

- **Direction of differentiation**: Where cells are heading
- **Dynamic genes**: Genes actively changing
- **State transitions**: Intermediate cell states
- **Terminal states**: End points of differentiation

### What Pseudotime Tells You

- **Ordering of cells**: Relative developmental time
- **Branching points**: Fate decision locations
- **Lineage relationships**: Which states lead to which
- **Driver genes**: Genes driving transitions

## Spatial Context

### Combining Velocity with Space

Questions to explore:
- Do velocity vectors align with spatial gradients?
- Are terminal states localized to specific regions?
- How does spatial organization relate to developmental time?

### Visualization Tips

- Overlay velocity arrows on spatial coordinates
- Color spatial plots by pseudotime
- Identify spatial zones of active differentiation

## Common Analysis Scenarios

### Tumor Heterogeneity
```
Focus on:
- Cancer stem cell identification
- Differentiation hierarchies within tumor
- EMT transitions
```

### Development
```
Focus on:
- Lineage specification
- Progenitor identification
- Temporal gene programs
```

### Regeneration/Wound Healing
```
Focus on:
- Activation states
- Re-differentiation paths
- Spatial gradients of repair
```

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| No velocity vectors | Low unspliced counts | Check data quality, try VeloVI |
| Random directions | Steady state cells | Focus on dynamic subpopulations |
| Wrong root assignment | Biological uncertainty | Try multiple roots, validate with markers |
| Disconnected trajectory | Batch effects | Integrate samples, increase neighbors |

## Output Summary

Successful dynamics analysis provides:
- **Pseudotime values**: Ordering of cells
- **Velocity vectors**: Direction of change
- **Fate probabilities**: Likelihood of terminal states
- **Driver genes**: Key regulators of transitions
- **Spatial mapping**: Where dynamics occur in tissue
