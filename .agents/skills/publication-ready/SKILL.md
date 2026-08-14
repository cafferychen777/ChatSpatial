---
name: publication-ready
description: |
  Generate publication-quality figures following journal standards for Nature, Cell, and other high-impact journals.
  Use when user wants to create figures for papers, presentations, or reports with professional formatting.
  Triggers: "publication figure", "paper figure", "Nature style", "journal figure",
  "export for publication", "figure panel", "high quality plot".
context: fork
agent: general-purpose
allowed-tools: Read, Grep, Glob
---

# Publication-Ready Figures

## Overview

This skill creates **publication-quality visualizations** following journal standards.

Key principles:
1. **Clarity over decoration**: Every element serves a purpose
2. **Consistent styling**: Uniform fonts, colors, and sizing
3. **Reproducibility**: Settings documented for revision cycles

## Journal Standards Quick Reference

### Figure Dimensions

| Journal | Single Column | 1.5 Column | Double Column |
|---------|--------------|------------|---------------|
| Nature | 89 mm | 120 mm | 183 mm |
| Cell | 85 mm | 114 mm | 174 mm |
| Science | 55 mm | 117 mm | 183 mm |

**Height**: Typically 1:1 to 1:1.5 aspect ratio; max ~230mm

### Resolution Requirements

| Format | Resolution | Use Case |
|--------|------------|----------|
| Print | 300-600 DPI | Final submission |
| Vector (PDF/SVG) | N/A | Preferred for line art |
| Web/Review | 150 DPI | Initial submission |

### Font Standards

```
Recommended fonts:
- Sans-serif: Arial, Helvetica (most common)
- Serif: Times New Roman (less common for figures)

Size guidelines:
- Axis labels: 7-8 pt
- Tick labels: 6-7 pt
- Panel labels (a, b, c): 8-10 pt, bold
- Annotations: 6-8 pt
```

## Color Guidelines

### Colorblind-Friendly Palettes

**Recommended categorical palette:**
```
Blue:    #0072B2
Orange:  #E69F00
Green:   #009E73
Yellow:  #F0E442
Sky:     #56B4E9
Red:     #D55E00
Pink:    #CC79A7
Gray:    #999999
```

**Sequential palettes:**
- Expression: viridis, magma, inferno
- Diverging: RdBu (centered at 0), coolwarm

### Color Usage Rules

1. **Maximum 7-8 colors** for categorical data
2. **Consistent color mapping** across all panels
3. **Avoid red-green combinations** (colorblind unfriendly)
4. **White background** for most journals

## Multi-Panel Figure Layout

### Panel Organization

```
Standard layouts:

2-panel: [a][b]

3-panel: [a][b]    or    [a][b][c]
         [c]

4-panel: [a][b]
         [c][d]

6-panel: [a][b][c]
         [d][e][f]
```

### Panel Labels

- Position: Top-left corner of each panel
- Style: Lowercase bold (a, b, c) or uppercase (A, B, C)
- Font: Same as figure, 8-10 pt bold

## Visualization Type Guidelines

### Spatial Plots

```
Essential elements:
- Scale bar (required for spatial)
- Clear spot/cell boundaries
- Consistent color scale across samples

Avoid:
- Overcrowded labels
- Low contrast colors
- Missing coordinate axes
```

### Heatmaps

```
Essential elements:
- Dendrogram if clustered
- Row/column labels (or clear groupings)
- Color scale bar with range

Settings:
- Center diverging scales at 0 or mean
- Use appropriate normalization (z-score for expression)
```

### UMAP/Embedding Plots

```
Essential elements:
- Legend (positioned outside if crowded)
- Point size appropriate for density
- Consistent coloring with spatial plots

Avoid:
- Axis labels (UMAP axes are arbitrary)
- Grid lines
- Cluttered legends
```

### Violin/Box Plots

```
Essential elements:
- Y-axis label with units
- Statistical comparisons if relevant
- Consistent colors with other panels

Good practices:
- Order groups logically
- Add individual points for small n
- Include median/quartile indicators
```

## Export Workflow

### Step 1: Configure Visualization

Use `visualize_data` with:
```
params:
  output_format: "pdf"  # Vector format preferred
  dpi: 300              # For raster elements
  figsize: [3.5, 3.5]   # In inches (89mm ≈ 3.5in)
```

### Step 2: Export Options

| Format | When to Use | File Size |
|--------|-------------|-----------|
| PDF | Vector graphics, line art | Small |
| SVG | Web, editable graphics | Small |
| PNG | Raster, screenshots | Medium |
| TIFF | Print submission (some journals) | Large |
| EPS | Legacy systems | Small |

### Step 3: Post-Processing (if needed)

For complex multi-panel figures:
1. Export individual panels as PDF/SVG
2. Assemble in vector software (Illustrator, Inkscape)
3. Add panel labels and annotations
4. Export final composite

## Common Figure Types for Spatial Analysis

### Figure 1: Data Overview
```
Typical panels:
a) Spatial scatter - clusters/domains
b) UMAP embedding
c) Violin plot - key markers
d) Spatial expression - signature genes
```

### Figure 2: Cell Composition
```
Typical panels:
a) Deconvolution - spatial proportions
b) Reference UMAP with cell types
c) Individual cell type maps
d) Proportion bar plots by region
```

### Figure 3: Cell Communication
```
Typical panels:
a) Circle plot - communication network
b) Dotplot - top interactions
c) Spatial L-R pairs
d) Pathway analysis
```

### Figure 4: Dynamics
```
Typical panels:
a) Velocity stream plot
b) Pseudotime trajectory
c) Gene trends along trajectory
d) Fate probability map
```

## Quality Checklist

Before submission:
- [ ] All text is legible at print size
- [ ] Colors are consistent across panels
- [ ] Scale bars present on spatial plots
- [ ] Legends are complete and clear
- [ ] File format matches journal requirements
- [ ] Resolution meets minimum standards
- [ ] Panel labels are properly positioned
- [ ] No overlapping elements

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Text too small | Increase figure size or font size |
| Colors muddy in print | Test grayscale conversion |
| File too large | Use vector format, reduce DPI for raster |
| Legend overlaps data | Move legend outside plot area |
| Inconsistent styling | Create style template, apply to all |
