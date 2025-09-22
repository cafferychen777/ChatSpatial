---
layout: default
title: Style Guide
nav_order: 99
nav_exclude: true
description: "ChatSpatial design system and component showcase"
---

# ChatSpatial Style Guide
{: .no_toc }

A showcase of the ChatSpatial documentation design system, components, and brand elements.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Brand Colors

Our color palette is inspired by spatial transcriptomics visualization and modern bioinformatics design:

<div class="features-grid">
  <div class="feature-card">
    <div style="background: linear-gradient(135deg, #0ea5e9 0%, #0369a1 100%); height: 60px; border-radius: 8px; margin-bottom: 1rem;"></div>
    <h4>Primary Blue</h4>
    <p>Main brand color used for buttons, links, and primary elements</p>
    <code>#0ea5e9 → #0369a1</code>
  </div>
  
  <div class="feature-card">
    <div style="background: linear-gradient(135deg, #d946ef 0%, #a21caf 100%); height: 60px; border-radius: 8px; margin-bottom: 1rem;"></div>
    <h4>Secondary Purple</h4>
    <p>Secondary brand color for accents and highlights</p>
    <code>#d946ef → #a21caf</code>
  </div>
  
  <div class="feature-card">
    <div style="background: linear-gradient(135deg, #22c55e 0%, #15803d 100%); height: 60px; border-radius: 8px; margin-bottom: 1rem;"></div>
    <h4>Accent Green</h4>
    <p>Success states and positive actions</p>
    <code>#22c55e → #15803d</code>
  </div>
</div>

---

## Typography

Our typography system uses Inter for body text and JetBrains Mono for code:

### Headings

# Heading 1 - Main page titles
## Heading 2 - Section headers  
### Heading 3 - Subsection headers
#### Heading 4 - Component headers
##### Heading 5 - Small headers
###### Heading 6 - Micro headers

### Body Text

This is a paragraph of body text. It uses the Inter font family for excellent readability across all devices and screen sizes. The line height and spacing are optimized for documentation reading.

**Bold text** and *italic text* are available for emphasis. You can also use `inline code` for technical terms.

### Code Blocks

```python
# ChatSpatial Python example
import chatspatial as cs

# Load spatial data
adata = cs.load_data("visium_data.h5ad")

# Analyze spatial patterns
results = cs.analyze_spatial_patterns(adata, method="moran")
```

---

## Buttons

Our button system provides consistent styling across different contexts:

<div style="display: flex; gap: 1rem; flex-wrap: wrap; margin: 2rem 0;">
  <button class="btn btn-primary">Primary Button</button>
  <button class="btn btn-outline">Outline Button</button>
  <button class="btn btn-secondary">Secondary Button</button>
  <button class="btn btn-success">Success Button</button>
  <button class="btn btn-ghost">Ghost Button</button>
</div>

<div style="display: flex; gap: 1rem; flex-wrap: wrap; margin: 2rem 0;">
  <button class="btn btn-primary btn-sm">Small Primary</button>
  <button class="btn btn-outline btn-sm">Small Outline</button>
  <button class="btn btn-primary btn-lg">Large Primary</button>
</div>

---

## Badges

Badges are used to display status, categories, or labels:

<div style="display: flex; gap: 0.5rem; flex-wrap: wrap; margin: 2rem 0;">
  <span class="badge badge-primary">Primary</span>
  <span class="badge badge-secondary">Secondary</span>
  <span class="badge badge-success">Success</span>
  <span class="badge badge-warning">Warning</span>
  <span class="badge badge-error">Error</span>
  <span class="badge badge-outline">Outline</span>
</div>

---

## Alerts

Our alert system provides contextual feedback to users:

<div class="alert alert-info">
  <div class="alert-content">
    <h4>Information</h4>
    <p>This is an informational alert. Use it to provide helpful context or additional information.</p>
  </div>
</div>

<div class="alert alert-success">
  <div class="alert-content">
    <h4>Success</h4>
    <p>This is a success alert. Use it to confirm that an action was completed successfully.</p>
  </div>
</div>

<div class="alert alert-warning">
  <div class="alert-content">
    <h4>Warning</h4>
    <p>This is a warning alert. Use it to draw attention to important information or potential issues.</p>
  </div>
</div>

<div class="alert alert-error">
  <div class="alert-content">
    <h4>Error</h4>
    <p>This is an error alert. Use it to indicate that something went wrong or requires immediate attention.</p>
  </div>
</div>

---

## Callouts

Just the Docs callouts with ChatSpatial styling:

{: .highlight }
This is a highlight callout. Use it to draw attention to important information.

{: .note }
This is a note callout. Use it for additional context or explanations.

{: .important }
This is an important callout. Use it for critical information that users must know.

{: .warning }
This is a warning callout. Use it to alert users about potential issues or limitations.

{: .new }
This is a new feature callout. Use it to highlight recently added functionality.

---

## Cards

Cards are used to group related content and create visual hierarchy:

<div class="features-grid">
  <div class="card">
    <div class="card-header">
      <h3>Card Title</h3>
    </div>
    <div class="card-body">
      <p>This is the card content. Cards are perfect for organizing related information and creating visual separation between different sections.</p>
    </div>
    <div class="card-footer">
      <button class="btn btn-primary btn-sm">Action</button>
      <button class="btn btn-outline btn-sm">Cancel</button>
    </div>
  </div>
  
  <div class="card">
    <div class="card-body">
      <h4>Simple Card</h4>
      <p>This card doesn't have a header or footer, just the main content area.</p>
    </div>
  </div>
</div>

---

## Progress Bars

Progress bars show completion status or loading states:

<div style="margin: 2rem 0;">
  <h4>Installation Progress</h4>
  <div class="progress">
    <div class="progress-bar" data-progress="75" style="width: 75%;"></div>
  </div>
  <p style="font-size: 0.875rem; color: var(--cs-gray-600); margin-top: 0.5rem;">75% complete</p>
</div>

<div style="margin: 2rem 0;">
  <h4>Analysis Progress</h4>
  <div class="progress">
    <div class="progress-bar" data-progress="45" style="width: 45%;"></div>
  </div>
  <p style="font-size: 0.875rem; color: var(--cs-gray-600); margin-top: 0.5rem;">45% complete</p>
</div>

---

## Tables

Tables with ChatSpatial styling:

| Method | Accuracy | Speed | Memory Usage |
|:-------|:---------|:------|:-------------|
| SpaGCN | 92.5% | Fast | Low |
| STAGATE | 89.3% | Medium | Medium |
| Leiden | 87.2% | Very Fast | Very Low |
| Louvain | 86.9% | Very Fast | Very Low |

---

## Code Examples

Enhanced code blocks with language labels and copy buttons:

```bash
# Install ChatSpatial
pip install chatspatial

# Run the MCP server
python -m chatspatial
```

```python
import chatspatial as cs
import scanpy as sc

# Load your spatial data
adata = cs.load_data("path/to/your/data.h5ad")

# Preprocess the data
cs.preprocess_data(adata, method="log")

# Identify spatial domains
domains = cs.identify_spatial_domains(adata, method="spagcn")

# Visualize results
cs.plot_spatial_domains(adata, color="spatial_domains")
```

---

## Responsive Design

Our design system is built mobile-first and works seamlessly across all device sizes:

- **Mobile** (< 500px): Single column, touch-friendly interactions
- **Tablet** (500px - 750px): Responsive grid, optimized for touch
- **Desktop** (750px - 1064px): Multi-column layout, hover states
- **Large** (> 1064px): Full layout with maximum content width

---

## Accessibility

ChatSpatial documentation follows WCAG 2.1 AA guidelines:

- ✅ Color contrast ratios meet AA standards
- ✅ Focus indicators for keyboard navigation
- ✅ Screen reader compatible markup
- ✅ Touch targets are at least 44px
- ✅ Reduced motion support
- ✅ Alternative text for images

---

## Logo Usage

The ChatSpatial logo combines spatial data visualization with AI communication:

<div style="display: flex; gap: 2rem; align-items: center; margin: 2rem 0; padding: 1.5rem; background: var(--cs-gray-50); border-radius: 12px;">
  <img src="/ChatSpatial/assets/images/logo.svg" alt="ChatSpatial Logo" style="height: 40px;">
  <img src="/ChatSpatial/assets/images/logo-icon.svg" alt="ChatSpatial Icon" style="height: 40px;">
</div>

### Logo Guidelines

- Use on light backgrounds with sufficient contrast
- Maintain minimum size of 32px height for the icon
- Preserve the aspect ratio and don't distort
- Allow clear space around the logo equal to the height

---

<div class="alert alert-info">
  <div class="alert-content">
    <h4>Using These Components</h4>
    <p>This style guide showcases the design system used throughout the ChatSpatial documentation. Most components are automatically styled by our CSS, while others can be added using the appropriate CSS classes shown in the examples above.</p>
  </div>
</div>