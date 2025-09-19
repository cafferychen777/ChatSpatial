# Navigation Structure Test

This document summarizes the Jekyll navigation structure we've configured for ChatSpatial documentation:

## Main Navigation (Level 1)

1. **Home** (nav_order: 1) - `/docs/index.md`
2. **What is MCP?** (nav_order: 2) - `/docs/resources/what_is_mcp.md`
3. **Getting Started** (nav_order: 3) - `/docs/getting-started/README.md` [has_children]
4. **Deployment** (nav_order: 4) - `/docs/deployment/README.md` [has_children]
5. **Tutorials** (nav_order: 5) - `/docs/tutorials/README.md` [has_children]
6. **Reference** (nav_order: 6) - `/docs/reference/README.md` [has_children]
7. **Examples** (nav_order: 7) - `/docs/examples/README.md` [has_children]

## Getting Started (Level 2)

- **Installation** (nav_order: 1) - `/docs/getting-started/installation.md`
- **Quick Start** (nav_order: 2) - `/docs/getting-started/quick-start.md`

## Deployment (Level 2)

- **Server Deployment** (nav_order: 1) - `/docs/deployment/server-deployment.md`

## Tutorials (Level 2)

1. **Core Tutorials** (nav_order: 1) - `/docs/tutorials/core/README.md` [has_children]
2. **Analysis Tutorials** (nav_order: 2) - `/docs/tutorials/analysis/README.md` [has_children]
3. **Advanced Tutorials** (nav_order: 3) - `/docs/tutorials/advanced/README.md` [has_children]
4. **Learning Paths** (nav_order: 4) - `/docs/tutorials/learning-paths/README.md` [has_children]

### Core Tutorials (Level 3)

- **Basic Spatial Analysis** (nav_order: 1) - `/docs/tutorials/core/basic_spatial_analysis.md`
- **Spatial Statistics** (nav_order: 2) - `/docs/tutorials/core/spatial_statistics.md`
- **Visualization Tutorial** (nav_order: 3) - `/docs/tutorials/core/visualization_tutorial.md`

### Analysis Tutorials (Level 3)

- **Cell Type Annotation** (nav_order: 1) - `/docs/tutorials/analysis/cell_type_annotation.md`
- **Cell Communication Analysis** (nav_order: 2) - `/docs/tutorials/analysis/cell_communication_analysis.md`
- **Spatial Enrichment** (nav_order: 3) - `/docs/tutorials/analysis/spatial_enrichment.md`

### Advanced Tutorials (Level 3)

- **GASTON Analysis** (nav_order: 1) - `/docs/tutorials/advanced/gaston_analysis.md`
- **Batch Integration** (nav_order: 2) - `/docs/tutorials/advanced/batch_integration.md`
- **Spatial Registration** (nav_order: 3) - `/docs/tutorials/advanced/spatial_registration.md`
- **Trajectory Analysis** (nav_order: 4) - `/docs/tutorials/advanced/trajectory_analysis.md`

### Learning Paths (Level 3)

- **Beginner Path** (nav_order: 1) - `/docs/tutorials/learning-paths/beginner.md`
- **Intermediate Path** (nav_order: 2) - `/docs/tutorials/learning-paths/intermediate.md`
- **Advanced Path** (nav_order: 3) - `/docs/tutorials/learning-paths/advanced.md`

## Reference (Level 2)

- **API Reference** (nav_order: 1) - `/docs/reference/api/README.md` [has_children]
- **Configuration** (nav_order: 2) - `/docs/reference/configuration.md`
- **Troubleshooting** (nav_order: 3) - `/docs/reference/troubleshooting/common_issues.md`

## Examples (Level 2)

- **Datasets** (nav_order: 1) - `/docs/examples/datasets/README.md`

## Features Configured

### Just the Docs Theme Features

- **Search**: Enabled with customized settings
- **Color Scheme**: Auto (light/dark mode)
- **Code Highlighting**: Rouge with line numbers
- **Mermaid Diagrams**: Enabled
- **SEO**: Jekyll SEO tag plugin
- **Callouts**: Multiple types (highlight, important, new, note, warning)
- **Back to Top**: Enabled
- **Edit on GitHub**: Enabled
- **Last Modified**: Enabled

### Navigation Features

- **Multi-level Navigation**: Up to 3 levels deep
- **Ordered Navigation**: Custom nav_order for all pages
- **Breadcrumbs**: Automatic via parent/grand_parent relationships
- **Collapsible Sections**: has_children for expandable navigation

### Additional Features

- **Social Links**: GitHub, Twitter in footer
- **Aux Links**: GitHub repo and PyPI package
- **Footer**: Copyright and license information
- **Plugins**: Feed, sitemap, SEO tag, mermaid

## Configuration Files

- `/docs/_config.yml` - Main Jekyll configuration
- `/docs/Gemfile` - Ruby dependencies
- `/.github/workflows/pages.yml` - GitHub Pages deployment

## Testing

To test locally:

```bash
cd docs
bundle install
bundle exec jekyll serve
```

Site should be available at `http://localhost:4000/ChatSpatial/`