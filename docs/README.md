# ChatSpatial Documentation

This directory contains the Jekyll-based documentation for ChatSpatial using the Just the Docs theme.

## Local Development

### Prerequisites

- Ruby 3.1+
- Bundler gem
- Git

### Setup

1. Install dependencies:
   ```bash
   cd docs
   bundle install
   ```

2. Serve locally:
   ```bash
   bundle exec jekyll serve
   ```

3. Open http://localhost:4000 in your browser

### Live Reload

For automatic rebuilds during development:
```bash
bundle exec jekyll serve --livereload
```

## Project Structure

```
docs/
├── _config.yml          # Jekyll configuration
├── Gemfile              # Ruby dependencies
├── index.md             # Homepage
├── getting-started/     # Getting started guides
│   ├── index.md
│   ├── installation.md
│   ├── quick-start.md
│   └── configuration.md
├── tutorials/           # Tutorial content
│   ├── index.md
│   ├── core/           # Core analysis tutorials
│   ├── analysis/       # Method-specific tutorials
│   ├── advanced/       # Advanced methods
│   └── learning-paths/ # Structured learning paths
├── reference/           # Reference documentation
│   ├── index.md
│   ├── api/            # API documentation
│   ├── quick-reference/ # Quick reference guides
│   └── troubleshooting/ # Problem-solving guides
└── examples/            # Example workflows
```

## Theme Customization

This site uses the [Just the Docs](https://just-the-docs.github.io/just-the-docs/index.md) theme.

### Custom Styling

Add custom CSS in `_sass/custom/custom.scss`:

```scss
// Custom variables
$body-font-family: -apple-system, BlinkMacSystemFont, "helvetica neue", helvetica, roboto, noto, "segoe ui", arial, sans-serif;
$mono-font-family: "SFMono-Regular", Menlo, Consolas, Monospace;

// Custom styles
.site-header {
  background-color: #f8f9fa;
}
```

### Navigation

Navigation is controlled by:
- `nav_order` in frontmatter
- `has_children: true` for parent pages
- `parent: Page Name` for child pages

### Callouts

Use Just the Docs callouts for highlights:

```markdown

Important information

Additional notes

Warning messages
```

## Content Guidelines

### File Naming
- Use lowercase with hyphens: `quick-start.md`
- Match navigation structure
- Include index.md for directories

### Frontmatter
Required for all pages:
```yaml
---
layout: default
title: Page Title
nav_order: 1
parent: Parent Page (if applicable)
---
```

### Images
- Store in `assets/images/`
- Use descriptive filenames
- Include alt text for accessibility

### Code Examples
- Use language-specific syntax highlighting
- Include working examples that users can copy
- Test all code examples before publishing

## Deployment

### GitHub Pages

The site is automatically deployed via GitHub Actions when changes are pushed to the main branch.

Configuration in `.github/workflows/deploy-docs.yml`:
- Builds on push to main
- Uses Ruby setup action
- Deploys to GitHub Pages

### Manual Deployment

To deploy manually:
```bash
# Build the site
bundle exec jekyll build

# Deploy _site directory to your hosting provider
```

## Contributing

### Adding New Content

1. Create new markdown files following the structure
2. Add appropriate frontmatter
3. Update navigation if needed
4. Test locally before committing

### Updating Existing Content

1. Edit the relevant markdown files
2. Check links and references
3. Test locally
4. Commit changes

### Style Guide

- Use clear, concise language
- Include code examples for all procedures
- Add troubleshooting sections where appropriate
- Link between related pages
- Use consistent formatting and terminology

## Troubleshooting

### Bundle Install Issues

```bash
# Clear bundle cache
bundle clean --force

# Reinstall dependencies
bundle install
```

### Jekyll Build Errors

```bash
# Clean Jekyll cache
bundle exec jekyll clean

# Rebuild with verbose output
bundle exec jekyll build --verbose
```

### Dependency Conflicts

```bash
# Update all gems
bundle update

# Check for security issues
bundle audit
```

## Links

- [Jekyll Documentation](https://jekyllrb.com/docs/index.md)
- [Just the Docs Theme](https://just-the-docs.github.io/just-the-docs/index.md)
- [GitHub Pages](https://docs.github.com/en/pages)