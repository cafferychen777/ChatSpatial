# Sphinx + PyData Theme Migration Guide

## Overview

ChatSpatial documentation has been migrated to **Sphinx with PyData Sphinx Theme** - the standard documentation system used by the scientific Python community.

## Why Sphinx + PyData Theme?

### Scientific Python Standard

PyData Sphinx Theme is the de facto standard for spatial transcriptomics and bioinformatics documentation:

- **Scanpy** - https://scanpy.readthedocs.io/
- **Squidpy** - https://squidpy.readthedocs.io/
- **AnnData** - https://anndata.readthedocs.io/
- **EnrichMap** - https://enrichmap.readthedocs.io/
- **NumPy**, **Pandas**, **SciPy**, and many more

### Key Advantages

1. **Community Alignment**: Users familiar with other spatial analysis tools will recognize the interface
2. **Rich Features**: Advanced cross-referencing, API documentation, bibliography support
3. **Professional Look**: Clean, academic-oriented design
4. **Better Integration**: Native intersphinx linking to other scientific packages
5. **Extensibility**: Huge ecosystem of Sphinx extensions

## What Changed

### Build System

- **Before**: MkDocs Material (modern, general-purpose)
- **After**: Sphinx + PyData Theme (scientific Python standard)

### File Format

- **Primary**: reStructuredText (.rst) for main structure pages
- **Secondary**: Markdown (.md) via MyST-Parser for content pages
- Both formats are fully supported side-by-side

### Configuration

- **Removed**: `mkdocs.yml`, `site/` directory
- **Added**: `docs/conf.py`, `docs/_build/` directory

### Theme Features

| Feature | PyData Theme | MkDocs Material |
|---------|--------------|-----------------|
| Sidebar Navigation | ✅ Yes | ❌ No (top tabs) |
| API Documentation | ✅ Native | ⚠️ Limited |
| Scientific Citations | ✅ Built-in | ❌ No |
| Code Documentation | ✅ Autodoc | ⚠️ Manual |
| Cross-references | ✅ Advanced | ⚠️ Basic |
| Community Standard | ✅ Scientific Python | ✅ General |

## New Features

### 1. Sphinx Extensions Enabled

- **sphinx.ext.autodoc**: Automatic API documentation from docstrings
- **sphinx.ext.napoleon**: Google/NumPy docstring support
- **sphinx.ext.intersphinx**: Link to other project docs (Scanpy, AnnData, etc.)
- **sphinx_copybutton**: Copy code blocks easily
- **sphinx_design**: Grid layouts, cards, tabs
- **sphinx_tabs**: Tabbed content
- **sphinx_togglebutton**: Collapsible sections
- **sphinxext.opengraph**: Social media preview cards
- **sphinxcontrib.bibtex**: Bibliography and citations

### 2. Grid Layouts

Beautiful card-based layouts on index pages:

```rst
.. grid:: 2
    :gutter: 3

    .. grid-item-card:: Getting Started
        :link: getting-started/installation
        :link-type: doc

        Install and configure ChatSpatial

    .. grid-item-card:: Tutorials
        :link: tutorials/index
        :link-type: doc

        Step-by-step guides
```

### 3. Intersphinx Linking

Reference other packages directly:

```rst
See :class:`anndata.AnnData` for data structure details.
Use :func:`scanpy.pp.normalize_total` for normalization.
```

### 4. API Documentation

Automatic documentation from Python docstrings:

```rst
.. automodule:: chatspatial.tools.annotation
   :members:
   :undoc-members:
```

## Building Documentation

### Local Development

```bash
# Build documentation
cd docs
sphinx-build -b html . _build/html

# Or use make (if Makefile exists)
make html

# Live reload with sphinx-autobuild
sphinx-autobuild . _build/html --port 8000
```

### Serve Locally

```bash
# Simple HTTP server
python -m http.server 8080 --directory _build/html

# Visit: http://127.0.0.1:8080
```

### Read the Docs

Documentation is automatically built on Read the Docs when:
- Commits pushed to main branch
- Pull requests created (preview builds)

Configuration: `.readthedocs.yaml`

## File Organization

```
docs/
├── conf.py                 # Sphinx configuration
├── index.rst              # Home page (reStructuredText)
├── _static/               # Static files (CSS, images, JS)
│   └── custom.css
├── _templates/            # Custom templates (optional)
├── references.bib         # Bibliography file
├── getting-started/
│   ├── index.rst         # Section index (reStructuredText)
│   ├── installation.md   # Content (Markdown)
│   ├── quick-start.md
│   └── configuration.md
├── tutorials/
│   ├── index.rst
│   └── ...
└── _build/               # Build output (git-ignored)
    └── html/
```

### When to Use .rst vs .md

**Use .rst for**:
- Index pages with complex structure
- Pages with advanced Sphinx directives
- API documentation
- Pages with bibliographies

**Use .md for**:
- Tutorial content
- Simple narrative documentation
- Existing content that doesn't need restructuring

## Writing Documentation

### reStructuredText Basics

```rst
Title
=====

Subtitle
--------

**Bold text** and *italic text*

- Bullet point
- Another point

#. Numbered item
#. Auto-numbered

Link to another page: :doc:`getting-started/installation`
External link: `ChatSpatial GitHub <https://github.com/cafferychen777/ChatSpatial>`_

.. code-block:: python

    def hello():
        print("Hello, World!")

.. note::
    This is a note callout

.. warning::
    This is a warning
```

### Markdown with MyST

ChatSpatial uses MyST-Parser for enhanced Markdown:

```markdown
# Title

## Subtitle

**Bold** and *italic*

- Bullet points work
- As expected

```python
# Code blocks
def hello():
    print("Hello!")
\```

:::{note}
MyST admonitions use ::: syntax
:::

[Link to docs](getting-started/installation.md)
```

## Customization

### Theme Options

Edit `docs/conf.py`:

```python
html_theme_options = {
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["navbar-icon-links", "theme-switcher"],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/...",
            "icon": "fab fa-github-square",
        },
    ],
    "logo": {
        "text": "ChatSpatial",
        "image_light": "_static/logo-light.svg",
        "image_dark": "_static/logo-dark.svg",
    },
}
```

### Custom CSS

Edit `docs/_static/custom.css`:

```css
:root {
    --pst-color-primary: #2196F3;
    --pst-color-secondary: #4CAF50;
}

.sd-card {
    border-radius: 8px;
}
```

## Common Tasks

### Add a New Page

1. Create markdown or .rst file in appropriate directory
2. Add to toctree in parent index.rst:

```rst
.. toctree::
   :maxdepth: 2

   existing-page
   new-page
```

### Add to Navigation

Edit section's `index.rst`:

```rst
.. toctree::
   :maxdepth: 2
   :caption: Section Name

   page1
   page2
   new-page
```

### Add Images

Place images in `docs/_static/` or `docs/images/`:

```rst
.. image:: _static/diagram.png
   :alt: Diagram description
   :width: 600px
```

Or in Markdown:

```markdown
![Diagram](/_static/diagram.png)
```

### Cross-Reference

Within ChatSpatial docs:

```rst
See :doc:`/getting-started/installation` for setup.
See :ref:`section-label` for details.
```

To other packages (via intersphinx):

```rst
:class:`scanpy.AnnData`
:func:`numpy.array`
```

## Troubleshooting

### Build Errors

**"Document may not begin with a transition"**
- Remove `---` at start of markdown files
- Run `python fix_markdown_transitions.py` (now removed)

**"Unknown directive type"**
- Check extension is enabled in `conf.py`
- Verify correct syntax for directive

**"WARNING: toctree contains reference to nonexisting document"**
- Check file exists
- Verify path is correct (no `.rst` or `.md` extension in toctree)

### Theme Issues

**Custom CSS not loading**
- Verify `_static/custom.css` exists
- Check `html_css_files` in `conf.py`
- Clear browser cache

**Navigation not showing**
- Check toctree structure in index.rst files
- Verify `:maxdepth:` settings
- Rebuild with `sphinx-build -E` (force rebuild)

## Resources

### Official Documentation

- [Sphinx](https://www.sphinx-doc.org/)
- [PyData Sphinx Theme](https://pydata-sphinx-theme.readthedocs.io/)
- [MyST-Parser](https://myst-parser.readthedocs.io/)
- [Read the Docs](https://docs.readthedocs.io/)

### Examples from Scientific Community

- [Scanpy Documentation](https://github.com/scverse/scanpy/tree/main/docs)
- [Squidpy Documentation](https://github.com/scverse/squidpy/tree/main/docs)
- [NumPy Documentation](https://github.com/numpy/numpy/tree/main/doc)

### Sphinx Extensions

- [sphinx-design](https://sphinx-design.readthedocs.io/) - Cards, grids, tabs
- [sphinx-copybutton](https://sphinx-copybutton.readthedocs.io/) - Copy buttons
- [sphinxcontrib-bibtex](https://sphinxcontrib-bibtex.readthedocs.io/) - Citations

## Migration Checklist

- [x] Install Sphinx and PyData theme
- [x] Create `conf.py` configuration
- [x] Convert main index pages to .rst
- [x] Fix markdown transition errors
- [x] Update Read the Docs configuration
- [x] Remove MkDocs files (mkdocs.yml, site/)
- [x] Remove Jekyll files (_config.yml, vendor/, etc.)
- [x] Update .gitignore for Sphinx
- [x] Test local build
- [x] Configure additional extensions
- [ ] Deploy to Read the Docs
- [ ] Update main README with new docs URL

## Contributing

When contributing documentation:

1. Keep narrative content in Markdown (.md)
2. Use reStructuredText (.rst) for structural pages
3. Test build locally before submitting PR
4. Follow existing formatting conventions
5. Add new pages to appropriate toctree

## Support

- [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
- [GitHub Discussions](https://github.com/cafferychen777/ChatSpatial/discussions)
- [Read the Docs Documentation](https://docs.readthedocs.io/)
