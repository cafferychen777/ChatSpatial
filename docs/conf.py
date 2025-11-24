# Configuration file for the Sphinx documentation builder.
# This file follows the PyData Sphinx Theme style used by scientific Python projects
#
# For the full list of built-in configuration values, see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime

# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "ChatSpatial"
copyright = f"{datetime.now().year}, ChatSpatial Project"
author = "ChatSpatial Contributors"
release = "1.0.0"
version = "1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # Sphinx built-in extensions
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    # Third-party extensions
    "myst_parser",  # Markdown support
    "sphinx_copybutton",  # Copy button for code blocks
    "sphinx_design",  # Design elements (cards, tabs, grids)
    "sphinxcontrib.bibtex",  # Bibliography support
    "sphinx_tabs.tabs",  # Tabbed content
    "sphinx_togglebutton",  # Toggle/collapsible buttons
    "sphinxext.opengraph",  # Open Graph meta tags for social media
]

# MyST-Parser configuration
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

# Support both .rst and .md files
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Patterns to exclude from source files
templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "_site",
    ".jekyll-cache",
    "Thumbs.db",
    ".DS_Store",
    "vendor",
    ".bundle",
    "Gemfile*",
    "_config.yml",
    "Rakefile",
    "convert_frontmatter.py",
    "MIGRATION_TO_READTHEDOCS.md",
    # Exclude Jekyll-specific directories
    "_includes",
    "_tutorials",
    "javascripts",
    # Exclude design and technical docs from main site
    "design/*",
    "technical/*",
    "reports/*",
    "deployment/*",
]

# The master toctree document
master_doc = "index"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]

# PyData theme options
# See: https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/layout.html
html_theme_options = {
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["navbar-icon-links", "theme-switcher"],
    "navbar_persistent": ["search-button"],
    "header_links_before_dropdown": 5,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/cafferychen777/ChatSpatial",
            "icon": "fab fa-github-square",
            "type": "fontawesome",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/chatspatial/",
            "icon": "fab fa-python",
            "type": "fontawesome",
        },
    ],
    "logo": {
        "text": "ChatSpatial",
        "image_light": "_static/logo.png",
        "image_dark": "_static/logo.png",
    },
    "use_edit_page_button": True,
    "show_toc_level": 2,
    "navigation_depth": 4,
    "show_nav_level": 1,
    "navigation_with_keys": True,
    "collapse_navigation": False,
    "search_bar_text": "Search documentation...",
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version", "theme-version"],
    "secondary_sidebar_items": ["page-toc", "edit-this-page", "sourcelink"],
}

# GitHub repository for edit button
html_context = {
    "github_user": "cafferychen777",
    "github_repo": "ChatSpatial",
    "github_version": "main",
    "doc_path": "docs",
}

# -- Options for intersphinx extension ---------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#configuration

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
    "anndata": ("https://anndata.readthedocs.io/en/latest/", None),
    "squidpy": ("https://squidpy.readthedocs.io/en/stable/", None),
}

# -- Options for autodoc extension -------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#configuration

autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
}

# -- Options for autosummary extension ---------------------------------------
autosummary_generate = True

# -- Options for napoleon extension ------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True

# -- Options for copybutton extension ----------------------------------------
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True
copybutton_remove_prompts = True
copybutton_line_continuation_character = "\\"

# -- Options for bibliography ------------------------------------------------
bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "unsrt"

# -- Additional configuration ------------------------------------------------

# The name of the Pygments (syntax highlighting) style to use
pygments_style = "sphinx"
pygments_dark_style = "monokai"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css"
html_css_files = [
    "custom.css",
]

# Custom sidebar templates, must be a dictionary that maps document names
# to template names
html_sidebars = {
    "**": ["sidebar-nav-bs"],
}

# If true, "Created using Sphinx" is shown in the HTML footer
html_show_sphinx = True

# Output file base name for HTML help builder
htmlhelp_basename = "ChatSpatialdoc"

# -- Options for LaTeX output ------------------------------------------------
latex_elements = {
    "papersize": "letterpaper",
    "pointsize": "10pt",
    "preamble": "",
    "figure_align": "htbp",
}

# Grouping the document tree into LaTeX files
latex_documents = [
    (
        master_doc,
        "ChatSpatial.tex",
        "ChatSpatial Documentation",
        "ChatSpatial Contributors",
        "manual",
    ),
]

# -- Options for manual page output ------------------------------------------
man_pages = [
    (
        master_doc,
        "chatspatial",
        "ChatSpatial Documentation",
        [author],
        1,
    )
]

# -- Options for Texinfo output ----------------------------------------------
texinfo_documents = [
    (
        master_doc,
        "ChatSpatial",
        "ChatSpatial Documentation",
        author,
        "ChatSpatial",
        "AI-Powered Spatial Transcriptomics Analysis via Model Context Protocol",
        "Miscellaneous",
    ),
]

# -- Options for Epub output -------------------------------------------------
epub_title = project
epub_exclude_files = ["search.html"]

# -- Options for Open Graph extension ----------------------------------------
ogp_site_url = "https://chatspatial.readthedocs.io/"
ogp_site_name = "ChatSpatial Documentation"
ogp_image = "https://github.com/cafferychen777/ChatSpatial/raw/main/figures/chatspatial_architecture/chatspatial_architecture_nature_v3.png"
ogp_description_length = 200
ogp_type = "website"

# -- Options for sphinx-tabs extension --------------------------------------
sphinx_tabs_disable_tab_closing = True
