# Configuration file for the Sphinx documentation builder.
# Simplified configuration for ChatSpatial documentation

import os
import sys
from datetime import datetime

# -- Path setup --------------------------------------------------------------
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
project = "ChatSpatial"
copyright = f"{datetime.now().year}, ChatSpatial Project"
author = "ChatSpatial Contributors"
release = "1.0.0"
version = "1.0"

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",  # Markdown support
    "sphinx_copybutton",  # Copy button for code blocks
    "sphinx_design",  # Design elements (cards, tabs, grids) for index.rst
]

# MyST-Parser configuration - enable useful Markdown extensions
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "linkify",
    "tasklist",
]

# Support both .rst and .md files
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Patterns to exclude from source files
exclude_patterns = [
    "_build",
    "_archive",  # Exclude archived documentation from indexing
    "Thumbs.db",
    ".DS_Store",
]

# The master toctree document
master_doc = "index"

# -- Options for HTML output -------------------------------------------------

# Use the default Read the Docs theme (alabaster is also clean and simple)
html_theme = "alabaster"

# Alabaster theme options - minimal and clean
html_theme_options = {
    "logo": "logo.png",
    "github_user": "cafferychen777",
    "github_repo": "ChatSpatial",
    "github_button": True,
    "github_banner": False,
    "description": "AI-Powered Spatial Transcriptomics Analysis",
    "fixed_sidebar": True,
}

# Path for static files (images, CSS, etc.)
html_static_path = ["_static"]

# Don't use custom CSS - let the theme handle everything
html_css_files = []

# Show "Created using Sphinx" in footer
html_show_sphinx = True

# Output file base name for HTML help builder
htmlhelp_basename = "ChatSpatialdoc"

# -- Options for copybutton extension ----------------------------------------
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
copybutton_remove_prompts = True
