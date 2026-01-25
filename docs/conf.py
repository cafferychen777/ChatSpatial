# Configuration file for the Sphinx documentation builder.
# ChatSpatial Documentation - Modern Furo Theme

import os
import sys
from datetime import datetime

# -- Path setup --------------------------------------------------------------
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
project = "ChatSpatial"
copyright = f"{datetime.now().year}, ChatSpatial Project"
author = "ChatSpatial Contributors"
release = "1.1.1"
version = "1.1"

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",  # Markdown support
    "sphinx_copybutton",  # Copy button for code blocks
    "sphinx_design",  # Design elements (cards, tabs, grids)
]

# MyST-Parser configuration
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "tasklist",
]

# Support both .rst and .md files
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Patterns to exclude
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The master toctree document
master_doc = "index"

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"

# Furo theme options
html_theme_options = {
    # Branding
    "light_logo": "logo.png",
    "dark_logo": "logo.png",

    # Sidebar
    "sidebar_hide_name": False,

    # Navigation
    "navigation_with_keys": True,

    # Source links
    "source_repository": "https://github.com/cafferychen777/ChatSpatial",
    "source_branch": "main",
    "source_directory": "docs/",

    # Light mode colors - Clean scientific palette
    "light_css_variables": {
        # Primary colors - Deep blue for trust and professionalism
        "color-brand-primary": "#1565C0",
        "color-brand-content": "#1976D2",

        # Background and surfaces
        "color-background-primary": "#ffffff",
        "color-background-secondary": "#f8fafc",
        "color-background-hover": "#f1f5f9",
        "color-background-border": "#e2e8f0",

        # Sidebar
        "color-sidebar-background": "#f8fafc",
        "color-sidebar-background-border": "#e2e8f0",
        "color-sidebar-brand-text": "#1e293b",
        "color-sidebar-link-text": "#475569",
        "color-sidebar-link-text--top-level": "#1e293b",

        # Admonitions
        "color-admonition-background": "#f0f9ff",

        # Code
        "color-code-background": "#f8fafc",
        "color-inline-code-background": "#f1f5f9",

        # Links
        "color-link": "#1565C0",
        "color-link--hover": "#1976D2",
    },

    # Dark mode colors - Comfortable night viewing
    "dark_css_variables": {
        # Primary colors
        "color-brand-primary": "#60a5fa",
        "color-brand-content": "#93c5fd",

        # Background and surfaces
        "color-background-primary": "#0f172a",
        "color-background-secondary": "#1e293b",
        "color-background-hover": "#334155",
        "color-background-border": "#334155",

        # Sidebar
        "color-sidebar-background": "#1e293b",
        "color-sidebar-background-border": "#334155",
        "color-sidebar-brand-text": "#f1f5f9",
        "color-sidebar-link-text": "#94a3b8",
        "color-sidebar-link-text--top-level": "#f1f5f9",

        # Admonitions
        "color-admonition-background": "#1e293b",

        # Code
        "color-code-background": "#1e293b",
        "color-inline-code-background": "#334155",

        # Links
        "color-link": "#60a5fa",
        "color-link--hover": "#93c5fd",
    },

    # Footer icons
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/cafferychen777/ChatSpatial",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
}

# Static files
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# Page title
html_title = "ChatSpatial"

# Favicon
# html_favicon = "_static/favicon.ico"

# Don't show Sphinx attribution
html_show_sphinx = False

# Last updated format
html_last_updated_fmt = "%b %d, %Y"

# -- Copybutton configuration ------------------------------------------------
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: "
copybutton_prompt_is_regexp = True
copybutton_remove_prompts = True

# -- Pygments syntax highlighting --------------------------------------------
pygments_style = "friendly"
pygments_dark_style = "monokai"
