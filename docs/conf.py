# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
from datetime import date

# Add the src directory to the path
sys.path.insert(0, os.path.abspath('../src'))

import pyrvt

project = "pyRVT"
copyright = f"2015-{date.today().year}, Albert R. Kottke"
author = "Albert R. Kottke"
release = pyrvt.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
    "sphinxext.opengraph",
    "sphinx_copybutton",
    "myst_parser",
    "nbsphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Bibliography configuration
bibtex_bibfiles = ["refs.bib"]
bibtex_reference_style = "author_year"

# Autosummary configuration
autosummary_generate = True

# Napoleon configuration for NumPy-style docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# Intersphinx configuration
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

# MyST configuration
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "deflist",
    "fieldlist",
    "html_admonition",
    "html_image",
    "colon_fence",
    "smartquotes",
    "replacements",
    "linkify",
    "strikethrough",
    "substitution",
    "tasklist",
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
pygments_style = "friendly"

html_theme = "furo"
html_static_path = ["_static"]
html_title = f"{project} v{release}"

# Furo theme options
html_theme_options = {
    "sidebar_hide_name": True,
    "navigation_with_keys": True,
    "top_of_page_button": "edit",
    "source_repository": "https://github.com/arkottke/pyrvt/",
    "source_branch": "main",
    "source_directory": "docs/",
}

# Show typehints as content of the function or method
autodoc_typehints = "description"
autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "inherited-members": True,
    "private-members": False,
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
    "no-index": False,
}

# Notebook execution configuration
nbsphinx_execute = "never"  # Don't execute notebooks during build
nbsphinx_allow_errors = True
