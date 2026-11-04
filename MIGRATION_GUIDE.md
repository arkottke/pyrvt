# Project Migration Guide: Hatch/MkDocs → uv/Sphinx

This guide provides a step-by-step process for migrating Python projects from Hatch+MkDocs to uv+Sphinx with standardized documentation structure.

## Overview

This migration covers:

- **Build system**: Hatch → uv + setuptools + setuptools-scm
- **Documentation**: MkDocs → Sphinx + Furo theme
- **Development workflow**: Hatch environments → scripts.sh + uv
- **Documentation structure**: Flat structure → organized sections

## Prerequisites

Before starting, ensure you have:

- Python 3.10+ installed
- Git repository with existing project
- Backup of your current configuration

## Step 1: Install uv

```bash
# On macOS and Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Or with pip
pip install uv
```

## Step 2: Update pyproject.toml

Replace your existing `pyproject.toml` with the new structure:

### 2.1 Build System Section

```toml
[build-system]
requires = ["setuptools>=64", "setuptools-scm>=8"]
build-backend = "setuptools.build_meta"
```

### 2.2 Project Section

```toml
[project]
name = "your-package-name"
dynamic = ["version"]
authors = [{ name = "Your Name", email = "your.email@example.com" }]
description = "Your package description."
readme = "README.md"
license = { text = "MIT" }  # or your license
requires-python = ">=3.10"
dependencies = [
    # Your core dependencies
    "numpy",
    "scipy",
    # Add others as needed
]

[project.optional-dependencies]
test = [
    "coverage[toml]",
    "pytest",
    "pytest-cov",
    # Add test-specific dependencies
]
docs = [
    "sphinx",
    "furo",
    "sphinx-copybutton",
    "sphinx-design",
    "sphinxcontrib-bibtex",
    "sphinxext-opengraph",
    "sphinx-autobuild",
    "matplotlib",
    "myst-parser",
    "nbsphinx",
    "pandoc",
    "ipython",
    "linkify-it-py",
    # Add your core dependencies here too for docs
]
style = [
    "black[colorama,jupyter]",
    "ruff",
]

[project.urls]
Homepage = "https://github.com/username/repository"
Documentation = "https://your-package.readthedocs.io"

[project.scripts]
your-cli-command = "your_package.module:main"  # if you have CLI

[tool.setuptools.packages.find]
where = ["src"]

[tool.setuptools_scm]
write_to = "src/your_package/_version.py"
```

### 2.3 Tool Configuration

```toml
[tool.ruff]
exclude = ["docs"]
line-length = 88
indent-width = 4
target-version = "py310"
extend-include = ["*.ipynb"]

[tool.ruff.format]
quote-style = "double"

[tool.ruff.lint]
select = ["E", "F", "W", "Q", "I"]
ignore = ["E203"]

[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_functions = ["test_*"]
addopts = ["--strict-markers", "--strict-config"]

[tool.coverage.run]
source = ["src/your_package"]
omit = ["*/tests/*", "*/test_*"]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "def __repr__",
    "raise AssertionError",
    "raise NotImplementedError",
]
```

## Step 3: Update Package Structure

### 3.1 Version Handling

Update `src/your_package/__init__.py`:

```python
"""Your package description."""

try:
    from ._version import __version__
except ImportError:
    # For development installs
    __version__ = "unknown"

__all__ = ["__version__"]
```

### 3.2 Create Development Scripts

Create `scripts.sh` in project root:

```bash
#!/usr/bin/env sh
# Development scripts for [PROJECT NAME] using uv
# Compatible with bash, zsh, and other POSIX shells

set -e

case "$1" in
  "completion")
    echo "Setting up tab completion for scripts.sh..."
    echo ""
    echo "Option 1 - Source the completion file (recommended):"
    echo "  source .scripts-completion.sh"
    echo ""
    echo "Option 2 - Add to your shell config permanently:"
    SCRIPT_PATH="$(realpath "$0")"

    if [ -n "$ZSH_VERSION" ]; then
      echo "  For zsh, add to ~/.zshrc:"
      echo "  compctl -k \"(install test test-cov format lint docs-build docs-serve docs-clean clean completion)\" $SCRIPT_PATH"
    elif [ -n "$BASH_VERSION" ]; then
      echo "  For bash, add to ~/.bashrc:"
      echo "  complete -W \"install test test-cov format lint docs-build docs-serve docs-clean clean completion\" $SCRIPT_PATH"
    else
      echo "  For zsh: compctl -k \"(install test test-cov format lint docs-build docs-serve docs-clean clean completion)\" $SCRIPT_PATH"
      echo "  For bash: complete -W \"install test test-cov format lint docs-build docs-serve docs-clean clean completion\" $SCRIPT_PATH"
    fi
    echo ""
    echo "Option 3 - Test completion now (temporary):"
    echo "  source .scripts-completion.sh && ./scripts.sh <TAB>"
    ;;
  "install")
    echo "Installing project in development mode..."
    uv sync --extra test --extra docs --extra style
    ;;
  "test")
    echo "Running tests..."
    uv run --extra test pytest
    ;;
  "test-cov")
    echo "Running tests with coverage..."
    uv run --extra test pytest --cov=src/your_package --cov-report=html --cov-report=term
    ;;
  "format")
    echo "Formatting code..."
    uv run --extra style black .
    uv run --extra style ruff check --fix .
    ;;
  "lint")
    echo "Checking code style..."
    uv run --extra style ruff check .
    ;;
  "docs-build")
    echo "Building documentation..."
    uv run --extra docs sphinx-build -b html docs docs/_build/html
    ;;
  "docs-serve")
    echo "Serving documentation..."
    uv run --extra docs sphinx-autobuild docs docs/_build/html --host localhost --port 8000
    ;;
  "docs-clean")
    echo "Cleaning documentation build..."
    rm -rf docs/_build
    ;;
  "clean")
    echo "Cleaning build artifacts..."
    rm -rf build/ dist/ *.egg-info/ .coverage htmlcov/ docs/_build/
    find . -name "*.pyc" -delete
    find . -name "__pycache__" -delete
    ;;
  *)
    echo "Usage: $0 {install|test|test-cov|format|lint|docs-build|docs-serve|docs-clean|clean|completion}"
    echo ""
    echo "Available commands:"
    echo "  install     - Install project in development mode"
    echo "  test        - Run tests"
    echo "  test-cov    - Run tests with coverage"
    echo "  format      - Format code with black and ruff"
    echo "  lint        - Check code style with ruff"
    echo "  docs-build  - Build documentation"
    echo "  docs-serve  - Serve documentation with auto-reload"
    echo "  docs-clean  - Clean documentation build"
    echo "  clean       - Clean all build artifacts"
    echo "  completion  - Set up tab completion for your shell"
    exit 1
    ;;
esac
```

Make it executable:

```bash
chmod +x scripts.sh
```

### 3.3 Create Tab Completion File

Create `.scripts-completion.sh`:

```bash
#!/usr/bin/env bash
# Tab completion for scripts.sh
# Source this file or add to your shell config

_scripts_completion() {
    local cur prev opts
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    opts="install test test-cov format lint docs-build docs-serve docs-clean clean completion"

    COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
    return 0
}

# For bash
if [ -n "$BASH_VERSION" ]; then
    complete -F _scripts_completion ./scripts.sh
    complete -F _scripts_completion scripts.sh
fi

# For zsh
if [ -n "$ZSH_VERSION" ]; then
    compctl -k "(install test test-cov format lint docs-build docs-serve docs-clean clean completion)" ./scripts.sh
    compctl -k "(install test test-cov format lint docs-build docs-serve docs-clean clean completion)" scripts.sh
fi
```

## Step 4: Set Up Sphinx Documentation

### 4.1 Create Documentation Directory Structure

```bash
mkdir -p docs/{getting-started,user-guide,developer-guide,about,api,_static}
```

### 4.2 Sphinx Configuration

Create `docs/conf.py`:

```python
# Configuration file for the Sphinx documentation builder.
import os
import sys
from datetime import date

# Add the src directory to the path
sys.path.insert(0, os.path.abspath("../src"))

import your_package

project = "Your Package Name"
copyright = f"2024-{date.today().year}, Your Name"
author = "Your Name"
release = your_package.__version__

# -- General configuration ---------------------------------------------------

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
    "sphinx_design",
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
pygments_style = "friendly"

html_theme = "furo"
html_static_path = ["_static"]
html_title = f"{project} v{release}"

# Furo theme options
html_theme_options = {
    "sidebar_hide_name": True,
    "navigation_with_keys": True,
    "top_of_page_button": "edit",
    "source_repository": "https://github.com/username/repository/",
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
```

## Step 5: Create Documentation Content

### 5.1 Main Index File

Create `docs/index.rst`:

```rst
Your Package Documentation
==========================

Brief description of your package and its purpose.

.. grid:: 2

   .. grid-item-card:: Getting Started
      :link: getting-started/index
      :link-type: doc

      New to [Package]? Start here for installation and basic usage.

   .. grid-item-card:: User Guide
      :link: user-guide/index
      :link-type: doc

      Learn how to use [Package] for your projects with tutorials and examples.

   .. grid-item-card:: API Reference
      :link: api/index
      :link-type: doc

      Complete API documentation for all modules and functions.

   .. grid-item-card:: Developer Guide
      :link: developer-guide/index
      :link-type: doc

      Contributing guidelines and development information.

Quick Start
-----------

Install with pip:

.. code-block:: bash

   $ python -m pip install your-package

Basic usage:

.. code-block:: python

   import your_package

   # Your quick example here

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   :hidden:

   getting-started/index
   getting-started/installation
   getting-started/quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   :hidden:

   user-guide/index
   user-guide/tutorials
   user-guide/examples
   user-guide/cli

.. toctree::
   :maxdepth: 2
   :caption: Reference
   :hidden:

   api/index
   references
   glossary

.. toctree::
   :maxdepth: 2
   :caption: Development
   :hidden:

   developer-guide/index
   developer-guide/contributing
   developer-guide/changelog
   developer-guide/release-notes

.. toctree::
   :maxdepth: 1
   :caption: About
   :hidden:

   about/authors
   about/license
   about/citation

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
```

### 5.2 Bibliography File

Create `docs/refs.bib`:

```bibtex
@article{example2024,
    title={Example Article},
    author={Author, A.},
    journal={Example Journal},
    year={2024},
    volume={1},
    pages={1--10}
}
```

## Step 6: Documentation Section Templates

I'll create template files for each section that you can customize for any project.

### Getting Started Section

Create `docs/getting-started/index.rst`, `installation.rst`, `quickstart.rst`

### User Guide Section

Create `docs/user-guide/index.rst`, `tutorials.rst`, `examples.rst`, `cli.rst`

### Developer Guide Section

Create `docs/developer-guide/index.rst`, `contributing.rst`, `changelog.rst`

### About Section

Create `docs/about/authors.rst`, `license.rst`, `citation.rst`

### Additional Files

Create `docs/references.rst`, `docs/glossary.rst`, `docs/api/index.rst`

## Step 7: Migration Checklist

- [ ] Update pyproject.toml with new build system
- [ ] Create scripts.sh and make executable
- [ ] Set up tab completion file
- [ ] Update package **init**.py for version handling
- [ ] Create docs directory structure
- [ ] Configure Sphinx with conf.py
- [ ] Create main documentation pages
- [ ] Set up bibliography (refs.bib)
- [ ] Test installation: `./scripts.sh install`
- [ ] Test documentation build: `./scripts.sh docs-build`
- [ ] Test development workflow: `./scripts.sh test`
- [ ] Update CI/CD workflows if needed
- [ ] Update README.md with new instructions

## Step 8: Testing the Migration

```bash
# Test the setup
./scripts.sh install
./scripts.sh test
./scripts.sh docs-build
./scripts.sh docs-serve  # Check in browser

# Test tab completion
source .scripts-completion.sh
./scripts.sh <TAB>  # Should show available commands
```

## Customization Notes

For each project, customize:

1. **Project metadata** in pyproject.toml (name, description, URLs)
2. **Dependencies** based on your project's needs
3. **Package structure** (update paths in conf.py and scripts.sh)
4. **Documentation content** to match your project's domain
5. **CLI commands** if your project has a command-line interface
6. **Bibliography** with relevant references for your field

## Common Issues and Solutions

**Import errors during docs build**: Ensure all dependencies are included in `docs` extra
**Version detection issues**: Verify setuptools-scm configuration and _version.py path
**Tab completion not working**: Source the completion file or add to shell config
**Tests failing**: Update test dependencies and pytest configuration

This process should work for most Python projects. The key is adapting the templates to your specific project's needs while maintaining the standardized structure and tooling.
