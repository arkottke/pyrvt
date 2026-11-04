#!/usr/bin/env python3
"""
Migration automation script for converting projects to uv + Sphinx structure.

This script automates the process of setting up a new project with the
standardized uv + Sphinx documentation structure.

Usage:
    python setup_project.py --name "My Package" --package "my_package" --author "Your Name" --email "your@email.com"

Or run interactively:
    python setup_project.py
"""

import argparse
import sys
from pathlib import Path


def get_template_dir():
    """Get the directory containing the templates."""
    script_dir = Path(__file__).parent
    return script_dir / "templates"


def prompt_for_info():
    """Prompt user for project information."""
    print("Setting up a new Python project with uv + Sphinx")
    print("=" * 50)

    project_name = input("Project display name (e.g., 'My Awesome Package'): ").strip()
    package_name = input("Package name (e.g., 'my_package'): ").strip()
    author_name = input("Author name: ").strip()
    author_email = input("Author email: ").strip()
    repo_url = input("GitHub repository URL (optional): ").strip()
    description = input("Short description: ").strip()

    # Extract username/repo from URL if provided
    username = ""
    repository = ""
    if repo_url:
        try:
            # Extract from https://github.com/username/repo format
            parts = repo_url.strip("/").split("/")
            if len(parts) >= 2:
                username = parts[-2]
                repository = parts[-1]
        except:
            pass

    return {
        "project_name": project_name,
        "package_name": package_name,
        "author_name": author_name,
        "author_email": author_email,
        "repo_url": repo_url,
        "username": username,
        "repository": repository,
        "description": description,
    }


def replace_placeholders(content, info):
    """Replace template placeholders with actual values."""
    replacements = {
        "PROJECT_NAME": info["project_name"],
        "YOUR_PACKAGE_NAME": info["package_name"],
        "YOUR_PACKAGE": info["package_name"],
        "YOUR_NAME": info["author_name"],
        "AUTHOR_NAME": info["author_name"],
        "YOUR_EMAIL": info["author_email"],
        "USERNAME": info["username"],
        "REPOSITORY": info["repository"],
        "YOUR_PACKAGE_DESCRIPTION": info["description"],
        "your-package-name": info["package_name"].replace("_", "-"),
        "your_package": info["package_name"],
    }

    result = content
    for placeholder, value in replacements.items():
        result = result.replace(placeholder, value)

    return result


def copy_and_process_template(template_path, target_path, info):
    """Copy a template file and replace placeholders."""
    with open(template_path, encoding="utf-8") as f:
        content = f.read()

    processed_content = replace_placeholders(content, info)

    target_path.parent.mkdir(parents=True, exist_ok=True)
    with open(target_path, "w", encoding="utf-8") as f:
        f.write(processed_content)

    print(f"Created: {target_path}")


def setup_project_structure(target_dir, info):
    """Set up the complete project structure."""
    target_path = Path(target_dir)
    template_dir = get_template_dir()

    print(f"Setting up project in: {target_path}")

    # Create package directory structure
    src_dir = target_path / "src" / info["package_name"]
    src_dir.mkdir(parents=True, exist_ok=True)

    # Create __init__.py with version import
    init_content = f'''"""{info["description"]}"""

try:
    from ._version import __version__
except ImportError:
    # For development installs
    __version__ = "unknown"

__all__ = ["__version__"]
'''
    with open(src_dir / "__init__.py", "w") as f:
        f.write(init_content)
    print(f"Created: {src_dir / '__init__.py'}")

    # Create tests directory
    tests_dir = target_path / "tests"
    tests_dir.mkdir(exist_ok=True)
    with open(tests_dir / "__init__.py", "w") as f:
        f.write("")
    with open(tests_dir / f"test_{info['package_name']}.py", "w") as f:
        f.write(f'''"""Tests for {info["package_name"]}."""

import pytest
import {info["package_name"]}


def test_import():
    """Test that the package can be imported."""
    assert {info["package_name"]}.__version__ is not None
''')
    print(f"Created: {tests_dir}")

    # Copy and process template files
    template_files = [
        ("pyproject.toml", target_path / "pyproject.toml"),
        ("scripts.sh", target_path / "scripts.sh"),
        (".scripts-completion.sh", target_path / ".scripts-completion.sh"),
        ("conf.py", target_path / "docs" / "conf.py"),
        ("docs/index.rst", target_path / "docs" / "index.rst"),
        (
            "docs/getting-started/index.rst",
            target_path / "docs" / "getting-started" / "index.rst",
        ),
        (
            "docs/getting-started/installation.rst",
            target_path / "docs" / "getting-started" / "installation.rst",
        ),
        (
            "docs/getting-started/quickstart.rst",
            target_path / "docs" / "getting-started" / "quickstart.rst",
        ),
    ]

    for template_file, target_file in template_files:
        template_path = template_dir / template_file
        if template_path.exists():
            copy_and_process_template(template_path, target_file, info)
        else:
            print(f"Warning: Template {template_path} not found")

    # Make scripts.sh executable
    scripts_path = target_path / "scripts.sh"
    if scripts_path.exists():
        scripts_path.chmod(0o755)
        print(f"Made executable: {scripts_path}")

    # Create additional documentation directories
    doc_dirs = [
        "docs/user-guide",
        "docs/developer-guide",
        "docs/about",
        "docs/api",
        "docs/_static",
    ]

    for doc_dir in doc_dirs:
        (target_path / doc_dir).mkdir(parents=True, exist_ok=True)
        print(f"Created directory: {target_path / doc_dir}")

    # Create basic files
    basic_files = {
        "docs/refs.bib": """@article{example2024,
    title={Example Article},
    author={Author, A.},
    journal={Example Journal},
    year={2024},
    volume={1},
    pages={1--10}
}""",
        "docs/glossary.rst": """Glossary
========

.. glossary::

   Term 1
      Definition of term 1.

   Term 2  
      Definition of term 2.
""",
        "README.md": f"""# {info["project_name"]}

{info["description"]}

## Installation

```bash
pip install {info["package_name"].replace("_", "-")}
```

## Quick Start

```python
import {info["package_name"]}

# Your example here
```

## Documentation

Full documentation is available at [Read the Docs](https://{info["package_name"].replace("_", "-")}.readthedocs.io/).

## Development

This project uses [uv](https://docs.astral.sh/uv/) for dependency management.

```bash
# Install development dependencies
./scripts.sh install

# Run tests
./scripts.sh test

# Build documentation
./scripts.sh docs-build
```

## License

MIT License - see LICENSE file for details.
""",
        "requirements.txt": "# This file is for compatibility only. Use pyproject.toml for dependencies.\n",
    }

    for file_path, content in basic_files.items():
        file_path = target_path / file_path
        content = replace_placeholders(content, info)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, "w") as f:
            f.write(content)
        print(f"Created: {file_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Set up a new Python project with uv + Sphinx"
    )
    parser.add_argument("--name", help="Project display name")
    parser.add_argument("--package", help="Package name")
    parser.add_argument("--author", help="Author name")
    parser.add_argument("--email", help="Author email")
    parser.add_argument("--repo", help="GitHub repository URL")
    parser.add_argument("--description", help="Short description")
    parser.add_argument(
        "--target", default=".", help="Target directory (default: current directory)"
    )

    args = parser.parse_args()

    # Get project information
    if all([args.name, args.package, args.author, args.email]):
        info = {
            "project_name": args.name,
            "package_name": args.package,
            "author_name": args.author,
            "author_email": args.email,
            "repo_url": args.repo or "",
            "description": args.description or "",
            "username": "",
            "repository": "",
        }

        # Extract username/repo from URL
        if info["repo_url"]:
            try:
                parts = info["repo_url"].strip("/").split("/")
                if len(parts) >= 2:
                    info["username"] = parts[-2]
                    info["repository"] = parts[-1]
            except:
                pass
    else:
        info = prompt_for_info()

    # Validate required fields
    if not all([info["project_name"], info["package_name"], info["author_name"]]):
        print("Error: Project name, package name, and author name are required")
        sys.exit(1)

    try:
        setup_project_structure(args.target, info)
        print("\n" + "=" * 50)
        print("Project setup complete!")
        print("\nNext steps:")
        print("1. cd to your project directory")
        print("2. Initialize git: git init")
        print("3. Install dependencies: ./scripts.sh install")
        print("4. Run tests: ./scripts.sh test")
        print("5. Build docs: ./scripts.sh docs-build")
        print("6. Set up tab completion: source .scripts-completion.sh")

    except Exception as e:
        print(f"Error setting up project: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
