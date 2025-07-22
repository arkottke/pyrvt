# Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

## Types of Contributions

### Report Bugs

Report bugs at https://github.com/arkottke/pyrvt/issues.

If you are reporting a bug, please include:

- Your operating system name and version.
- Any details about your local setup that might be helpful in troubleshooting.
- Detailed steps to reproduce the bug.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

### Write Documentation

pyRVT could always use more documentation, whether as part of the
official pyRVT docs, in docstrings, or even on the web in blog posts,
articles, and such. Docstrings should be formatted using the
[NumPy conventions](https://numpydoc.readthedocs.io/en/latest/format.html)

### Submit Feedback

The best way to send feedback is to file an issue at https://github.com/arkottke/pyrvt/issues.

If you are proposing a feature:

- Explain in detail how it would work.
- Keep the scope as narrow as possible, to make it easier to implement.
- Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

## Get Started!

Ready to contribute? Here's how to set up `pyRVT` for local development.

### Prerequisites

This project uses [uv](https://docs.astral.sh/uv/) for dependency management. Install uv first:

```bash
# On macOS and Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Or with pip
pip install uv
```

### Setup

1. Fork the `pyRVT` repo on GitHub.
2. Clone your fork locally

```
$ git clone git@github.com:your_name_here/pyrvt.git
```

3. Create a branch for local development::

```
$ git checkout -b name-of-your-bugfix-or-feature
```

4. Install the development environment using uv:

```
$ ./scripts.sh install
```

Now you can make your changes locally.

5. When you're done making changes, check that your changes pass formatting and the
   tests.

```
$ ./scripts.sh format
$ ./scripts.sh test
```

The documentation can be built and served with:

```
$ ./scripts.sh docs-build
$ ./scripts.sh docs-serve
```

6. Commit your changes and push your branch to GitHub::

```
$ git add .
$ git commit -m "Your detailed description of your changes."
$ git push origin name-of-your-bugfix-or-feature
```

7. Submit a pull request through the GitHub website.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.md.
3. The pull request should work for Python 3.7 and later.

## Tips

This project uses [uv](https://docs.astral.sh/uv/) for dependency management and a custom `scripts.sh` file for common development tasks.

Available development commands:

```
$ ./scripts.sh install      # Install project in development mode
$ ./scripts.sh test         # Run tests
$ ./scripts.sh test-cov     # Run tests with coverage
$ ./scripts.sh format       # Format code with black and ruff
$ ./scripts.sh lint         # Check code style with ruff
$ ./scripts.sh docs-build   # Build documentation
$ ./scripts.sh docs-serve   # Serve documentation with auto-reload
$ ./scripts.sh docs-clean   # Clean documentation build
$ ./scripts.sh clean        # Clean all build artifacts
```

To run a subset of tests:

```
$ uv run pytest tests/test_specific_module.py
```
