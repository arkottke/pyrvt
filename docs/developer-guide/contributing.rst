Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/arkottke/pyrvt/issues.

If you are reporting a bug, please include:

- Your operating system name and version.
- Any details about your local setup that might be helpful in troubleshooting.
- Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

pyRVT could always use more documentation, whether as part of the
official pyRVT docs, in docstrings, or even on the web in blog posts,
articles, and such. Docstrings should be formatted using the
`NumPy conventions <https://numpydoc.readthedocs.io/en/latest/format.html>`_

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/arkottke/pyrvt/issues.

If you are proposing a feature:

- Explain in detail how it would work.
- Keep the scope as narrow as possible, to make it easier to implement.
- Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up ``pyRVT`` for local development.

Prerequisites
~~~~~~~~~~~~~

This project uses `uv <https://docs.astral.sh/uv/>`_ for dependency management. Install uv first:

.. code-block:: bash

   # On macOS and Linux
   curl -LsSf https://astral.sh/uv/install.sh | sh

   # Or with pip
   pip install uv

Setup
~~~~~

1. Fork the ``pyRVT`` repo on GitHub.
2. Clone your fork locally:

   .. code-block:: bash

      $ git clone git@github.com:your_name_here/pyrvt.git

3. Create a branch for local development:

   .. code-block:: bash

      $ git checkout -b name-of-your-bugfix-or-feature

4. Install the development environment using uv:

   .. code-block:: bash

      $ ./scripts.sh install

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass formatting and the
   tests:

   .. code-block:: bash

      $ ./scripts.sh format
      $ ./scripts.sh test

   The documentation can be built and served with:

   .. code-block:: bash

      $ ./scripts.sh docs-build
      $ ./scripts.sh docs-serve

6. Commit your changes and push your branch to GitHub:

   .. code-block:: bash

      $ git add .
      $ git commit -m "Your detailed description of your changes."
      $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.md.
3. The pull request should work for Python 3.10 and later.

Development Commands
--------------------

This project uses [uv](https://docs.astral.sh/uv/) for dependency management and a custom ``scripts.sh`` file for common development tasks.

Available development commands:

.. code-block:: bash

   $ ./scripts.sh install      # Install project in development mode
   $ ./scripts.sh test         # Run tests
   $ ./scripts.sh test-cov     # Run tests with coverage
   $ ./scripts.sh format       # Format code with black and ruff
   $ ./scripts.sh lint         # Check code style with ruff
   $ ./scripts.sh docs-build   # Build documentation
   $ ./scripts.sh docs-serve   # Serve documentation with auto-reload
   $ ./scripts.sh docs-clean   # Clean documentation build
   $ ./scripts.sh clean        # Clean all build artifacts

To run a subset of tests:

.. code-block:: bash

   $ uv run pytest tests/test_specific_module.py

Code Style and Quality
----------------------

**Formatting**: Code is automatically formatted using Black and Ruff.

**Type Hints**: Use type hints for public APIs and complex functions.

**Documentation**: All public functions and classes should have NumPy-style docstrings.

**Testing**: Write tests for new functionality. Aim for high test coverage.

**Performance**: Profile code before optimizing. Use NumPy vectorization and consider Numba for hot loops.

Release Process
---------------

Releases are handled by the maintainers. The process includes:

1. Update version numbers
2. Update changelog
3. Create GitHub release
4. Publish to PyPI
5. Update documentation

Community Guidelines
--------------------

- Be respectful and inclusive
- Follow the code of conduct
- Help newcomers get started
- Provide constructive feedback
- Celebrate contributions from all skill levels

Getting Help
------------

If you need help with development:

1. Check the existing documentation
2. Search GitHub issues for similar problems
3. Ask questions in GitHub Discussions
4. Contact the maintainers

Thank you for contributing to pyRVT!
