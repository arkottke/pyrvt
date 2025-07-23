Developer Guide
===============

Welcome to the pyRVT developer documentation! This section contains information
for contributors and developers who want to extend or modify pyRVT.

.. toctree::
   :maxdepth: 2

   contributing
   changelog
   release-notes

Getting Started with Development
--------------------------------

If you're interested in contributing to pyRVT, this guide will help you:

- Set up a development environment
- Understand the code structure
- Learn the contribution workflow
- Follow coding standards and best practices
- Run tests and build documentation

Development Environment Setup
-----------------------------

pyRVT uses modern Python tooling for development:

**Package Management**: `uv <https://docs.astral.sh/uv/>`_ for fast dependency resolution and environment management

**Code Quality**: Black for formatting, Ruff for linting

**Testing**: pytest with coverage reporting

**Documentation**: Sphinx with Furo theme

**Version Control**: Git with conventional commit messages

Quick Development Setup
~~~~~~~~~~~~~~~~~~~~~~~

1. **Clone and enter the repository**:

   .. code-block:: bash

      $ git clone https://github.com/arkottke/pyrvt.git
      $ cd pyrvt

2. **Install uv** (if not already installed):

   .. code-block:: bash

      # On macOS and Linux
      $ curl -LsSf https://astral.sh/uv/install.sh | sh

3. **Set up development environment**:

   .. code-block:: bash

      $ ./scripts.sh install

4. **Verify the setup**:

   .. code-block:: bash

      $ ./scripts.sh test
      $ ./scripts.sh lint

Development Commands
~~~~~~~~~~~~~~~~~~~~

The project includes a ``scripts.sh`` file with common development tasks:

.. code-block:: bash

   $ ./scripts.sh install      # Install in development mode with all extras
   $ ./scripts.sh test         # Run the test suite
   $ ./scripts.sh test-cov     # Run tests with coverage reporting
   $ ./scripts.sh format       # Format code with Black and fix linting issues
   $ ./scripts.sh lint         # Check code style and quality
   $ ./scripts.sh docs-build   # Build documentation
   $ ./scripts.sh docs-serve   # Serve documentation with live reload
   $ ./scripts.sh clean        # Clean build artifacts

Code Structure
--------------

pyRVT is organized into several main modules:

``pyrvt.motions``
   Contains classes for representing ground motions and their transformations.
   Key classes include ``RvtMotion``, ``CompatibleRvtMotion``, and seismological models.

``pyrvt.peak_calculators``
   Implements various peak factor models for random vibration theory calculations.
   Each model is implemented as a separate class with a common interface.

``pyrvt.runner``
   Provides the command-line interface and batch processing functionality.

``pyrvt.tools``
   Contains utility functions for signal processing, numerical integration,
   and other supporting calculations.

Key Design Principles
---------------------

**Modularity**: Each peak factor model is implemented as a separate class with a consistent interface.

**Extensibility**: New peak factor models can be easily added by inheriting from the base class.

**Performance**: Numerical computations use NumPy vectorization and Numba JIT compilation where beneficial.

**Testing**: Comprehensive test suite ensures reliability and prevents regressions.

**Documentation**: All public APIs are documented with examples and mathematical background.

Contributing Workflow
---------------------

We welcome contributions! Please see the :doc:`contributing` guide for detailed information about:

- How to report bugs and request features
- Code style and testing requirements  
- Pull request process
- Documentation requirements

For a quick overview:

1. Fork the repository on GitHub
2. Create a new branch for your changes
3. Make your changes with tests and documentation
4. Run ``./scripts.sh format`` and ``./scripts.sh test``
5. Submit a pull request

Development Best Practices
---------------------------

**Code Style**
   - Use Black for formatting (enforced by CI)
   - Follow PEP 8 guidelines
   - Use type hints for public APIs
   - Write docstrings in NumPy format

**Testing**
   - Write tests for all new functionality
   - Maintain high test coverage (>90%)
   - Test edge cases and error conditions
   - Use pytest fixtures for common setup

**Documentation**
   - Document all public functions and classes
   - Include mathematical background where relevant
   - Provide usage examples
   - Update this documentation for significant changes

**Performance**
   - Profile code before optimizing
   - Use NumPy vectorization for array operations
   - Consider Numba JIT for hot loops
   - Benchmark performance-critical changes

Release Process
---------------

pyRVT follows semantic versioning (MAJOR.MINOR.PATCH):

- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality in a backward-compatible manner  
- **PATCH**: Backward-compatible bug fixes

See :doc:`release-notes` for the detailed release history.
