Changelog
=========

This changelog documents all notable changes to pyRVT.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

[Unreleased]
------------

- Documentation reorganization with improved structure
- Migration from Hatch to uv for dependency management
- Enhanced Sphinx documentation with Furo theme

[0.7.3] - 2024-XX-XX
---------------------

Added
~~~~~
- Enhanced peak factor models
- Improved error handling and validation
- Better support for modern Python versions

Changed
~~~~~~~
- Updated dependencies to latest versions
- Improved performance of numerical computations

Fixed
~~~~~
- Various bug fixes and stability improvements

[0.7.2] - 2024-XX-XX
---------------------

Added
~~~~~
- Additional peak calculator implementations
- Extended test coverage

Fixed
~~~~~
- Edge case handling in peak factor calculations

[0.7.1] - 2024-XX-XX  
---------------------

Changed
~~~~~~~
- Documentation improvements
- Code quality enhancements

Fixed
~~~~~
- Minor bug fixes

[0.7.0] - 2024-XX-XX
---------------------

Added
~~~~~
- Support for Python 3.10+
- New peak factor models
- Enhanced command-line interface

Changed
~~~~~~~
- **Breaking**: Minimum Python version increased to 3.10
- Updated API for improved consistency
- Performance optimizations

Removed
~~~~~~~
- Support for Python 3.9 and earlier

Migration Guide
---------------

**From 0.6.x to 0.7.x**

1. **Python Version**: Upgrade to Python 3.10 or later
2. **Dependencies**: Update all dependencies with ``pip install --upgrade pyrvt``
3. **API Changes**: Most APIs remain the same, check deprecation warnings
4. **Testing**: Run your existing code to verify compatibility

**From 0.5.x to 0.6.x**

Major changes included enhanced peak factor models and improved documentation.
Most user-facing APIs remained stable.

See the complete changelog in the `GitHub repository <https://github.com/arkottke/pyrvt/releases>`_
for detailed release notes and technical changes.
