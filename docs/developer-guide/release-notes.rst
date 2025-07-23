Release Notes
=============

This page contains detailed release notes for pyRVT versions.

For a complete changelog, see :doc:`changelog`.

Major Releases
--------------

Version 0.7.x Series
~~~~~~~~~~~~~~~~~~~~~

**New Features:**
- Enhanced peak factor models with improved accuracy
- Support for modern Python versions (3.10+)
- Improved command-line interface
- Better error handling and validation

**Breaking Changes:**
- Minimum Python version increased to 3.10
- Some API changes in peak calculator interface

**Performance Improvements:**
- Faster numerical computations using updated NumPy
- Optimized memory usage for large datasets

Version 0.6.x Series
~~~~~~~~~~~~~~~~~~~~~

**New Features:**
- Additional peak factor models
- Enhanced documentation with examples
- Improved test coverage

**Bug Fixes:**
- Fixed edge cases in peak factor calculations
- Resolved issues with very short or long durations

Migration Guide
---------------

Upgrading from 0.6.x to 0.7.x
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Python Version Requirements:**
pyRVT 0.7.x requires Python 3.10 or later. If you're using an older Python version,
you'll need to upgrade your Python installation.

**API Changes:**
Most user-facing APIs remain the same. The main changes are:

- Some internal peak calculator methods have changed signatures
- Improved error messages and validation

**Installation:**
Simply update your installation:

.. code-block:: bash

   $ pip install --upgrade pyrvt

**Testing Your Code:**
After upgrading, run your existing code and tests to ensure compatibility.
Most code should work without changes.

Deprecation Policy
------------------

pyRVT follows these deprecation practices:

1. **Advance Notice**: Deprecated features are marked as deprecated for at least one minor release before removal
2. **Clear Warnings**: Deprecation warnings are issued when deprecated features are used
3. **Migration Path**: Alternative approaches are always provided before deprecating features
4. **Documentation**: Deprecations are clearly documented in release notes and API documentation

Security Updates
-----------------

Security-related updates are handled with high priority:

- Critical security issues are addressed in patch releases
- Security advisories are published through GitHub Security Advisories
- Users are notified through multiple channels (GitHub, PyPI, documentation)

Getting Support
---------------

For questions about specific releases:

1. Check the :doc:`changelog` for detailed changes
2. Review the API documentation for updated interfaces
3. Search or create issues on `GitHub <https://github.com/arkottke/pyrvt/issues>`_
4. Join discussions in the GitHub Discussions section

Planned Features
----------------

**Upcoming in 0.8.x:**
- Additional ground motion models
- Enhanced plotting capabilities  
- Performance optimizations
- Extended file format support

**Future Releases:**
- Integration with other seismological tools
- Web-based interface
- Extended validation against published results

Note: Feature plans are subject to change based on community feedback and development priorities.
