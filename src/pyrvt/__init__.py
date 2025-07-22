#!/usr/bin/python
"""pyRVT module. Provides random vibration theory models for earthquake ground motion
calculations.
"""

from importlib.metadata import version

from . import motions, peak_calculators, tools

__all__ = [
    "tools",
    "motions",
    "peak_calculators",
]

__author__ = "Albert Kottke"
__copyright__ = "Copyright 2016-2025 Albert Kottke"
__license__ = "MIT"
__title__ = "pyRVT"

try:
    from ._version import version as __version__
except ImportError:
    # Fallback for development installs without setuptools-scm
    try:
        from importlib.metadata import version

        __version__ = version("pyRVT")
    except ImportError:
        __version__ = "unknown"
