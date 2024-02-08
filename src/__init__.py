"""pyRVT module."""

from pkg_resources import get_distribution

from . import motions
from . import peak_calculators
from . import tools

__all__ = [
    "tools",
    "motions",
    "peak_calculators",
]

__author__ = "Albert Kottke"
__copyright__ = "Copyright 2016-2024 Albert Kottke"
__license__ = "MIT"
__title__ = "pyRVT"
__version__ = get_distribution("pyRVT").version
