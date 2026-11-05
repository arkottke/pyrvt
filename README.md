# pyRVT

[![PyPI](https://img.shields.io/pypi/v/pyrvt.svg)](https://pypi.org/project/pyrvt/)
[![Python](https://img.shields.io/python/required-version-toml?tomlFilePath=https://raw.githubusercontent.com/arkottke/pyrvt/main/pyproject.toml)](https://pypi.org/project/pyrvt/)
[![Documentation](https://readthedocs.org/projects/pyrvt/badge/?version=latest)](https://pyrvt.readthedocs.io/?badge=latest)
[![Build Status](https://github.com/arkottke/pyrvt/actions/workflows/python-app.yml/badge.svg)](https://github.com/arkottke/pyrvt/actions/workflows/python-app.yml)
[![Code Quality](https://api.codacy.com/project/badge/Grade/4f1fe64804bc45f89b6386666ae20696)](https://app.codacy.com/gh/arkottke/pyrvt/dashboard)
[![Test Coverage](https://api.codacy.com/project/badge/Coverage/4f1fe64804bc45f89b6386666ae20696)](https://app.codacy.com/gh/arkottke/pyrvt/dashboard)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/arkottke/pyrvt/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/5086299.svg)](https://zenodo.org/badge/latestdoi/5086299)

Random vibration theory implemented in Python.

A Python library and command-line application for using random vibration
theory (RVT) to transform between the acceleration Fourier amplitude spectrum
and the acceleration response spectrum.

pyRVT provides multiple peak factor models in a common framework such that they
can be compared and tested. Additionally, it provides an interface to define
RVT-based ground motion models through specification of the Fourier amplitude
spectrum, the acceleration response spectrum, or a seismological model.

## Installation

pyRVT requires Python 3.13 or later and can be installed with
[pip](https://pip.pypa.io):

```bash
python -m pip install pyrvt
```

The command-line application reads and writes spreadsheets through an optional
extra:

```bash
python -m pip install "pyrvt[cli]"
```

To work from the latest source code instead:

```bash
git clone https://github.com/arkottke/pyrvt.git
cd pyrvt
python -m pip install .
```

## Documentation

Information on the installation and usage can be found in the
[documentation](https://pyrvt.readthedocs.io/). The
[changelog](https://github.com/arkottke/pyrvt/blob/main/CHANGELOG.md) lists the
changes in each release.

## Related packages

pyRVT is developed together with:

- [pyGMM](https://github.com/arkottke/pygmm): ground motion, Fourier amplitude
  spectrum and soil curve models
- [pyRotd](https://github.com/arkottke/pyrotd): response spectra and rotated
  response spectra from acceleration time series
- [pyStrata](https://github.com/arkottke/pystrata): site response analysis

## Citation

When citing the software, please reference the
[DOI](https://zenodo.org/badge/latestdoi/5086299).

## License

pyRVT is made available under the MIT License. See
[LICENSE](https://github.com/arkottke/pyrvt/blob/main/LICENSE).
