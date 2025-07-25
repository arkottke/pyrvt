# pyRVT

[![Docs](https://readthedocs.org/projects/pyrvt/badge/?version=latest)](https://pyrvt.readthedocs.io/)
[![PyPi Cheese
Shop](https://img.shields.io/pypi/v/pyrvt.svg)](https://img.shields.io/pypi/v/pyrvt.svg)
[![Documentation](https://readthedocs.org/projects/pyrvt/badge/?version=latest)](https://pyrvt.readthedocs.io/?badge=latest)
[![Build status](https://github.com/arkottke/pyrvt/actions/workflows/python-app.yml/badge.svg)](https://github.com/arkottke/pyrvt/actions/workflows/python-app.yml)
[![Code
Quality](https://api.codacy.com/project/badge/Grade/4f1fe64804bc45f89b6386666ae20696)](https://www.codacy.com/manual/arkottke/pyrvt)
[![Test
Coverage](https://api.codacy.com/project/badge/Coverage/4f1fe64804bc45f89b6386666ae20696)](https://www.codacy.com/manual/arkottke/pyrvt)
![License](https://img.shields.io/badge/license-MIT-blue.svg)
[![DOI](https://zenodo.org/badge/5086299.svg)](https://zenodo.org/badge/latestdoi/5086299)

A Python library and command-line application for using random vibration
theory to transform between acceleration Fourier amplitude spectrum and
acceleration response spectrum.

Information on the installation and usage can be found in the
[documentation](https://pyrvt.readthedocs.io/).

`pyRVT` provides random vibration theory (RVT) models for use in earthquake
ground motion models. It provides multiple peak factor models in a common
framework such that they can be compared and tested. Additionally, it provides
an interface to define RVT based ground motion models through specification of
the Fourier amplitude spectrum, acceleration response spectrum, or calculated by
seismological models.

## Installation

`pyRVT` can be installed with [pip](https://pip.pypa.io):

```bash
$ python -m pip install pyrvt
```

Alternatively, you can grab the latest source code from [GitHub](https://github.com/arkottke/pyrvt):

```
$ git clone https://github.com/arkottke/pyrvt.git
$ cd pyrvt
$ pip install .
```

## Citation

When citing the software reference the [DOI](https://zenodo.org/records/3630729).

## License

`pyRVT` is made available under the MIT License.
