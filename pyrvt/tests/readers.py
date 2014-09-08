#!/usr/bin/env python3
# encoding: utf-8

# pyRVT: Seismological random vibration theory implemented with Python
# Copyright (C) 2013-2014 Albert R. Kottke albert.kottke@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Scripts for loading test case data."""

import numpy as np


def load_fourier_spectrum(fname):
    """Load Fourier amplitude spectrum file created by fas_drvr.exe.

    Inputs
    ------
    fname : string
        File name of the data file

    Returns
    -------
    event : dict
        Dictionary containing the event
    """
    assert fname.endswith('_fs.col')
    rows = np.loadtxt(fname, skiprows=2)

    return dict(
        mag=rows[0, 0],
        dist=rows[0, 1],
        freqs=rows[:, 4],
        fourier_amps=rows[:, 8],
    )


def load_rvt_response_spectrum(fname):
    """Load response spectrum file created by fas_drvr.exe.

    Inputs
    ------
    fname : string
        File name of the data file

    Returns
    -------
    event : dict
        Dictionary containing the event
    """
    assert fname.endswith('_rs.rv.col')
    rows = np.loadtxt(fname, skiprows=2)

    return dict(
        damping=rows[0, 0],
        mag=rows[0, 3],
        dist=rows[0, 4],
        duration=rows[0, 16],
        freqs=rows[:, 2],
        spec_accels=rows[:, 11]
    )
