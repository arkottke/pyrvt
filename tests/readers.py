#!/usr/bin/python
# -*- coding: utf-8 -*-
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
        fourier_amps=rows[:, 8], )


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
        spec_accels=rows[:, 11])
