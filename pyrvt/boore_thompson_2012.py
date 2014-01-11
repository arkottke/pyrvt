#!/usr/bin/env python
# encoding: utf-8

"""Oscillator duration correction for Boore and Thompson (2012).

Parameters files downloaded from D. Boore's website on Jan. 12, 2013 [1]_.

References
----------
Boore, D. M., & Thompson, E. M. (2012). Empirical Improvements for Estimating
Earthquake Response Spectra with Random‐Vibration Theory. Bulletin of the
Seismological Society of America, 102(2), 761-772.

.. [1] http://daveboore.com/rv_improvements_esupp/rv_improvements_esupp.html
"""

import os

import numpy as np
from scipy.interpolate import LinearNDInterpolator

def load_bt12_data(region):
    """Load data from the Boore and Thompson (2012) parameter files.

    Input
    -----
    region : str
        Region for which the parameters were developed. Possible options are:
            'wna' - Western North America (active tectonic)
            'ena' - Eastern North America (stable tectonic)

    Returns
    -------
    np.recarray
        Parameters for the region
    """
    fname = os.path.join(
        os.path.dirname(__file__), 'data',
        region + '_bt12_trms4osc.pars')

    return np.rec.fromrecords(
        np.loadtxt(fname, skiprows=4, usecols=range(9)),
        names='mag,dist,c1,c2,c3,c4,c5,c6,c7')

# Create interpolators for each of the regions
INTERPS = {}
for region in ['wna', 'ena']:
    d = load_bt12_data(region)
    i = LinearNDInterpolator(
        np.c_[d.mag, np.log(d.dist)],
        np.c_[d.c1, d.c2, d.c3, d.c4, d.c5, d.c6, d.c7])
    INTERPS[region] = i

def duration_osc(region, mag, dist, osc_freq, osc_damping, duration_gm):
    """Compute the duration ratio.

    The duration ratio is defined by Equation (10) in Boore and Thompson
    (2012). Magnitude and distance is interpolated.

    Parameters
    ----------
    region : str
        Region for which the parameters were developed. Possible options are:
            'wna' - Western North America (active tectonic)
            'ena' - Eastern North America (stable tectonic)
    mag : float
        Magnitude of the event
    dist : float
        Distance of the event in [km].
    osc_freq : float
        Oscillator frequency in [Hz] -- the reciprocal of the oscillator period.
    osc_damping : float
        Oscillator damping in decimal.
    duration_gm : float
        Duration of the ground motion in [sec].

    Returns
    -------
    duration_ratio : float
        ratio (D_rms/D_ex) of the root-mean-squared duration (D_rms) to the ground motion
        duration (D_ex).

    Notes
    -----
    The interpolant is constructed by triangulating the input data
    with Qhull [1]_, and on each triangle performing linear
    barycentric interpolation.

    References
    ----------
    .. [1] http://www.qhull.org/
    """
    assert region in ['wna', 'ena']

    c1, c2, c3, c4, c5, c6, c7 = INTERPS[region](mag, np.log(dist))

    # Ratio of the oscillator period and the ground motion duration.
    foo = 1 / (osc_freq * duration_gm)

    dur_ratio = ((c1 + c2 * (1 - foo ** c3) / (1 + foo ** c3))
        * (1 + c4 / (2 * np.pi * osc_damping) *
           (foo / (1 + c5 * foo ** c6)) ** c7))

    return duration_gm * dur_ratio
