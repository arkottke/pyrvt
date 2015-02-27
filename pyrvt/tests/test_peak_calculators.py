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

import os

import numpy as np
from numpy.testing import assert_allclose

from .. import motions
from .. import peak_calculators
from ..tests import readers


def check_osc_accels(peak_calculator, fname_fs, fname_rs):
    path = os.path.join(os.path.dirname(__file__), 'data')
    fs = readers.load_fourier_spectrum(
        os.path.join(path, fname_fs))
    rs = readers.load_rvt_response_spectrum(os.path.join(path, fname_rs))

    rvt_motion = motions.RvtMotion(
        freqs=fs['freqs'],
        fourier_amps=fs['fourier_amps'],
        duration=rs['duration'],
        peak_calculator=peak_calculator)

    spec_accels = rvt_motion.compute_osc_accels(rs['freqs'], rs['damping'])

    # The SMSIM output file only provides 4 significant digits
    assert_allclose(spec_accels, rs['spec_accels'], rtol=0.001, atol=0.05)


def test_boore_joyner_1984():
    check_osc_accels(
        peak_calculators.BooreJoyner1984(),
        'test-bj84.m6.00r0020.0_fs.col',
        'test-bj84.m6.00r020.0_rs.rv.col',
    )


def test_liu_pezeshk_1999():
    check_osc_accels(
        peak_calculators.LiuPezeshk1999(),
        'test-lp99.m6.00r0020.0_fs.col',
        'test-lp99.m6.00r020.0_rs.rv.col',
    )


def test_boore_thompson_2012_wna():
    check_osc_accels(
        peak_calculators.BooreThompson2012('wna', 6, 20.),
        'test-bt12_wna.m6.00r0020.0_fs.col',
        'test-bt12_wna.m6.00r020.0_rs.rv.col',
    )


def test_boore_thompson_2012_ena():
    check_osc_accels(
        peak_calculators.BooreThompson2012('ena', 6, 20.),
        'test-bt12_ena.m6.00r0020.0_fs.col',
        'test-bt12_ena.m6.00r020.0_rs.rv.col',
    )


def check_formulation(method):
    mag = 6.5
    dist = 20
    region = 'cena'

    m = motions.SourceTheoryMotion(
        mag, dist, region,
        peak_calculator=peak_calculators.get_peak_calculator(
            method, dict(mag=mag, dist=dist, region=region))
    )
    m.compute_fourier_amps(np.logspace(-1.5, 2, 1024))

    assert np.any(np.isreal(m.fourier_amps))


def test_peak_factor_formulations():
    methods = ['V75', 'D64', 'DK85', 'TM87']
    for method in methods:
        yield check_formulation, method
