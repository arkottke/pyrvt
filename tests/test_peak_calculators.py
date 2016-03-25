#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

import numpy as np
from numpy.testing import assert_allclose, assert_equal
import pytest

import pyrvt

from . import readers


def check_osc_accels(peak_calculator, fname_fs, fname_rs):
    path = os.path.join(os.path.dirname(__file__), 'data')
    fs = readers.load_fourier_spectrum(
        os.path.join(path, fname_fs))
    rs = readers.load_rvt_response_spectrum(os.path.join(path, fname_rs))

    rvt_motion = pyrvt.motions.RvtMotion(
        freqs=fs['freqs'],
        fourier_amps=fs['fourier_amps'],
        duration=rs['duration'],
        peak_calculator=peak_calculator)

    spec_accels = rvt_motion.compute_osc_accels(rs['freqs'], rs['damping'])

    # The SMSIM output file only provides 4 significant digits
    assert_allclose(spec_accels, rs['spec_accels'], rtol=0.001, atol=0.05)


def test_boore_joyner_1984():
    bj84 = pyrvt.peak_calculators.BooreJoyner1984()
    check_osc_accels(
        bj84,
        'test-bj84.m6.00r0020.0_fs.col',
        'test-bj84.m6.00r020.0_rs.rv.col',
    )

    assert_equal(bj84.name, 'Boore & Joyner (1984)')
    assert_equal(bj84.abbrev, 'BJ84')


def test_liu_pezeshk_1999():
    check_osc_accels(
        pyrvt.peak_calculators.LiuPezeshk1999(),
        'test-lp99.m6.00r0020.0_fs.col',
        'test-lp99.m6.00r020.0_rs.rv.col',
    )


def test_boore_thompson_2012_wna():
    check_osc_accels(
        pyrvt.peak_calculators.BooreThompson2012('wna', 6, 20.),
        'test-bt12_wna.m6.00r0020.0_fs.col',
        'test-bt12_wna.m6.00r020.0_rs.rv.col',
    )


def test_boore_thompson_2012_ena():
    check_osc_accels(
        pyrvt.peak_calculators.BooreThompson2012('ena', 6, 20.),
        'test-bt12_ena.m6.00r0020.0_fs.col',
        'test-bt12_ena.m6.00r020.0_rs.rv.col',
    )


@pytest.mark.parametrize('method', ['V75', 'D64', 'DK85', 'TM87'])
def test_formulations(method):
    mag = 6.5
    dist = 20
    region = 'cena'

    m = pyrvt.motions.SourceTheoryMotion(
        mag, dist, region,
        peak_calculator=pyrvt.peak_calculators.get_peak_calculator(
            method, dict(mag=mag, dist=dist, region=region))
    )
    m.calc_fourier_amps(np.logspace(-1.5, 2, 1024))

    osc_freqs = np.logspace(-1, 2, num=50)
    osc_accels = m.compute_osc_accels(osc_freqs, 0.05)

    assert np.all(np.isreal(osc_accels))
