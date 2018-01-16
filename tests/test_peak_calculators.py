#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib

import numpy as np
from numpy.testing import assert_allclose, assert_string_equal
import pytest

import pyrvt

from . import readers


@pytest.mark.parametrize('peak_calculator,abbrev,suffix', [
    (pyrvt.peak_calculators.BooreJoyner1984(), 'bj84', 'm6.00r020.0'),
    (pyrvt.peak_calculators.LiuPezeshk1999(), 'lp99', 'm6.00r020.0'),
    (pyrvt.peak_calculators.BooreThompson2012('wna', 6, 20.),
     'bt12_wna', 'm6.00r020.0'),
    (pyrvt.peak_calculators.BooreThompson2012('ena', 6, 20.),
     'bt12_ena', 'm6.00r020.0'),
    (pyrvt.peak_calculators.BooreThompson2015('wna', 6, 21.3),
     'bt15_wna', 'm6.00rps021.3'),
    (pyrvt.peak_calculators.BooreThompson2015('ena', 6, 20.8),
     'bt15_ena', 'm6.00rps020.8'),
])
def test_osc_accels(peak_calculator, abbrev, suffix):
    fpath = pathlib.Path(__file__).parent / 'data'
    fs = readers.load_fourier_spectrum(
        str(fpath / f'test-{abbrev}.{suffix}_fs.col')
    )
    rs = readers.load_rvt_response_spectrum(
        str(fpath / f'test-{abbrev}.{suffix}_rs.rv.col')
    )

    rvt_motion = pyrvt.motions.RvtMotion(
        freqs=fs['freqs'],
        fourier_amps=fs['fourier_amps'],
        duration=rs['duration'],
        peak_calculator=peak_calculator)

    spec_accels = rvt_motion.calc_osc_accels(rs['freqs'], rs['damping'])

    # The SMSIM output file only provides 4 significant digits
    assert_allclose(spec_accels, rs['spec_accels'], rtol=0.05, atol=0.01)


@pytest.fixture
def bj84_pc():
    return pyrvt.peak_calculators.BooreJoyner1984()


def test_name(bj84_pc):
    assert_string_equal(bj84_pc.name, 'Boore & Joyner (1984)')


def test_abbrev(bj84_pc):
    assert_string_equal(bj84_pc.abbrev, 'BJ84')


@pytest.mark.parametrize('method', ['V75', 'D64', 'DK85', 'TM87', 'BT12',
                                    'BT15', 'WR18'])
def test_formulations(method):
    mag = 6.5
    dist = 20
    region = 'cena'

    m = pyrvt.motions.SourceTheoryMotion(
        mag,
        dist,
        region,
        peak_calculator=pyrvt.peak_calculators.get_peak_calculator(
            method, dict(mag=mag, dist=dist, region=region)))
    m.calc_fourier_amps(np.logspace(-1.5, 2, 1024))

    osc_freqs = np.logspace(-1, 2, num=50)
    osc_accels = m.calc_osc_accels(osc_freqs, 0.05)
    assert np.all(np.isreal(osc_accels))
