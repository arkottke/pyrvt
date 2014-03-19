#!/usr/bin/env python3
# encoding: utf-8

import os

from numpy.testing import assert_almost_equal, assert_allclose, assert_array_max_ulp

import pyrvt.motions as motions
import pyrvt.peak_calculators as peak_calculators
import pyrvt.tests.readers as readers

def check_osc_resp(peak_calculator, fname_fs, fname_rs):
    path = os.path.join(os.path.dirname(__file__), 'data')
    fs = readers.load_fourier_spectrum(
        os.path.join(path, fname_fs))
    rs = readers.load_rvt_response_spectrum(
            os.path.join(path, fname_rs))

    rvt_motion = motions.RvtMotion(
        freqs = fs['freqs'],
        fourier_amps = fs['fourier_amps'],
        duration = rs['duration'],
        peak_calculator = peak_calculator)

    spec_accel = rvt_motion.compute_osc_resp(
        rs['freqs'], rs['damping'])

    # The SMSIM output file only provides 4 significant digits
    assert_allclose(spec_accel, rs['spec_accels'], rtol = 0.001, atol=0.05)

def test_boore_joyner_1984():
    check_osc_resp(
        peak_calculators.BooreJoyner1984(),
        'test-bj84.m6.00r0020.0_fs.col',
        'test-bj84.m6.00r020.0_rs.rv.col',
    )

def test_liu_pezeshk_1999():
    check_osc_resp(
        peak_calculators.LiuPezeshk1999(),
        'test-lp99.m6.00r0020.0_fs.col',
        'test-lp99.m6.00r020.0_rs.rv.col',
    )

def test_boore_thompson_2012_wna():
    check_osc_resp(
        peak_calculators.BooreThompson2012('wna', 6, 20.),
        'test-bt12_wna.m6.00r0020.0_fs.col',
        'test-bt12_wna.m6.00r020.0_rs.rv.col',
    )

def test_boore_thompson_2012_ena():
    check_osc_resp(
        peak_calculators.BooreThompson2012('ena', 6, 20.),
        'test-bt12_ena.m6.00r0020.0_fs.col',
        'test-bt12_ena.m6.00r020.0_rs.rv.col',
    )
