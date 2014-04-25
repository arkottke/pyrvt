#!/usr/bin/env python3
# encoding: utf-8

import os

import numpy as np
from numpy.testing import assert_allclose

import pyrvt.motions as motions
import pyrvt.peak_calculators as peak_calculators
import pyrvt.tests.readers as readers


def check_osc_resp(peak_calculator, fname_fs, fname_rs):
    path = os.path.join(os.path.dirname(__file__), 'data')
    fs = readers.load_fourier_spectrum(
        os.path.join(path, fname_fs))
    rs = readers.load_rvt_response_spectrum(os.path.join(path, fname_rs))

    rvt_motion = motions.RvtMotion(
        freqs=fs['freqs'],
        fourier_amps=fs['fourier_amps'],
        duration=rs['duration'],
        peak_calculator=peak_calculator)

    spec_accel = rvt_motion.compute_osc_resp(
        rs['freqs'], rs['damping'])

    # The SMSIM output file only provides 4 significant digits
    assert_allclose(spec_accel, rs['spec_accels'], rtol=0.001, atol=0.05)


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


def test_compatible_rvt_motion():
    # Compute the target from the point source model.
    target = motions.SourceTheoryMotion(
        6., 20., 'wna',
        peak_calculator=peak_calculators.DerKiureghian1985())
    target.compute_fourier_amps(np.logspace(-1.5, 2, 1024))

    osc_freqs = np.logspace(-1, 2)
    osc_resp_target = target.compute_osc_resp(osc_freqs, damping=0.05)

    print(osc_resp_target)

    compat = motions.CompatibleRvtMotion(
        osc_freqs, osc_resp_target,
        duration=target.duration, damping=0.05,
        peak_calculator=peak_calculators.DerKiureghian1985())

    osc_resp_compat = compat.compute_osc_resp(osc_freqs, damping=0.05)

    # Might be off by a few percent because of difficulties with the inversion.
    assert_allclose(osc_resp_target, osc_resp_compat, rtol=0.03, atol=0.05)
