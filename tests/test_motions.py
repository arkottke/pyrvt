#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from numpy.testing import assert_allclose

import pyrvt


def test_calc_attenuation():
    m = pyrvt.motions.SourceTheoryMotion(5.5, 0, 'cena', depth=1)
    m.calc_fourier_amps()

    atten, r_value, freqs, fitted = m.calc_attenuation(50)

    assert_allclose(0.006, atten, rtol=0.01)
    assert_allclose(1.0, r_value, rtol=0.01)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(m.freqs, m.fourier_amps, 'b-')
    # ax.set_xlabel('Frequency (Hz)')
    # ax.set_xscale('log')
    # ax.set_ylabel('Amplitude')
    # ax.set_yscale('log')
    #
    # fig.savefig('test')


def test_compatible_rvt_motion():
    # Compute the target from the point source model.
    target = pyrvt.motions.SourceTheoryMotion(
        6.,
        20.,
        'wna',
        peak_calculator=pyrvt.peak_calculators.DerKiureghian1985())
    target.calc_fourier_amps(np.logspace(-1.5, 2, 1024))

    osc_freqs = np.logspace(-1, 2, num=50)
    osc_accels_target = target.calc_osc_accels(osc_freqs, 0.05)

    compat = pyrvt.motions.CompatibleRvtMotion(
        osc_freqs,
        osc_accels_target,
        duration=target.duration,
        osc_damping=0.05,
        peak_calculator=pyrvt.peak_calculators.DerKiureghian1985())

    osc_accels_compat = compat.calc_osc_accels(osc_freqs, 0.05)

    # Might be off by a few percent because of difficulties with the inversion.
    assert_allclose(osc_accels_target, osc_accels_compat, rtol=0.03, atol=0.05)

    # fig, axes = plt.subplots(2, 1)
    #
    # ax = axes.flat[0]
    #
    # ax.plot(target.freqs, target.fourier_amps, 'b-', label='Target')
    # ax.plot(compat.freqs, compat.fourier_amps, 'r--', label='Compatible')
    #
    # ax.set_xlabel('Frequency (Hz)')
    # ax.set_xscale('log')
    # ax.set_ylabel('Fourier Ampl. (cm/s)')
    # ax.set_yscale('log')
    #
    # ax = axes.flat[1]
    #
    # ax.plot(osc_freqs, osc_resp_target, 'b-', label='Target')
    # ax.plot(osc_freqs, osc_resp_compat, 'r--', label='Compatible')
    #
    # ax.set_xlabel('Frequency (Hz)')
    # ax.set_xscale('log')
    # ax.set_ylabel('Spectral Accel. (cm/s²)')
    # ax.set_yscale('log')
    #
    # fig.savefig('compatible_fas.png', dpi=300)
