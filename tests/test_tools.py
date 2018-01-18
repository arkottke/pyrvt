#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import tempfile
import shutil

from numpy.testing import assert_allclose
import pytest

import pyrvt

from . import make_relpath

exts = ['.csv', '.xls', '.xlsx', '.ods']


@pytest.mark.parametrize('ext', exts)
def check_read_events(ext):
    fpath = make_relpath('data', 'test_sa' + ext)
    ext, periods, events = pyrvt.tools.read_events(fpath, 'psa')

    # Test the periods
    assert_allclose(periods, [
        0.010000, 0.010723, 0.011498, 0.012328, 0.013219, 0.014175, 0.015199,
        0.016298, 0.017475, 0.018738, 0.020092, 0.021544, 0.023101, 0.024771,
        0.026561, 0.028480, 0.030539, 0.032745, 0.035112, 0.037649, 0.040370,
        0.043288, 0.046416, 0.049770, 0.053367, 0.057224, 0.061359, 0.065793,
        0.070548, 0.075646, 0.081113, 0.086975, 0.093260, 0.100000, 0.107227,
        0.114976, 0.123285, 0.132194, 0.141747, 0.151991, 0.162975, 0.174753,
        0.187382, 0.200923, 0.215443, 0.231013, 0.247708, 0.265609, 0.284804,
        0.305386, 0.327455, 0.351119, 0.376494, 0.403702, 0.432876, 0.464159,
        0.497702, 0.533670, 0.572237, 0.613591, 0.657933, 0.705480, 0.756463,
        0.811131, 0.869749, 0.932603, 1.000000, 1.072267, 1.149757, 1.232847,
        1.321941, 1.417474, 1.519911, 1.629751, 1.747528, 1.873817, 2.009233,
        2.154435, 2.310130, 2.477076, 2.656088, 2.848036, 3.053856, 3.274549,
        3.511192, 3.764936, 4.037017, 4.328761, 4.641589, 4.977024, 5.336699,
        5.722368, 6.135907, 6.579332, 7.054802, 7.564633, 8.111308, 8.697490,
        9.326033, 10.000000
    ])

    # Test the characteristics of an event
    e = events[0]
    assert e['magnitude'] == 5
    assert e['distance'] == 5
    assert e['vs30'] == 760
    assert_allclose(e['kappa'], 0.039447)
    assert_allclose(e['duration'], 1.361042)
    assert_allclose(e['region'], 'wna')

    # Test the response values of the last event
    e = events[-1]
    assert_allclose(
        e['psa'],
        [
            0.072912, 0.072902, 0.072895, 0.072891, 0.072895, 0.072907,
            0.072921, 0.072937, 0.072954, 0.072973, 0.073151, 0.073831,
            0.074725, 0.075723, 0.076905, 0.078207, 0.079820, 0.081709,
            0.083855, 0.086387, 0.089205, 0.092482, 0.096168, 0.100326,
            0.104895, 0.109822, 0.115028, 0.120412, 0.125902, 0.131879,
            0.137615, 0.142984, 0.148044, 0.152837, 0.157048, 0.160672,
            0.163975, 0.166902, 0.169116, 0.170905, 0.172079, 0.172499,
            0.172268, 0.171136, 0.169345, 0.166969, 0.163935, 0.160346,
            0.156536, 0.152141, 0.147470, 0.142458, 0.137136, 0.131678,
            0.126159, 0.120542, 0.114820, 0.109102, 0.103373, 0.097648,
            0.091909, 0.086209, 0.080628, 0.075285, 0.070170, 0.065238,
            0.060512, 0.056171, 0.052115, 0.048319, 0.044796, 0.041584,
            0.038665, 0.035968, 0.033422, 0.031031, 0.028812, 0.026757,
            0.024839, 0.023024, 0.021317, 0.019658, 0.018069, 0.016560,
            0.015141, 0.013817, 0.012545, 0.011336, 0.010203, 0.009155,
            0.008209, 0.007338, 0.006551, 0.005838, 0.005193, 0.004609,
            0.004082, 0.003605, 0.003173, 0.002780
        ], )


@pytest.mark.parametrize('ext', exts)
def check_write_events(ext):
    # Load the original data
    fpath = make_relpath('data', 'test_sa' + ext)
    ext, periods, events = pyrvt.tools.read_events(fpath, 'psa')

    # Write the data
    handle, dst_fname = tempfile.mkstemp(suffix=ext)
    pyrvt.tools.write_events(dst_fname, periods, 'Period (s)', 'psa', 'SA (g)',
                             events)
    os.close(handle)

    # Reload the data
    _ext, _periods, _events = pyrvt.tools.read_events(dst_fname, 'psa')

    # Delete the temporary file
    os.unlink(dst_fname)

    assert ext == _ext
    assert_allclose(periods, _periods)

    for event, _event in zip(events, _events):
        for key in event:
            if key == 'region':
                assert event[key] == _event[key]
            else:
                assert_allclose(event[key], _event[key])


def test_calc_compatible_spectra():
    src_fname = os.path.join(os.path.dirname(__file__), 'data', 'test_sa.csv')
    ext, periods, events = pyrvt.tools.read_events(src_fname, 'psa')

    pyrvt.tools.calc_compatible_spectra('LP99', periods, events[:1], 0.05)

    # Test that the fit is within 2% of the target
    assert_allclose(events[0]['psa'], events[0]['psa_calc'], rtol=0.02)


@pytest.mark.parametrize('fixed_spacing', [[True, False]])
def test_operation_psa2fa(fixed_spacing):
    src_fname = os.path.join(os.path.dirname(__file__), 'data', 'test_sa.csv')
    dest_dirname = tempfile.mkdtemp()

    # Do not need to check the output as it is checked in
    # test_calc_compatible_spectra
    pyrvt.tools.operation_psa2fa(src_fname, dest_dirname, 0.05, 'LP99',
                                 fixed_spacing)
    shutil.rmtree(dest_dirname)


@pytest.mark.parametrize('fixed_spacing', [[True, False]])
def test_operation_fa2psa(fixed_spacing):
    src_fname = os.path.join(os.path.dirname(__file__), 'data', 'test_fa.csv')
    dest_dirname = tempfile.mkdtemp()

    # Do not need to check the output as it is checked in
    # test_calc_compatible_spectra
    pyrvt.tools.operation_fa2psa(src_fname, dest_dirname, 0.05, 'LP99',
                                 fixed_spacing)
    shutil.rmtree(dest_dirname)
