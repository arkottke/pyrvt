#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib

import numpy as np
import pyexcel
import pytest

import pyrvt
import pysra

from numpy.testing import assert_allclose, assert_string_equal

from . import readers

fpath_data = pathlib.Path(__file__).parent / 'data'


@pytest.mark.parametrize('peak_calculator,abbrev,suffix', [
    (pyrvt.peak_calculators.BooreJoyner1984(), 'bj84', 'm6.00r020.0'),
    (pyrvt.peak_calculators.LiuPezeshk1999(), 'lp99', 'm6.00r020.0'),
    (pyrvt.peak_calculators.BooreThompson2012('wna', 6, 20.), 'bt12_wna',
     'm6.00r020.0'),
    (pyrvt.peak_calculators.BooreThompson2012('ena', 6, 20.), 'bt12_ena',
     'm6.00r020.0'),
    (pyrvt.peak_calculators.BooreThompson2015('wna', 6, 21.3), 'bt15_wna',
     'm6.00rps021.3'),
    (pyrvt.peak_calculators.BooreThompson2015('ena', 6, 20.8), 'bt15_ena',
     'm6.00rps020.8'),
])
def test_osc_accels(peak_calculator, abbrev, suffix):
    fs = readers.load_fourier_spectrum(
        str(fpath_data / f'test-{abbrev}.{suffix}_fs.col'))
    rs = readers.load_rvt_response_spectrum(
        str(fpath_data / f'test-{abbrev}.{suffix}_rs.rv.col'))

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


@pytest.mark.parametrize(
    'method', ['V75', 'D64', 'DK85', 'TM87', 'BT12', 'BT15', 'WR18'])
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


def read_wang_rathje_18_data(motion_id):
    mag = {0: 5, 1: 6.5, 2: 8}[motion_id]

    wb = pyexcel.get_book(file_name=str(fpath_data / 'wang_rathje_2018.xlsx'))
    # Input Fourier amplitude
    ws = wb['FAS and Dgm']
    freqs = np.array(ws.column[0][2:1002])
    fourier_ampls = np.array(ws.column[1 + motion_id][2:1002])
    duration = ws[2, 5 + motion_id]
    motion = pysra.motion.RvtMotion(
        freqs, fourier_ampls, duration,
        pyrvt.peak_calculators.WangRathje2018('cena', mag, 20)
    )

    expected = {}
    for key, sheetname in [('rock', 'SaRock'),
                           ('surface', 'SaSurf (Modified)')]:
        ws = wb[sheetname]
        osc_freq = np.array(ws.column[0][2:303])
        spec_acc = np.array(ws.column[1 + motion_id][2:303])
        expected[key] = np.rec.fromarrays(
            [osc_freq, spec_acc], names='osc_freq,spec_acc')

    return motion, expected


@pytest.mark.parametrize('motion_id', [0, 1, 2])
@pytest.mark.parametrize('location', ['rock', 'surface'])
def test_wang_rathje(motion_id, location):
    motion, expected = read_wang_rathje_18_data(motion_id)

    if location == 'rock':
        actual = motion.calc_osc_accels(expected[location].osc_freq)
    elif location == 'surface':
        calc = pysra.propagation.LinearElasticCalculator()
        profile = pysra.site.Profile([
            pysra.site.Layer(
                pysra.site.SoilType('Soil', 18., mod_reduc=None, damping=0.01),
                316, 400),
            pysra.site.Layer(
                pysra.site.SoilType('Rock', 22., mod_reduc=None, damping=0.01),
                0, 3000),
        ])
        calc(motion, profile, profile.location('outcrop', index=-1))
        rso = pysra.output.ResponseSpectrumOutput(
            expected[location].osc_freq,
            pysra.output.OutputLocation('outcrop', depth=0), 0.05)
        # Compute the response spectrum at the surface
        rso(calc)
        actual = rso.values

    else:
        raise NotImplementedError

    assert_allclose(actual, expected[location].spec_acc, rtol=0.02)
