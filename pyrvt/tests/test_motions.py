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
