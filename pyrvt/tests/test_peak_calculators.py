
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

import numpy as np

from .. import motions, peak_calculators


def test_peak_factor_formulations():
    methods = ['V75', 'D64', 'DK85', 'TM87', 'CLH56', 'BJ84', 'LP99', 'BT12']
    for method in methods:
        yield check_formulation, method


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