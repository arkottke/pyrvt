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

import sys

from .. import runner

from numpy.testing import assert_almost_equal, assert_string_equal, assert_equal

def test_args():
    tests = [
        (
            'test psa2fa -i input -o output -d 0.05 -f -m BJ84',
            [
                ('operation', 'psa2fa'),
                ('input', 'input'),
                ('output', 'output'),
                ('damping', 0.05),
                ('fixed_spacing', True),
                ('method', 'BJ84'),
            ]
        ),
        (
            'test fa2psa --input input --output output --method BT12',
            [
                ('operation', 'fa2psa'),
                ('input', 'input'),
                ('output', 'output'),
                ('damping', 0.05),
                ('fixed_spacing', False),
                ('method', 'BT12'),
            ]
        ),

    ]

    for cmd_line, values in tests:
        yield check_args, cmd_line, values

def check_args(cmd_line, values):
    args = 'test psa2fa -i test -o test -d 0.05 -f -m V75'
    sys.argv = cmd_line.split()

    args = runner.parser.parse_args()

    for name, value in values:
        if name == 'damping':
            f = assert_almost_equal
        elif name == 'fixed_spacing':
            f = assert_equal
        else:
            f = assert_string_equal

        f(value, getattr(args, name))
