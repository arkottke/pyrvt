#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys

from numpy.testing import (assert_almost_equal, assert_string_equal,
                           assert_equal)

import pyrvt.runner


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

    args = pyrvt.runner.parser.parse_args()

    for name, value in values:
        if name == 'damping':
            f = assert_almost_equal
        elif name == 'fixed_spacing':
            f = assert_equal
        else:
            f = assert_string_equal

        f(value, getattr(args, name))
