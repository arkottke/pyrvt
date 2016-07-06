#!/usr/bin/python
# -*- coding: utf-8 -*-

import pyrvt.runner


def test_args_short():
    cmd_line = 'psa2fa -i input -o output -d 0.05 -f -m BJ84'
    args = pyrvt.runner.parser.parse_args(args=cmd_line.split())

    assert args.operation == 'psa2fa'
    assert args.input == 'input'
    assert args.output == 'output'
    assert args.damping == 0.05
    assert args.fixed_spacing is True
    assert args.method == 'BJ84'


def test_args_long():
    cmd_line = 'fa2psa --input input --output output --method BT12'
    args = pyrvt.runner.parser.parse_args(args=cmd_line.split())

    assert args.operation == 'fa2psa'
    assert args.input == 'input'
    assert args.output == 'output'
    assert args.damping == 0.05
    assert args.method == 'BT12'
