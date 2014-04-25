#!/usr/bin/env python
# encoding: utf-8

# This file is part of IRVT.
#
# Foobar is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# IRVT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with IRVT.  If not, see <http://www.gnu.org/licenses/>.

'''
File: rvt_operator.py
Author: Albert Kottke
Description: Provides a command line interface for performing RVT calculations.
'''


import argparse

import pyrvt


parser = argparse.ArgumentParser(
    description='Compute response or Fourier amplitude spectra using RVT.')
parser.add_argument(
    'operation',
    help='''Operation to be performed. [psa2fa] converts from
    (pseudo)-spectral acceleration to Fourier amplitude.  [fa2psa] converts
    from Fourier amplitude to (pseudo)-spectral acceleration.''',
    choices=['psa2fa', 'fa2psa'])
parser.add_argument(
    '-i', '--input', dest='src',
    help='''Path containing the input file(s). Supported file types are
    csv, xls, and xlsx -- provided the required packages have been
    installed. A single file or glob can be specified. An example of a
    glob would be "input/*_sa.xls" for all files within directory "input"
    ending in "_sa.xls".''',
    required=True)
parser.add_argument(
    '-o', '--output', dest='dst',
    help='''Path where the output files should be created. If this
    directory does not exist it will be created. Default: ./output''',
    default='./output')
parser.add_argument(
    '-d', '--damping', dest='damping', default=0.05,
    help='''Oscillator damping in decimal.  Default: 0.05.''')
parser.add_argument(
    '-f', '--fixed-spacing', dest='fixed_spacing',
    action='store_true', help='''Fixed spacing of the oscillator period of
    0.01 to 10 sec log-spaced with 100 points. Target SA values will be
    interpolated if needed''')
parser.add_argument(
    '-m', '--method', dest='method',
    help='''Specify the peak factor calculation method. Possible options
    are:
        DK85: Der Kiureghian (1985)
        BJ84: Boore and Joyner (1984)
        LP99: Liu and Pezeshk (1999) [default]
        BT12: Boore and Thompson (2012)
    If the BT12 method is used, then the magnitude, distance and region
    must be provided.'''
)

args = parser.parse_args()

if args.operation == 'sa2fa':
    pyrvt.tools.operation_sa2fa(args.src, args.dst, args.damping,
                                args.method, args.fixed_spacing)
elif args.operation == 'fa2sa':
    pyrvt.tools.operation_fa2sa(args.src, args.dst, args.damping,
                                args.method, args.fixed_spacing)
else:
    raise NotImplementedError
