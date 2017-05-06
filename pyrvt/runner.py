#!/usr/bin/python
# -*- coding: utf-8 -*-
"""The command line interface for calling pyRVT.

See :doc:`usage` for more details.
"""

import argparse

from . import __version__
from .tools import operation_psa2fa, operation_fa2psa
from .motions import DEFAULT_CALC

parser = argparse.ArgumentParser(
    prog='pyrvt',
    description='Compute response or Fourier amplitude spectra using RVT.')
parser.add_argument(
    '--version',
    action='version',
    version='%(prog)s version ' + str(__version__))
parser.add_argument(
    'operation',
    help='''Operation to be performed: [psa2fa] converts from
    pseudo-spectral acceleration to Fourier amplitude, and [fa2psa] converts
    from Fourier amplitude to pseudo-spectral acceleration.''',
    choices=['psa2fa', 'fa2psa'])
parser.add_argument(
    '-i',
    '--input',
    help='''Path containing the input file(s). Supported file types are
    csv, xls, and xlsx -- provided the required packages have been
    installed. A single file or glob can be specified. An example of a
    glob would be "input/*_sa.xls" for all files within directory "input"
    ending in "_sa.xls".''',
    required=True)
parser.add_argument(
    '-o',
    '--output',
    help='''Path where the output files should be created. If this
    directory does not exist it will be created. Default: ./output''',
    default='./output')
parser.add_argument(
    '-d',
    '--damping',
    default=0.05,
    type=float,
    help='''Oscillator damping in decimal.  Default: 0.05.''')
parser.add_argument(
    '-f',
    '--fixed-spacing',
    action='store_true',
    help='''Fixed spacing of the oscillator period of
    0.01 to 10 sec log-spaced with 100 points. Target SA values will be
    interpolated if needed''')
parser.add_argument(
    '-m',
    '--method',
    default=DEFAULT_CALC,
    choices=['BJ84', 'BT12', 'DK85', 'LP99', 'TM87', 'V75'],
    help='''Specify the peak factor calculation method. Possible options
    are: [BJ84] Boore and Joyner (1984), [BT12] Boore and Thompson (2012),
    [DK85] Der Kiureghian (1985), [LP99] Liu and Pezeshk (1999), [TM87] Toro
    and McGuire (1987), and [V75] Vanmarcke (1975). If the BT12 method is used,
    then the magnitude, distance and region must be provided by the input
    files. If no value is provided, then '%(default)s' is used as the
    default.''')


def main():
    """Perform the command line operations."""
    args = parser.parse_args()

    if args.operation == 'psa2fa':
        operation_psa2fa(args.input, args.output, args.damping, args.method,
                         args.fixed_spacing)
    elif args.operation == 'fa2psa':
        operation_fa2psa(args.input, args.output, args.damping, args.method,
                         args.fixed_spacing)
    else:
        raise NotImplementedError


if __name__ == '__main__':
    main()
