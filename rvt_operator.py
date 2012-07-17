#!/usr/bin/python

import argparse
import glob
import os

import numpy as np
import xlrd
import xlwt

import rvt

damping = 0.05
parameter_names = [
    ('magnitude', 'Magnitude'),
    ('distance', 'Distance (km)'),
    ('vs30', 'Vs30 (m/s)'),
    ('kappa', 'Kappa0 (sec)'),
    ('duration', 'Duration (sec)'),
        ]


def load_events(filename, response_key='sa'):
    '''Load an Excel work book'''

    start_rowx = len(parameter_names) + 1

    wb = xlrd.open_workbook(filename)
    ws = wb.sheet_by_index(0)

    parameters = {key: ws.row_values(i, start_colx=1)
            for i, (key, label) in enumerate(parameter_names)}

    reference = np.array(ws.col_values(0, start_rowx))

    events = []
    for i in range(ws.ncols - 1):
        sa = np.array(ws.col_values(i + 1, start_rowx))

        # Extract the appropriate attributes
        e = {k: v[i] for k, v in parameters.iteritems() if not k == 'duration'}
        e[response_key] = sa

        if parameters['duration'][i] in ['wus', 'ceus']:
            e['region'] = parameters['duration'][i]
            e['duration'] = None
        else:
            e['duration'] = float(parameters['duration'][i])
            e['region'] = None

        events.append(e)

    return reference, events


def export_events(filename, reference, reference_label, response_key,
        response_label, events):
    wb = xlwt.Workbook()
    ws = wb.add_sheet('Sheet 1')
    row = 0

    for key, label in parameter_names:
        ws.write(row, 0, label)

        for column, e in enumerate(events):
            ws.write(row, column + 1, e[key])

        row += 1

    ws.write(row, 0, reference_label)

    for column in range(len(events)):
        ws.write(row, column + 1, response_label)

    row += 1
    for i in range(len(reference)):
        ws.write(row + i, 0, reference[i])

        for column, e in enumerate(events):
            ws.write(row + i, 1 + column, e[response_key][i])

    wb.save(filename)


def compute_compatible_spectra(period, events, damping=0.05):
    ''' Compute the response spectrum compatible motions. '''
    target_freq = 1. / period

    for e in events:
        crm = rvt.CompatibleRvtMotion(target_freq, e['sa_target'],
                damping=damping, magnitude=e['magnitude'],
                distance=e['distance'], duration=e['duration'],
                region=e['region'])

        freq = crm.freq

        if not e['duration']:
            e['duration'] = crm.duration

        e['fa'] = crm.fourier_amp
        e['sa_calc'] = crm.compute_osc_resp(target_freq, damping)

    return freq


def operation_sa2fa(src, dest, damping, fixed_spacing):
    '''Compute the acceleration response spectrum from a Fourier amplitude
    spectrum.'''

    for filename_src in glob.iglob(src):
        if filename_src[-3:] != 'xls':
            continue

        period, events = load_events(filename_src, 'sa_target')

        if fixed_spacing:
            # Interpolate the periods to a smaller range
            _period = np.logspace(-2, 1, 100)

            for e in events:
                e['sa_target'] = np.exp(np.interp(
                    _period, period, np.log(e['sa_target'])))

            period = _period

        # Compute the FA from the PSA
        freq = compute_compatible_spectra(period, events, damping=damping)

        if not os.path.exists(dest):
            os.makedirs(dest)

        basename = os.path.basename(filename_src)
        pathname_dest = os.path.join(dest, basename.rsplit('_', 1)[0])

        export_events(pathname_dest + '_sa.xls', period, 'Period (s)',
                'sa_calc', 'Sa (g)', events)
        export_events(pathname_dest + '_fa.xls', freq, 'Frequency (Hz)', 'fa',
                'FA (g-s)', events)


def operation_fa2sa(src, dest, damping, fixed_spacing):
    '''Compute the Fourier amplitude spectrum from a acceleration response
    spectrum.'''

    if fixed_spacing:
        period = np.logspace(-2, 1, 100)
        osc_freq = 1. / period

    for filename_src in glob.iglob(src):
        if filename_src[-3:] != 'xls':
            continue

        freq, events = load_events(filename_src, 'fa')

        if not fixed_spacing:
            osc_freq = freq
            period = 1. / osc_freq

        for e in events:
            m = rvt.RvtMotion(freq=freq, fourier_amp=e['fa'],
                    duration=e['duration'])
            e['sa'] = m.compute_osc_resp(osc_freq, damping)

        if not os.path.exists(dest):
            os.makedirs(dest)

        basename = os.path.basename(filename_src)

        pathname_dest = os.path.join(dest, basename.rsplit('_', 1)[0])

        export_events(pathname_dest + '_sa.xls', period, 'Period (s)', 'sa',
                'Sa (g)', events)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute reponse or Fourier amplitude spectra using RVT.')
    parser.add_argument('operation',
            help='Operation to be performed. '
            '[sa2fa] converts from (psuedo)-spectral acceleration to Fourier amplitude. '
            '[fa2sa] converts from Fourier amplitude to (psuedo)-spectral acceleration.',
            choices=['sa2fa', 'fa2sa'])
    parser.add_argument('-i', '--input', dest='src',
                    help='Path containing the input files. Single file or glob'
                    'can be specified. An example of a glob would be'
                    '"input/*_sa.xls" for all files within directory "input"'
                    'ending in "_sa.xls".', required=True)
    parser.add_argument('-o', '--output', dest='dest',
                    help='Path where the output files should be created. If this'
                    'directory does not exist it will be created. Default: ./output',
                    default='./output')
    parser.add_argument('-d', '--damping', dest='damping', default=0.05,
            help='Oscillator damping in decimal.  Default: 0.05.')
    parser.add_argument('-f', '--fixed-spacing', dest='fixed_spacing', action='store_true',
            help='Fixed spacing of the oscillator period '
            ' of 0.01 to 10 sec log-spaced with 100 points. Target SA values'
            ' will be interpolated if needed')

    args = parser.parse_args()

    if args.operation == 'sa2fa':
        operation_sa2fa(args.src, args.dest, args.damping, args.fixed_spacing)
    elif args.operation == 'fa2sa':
        operation_fa2sa(args.src, args.dest, args.damping, args.fixed_spacing)
