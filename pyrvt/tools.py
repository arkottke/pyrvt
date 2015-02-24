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

"""
File: tools.py
Author: Albert Kottke
Description: Tools for reading/writing of files and performing operations.
"""

import csv
import glob
import os
import sys

import numpy as np

from pyrvt.peak_calculators import get_peak_calculator, get_region
from pyrvt import motions

# Try to load the modules required for reading/writing files
try:
    import xlrd
except ImportError:
    xlrd = None

try:
    import xlwt
except ImportError:
    xlwt = None

try:
    import openpyxl
except ImportError:
    openpyxl = None

PARAMETER_NAMES = [
    ('magnitude', 'Magnitude'),
    ('distance', 'Distance (km)'),
    ('vs30', 'Vs30 (m/s)'),
    ('kappa', 'Site Atten., Kappa0 (sec)'),  # Site Atten., κ₀
    ('duration', 'Duration (sec)'),
    ('region', 'Region'),
]


def read_events(fname, response_type='psa'):
    """Read data from the file an Excel work book.

    Parameters
    ----------
    fname : str
        Filename of the input file.
    response_type : str
        Type of response. Possible options:
            psa - (psuedo) spectral accleration [default]
            fa  - Fourier amplitude

    Returns
    -------
    ext : str
        Extension of input file
    reference : numpy.array
        Reference of the response. This is either period [sec] for
        repsonse_type 'psa' or frequency [Hz] for response_type 'fa'
    events : list
        List of events read from the file
    """
    assert response_type in ['psa_target', 'fa']

    ext = os.path.splitext(fname)[1].lower()
    # Load the file depending on the format
    if ext == '.csv':
        def parse(s):
            try:
                return float(s)
            except ValueError:
                return s

        with open(fname) as fp:
            reader = csv.reader(fp)
            rows = [[parse(r) for r in row] for row in reader]
    elif ext == '.xls':
        if xlrd is None:
            raise RuntimeError('xlrd is required to open an xls file')
        wb = xlrd.open_workbook(fname)
        ws = wb.sheet_by_index(0)
        rows = [ws.row_values(i) for i in range(ws.nrows)]
    elif ext == '.xlsx':
        if openpyxl is None:
            raise RuntimeError('openpyxl is required to open an xlsx file')
        wb = openpyxl.load_workbook(fname)
        ws = wb.worksheets[0]
        rows = [[r.value for r in row] for row in ws.rows]
    else:
        raise NotImplementedError

    parameters = {key: rows[i][1:]
                  for i, (key, label) in enumerate(PARAMETER_NAMES)}

    event_row = len(parameters) + 1
    event_count = len(rows[0]) - 1

    reference = np.array([row[0] for row in rows[event_row:]])

    events = []
    for i in range(event_count):
        resps = np.array([row[i + 1] for row in rows[event_row:]])
        # Extract the appropriate attributes
        e = {k: v[i] for k, v in parameters.items()}
        e[response_type] = resps
        events.append(e)

    return ext, reference, events


def write_events(fname, reference, reference_label, response_type,
                 response_label, events):
    """Write the events to a file.

    Parameters
    ----------
    fname : str
        Filename
    reference_target : numpy.array
        Value of the reference (e.g., frequency values)
    reference_label : str
        Label of the reference (e.g., Frequency [Hz])
    response_type : str
        Type of response. Possible options:
            psa - Psuedo spectral accleration [default]
            fa  - Fourier amplitude
    events : list
        Events to write to file
    """

    # Create the rows of output
    rows = []
    # Output the parameters
    for key, label in PARAMETER_NAMES:
        rows.append([label] + [e[key] for e in events])

    rows.append(
        [reference_label] + len(events) * [response_label])

    # Output the response spectra
    for i in range(len(reference)):
        rows.append(
            [reference[i]] + [e[response_type][i] for e in events])

    # Create the directory
    dirname = os.path.dirname(fname)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    # Write the file
    ext = os.path.splitext(fname)[1].lower()

    if ext == '.csv':
        if sys.version_info < (3, 1):
            fp = open(fname, 'wt')
        else:
            fp = open(fname, 'wt', newline='')

        writer = csv.writer(fp)
        writer.writerows(rows)

        fp.close()
    elif ext == '.xls':
        if xlwt is None:
            raise RuntimeError('xlwt is required to open an xls file')

        wb = xlwt.Workbook()
        ws = wb.add_sheet('Sheet 1')

        for i, row in enumerate(rows):
            for j, cell in enumerate(row):
                ws.write(i, j, cell)

        wb.save(fname)
    elif ext == '.xlsx':
        if openpyxl is None:
            raise RuntimeError('openpyxl is required to open an xlsx file')

        wb = openpyxl.Workbook()
        ws = wb.worksheets[0]

        for row in rows:
            ws.append(row)

        wb.save(fname)
    else:
        raise NotImplementedError


def compute_compatible_spectra(method, periods, events, damping):
    """ Compute the response spectrum compatible motions.
    """
    target_freqs = 1. / periods

    for e in events:
        crm = motions.CompatibleRvtMotion(
            target_freqs,
            e['psa_target'],
            damping=damping,
            magnitude=e['magnitude'],
            distance=e['distance'],
            duration=e['duration'],
            region=e['region'],
            peak_calculator=get_peak_calculator(
                method, dict(region=e['region'], mag=e['magnitude'],
                             dist=e['distance']))
        )

        freqs = crm.freqs

        if not e['duration']:
            e['duration'] = crm.duration

        e['fa'] = crm.fourier_amps
        e['psa_calc'] = crm.compute_osc_resp(target_freqs, damping)

    return freqs


def operation_psa2fa(src, dst, damping, method='LP99', fixed_spacing=True):
    '''Compute the acceleration response spectrum from a Fourier amplitude
    spectrum.'''

    for filename_src in glob.iglob(src):
        ext, periods, events = read_events(filename_src, 'psa_target')

        if fixed_spacing:
            # Interpolate the periods to a smaller range
            _periods = np.logspace(-2, 1, 100)

            for e in events:
                e['psa_target'] = np.exp(np.interp(
                    _periods, periods, np.log(e['psa_target'])))

            periods = _periods

        # Compute the FA from the PSA
        freqs = compute_compatible_spectra(method, periods, events,
                                           damping=damping)

        if not os.path.exists(dst):
            os.makedirs(dst)

        basename = os.path.basename(filename_src)
        pathname_dst = os.path.join(dst, basename.rsplit('_', 1)[0])

        write_events(pathname_dst + '_sa' + ext, periods, 'Period (s)',
                     'psa_calc', 'Sa (g)', events)
        write_events(pathname_dst + '_fa' + ext, freqs, 'Frequency (Hz)',
                     'fa', 'FA (g-s)', events)


def operation_fa2psa(src, dst, damping, method='LP99', fixed_spacing=True):
    '''Compute the Fourier amplitude spectrum from a acceleration response
    spectrum.'''

    if fixed_spacing:
        periods = np.logspace(-2, 1, 100)
        osc_freqs = 1. / periods

    for filename_src in glob.iglob(src):
        ext, freqs, events = read_events(filename_src, 'fa')

        if not fixed_spacing:
            osc_freq = freqs
            periods = 1. / osc_freq

        for e in events:
            m = motions.RvtMotion(
                freqs=freqs,
                fourier_amps=e['fa'],
                duration=e['duration'],
                peak_calculator=get_peak_calculator(
                    method, dict(region=e['region'], mag=e['magnitude'],
                                 dist=e['distance']))
            )
            e['sa'] = m.compute_osc_resp(osc_freqs, damping)

        if not os.path.exists(dst):
            os.makedirs(dst)

        basename = os.path.basename(filename_src)

        pathname_dst = os.path.join(dst, basename.rsplit('_', 1)[0])

        write_events(pathname_dst + '_sa' + ext, periods, 'Period (s)', 'sa',
                     'PSA (g)', events)
