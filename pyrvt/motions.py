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
File: motions.py
Author: Albert Kottke
Description: Random vibration theory motions.
"""

import enum

import numpy as np

from scipy.interpolate import interp1d

from . import peak_calculators

DEFAULT_CALC = peak_calculators.Vanmarcke1975()


def compute_sdof_tf(freqs, osc_freq, osc_damping):
    """Compute the single-degree-of-freedom transfer function.

    Parameters
    ----------
    freqs : numpy.array
        Frequencies [Hz] at which the transfer function should be calculated.
    osc_freq : float
        Frequency of the oscillator [Hz]
    osc_damping : float
        Damping of the oscillator [decimal].
    Returns
    -------
    numpy.array
        Complex valued ransfer function

    """
    return (-freqs ** 2. /
            (freqs ** 2 - osc_freq ** 2
             - 2.j * osc_damping * osc_freq * freqs))


def compute_stress_drop(magnitude):
    """Compute the stress drop using Atkinson and Boore (2011) model.

    Parameters
    ----------
    magnitude : float
        moment magnitude of the stress drop.

    Returns
    -------
    stress_drop :
        stress drop in bars.

    """
    return 10 ** (3.45 - 0.2 * max(magnitude, 5.))


def compute_geometric_spreading(dist, coefs):
    """Compute the geometric spreading defined by piecewise linear model.

    Parameters
    ----------
    dist : float
        closest distance to the rupture surface [km]

    coefs : list
        list of (slope, limit) tuples that define the attenuation. For an
        inifinite distance use None. For example, [(1, None)] would provide for
        1/R geometric spreading to an infinite distance.

    Returns
    -------
    geometric_spreading : float
        geometric spreading function, Z(R)

    """
    initial = 1
    gs = 1
    for slope, limit in coefs:
        # Compute the distance limited by the maximum distance of the slope.
        _dist = min(dist, limit) if limit else dist
        gs *= (initial / _dist) ** slope

        if _dist < dist:
            initial = _dist
        else:
            break

    return gs


class RvtMotion(object):
    def __init__(self, freqs=None, fourier_amps=None, duration=None,
                 peak_calculator=DEFAULT_CALC):
        self.freqs = freqs
        self.fourier_amps = fourier_amps
        self.duration = duration
        self.peak_calculator = peak_calculator

    def compute_osc_resp(self, osc_freqs, damping=0.05):
        """Compute the psuedo-acceleration spectral response of an oscillator
        with a specific frequency and damping.

        Parameters
        ----------
        osc_freq : numpy.array
            Natural frequency of the oscillator [Hz]
        damping : float (optional)
            Fractional damping of the oscillator.
        resp_type : OscResponse
            Type of osciallator response. Default is
            OscResponse.pseudo_acc
        Returns
        -------
        psa : numpy.array
            peak psuedo spectral acceleration of the oscillator

        """
        # Conversion from acceleration to displacement
        tf_gm = (2.j * np.pi * self.freqs) ** -2

        def compute_spec_accel(fn):
            return self.compute_peak(
                tf_gm * compute_sdof_tf(self.freqs, fn, damping),
                osc_freq=fn, osc_damping=damping)

        resp = np.array([compute_spec_accel(f) for f in osc_freqs])

        # Conversion between spectral displacement and pseudo acceleration.
        resp *= (2 * np.pi * osc_freqs) ** 2

        return resp

    def compute_peak(self, transfer_func=None, osc_freq=None,
                     osc_damping=None):
        """Compute the peak response.

        """
        if transfer_func is None:
            fourier_amps = self.fourier_amps
        else:
            fourier_amps = np.abs(transfer_func) * self.fourier_amps

        return self.peak_calculator(self.duration, self.freqs, fourier_amps,
                                    osc_freq=osc_freq, osc_damping=osc_damping)


class SourceTheoryMotion(RvtMotion):
    """Single-corner source theory model."""
    def __init__(self, magnitude, distance, region,
                 peak_calculator=DEFAULT_CALC,
                 stress_drop=None, depth=8):
        """Compute the duration using the Atkinson and Boore (1995) model.

        Parameters
        ----------
        magnitude : float
            moment magnitude of the event
        distance : float
            distance in km
        stress_drop : float or None
            stress drop of the event in bar. If stress_drop is None, then it
            will be computed using the Atkinson and Boore (2011) model.
        depth : float
            hypocenter depth [km]. Default is 8 km.

        Returns
        -------
        duration : float
            ground motion duration

        """
        super(SourceTheoryMotion, self).__init__(
            peak_calculator=peak_calculator)

        self.magnitude = magnitude
        self.distance = distance
        self.region = peak_calculators.get_region(region)

        if self.region == 'wna':
            # Default parameters for the WUS from Campbell (2003)
            self.shear_velocity = 3.5
            self.path_atten_coeff = 180.
            self.path_atten_power = 0.45
            self.density = 2.8
            self.site_atten = 0.04

            self.geometric_spreading = [(1, 40), (0.5, None)]

            if stress_drop:
                self.stress_drop = stress_drop
            else:
                self.stress_drop = 100.

            # Crustal amplification from Campbell (2003) using the
            # log-frequency and the amplification based on a quarter-wave
            # length approximation
            self.site_amp = interp1d(
                np.log([0.01, 0.09, 0.16, 0.51, 0.84, 1.25, 2.26, 3.17, 6.05,
                        16.60, 61.20, 100.00]),
                [1.00, 1.10, 1.18, 1.42, 1.58, 1.74, 2.06, 2.25, 2.58, 3.13,
                 4.00, 4.40])
        elif self.region == 'cena':
            # Default parameters for the CEUS from Campbell (2003)
            self.shear_velocity = 3.6
            self.density = 2.8
            self.path_atten_coeff = 680.
            self.path_atten_power = 0.36
            self.site_atten = 0.006

            self.geometric_spreading = [(1, 70), (0, 130), (0.5, None)]

            if stress_drop:
                self.stress_drop = stress_drop
            else:
                self.stress_drop = 150.

            # Crustal amplification from Campbell (2003) using the
            # log-frequency and the amplification based on a quarter-wave
            # length approximation
            self.site_amp = interp1d(
                np.log([0.01, 0.10, 0.20, 0.30, 0.50, 0.90, 1.25, 1.80, 3.00,
                        5.30, 8.00, 14.00, 30.00, 60.00, 100.00]),
                [1.00, 1.02, 1.03, 1.05, 1.07, 1.09, 1.11, 1.12, 1.13, 1.14,
                 1.15, 1.15, 1.15, 1.15, 1.15])

        else:
            raise NotImplementedError

        # Depth to rupture
        self.depth = 8.
        self.hypo_distance = np.sqrt(self.distance ** 2. + self.depth ** 2.)

        # Constants
        self.seismic_moment = 10. ** (1.5 * (self.magnitude + 10.7))
        self.corner_freq = (4.9e6 * self.shear_velocity
                            * (self.stress_drop
                               / self.seismic_moment) ** (1./3.))

    def compute_duration(self):
        # Source component
        duration_source = 1. / self.corner_freq

        # Path component
        if self.region == 'wna':
            duration_path = 0.05 * self.hypo_distance
        elif self.region == 'cena':
            duration_path = 0.

            if 10 < self.hypo_distance:
                # 10 < R <= 70 km
                duration_path += 0.16 * (min(self.hypo_distance, 70) - 10.)

            if 70 < self.hypo_distance:
                # 70 < R <= 130 km
                duration_path += -0.03 * (min(self.hypo_distance, 130) - 70.)

            if 130 < self.hypo_distance:
                # 130 km < R
                duration_path += 0.04 * (self.hypo_distance - 130.)
        else:
            raise NotImplementedError

        return duration_source + duration_path

    def compute_fourier_amps(self, freqs):
        self.freqs = np.asarray(freqs)
        self.duration = self.compute_duration()

        # Model component
        const = (0.55 * 2.) / (np.sqrt(2.) * 4. * np.pi * self.density
                               * self.shear_velocity ** 3.)
        source_comp = (const * self.seismic_moment
                       / (1. + (self.freqs / self.corner_freq) ** 2.))

        # Path component
        path_atten =\
            self.path_atten_coeff * self.freqs ** self.path_atten_power
        geo_atten = compute_geometric_spreading(self.hypo_distance,
                                                self.geometric_spreading)

        path_comp = geo_atten * np.exp(
            (-np.pi * self.freqs * self.hypo_distance)
            / (path_atten * self.shear_velocity))

        # Site component
        site_dim = np.exp(-np.pi * self.site_atten * self.freqs)
        site_comp = self.site_amp(np.log(self.freqs)) * site_dim

        # Conversion factor to convert from dyne-cm into gravity-sec
        conv = 1.e-20 / 980.7
        # Combine the three components and convert from displacement to
        # acceleleration
        self.fourier_amps = (conv * (2. * np.pi * self.freqs) ** 2.
                             * source_comp * path_comp * site_comp)


class CompatibleRvtMotion(RvtMotion):
    def __init__(self, osc_freqs, osc_resp_target, duration=None, damping=0.05,
                 magnitude=None, distance=None, stress_drop=None, region=None,
                 window_len=None, peak_calculator=DEFAULT_CALC):
        """Compute a Fourier amplitude spectrum that is compatible with a
        target response spectrum."""
        super(CompatibleRvtMotion, self).__init__(
            peak_calculator=peak_calculator)

        osc_freqs = np.asarray(osc_freqs)
        osc_resp_target = np.asarray(osc_resp_target)

        # Order by increasing frequency
        ind = osc_freqs.argsort()
        osc_freqs = osc_freqs[ind]
        osc_resp_target = osc_resp_target[ind]

        if duration:
            self.duration = duration
        else:
            stm = SourceTheoryMotion(magnitude, distance, region, stress_drop)
            self.duration = stm.compute_duration()

        # Compute initial value using Vanmarcke methodology. The response is
        # first computed at the lowest frequency and then subsequently comptued
        # at higher frequencies.

        peak_factor = 2.5
        fa_sqr_cur = None
        fa_sqr_prev = 0.
        total = 0.

        sdof_factor = np.pi / (4. * damping) - 1.

        fourier_amps = np.empty_like(osc_freqs)

        for i, (osc_freq, osc_resp) in enumerate(
                zip(osc_freqs, osc_resp_target)):
            fa_sqr_cur = (
                ((self.duration * osc_resp ** 2) / (2 * peak_factor ** 2) -
                 total) / (osc_freq * sdof_factor))

            if fa_sqr_cur < 0:
                fourier_amps[i] = fourier_amps[i-1]
                fa_sqr_cur = fourier_amps[i] ** 2
            else:
                fourier_amps[i] = np.sqrt(fa_sqr_cur)

            if i == 0:
                total = fa_sqr_cur * osc_freq / 2.
            else:
                total += ((fa_sqr_cur - fa_sqr_prev)
                          / 2 * (osc_freq - osc_freqs[i-1]))

        # The frequency needs to be extended to account for the fact that the
        # osciallator transfer function has a width.
        self.freqs = np.logspace(
            np.log10(osc_freqs[0] / 2.),
            np.log10(2 * osc_freqs[-1]), 1024)
        self.fourier_amps = np.empty_like(self.freqs)

        # Indices of the first and last point with the range of the provided
        # response spectra
        indices = np.argwhere(
            (osc_freqs[0] < self.freqs) & (self.freqs < osc_freqs[-1]))
        first = indices[0, 0]
        # last is extend one past the usable range to allow use of first:last
        # notation
        last = indices[-1, 0] + 1

        log_freqs = np.log(self.freqs)
        log_osc_freqs = np.log(osc_freqs)

        self.fourier_amps[first:last] = np.exp(np.interp(
            log_freqs[first:last], log_osc_freqs, np.log(fourier_amps)))

        def extrapolate():
            # Extrapolate the first and last value of Fourier amplitude
            # spectrum.
            def _extrap(freq, freqs, fourier_amps, max_slope=None):
                # Extrapolation is performed in log-space using the first and
                # last two points
                xi = np.log(freq)
                x = np.log(freqs)
                y = np.log(fourier_amps)

                slope = (y[1] - y[0]) / (x[1] - x[0])

                if max_slope:
                    slope = min(slope, max_slope)

                return np.exp(slope * (xi - x[0]) + y[0])

            # Update the first point using the second and third points
            self.fourier_amps[0:first] = _extrap(
                self.freqs[0:first],
                self.freqs[first:first+2],
                self.fourier_amps[first:first+2], None)
            # Update the last point using the third- and second-to-last points
            self.fourier_amps[last:] = _extrap(
                self.freqs[last:],
                self.freqs[last-2:last],
                self.fourier_amps[last-2:last], None)

        extrapolate()

        # Apply a ratio correction between the computed at target response
        # spectra
        self.iterations = 0
        self.rmse = 1.

        max_iterations = 30
        tolerance = 5e-6

        osc_resps = self.compute_osc_resp(osc_freqs, damping)

        if window_len:
            window = np.ones(window_len, 'd')
            window /= window.sum()

        while self.iterations < max_iterations and tolerance < self.rmse:
            # Correct the FAS by the ratio of the target to computed
            # osciallator response. The ratio is applied over the same
            # frequency range. The first and last points in the FAS are
            # determined through extrapolation.

            self.fourier_amps[first:last] *= np.exp(np.interp(
                log_freqs[first:last], log_osc_freqs,
                np.log((osc_resp_target / osc_resps))))

            extrapolate()

            # Apply a running average to smooth the signal
            if window_len:
                self.fourier_amps = np.convolve(
                    window, self.fourier_amps, 'same')

            # Compute the fit between the target and computed oscillator
            # response
            self.rmse = np.sqrt(np.mean(
                (osc_resp_target
                 - self.compute_osc_resp(osc_freqs, damping)) ** 2))

            self.iterations += 1

            # Recompute the response spectrum
            osc_resps = self.compute_osc_resp(osc_freqs, damping)
