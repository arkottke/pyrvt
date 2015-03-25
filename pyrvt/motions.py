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
Classes and functions used to define random vibration theory (RVT) based
motions.

References
----------
.. [AH84] Anderson, J. G., & Hough, S. E. (1984). A model for the shape of the
    Fourier amplitude spectrum of acceleration at high frequencies. Bulletin of
    the Seismological Society of America, 74(5), 1969-1993.

.. [AB11] Atkinson, G. M., & Boore, D. M. (2011). Modifications to existing
    ground-motion prediction equations in light of new data. Bulletin of the
    Seismological Society of America, 101(3), 1121-1135.

.. [C03] Campbell, K. W. (2003). Prediction of strong ground motion using
    the hybrid empirical method and its use in the development of ground-motion
    (attenuation) relations in eastern North America.  Bulletin of the
    Seismological Society of America, 93(3), 1012-1033.

.. [GV76] Gasparini, D. A., & Vanmarcke, E. H. (1976). Simulated earthquake
    motions compatible with prescribed response spectra.  Massachusetts
    Institute of Technology, Department of Civil Engineering, Constructed
    Facilities Division.

"""

import numpy as np

from scipy.stats import linregress
from scipy.interpolate import interp1d

from . import peak_calculators

DEFAULT_CALC = 'V75'


def increasing_x(x, y=None):
    """Check if *x* is monotonically increasing. If not, reverse it.

    Parameters
    ----------
    x : :class:`numpy.array`
        X values, which are checked.

    y : :class:`numpy.array` or ``None``, default: ``None``
        Y values, which are reversed if *x* is reversed

    Returns
    -------
    If *y* is not ``None``, then returns (*x*, *y*), otherwise only returns *x*.
    x : :class:`numpy.array`
        X values in monotonically increasing order.

    y : :class:`numpy.array`
        Y values in same order as *x*.

    Raises
    ------
    :class:`NotImplementedError`
        If *x* is not monotonic.
    """

    diffs = np.diff(x)
    if np.all(0 <= diffs):
        # All increasing, do nothing
        pass
    elif np.all(diffs <= 0):
        # All decreasing, reverse
        x = x[::-1]
        if y is not None:
            y = y[::-1]
    else:
        raise NotImplementedError('Values are not regularly ordered.')

    if y is None:
        return x
    else:
        return x, y


def log_spaced_values(min, max):
    """Generate values with constant log-spacing.

    Values are generated with 512 points per decade.

    Parameters
    ----------
    min : float
        Minimum value of the range.

    max : float
        Maximum value of the range.

    Returns
    -------
    value : :class:`numpy.array`
        Log-spaced values

    """
    lower = np.log10(min)
    upper = np.log10(max)
    count = np.ceil(512 * (upper - lower))
    return np.logspace(lower, upper, count)


def compute_sdof_tf(freqs, osc_freq, osc_damping):
    """Compute the single-degree-of-freedom transfer function. When applied on
    the acceleration Fourier amplitude spectrum, it provides the
    psuedo-spectral acceleration.

    Parameters
    ----------
    freqs : :class:`numpy.array`
        Frequencies [Hz] at which the transfer function should be calculated

    osc_freq : float
        Frequency of the oscillator [Hz]

    osc_damping : float
        Damping ratio of the oscillator [decimal]

    Returns
    -------
    tf : numpy.array
        Complex valued transfer function

    """
    return (-osc_freq ** 2. /
            (freqs ** 2 - osc_freq ** 2 -
             2.j * osc_damping * osc_freq * freqs))


def compute_stress_drop(magnitude):
    """Compute the stress drop using [AB11]_ model.

    Parameters
    ----------
    magnitude : float
        Moment magnitude of the stress drop

    Returns
    -------
    stress_drop : float
        Stress drop [bars]

    """
    return 10 ** (3.45 - 0.2 * max(magnitude, 5.))


def compute_geometric_spreading(dist, coefs):
    """Compute the geometric spreading defined by piecewise linear model.

    Parameters
    ----------
    dist : float
        Closest distance to the rupture surface [km]

    coefs : list
        List of (slope, limit) tuples that define the attenuation. For an
        infinite distance use None. For example, [(1, None)] would provide for
        1/R geometric spreading to an infinite distance.

    Returns
    -------
    geometric_spreading : float
        Geometric spreading coefficient

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
    """A :class:`~.motions.RvtMotion` object used to store the properties
    of an RVT motion.

    Parameters
    ----------
    freqs : :class:`numpy.array` or ``None``, default: ``None``
        Frequency array [Hz]

    fourier_amps : :class:`numpy.array` or ``None``, default: ``None``
        Absolute value of acceleration Fourier amplitudes.

    duration : float or ``None``, default: ``None``
        Ground motion duration [dec].

    peak_calculator : str or :class:`~.peak_calculators.Calculator`, default: ``None``
        Peak calculator to use. If ``None``, then the default peak
        calculator is used. The peak calculator may either be specified by a
        :class:`~.peak_calculators.Calculator` object, or by the initials of
        the calculator.

    calc_kwds : dict or ``None``, default: ``None``
        Keywords to be passed during the creation the peak calculator. These
        keywords are only required for some peak calculators.

    """
    def __init__(self, freqs=None, fourier_amps=None, duration=None,
                 peak_calculator=None, calc_kwds=None):
        self._freqs = freqs
        self._fourier_amps = fourier_amps
        self._duration = duration

        if self._freqs is not None:
            self._freqs, self._fourier_amps = increasing_x(
                self._freqs, self._fourier_amps)

        if isinstance(peak_calculator, peak_calculators.Calculator):
            self.peak_calculator = peak_calculator
        else:
            self.peak_calculator = peak_calculators.get_peak_calculator(
                peak_calculator or DEFAULT_CALC, calc_kwds)

    @property
    def freqs(self):
        """Frequency values."""
        return self._freqs

    @property
    def fourier_amps(self):
        """Acceleration Fourier amplitude values."""
        return self._fourier_amps

    @property
    def duration(self):
        """Duration of the ground motion for RVT analysis."""
        return self._duration

    def compute_osc_accels(self, osc_freqs, osc_damping=0.05):
        """Compute the pseudo-acceleration spectral response of an oscillator
        with a specific frequency and damping.

        Parameters
        ----------
        osc_freqs : :class:`numpy.array`
            Natural frequencies of the oscillator [Hz]

        osc_damping : float, default: 0.05
            Damping ratio of the oscillator.

        Returns
        -------
        psa : :class:`numpy.array`
            Peak pseudo-spectral acceleration of the oscillator

        """
        def compute_spec_accel(osc_freq):
            return self.compute_peak(
                compute_sdof_tf(self._freqs, osc_freq, osc_damping),
                osc_freq, osc_damping)

        resp = np.array([compute_spec_accel(of) for of in osc_freqs])
        return resp

    def compute_peak(self, transfer_func=None, osc_freq=None,
                     osc_damping=None):
        """Compute the peak response.

        Parameters
        ----------
        transfer_func : :class:`numpy.array` or ``None``, default: ``None``
            Transfer function to apply to the motion. If ``None``, then no
            transfer function is applied.

        osc_freq : float or ``None``, default: ``None``
            Oscillator frequency for correction of RVT peak calculation.
            Typically, both the *osc_freq* and *osc_damping* need to be
            specified for a correction to be computed.

        osc_damping : float or ``None``, default: ``None``
            Oscillator frequency for correction of RVT peak calculation.

        """
        if transfer_func is None:
            fourier_amps = self._fourier_amps
        else:
            fourier_amps = np.abs(transfer_func) * self._fourier_amps

        return self.peak_calculator(self._duration, self._freqs, fourier_amps,
                                    osc_freq=osc_freq, osc_damping=osc_damping,
                                    full_output=False)

    def compute_attenuation(self, min_freq, max_freq=None, full=False):
        """Compute the site attenuation (κ) based on a log-linear fit.

        This function computes the site attenuation defined by [AH84]_ as:

        .. math::
            a(f) = A_0 \exp(-\pi \kappa f) \\text( for ) f > f_E

        for a single Fourier amplitude spectrum


        Parameters
        ----------
        min_freq : float
            minimum frequency of the fit

        max_freq : float or ``None``, default: ``None``
            maximum frequency of the fit. If ``None``, then the maximum
            frequency range is used.

        full : bool, default: ``False``
            If the complete output should be returned

        Returns
        -------
        If *full* is ``False``, then only *atten* is returned. Otherwise,
        the tuple (*atten*, *r_value*, *freqs*, *fitted*) is return.

        atten : float
            attenuation parameter

        r_sqrt : float
            sqared correlation coefficient of the fit (R²). See
            `scipy.stats.linregress`.

        freqs : :class:`numpy.array`
            selected frequencies

        fitted : :class:`numpy.array`
            fitted values

        """
        max_freq = max_freq or self.freqs[-1]
        mask = (min_freq <= self.freqs) & (self.freqs <= max_freq)

        slope, intercept, r_value, p_value, stderr = linregress(
            self.freqs[mask], np.log(self.fourier_amps[mask]))

        atten = slope / -np.pi

        if full:
            freqs = self.freqs[mask]
            fitted = np.exp(intercept + slope * freqs)
            return atten, r_value ** 2, freqs, fitted
        else:
            return atten


class SourceTheoryMotion(RvtMotion):
    """Single-corner source theory model with default parameters from [C03]_.

    Parameters
    ----------
    magnitude : float
        Moment magnitude of the event

    distance : float
        Epicentral distance [km]

    region : {'cena', 'wna'}, str
        Region for the parameters. Either 'cena' for Central and Eastern
        North America, or 'wna' for Western North America.

    stress_drop : float or None, default: ``None``
        Stress drop of the event [bars]. If ``None``, then the default value is
        used. For *region* = 'cena', the default value is computed by the
        [AB11]_ model, while for *region* = 'wna' the default value is 100
        bars.

    depth : float, default: 8
        Hypocenter depth [km]. The *depth* is combined with the
        *distance* to compute the hypocentral distance.

    peak_calculator : str or :class:`~.peak_calculators.Calculator`, default: ``None``
        Peak calculator to use. If ``None``, then the default peak
        calculator is used. The peak calculator may either be specified by a
        :class:`~.peak_calculators.Calculator` object, or by the initials of
        the calculator.

    calc_kwds : dict or ``None``, default: ``None``
        Keywords to be passed during the creation the peak calculator. These
        keywords are only required for some peak calculators.

    """

    def __init__(self, magnitude, distance, region, stress_drop=None, depth=8,
                 peak_calculator=None, calc_kwds=None):
        super(SourceTheoryMotion, self).__init__(
            peak_calculator=peak_calculator, calc_kwds=calc_kwds)

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
                 4.00, 4.40],
                bounds_error=False)
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
                self.stress_drop = compute_stress_drop(magnitude)

            # Crustal amplification from Campbell (2003) using the
            # log-frequency and the amplification based on a quarter-wave
            # length approximation
            self.site_amp = interp1d(
                np.log([0.01, 0.10, 0.20, 0.30, 0.50, 0.90, 1.25, 1.80, 3.00,
                        5.30, 8.00, 14.00, 30.00, 60.00, 100.00]),
                [1.00, 1.02, 1.03, 1.05, 1.07, 1.09, 1.11, 1.12, 1.13, 1.14,
                 1.15, 1.15, 1.15, 1.15, 1.15],
                bounds_error=False)

        else:
            raise NotImplementedError

        # Depth to rupture
        self.depth = depth
        self.hypo_distance = np.sqrt(self.distance ** 2. + self.depth ** 2.)

        # Constants
        self.seismic_moment = 10. ** (1.5 * (self.magnitude + 10.7))
        self.corner_freq = (4.9e6 * self.shear_velocity *
                            (self.stress_drop / self.seismic_moment) **
                            (1. / 3.))

    def compute_duration(self):
        """Compute the duration by combination of source and path.

        """
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

    def compute_fourier_amps(self, freqs=None):
        """Compute the acceleration Fourier amplitudes for a frequency range.

        Parameters
        ----------
        freqs : :class:`numpy.array` or ``None``, default: ``None``
            Frequency range. If no frequency range is specified then
            `log_spaced_values(0.05, 200.)` is used.

        """
        if freqs is None:
            self._freqs = log_spaced_values(0.05, 200.)
        else:
            self._freqs = increasing_x(np.asarray(freqs))

        self._duration = self.compute_duration()

        # Model component
        const = (0.55 * 2.) / (np.sqrt(2.) * 4. * np.pi * self.density *
                               self.shear_velocity ** 3.)
        source_comp = (const * self.seismic_moment /
                       (1. + (self._freqs / self.corner_freq) ** 2.))

        # Path component
        path_atten = (self.path_atten_coeff * self._freqs **
                      self.path_atten_power)
        geo_atten = compute_geometric_spreading(self.hypo_distance,
                                                self.geometric_spreading)

        path_comp = geo_atten * np.exp(
            (-np.pi * self._freqs * self.hypo_distance) /
            (path_atten * self.shear_velocity))

        # Site component
        site_dim = np.exp(-np.pi * self.site_atten * self._freqs)

        ln_freqs = np.log(self._freqs)
        site_amp = self.site_amp(ln_freqs)
        if np.any(np.isnan(site_amp)):
            # Need to extrapolate
            mask = ln_freqs < self.site_amp.x[0]
            site_amp[mask] = self.site_amp.y[0]

            mask = self.site_amp.x[-1] < ln_freqs
            site_amp[mask] = self.site_amp.y[-1]

        site_comp = site_amp * site_dim

        # Conversion factor to convert from dyne-cm into gravity-sec
        conv = 1.e-20 / 980.7
        # Combine the three components and convert from displacement to
        # acceleration
        self._fourier_amps = (conv * (2. * np.pi * self._freqs) ** 2. *
                              source_comp * path_comp * site_comp)


class CompatibleRvtMotion(RvtMotion):
    """A :class:`~.motions.CompatibleRvtMotion` object is used to compute a
    Fourier amplitude spectrum that is compatible with a target response
    spectrum.

    Parameters
    ----------
    osc_freqs : :class:`numpy.array`
        Frequencies of the oscillator response [Hz].

    osc_accels_target : :class:`numpy.array`
        Spectral acceleration of the oscillator at the specified frequencies
        [g].

    duration : float or None, default: None
        Duration of the ground motion [sec]. If ``None``, then the duration is
        computed using

    osc_damping : float, default: 0.05
        Damping ratio of the oscillator [dec].

    event_kwds : dict or ``None``, default: ``None``
        Keywords passed to :class:`~.motions.SourceTheoryMotion` and used
        to compute the duration of the motion. Only *duration* or
        *event_kwds* should be specified.

    window_len : int or ``None``, default: ``None``
        Window length used for smoothing the computed Fourier amplitude
        spectrum. If ``None``, then no smoothing is applied. The smoothing is
        applied as a moving average with a width of ``window_len``.

    peak_calculator : str or :class:`~.peak_calculators.Calculator`, default: ``None``
        Peak calculator to use. If ``None``, then the default peak
        calculator is used. The peak calculator may either be specified by a
        :class:`~.peak_calculators.Calculator` object, or by the initials of
        the calculator.

    calc_kwds : dict or ``None``, default: ``None``
        Keywords to be passed during the creation the peak calculator. These
        keywords are only required for some peak calculators.

    """
    def __init__(self, osc_freqs, osc_accels_target, duration=None,
                 osc_damping=0.05, event_kwds=None, window_len=None,
                 peak_calculator=None, calc_kwds=None):
        super(CompatibleRvtMotion, self).__init__(
            peak_calculator=peak_calculator)

        osc_freqs, osc_accels_target = increasing_x(
            np.asarray(osc_freqs), np.asarray(osc_accels_target))

        if duration:
            self._duration = duration
        else:
            stm = SourceTheoryMotion(**event_kwds)
            self._duration = stm.compute_duration()

        fourier_amps = self._estimate_fourier_amps(
            osc_freqs, osc_accels_target, osc_damping)

        # The frequency needs to be extended to account for the fact that the
        # oscillator transfer function has a width. The number of frequencies
        # depends on the range of frequencies provided.
        self._freqs = log_spaced_values(osc_freqs[0] / 2.,
                                        2. * osc_freqs[-1])
        self._fourier_amps = np.empty_like(self._freqs)

        # Indices of the first and last point with the range of the provided
        # response spectra
        indices = np.argwhere(
            (osc_freqs[0] < self._freqs) & (self._freqs < osc_freqs[-1]))
        first = indices[0, 0]
        # last is extend one past the usable range to allow use of first:last
        # notation
        last = indices[-1, 0] + 1
        log_freqs = np.log(self._freqs)
        log_osc_freqs = np.log(osc_freqs)

        self._fourier_amps[first:last] = np.exp(np.interp(
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
            self._fourier_amps[0:first] = _extrap(
                self._freqs[0:first],
                self._freqs[first:first+2],
                self._fourier_amps[first:first+2], None)
            # Update the last point using the third- and second-to-last points
            self._fourier_amps[last:] = _extrap(
                self._freqs[last:],
                self._freqs[last-2:last],
                self._fourier_amps[last-2:last], None)

        extrapolate()

        # Apply a ratio correction between the computed at target response
        # spectra
        self.iterations = 0
        self.rmse = 1.

        max_iterations = 30
        tolerance = 5e-6

        osc_accels = self.compute_osc_accels(osc_freqs, osc_damping)

        # Smoothing operator
        if window_len:
            window = np.ones(window_len, 'd')
            window /= window.sum()

        while self.iterations < max_iterations and tolerance < self.rmse:
            # Correct the FAS by the ratio of the target to computed
            # osciallator response. The ratio is applied over the same
            # frequency range. The first and last points in the FAS are
            # determined through extrapolation.

            self._fourier_amps[first:last] *= np.exp(np.interp(
                log_freqs[first:last], log_osc_freqs,
                np.log((osc_accels_target / osc_accels))))

            extrapolate()

            # Apply a running average to smooth the signal
            if window_len:
                self._fourier_amps = np.convolve(
                    window, self._fourier_amps, 'same')

            # Recompute the response spectrum
            osc_accels = self.compute_osc_accels(osc_freqs, osc_damping)

            # Compute the fit between the target and computed oscillator
            # response
            self.rmse = np.sqrt(np.mean(
                (osc_accels_target - osc_accels) ** 2))

            self.iterations += 1

    def _estimate_fourier_amps(self, osc_freqs, osc_accels, osc_damping):
        """Compute an estimate of the FAS using the [GV76]_ methodology.

        Parameters
        ----------
        osc_freqs : :class:`numpy.array`
            Oscillator frequencies in increasing order [Hz]

        osc_accels : :class:`numpy.array`
            Psuedo-spectral accelerations of the oscillator [g].

        osc_damping : float
            Damping ratio of the oscillator.

        Returns
        -------
        fourier_amps : :class:`numpy.array`
            acceleration Fourier amplitude values at the specified
            frequencies specifed by *osc_freqs*.

        """

        # Compute initial value using Vanmarcke methodology. The response is
        # first computed at the lowest frequency and then subsequently computed
        # at higher frequencies.

        peak_factor = 2.5
        fa_sqr_cur = None
        fa_sqr_prev = 0.
        total = 0.

        sdof_factor = np.pi / (4. * osc_damping) - 1.

        fourier_amps = np.empty_like(osc_freqs)

        for i, (osc_freq, osc_accel) in enumerate(
                zip(osc_freqs, osc_accels)):
            # TODO simplify equation and remove duration
            fa_sqr_cur = (
                ((self.duration * osc_accel ** 2) / (2 * peak_factor ** 2) -
                 total) / (osc_freq * sdof_factor))

            if fa_sqr_cur < 0:
                fourier_amps[i] = fourier_amps[i-1]
                fa_sqr_cur = fourier_amps[i] ** 2
            else:
                fourier_amps[i] = np.sqrt(fa_sqr_cur)

            if i == 0:
                total = fa_sqr_cur * osc_freq / 2.
            else:
                total += ((fa_sqr_cur - fa_sqr_prev) /
                          2 * (osc_freq - osc_freqs[i-1]))

        return fourier_amps
