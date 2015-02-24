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
File: peak_calculators.py
Author: Albert Kottke
Description: Published peak factor models, which compute the expected peak
ground motion. A specific model may include oscillator duration correction.
"""

import os

import numpy as np

from scipy.integrate import quad
from scipy.interpolate import LinearNDInterpolator


def compute_moments(freqs, fourier_amps, orders):
    """Compute the spectral moments.

    The spectral moment is computed using the squared Fourier amplitude
    spectrum.

    Parameters
    ----------
    freqs : numpy.array
        Frequency of the Fourier amplitude spectrum [Hz]

    fourier_amps : numpy.array
        Amplitude of the Fourier amplitude spectrum [g-s]

    Returns
    -------
    moments : list
        Spectral moments.

    """
    squared_fa = np.square(fourier_amps)

    # Use trapzoid integration to compute the requested moments.
    moments = [2. * np.trapz(
        np.power(2 * np.pi * freqs, o) * squared_fa, freqs)
        for o in orders]

    return moments


class Calculator(object):
    NAME = ''
    ABBREV = ''

    def __init__(self, **kwargs):
        pass

    @property
    def name(self):
        """Name of the calculator."""
        return self.NAME

    @property
    def abbrev(self):
        """Abbreviated name of the calculator."""
        return self.ABBREV

    def __call__(self, peak_factor, resp_rms, full_output):
        """Return the computed peak response."""
        resp_peak = resp_rms * peak_factor

        if full_output:
            return resp_peak, peak_factor
        else:
            return resp_peak

    @classmethod
    def limited_num_zero_crossings(cls, num_zero_crossings):
        return max(1.33, num_zero_crossings)


class Vanmarcke1975(Calculator):
    """Vanmarcke (1975) [#]_ peak factor which includes the effects of clumping.

    The peak factor equation is from equation (2) in (Der Kiureghian, 1980)
    [#]_ , which is based on equation (29) in (Vanmarcke, 1975) .

    The cumulative density function (CDF) of the peak is defined as:

    .. math::
        F_x(x) = \\left[1 - \\exp\\left(-x^2/2\\right)\\right]
        \\exp\\left[-N_z \\frac{1 -
            \\exp\\left(-\\sqrt{\\pi/2} \\delta_e x\\right)}{\\exp(x^2 / 2) -
            1 }\\right]

    where :math:`N_z` is the number of zero crossings, :math:`delta_e` is the
    effective bandwidth (:math:`\delta ^ 1.2`).

    Typically, the expected value of the peak factor is calculated using the
    probability density function (:math:`f_x(x) = \\frac{d}{dx} F_x(x)`) by:

    .. math::
        E[x] = \\int_0^\\infty x f_x(x) dx

    However, because of the properties of :math:`F_x(x)`, specifically that it
    has non-zero probablities for only positive values, :math:`E[x]` can be
    computed directly from :math:`F_x(x)`.

    .. math::
        E[x] = \\int_0^\\infty 1 - F_x(x) dx.

    This is based on the following references [#]_ and [#]_.


    References
    ----------
    .. [#] Der Kiureghian, A. (1980). Structural response to stationary
        excitation. Journal of the Engineering Mechanics Division, 106(6),
        1195-1213.

    .. [#] Vanmarcke, E. H. (1975). On the distribution of the first-passage
        time for normal stationary random processes. Journal of applied
        mechanics, 42(1), 215-220.

    .. [#] http://en.wikipedia.org/wiki/Expected_value#Formulas_for_special_cases

    .. [#] http://stats.stackexchange.com/a/13377/48461

    """
    NAME = 'Vanmarcke (1975)'
    ABBREV = 'V75'

    def __init__(self, **kwargs):
        super(Vanmarcke1975, self).__init__(**kwargs)

    def __call__(self, gm_duration, freqs, fourier_amps, osc_freq, osc_damping,
                 full_output=False, **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        gm_duration : float
            Duration of the strong-motion phase of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Aris
            intensity [sec]
        freqs : numpy.array
            Frequency of the Fourier amplitude spectrum [Hz]
        fourier_amps: numpy.array
            Amplitude of the Fourier amplitude spectrum with a single
            degree of freedom oscillator already applied if being used. Units
            are not important.
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Damping of the oscillator [decimal]. For example, 0.05 for 5%.
        full_output : bool
            (optional) If the full output should be returned. Default is
            ``False``.

        Returns
        -------
        max_resp : float
            Expected maximum response
        peak_factor : float
            Associated peak factor. Only provided in `full_output` is True.

        """

        m0, m1, m2 = compute_moments(freqs, fourier_amps, [0, 1, 2])

        # Compute the root-mean-squared response
        resp_rms = np.sqrt(m0 / gm_duration)

        bandwidth = np.sqrt(1 - (m1 * m1) / (m0 * m2))
        bandwidth_eff = bandwidth ** 1.2

        num_zero_crossings = self.limited_num_zero_crossings(
            gm_duration * np.sqrt(m2 / m0) / np.pi)

        def ccdf(x):
            """ The expected peak factor is computed as the integral of the
            complementary CDF (1 - CDF(x)).

            """
            return (1 - (1 - np.exp(-x ** 2 / 2)) *
                    np.exp(-1 * num_zero_crossings *
                           (1 - np.exp(-1 * np.sqrt(np.pi / 2) *
                                       bandwidth_eff * x)) /
                           (np.exp(x ** 2 / 2) - 1)))

        peak_factor = quad(ccdf, 0, np.inf)[0]

        if osc_freq and osc_damping:
            peak_factor *= self.nonstationarity_factor(
                osc_damping, osc_freq, gm_duration)

        return super(Vanmarcke1975, self).__call__(
            peak_factor, resp_rms, full_output)

    @classmethod
    def nonstationarity_factor(cls, osc_damping, osc_freq,  gm_duration):
        return np.sqrt(
            1 - np.exp(-4 * np.pi * osc_damping * osc_freq * gm_duration))


class Davenport1964(Calculator):
    """RVT calculation using the asymptotic solution proposed by Davenport
    (1964) [#]_.

    References
    ----------
    .. [#] Davenport, A. G. (1964). Note on the distribution of the largest
        value of a random function with application to gust loading. In
        Institute of Civil Engineering Proceedings, 28(2), 187-196.

    """

    NAME = 'Davenport (1964)'
    ABBREV = 'D64'

    def __init__(self, **kwargs):
        super(Davenport1964, self).__init__(**kwargs)

    def __call__(self, gm_duration, freqs, fourier_amps, full_output=False,
                 **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        gm_duration : float
            Duration of the strong-motion phase of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity [sec]
        freqs : numpy.array
            Frequency of the Fourier amplitude spectrum [Hz]
        fourier_amps: numpy.array
            Amplitude of the Fourier amplitude spectrum with a single
            degree of freedom oscillator already applied if being used. Units
            are not important.
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Damping of the oscillator [decimal]. For example, 0.05 for 5%.
        full_output : bool
            (optional) If the full output should be returned

        Returns
        -------
        max_resp : float
            Expected maximum response
        peak_factor : float
            Associated peak factor. Only provided in `full_output` is True.

        """

        m0, m2 = compute_moments(freqs, fourier_amps, [0, 2])

        # Compute the root-mean-squared response
        resp_rms = np.sqrt(m0 / gm_duration)

        # Compute the number of zero crossings
        num_zero_crossings = self.limited_num_zero_crossings(
            gm_duration * np.sqrt(m2 / m0) / np.pi)

        # Compute the peak factor
        bar = np.sqrt(2 * np.log(num_zero_crossings))
        peak_factor = bar + 0.5772 / bar

        return super(Davenport1964, self).__call__(
            peak_factor, resp_rms, full_output)


class DerKiureghian1985(Calculator):
    """RVT calculation using peak factor derived by Davenport (1964) [#]_ with
    limits suggested by Der Kiureghian and Igusa [#]_.

    References
    ----------
    .. [#] Davenport, A. G. (1964). Note on the distribution of the largest
        value of a random function with application to gust loading. In
        Institute of Civil Engineering Proceedings, 28(2), 187-196.
    .. [#] Igusa, T., & Der Kiureghian, A. (1985). Dynamic response of
        multiply supported secondary systems. Journal of Engineering Mechanics,
        111(1), 20-41.

    """

    NAME = 'Der Kiureghian (1985)'
    ABBREV = 'DK85'

    def __init__(self, **kwargs):
        super(DerKiureghian1985, self).__init__(**kwargs)

    def __call__(self, gm_duration, freqs, fourier_amps, full_output=False,
                 **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        gm_duration : float
            Duration of the strong-motion phase of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Aris
            intensity [sec]
        freqs : numpy.array
            Frequency of the Fourier amplitude spectrum [Hz]
        fourier_amps: numpy.array
            Amplitude of the Fourier amplitude spectrum with a single
            degree of freedom oscillator already applied if being used. Units
            are not important.
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Damping of the oscillator [decimal]. For example, 0.05 for 5%.
        full_output : bool
            (optional) If the full output should be returned. Default is
            ``False``.

        Returns
        -------
        max_resp : float
            Expected maximum response
        peak_factor : float
            Associated peak factor. Only provided in `full_output` is ``True``.

        """

        m0, m1, m2 = compute_moments(freqs, fourier_amps, [0, 1, 2])

        # Compute the root-mean-squared response
        resp_rms = np.sqrt(m0 / gm_duration)

        # Compute the number of zero crossings
        num_zero_crossings = gm_duration * np.sqrt(m2 / m0) / np.pi

        # Reduce the rate of zero crossings based on the bandwidth
        bandwidth = np.sqrt(1 - (m1 * m1) / (m0 * m2))
        if bandwidth <= 0.1:
            eff_crossings = max(2.1, 2 * bandwidth * num_zero_crossings)
        elif 0.1 < bandwidth <= 0.69:
            eff_crossings = \
                (1.63 * bandwidth ** 0.45 - 0.38) * num_zero_crossings
        else:
            eff_crossings = num_zero_crossings

        eff_crossings = self.limited_num_zero_crossings(eff_crossings)

        # Compute the peak factor
        bar = np.sqrt(2 * np.log(eff_crossings))
        peak_factor = bar + 0.5772 / bar

        return super(DerKiureghian1985, self).__call__(
            peak_factor, resp_rms, full_output)


class ToroMcGuire1987(Calculator):
    """RVT calculation using peak factor derived by Davenport (1964) with
    modifications proposed by Toro and McGuire [#]_.

    References
    ----------
    .. [#] Davenport, A. G. (1964). Note on the distribution of the largest
        value of a random function with application to gust loading. In
        Institute of Civil Engineering Proceedings, 28(2), 187-196.
    .. [#] Toro, G. R., & McGuire, R. K. (1987). An investigation into
        earthquake ground motion characteristics in eastern North America.
        Bulletin of the Seismological Society of America, 77(2), 468-489.

    """

    NAME = 'Toro & McGuire (1987)'
    ABBREV = 'TM87'

    def __init__(self, **kwargs):
        super(ToroMcGuire1987, self).__init__(**kwargs)

    def __call__(self, gm_duration, freqs, fourier_amps, osc_freq=None,
                 osc_damping=None, full_output=False, **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        gm_duration : float
            Duration of the strong-motion phase of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Aris
            intensity [sec]
        freqs : numpy.array
            Frequency of the Fourier amplitude spectrum [Hz]
        fourier_amps: numpy.array
            Amplitude of the Fourier amplitude spectrum with a single
            degree of freedom oscillator already applied if being used. Units
            are not important.
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Damping of the oscillator [decimal]. For example, 0.05 for 5%.
        full_output : bool
            (optional) If the full output should be returned. Default is
            ``False``.

        Returns
        -------
        max_resp : float
            Expected maximum response
        peak_factor : float
            Associated peak factor. Only provided in `full_output` is ``True``.

        """

        m0, m1, m2 = compute_moments(freqs, fourier_amps, [0, 1, 2])

        # Vanmarcke's (1976) bandwidth measure and central frequency
        bandwidth = np.sqrt(1 - (m1 * m1) / (m0 * m2))
        freq_cent = np.sqrt(m2 / m0) / (2 * np.pi)

        num_zero_crossings = self.limited_num_zero_crossings(
            2 * freq_cent * gm_duration * (1.63 * bandwidth ** 0.45 - 0.38))

        foo = np.sqrt(2 * np.log(num_zero_crossings))
        peak_factor = (foo + 0.5772 / foo)

        if osc_freq and osc_damping:
            peak_factor *= Vanmarcke1975.nonstationarity_factor(
                osc_damping, osc_freq, gm_duration)

        # Compute the root-mean-squared response
        resp_rms = np.sqrt(m0 / gm_duration)

        return super(ToroMcGuire1987, self).__call__(
            peak_factor, resp_rms, full_output)


class CartwrightLonguetHiggins1956(Calculator):
    """RVT calculation based on the peak factor definition by Cartwright and
    Longuet-Higgins (1956) [#]_. No modification for stationarity.

    References
    ----------
    .. [#] Cartwright, D. E., & Longuet-Higgins, M. S. (1956). The
        statistical distribution of the maxima of a random function.
        Proceedings of the Royal Society of London. Series A. Mathematical and
        Physical Sciences, 237(1209), 212-232.
    """

    NAME = 'Cartwright & Longuet-Higgins (1956)'
    ABBREV = 'CLH56'

    def __init__(self, **kwargs):
        super(CartwrightLonguetHiggins1956, self).__init__(**kwargs)

    def __call__(self, gm_duration, freqs, fourier_amps, osc_freq=None,
                 osc_damping=None, full_output=False, **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        gm_duration : float
            Duration of the strong-motion phase of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Aris
            intensity [sec]
        freqs : numpy.array
            Frequency of the Fourier amplitude spectrum [Hz]
        fourier_amps: numpy.array
            Amplitude of the Fourier amplitude spectrum with a single
            degree of freedom oscillator already applied if being used. Units
            are not important.
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Fractional damping of the oscillator. For example, 0.05 for 5%
            damping.
        full_output : bool
            (optional) If the full output should be returned. Default is
            ``False``.

        Returns
        -------
        max_resp : float
            Expected maximum response
        peak_factor : float
            Associated peak factor. Only provided in `full_output` is ``True``.

        """
        m0, m1, m2, m4 = compute_moments(freqs, fourier_amps, [0, 1, 2, 4])

        bandwidth = np.sqrt((m2 * m2) / (m0 * m4))
        num_extrema = max(2., np.sqrt(m4 / m2) * gm_duration / np.pi)

        # Compute the peak factor by the indefinite integral.
        peak_factor = np.sqrt(2.) * quad(
            lambda z: 1. - (1. - bandwidth * np.exp(-z * z)) ** num_extrema,
            0, np.inf)[0]

        # Compute the root-mean-squared response -- correcting for the RMS
        # duration.
        if osc_freq and osc_damping:
            rms_duration = self.compute_duration_rms(
                gm_duration, osc_freq, osc_damping, m0, m1, m2)
        else:
            rms_duration = gm_duration

        resp_rms = np.sqrt(m0 / rms_duration)

        return super(CartwrightLonguetHiggins1956, self).__call__(
            peak_factor, resp_rms, full_output)

    def compute_duration_rms(self, gm_duration, osc_freq, osc_damping, m0, m1,
                             m2):
        """Compute the RMS duration if needed."""
        return gm_duration


class BooreJoyner1984(CartwrightLonguetHiggins1956):
    """RVT calculation based on the peak factor definition by Cartwright and
    Longuet-Higgins (1956) [#]_ along with the root-mean-squared duration
    correction proposed by Boore and Joyner (1984) [#]_.

    This RVT calculation is used by SMSIM and is described in Boore (2003)
    [#]_.

    References
    ----------
    .. [#] Cartwright, D. E., & Longuet-Higgins, M. S. (1956). The
        statistical distribution of the maxima of a random function.
        Proceedings of the Royal Society of London. Series A. Mathematical and
        Physical Sciences, 237(1209), 212-232.
    .. [#] Boore, D. M., & Joyner, W. B. (1984). A note on the use of random
        vibration theory to predict peak amplitudes of transient signals.
        Bulletin of the Seismological Society of America, 74(5), 2035-2039.
    .. [#] Boore, D. M. (2003). Simulation of ground motion using the
        stochastic method. Pure and applied geophysics, 160(3-4), 635-676.

    """

    NAME = 'Boore & Joyner (1984)'
    ABBREV = 'BJ84'

    def __init__(self, **kwargs):
        super(BooreJoyner1984, self).__init__(**kwargs)

    def compute_duration_rms(self, gm_duration, osc_freq, osc_damping, m0, m1,
                             m2):
        """Compute the oscillator duration used in the calculation of the
        root-mean-squared response.

        Based on Boore and Joyner (1984).

        Parameters
        ----------
        gm_duration : float
            Duration of the strong-motion phase of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Aris
            intensity [sec]
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Fractional damping of the oscillator. For example, 0.05 for 5%
            damping.

        Returns
        -------
        duration_rms : float
            Duration of the root-mean-squared oscillator response [sec]

        """
        power = 3.
        coef = 1. / 3.

        # This equation was rewritten in Boore and Thompson (2012).
        foo = 1. / (osc_freq * gm_duration)
        dur_ratio = (1 + 1. / (2 * np.pi * osc_damping) *
                     (foo / (1 + coef * foo ** power)))

        return gm_duration * dur_ratio


class LiuPezeshk1999(BooreJoyner1984):
    """ RVT calculation based on the peak factor definition by Cartwright and
    Longuet-Higgins (1956) [#]_ along with the root-mean-squared duration
    correction proposed by Liu and Pezeshk (1999) [#]_.

    References
    ----------
    .. [#] Cartwright, D. E., & Longuet-Higgins, M. S. (1956). The statistical
        distribution of the maxima of a random function. Proceedings of the
        Royal Society of London. Series A. Mathematical and Physical Sciences,
        237(1209), 212-232.
    .. [#] Liu, L., & Pezeshk, S. (1999). An improvement on the estimation of
        pseudoresponse spectral velocity using RVT method. Bulletin of the
        Seismological Society of America, 89(5), 1384-1389.

    """

    NAME = 'Liu & Pezeshk (1999)'
    ABBREV = 'LP99'

    def __init__(self, **kwargs):
        super(LiuPezeshk1999, self).__init__(**kwargs)

    def compute_duration_rms(self, gm_duration, osc_freq, osc_damping, m0, m1,
                             m2):
        """Compute the oscillator duration used in the calculation of the
        root-mean-squared response.

        Based on Liu and Pezeshk (1999).

        Parameters
        ----------
        gm_duration : float
            Duration of the strong-motion phase of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Aris
            intensity [sec]
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Fractional damping of the oscillator. For example, 0.05 for 5%
            damping.
        m0 : float
            Zero-th moment of the Fourier amplitude spectrum
        m1 : float
            First moment of the Fourier amplitude spectrum
        m2 : float
            Second moment of the Fourier amplitude spectrum

        Returns
        -------
        duration_rms : float
            Duration of the root-mean-squared oscillator response [sec]

        """

        power = 2.
        coef = np.sqrt(2 * np.pi * (1. - (m1 * m1) / (m0 * m2)))

        # Same model as used in Boore and Joyner (1984). This equation was
        # rewritten in Boore and Thompson (2012).
        foo = 1. / (osc_freq * gm_duration)
        dur_ratio = (1 + 1. / (2 * np.pi * osc_damping) *
                     (foo / (1 + coef * foo ** power)))

        return gm_duration * dur_ratio


def _load_bt12_data(region):
    """Load data from the Boore and Thompson (2012) parameter files.

    Parameters
    ----------
    region : str
        Region for which the parameters were developed. Possible options are:
            'wna' - Western North America (active tectonic)
            'ena' - Eastern North America (stable tectonic)

    Returns
    -------
    parameters : np.recarray
        Parameters for the region

    """
    fname = os.path.join(
        os.path.dirname(__file__), 'data',
        region + '_bt12_trms4osc.pars')

    return np.rec.fromrecords(
        np.loadtxt(fname, skiprows=4, usecols=range(9)),
        names='mag,dist,c1,c2,c3,c4,c5,c6,c7')

# Load coefficient interpolators for Boore and Thompson (2012)
_BT12_INTERPS = {}
for region in ['wna', 'cena']:
    d = _load_bt12_data(region)
    i = LinearNDInterpolator(
        np.c_[d.mag, np.log(d.dist)],
        np.c_[d.c1, d.c2, d.c3, d.c4, d.c5, d.c6, d.c7])
    _BT12_INTERPS[region] = i


class BooreThompson2012(BooreJoyner1984):
    """RVT calculation based on the peak factor definition by Cartwright and
    Longuet-Higgins (1956) [#]_ along with the root-mean-squared duration
    correction proposed by Boore and Thompson (2012) [#]_.

    References
    ----------
    .. [#] Cartwright, D. E., & Longuet-Higgins, M. S. (1956). The statistical
        distribution of the maxima of a random function.  Proceedings of the
        Royal Society of London. Series A. Mathematical and Physical Sciences,
        237(1209), 212-232.
    .. [#] Boore, D. M., & Thompson, E. M. (2012). Empirical Improvements for
        Estimating Earthquake Response Spectra with Random Vibration Theory.
        Bulletin of the Seismological Society of America, 102(2), 761-772.

    """

    NAME = 'Boore & Thompson (2012)'
    ABBREV = 'BT12'

    def __init__(self, region, mag, dist, **kwargs):
        """Initialize the RVT calculator.

        The duration ratio is defined by Equation (10) in Boore and Thompson
        (2012). Magnitude and distance is interpolated using Qhulls.

        Parameters
        ----------
        region : str
            Region for which the parameters were developed. Possible options
            are:
                'wna' - Western North America (active tectonic)
                'cena' - Central and Eastern North America (stable tectonic)
        mag : float
            Magnitude of the event
        dist : float
            Distance of the event in [km].

        Notes
        -----
        The interpolant is constructed by triangulating the input data
        with Qhull [#r1]_, and on each triangle performing linear
        barycentric interpolation.

        References
        ----------
        .. _[#r1] http://www.qhull.org/

        """
        super(BooreThompson2012, self).__init__(**kwargs)

        region = get_region(region)
        self._COEFS = _BT12_INTERPS[region](mag, np.log(dist))

    def compute_duration_rms(self, gm_duration, osc_freq, osc_damping, *args,
                             **kwargs):
        """Compute the oscillator duration used in the calculation of the
        root-mean-squared response.

        Based on Boore and Joyner (1984).

        Parameters
        ----------
        gm_duration : float
            Duration of the strong-motion phase of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Aris
            intensity [sec]
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Fractional damping of the oscillator. For example, 0.05 for 5%
            damping.

        Returns
        -------
        duration_rms : float
            Duration of the root-mean-squared oscillator response [sec]

        """
        c1, c2, c3, c4, c5, c6, c7 = self._COEFS

        foo = 1 / (osc_freq * gm_duration)
        dur_ratio = ((c1 + c2 * (1 - foo ** c3) / (1 + foo ** c3)) *
                     (1 + c4 / (2 * np.pi * osc_damping) *
                      (foo / (1 + c5 * foo ** c6)) ** c7))

        return gm_duration * dur_ratio


def get_peak_calculator(method, calc_kwds):
    """Select a peak calculator based on a string.

    Parameters
    ----------
    method : str
        Name of the peak calculation method
    calc_kwds : dict
        Keywords passed to the calculator

    Returns
    -------
    calculator : pyrvt.peak_calculators.Calculator

    """

    calculators = [
        BooreJoyner1984,
        BooreThompson2012,
        CartwrightLonguetHiggins1956,
        Davenport1964,
        DerKiureghian1985,
        LiuPezeshk1999,
        ToroMcGuire1987,
        Vanmarcke1975,
    ]

    for calculator in calculators:
        if method in [calculator.NAME, calculator.ABBREV]:
            return calculator(**calc_kwds)

    raise NotImplementedError('No calculator for: %s' % method)


def get_region(region):
    """Return the region naming used in this package.

    Parameters
    ----------
    region : str
        Name for region

    Returns
    -------
    region : str

    """
    region = region.lower()
    if region in ['cena', 'ena', 'ceus', 'eus']:
        return 'cena'
    elif region in ['wna', 'wus']:
        return 'wna'
    else:
        raise NotImplementedError
