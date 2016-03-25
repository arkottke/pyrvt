#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Published peak factor models, which compute the expected peak ground motion. A
specific model may include oscillator duration correction.
"""

import os

import numpy as np

from scipy.integrate import quad
from scipy.interpolate import LinearNDInterpolator


def compute_moments(freqs, fourier_amps, orders):
    """Compute the spectral moments.

    The spectral moment is computed using the squared Fourier amplitude
    spectrum.

    Args:
        freqs (:class:`np.ndarray`): Frequency of the Fourier amplitude
            spectrum (Hz)
        fourier_amps (:class:`np.ndarray`): Amplitude of the Fourier
            amplitude spectrum.
    Returns:
        List[float]: List of computed spectral moments.
    """
    squared_fa = np.square(fourier_amps)

    # Use trapzoidal integration to compute the requested moments.
    moments = [2. * np.trapz(
        np.power(2 * np.pi * freqs, o) * squared_fa, freqs)
        for o in orders]

    return moments


class Calculator(object):
    """:class:`Calculator` is a base class used for all peak calculator
    classes.
    """

    NAME = ''
    ABBREV = ''

    _MIN_ZERO_CROSSINGS = 1.33

    def __init__(self, **kwds):
        pass

    @property
    def name(self):
        """Name of the calculator."""
        return self.NAME

    @property
    def abbrev(self):
        """Abbreviated name of the calculator."""
        return self.ABBREV

    @property
    def min_zero_crossings(self):
        """Minimum number of zero crossings."""
        return self._MIN_ZERO_CROSSINGS

    @classmethod
    def limited_num_zero_crossings(cls, num_zero_crossings):
        """Limit the number of zero crossing to a static limit."""
        return max(cls._MIN_ZERO_CROSSINGS, num_zero_crossings)


class Vanmarcke1975(Calculator):
    """Vanmarcke (1975, :cite:`vanmarcke75`) peak factor which includes the
    effects of clumping.

    The peak factor equation is from Equation (2) in Der Kiureghian (1980,
    :cite:`derkiureghian80`), which is based on Equation (29) in
    :cite:`vanmarcke75`.

    The cumulative density function (CDF) of the peak is defined as:

    .. math::
        F_x(x) = \\left[1 - \\exp\\left(-x^2/2\\right)\\right]
        \\exp\\left[-N_z \\frac{1 -
            \\exp\\left(-\\sqrt{\\pi/2} \\delta_e x\\right)}{\\exp(x^2 / 2) -
            1 }\\right]

    where :math:`N_z` is the number of zero crossings, :math:`\delta_e` is the
    effective bandwidth (:math:`\delta^{1.2}`).

    Typically, the expected value of the peak factor is calculated by
    integrating over the probability density function (i.e., :math:`f_x(x) =
    \\frac{d}{dx} F_x( x)`):

    .. math::
        E[x] = \\int_0^\\infty x f_x(x) dx

    However, because of the properties of :math:`F_x(x)`, specifically that it
    has non-zero probabilities for only positive values, :math:`E[x]` can be
    computed directly from :math:`F_x(x)`.

    .. math::
        E[x] = \\int_0^\\infty 1 - F_x(x) dx.

    This is based on the following sources [#]_ and [#]_.

    .. [#] http://en.wikipedia.org/wiki/Expected_value#Formulas_for_special_cases
    .. [#] http://stats.stackexchange.com/a/13377/48461

    """
    NAME = 'Vanmarcke (1975)'
    ABBREV = 'V75'

    def __init__(self, **kwargs):
        super(Vanmarcke1975, self).__init__(**kwargs)

    def __call__(self, duration, freqs, fourier_amps, osc_freq, osc_damping,
                 **kwargs):
        """Compute the peak response.

        Args:
            duration (float): Duration of the stationary portion of the
                ground motion. Typically defined as the duration between the 5%
                and 75% normalized AriasA intensity (sec).
            freqs (:class:`np.ndarray`): Frequency of the Fourier amplitude
                spectrum (Hz).
            fourier_amps (:class:`np.ndarray`):  Amplitude of the Fourier
                amplitude spectrum with a single degree of freedom oscillator
                already applied if being used. Units are not important.
            osc_freq (float): Frequency of the oscillator (Hz).
            osc_damping (float): Fractional damping of the oscillator (dec). For
                example, 0.05 for a damping ratio of 5%.
        Returns:
            (tuple): tuple containing:

                - max_resp (float): expected maximum response.
                - peak_factor (float): associated peak factor.
        """

        m0, m1, m2 = compute_moments(freqs, fourier_amps, [0, 1, 2])

        # Compute the root-mean-squared response
        resp_rms = np.sqrt(m0 / duration)

        bandwidth = np.sqrt(1 - (m1 * m1) / (m0 * m2))
        bandwidth_eff = bandwidth ** 1.2

        num_zero_crossings = self.limited_num_zero_crossings(
            duration * np.sqrt(m2 / m0) / np.pi)

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
                osc_damping, osc_freq, duration)

        return peak_factor * resp_rms, peak_factor

    @classmethod
    def nonstationarity_factor(cls, osc_damping, osc_freq,  duration):
        return np.sqrt(
            1 - np.exp(-4 * np.pi * osc_damping * osc_freq * duration))


class Davenport1964(Calculator):
    """RVT calculation using the asymptotic solution proposed by
    Davenport (1964, :cite:`davenport64`).
    """

    NAME = 'Davenport (1964)'
    ABBREV = 'D64'

    def __init__(self, **kwargs):
        super(Davenport1964, self).__init__(**kwargs)

    def __call__(self, duration, freqs, fourier_amps, **kwargs):
        """Compute the peak response.

        Args:
            duration (float): Duration of the stationary portion of the
                ground motion. Typically defined as the duration between the 5%
                and 75% normalized AriasA intensity (sec).
            freqs (:class:`np.ndarray`): Frequency of the Fourier amplitude
                spectrum (Hz).
            fourier_amps (:class:`np.ndarray`):  Amplitude of the Fourier
                amplitude spectrum with a single degree of freedom oscillator
                already applied if being used. Units are not important.
        Returns:
            tuple: tuple containing:

                - max_resp (float): expected maximum response.
                - peak_factor (float): associated peak factor.
        """
        m0, m2 = compute_moments(freqs, fourier_amps, [0, 2])

        # Compute the root-mean-squared response
        resp_rms = np.sqrt(m0 / duration)

        # Compute the number of zero crossings
        num_zero_crossings = self.limited_num_zero_crossings(
            duration * np.sqrt(m2 / m0) / np.pi)

        peak_factor = self.asymtotic_approx(num_zero_crossings)

        return peak_factor * resp_rms, peak_factor

    @classmethod
    def asymtotic_approx(self, zero_crossings):
        """Compute the peak factor from the asymptotic approximation.

        Args:
            zero_crossings (float): Number of zero crossing.
        Returns:
            (float): calculated peak factor
        """
        x = np.sqrt(2 * np.log(zero_crossings))
        return x + 0.5772 / x


class DerKiureghian1985(Davenport1964):
    """RVT calculation using peak factor derived by Davenport (1964,
    :cite:`davenport64`) with limits suggested by Igusa & Der Kiureghian
    (1985, :cite:`igusa85`).
    """

    NAME = 'Der Kiureghian (1985)'
    ABBREV = 'DK85'

    def __init__(self, **kwargs):
        super(DerKiureghian1985, self).__init__(**kwargs)

    def __call__(self, duration, freqs, fourier_amps, **kwargs):
        """Compute the peak response.

        Args:
            duration (float): Duration of the stationary portion of the
                ground motion. Typically defined as the duration between the 5%
                and 75% normalized AriasA intensity (sec).
            freqs (:class:`np.ndarray`): Frequency of the Fourier amplitude
                spectrum (Hz).
            fourier_amps (:class:`np.ndarray`):  Amplitude of the Fourier
                amplitude spectrum with a single degree of freedom oscillator
                already applied if being used. Units are not important.
        Returns:
            tuple: tuple containing:

                - max_resp (float): expected maximum response.
                - peak_factor (float): associated peak factor.
        """

        m0, m1, m2 = compute_moments(freqs, fourier_amps, [0, 1, 2])

        # Compute the root-mean-squared response
        resp_rms = np.sqrt(m0 / duration)

        # Compute the number of zero crossings
        num_zero_crossings = duration * np.sqrt(m2 / m0) / np.pi

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
        peak_factor = self.asymtotic_approx(eff_crossings)

        return peak_factor * resp_rms, peak_factor


class ToroMcGuire1987(Davenport1964):
    """RVT calculation using asymptotic solution proposed by Davenport
    (1964, :cite:`davenport64`) with modifications proposed by Toro & McGuire
    (1987, :cite:`toro87`).
    """

    NAME = 'Toro & McGuire (1987)'
    ABBREV = 'TM87'

    def __init__(self, **kwargs):
        super(ToroMcGuire1987, self).__init__(**kwargs)

    def __call__(self, duration, freqs, fourier_amps, osc_freq=None,
                 osc_damping=None, **kwargs):
        """Compute the peak response.

        Args:
            duration (float): Duration of the stationary portion of the
                ground motion. Typically defined as the duration between the 5%
                and 75% normalized AriasA intensity (sec).
            freqs (:class:`np.ndarray`): Frequency of the Fourier amplitude
                spectrum (Hz).
            fourier_amps (:class:`np.ndarray`):  Amplitude of the Fourier
                amplitude spectrum with a single degree of freedom oscillator
                already applied if being used. Units are not important.
            osc_freq (float): Frequency of the oscillator (Hz).
            osc_damping (float): Fractional damping of the oscillator (dec). For
                example, 0.05 for a damping ratio of 5%.
        Returns:
            tuple: tuple containing:

                - max_resp (float): expected maximum response.
                - peak_factor (float): associated peak factor.
        """

        m0, m1, m2 = compute_moments(freqs, fourier_amps, [0, 1, 2])

        # Vanmarcke's (1976) bandwidth measure and central frequency
        bandwidth = np.sqrt(1 - (m1 * m1) / (m0 * m2))
        freq_cent = np.sqrt(m2 / m0) / (2 * np.pi)

        num_zero_crossings = self.limited_num_zero_crossings(
            2 * freq_cent * duration * (1.63 * bandwidth ** 0.45 - 0.38))

        peak_factor = self.asymtotic_approx(num_zero_crossings)

        if osc_freq and osc_damping:
            peak_factor *= Vanmarcke1975.nonstationarity_factor(
                osc_damping, osc_freq, duration)

        # Compute the root-mean-squared response
        resp_rms = np.sqrt(m0 / duration)

        return peak_factor * resp_rms, peak_factor


class CartwrightLonguetHiggins1956(Calculator):
    """RVT calculation based on the peak factor definition by Cartwright and
    Longuet-Higgins (1956, :cite:`cartwright56`) using the
    integral provided by Boore (2003, :cite:`boore03`).
    """

    NAME = 'Cartwright & Longuet-Higgins (1956)'
    ABBREV = 'CLH56'

    def __init__(self, **kwargs):
        super(CartwrightLonguetHiggins1956, self).__init__(**kwargs)

    def __call__(self, duration, freqs, fourier_amps, osc_freq=None,
                 osc_damping=None, **kwargs):
        """Compute the peak response.

        Args:
            duration (float): Duration of the stationary portion of the
                ground motion. Typically defined as the duration between the 5%
                and 75% normalized AriasA intensity (sec).
            freqs (:class:`np.ndarray`): Frequency of the Fourier amplitude
                spectrum (Hz).
            fourier_amps (:class:`np.ndarray`):  Amplitude of the Fourier
                amplitude spectrum with a single degree of freedom oscillator
                already applied if being used. Units are not important.
            osc_freq (float): Frequency of the oscillator (Hz).
            osc_damping (float): Fractional damping of the oscillator (dec). For
                example, 0.05 for a damping ratio of 5%.
        Returns:
            tuple: tuple containing:

                - max_resp (float): expected maximum response.
                - peak_factor (float): associated peak factor.
        """
        m0, m1, m2, m4 = compute_moments(freqs, fourier_amps, [0, 1, 2, 4])

        bandwidth = np.sqrt((m2 * m2) / (m0 * m4))
        num_extrema = max(2., np.sqrt(m4 / m2) * duration / np.pi)

        # Compute the peak factor by the indefinite integral.
        peak_factor = np.sqrt(2.) * quad(
            lambda z: 1. - (1. - bandwidth * np.exp(-z * z)) ** num_extrema,
            0, np.inf)[0]

        # Compute the root-mean-squared response -- correcting for the RMS
        # duration.
        if osc_freq and osc_damping:
            rms_duration = self.compute_duration_rms(
                duration, osc_freq, osc_damping, m0, m1, m2)
        else:
            rms_duration = duration

        resp_rms = np.sqrt(m0 / rms_duration)

        return peak_factor * resp_rms, peak_factor

    def compute_duration_rms(self, duration, osc_freq, osc_damping, m0, m1, m2):
        """Compute the RMS duration.

        Not used by :class:`.CartwrightLonguetHiggins1956`.

        Args:
            duration (float): Duration of the stationary portion of the
                ground motion. Typically defined as the duration between the 5%
                and 75% normalized AriasA intensity (sec).
            osc_freq (float): Frequency of the oscillator (Hz).
            osc_damping (float): Fractional damping of the oscillator (dec). For
                example, 0.05 for a damping ratio of 5%.
            m0 (float): Zero-th moment of the Fourier amplitude spectrum.
            m1 (float): First moment of the Fourier amplitude spectrum.
            m2 (float): Second moment of the Fourier amplitude spectrum
        Returns:
            (float): Duration of the root-mean-squared oscillator response
                (sec).
        """
        del (osc_freq, osc_damping, m0, m1, m2)
        return duration


class BooreJoyner1984(CartwrightLonguetHiggins1956):
    """RVT calculation based on the peak factor definition by Cartwright &
    Longuet-Higgins (1956, :cite:`cartwright56`) and along with the
    root-mean-squared duration correction proposed by Boore & Joyner (1984,
    :cite:`boore84`).

    This RVT calculation is used by SMSIM and is described in Boore
    (2003, :cite:`boore03`).
    """

    NAME = 'Boore & Joyner (1984)'
    ABBREV = 'BJ84'

    def __init__(self, **kwargs):
        super(BooreJoyner1984, self).__init__(**kwargs)

    def compute_duration_rms(self, duration, osc_freq, osc_damping, m0, m1, m2):
        """Compute the oscillator duration used in the calculation of the
        root-mean-squared response.

        Based on  :cite:`boore84`.

        Args:
            duration (float): Duration of the stationary portion of the
                ground motion. Typically defined as the duration between the 5%
                and 75% normalized AriasA intensity (sec).
            osc_freq (float): Frequency of the oscillator (Hz).
            osc_damping (float): Fractional damping of the oscillator (dec). For
                example, 0.05 for a damping ratio of 5%.
            m0 (float): Zero-th moment of the Fourier amplitude spectrum.
            m1 (float): First moment of the Fourier amplitude spectrum.
            m2 (float): Second moment of the Fourier amplitude spectrum
        Returns:
            (float): Duration of the root-mean-squared oscillator response
                (sec).
        """
        del (m0, m1, m2)

        power = 3.
        coef = 1. / 3.

        # This equation was rewritten in Boore and Thompson (2012).
        foo = 1. / (osc_freq * duration)
        dur_ratio = (1 + 1. / (2 * np.pi * osc_damping) *
                     (foo / (1 + coef * foo ** power)))

        return duration * dur_ratio


class LiuPezeshk1999(BooreJoyner1984):
    """RVT calculation based on the peak factor definition by Cartwright &
    Longuet-Higgins (1956, :cite:`cartwright56`) along with the
    root-mean-squared duration correction proposed by Liu & Pezeshk
    (1999, :cite:`liu99`).
    """

    NAME = 'Liu & Pezeshk (1999)'
    ABBREV = 'LP99'

    def __init__(self, **kwargs):
        super(LiuPezeshk1999, self).__init__(**kwargs)

    def compute_duration_rms(self, duration, osc_freq, osc_damping, m0, m1, m2):
        """Compute the oscillator duration used in the calculation of the
        root-mean-squared response.

        Based on :cite:`liu99`.

        Args:
            duration (float): Duration of the stationary portion of the
                ground motion. Typically defined as the duration between the 5%
                and 75% normalized AriasA intensity (sec).
            osc_freq (float): Frequency of the oscillator (Hz).
            osc_damping (float): Fractional damping of the oscillator (dec). For
                example, 0.05 for a damping ratio of 5%.
            m0 (float): Zero-th moment of the Fourier amplitude spectrum.
            m1 (float): First moment of the Fourier amplitude spectrum.
            m2 (float): Second moment of the Fourier amplitude spectrum
        Returns:
            (float): Duration of the root-mean-squared oscillator response
                (sec).
        """
        power = 2.
        coef = np.sqrt(2 * np.pi * (1. - (m1 * m1) / (m0 * m2)))

        # Same model as used in Boore and Joyner (1984). This equation was
        # rewritten in Boore and Thompson (2012).
        foo = 1. / (osc_freq * duration)
        dur_ratio = (1 + 1. / (2 * np.pi * osc_damping) *
                     (foo / (1 + coef * foo ** power)))

        return duration * dur_ratio


def _load_bt12_data(region):
    """Load data from the Boore & Thompson (2012, :cite:`boore12`) parameter
    files.

    Args:
        region (str): Region for which the parameters were developed. Valid
            options: 'wna' for Western North America (active tectonic), or
            'cena' for Eastern North America (stable tectonic).
    Returns:
        (np.recarray): Parameters for the region.
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
    """:class:`BooreThompson2012` is an RVT calculator based on the peak factor
    definition by Cartwright & Longuet-Higgins (1956, :cite:`cartwright56`
    along with the root-mean-squared duration correction proposed by Boore &
    Thompson (2012, :cite:`boore12`).

    The duration ratio is defined by Equation (10) in :cite:`boore12`.
    Magnitude and distance is interpolated using Qhulls.

    Notes:
        The interpolant is constructed by triangulating the input data with
        Qhull [#]_, and on each triangle performing linear barycentric
        interpolation.

    .. [#] http://www.qhull.org/
    """
    NAME = 'Boore & Thompson (2012)'
    ABBREV = 'BT12'

    def __init__(self, region, mag, dist, **kwargs):
        """Initialize the model and interpolate the coefficients.

        Args:
            region (str): Region for which the parameters were developed.
                Valid options are: 'wna' for Western North America (active
                tectonic), and 'cena' for Central and Eastern North America (
                stable tectonic).
            mag (float): Magnitude of the event.
            dist (float): Distance of the event in (km).
        """
        super(BooreThompson2012, self).__init__(**kwargs)

        region = get_region(region)
        self._COEFS = _BT12_INTERPS[region](mag, np.log(dist))

    def compute_duration_rms(self, duration, osc_freq, osc_damping, m0, m1, m2):
        """Compute the RMS duration.

        Based on :cite:`boore12`.

        Args:
            duration (float): Duration of the stationary portion of the
                ground motion. Typically defined as the duration between the 5%
                and 75% normalized AriasA intensity (sec).
            osc_freq (float): Frequency of the oscillator (Hz).
            osc_damping (float): Fractional damping of the oscillator (dec). For
                example, 0.05 for a damping ratio of 5%.
            m0 (float): Zero-th moment of the Fourier amplitude spectrum.
            m1 (float): First moment of the Fourier amplitude spectrum.
            m2 (float): Second moment of the Fourier amplitude spectrum
        Returns:
            (float): Duration of the root-mean-squared oscillator response
                (sec).
        """
        del (m0, m1, m2)

        c1, c2, c3, c4, c5, c6, c7 = self._COEFS

        foo = 1 / (osc_freq * duration)
        dur_ratio = ((c1 + c2 * (1 - foo ** c3) / (1 + foo ** c3)) *
                     (1 + c4 / (2 * np.pi * osc_damping) *
                      (foo / (1 + c5 * foo ** c6)) ** c7))

        return duration * dur_ratio


def get_peak_calculator(method, calc_kwds):
    """Select a peak calculator based on a string.

    Args:
        method (str): Name of the peak calculation method
        calc_kwds (Dict): Keywords passed to the calculator
    Returns:
        (calculator): :class:`.Calculator`
    """
    calc_kwds = calc_kwds or dict()

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

    raise NotImplementedError('No calculator for: %s', method)


def get_region(region):
    """Return the region naming used in this package.

    Args:
        region (str): Regional synonym.
    Returns:
        (str): Region either 'cena' or 'wna'.
    """
    region = region.lower()
    if region in ['cena', 'ena', 'ceus', 'eus']:
        return 'cena'
    elif region in ['wna', 'wus']:
        return 'wna'
    else:
        raise NotImplementedError
