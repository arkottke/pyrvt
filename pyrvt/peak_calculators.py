#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Peak factor models.

Published peak factor models, which compute the expected peak ground motion. A
specific model may include oscillator duration correction.
"""

import ctypes
import itertools
import pathlib

import numpy as np
import numba

from scipy.integrate import quad
from scipy.interpolate import LinearNDInterpolator
from scipy.signal import argrelmax


@numba.jit
def trapz(x, y):
    """Trapezoidal integration written in numba.

    Parameters
    ----------
    x : array_like
        sample points to corresponding to the `y` values.
    y : array_like
        Input array to integrate

    Returns
    -------
    total : float
        Definite integral as approximated by the trapezoidal rule.

    """
    n = x.shape[0]
    total = 0
    for i in range(n - 1):
        total += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])

    return total


# From: https://github.com/scipy/scipy/issues/4831#issuecomment-258501648
# Set the signature of the future types function
# The second argument must be a pointer or it won't work right,
# but this will cause problems because of the bug in scipy
# because scipy looks for a double instead of a pointer
c_sig = numba.types.double(numba.types.intc,
                           numba.types.CPointer(numba.types.double))


# Turn the integrand function into a ctypes function
@numba.cfunc(c_sig)
def _calc_vanmarcke1975_ccdf(n, a):
    """Calculate the Vanmarcke (1975) complementary CDF.

    Parameters
    ----------
    n : int
        Length of arguments
    a : list of floats
        Arguments specifying:
            - function value `x`
            - number of zero crossing
            - effective bandwdith

    Returns
    -------
    ccdf : float
        Complementary CDF value

    """
    args = numba.carray(a, n)
    x = args[0]
    num_zero_crossings = args[1]
    bandwidth_eff = args[2]

    return (1 - (1 - np.exp(-x ** 2 / 2)) * np.exp(-1 * num_zero_crossings * (
        1 - np.exp(-1 * np.sqrt(np.pi / 2) * bandwidth_eff * x)) /
                                                   (np.exp(x ** 2 / 2) - 1)))


# Force the argtypes to be what quad expects
_calc_vanmarcke1975_ccdf.ctypes.argtypes = (ctypes.c_int, ctypes.c_double)


@numba.cfunc(c_sig)
def _calc_cartwright_pf(n, a):
    """Integrand for the Cartwright and Longuet-Higgins peak factor.

    Parameters
    ----------
    n : int
        Length of arguments
    a : list of floats
        Arguments specifying:
            - function value `x`
            - number of extrema
            - bandwdith

    Returns
    -------
    dpf : float
        Portion of the peak factor

    """
    args = numba.carray(a, n)
    x = args[0]
    num_extrema = args[1]
    bandwidth = args[2]
    return 1. - (1. - bandwidth * np.exp(-x * x)) ** num_extrema


# Force the argtypes to be what quad expects
_calc_cartwright_pf.ctypes.argtypes = (ctypes.c_int, ctypes.c_double)


def calc_moments(freqs, fourier_amps, orders):

    squared_fa = np.square(fourier_amps)

    # Use trapzoidal integration to compute the requested moments.
    moments = [
        2. * trapz(freqs, np.power(2 * np.pi * freqs, o) * squared_fa)
        for o in orders
    ]

    return moments


class SquaredSpectrum(object):
    """Squared Fourier amplitude spectrum.

    Used to store calculated spectral moments during calculations.

    Parameters
    ----------
    freqs : array_like
        Frequency of the Fourier amplitude spectrum (Hz)
    fourier_amps : array_like
        Amplitude of the Fourier amplitude spectrum.
    """

    def __init__(self, freqs, fourier_amps):
        self._freqs = freqs
        self._squared_fa = np.square(fourier_amps)
        self._moments = {}

    def moment(self, num):
        """Compute the spectral moments.

        The spectral moment is computed using the squared Fourier amplitude
        spectrum.

        Returns
        -------
        moment : float
            Computed spectral moments.
        """
        num = int(num)
        try:
            moment = self._moments[num]
        except KeyError:
            moment = 2. * trapz(self._freqs,
                                np.power(2 * np.pi * self._freqs, num) *
                                self._squared_fa)
            self._moments[num] = moment

        return moment

    def moments(self, *nums):
        return [self.moment(n) for n in nums]


class Calculator(object):
    """Base class used for all peak calculator classes."""

    NAME = ''
    ABBREV = ''

    _MIN_ZERO_CROSSINGS = 1.33

    def __init__(self, **kwds):
        """Initialize the object."""
        super().__init__()
        self._spectrum = None

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

    def __call__(self, duration, freqs, fourier_amps, **kwargs):
        """Compute the peak response.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).
        freqs : array_like
            Frequency of the Fourier amplitude spectrum (Hz).
        fourier_amps : array_like
             Amplitude of the Fourier amplitude spectrum with a single degree
             of freedom oscillator already applied if being used. Units are
             not important.

        Returns
        -------
        max_resp : float
            expected maximum response.
        peak_factor : float
            associated peak factor.

        """
        self._spectrum = SquaredSpectrum(freqs, fourier_amps)

        peak_factor = self._calc_peak_factor(duration, **kwargs)

        duration_rms = self._calc_duration_rms(duration, freqs=freqs, **kwargs)
        # Compute the root-mean-squared response.
        resp_rms = np.sqrt(self._spectrum.moment(0) / duration_rms)

        self._spectrum = None
        return peak_factor * resp_rms, peak_factor

    def _calc_peak_factor(self, duration, **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).
        freqs : array_like
            Frequency of the Fourier amplitude spectrum (Hz).
        fourier_amps : array_like
             Amplitude of the Fourier amplitude spectrum with a single degree
             of freedom oscillator already applied if being used. Units are
             not important.

        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        raise NotImplementedError

    def _calc_duration_rms(self, duration, **kwargs):
        """Modify a duration to correct for stationarity.

        Default implemenation does nothing.

        Returns
        -------
        duration : float
            Modified duration.
        """
        return duration


class Vanmarcke1975(Calculator):
    r"""Vanmarcke (1975) peak factor.

    The Vanmarcke (1975, :cite:`vanmarcke75`) peak factor, which includes the
    effects of clumping.  The peak factor equation is from Equation (2) in Der
    Kiureghian (1980, :cite:`derkiureghian80`), which is based on Equation (29)
    in :cite:`vanmarcke75`.

    The cumulative density function (CDF) of the peak is defined as:

    .. math::
        F_x(x) = \left[1 - \exp\left(-x^2/2\right)\right]
        \exp\left[-N_z \frac{1 -
            \exp\left(-\sqrt{\pi/2} \delta_e x\right)}{\exp(x^2 / 2) -
            1 }\right]

    where :math:`N_z` is the number of zero crossings, :math:`\delta_e` is the
    effective bandwidth (:math:`\delta^{1.2}`).

    Typically, the expected value of the peak factor is calculated by
    integrating over the probability density function (i.e., :math:`f_x(x) =
    \frac{d}{dx} F_x( x)`):

    .. math::
        E[x] = \int_0^\infty x f_x(x) dx

    However, because of the properties of :math:`F_x(x)`, specifically that it
    has non-zero probabilities for only positive values, :math:`E[x]` can be
    computed directly from :math:`F_x(x)`.

    .. math::
        E[x] = \int_0^\infty 1 - F_x(x) dx.

    This is based on the following sources [#]_ and [#]_.

    .. # noqa
    .. [#] http://en.wikipedia.org/wiki/Expected_value#Formulas_for_special_cases
    .. [#] http://stats.stackexchange.com/a/13377/48461

    Parameters
    ----------
    use_nonstationarity_factor : bool
        If the non-stationarity factor should be applied.
    """

    NAME = 'Vanmarcke (1975)'
    ABBREV = 'V75'

    def __init__(self, use_nonstationarity_factor=True, **kwargs):
        """Initialize the class."""
        super().__init__(**kwargs)
        self._use_nonstationarity_factor = use_nonstationarity_factor

    def _calc_peak_factor(self, duration, **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).
        osc_freq : float
            Frequency of the oscillator (Hz).
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for
            a damping ratio of 5%.
        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        m0, m1, m2 = self._spectrum.moments(0, 1, 2)

        bandwidth = np.sqrt(1 - (m1 * m1) / (m0 * m2))
        bandwidth_eff = bandwidth ** 1.2

        num_zero_crossings = self.limited_num_zero_crossings(
            duration * np.sqrt(m2 / m0) / np.pi)
        # The expected peak factor is computed as the integral of the
        # complementary CDF (1 - CDF(x)).
        peak_factor = quad(
            _calc_vanmarcke1975_ccdf.ctypes,
            0,
            np.inf,
            args=(num_zero_crossings, bandwidth_eff))[0]

        osc_freq = kwargs.get('osc_freq', None)
        osc_damping = kwargs.get('osc_damping', None)
        if (osc_freq and osc_damping) and self._use_nonstationarity_factor:
            peak_factor *= self.nonstationarity_factor(osc_damping, osc_freq,
                                                       duration)

        return peak_factor

    @classmethod
    def nonstationarity_factor(cls, osc_damping, osc_freq, duration):
        """Compute nonstationarity factor to modify duration.

        Parameters
        ----------
        osc_damping : float
            Oscillator damping (decimal).
        osc_freq : float
            Oscillator frequency (Hz).
        duration : float
            Duration of the stationary portion of the ground motion

        Returns
        -------
        float
            Nonstationarity factor.

        """
        return np.sqrt(1 - np.exp(-4 * np.pi * osc_damping * osc_freq *
                                  duration))


class Davenport1964(Calculator):
    """Davenport (1964) peak factor.

    RVT calculation using the asymptotic solution proposed by
    Davenport (1964, :cite:`davenport64`).
    """

    NAME = 'Davenport (1964)'
    ABBREV = 'D64'

    def __init__(self, **kwargs):
        """Initialize the class."""
        super().__init__(**kwargs)

    def _calc_peak_factor(self, duration, **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).

        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        m0, m2 = self._spectrum.moments(0, 2)

        # Compute the number of zero crossings
        num_zero_crossings = self.limited_num_zero_crossings(
            duration * np.sqrt(m2 / m0) / np.pi)

        peak_factor = self.asymtotic_approx(num_zero_crossings)

        return peak_factor

    @classmethod
    def asymtotic_approx(self, zero_crossings):
        """Compute the peak factor from the asymptotic approximation.

        Parameters
        ----------
        zero_crossings : float
            Number of zero crossing.

        Returns
        -------
        approx_peak_factor : float
            Calculated peak factor.

        """
        x = np.sqrt(2 * np.log(zero_crossings))
        return x + 0.5772 / x


class DerKiureghian1985(Davenport1964):
    """Der Kiureghian (1985) peak factor.

    RVT calculation using peak factor derived by Davenport (1964,
    :cite:`davenport64`) with limits suggested by Igusa & Der Kiureghian
    (1985, :cite:`igusa85`).
    """

    NAME = 'Der Kiureghian (1985)'
    ABBREV = 'DK85'

    def __init__(self, **kwargs):
        """Initialize the class."""
        super().__init__(**kwargs)

    def _calc_peak_factor(self, duration, **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).
        freqs : array_like
            Frequency of the Fourier amplitude spectrum (Hz).
        fourier_amps : array_like
             Amplitude of the Fourier amplitude spectrum with a single degree
             of freedom oscillator already applied if being used. Units are
             not important.

        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        m0, m1, m2 = self._spectrum.moments(0, 1, 2)

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

        return peak_factor


class ToroMcGuire1987(Davenport1964):
    """Toro and McGuire (1987) peak factor.

    Peak factor equation using asymptotic solution proposed by Davenport (1964,
    :cite:`davenport64`) with modifications proposed by Toro & McGuire (1987,
    :cite:`toro87`).
    """

    NAME = 'Toro & McGuire (1987)'
    ABBREV = 'TM87'

    def __init__(self, **kwargs):
        """Initialize the class."""
        super().__init__(**kwargs)

    def _calc_peak_factor(self, duration, **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).
        freqs : array_like
            Frequency of the Fourier amplitude spectrum (Hz).
        fourier_amps : array_like
             Amplitude of the Fourier amplitude spectrum with a single degree
             of freedom oscillator already applied if being used. Units are
             not important.

        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        m0, m1, m2 = self._spectrum.moments(0, 1, 2)

        # Vanmarcke's (1976) bandwidth measure and central frequency
        bandwidth = np.sqrt(1 - (m1 * m1) / (m0 * m2))
        freq_cent = np.sqrt(m2 / m0) / (2 * np.pi)

        num_zero_crossings = self.limited_num_zero_crossings(
            2 * freq_cent * duration * (1.63 * bandwidth ** 0.45 - 0.38))

        peak_factor = self.asymtotic_approx(num_zero_crossings)

        osc_freq = kwargs.get('osc_freq', None)
        osc_damping = kwargs.get('osc_damping', None)
        if osc_freq and osc_damping:
            peak_factor *= Vanmarcke1975.nonstationarity_factor(
                osc_damping, osc_freq, duration)

        return peak_factor


class CartwrightLonguetHiggins1956(Calculator):
    """Cartwight and Longuet-Higgins (1956) peak factor.

    RVT calculation based on the peak factor definition by Cartwright and
    Longuet-Higgins (1956, :cite:`cartwright56`) using the
    integral provided by Boore (2003, :cite:`boore03`).
    """

    NAME = 'Cartwright & Longuet-Higgins (1956)'
    ABBREV = 'CLH56'

    def __init__(self, **kwargs):
        """Initialize the class."""
        super().__init__(**kwargs)

    def _calc_peak_factor(self, duration, **kwargs):
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).
        freqs : array_like
            Frequency of the Fourier amplitude spectrum (Hz).
        fourier_amps : array_like
             Amplitude of the Fourier amplitude spectrum with a single degree
             of freedom oscillator already applied if being used. Units are
             not important.

        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        m0, m1, m2, m4 = self._spectrum.moments(0, 1, 2, 4)

        bandwidth = np.sqrt((m2 * m2) / (m0 * m4))
        num_extrema = max(2., np.sqrt(m4 / m2) * duration / np.pi)
        # Compute the peak factor by the indefinite integral.
        peak_factor = np.sqrt(2.) * quad(
            _calc_cartwright_pf.ctypes,
            0,
            np.inf,
            args=(num_extrema, bandwidth))[0]

        return peak_factor


class BooreJoyner1984(CartwrightLonguetHiggins1956):
    """Boore and Joyner (1984) peak factor.

    RVT calculation based on the peak factor definition by Cartwright &
    Longuet-Higgins (1956, :cite:`cartwright56`) and along with the
    root-mean-squared duration correction proposed by Boore & Joyner (1984,
    :cite:`boore84`).

    This RVT calculation is used by SMSIM and is described in Boore
    (2003, :cite:`boore03`).
    """

    NAME = 'Boore & Joyner (1984)'
    ABBREV = 'BJ84'

    def __init__(self, **kwargs):
        """Initialize the class."""
        super().__init__(**kwargs)

    def _calc_duration_rms(self, duration, **kwargs):
        """Compute the oscillator duration.

        Oscillator duration is used in the calculation of the root-mean-squared
        response and based on  :cite:`boore84`.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).
        osc_freq : float
            Frequency of the oscillator (Hz).
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for
            a damping ratio of 5%.

        Returns
        -------
        duration_rms : float
            Duration of the root-mean-squared oscillator response (sec).

        """
        osc_freq = kwargs.get('osc_freq', None)
        osc_damping = kwargs.get('osc_damping', None)

        if osc_damping and osc_freq:
            power = 3.
            coef = 1. / 3.
            # This equation was rewritten in Boore and Thompson (2012).
            foo = 1. / (osc_freq * duration)
            dur_ratio = (1 + 1. / (2 * np.pi * osc_damping) *
                         (foo / (1 + coef * foo ** power)))
            duration *= dur_ratio

        return duration


class LiuPezeshk1999(BooreJoyner1984):
    """Liu and Pezeshk (1999) peak factor.

    RVT calculation based on the peak factor definition by Cartwright &
    Longuet-Higgins (1956, :cite:`cartwright56`) along with the
    root-mean-squared duration correction proposed by Liu & Pezeshk
    (1999, :cite:`liu99`).
    """

    NAME = 'Liu & Pezeshk (1999)'
    ABBREV = 'LP99'

    def __init__(self, **kwargs):
        """Initialize the class."""
        super().__init__(**kwargs)

    def _calc_duration_rms(self, duration, **kwargs):
        """Compute the oscillator duration.

        Oscillator duration is used in the calculation of the root-mean-squared
        response and based on :cite:`liu99`.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).
        osc_freq : float
            Frequency of the oscillator (Hz).
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for a
            damping ratio of 5%.

        Returns
        -------
        duration_rms : float
            Duration of the root-mean-squared oscillator response (sec).

        """
        osc_freq = kwargs.get('osc_freq', None)
        osc_damping = kwargs.get('osc_damping', None)
        if osc_freq and osc_damping:
            m0, m1, m2 = self._spectrum.moments(0, 1, 2)

            power = 2.
            coef = np.sqrt(2 * np.pi * (1. - (m1 * m1) / (m0 * m2)))

            # Same model as used in Boore and Joyner (1984). This equation was
            # rewritten in Boore and Thompson (2012).
            foo = 1. / (osc_freq * duration)
            dur_ratio = (1 + 1. / (2 * np.pi * osc_damping) *
                         (foo / (1 + coef * foo ** power)))
            duration *= dur_ratio

        return duration


def _make_bt_interpolator(region, ref):
    """Load data from the Boore & Thompson (2012) parameter files.

    Parameters
    ----------
    region : str
        Region for which the parameters were developed. Valid options: 'wna'
        for Western North America (active tectonic), or 'cena' for Eastern
        North America (stable tectonic).
    ref : str
        Reference document. Either: bt12 or bt15 for Boore & Thompson (2012) or
         (2015), respectively.

    Returns
    -------
    interpolator : :class:`scipy.interpolate.LinearNDInterpolator`
        Interpolator for the data.

    """

    fpath = pathlib.Path(__file__).parent.joinpath(
        'data', f'{region}_{ref}_trms4osc.pars.gz')
    d = np.rec.fromrecords(
        np.loadtxt(str(fpath), skiprows=4, usecols=range(9)),
        names='mag,dist,c1,c2,c3,c4,c5,c6,c7')

    return LinearNDInterpolator(
        np.c_[d.mag, np.log(d.dist)],
        np.c_[d.c1, d.c2, d.c3, d.c4, d.c5, d.c6, d.c7])


# Load coefficient interpolators for Boore and Thompson (2012)
_BT_INTERPS = {
    (region, ref): _make_bt_interpolator(region, ref)
    for region, ref in itertools.product(['wna', 'cena'], ['bt12', 'bt15'])
}


class BooreThompson(object):
    """Abstract class for the Boore & Thompson duration correction.

    The duration ratio is defined by Equation (10) in :cite:`boore12`.
    Magnitude and distance is interpolated using
    `scipy.interpolate.LinearNDInterpolator` on the natural log of the
    distance.

    Parameters
    ----------
    region : str
        Region for which the parameters were developed.  Valid options
        are: 'wna' for Western North America (active tectonic), and 'cena'
        for Central and Eastern North America ( stable tectonic).
    mag : float
        Magnitude of the event.
    dist : float
        Distance of the event in (km).
    ref : str
        Reference for coefficients, either: 'bt12' or 'bt15'
    """

    def __init__(self, region, mag, dist, ref, **kwargs):
        """Initialize the class."""
        super().__init__(**kwargs)
        region = get_region(region)
        self._COEFS = _BT_INTERPS[(region, ref)](mag, np.log(dist))

    def _calc_duration_rms(self, duration, **kwargs):
        """Compute the RMS duration.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).
        osc_freq : float
            Frequency of the oscillator (Hz).
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for
            a damping ratio of 5%.

        Returns
        -------
        duration_rms : float
            Duration of the root-mean-squared oscillator response (sec).

        """
        osc_freq = kwargs.get('osc_freq', None)
        osc_damping = kwargs.get('osc_damping', None)
        if osc_freq and osc_damping:
            c1, c2, c3, c4, c5, c6, c7 = self._COEFS

            foo = 1 / (osc_freq * duration)
            dur_ratio = ((c1 + c2 * (1 - foo ** c3) / (1 + foo ** c3)) *
                         (1 + c4 / (2 * np.pi * osc_damping) *
                          (foo / (1 + c5 * foo ** c6)) ** c7))
            duration *= dur_ratio

        return duration


class BooreThompson2012(BooreThompson, BooreJoyner1984):
    """Boore and Thompson (2012) peak factor.

    Peak calculation based on the peak factor definition by Cartwright &
    Longuet-Higgins (1956, :cite:`cartwright56` along with the
    root-mean-squared duration correction proposed by Boore & Thompson (2012,
    :cite:`boore12`).


    Parameters
    ----------
    region : str
        Region for which the parameters were developed.  Valid options
        are: 'wna' for Western North America (active tectonic), and 'cena'
        for Central and Eastern North America ( stable tectonic).
    mag : float
        Magnitude of the event.
    dist : float
        Distance of the event in (km).

    """

    NAME = 'Boore & Thompson (2012)'
    ABBREV = 'BT12'

    def __init__(self, region, mag, dist, **kwargs):
        """Initialize the class."""
        BooreThompson.__init__(self, region, mag, dist, 'bt12', **kwargs)
        BooreJoyner1984.__init__(self, **kwargs)


class BooreThompson2015(BooreThompson, Vanmarcke1975):
    """Boore and Thompson (2015) peak factor.

    Peak calculation based on the peak factor definition by Vanmarcke
    (1975, :cite:`vanmarcke75`) along with the root-mean-squared duration
    correction proposed by Boore & Thompson (2015, :cite:`boore15`).


    Parameters
    ----------
    region : str
        Region for which the parameters were developed.  Valid options
        are: 'wna' for Western North America (active tectonic), and 'cena'
        for Central and Eastern North America ( stable tectonic).
    mag : float
        Magnitude of the event.
    dist : float
        Distance of the event in (km).

    """

    NAME = 'Boore & Thompson (2015)'
    ABBREV = 'BT15'

    def __init__(self, region, mag, dist, **kwargs):
        """Initialize the class."""
        BooreThompson.__init__(self, region, mag, dist, 'bt15', **kwargs)
        Vanmarcke1975.__init__(
            self, use_nonstationarity_factor=False, **kwargs)


class WangRathje2018(BooreThompson2015):
    """Wang & Rathje (2018) peak factor.

    Peak calculation based on the peak factor definition by Vanmarcke (1975,
    :cite:`vanmarcke75`) along with correction for oscillator duration and site
    amplification as described in Wang & Rathje (2018, :cite:`rathje18`).
    """

    NAME = 'Wang & Rathje (2018) '
    ABBREV = 'WR18'

    # Coefficients from Table 2, and paragraph after Equation (8)
    COEFS = np.rec.fromrecords(
        [(1, 0.2688, 0.0030, 1.8380, -0.0198, 0.091),
         (2, 0.2555, -0.0002, 1.2154, -0.0183, 0.081),
         (3, 0.2287, -0.0014, 0.9404, -0.0130, 0.056)],
        names='mode,a,b,d,e,sd', )

    def __init__(self, region, mag, dist, **kwargs):
        """Initialize the class."""
        BooreThompson2015.__init__(self, region, mag, dist, **kwargs)

    def _calc_duration_rms(self, duration, **kwargs):
        """Compute the RMS duration.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity (sec).
        osc_freq : float
            Frequency of the oscillator (Hz).
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for
            a damping ratio of 5%.
        site_tf : array_like
            Transfer function for applied to compute site effects.

        Returns
        -------
        duration_rms : float
            Duration of the root-mean-squared oscillator response (sec).

        """
        duration_rms = BooreThompson2015._calc_duration_rms(self, duration,
                                                            **kwargs)
        osc_freq = kwargs.get('osc_freq', None)

        site_tf = np.abs(kwargs.get('site_tf', []))
        if np.any(site_tf > 1):
            # Modify duration for site effects

            # Compute the expected rock oscillator duration

            # Equation 4a
            f_lim = 5.274 * duration ** -0.640
            ratio = 1
            if 0.1 <= osc_freq < f_lim:
                # Equation 4b
                dur_0 = 31.858 * duration ** -0.849
                # Equation 4c
                dur_min = 1.009 * duration / (3.583 + duration)
                # Equation 3b
                b = 1 / (dur_0 - dur_min)
                # Equation 3a
                a = (1 / (dur_0 - 1) - b) * (f_lim - 0.1)
                # Equation 2
                ratio = (dur_0 - (osc_freq - 0.1) / (a + b * (osc_freq - 0.1)))

            dur_osc_rock = ratio * duration

            # Compute the expected soil oscillator duration

            # Peaks in the transfer function
            indices = argrelmax(site_tf)[0][:3]

            freqs = kwargs['freqs']
            modes_f = freqs[indices]
            modes_a = site_tf[indices]

            # Amplitude / frequency ratio of the first mode
            af_ratio = modes_a[0] / modes_f[0]

            c = self.COEFS.a * af_ratio + self.COEFS.b * af_ratio ** 2
            m = self.COEFS.d * af_ratio + self.COEFS.e * af_ratio ** 2
            incr_max = c * np.exp(-duration / m)

            incr = incr_max * np.exp(
                -np.log(osc_freq / modes_f) ** 2 /
                (2 * self.COEFS.sd ** 2)
            )

            dur_osc_soil = dur_osc_rock + incr.sum()

            # Scale the RMS duration by the ratio in soil to rock durations
            duration_rms *= (dur_osc_soil / dur_osc_rock)

        return duration_rms


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
    calc : calculator
        :class:`.Calculator`

    """
    calc_kwds = calc_kwds or dict()

    calculators = [
        BooreJoyner1984,
        BooreThompson2012,
        BooreThompson2015,
        CartwrightLonguetHiggins1956,
        Davenport1964,
        DerKiureghian1985,
        LiuPezeshk1999,
        ToroMcGuire1987,
        Vanmarcke1975,
        WangRathje2018,
    ]

    for calculator in calculators:
        if method in [calculator.NAME, calculator.ABBREV]:
            return calculator(**calc_kwds)
    else:
        raise NotImplementedError('No calculator for: %s', method)


def get_region(region):
    """Return the region naming used in this package.

    Parameters
    ----------
    region : str
        Regional synonym.

    Returns
    -------
    region : str
        Region either 'cena' or 'wna'.

    """
    region = region.lower()
    if region in ['cena', 'ena', 'ceus', 'eus']:
        return 'cena'
    elif region in ['wna', 'wus']:
        return 'wna'
    else:
        raise NotImplementedError('No recognized region for: %s', region)
