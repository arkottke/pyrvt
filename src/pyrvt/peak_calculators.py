#!/usr/bin/python
"""
Peak factor models.

Published peak factor models, which compute the expected peak ground motion. A
specific model may include non-stationarity adjustments such as a oscillator
duration correction.
"""

import ctypes
import itertools
import pathlib
import warnings
from abc import ABC, abstractmethod
from typing import Any

import numba
import numpy as np
import numpy.typing as npt
from scipy.integrate import quad
from scipy.interpolate import LinearNDInterpolator
from scipy.signal import argrelmax

try:
    _trapz = np.trapezoid
except AttributeError:
    _trapz = np.trapz


# FIXME: Is this needed?
@numba.jit(nopython=True)
def trapz(x: npt.ArrayLike, y: npt.ArrayLike) -> float:
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
c_sig = numba.types.double(numba.types.intc, numba.types.CPointer(numba.types.double))


# Turn the integrand function into a ctypes function
@numba.cfunc(c_sig)
def _calc_vanmarcke1975_ccdf(n: int, a) -> float:
    """Calculate the Vanmarcke (1975) complementary CDF.

    Parameters
    ----------
    n : int
        Length of arguments
    a : list of floats
        Arguments specifying:
            - function value `x`
            - number of zero crossings
            - effective bandwdith

    Returns
    -------
    ccdf : float
        Complementary CDF value

    """
    args = numba.carray(a, n)
    x, num_z, bandwidth_eff = args

    fact = np.exp(-(x**2) / 2)

    return 1 - (
        (1 - fact)
        * np.exp(
            -num_z
            * fact
            * (1 - np.exp(-np.sqrt(np.pi / 2) * bandwidth_eff * x))
            / (1 - fact)
        )
    )


# Force the argtypes to be what quad expects
_calc_vanmarcke1975_ccdf.ctypes.argtypes = (ctypes.c_int, ctypes.c_double)


# Turn the integrand function into a ctypes function
@numba.cfunc(c_sig)
def _calc_log_vanmarcke1975_ccdf(n: int, a) -> float:
    """Calculate 1 - Vanmarcke (1975) complementary CDF.

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
    # Convert from log to natural
    a[0] = np.exp(a[0])

    value = _calc_vanmarcke1975_ccdf(n, a)
    if a[0] < 1:
        value = 1 - value
    return value


_calc_log_vanmarcke1975_ccdf.ctypes.argtypes = (ctypes.c_int, ctypes.c_double)


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
    return 1.0 - (1.0 - bandwidth * np.exp(-x * x)) ** num_extrema


# Force the argtypes to be what quad expects
_calc_cartwright_pf.ctypes.argtypes = (ctypes.c_int, ctypes.c_double)


def calc_moments(
    freqs: npt.ArrayLike, fourier_amps: npt.ArrayLike, orders: list[int]
) -> list[float]:
    """Compute the moments of a Fourier amplitude spectrumself.

    Parameters
    ----------
    freqs : array_like
        Frequency of the Fourier amplitude spectrum (Hz)
    fourier_amps : array_like
        Amplitude of the Fourier amplitude spectrum.
    orders : list
        Moments to consider

    Returns
    -------
    moment : list
        Computed spectral moments.
    """
    squared_fa = np.square(fourier_amps)

    # Use trapzoidal integration to compute the requested moments.
    moments = [
        2.0 * _trapz(np.power(2 * np.pi * freqs, o) * squared_fa, x=freqs)
        for o in orders
    ]

    return moments


class SquaredSpectrum:
    """Squared Fourier amplitude spectrum.

    Used to store calculated spectral moments during calculations.

    Parameters
    ----------
    freqs : array_like
        Frequency of the Fourier amplitude spectrum (Hz)
    fourier_amps : array_like
        Amplitude of the Fourier amplitude spectrum.
    """

    def __init__(self, freqs: npt.ArrayLike, fourier_amps: npt.ArrayLike):
        """Initialize SquaredSpectrum."""
        self._freqs = np.asarray(freqs)
        self._squared_fa = np.square(fourier_amps)
        self._moments = {}

    def moment(self, num: int) -> float:
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
            moment = 2.0 * _trapz(
                np.power(2 * np.pi * self._freqs, num) * self._squared_fa, x=self._freqs
            )
            self._moments[num] = moment

        return moment

    def moments(self, *nums) -> list[float]:
        """Return the computed moments.

        Returns
        -------
        moments : list[float]
            Computed spectral moments.
        """
        return [self.moment(n) for n in nums]

    @property
    def freqs(self) -> np.ndarray:
        return self._freqs

    @property
    def squared_fa(self) -> np.ndarray:
        return self._squared_fa


class Calculator(ABC):
    """Base class used for all peak calculator classes.

    Provides the interface that is used by all the other classes. Specific peak
    calculators modify the `_calc_peak_factor` method and potentially
    `_calc_duration_rms` method. Using any calculator is done through calling of
    `__call__`.

    Attributes
    ----------
    NAME : str
        Complete reference of the peak calculator

    ABBREV : str
        Abbreviation of the reference

    _MIN_ZERO_CROSSINGS : float
        Minimum number of zero crossings.
    """

    NAME: str = ""
    ABBREV: str = ""
    # FIXME: Provide justification for this number
    _MIN_ZERO_CROSSINGS: float = 1.33

    def __init__(self, **kwds):
        """Initialize the object."""
        super().__init__()
        self._spectrum = None

    @classmethod
    def limited_num_zero_crossings(cls, num_zero_crossings: float) -> float:
        """Limit the number of zero crossings to a static limit.

        The minimum is configured through adjusting `cls._MIN_ZERO_CROSSINGS`.

        Parameters
        ----------
        num_zero_crossings : float
            Calculated number of zero crossings

        Returns
        -------
        float
            Limited number of zero crossings
        """
        return max(cls._MIN_ZERO_CROSSINGS, num_zero_crossings)

    def __call__(
        self,
        duration: float,
        freqs: npt.ArrayLike,
        fourier_amps: npt.ArrayLike,
        **kwds,
    ) -> tuple[float]:
        """Compute the peak response.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
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
        sspectrum = SquaredSpectrum(freqs, fourier_amps)
        peak_factor = self._calc_peak_factor(duration, sspectrum, **kwds)
        duration_rms = self._calc_duration_rms(duration, sspectrum, **kwds)
        # Compute the root-mean-squared response.
        resp_rms = np.sqrt(sspectrum.moment(0) / duration_rms)

        return peak_factor * resp_rms, peak_factor

    @abstractmethod
    def _calc_peak_factor(self, duration: float, sspectrum: SquaredSpectrum) -> float:
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion. Typically
            defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.

        Returns
        -------
        peak_factor : float
            associated peak factor.

        """

    def _calc_duration_rms(
        self, duration: float, sspectrum: SquaredSpectrum, **kwds
    ) -> float:
        """Modify root-mean-squared duration to account for nonstationarity.

        Default implemenation does nothing.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the  ground motion [sec].
            Typically defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.
        Returns
        -------
        duration : float
            Modified duration.
        """
        return duration


class Vanmarcke1975(Calculator):
    r"""Vanmarcke (1975) peak factor.

    The Vanmarcke (1975, Vanmarcke (1975)) peak factor, which includes the effects
    of clumping.  The peak factor equation is from Equation (2) in
    Der Kiureghian (1980), which is based on Equation (29) in Vanmarcke (1975).

    The cumulative density function (CDF) of the peak is defined as:

    .. math::
        F_x(x) = \left[1 - \exp\left(-x^2/2\right)\right]
        \exp\left[-N_z \frac{1 -
            \exp\left(-\sqrt{\pi/2} \delta_e x\right)}{\exp(x^2 / 2) -
            1 }\right]

    where $N_z$ is the number of zero crossings, $\delta_e$ is the effective
    bandwidth ($\delta^{1.2}$).

    Typically, the expected value of the peak factor is calculated by integrating over
    the probability density function (i.e., $f_x(x) = \frac{d}{dx} F_x( x)$):

    .. math::
        E[x] = \int_0^\infty x f_x(x) dx

    However, because of the properties of $F_x(x)$, specifically that it has
    non-zero probabilities for only positive values, $E[x]$ can be computed
    directly from $F_x(x)$.

    .. math::
        E[x] = \int_0^\infty 1 - F_x(x) dx.

    This is based on the following sources [1] and [2].

    [1]: http://en.wikipedia.org/wiki/Expected_value#Formulas_for_special_cases
    [2]: http://stats.stackexchange.com/a/13377/48461

    Attributes
    ----------
    NAME : str
        Complete reference of the peak calculator

    ABBREV : str
        Abbreviation of the reference

    _MIN_ZERO_CROSSINGS : float
        Minimum number of zero crossings.
    """

    NAME: str = "Vanmarcke (1975)"
    ABBREV: str = "V75"

    def __init__(self, **kwds):
        """Initialize the class."""
        super().__init__(**kwds)

    def _calc_peak_factor(
        self, duration: float, sspectrum: SquaredSpectrum, **kwds
    ) -> float:
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the  ground motion [sec].
            Typically defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.
        osc_freq : float
            Frequency of the oscillator (Hz).
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for a damping
            ratio of 5%.
        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        m0, m1, m2 = sspectrum.moments(0, 1, 2)

        bandwidth = np.sqrt(1 - (m1 * m1) / (m0 * m2))
        bandwidth_eff = bandwidth**1.2

        num_zero_crossings = self.limited_num_zero_crossings(
            duration * np.sqrt(m2 / m0) / np.pi
        )
        # The expected peak factor is computed as the integral of the
        # complementary CDF (1 - CDF(x)).
        peak_factor = quad(
            _calc_vanmarcke1975_ccdf.ctypes,
            0,
            np.inf,
            args=(num_zero_crossings, bandwidth_eff),
        )[0]

        return peak_factor


class Davenport1964(Calculator):
    """Davenport (1964) peak factor.

    RVT calculation using the asymptotic solution proposed by Davenport (1964).

    Attributes
    ----------
    NAME : str
        Complete reference of the peak calculator

    ABBREV : str
        Abbreviation of the reference

    _MIN_ZERO_CROSSINGS : float
        Minimum number of zero crossings.
    """

    #:  Name of the calculator
    NAME: str = "Davenport (1964)"
    #: Abbreviation of the calculator
    ABBREV: str = "D64"

    def __init__(self, **kwds):
        """Initialize the class."""
        super().__init__(**kwds)

    def _calc_peak_factor(
        self, duration: float, sspectrum: SquaredSpectrum, **kwds
    ) -> float:
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the  ground motion [sec].
            Typically defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.
        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        m0, m2 = sspectrum.moments(0, 2)

        # Compute the number of zero crossings
        num_zero_crossings = self.limited_num_zero_crossings(
            duration * np.sqrt(m2 / m0) / np.pi
        )

        peak_factor = self.asymtotic_approx(num_zero_crossings)

        return peak_factor

    @classmethod
    def asymtotic_approx(cls, zero_crossings: float) -> float:
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

    RVT calculation using peak factor derived by Davenport (1964) with limits
    suggested by :cite:t:`igusa85`.

    Attributes
    ----------
    NAME : str
        Complete reference of the peak calculator

    ABBREV : str
        Abbreviation of the reference

    _MIN_ZERO_CROSSINGS : float
        Minimum number of zero crossings.
    """

    #:  Name of the calculator
    NAME: str = "Der Kiureghian (1985)"
    #: Abbreviation of the calculator
    ABBREV: str = "DK85"

    def __init__(self, **kwds):
        """Initialize the class."""
        super().__init__(**kwds)

    def _calc_peak_factor(
        self, duration: float, sspectrum: SquaredSpectrum, **kwds
    ) -> float:
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the  ground motion [sec].
            Typically defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.

        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        m0, m1, m2 = sspectrum.moments(0, 1, 2)

        # Compute the number of zero crossings
        num_zero_crossings = duration * np.sqrt(m2 / m0) / np.pi

        # Reduce the rate of zero crossings based on the bandwidth
        bandwidth = np.sqrt(1 - (m1 * m1) / (m0 * m2))
        if bandwidth <= 0.1:
            eff_crossings = max(2.1, 2 * bandwidth * num_zero_crossings)
        elif 0.1 < bandwidth <= 0.69:
            eff_crossings = (1.63 * bandwidth**0.45 - 0.38) * num_zero_crossings
        else:
            eff_crossings = num_zero_crossings

        eff_crossings = self.limited_num_zero_crossings(eff_crossings)
        peak_factor = self.asymtotic_approx(eff_crossings)

        return peak_factor


class ToroMcGuire1987(Davenport1964):
    """Toro and McGuire (1987) peak factor.

    Peak factor equation using asymptotic solution proposed by Davenport (1964)
    with modifications proposed by Toro & McGuire (1987).

    Parameters
    ----------
    use_nonstationarity_factor : bool
        If the non-stationarity factor should be applied.

    Attributes
    ----------
    NAME : str
        Complete reference of the peak calculator

    ABBREV : str
        Abbreviation of the reference

    _MIN_ZERO_CROSSINGS : float
        Minimum number of zero crossings.
    """

    #:  Name of the calculator
    NAME: str = "Toro & McGuire (1987)"
    #: Abbreviation of the calculator
    ABBREV: str = "TM87"

    def __init__(self, use_nonstationarity_factor: bool = True, **kwds):
        """Initialize the class."""
        super().__init__(**kwds)

        self._use_nonstationarity_factor = use_nonstationarity_factor

    def _calc_peak_factor(
        self, duration: float, sspectrum: SquaredSpectrum, **kwds
    ) -> float:
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the  ground motion [sec].
            Typically defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.

        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        m0, m1, m2 = sspectrum.moments(0, 1, 2)

        # Vanmarcke's (1976) bandwidth measure and central frequency
        bandwidth = np.sqrt(1 - (m1 * m1) / (m0 * m2))
        freq_cent = np.sqrt(m2 / m0) / (2 * np.pi)

        num_zero_crossings = self.limited_num_zero_crossings(
            2 * freq_cent * duration * (1.63 * bandwidth**0.45 - 0.38)
        )

        peak_factor = self.asymtotic_approx(num_zero_crossings)

        osc_freq = kwds.get("osc_freq", None)
        osc_damping = kwds.get("osc_damping", None)
        if osc_freq and osc_damping:
            peak_factor *= self.nonstationarity_factor(osc_damping, osc_freq, duration)

        return peak_factor

    @classmethod
    def nonstationarity_factor(
        cls, osc_damping: float, osc_freq: float, duration: float
    ) -> float:
        """Compute nonstationarity factor to the peak response.

        Vanmarcke (1975) provides a recommendation for scaling the damping of the SDOF
        oscillator to account for nonstationarity in Equation 8.30:
        $$
        \\xi_t = \\frac{\\xi}{1 - \\exp\\left(-4 \\pi \\xi f_n t\\right)}
        $$
        Toro and McGuire (1987) simplified this to scale the $X_{rms}$ by:
        $$
        n_f = \\sqrt{1 - \\exp\\left(-4 \\pi \\xi f_n T\\right)}
        $$

        The simplified model from Toro and McGuire (1987) is used here.


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
        fact = -4 * np.pi * osc_damping * osc_freq * duration

        # Here for some conditions the nonstationarity factor gets very small (-700)
        # which causes and underflow the in exponent below.
        if fact < -10:
            nsf = 1
        else:
            nsf = np.sqrt(1 - np.exp(fact))

        return nsf


class CartwrightLonguetHiggins1956(Calculator):
    """Cartwight and Longuet-Higgins (1956) peak factor.

    RVT calculation based on the peak factor definition by :cite:t:`cartwright56` using
    the integral provided by :cite:t:`boore03`.

    Attributes
    ----------
    NAME : str
        Complete reference of the peak calculator

    ABBREV : str
        Abbreviation of the reference

    _MIN_ZERO_CROSSINGS : float
        Minimum number of zero crossings.
    """

    #: Name of the calculator
    NAME: str = "Cartwright & Longuet-Higgins (1956)"
    #: Abbreviation of the calculator
    ABBREV: str = "CLH56"

    def __init__(self, **kwds):
        """Initialize the class."""
        super().__init__(**kwds)

    def _calc_peak_factor(
        self, duration: float, sspectrum: SquaredSpectrum, **kwds
    ) -> float:
        """Compute the peak factor.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the  ground motion [sec].
            Typically defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.

        Returns
        -------
        peak_factor : float
            associated peak factor.

        """
        m0, m2, m4 = sspectrum.moments(0, 2, 4)

        bandwidth = np.sqrt((m2 * m2) / (m0 * m4))
        num_extrema = max(2.0, np.sqrt(m4 / m2) * duration / np.pi)
        # Compute the peak factor by the indefinite integral.
        peak_factor = (
            np.sqrt(2.0)
            * quad(
                _calc_cartwright_pf.ctypes, 0, np.inf, args=(num_extrema, bandwidth)
            )[0]
        )

        return peak_factor


class BooreJoyner1984(CartwrightLonguetHiggins1956):
    """Boore and Joyner (1984) peak factor.

    RVT calculation based on the peak factor definition by :cite:t:`cartwright56` and
    along with the root-mean-squared duration correction proposed by :cite:t:`boore84`.

    This RVT calculation is used by SMSIM and is described in :cite:t:`boore03`.

    Attributes
    ----------
    NAME : str
        Complete reference of the peak calculator

    ABBREV : str
        Abbreviation of the reference

    _MIN_ZERO_CROSSINGS : float
        Minimum number of zero crossings.
    """

    #:  Name of the calculator
    NAME: str = "Boore & Joyner (1984)"
    #: Abbreviation of the calculator
    ABBREV: str = "BJ84"

    def __init__(self, **kwds):
        """Initialize the class."""
        super().__init__(**kwds)

    def _calc_duration_rms(
        self,
        duration: float,
        sspectrum: SquaredSpectrum,
        *,
        osc_damping: float = 0.05,
        osc_freq: float | None = None,
        **kwds,
    ) -> float:
        """Modify root-mean-squared duration to account for nonstationarity.

        Default implemenation does nothing.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the  ground motion [sec].
            Typically defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for a damping
            ratio of 5%.
        osc_freq : float
            Oscillator frequency [Hz].
        Returns
        -------
        duration : float
            Modified duration.
        """
        if osc_damping and osc_freq:
            power = 3.0
            coef = 1.0 / 3.0
            # This equation was rewritten in Boore and Thompson (2012).
            foo = 1.0 / (osc_freq * duration)
            dur_ratio = 1 + 1.0 / (2 * np.pi * osc_damping) * (
                foo / (1 + coef * foo**power)
            )
            duration *= dur_ratio

        return duration


class LiuPezeshk1999(BooreJoyner1984):
    """Liu and Pezeshk (1999) peak factor.

    RVT calculation based on the peak factor definition by :cite:t:`cartwright56` along
    with the root-mean-squared duration correction proposed by :cite:t:`liu99`.
    """

    #:  Name of the calculator
    NAME: str = "Liu & Pezeshk (1999)"
    #: Abbreviation of the calculator
    ABBREV: str = "LP99"

    def __init__(self, **kwds):
        """Initialize the class."""
        super().__init__(**kwds)

    def _calc_duration_rms(
        self,
        duration: float,
        sspectrum: SquaredSpectrum,
        *,
        osc_damping: float = 0.05,
        osc_freq: float | None = None,
        **kwds,
    ) -> float:
        """Modify root-mean-squared duration to account for nonstationarity.

        Default implemenation does nothing.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the  ground motion [sec].
            Typically defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for a damping
            ratio of 5%.
        osc_freq : float
            Oscillator frequency [Hz].
        Returns
        -------
        duration : float
            Modified duration.
        """
        if osc_freq and osc_damping:
            m0, m1, m2 = sspectrum.moments(0, 1, 2)

            power = 2.0
            coef = np.sqrt(2 * np.pi * (1.0 - (m1 * m1) / (m0 * m2)))

            # Same model as used in Boore and Joyner (1984). This equation was
            # rewritten in Boore and Thompson (2012).
            foo = 1.0 / (osc_freq * duration)
            dur_ratio = 1 + 1.0 / (2 * np.pi * osc_damping) * (
                foo / (1 + coef * foo**power)
            )
            duration *= dur_ratio

        return duration


def _make_bt_interpolator(region: str, ref: str) -> LinearNDInterpolator:
    """Load data from the :cite:t:`boore12` and Boore & Thompson (2015) parameter files.

    Parameters
    ----------
    region : str
        Region for which the parameters were developed. Valid options: 'wna' for Western
        North America (active tectonic), or 'cena' for Central and Eastern North America
        (stable tectonic).
    ref : str
        Reference document. Either: 'bt12' or 'bt15'.

    Returns
    -------
    LinearNDInterpolator
        Interpolator for the data.
    """
    fpath = pathlib.Path(__file__).parent.joinpath(
        "data", f"{region}_{ref}_trms4osc.pars.gz"
    )
    d = np.rec.fromrecords(
        np.loadtxt(str(fpath), skiprows=4, usecols=range(9)),
        names="mag,dist,c1,c2,c3,c4,c5,c6,c7",
    )

    return LinearNDInterpolator(
        np.c_[d.mag, np.log(d.dist)], np.c_[d.c1, d.c2, d.c3, d.c4, d.c5, d.c6, d.c7]
    )


# Load coefficient interpolators for Boore and Thompson (2012)
_BT_INTERPS = {
    (region, ref): _make_bt_interpolator(region, ref)
    for region, ref in itertools.product(["wna", "cena"], ["bt12", "bt15"])
}


class BooreThompson(Calculator):
    """Abstract class for the Boore & Thompson duration correction.

    The duration ratio is defined by Equation (10) in :cite:t:`boore12`. Magnitude and
    distance is interpolated using `scipy.interpolate.LinearNDInterpolator` on the
    natural log of the distance.

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

    def __init__(self, region, mag, dist, ref, **kwds):
        """Initialize the class."""
        super().__init__(**kwds)
        region = get_region(region)
        self._COEFS = _BT_INTERPS[(region, ref)](mag, np.log(dist))

    def _calc_duration_rms(
        self,
        duration: float,
        sspectrum: SquaredSpectrum,
        osc_damping: float = 0.05,
        osc_freq: float | None = None,
        **kwds,
    ) -> float:
        """Modify root-mean-squared duration to account for nonstationarity.

        Default implemenation does nothing.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the  ground motion [sec].
            Typically defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for a damping
            ratio of 5%.
        osc_freq : float
            Oscillator frequency [Hz].
        Returns
        -------
        duration : float
            Modified duration.
        """
        if osc_freq and osc_damping:
            c1, c2, c3, c4, c5, c6, c7 = self._COEFS

            foo = 1 / (osc_freq * duration)
            dur_ratio = (c1 + c2 * (1 - foo**c3) / (1 + foo**c3)) * (
                1 + c4 / (2 * np.pi * osc_damping) * (foo / (1 + c5 * foo**c6)) ** c7
            )
            duration *= dur_ratio

        return duration


class BooreThompson2012(BooreThompson, BooreJoyner1984):
    """Boore and Thompson (2012) peak factor.

    Peak calculation based on the peak factor definition by :cite:t:`cartwright56`
    along with the root-mean-squared duration correction proposed by :cite:t:`boore12`.


    Parameters
    ----------
    region : str
        Region for which the parameters were developed.  Valid options are: 'wna' for
        Western North America (active tectonic), and 'cena' for Central and Eastern
        North America ( stable tectonic).
    mag : float
        Magnitude of the event.
    dist : float
        Distance of the event in (km).

    Attributes
    ----------
    NAME : str
        Complete reference of the peak calculator

    ABBREV : str
        Abbreviation of the reference

    _MIN_ZERO_CROSSINGS : float
        Minimum number of zero crossings.

    """

    #:  Name of the calculator
    NAME: str = "Boore & Thompson (2012)"
    #: Abbreviation of the calculator
    ABBREV: str = "BT12"

    def __init__(self, region, mag, dist, **kwds):
        """Initialize the class."""
        BooreThompson.__init__(self, region, mag, dist, "bt12", **kwds)
        BooreJoyner1984.__init__(self, **kwds)


class BooreThompson2015(BooreThompson, Vanmarcke1975):
    """Boore and Thompson (2015) peak factor.

    Peak calculation based on the peak factor definition by Vanmarcke (1975) along
    with the root-mean-squared duration correction proposed by Boore & Thompson (2015).


    Parameters
    ----------
    region : str
        Region for which the parameters were developed.  Valid options are: 'wna' for
        Western North America (active tectonic), and 'cena' for Central and Eastern
        North America ( stable tectonic).
    mag : float
        Magnitude of the event.
    dist : float
        Distance of the event in (km).

    Attributes
    ----------
    NAME : str
        Complete reference of the peak calculator

    ABBREV : str
        Abbreviation of the reference

    _MIN_ZERO_CROSSINGS : float
        Minimum number of zero crossings.

    """

    #:  Name of the calculator
    NAME: str = "Boore & Thompson (2015)"
    #: Abbreviation of the calculator
    ABBREV: str = "BT15"

    def __init__(self, region, mag, dist, **kwds):
        """Initialize the class."""
        BooreThompson.__init__(self, region, mag, dist, "bt15", **kwds)
        Vanmarcke1975.__init__(self, **kwds)


class WangRathje2018(BooreThompson2015):
    """Wang & Rathje (2018) peak factor.

    Peak calculation based on the peak factor definition by Vanmarcke (1975) along
    with correction for oscillator duration by Boore & Thompson (2015) and site
    amplification as described in Wang & Rathje (2018).

    Parameters
    ----------
    region : str
        Region for which the parameters were developed.  Valid options are: 'wna' for
        Western North America (active tectonic), and 'cena' for Central and Eastern
        North America ( stable tectonic).
    mag : float
        Magnitude of the event.
    dist : float
        Distance of the event in (km).

    Attributes
    ----------
    NAME : str
        Complete reference of the peak calculator

    ABBREV : str
        Abbreviation of the reference

    _MIN_ZERO_CROSSINGS : float
        Minimum number of zero crossings.
    """

    #:  Name of the calculator
    NAME: str = "Wang & Rathje (2018) "
    #: Abbreviation of the calculator
    ABBREV: str = "WR18"

    # Coefficients from Table 2, and paragraph after Equation (8)
    COEFS = np.rec.fromrecords(
        [
            (1, 0.2688, 0.0030, 1.8380, -0.0198, 0.091),
            (2, 0.2555, -0.0002, 1.2154, -0.0183, 0.081),
            (3, 0.2287, -0.0014, 0.9404, -0.0130, 0.056),
        ],
        names="mode,a,b,d,e,sd",
    )

    def __init__(self, region, mag, dist, **kwds):
        """Initialize the class."""
        BooreThompson2015.__init__(self, region, mag, dist, **kwds)

    def _calc_duration_rms(
        self,
        duration: float,
        sspectrum: SquaredSpectrum,
        *,
        osc_damping: float = 0.05,
        osc_freq: float | None = None,
        site_tf: npt.ArrayLike | None = None,
        **kwds,
    ) -> float:
        """Modify root-mean-squared duration to account for nonstationarity.

        Default implemenation does nothing.

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the  ground motion [sec].
            Typically defined as the duration between the 5% and 75% normalized Arias
            intensity [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` that defines the frequency content of the
            motion.
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for a damping
            ratio of 5%.
        osc_freq : float
            Oscillator frequency [Hz].
        site_tf : array_like
            Transfer function for applied to compute site effects.
        Returns
        -------
        duration : float
            Modified duration.
        """
        duration_rms = BooreThompson2015._calc_duration_rms(
            self, duration, sspectrum, osc_damping=osc_damping, osc_freq=osc_freq
        )

        if osc_freq and osc_damping and site_tf is not None:
            # Modify duration for site effects

            # Work only on the absolute value
            site_tf = np.abs(site_tf)
            # Peaks in the transfer function
            indices = argrelmax(site_tf)[0][:3]

            # Some transfer functions might have have peaks
            if len(indices):
                # In some instances higher-modes are not found -- especially due to high
                # damping values. In this case, we truncate the coefficients
                C = self.COEFS[: len(indices)]

                # if len(indices) != 3:
                #     indices = np.r_[
                #         indices,
                #         (3 - len(indices))
                #         * [
                #             10,
                #         ],
                #     ]
                #     valid_modes = len(indices)
                # else:
                #     valid_modes = 3

                modes_f = sspectrum.freqs[indices]
                modes_a = site_tf[indices]

                # Amplitude / frequency ratio of the first mode
                af_ratio = modes_a[0] / modes_f[0]

                c = C.a * af_ratio + C.b * af_ratio**2
                m = C.d * af_ratio + C.e * af_ratio**2
                incr_max = c * np.exp(-duration / m)

                incr = incr_max * np.exp(
                    -((np.log(osc_freq / modes_f)) ** 2) / (2 * C.sd**2)
                )
                duration_rms += incr.sum()
                # duration_rms += incr[0:valid_modes].sum()

        return duration_rms


class SeifriedEtAl2025(Calculator):
    """Seifried et al. (2025) peak factor calculator."""

    #:  Name of the calculator

    NAME: str = "Seifried et al. (2025)"
    #: Abbreviation of the calculator
    ABBREV: str = "Sea25"

    # Coefficients from M. Bahrampouri (7/19/2025)
    _COEF_A: float = 0.541
    _COEF_B: float = 2.456

    def __init__(
        self,
        use_nonstationarity_factor: bool = True,
        mean_calc="arithmetic",
        **kwds: Any,
    ) -> None:
        """Initialize SeifriedEtAl2025.

        Parameters
        ----------
        use_nonstationarity_factor : bool, optional
            If True, the nonstationarity adjustment is applied, by default True.
        **kwds : Any
            Additional keyword arguments.
        """
        super().__init__(**kwds)
        self._use_nonstationarity_factor = use_nonstationarity_factor
        self._mean_calc = mean_calc

        if self._mean_calc not in ["arithmetic", "geometric"]:
            raise ValueError(
                f"Mean calculation method '{self._mean_calc}' is not supported. "
                "Use 'arithmetic' or 'geometric'."
            )

        if self._mean_calc == "geometric":
            warnings.warn(
                "Coefficients A and Be developed for airithmetic mean calculation."
            )

    def _calc_peak_factor(
        self, duration: float, sspectrum: SquaredSpectrum, **kwds: Any
    ) -> float:
        """Compute the peak factor for Seifried et al. (2025).

        Parameters
        ----------
        duration : float
            Duration of the stationary portion of the ground motion [sec].
        sspectrum : SquaredSpectrum
            Instance of `SquaredSpectrum` defining frequency content.
        **kwds : Any
            May contain 'osc_freq' and 'osc_damping' parameters if relevant.

        Returns
        -------
        float
            Computed peak factor.
        """
        m0, m2 = sspectrum.moments(0, 2)

        tE = 4 * _trapz(np.square(sspectrum.squared_fa), x=sspectrum.freqs) / m0**2

        # Central frequency
        freq_cent = np.sqrt(m2 / m0) / (2 * np.pi)
        # Effective bandwidth from Winterstein
        bandwidth_eff = np.sqrt(2 / (freq_cent * tE)) / np.pi
        # Effective damping
        damp_eff = 1 / (2 * np.pi * freq_cent * tE)
        # Number of zero crossings
        num_zero_crossings = self.limited_num_zero_crossings(2 * duration * freq_cent)

        # Integration done in log-space. Need to seperate the parts below
        # and above ln(1)
        if self._mean_calc == "geometric":
            left = quad(
                _calc_log_vanmarcke1975_ccdf.ctypes,
                # FIXME Better bounds?
                -5,
                0,
                args=(num_zero_crossings, bandwidth_eff),
            )[0]
            right = quad(
                _calc_log_vanmarcke1975_ccdf.ctypes,
                0,
                # FIXME Better bounds?
                5,
                args=(num_zero_crossings, bandwidth_eff),
            )[0]
            peak_factor = np.exp(-left + right)
        elif self._mean_calc == "arithmetic":
            # The expected peak factor is computed as the integral of the
            # complementary CDF (1 - CDF(x)).
            peak_factor = quad(
                _calc_vanmarcke1975_ccdf.ctypes,
                0,
                np.inf,
                args=(num_zero_crossings, bandwidth_eff),
            )[0]
        else:
            raise NotImplementedError(
                f"Mean calculation method '{self._mean_calc}' is not implemented."
            )

        if self._use_nonstationarity_factor:
            peak_factor *= np.sqrt(
                1
                - np.exp(
                    (-1 * (4 * np.pi * damp_eff * freq_cent * duration) ** self._COEF_A)
                    / self._COEF_B
                )
            )

        return peak_factor


def get_peak_calculator(method: str, calc_kwds: dict[str, Any] | None) -> Calculator:
    """Select a peak calculator based on a XXDD string.

    The format of the string is XX for author initials, and then DD for the last two
    years of the date published (e.g., 'BJ84' for Boore & Joyner 1984).

    Parameters
    ----------
    method : str
        Name or abbreviation of the peak calculation method.
    calc_kwds : dict[str, Any] | None
        Additional keywords passed to the calculator constructor.

    Returns
    -------
    Calculator
        A matching peak calculator instance.

    Raises
    ------
    NotImplementedError
        If no matching calculator is found.
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
        SeifriedEtAl2025,
    ]

    for calculator in calculators:
        if method in [calculator.NAME, calculator.ABBREV]:
            return calculator(**calc_kwds)
    else:
        raise NotImplementedError("No calculator for: %s", method)


def get_region(region: str) -> str:
    """Return the internal region naming used in this package.

    Parameters
    ----------
    region : str
        Regional synonym.

    Returns
    -------
    str
        Region either 'cena' or 'wna'.

    Raises
    ------
    NotImplementedError
        If the region is unknown.
    """
    region = region.lower()
    if region in ["cena", "ena", "ceus", "eus"]:
        return "cena"
    elif region in ["wna", "wus"]:
        return "wna"
    else:
        raise NotImplementedError("No recognized region for: %s", region)
