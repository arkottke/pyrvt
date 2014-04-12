#!/usr/bin/env python
# encoding: utf-8

"""Published peak factor models, which compute the expected peak ground
motion. A specific model may include oscillator duration correction."""

import os

import numpy as np

from scipy import Inf
from scipy.integrate import quad
from scipy.interpolate import LinearNDInterpolator


def compute_moments(freqs, fourier_amps, orders):
    """Compute the spectral moments.

    The spectral moment is computed using the squared Fourier amplitude
    spectrum.

    Parameters
    ----------
    freqs : numpy.array
        Frequency of the Fouier amplitude spectrum [Hz]

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


class DerKiureghian1983(object):
    """RVT calculation using peak factor derived by Davenport (1964) with
    limits suggested by Kiureghian and Neuenhofer [1]_.

    References
    ----------
    .. [1] Kiureghian, A. D., & Neuenhofer, A. (1992). Response spectrum method
    for multi‐support seismic excitations. Earthquake Engineering & Structural
    Dynamics, 21(8), 713-740.

    """
    def __init__(self, **kwds):
        # No class variables are required
        pass

    def __call__(self, gm_duration, freqs, fourier_amps, **kwds):
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
            Amplitude of the Fourie amplitude spectrum with a single
            degree of freedom oscillator already applied if being used. Units
            are not important.
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Damping of the oscillator [decimal]. For example, 0.05 for 5%.

        Returns
        -------
        max_resp : float
            Expected maximum response
        """

        m0, m2 = compute_moments(freqs, fourier_amps, [0, 2])

        # Compute the peak factor
        foo = gm_duration * np.sqrt(m2 / m0) / np.pi
        if foo < 1.6:
            # To avoid unreasonable values of the peak factor for small
            # frequencies, peak_factor of 1.56 is assumed for foo less than
            # 1.56. This is based on the recommendation from Der Kiureghian
            # (1992).
            peak_factor = 1.56
        else:
            bar = np.sqrt(2 * np.log(foo))
            peak_factor = bar + 0.577 / bar

        # Compute the root-mean-squared response
        resp_rms = np.sqrt(m0 / gm_duration)

        return peak_factor * resp_rms


class ToroMcGuire1987(object):
    """RVT calculation using peak factor derived by Davenport (1964) with
    modifications proposed by Toro and McGuire [1]_.

    References
    ----------
    .. [1] Toro, G. R., & McGuire, R. K. (1987). An investigation into
    earthquake ground motion characteristics in eastern North America. Bulletin
    of the Seismological Society of America, 77(2), 468-489.
    """
    def __init__(self, **kwds):
        # No class variables are required
        pass

    def __call__(self, gm_duration, freqs, fourier_amps, osc_freq=None,
                 osc_damping=None):
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
            Amplitude of the Fourie amplitude spectrum with a single
            degree of freedom oscillator already applied if being used. Units
            are not important.
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Damping of the oscillator [decimal]. For example, 0.05 for 5%.

        Returns
        -------
        max_resp : float
            Expected maximum response
        """

        m0, m1, m2 = compute_moments(freqs, fourier_amps, [0, 1, 2])

        # Vanmarcke's (1976) bandwidth measure and central frequency
        bandwidth = np.sqrt(1 - m1 ** 2 / (m0 * m2))
        freq_cent = np.sqrt(m2 / m0) / (2 * np.pi)

        zero_crossings = max(
            2 * freq_cent * gm_duration * (1.63 * bandwidth ** 0.45 - 0.38),
            1.33)

        foo = np.sqrt(2 * np.log(zero_crossings))
        peak_factor = (foo + 0.577 / foo)

        if osc_freq and osc_damping:
            peak_factor *= np.sqrt(
                1 - np.exp(-2 * osc_damping * osc_freq * gm_duration))

        # Compute the root-mean-squared response
        resp_rms = np.sqrt(m0 / gm_duration)

        return peak_factor * resp_rms


class BooreJoyner1984(object):
    """RVT calculation based on the peak factor definition by Cartwright and
    Longuet-Higgins (1956) [1]_ along with the root-mean-squared duration
    correction proposed by Boore and Joyner (1984) [2]_.

    This RVT calculation is used by SMSIM and is described in Boore(2003) [3]_.

    References
    ----------
    .. [1] Cartwright, D. E., & Longuet-Higgins, M. S. (1956). The statistical
    distribution of the maxima of a random function. Proceedings of the Royal
    Society of London. Series A. Mathematical and Physical Sciences, 237(1209),
    212-232.
    .. [2] Boore, D. M., & Joyner, W. B. (1984). A note on the use of random
    vibration theory to predict peak amplitudes of transient signals. Bulletin
    of the Seismological Society of America, 74(5), 2035-2039.
    .. [3] Boore, D. M. (2003). Simulation of ground motion using the
    stochastic method. Pure and applied geophysics, 160(3-4), 635-676.

    """
    def __call__(self, gm_duration, freqs, fourier_amps, osc_freq=None,
                 osc_damping=None, **kwds):
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
            Amplitude of the Fourie amplitude spectrum with a single
            degree of freedom oscillator already applied if being used. Units
            are not important.
        osc_freq : float
            Frequency of the oscillator [Hz]
        osc_damping : float
            Fractional damping of the oscillator. For example, 0.05 for 5%
            damping.

        Returns
        -------
        max_resp : float
            Expected maximum response

        """

        m0, m1, m2, m4 = compute_moments(freqs, fourier_amps, [0, 1, 2, 4])

        bandwidth = np.sqrt((m2 * m2) / (m0 * m4))
        num_extrema = max(2., np.sqrt(m4 / m2) * gm_duration / np.pi)

        # Compute the peak factor by the indefinite integral
        peak_factor = np.sqrt(2.) * quad(
            lambda z: 1. - (1. - bandwidth * np.exp(-z * z)) ** num_extrema,
            0, Inf)[0]

        # Compute the root-mean-squared response
        if osc_freq is None:
            rms_duration = gm_duration
        else:
            rms_duration = self.compute_duration_rms(
                gm_duration, osc_freq, osc_damping, m0, m1, m2)

        resp_rms = np.sqrt(m0 / rms_duration)

        return peak_factor * resp_rms

    def compute_duration_rms(self, gm_duration, osc_freq, osc_damping, *args,
                             **kwds):
        """Compute the oscillator duration used in the calculation of the
        root-mean-squared response.

        Based on Boore and Joyner (1984).

        Paramaters
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
        dur_ratio = (1 + 1. / (2 * np.pi * osc_damping)
                     * (foo / (1 + coef * foo ** power)))

        return gm_duration * dur_ratio


class LiuPezeshk1999(BooreJoyner1984):
    """ RVT calculation based on the peak factor definition by Cartwright and
    Longuet-Higgins (1956) [1]_ along with the root-mean-squared duration
    correction proposed by Liu and Pezeshk (1999) [2]_.

    References
    ----------
    .. [1] Cartwright, D. E., & Longuet-Higgins, M. S. (1956). The statistical
    distribution of the maxima of a random function. Proceedings of the Royal
    Society of London. Series A. Mathematical and Physical Sciences, 237(1209),
    212-232.
    .. [2] Liu, L., & Pezeshk, S. (1999). An improvement on the estimation of
        pseudoresponse spectral velocity using RVT method. Bulletin of the
        Seismological Society of America, 89(5), 1384-1389.

    """

    def compute_duration_rms(self, gm_duration, osc_freq, osc_damping,
                             m0, m1, m2, *args, **kwds):
        """Compute the oscillator duration used in the calculation of the
        root-mean-squared response.

        Based on Liu and Pezeshk (1999).

        Paramaters
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
        dur_ratio = (1 + 1. / (2 * np.pi * osc_damping)
                     * (foo / (1 + coef * foo ** power)))

        rms_duration = gm_duration * dur_ratio

        return rms_duration


def _load_bt12_data(region):
    """Load data from the Boore and Thompson (2012) parameter files.

    Input
    -----
    region : str
        Region for which the parameters were developed. Possible options are:
            'wna' - Western North America (active tectonic)
            'ena' - Eastern North America (stable tectonic)

    Returns
    -------
    np.recarray
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
    Longuet-Higgins (1956) [1]_ along with the root-mean-squared duration
    correction proposed by Boore and Thompson (2012) [2]_.

    References
    ----------
     .. [1] Cartwright, D. E., & Longuet-Higgins, M. S. (1956). The statistical
    distribution of the maxima of a random function. Proceedings of the Royal
    Society of London. Series A. Mathematical and Physical Sciences, 237(1209),
    212-232.
     .. [2] Boore, D. M., & Thompson, E. M. (2012). Empirical Improvements for
        Estimating Earthquake Response Spectra with Random‐Vibration Theory.
        Bulletin of the Seismological Society of America, 102(2), 761-772.

    """

    def __init__(self, region, mag, dist, **kwds):
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
        with Qhull [1]_, and on each triangle performing linear
        barycentric interpolation.

        References
        ----------
        .. [1] http://www.qhull.org/

        """
        super(BooreThompson2012, self).__init__(**kwds)
        region = get_region(region)
        self._CEOFS = _BT12_INTERPS[region](mag, np.log(dist))

    def compute_duration_rms(self, gm_duration, osc_freq, osc_damping, *args,
                             **kwds):
        """Compute the oscillator duration used in the calculation of the
        root-mean-squared response.

        Based on Boore and Joyner (1984).

        Paramaters
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
        dur_ratio = ((c1 + c2 * (1 - foo ** c3) / (1 + foo ** c3))
                     * (1 + c4 / (2 * np.pi * osc_damping) *
                        (foo / (1 + c5 * foo ** c6)) ** c7))

        return gm_duration * dur_ratio


def get_peak_calculator(method):
    """Select a peak calculator based on a string.

    """
    if method in ['DK83', 'DerKiureghian1983']:
        return DerKiureghian1983
    if method in ['TM87', 'ToroMcGuire1987']:
        return ToroMcGuire1987
    elif method in ['BJ84', 'BooreJoyner1984']:
        return BooreJoyner1984
    elif method in ['LP99', 'LiuPezeshk1999']:
        return LiuPezeshk1999
    elif method in ['BT12', 'BooreThompson2012']:
        return BooreThompson2012
    else:
        raise NotImplementedError


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
