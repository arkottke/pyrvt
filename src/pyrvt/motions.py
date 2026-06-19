"""Random vibration theory (RVT) based motions.


The module attribute `DEFAULT_CALC` is used to control the default peak factor
calculator is one is not provided when a class is initialized.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import numpy.typing as npt
from scipy.constants import g as gravity
from scipy.interpolate import interp1d
from scipy.stats import linregress

from . import peak_calculators

DEFAULT_CALC = "V75"


def sort_increasing(*args):
    """Sort arrays such that they are increasing.

    Check if the first array is is increasing, if not reverse the order. Same
    operation is applied to additional arrays.

    Args:
    args : array_like
        arrays to be re-ordered.

    Returns
    -------
    tuple
        tuple containing sorted :class:`numpy.ndarray`'s.

    Raises
    ------
    :class:`NotImplementedError`
        If first array is not monotonic.

    """
    diffs = np.diff(args[0])
    if np.all(diffs >= 0):
        # All increasing, do nothing
        pass
    elif np.all(diffs <= 0):
        # All decreasing, reverse
        args = [a[::-1] for a in args]
    else:
        raise NotImplementedError("Values are not regularly ordered.")

    return args


def log_spaced_values(lower: float, upper: float, per_decade: int = 512) -> np.ndarray:
    """Generate values with constant log-spacing.

    Parameters
    ----------
    lower : float
        lower end of the range.
    upper : float
        upper end of the range.
    per_decade : int, optional
        number of points per decade. Default is 512 points per decade.

    Returns
    -------
    values : :class:`numpy.ndarray`
        Log-spaced values.

    """
    lower = np.log10(lower)
    upper = np.log10(upper)
    count = int(np.ceil(per_decade * (upper - lower)))
    return np.logspace(lower, upper, num=count)


def calc_sdof_tf(
    freqs: npt.ArrayLike, osc_freq: float, osc_damping: float
) -> np.ndarray:
    """Single-degree-of-freedom transfer function.

    When applied on the acceleration Fourier amplitude spectrum, it provides
    the pseudo-spectral acceleration.

    Parameters
    ----------
    freqs : array_like
        Frequencies at which the transfer function should be calculated  (Hz).
    osc_freq : float
        Frequency of the oscillator (Hz).
    osc_damping : float
        Fractional damping of the oscillator (decimal).

    Returns
    -------
    :class:`numpy.ndarray`
        Complex valued transfer function.

    """
    freqs = np.asarray(freqs)
    return -(osc_freq**2.0) / (
        freqs**2 - osc_freq**2 - 2.0j * osc_damping * osc_freq * freqs
    )


class RvtMotion:
    """Random vibration theory motion.

    Parameters
    ----------
    freqs : array_like, optional
        Frequency array (Hz).
    fourier_amps : array_like, optional
        Absolute value of acceleration Fourier amplitudes.
    duration : float, optional
        Ground motion duration (sec).
    peak_calculator : `Calculator`, optional
        Peak calculator to use. If `None`, then the default peak
        calculator is used.  The peak calculator may either be specified
        by a [pyrvt.peak_calculators.Calculator][] instance, or created by the
        abbreviation of the calculator using
        [pyrvt.peak_calculators.get_peak_calculator][].
    calc_kwds : dict, optional
        Keywords to be passed during the creation the peak calculator.
        These keywords are only required for some peak calculators.
    """

    def __init__(
        self,
        freqs: npt.ArrayLike | None = None,
        fourier_amps: npt.ArrayLike | None = None,
        duration: float | None = None,
        peak_calculator: str | peak_calculators.Calculator | None = None,
        calc_kwds: dict | None = None,
    ):
        """Initialize the class."""
        self._freqs = freqs
        self._fourier_amps = fourier_amps
        self._duration = duration

        self._pga = None
        self._pgv = None
        self._arias_intensity = None
        self._cav = None

        if self._freqs is not None:
            self._freqs, self._fourier_amps = sort_increasing(
                self._freqs, self._fourier_amps
            )

        if isinstance(peak_calculator, peak_calculators.Calculator):
            self.peak_calculator = peak_calculator
        else:
            self.peak_calculator = peak_calculators.get_peak_calculator(
                peak_calculator or DEFAULT_CALC, calc_kwds
            )

    @property
    def freqs(self) -> np.ndarray:
        """Frequency values (Hz)."""
        return self._freqs

    @property
    def angular_freqs(self) -> np.ndarray:
        """Angular frequency values (rad/sec)."""
        return 2 * np.pi * self._freqs

    @property
    def fourier_amps(self) -> np.ndarray:
        """Acceleration Fourier amplitude values (g-sec)."""
        return self._fourier_amps

    @property
    def duration(self) -> float:
        """Duration of the ground motion for RVT analysis."""
        return self._duration

    @property
    def pga(self) -> float:
        """Peak ground acceleration (g)."""
        if self._pga is None:
            self._pga = self.calc_pga()
        return self._pga

    @property
    def pgv(self) -> float:
        """Peak ground velocity (cm/sec)."""
        if self._pgv is None:
            self._pgv = self.calc_pgv()
        return self._pgv

    @property
    def arias_intensity(self) -> float:
        """Arias intensity (m/s)."""
        if self._arias_intensity is None:
            self._arias_intensity = self.calc_arias_intensity()
        return self._arias_intensity

    @property
    def cav(self) -> float:
        """Cumulative absolute velocity (m/s)."""
        if self._cav is None:
            self._cav = self.calc_cav()
        return self._cav

    @classmethod
    def from_fas(
        cls,
        fas,
        peak_calculator: "str | peak_calculators.Calculator | None" = None,
        calc_kwds: dict | None = None,
    ) -> "RvtMotion":
        """Build an :class:`RvtMotion` from any object exposing a Fourier-spectrum shape.

        Parameters
        ----------
        fas : object
            Any object exposing ``freqs`` [Hz], ``fourier_amps`` [g-sec], and
            ``duration`` [sec] attributes (e.g. an instance of a
            ``pygmm.fourier_spectrum`` model, or a
            :class:`pygmm.contracts.FourierSpectrum` dataclass).
        peak_calculator, calc_kwds
            Forwarded to :class:`RvtMotion`.
        """
        return cls(
            freqs=np.asarray(fas.freqs),
            fourier_amps=np.asarray(fas.fourier_amps),
            duration=float(fas.duration),
            peak_calculator=peak_calculator,
            calc_kwds=calc_kwds,
        )

    def calc_pga(self, transfer_func: npt.ArrayLike | None = None) -> float:
        """Peak ground acceleration [g] via RVT.

        Parameters
        ----------
        transfer_func : array_like, optional
            Additional transfer function applied to the acceleration FAS prior
            to the peak calculation.
        """
        return self.calc_peak(transfer_func)

    def calc_pgv(self, transfer_func: npt.ArrayLike | None = None) -> float:
        """Peak ground velocity [cm/sec] via RVT.

        Computed by integrating the acceleration FAS in the frequency domain
        (multiplication by :math:`1 / (i\\omega)`) and then applying the peak
        calculator. The result is scaled from g-sec to cm/sec.
        """
        omega = self.angular_freqs
        mask = ~np.isclose(omega, 0)
        tf_av = np.zeros_like(omega, dtype=complex)
        tf_av[mask] = 1 / (omega[mask] * 1j)
        if transfer_func is not None:
            tf_av = tf_av * np.asarray(transfer_func)
        # g-sec * (1/rad-sec) -> g-sec * sec = g; multiply by gravity (m/s^2)
        # then by 100 to convert m/s -> cm/s.
        return gravity * 100 * self.calc_peak(tf_av)

    def calc_arias_intensity(
        self, transfer_func: npt.ArrayLike | None = None
    ) -> float:
        """Compute the Arias intensity.

        Parameters
        ----------
        transfer_func : array_like, optional
            Transfer function to apply to the motion. If ``None``, no
            transfer function is applied.

        Returns
        -------
        arias_intensity : float
            Arias intensity (m/s).

        """
        tf = 1 if transfer_func is None else np.abs(np.asarray(transfer_func))
        fa = tf * self._fourier_amps
        m0 = np.trapezoid(fa**2, self._freqs)
        return np.pi * gravity / 2 * m0

    def calc_cav(self, transfer_func: npt.ArrayLike | None = None) -> float:
        """Compute the cumulative absolute velocity (CAV).

        Uses an empirical regression on Arias intensity and duration based on
        observed ground motions.

        Parameters
        ----------
        transfer_func : array_like, optional
            Transfer function to apply to the motion. If ``None``, no
            transfer function is applied.

        Returns
        -------
        cav : float
            Cumulative absolute velocity (m/s).

        """
        return np.exp(
            1.553
            + 0.496 * np.log(self.calc_arias_intensity(transfer_func))
            + 0.356 * np.log(self.duration)
        )

    def calc_osc_accels(
        self,
        osc_freqs: npt.ArrayLike,
        osc_damping: float = 0.05,
        trans_func: npt.ArrayLike | None = None,
    ) -> np.ndarray:
        """Pseudo-acceleration spectral response of an oscillator.

        Parameters
        ----------
        osc_freqs : float
            Frequency of the oscillator (Hz).
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for a
            damping ratio of 5%.
        trans_func : array_like, optional
            Transfer function to be applied to motion prior calculation of the
            oscillator response.

        Returns
        -------
        spec_accels : `numpy.ndarray`
            Peak pseudo-spectral acceleration of the oscillator

        """
        has_trans_func = trans_func is not None
        trans_func = 1 if trans_func is None else np.asarray(trans_func)

        resp = np.array(
            [
                self.calc_peak(
                    trans_func * calc_sdof_tf(self.freqs, of, osc_damping),
                    osc_freq=of,
                    osc_damping=osc_damping,
                    site_tf=trans_func if has_trans_func else None,
                )
                for of in osc_freqs
            ]
        )

        return resp

    def calc_peak(self, transfer_func: npt.ArrayLike | None = None, **kwds) -> float:
        """Compute the peak response.

        Parameters
        ----------
        transfer_func : array_like, optional
            Transfer function to apply to the motion. If `None`, then no
            transfer function is applied.

        Returns
        -------
        peak : float
            Calculated peak

        """
        if transfer_func is None:
            fourier_amps = self._fourier_amps
        else:
            fourier_amps = np.abs(transfer_func) * self._fourier_amps

        return self.peak_calculator(self._duration, self._freqs, fourier_amps, **kwds)[
            0
        ]

    def calc_attenuation(
        self, min_freq: float, max_freq: float | None = None
    ) -> tuple[float, float, np.ndarray, np.ndarray]:
        r"""Compute the site attenuation (κ) based on a log-linear fit.

        Parameters
        ----------
        min_freq : float
            minimum frequency of the fit (Hz).
        max_freq : float, optional
            maximum frequency of the fit. If `None`, then the maximum frequency range is
            used.

        Returns
        -------
        atten : float
            attenuation parameter.
        r_sqr : float
            squared correlation coefficient of the fit (R²). See
            `scipy.stats.linregress`.
        freqs : array_like
            selected frequencies
        fitted : array_like
            fitted values

        Notes
        -----
        This function computes the site attenuation defined by Anderson & Hough (1984)
        [@anderson84] as:

        $$
        a(f) = A_0 \exp(-\pi \kappa f) \text( for ) f > f_E
        $$

        for a single Fourier amplitude spectrum

        """
        max_freq = max_freq or self.freqs[-1]
        mask = (min_freq <= self.freqs) & (self.freqs <= max_freq)

        slope, intercept, r_value, p_value, stderr = linregress(
            self.freqs[mask], np.log(self.fourier_amps[mask])
        )

        atten = slope / -np.pi
        freqs = self.freqs[mask]
        fitted = np.exp(intercept + slope * freqs)
        return atten, r_value**2, freqs, fitted



class CompatibleRvtMotion(RvtMotion):
    """Response spectrum compatible RVT motion.

    A [`CompatibleRvtMotion`][pyrvt.motions.CompatibleRvtMotion] object is used to
    compute a Fourier amplitude spectrum that is compatible with a target response
    spectrum.

    """

    def __init__(
        self,
        osc_freqs: npt.ArrayLike,
        osc_accels_target: npt.ArrayLike,
        duration: float,
        osc_damping: float = 0.05,
        window_len: int | None = None,
        peak_calculator: str | peak_calculators.Calculator | None = None,
        calc_kwds: dict | None = None,
    ):
        """Initialize the motion.

        Parameters
        ----------
        osc_freqs : array_like
            Frequencies of the oscillator response (Hz).
        osc_accels_target : array_like
            Spectral acceleration of the oscillator at the specified
            frequencies (g).
        duration : float
            Duration of the ground motion (sec).
        osc_damping : float, optional
            Fractional damping of the oscillator (dec). Default value is 0.05
            for a damping ratio of 5%.
        window_len : int, optional
            Window length used for smoothing the computed Fourier amplitude
            spectrum. If `None`, then no smoothing is applied. The smoothing
            is applied as a moving average with a width of `window_len`.
        peak_calculator : `Calculator`, optional
            Peak calculator to use. If `None`, then the default peak
            calculator is used. The peak calculator may either be specified by
            a [pyrvt.peak_calculators.Calculator][] object, or by the
            initials of the calculator using
            [pyrvt.peak_calculators.get_peak_calculator][].
        calc_kwds : dict, optional
            Keywords to be passed during the creation the peak calculator.
            These keywords are only required for some peak calculators.

        """
        super().__init__(peak_calculator=peak_calculator, calc_kwds=calc_kwds)

        osc_freqs, osc_accels_target = sort_increasing(
            np.asarray(osc_freqs), np.asarray(osc_accels_target)
        )

        self._duration = duration

        fourier_amps = self._estimate_fourier_amps(
            osc_freqs, osc_accels_target, osc_damping
        )

        # The frequency needs to be extended to account for the fact that the
        # oscillator transfer function has a width. The number of frequencies
        # depends on the range of frequencies provided.
        self._freqs = log_spaced_values(osc_freqs[0] / 2.0, 2.0 * osc_freqs[-1])
        self._fourier_amps = np.empty_like(self._freqs)

        # Indices of the first and last point with the range of the provided
        # response spectra
        indices = np.argwhere(
            (osc_freqs[0] < self._freqs) & (self._freqs < osc_freqs[-1])
        )
        first = indices[0, 0]
        # last is extend one past the usable range to allow use of first:last
        # notation
        last = indices[-1, 0] + 1
        log_freqs = np.log(self._freqs)
        log_osc_freqs = np.log(osc_freqs)

        self._fourier_amps[first:last] = np.exp(
            np.interp(log_freqs[first:last], log_osc_freqs, np.log(fourier_amps))
        )

        def extrapolate():
            """Extrapolate the first and last value of FAS."""

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
                self._freqs[first : first + 2],
                self._fourier_amps[first : first + 2],
                None,
            )
            # Update the last point using the third- and second-to-last points
            self._fourier_amps[last:] = _extrap(
                self._freqs[last:],
                self._freqs[last - 2 : last],
                self._fourier_amps[last - 2 : last],
                None,
            )

        extrapolate()

        # Apply a ratio correction between the computed at target response
        # spectra
        self.iterations = 0
        self.rmse = 1.0

        max_iterations = 30
        tolerance = 5e-6

        osc_accels = self.calc_osc_accels(osc_freqs, osc_damping)

        # Smoothing operator
        if window_len:
            window = np.ones(window_len, "d")
            window /= window.sum()

        while self.iterations < max_iterations and tolerance < self.rmse:
            # Correct the FAS by the ratio of the target to computed
            # oscillator response. The ratio is applied over the same
            # frequency range. The first and last points in the FAS are
            # determined through extrapolation.

            self._fourier_amps[first:last] *= np.exp(
                np.interp(
                    log_freqs[first:last],
                    log_osc_freqs,
                    np.log(osc_accels_target / osc_accels),
                )
            )

            extrapolate()

            # Apply a running average to smooth the signal
            if window_len:
                self._fourier_amps = np.convolve(window, self._fourier_amps, "same")

            # Recompute the response spectrum
            osc_accels = self.calc_osc_accels(osc_freqs, osc_damping)

            # Compute the fit between the target and computed oscillator
            # response
            self.rmse = np.sqrt(np.mean((osc_accels_target - osc_accels) ** 2))

            self.iterations += 1

    @classmethod
    def from_response_spectrum(
        cls,
        rs,
        duration: float,
        **kw,
    ) -> "CompatibleRvtMotion":
        """Create from any object with .periods, .spec_accels, .damping.

        Duck-typed: accepts ``pygmm.contracts.ResponseSpectrum`` or any object
        with the same attributes.

        Parameters
        ----------
        rs :
            Response spectrum with ``.periods`` [s], ``.spec_accels`` [g],
            and ``.damping`` [decimal].
        duration : float
            Ground-motion duration [s].
        **kw
            Forwarded to :class:`CompatibleRvtMotion` (e.g. ``peak_calculator``).
        """
        return cls(
            osc_freqs=1.0 / np.asarray(rs.periods),
            osc_accels_target=np.asarray(rs.spec_accels),
            duration=duration,
            osc_damping=rs.damping,
            **kw,
        )

    def _estimate_fourier_amps(
        self, osc_freqs: npt.ArrayLike, osc_accels: npt.ArrayLike, osc_damping: float
    ) -> np.ndarray:
        """Estimate the Fourier amplitudes.

        Compute an estimate of the FAS using the Gasparini & Vanmarcke (1976)
        methodology.  The response is first computed at the lowest frequency and then
        subsequently computed at higher frequencies.

        Parameters
        ----------
        osc_freqs : array_like
            Oscillator frequencies in increasing order (Hz).
        osc_accels : array_like
            Psuedo-spectral accelerations of the oscillator (g).
        osc_damping : float
            Fractional damping of the oscillator (dec). For example, 0.05 for a
            damping ratio of 5%.

        Returns
        -------
        :class:`numpy.ndarray`
            acceleration Fourier amplitude values at the specified frequencies
            specifed by `osc_freqs`.

        """
        # Compute initial value using Vanmarcke methodology.
        peak_factor = 2.5
        fa_sqr_prev = 0.0
        total = 0.0
        sdof_factor = np.pi / (4.0 * osc_damping) - 1.0
        fourier_amps = np.empty_like(osc_freqs)
        for i, (osc_freq, osc_accel) in enumerate(zip(osc_freqs, osc_accels)):
            # TODO: simplify equation and remove duration
            fa_sqr_cur = (
                (self.duration * osc_accel**2) / (2 * peak_factor**2) - total
            ) / (osc_freq * sdof_factor)

            if fa_sqr_cur < 0:
                fourier_amps[i] = fourier_amps[i - 1]
                fa_sqr_cur = fourier_amps[i] ** 2
            else:
                fourier_amps[i] = np.sqrt(fa_sqr_cur)

            if i == 0:
                total = fa_sqr_cur * osc_freq / 2.0
            else:
                total += (fa_sqr_cur - fa_sqr_prev) / 2 * (osc_freq - osc_freqs[i - 1])
        return fourier_amps
