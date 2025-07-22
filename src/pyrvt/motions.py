"""Random vibration theory (RVT) based motions.


The module attribute `DEFAULT_CALC` is used to control the default peak factor
calculator is one is not provided when a class is initialized.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import numpy.typing as npt
from scipy.constants import g as gravity
from scipy.interpolate import interp1d
from scipy.stats import linregress

from . import peak_calculators

DEFAULT_CALC = "V75"


def sort_increasing(*args: npt.ArrayLike) -> tuple[np.ndarray, ...]:
    """Sort arrays such that they are increasing.

    Check if the first array is is increasing, if not reverse the order. Same
    operation is applied to additional arrays.

    Parameters
    ----------
    *args : array_like
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


def calc_stress_drop(magnitude: float) -> float:
    """Stress drop using Atkinson & Boore (2011) model.

    Parameters
    ----------
    magnitude : float
         Moment magnitude of the stress drop.

    Returns
    -------
    stress_drop : float
        Stress drop (bars).

    """
    return 10 ** (3.45 - 0.2 * max(magnitude, 5.0))


def calc_geometric_spreading(
    dist: float, params: list[tuple[float, float | None]]
) -> float:
    """Geometric spreading defined by piece-wise linear model.

    Parameters
    ----------
    dist : float
        Closest distance to the rupture surface (km).
    params : List[(float,Optional[float])]
        List of (slope, limit) tuples that define the attenuation. For an
        infinite distance use `None`.  For example, [(1, `None`)] would provide
        for 1/R geometric spreading to an infinite distance.

    Returns
    -------
    coeff : float
        Geometric spreading coefficient.

    """
    initial = 1
    coeff = 1
    for slope, limit in params:
        # Compute the distance limited by the maximum distance of the slope.
        _dist = min(dist, limit) if limit else dist
        coeff *= (initial / _dist) ** slope
        if _dist < dist:
            initial = _dist
        else:
            break

    return coeff


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
    ) -> None:
        """Initialize the class."""
        self._freqs = freqs
        self._fourier_amps = fourier_amps
        self._duration = duration

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
    def fourier_amps(self) -> np.ndarray:
        """Acceleration Fourier amplitude values (g-sec)."""
        return self._fourier_amps

    @property
    def duration(self) -> float:
        """Duration of the ground motion for RVT analysis."""
        return self._duration

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

        # Need to perserve the site_tf for Wang and Rathje. It expects None
        if trans_func is None:
            trans_func = 1
            site_tf = None
        else:
            site_tf = trans_func = np.asarray(trans_func)

        resp = np.array(
            [
                self.calc_peak(
                    trans_func * calc_sdof_tf(self.freqs, of, osc_damping),
                    osc_freq=of,
                    osc_damping=osc_damping,
                    site_tf=site_tf,
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


class SourceTheoryMotion(RvtMotion):
    """Single-corner source theory model.

    The single-corner source theory model uses default parameters from Campbell (2003).
    """

    def __init__(
        self,
        magnitude: float,
        distance: float,
        region: str,
        stress_drop: float | None = None,
        depth: float | None = 8,
        peak_calculator: str | peak_calculators.Calculator | None = None,
        calc_kwds: dict | None = None,
        freqs: npt.ArrayLike | None = None,
        disable_site_amp: bool = False,
    ) -> None:
        """Initialize the motion.

        Parameters
        ----------
        magnitude : float
            Moment magnitude of the event.
        distance : float
            Epicentral distance (km).
        region : str
            Region for the parameters. Either 'cena' for Central and Eastern North
            America, or 'wna' for Western North America.
        stress_drop : float, optional
            Stress drop of the event (bars).  If `None`, then the default value
            is used. For `region` is 'cena', the default value is computed by
            the  model, while for `region` is 'wna' the
            default value is 100 bars.
        depth : float, optional
            Hypocenter depth (km). The `depth` is combined with the `distance`
            to compute the hypocentral distance.
        peak_calculator : `Calculator`, optional
            Peak calculator to use. If `None`, then the default peak
            calculator is used. The peak calculator may either be specified by
            a [pyrvt.peak_calculators.Calculator][] object, or by the
            initials of the calculator using
            [pyrvt.peak_calculators.get_peak_calculator][].
        calc_kwds : dict, optional
            Keywords to be passed during the creation the peak calculator.
            These keywords are only required for some peak calculators.
        freqs : array_like
            frequencies for which the Fourier amplitude spectrum should be computed.
            Defaults to `np.geomspace(0.05, 200, 512)`
        disable_site_amp: bool, optional
            if the crustal site amplification should be disable. Defaults to *False*.

        """
        super().__init__(peak_calculator=peak_calculator, calc_kwds=calc_kwds)

        self._disable_site_amp = disable_site_amp
        self.magnitude = magnitude
        self.distance = distance
        self.region = peak_calculators.get_region(region)

        if self.region == "wna":
            # Default parameters for the WUS from Campbell (2003)
            self.shear_velocity = 3.5
            self.path_atten_coeff = 180.0
            self.path_atten_power = 0.45
            self.density = 2.8
            self.site_atten = 0.04

            self.geometric_spreading = [(1, 40), (0.5, None)]

            if stress_drop:
                self.stress_drop = stress_drop
            else:
                self.stress_drop = 100.0

            # Crustal amplification from Campbell (2003) using the
            # log-frequency and the amplification based on a quarter-wave
            # length approximation
            self.site_amp = interp1d(
                np.log(
                    [
                        0.01,
                        0.09,
                        0.16,
                        0.51,
                        0.84,
                        1.25,
                        2.26,
                        3.17,
                        6.05,
                        16.60,
                        61.20,
                        100.00,
                    ]
                ),
                [
                    1.00,
                    1.10,
                    1.18,
                    1.42,
                    1.58,
                    1.74,
                    2.06,
                    2.25,
                    2.58,
                    3.13,
                    4.00,
                    4.40,
                ],
                bounds_error=False,
            )
        elif self.region == "cena":
            # Default parameters for the CEUS from Campbell (2003)
            self.shear_velocity = 3.6
            self.density = 2.8
            self.path_atten_coeff = 680.0
            self.path_atten_power = 0.36
            self.site_atten = 0.006

            self.geometric_spreading = [(1, 70), (0, 130), (0.5, None)]

            if stress_drop:
                self.stress_drop = stress_drop
            else:
                self.stress_drop = calc_stress_drop(magnitude)

            # Crustal amplification from Campbell (2003) using the
            # log-frequency and the amplification based on a quarter-wave
            # length approximation
            self.site_amp = interp1d(
                np.log(
                    [
                        0.01,
                        0.10,
                        0.20,
                        0.30,
                        0.50,
                        0.90,
                        1.25,
                        1.80,
                        3.00,
                        5.30,
                        8.00,
                        14.00,
                        30.00,
                        60.00,
                        100.00,
                    ]
                ),
                [
                    1.00,
                    1.02,
                    1.03,
                    1.05,
                    1.07,
                    1.09,
                    1.11,
                    1.12,
                    1.13,
                    1.14,
                    1.15,
                    1.15,
                    1.15,
                    1.15,
                    1.15,
                ],
                bounds_error=False,
                fill_value=(1.0, 1.15),
            )

        else:
            raise NotImplementedError

        # Depth to rupture
        self.depth = depth
        self.hypo_distance = np.sqrt(self.distance**2.0 + self.depth**2.0)

        # Constants
        self.seismic_moment = 10.0 ** (1.5 * (self.magnitude + 10.7))
        self.corner_freq = (
            4.9e6
            * self.shear_velocity
            * (self.stress_drop / self.seismic_moment) ** (1.0 / 3.0)
        )

        # Combine the three components and convert from displacement to acceleration
        self.calc_fourier_amps(freqs)
        self._duration = self.calc_duration()

    def calc_duration(self) -> float:
        """Compute the duration by combination of source and path.

        Returns
        -------
        duration : float
            Computed duration

        """
        # Source component
        duration_source = 1.0 / self.corner_freq

        # Path component
        if self.region == "wna":
            duration_path = 0.05 * self.hypo_distance
        elif self.region == "cena":
            duration_path = 0.0
            if self.hypo_distance > 10:
                # 10 < R <= 70 km
                duration_path += 0.16 * (min(self.hypo_distance, 70) - 10.0)
            if self.hypo_distance > 70:
                # 70 < R <= 130 km
                duration_path += -0.03 * (min(self.hypo_distance, 130) - 70.0)
            if self.hypo_distance > 130:
                # 130 km < R
                duration_path += 0.04 * (self.hypo_distance - 130.0)
        else:
            raise NotImplementedError

        return duration_source + duration_path

    def calc_fourier_amps(self, freqs: npt.ArrayLike | None = None) -> np.ndarray:
        """Compute the acceleration Fourier amplitudes for a frequency range.

        Parameters
        ----------
        freqs : array_like, optional
            Frequency range. If no frequency range is specified then
            :func:`log_spaced_values(0.05, 200.)` is used.

        Returns
        -------
        fourier_amps : :class:`np.ndarray`
            acceleration Fourier amplitudes

        """
        if freqs is None:
            self._freqs = log_spaced_values(0.05, 200.0)
        else:
            (self._freqs,) = sort_increasing(np.asarray(freqs))

        self._duration = self.calc_duration()

        # Model component
        const = (0.55 * 2.0) / (
            np.sqrt(2.0) * 4.0 * np.pi * self.density * self.shear_velocity**3.0
        )
        source_comp = (
            const
            * self.seismic_moment
            / (1.0 + (self._freqs / self.corner_freq) ** 2.0)
        )

        # Path component
        path_atten = self.path_atten_coeff * self._freqs**self.path_atten_power
        geo_atten = calc_geometric_spreading(
            self.hypo_distance, self.geometric_spreading
        )

        path_comp = geo_atten * np.exp(
            (-np.pi * self._freqs * self.hypo_distance)
            / (path_atten * self.shear_velocity)
        )

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

        site_comp = 1 if self._disable_site_amp else (site_amp * site_dim)

        # Conversion factor to convert from dyne-cm into gravity-sec
        conv = 1.0e-20 / (100 * gravity)
        # Combine the three components and convert from displacement to
        # acceleration

        self._fourier_amps = (
            conv
            * (2.0 * np.pi * self._freqs) ** 2.0
            * source_comp
            * path_comp
            * site_comp
        )


class StaffordEtAl22Motion(RvtMotion):
    # Site amplification from Al Atik & Abrahamson (2021) generic rock amplification
    # function for Vs30 of 760 m/s obtained from inversion of the Chiou & Youngs (2014)
    # response spectral model note that the ordinate at (0.01, 1.0) has been added to
    # the amp function that was provided. These values were taken from:
    # https://github.com/pstafford/StochasticGroundMotionSimulation.jl/blob/master/src/fourier/PJSsite.jl
    #
    # In the original code, the interpolation is done on log(amplitdue) and linear
    # frequency. I would typically do this on log frequency and linear amplitude.
    _ln_site_amp_interpolator = None

    def __init__(
        self,
        mag: float,
        dist_rup: float | None = None,
        dist_jb: float | None = None,
        mechanism: str = "U",
        method: str = "continuous",
        delta_ztor: float = 0,
        freqs: npt.ArrayLike | None = None,
        disable_site_amp: bool = False,
    ) -> None:
        """Point source model developed by Stafford (2021) for a Vs30 of 760 m/s.

        Use an RVT framework, and assume the following duration/peak factor models:
            - Boore & Thompson (2014) for excitation duration
            - Boore & Thompson (2015) for RMS duration
            - Vanmarke (1975)/Der Kiureghian (1980) peak factor expression;
              per Boore & Thompson (2015)

        Note that other peak factors models are not permitted to be consistent with the
        Sea22 model.

        Parameters
        ----------
        mag : float
            moment magnitude of the event.
        dist_rup : float, optional
            closest distance to the rupture (km). Either *dist_rup* or *dist_jb* must be
            provided.
        dist_jb : float, optional
            Joyner-Boore distance (km). Either *dist_rup* or *dist_jb* must be
            provided.
        mechanism : str, optional
            earthquake mechansim. Options are: "U", "SS", "NS", "RS".  Defaults to 'U'.
        method: str, optional
            geometric spreading model. Options are: "continuous" or "trilinear".
            Defaults to "continuous".
        delta_ztor : float, optional
            difference in the top of the rutpure (km)
        freqs : array_like
            frequencies for which the Fourier amplitude spectrum should be computed.
            Defaults to `np.geomspace(0.05, 200, 512)`
        disable_site_amp: bool, optional
            if the crustal site amplification should be disable. Defaults to *False*.
        """

        if dist_rup is None and dist_jb is None:
            raise NotImplementedError
        elif dist_rup is None:
            # Compute rupture distance for dist_jb
            depth_tor = StaffordEtAl22Motion.calc_depth_tor(mag, mechanism)
            dist_rup = np.sqrt(depth_tor**2 + dist_jb**2)

        super().__init__(
            peak_calculator="BT15",
            calc_kwds={"mag": mag, "dist": dist_rup, "region": "wus"},
        )

        if freqs is None:
            self._freqs = np.geomspace(0.05, 200, 512)
        else:
            self._freqs = np.asarray(freqs)

        # Constants
        shear_vel = 3.5
        density = 2.75
        site_atten = 0.039

        # Parameters
        if method == "continuous":
            # Stress parameter components
            ln_ds0 = 4.599
            dln_ds0 = 0.4624
            dz_a = 0.0453
            dz_b = 0.109
            # Geometric spreading
            y_1 = 1.1611
            y_f = 0.5
            r_t = 50
            # Finite fault
            h_a = -0.8712
            h_b = 0.4451
            h_c = 1.1513
            h_d = 5.0948
            h_e = 7.2725
            # Anelastic attenuation
            Q_0 = 205.4
            n_a = 0.6884
            # Magnitude scaling of eta
            n_b = 0.1354
            n_c = 5.1278
        elif method == "trilinear":
            # Stress parameter components
            ln_ds0 = 5.07
            dln_ds0 = 0.6451
            dz_a = 0.4077
            dz_b = 0.117
            # Geometric spreading
            y_1 = 1.1680
            y_2 = 0.9293
            y_f = 0.5
            # Finite fault
            h_a = -0.7771
            h_b = 0.4768
            h_c = 1.1513
            h_d = 3.418
            h_e = 7.088
            # Anelastic attenuation
            Q_0 = 183.7
            n = 0.7077
        else:
            raise NotImplementedError

        # Stress drop in bars
        stress_drop = np.exp(
            ln_ds0
            + dln_ds0 * np.minimum(mag - 5, 0)
            + delta_ztor * (dz_a + dz_b / np.cosh(2 * np.maximum(mag - 4.5, 0)))
        )

        # Source spectrum
        const = (0.55 / np.sqrt(2) * 2) / (4 * np.pi * density * shear_vel**3) * 1e-20
        seismic_moment = 10 ** (1.5 * (mag + 10.7))

        corner_freq = 4.9058e6 * shear_vel * (stress_drop / seismic_moment) ** (1 / 3)
        source_comp = (const * seismic_moment) / (1 + (self._freqs / corner_freq) ** 2)

        # Path scaling

        # Finite fault factor h(m)
        fault_fact = np.exp(
            h_a
            + h_b * mag
            + ((h_b - h_c) / h_d) * np.log(1 + np.exp(-h_d * (mag - h_e)))
        )
        dist_ps = dist_rup + fault_fact
        # Geometric spreading term
        if method == "continuous":
            geom_spread = np.exp(
                -y_1 * np.log(dist_ps)
                + (y_1 - y_f) / 2 * np.log((dist_rup**2 + r_t**2) / (1**2 + r_t**2))
            )
            n = n_a + n_b * np.tanh(mag - n_c)
            dist_ae = dist_rup
        elif method == "trilinear":
            geom_spread = calc_geometric_spreading(
                dist_ps, [(y_1, 25), (y_2, 85), (y_f, None)]
            )
            dist_ae = dist_ps
        else:
            raise NotImplementedError

        # Distance metric is different between the two forms
        anelastic_atten = np.exp(
            -(np.pi * self._freqs ** (1 - n) * dist_ae) / (Q_0 * shear_vel)
        )

        path_comp = geom_spread * anelastic_atten

        # Convert to acceleration (cm/s) and then into g-se
        conv = (2 * np.pi * self._freqs) ** 2 / (gravity * 100)

        if disable_site_amp:
            site_tf = 1.0
        else:
            site_tf = StaffordEtAl22Motion.site_amp(self._freqs, site_atten)

        # Combine the three components and convert from displacement to acceleration
        self._fourier_amps = conv * source_comp * path_comp * site_tf
        self._dist_ps = dist_ps

        self._duration = StaffordEtAl22Motion.calc_duration(corner_freq, dist_ps)

    @classmethod
    def site_amp(cls, freqs: npt.ArrayLike, site_atten: float) -> np.ndarray:
        if cls._ln_site_amp_interpolator is None:
            # Load the data
            data = np.genfromtxt(
                gzip.open(Path(__file__).parent / "data" / "sea22-site_amp.csv.gz"),
                delimiter=",",
                names=True,
                skip_header=1,
            ).view(np.recarray)

            _ln_site_amp = np.log(data["site_amp"])
            cls._ln_site_amp_interpolator = interp1d(
                data["freq"],
                _ln_site_amp,
                kind="linear",
                bounds_error=False,
                fill_value=(_ln_site_amp[0], _ln_site_amp[-1]),
            )

        return np.exp(cls._ln_site_amp_interpolator(freqs)) * np.exp(
            -np.pi * site_atten * freqs
        )

    @property
    def dist_ps(self) -> float:
        """Equivalent point source distance (km)."""
        return self._dist_ps

    @staticmethod
    def calc_depth_tor(mag: float, mechanism: str) -> float:
        """Top of rupture model from Chiou and Youngs (2014)

        Parameters
        ----------
        mag : float
            moment magnitude of the event (:math:`M_w`)
        mechanism : str
            fault mechanism. Valid options: "U", "SS", "NS",
            "RS".

        Returns
        -------
        depth_tor : float
            estimated depth to top of rupture (km)

        """
        if mechanism == "RS":
            # Reverse and reverse-oblique faulting
            fact = 2.704 - 1.226 * max(mag - 5.849, 0)
        else:
            # Combined strike-slip and normal faulting
            fact = 2.673 - 1.136 * max(mag - 4.970, 0)

        return max(fact, 0) ** 2

    @staticmethod
    def calc_duration(corner_freq: float, dist_ps: float) -> float:
        """Boore & Thomspson (2014) duration model."""

        # Source component. Equation 2
        d_s = 1.0 / corner_freq

        # Path component. Table 1
        # Maximum distance set to be 10 km
        DISTS = [0, 7, 45, 125, 175, 270]
        D_P = [0.0, 2.4, 8.4, 10.9, 17.4, 34.2]

        if dist_ps < DISTS[-1]:
            d_p = np.interp(dist_ps, DISTS, D_P)
        else:
            d_p = D_P[-1] + 0.156 * (dist_ps - DISTS[-1])

        return d_s + d_p


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
        duration: float | None = None,
        osc_damping: float | None = 0.05,
        event_kwds: dict | None = None,
        window_len: int | None = None,
        peak_calculator: str | peak_calculators.Calculator | None = None,
        calc_kwds: dict | None = None,
    ) -> None:
        """Initialize the motion.

        Parameters
        ----------
        osc_freqs : array_like
            Frequencies of the oscillator response (Hz).
        osc_accels_target : array_like
            Spectral acceleration of the oscillator at the specified
            frequencies (g).
        duration : float, optional
            Duration of the ground motion (sec). If `None`, then the duration
            is computed using the `event_kwds`.
        osc_damping : float, optional
            Fractional damping of the oscillator (dec). Default value is 0.05
            for a damping ratio of 5%.
        event_kwds : Dict, optional
            Keywords passed to :class:`~.motions.SourceTheoryMotion` and used
            to compute the duration of the motion. Either `duration` or
            `event_kwds` should be specified.
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

        if duration:
            self._duration = duration
        else:
            stm = SourceTheoryMotion(**event_kwds)
            self._duration = stm.calc_duration()

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
