import numpy as np
import pytest
from numpy.testing import assert_allclose

import pyrvt


@pytest.fixture
def motion():
    """Synthetic single-corner motion in g-sec units."""
    freqs = pyrvt.motions.log_spaced_values(0.05, 100.0, per_decade=64)
    fourier_amps = (
        0.05
        * (freqs**2 / (1 + (freqs / 0.5) ** 2))
        * np.exp(-np.pi * 0.035 * freqs)
        / freqs**2
    )
    return pyrvt.motions.RvtMotion(
        freqs=freqs, fourier_amps=fourier_amps, duration=8.0, peak_calculator="V75"
    )


def test_peak_params_invariants(motion):
    """Properties, caching, transfer-function plumbing, and scaling."""
    assert_allclose(motion.angular_freqs, 2 * np.pi * motion.freqs, rtol=0)

    # Each cached property agrees with the method that fills it
    assert motion.pga == motion.calc_pga()
    assert motion.pgv == motion.calc_pgv()
    assert motion.arias_intensity == motion.calc_arias_intensity()
    assert motion.cav == motion.calc_cav()

    # `tf` is the public keyword; this pins the signature
    ones = np.ones_like(motion.freqs)
    assert motion.calc_pga(tf=ones) == motion.pga
    assert motion.calc_pgv(tf=ones) == motion.pgv
    assert motion.calc_arias_intensity(tf=ones) == motion.arias_intensity

    # The peak scales linearly with the transfer function, Arias quadratically
    assert motion.calc_pga(2 * ones) == pytest.approx(2 * motion.pga, rel=1e-12)
    assert motion.calc_pgv(2 * ones) == pytest.approx(2 * motion.pgv, rel=1e-12)
    assert motion.calc_arias_intensity(2 * ones) == pytest.approx(
        4 * motion.arias_intensity, rel=1e-12
    )

    # CAV is an empirical function of Arias intensity and duration
    assert motion.cav == pytest.approx(
        np.exp(
            1.553
            + 0.496 * np.log(motion.arias_intensity)
            + 0.356 * np.log(motion.duration)
        ),
        rel=1e-12,
    )


def test_peak_params_regression(motion):
    """Locked values, guarding against silent numerical drift."""
    assert motion.pga == pytest.approx(0.031619466772979048, rel=1e-10)
    assert motion.pgv == pytest.approx(27.562834894066476, rel=1e-10)
    assert motion.arias_intensity == pytest.approx(0.0244873720524409, rel=1e-10)
    assert motion.cav == pytest.approx(1.5735287667201079, rel=1e-10)


def test_compatible_rvt_motion():
    """CompatibleRvtMotion recovers its target spectrum within tolerance."""
    pygmm = pytest.importorskip("pygmm")

    src = pygmm.fourier_spectrum.SourceTheoryModel(
        magnitude=6.0, distance=20.0, region="wna"
    )
    motion = pyrvt.motions.RvtMotion.from_fas(
        src, peak_calculator=pyrvt.peak_calculators.DerKiureghian1985()
    )

    osc_freqs = np.logspace(-1, 2, num=50)
    osc_accels_target = motion.calc_osc_accels(osc_freqs, 0.05)

    compat = pyrvt.motions.CompatibleRvtMotion(
        osc_freqs,
        osc_accels_target,
        duration=src.duration,
        osc_damping=0.05,
        peak_calculator=pyrvt.peak_calculators.DerKiureghian1985(),
    )

    osc_accels_compat = compat.calc_osc_accels(osc_freqs, 0.05)

    # Might be off by a few percent because of difficulties with the inversion.
    assert_allclose(osc_accels_target, osc_accels_compat, rtol=0.03, atol=0.05)


def test_arias_intensity_matches_time_domain():
    """Arias intensity from the one-sided FAS equals the time-domain integral."""
    rng = np.random.default_rng(0)
    time_step, count = 0.01, 8192
    times = time_step * np.arange(count)
    accels = 0.1 * rng.standard_normal(count) * np.exp(-(((times - 30) / 10) ** 2))

    freqs = np.fft.rfftfreq(count, time_step)[1:]
    fourier_amps = np.abs(np.fft.rfft(accels))[1:] * time_step
    motion = pyrvt.motions.RvtMotion(freqs, fourier_amps, duration=60.0)

    gravity = 9.80665
    expected = np.pi / (2 * gravity) * np.sum((gravity * accels) ** 2) * time_step
    assert motion.calc_arias_intensity() == pytest.approx(expected, rel=1e-3)
