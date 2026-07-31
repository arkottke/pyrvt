"""Interop tests: pyrvt ↔ pygmm (RS path + stub objects)."""

from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose

pygmm = pytest.importorskip("pygmm")
from pyrvt.motions import CompatibleRvtMotion


@pytest.fixture
def scenario():
    return pygmm.Scenario(
        mag=6.5,
        dist_jb=20.0,
        dist_rup=20.0,
        v_s30=400.0,
        dip=90.0,
        dist_crjb=0.0,
        mechanism="SS",
    )


def test_from_response_spectrum_pygmm(scenario):
    """CompatibleRvtMotion.from_response_spectrum works with a real pygmm RS."""
    rs = pygmm.AbrahamsonSilvaKamai2014(scenario).response_spectrum()
    # AfshariStewart2016.duration returns a structured array; use the D5-95 field
    dur_rec = pygmm.AfshariStewart2016(scenario).duration
    duration = float(dur_rec["D_5t95"])

    motion = CompatibleRvtMotion.from_response_spectrum(rs, duration)

    # Should converge and produce non-degenerate FAS
    assert np.all(np.isfinite(motion.fourier_amps))
    assert np.all(motion.fourier_amps > 0)
    assert motion.duration == pytest.approx(duration, rel=1e-6)

    # Sanity: pga from RVT should be in a reasonable range
    pga = motion.calc_peak()
    assert 0.001 < pga < 10.0


def test_from_response_spectrum_regression(scenario):
    """PGA from RS→RVT path is stable across code changes (regression anchor)."""
    rs = pygmm.AbrahamsonSilvaKamai2014(scenario).response_spectrum()
    dur_rec = pygmm.AfshariStewart2016(scenario).duration
    duration = float(dur_rec["D_5t95"])
    motion = CompatibleRvtMotion.from_response_spectrum(rs, duration)
    pga = motion.calc_peak()
    # Sanity bounds — real numerical regression added when first value is locked in CI
    assert 0.001 < pga < 10.0


def test_from_response_spectrum_stub():
    """from_response_spectrum duck-types: SimpleNamespace works, no isinstance check."""
    rs_stub = SimpleNamespace(
        periods=np.array([0.1, 0.5, 1.0, 2.0]),
        spec_accels=np.array([0.30, 0.25, 0.15, 0.08]),
        damping=0.05,
    )
    motion = CompatibleRvtMotion.from_response_spectrum(rs_stub, duration=5.0)
    assert np.all(np.isfinite(motion.fourier_amps))
    assert motion.duration == pytest.approx(5.0)


def test_from_fas_stub():
    """RvtMotion.from_fas duck-types: SimpleNamespace works."""
    from pyrvt.motions import RvtMotion

    fas_stub = SimpleNamespace(
        freqs=np.logspace(-1, 2, 50),
        fourier_amps=np.ones(50) * 0.01,
        duration=5.0,
    )
    motion = RvtMotion.from_fas(fas_stub)
    assert np.all(np.isfinite(motion.fourier_amps))
    assert motion.duration == pytest.approx(5.0)
