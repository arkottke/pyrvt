import numpy as np
import pyrvt
import pytest
from numpy.testing import assert_allclose


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
