import re
from pathlib import Path

import numpy as np
import pyrvt
import pytest
from numpy.testing import assert_allclose


def test_calc_attenuation():
    m = pyrvt.motions.SourceTheoryMotion(5.5, 0, "cena", depth=1)
    m.calc_fourier_amps()

    atten, r_value, freqs, fitted = m.calc_attenuation(50)

    assert_allclose(0.006, atten, rtol=0.01)
    assert_allclose(1.0, r_value, rtol=0.01)


def test_compatible_rvt_motion():
    # Compute the target from the point source model.
    target = pyrvt.motions.SourceTheoryMotion(
        6.0, 20.0, "wna", peak_calculator=pyrvt.peak_calculators.DerKiureghian1985()
    )
    target.calc_fourier_amps(np.logspace(-1.5, 2, 1024))

    osc_freqs = np.logspace(-1, 2, num=50)
    osc_accels_target = target.calc_osc_accels(osc_freqs, 0.05)

    compat = pyrvt.motions.CompatibleRvtMotion(
        osc_freqs,
        osc_accels_target,
        duration=target.duration,
        osc_damping=0.05,
        peak_calculator=pyrvt.peak_calculators.DerKiureghian1985(),
    )

    osc_accels_compat = compat.calc_osc_accels(osc_freqs, 0.05)

    # Might be off by a few percent because of difficulties with the inversion.
    assert_allclose(osc_accels_target, osc_accels_compat, rtol=0.03, atol=0.05)


def iter_fas_test_cases():
    """Iterate through the FAS examples provided by P. Stafford."""
    srcpath = Path(__file__).parent / "data"

    for fpath in srcpath.glob("PJSfasSpectraDurationScaling_*.csv"):
        tests = np.genfromtxt(fpath, delimiter=",", names=True)

        freqs = tests["Frequency_Hz"]
        method = "continuous" if "Optimal" in fpath.name else "trilinear"

        for n in tests.dtype.names:
            if "FAS" not in n:
                continue
            key = n.split("_", 1)[1]
            m = re.search(r"M(?P<mag>\d+)_R(?P<dist>\d+)", key)
            mag, dist = (float(p) for p in m.groups())

            yield {
                "mag": mag,
                "dist": dist,
                "dur_ex": tests[f"Dex_{key}"][0],
                "dur_rms": tests[f"Drms_{key}"],
                "method": method,
                "freqs": freqs,
                # Convert from m/s to cm/s
                "fourier_amps": tests[n] * 100,
            }


@pytest.fixture(params=iter_fas_test_cases())
def stafford_fas(request):
    case = request.param
    # Compute the values

    mot = pyrvt.motions.StaffordEtAl22Motion(
        case["mag"],
        dist_jb=case["dist"],
        disable_site_amp=False,
        method=case["method"],
        freqs=case["freqs"],
    )

    return {
        "freqs": case["freqs"],
        "scenario": {k: case[k] for k in ["mag", "dist", "method"]},
        "desired": {
            "fourier_amps": case["fourier_amps"],
            "duration": case["dur_ex"],
        },
        "actual": {
            "fourier_amps": mot.fourier_amps,
            "duration": mot.duration,
            "dist_ps": mot.dist_ps,
        },
    }


# @pytest.mark.parametrize("key", ["fourier_amps", "duration"])
def test_stafford_fas(stafford_fas):
    assert_allclose(
        stafford_fas["actual"]["fourier_amps"],
        stafford_fas["desired"]["fourier_amps"],
        rtol=0.002,
    )

    assert_allclose(
        stafford_fas["actual"]["duration"],
        stafford_fas["desired"]["duration"],
        rtol=0.002,
    )
