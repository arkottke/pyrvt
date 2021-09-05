import datetime
import functools
import multiprocessing

import numpy as np

import pyrvt


def calc_crm(freqs, osc_damping, peak_calculator, event):
    crm = pyrvt.motions.CompatibleRvtMotion(
        freqs,
        event["psa"],
        duration=event["duration"],
        osc_damping=osc_damping,
        peak_calculator=peak_calculator,
    )
    return crm


if __name__ == "__main__":
    # Pool of workers
    POOL = multiprocessing.Pool()

    osc_damping = 0.05
    peak_calculator = pyrvt.peak_calculators.get_peak_calculator("V75", {})

    ext, freqs, events = pyrvt.tools.read_events(
        "../examples/example_targetSa.csv", "psa"
    )

    start = datetime.datetime.now()
    crms = POOL.map(
        functools.partial(calc_crm, freqs, osc_damping, peak_calculator), events
    )
    end = datetime.datetime.now()
    print("Duration:", end - start)
