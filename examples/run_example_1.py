"""Simple example of calculating a Fourier amplitude spectrum"""

import datetime
import os

import pyrvt

ext, periods, events = pyrvt.tools.read_events(
    os.path.join(os.path.dirname(__file__), "example_targetSa.csv"), "psa"
)
event = events[0]

target_freqs = 1.0 / periods
damping = 0.05
method = "BJ84"

calc_kwds = {key: event[key] for key in ["magnitude", "distance", "region"]}

for method in ["BJ84", "CLH56", "DK85", "LP99", "TM87", "V75"]:
    start = datetime.datetime.now()
    crm = pyrvt.motions.CompatibleRvtMotion(
        target_freqs,
        event["psa"],
        duration=event["duration"],
        osc_damping=damping,
        peak_calculator=pyrvt.tools.get_peak_calculator(method, calc_kwds),
    )
    psa_calc = crm.calc_osc_accels(target_freqs, damping)
    duration = datetime.datetime.now() - start
    print(f"{method:10}{duration}")
