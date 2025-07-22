import datetime

import dask

import pyrvt

# def inc(a):
#     return a + 1
#
# def add(a, b):
#     return a + b
#
# x = dask.delayed(inc)(1)
# y = dask.delayed(inc)(2)
# z = dask.delayed(add)(x, y)

ext, freqs, events = pyrvt.tools.read_events("../examples/example_targetSa.csv", "psa")

osc_damping = 0.05
peak_calculator = pyrvt.peak_calculators.get_peak_calculator("BJ84", {})

start = datetime.datetime.now()
crms = []
for event in events:
    crm = dask.delayed(pyrvt.motions.CompatibleRvtMotion)(
        freqs,
        event["psa"],
        duration=event["duration"],
        osc_damping=osc_damping,
        peak_calculator=peak_calculator,
    )
    crms.append(crm)

results = dask.compute(*crms)
end = datetime.datetime.now()
print("Duration:", end - start)
