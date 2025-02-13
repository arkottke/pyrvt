---
title: Use as a library
---

# Use as a library

The intent in developing `pyRVT` the was for use in site response calculations
performed by [pyStrata](https://github.com/arkottke/pystrata). However, it can
also be used for independent calculations through [peak
calculation](#peak-calculators) or through [ground motion
models](#ground-motion-models).

## Peak calculators

The `Calculator` class in **pyRVT** is the base class for implementations of
peak factor models. Each derived calculator implements methods to:

1. Compute peak factors for given duration and frequency content.
2. Optionally handle hardware- or site-specific parameters.

Typically, you get an instance by calling:

```python
from pyrvt.peak_calculators import get_peak_calculator

calc = get_peak_calculator("BJ84", {})
peak_factor = calc.calc_peak(duration=30.0, sspectrum=my_spectrum)
```

In this example, `"BJ84"` refers to a specific peak calculator (Boore and Joyner).
You can choose other calculators, such as `"BT15"` or `"V75"`, as needed.

## Ground motion models

### Overview of Classes in motions.py

- **RvtMotion**
  The base class that stores frequency arrays, Fourier amplitudes, and durations. It also provides methods for calculating oscillator accelerations and peaks.

- **SourceTheoryMotion**
  Extends RvtMotion for source-theory-based ground motions. It handles magnitude, distance, and region-specific parameters.

- **StaffordEtAl22Motion**
  Another RvtMotion subclass modeling motions via Stafford et al. (2022) relationships. It computes site amplification, duration, and Fourier amplitudes based on magnitude, distance, and other inputs.

- **CompatibleRvtMotion**
  Creates an RVT-compatible target motion from a set of oscillator frequencies and corresponding spectral accelerations.

### Example Usage (RvtMotion)

```python
from pyrvt.motions import RvtMotion
import numpy as np

# Frequency array in Hz
freqs = np.array([0.5, 1.0, 2.0])
# Sample Fourier amplitudes
fourier_amps = np.array([1e-3, 2e-3, 1.5e-3])

# Create an RVT motion instance
motion = RvtMotion(
    freqs=freqs,
    fourier_amps=fourier_amps,
    duration=30.0,
    peak_calculator="BJ84"
)

# Calculate peak amplitude
peak_value = motion.calc_peak()
print("Peak Value:", peak_value)
```

## Example Usage (SourceTheoryMotion)

```python
import numpy as np
from pyrvt.motions import SourceTheoryMotion

# Magnitude and distance
magnitude = 6.0
distance = 10.0

# Frequency array in Hz
freqs = np.logspace(-1, 2, 50)  # 0.1 to 100 Hz
# Create a SourceTheoryMotion instance
source_motion = SourceTheoryMotion(
    magnitude=magnitude,
    distance=distance,
    region="wna",          # Western North America
    stress_drop=100.0,     # Example stress drop
    freqs=freqs,
    duration=30.0,
    peak_calculator="BT15" # Boore & Thompson 2015
)

# Calculate Fourier amplitudes
fourier_amps = source_motion.calc_fourier_amps()
print("Fourier Amplitudes:", fourier_amps)

# Calculate peak amplitude
peak_value = source_motion.calc_peak()
print("Peak Value:", peak_value)
```
