Quick Start
===========

The three operations of pyRVT: the response spectrum of a Fourier amplitude
spectrum, the Fourier amplitude spectrum of a response spectrum, and the
comparison of peak factor models. :doc:`../user-guide/index` describes how the
parts fit together.

Response spectrum of a Fourier amplitude spectrum
-------------------------------------------------

.. code-block:: python

   import numpy as np
   import pyrvt

   freqs = np.geomspace(0.1, 100.0, 100)       # Hz
   fourier_amps = 0.1 * freqs / (1 + (freqs / 5.0) ** 2) ** 2   # g-s
   duration = 10.0                            # s

   motion = pyrvt.motions.RvtMotion(freqs, fourier_amps, duration)

   osc_freqs = np.geomspace(0.1, 50.0, 50)     # Hz
   spec_accels = motion.calc_osc_accels(osc_freqs, osc_damping=0.05)   # g

Fourier amplitude spectrum of a response spectrum
-------------------------------------------------

.. code-block:: python

   compatible = pyrvt.motions.CompatibleRvtMotion(
       osc_freqs, spec_accels, duration=duration, osc_damping=0.05
   )
   compatible.freqs, compatible.fourier_amps

Peak factor models
------------------

A motion is created with the abbreviation of a peak calculator. The calculators
that depend on the earthquake take its region, magnitude and distance.

.. code-block:: python

   for name in ["V75", "BJ84", "CLH56"]:
       m = pyrvt.motions.RvtMotion(freqs, fourier_amps, duration, peak_calculator=name)
       print(name, m.calc_osc_accels(osc_freqs, 0.05).max())

   calc_kwds = {"region": "wna", "mag": 6.5, "dist": 20.0}
   m = pyrvt.motions.RvtMotion(
       freqs, fourier_amps, duration, peak_calculator="BT15", calc_kwds=calc_kwds
   )

From a seismological model
--------------------------

.. code-block:: python

   from pygmm import fourier_spectrum

   fas = fourier_spectrum.SourceTheoryModel(magnitude=6.5, distance=20.0, region="wna")
   motion = pyrvt.motions.RvtMotion.from_fas(fas)

Command line
------------

.. code-block:: bash

   $ pyrvt psa2fa -i events.csv -o fourier_spectra.csv
   $ pyrvt fa2psa -i fourier_spectra.csv -o response_spectra.csv

See :doc:`../user-guide/cli`.

Units
-----

Frequencies are in Hz, periods and durations in s, Fourier amplitudes in g-s,
response spectra in g, and damping is a fraction.

Next steps
----------

- :doc:`../user-guide/index`
- :doc:`../examples/index`
- :doc:`../api/index`
