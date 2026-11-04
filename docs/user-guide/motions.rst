Motions
=======

Direct specification
--------------------

An :class:`~pyrvt.motions.RvtMotion` is defined by the frequencies and amplitudes
of a Fourier amplitude spectrum (in g-s), a duration (in s), and optionally a
peak calculator.

.. code-block:: python

   import numpy as np
   import pyrvt

   freqs = np.geomspace(0.1, 100.0, 100)
   fourier_amps = 0.1 * freqs / (1 + (freqs / 5.0) ** 2) ** 2
   motion = pyrvt.motions.RvtMotion(freqs, fourier_amps, duration=10.0)

   osc_freqs = np.geomspace(0.1, 50.0, 50)
   spec_accels = motion.calc_osc_accels(osc_freqs, 0.05)   # osc_damping = 5%

``calc_osc_accels(osc_freqs, osc_damping=0.05, trans_func=None)`` returns the
pseudo-spectral accelerations, in g. A ``trans_func`` evaluated at the
frequencies of the motion filters the spectrum first, which is the response at
another location. ``calc_peak(transfer_func=None)`` returns the peak of the motion,
``pga`` and ``pgv`` are the peaks of the acceleration and velocity,
and ``arias_intensity`` and ``cav`` are the intensity measures. When no peak
calculator is given, the default of the library is used.

From another model
------------------

``RvtMotion.from_fas(fas, peak_calculator=None, calc_kwds=None)`` takes any object
that has ``freqs``, ``fourier_amps`` and ``duration``, such as a pyGMM
Fourier amplitude spectrum model.

.. code-block:: python

   from pygmm import fourier_spectrum

   fas = fourier_spectrum.SourceTheoryModel(magnitude=6.5, distance=20.0, region="wna")
   motion = pyrvt.motions.RvtMotion.from_fas(fas, peak_calculator="V75")

Compatible with a response spectrum
-----------------------------------

A :class:`~pyrvt.motions.CompatibleRvtMotion` finds the Fourier amplitude spectrum
whose response spectrum matches a target, for a duration. It has the attributes of
an ``RvtMotion``, with ``freqs`` and ``fourier_amps`` the compatible spectrum.

.. code-block:: python

   compatible = pyrvt.motions.CompatibleRvtMotion(
       osc_freqs,
       spec_accels,
       duration=10.0,
       osc_damping=0.05,
       peak_calculator="V75",
   )

``window_len`` is the width of a moving average that smooths the computed Fourier
amplitude spectrum. The attenuation that the
spectrum implies above a frequency is ``calc_attenuation(min_freq, max_freq=None)``.

Examples
--------

- :doc:`../examples/example-01`: a motion from a source-theory spectrum.
- :doc:`../examples/example-02`: compatible motions for the events of a file.
