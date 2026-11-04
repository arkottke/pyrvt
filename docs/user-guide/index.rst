.. _usage:

User Guide
==========

pyRVT relates the Fourier amplitude spectrum of a ground motion to its response
spectrum. The relationship needs a peak factor model and a duration. This page
describes the pieces and how they are combined; the pages that follow describe
each of them.

.. toctree::
   :maxdepth: 2

   motions
   peak-calculators
   cli

What it computes
----------------

The response of an oscillator to a motion with a Fourier amplitude spectrum
:math:`|A(f)|` has a root-mean-square (RMS) value that follows from Parseval's
theorem, and the peak response is that RMS multiplied by a peak factor. The peak
factor depends on the duration of the motion and on the shape of the spectrum, and
there are several models of it. The glossary defines the terms.

The pieces
----------

**Peak calculator** (:doc:`peak-calculators`)
    A :class:`~pyrvt.peak_calculators.Calculator` computes the peak factor. A
    calculator is selected by an abbreviation of the authors and year,
    such as ``"V75"`` for Vanmarcke (1975), with
    :func:`~pyrvt.peak_calculators.get_peak_calculator`, or by passing its name to
    a motion.

**Motion** (:doc:`motions`)
    A :class:`~pyrvt.motions.RvtMotion` has a Fourier amplitude spectrum
    (``freqs`` and ``fourier_amps``), a ``duration`` and a peak calculator. Its
    methods compute the response of oscillators (``calc_osc_accels``), the peak of
    the motion or of its transfer function to another location (``calc_peak``), and
    the peak ground motion (``pga``, ``pgv``). A
    :class:`~pyrvt.motions.CompatibleRvtMotion` is built the other way around: from
    a target response spectrum, it finds the Fourier amplitude spectrum whose
    response spectrum matches it.

**Tools and the command line** (:doc:`cli`)
    :mod:`pyrvt.tools` reads and writes events (response spectra with their
    magnitude, distance, region and duration) and reads accelerograms.
    ``pyrvt`` converts between response spectra and Fourier amplitude spectra
    for the events of a file.

How they fit together
---------------------

.. code-block:: text

   Fourier spectrum (freqs, fourier_amps) ─┐
   duration ───────────────────────────────┼──► RvtMotion ──► calc_osc_accels ──► response spectrum
   peak calculator ────────────────────────┘        │
                                                    └────────► calc_peak(transfer_func)

   target response spectrum ─┐
   duration                  ├──► CompatibleRvtMotion ──► fourier_amps (compatible)
   peak calculator ──────────┘

.. code-block:: python

   import numpy as np
   import pyrvt

   freqs = np.geomspace(0.1, 100.0, 100)
   fourier_amps = 0.1 * freqs / (1 + (freqs / 5.0) ** 2) ** 2

   motion = pyrvt.motions.RvtMotion(
       freqs, fourier_amps, duration=10.0, peak_calculator="V75"
   )
   osc_freqs = np.geomspace(0.1, 50.0, 50)
   spec_accels = motion.calc_osc_accels(osc_freqs, 0.05)

   target = motion.calc_osc_accels(osc_freqs, 0.05)
   compatible = pyrvt.motions.CompatibleRvtMotion(
       osc_freqs, target, duration=10.0, peak_calculator="V75"
   )
   compatible.freqs, compatible.fourier_amps

Other libraries
---------------

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Library
     - Connection
     - Where
   * - pyGMM
     - Any object with ``freqs``, ``fourier_amps`` and ``duration``, such as a pyGMM
       Fourier amplitude spectrum model, is the input of ``RvtMotion.from_fas``
     - :doc:`motions`
   * - pyStrata
     - Its RVT motions are pyRVT motions, and it uses ``calc_peak`` with the
       transfer function of a site response calculation to get the peak at the
       surface
     - :doc:`motions`

Examples
--------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Topic
     - Example
   * - Response spectra from the same Fourier spectrum with each peak calculator
     - :doc:`../examples/example-01`
   * - A Fourier spectrum compatible with a target response spectrum
     - :doc:`../examples/example-02`
