Use as a library
================

The intent in developing *pyRVT* the was for use in site response calculations
performed by `pyStrata <https://github.com/arkottke/pystrata>`_. However, it can
also be used for independent calculations through :ref:`peak-calculators` or through :ref:`ground-motion-models`.

.. _peak-calculators:

Peak calculators
----------------

The :class:`~pyrvt.peak_calculators.Calculator` class in **pyRVT** is the base class for implementations of
peak factor models. Each derived calculator implements methods to:

1. Compute peak factors for given duration and frequency content.
2. Optionally handle hardware- or site-specific parameters.

Typically, you get an instance by calling:

.. code-block:: python

   from pyrvt.peak_calculators import get_peak_calculator

   calc = get_peak_calculator("BJ84", {})
   peak_factor = calc.calc_peak(duration=30.0, sspectrum=my_spectrum)

In this example, ``"BJ84"`` refers to a specific peak calculator (Boore and Joyner).

Available Peak Calculators
~~~~~~~~~~~~~~~~~~~~~~~~~~

The following peak calculators are available in pyRVT:

- :class:`~pyrvt.peak_calculators.SeifriedEtAl2025`: Seifried et al. (2025) model
- :class:`~pyrvt.peak_calculators.WangRathje2018`: Wang and Rathje (2018) model  
- :class:`~pyrvt.peak_calculators.BooreThompson2015`: Boore and Thompson (2015) model
- :class:`~pyrvt.peak_calculators.Vanmarcke1975`: Vanmarcke (1975) model
- :class:`~pyrvt.peak_calculators.Davenport1964`: Davenport (1964) model
- :class:`~pyrvt.peak_calculators.DerKiureghian1985`: Der Kiureghian (1985) model

.. _ground-motion-models:

Ground motion models
--------------------

Motion models in pyRVT provide different ways to specify the frequency content and characteristics
of earthquake ground motions for RVT calculations.

The main motion classes are:

- :class:`~pyrvt.motions.RvtMotion`: Direct specification of Fourier amplitude spectrum
- :class:`~pyrvt.motions.CompatibleRvtMotion`: Derived from acceleration response spectrum
- :class:`~pyrvt.motions.SourceTheoryMotion`: Seismological ω²-model based motion
- :class:`~pyrvt.motions.StaffordEtAl22Motion`: Optimized functional form motion model

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   from pyrvt import motions, peak_calculators
   
   # Create frequency array
   freqs = np.logspace(-1, 2, 100)
   
   # Define Fourier amplitude spectrum
   fourier_amps = np.exp(-freqs/10) * 100
   
   # Create RVT motion
   motion = motions.RvtMotion(
       freqs=freqs,
       fourier_amps=fourier_amps,
       duration=20.0
   )
   
   # Get peak calculator
   calc = peak_calculators.get_peak_calculator("V75")
   
   # Calculate response spectrum
   periods = np.logspace(-1, 1, 50) 
   damping = 0.05
   
   resp_spec = motion.calc_response_spectrum(
       periods=periods,
       damping=damping,
       peak_calculator=calc
   )
