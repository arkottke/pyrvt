Quick Start
===========

This guide will get you up and running with pyRVT in just a few minutes.

Basic Usage
-----------

pyRVT provides several ways to work with ground motion data. Here are the most common use cases:

Creating a Motion from Fourier Amplitude Spectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pyrvt
   import numpy as np
   import matplotlib.pyplot as plt
   
   # Define frequency array (Hz)
   freqs = np.logspace(-1, 2, 100)  # 0.1 to 100 Hz
   
   # Define Fourier amplitude spectrum (g-s)
   # This example uses a simple flat spectrum
   fourier_amps = np.ones_like(freqs) * 0.1
   
   # Define motion duration (s)
   duration = 10.0
   
   # Create RVT motion
   motion = pyrvt.motions.RvtMotion(freqs, fourier_amps, duration)
   
   # Calculate response spectrum
   osc_freqs = np.logspace(-1, 1.5, 50)  # Oscillator frequencies
   damping = 0.05  # 5% damping
   resp_spec = motion.calc_osc_accels(osc_freqs, damping)
   
   # Plot results
   plt.figure(figsize=(10, 4))
   
   plt.subplot(1, 2, 1)
   plt.loglog(freqs, fourier_amps)
   plt.xlabel('Frequency (Hz)')
   plt.ylabel('Fourier Amplitude (g-s)')
   plt.title('Input FAS')
   plt.grid(True)
   
   plt.subplot(1, 2, 2)
   plt.loglog(1/osc_freqs, resp_spec)
   plt.xlabel('Period (s)')
   plt.ylabel('Response Acceleration (g)')
   plt.title('Calculated Response Spectrum')
   plt.grid(True)
   
   plt.tight_layout()
   plt.show()

Converting Response Spectrum to Fourier Amplitude Spectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pyrvt
   import numpy as np
   
   # Define response spectrum
   periods = np.logspace(-1, 1, 50)  # 0.1 to 10 seconds
   resp_spec = 0.1 * np.ones_like(periods)  # Flat response spectrum at 0.1g
   
   # Define motion duration
   duration = 15.0
   
   # Create compatible motion
   motion = pyrvt.motions.CompatibleRvtMotion(
       osc_freqs=1/periods,
       osc_accels=resp_spec,
       osc_damping=0.05,
       duration=duration
   )
   
   # Get the compatible Fourier amplitude spectrum
   freqs = motion.freqs
   fourier_amps = motion.fourier_amps
   
   print(f"Generated FAS with {len(freqs)} frequency points")
   print(f"Frequency range: {freqs[0]:.3f} to {freqs[-1]:.1f} Hz")

Using Different Peak Factor Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pyRVT includes several peak factor models. Here's how to use different ones:

.. code-block:: python

   import pyrvt
   import numpy as np
   
   # Create test motion
   freqs = np.logspace(-1, 2, 100)
   fourier_amps = np.ones_like(freqs) * 0.1
   duration = 10.0
   
   # Compare different peak factor models
   peak_models = [
       pyrvt.peak_calculators.Vanmarcke1975(),
       pyrvt.peak_calculators.BooreThompson2015(),
       pyrvt.peak_calculators.WangRathje2018(),
   ]
   
   osc_freqs = np.logspace(0, 1, 20)  # 1 to 10 Hz
   damping = 0.05
   
   results = {}
   for peak_calc in peak_models:
       motion = pyrvt.motions.RvtMotion(freqs, fourier_amps, duration, 
                                       peak_calculator=peak_calc)
       resp_spec = motion.calc_osc_accels(osc_freqs, damping)
       results[peak_calc.__class__.__name__] = resp_spec
   
   # Compare results
   for name, resp_spec in results.items():
       print(f"{name}: Max response = {np.max(resp_spec):.3f} g")

Command-Line Interface
----------------------

pyRVT also provides a command-line interface for quick conversions:

Convert response spectrum to Fourier amplitude spectrum:

.. code-block:: bash

   $ pyrvt rv2fa --input-file response_spectrum.csv --output-file fourier_spectrum.csv

Convert Fourier amplitude spectrum to response spectrum:

.. code-block:: bash

   $ pyrvt fa2rv --input-file fourier_spectrum.csv --output-file response_spectrum.csv

See the :doc:`../user-guide/cli` for more details on command-line usage.

Working with Seismological Models
----------------------------------

pyRVT includes seismological models for generating Fourier amplitude spectra:

.. code-block:: python

   import pyrvt
   
   # Define earthquake parameters
   magnitude = 6.5
   distance = 20.0  # km
   duration = 10.0
   
   # Create motion using source theory
   freqs = np.logspace(-1, 2, 100)
   motion = pyrvt.motions.SourceTheoryMotion(
       mag=magnitude,
       dist=distance,
       duration=duration,
       freqs=freqs
   )
   
   # Calculate response spectrum
   osc_freqs = np.logspace(-1, 1, 30)
   resp_spec = motion.calc_osc_accels(osc_freqs, damping=0.05)
   
   print(f"Peak response: {np.max(resp_spec):.3f} g")

Next Steps
----------

Now that you've seen the basics, you can:

1. Explore the :doc:`../user-guide/tutorials` for more detailed examples
2. Learn about the :doc:`../user-guide/cli` for batch processing
3. Check out the :doc:`../api/index` for complete function documentation
4. See the :doc:`../user-guide/examples` for real-world applications

Key Points to Remember
----------------------

- **Frequency units**: pyRVT uses Hz for frequency and periods in seconds
- **Amplitude units**: Fourier amplitudes are typically in g-s, response spectra in g
- **Damping**: Response spectra calculations require a damping ratio (typically 0.05 for 5%)
- **Duration**: Motion duration significantly affects peak factor calculations
- **Peak factor models**: Different models may give different results - choose based on your application
