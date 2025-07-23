Getting Started
===============

Welcome to pyRVT! This section will help you get up and running with the library.

pyRVT is a Python library and command-line application for using random vibration
theory to transform between acceleration Fourier amplitude spectrum and
acceleration response spectrum.

.. toctree::
   :maxdepth: 2

   installation
   quickstart

What is Random Vibration Theory?
--------------------------------

Random vibration theory (RVT) is a statistical method used to relate the frequency content
of earthquake ground motions (described by Fourier amplitude spectra) to their peak
response characteristics (described by response spectra). This relationship is particularly
useful in earthquake engineering for:

- Converting between Fourier amplitude spectra and response spectra
- Ground motion prediction and simulation
- Seismic hazard analysis
- Structural response analysis

Key Concepts
------------

**Fourier Amplitude Spectrum (FAS)**
   The frequency-domain representation of ground motion acceleration, showing
   the amplitude of different frequency components.

**Response Spectrum**
   A plot showing the maximum response of single-degree-of-freedom oscillators
   with different natural periods to a given ground motion.

**Peak Factor**
   A statistical parameter that relates the RMS (root-mean-square) response
   to the peak response of a system.

**Duration**
   The effective duration of strong ground motion, which affects the peak factor
   calculation.

Why Use pyRVT?
--------------

- **Multiple Peak Factor Models**: Choose from several well-established peak factor models
- **Flexible Interface**: Work with Fourier amplitude spectra, response spectra, or seismological models
- **Command-line Interface**: Use pyRVT directly from the command line for quick conversions
- **Well-tested**: Extensively tested against published results and other implementations
- **Open Source**: MIT licensed and actively maintained

Next Steps
----------

1. :doc:`installation` - Install pyRVT on your system
2. :doc:`quickstart` - Learn basic usage with simple examples
3. :doc:`../user-guide/tutorials` - Dive deeper with comprehensive tutorials
