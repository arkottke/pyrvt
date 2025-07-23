Glossary
========

This glossary defines key terms used in pyRVT and random vibration theory.

.. glossary::

   Acceleration Response Spectrum
      A plot showing the maximum acceleration response of single-degree-of-freedom
      oscillators with different natural periods when subjected to a specific ground motion.
      Typically computed for a fixed damping ratio (commonly 5%).

   Damping
      A measure of energy dissipation in a vibrating system. In structural dynamics,
      typically expressed as a fraction of critical damping (e.g., 0.05 for 5% damping).

   Duration
      The effective duration of strong ground motion. Different definitions exist,
      including significant duration (time between 5% and 95% of cumulative energy)
      and RMS duration used in random vibration theory calculations.

   Fourier Amplitude Spectrum (FAS)
      The frequency-domain representation of a time series, showing the amplitude
      of different frequency components. For ground motions, typically expressed
      in units of g-s (acceleration times time).

   Ground Motion Prediction Equation (GMPE)
      An empirical or semi-empirical equation that predicts ground motion intensity
      measures (such as peak ground acceleration or response spectral values) based
      on earthquake source, path, and site parameters.

   Kappa (κ)
      A parameter describing high-frequency attenuation of seismic waves, particularly
      important for site effects. Often modeled as an exponential decay factor
      in the form exp(-πκf) where f is frequency.

   Moment Magnitude
      A logarithmic scale used to measure earthquake size, based on the seismic moment.
      Denoted as Mw, it provides a physically meaningful measure of earthquake strength.

   Peak Factor
      A statistical parameter relating the peak value of a random process to its
      root-mean-square (RMS) value. Central to random vibration theory calculations.

   Power Spectral Density (PSD)
      The frequency-domain representation of a random process showing how power
      is distributed across different frequencies. Related to the square of the
      Fourier amplitude spectrum.

   Random Vibration Theory (RVT)
      A statistical approach for relating frequency-domain characteristics of
      ground motion (Fourier amplitude spectrum) to peak response characteristics
      (response spectrum) using probabilistic methods.

   Response Spectrum
      See :term:`Acceleration Response Spectrum`. May also refer to velocity or
      displacement response spectra computed similarly.

   RMS (Root Mean Square)
      A statistical measure of the magnitude of a varying quantity. For random
      processes, it represents the square root of the mean of the squared values.

   Seismic Moment
      A measure of earthquake size based on the area of fault rupture, average
      slip, and rock rigidity. Used to calculate moment magnitude.

   Single-Degree-of-Freedom (SDOF)
      A simplified structural model with one degree of freedom, characterized by
      mass, damping, and stiffness. Used as the basis for response spectrum calculations.

   Spectral Acceleration
      The maximum acceleration response of a single-degree-of-freedom oscillator
      with a specific natural period and damping ratio when subjected to ground motion.

   Stochastic Simulation
      A method for generating synthetic ground motions using random processes
      constrained by statistical properties (such as Fourier amplitude spectra
      and durations) derived from empirical data or theoretical models.

   Strong Motion
      The portion of earthquake ground shaking that is of primary engineering
      interest, typically characterized by significant amplitude and frequency content
      that can cause structural damage.

   VS30
      The time-averaged shear-wave velocity in the upper 30 meters of the subsurface.
      A commonly used parameter for site characterization in earthquake engineering.

Mathematical Terms
------------------

.. glossary::

   Autocorrelation Function
      A mathematical function that describes the correlation of a signal with a delayed
      copy of itself as a function of delay time.

   First-Passage Time
      In probability theory, the first time a stochastic process reaches a specified
      threshold level. Important for calculating peak factor distributions.

   Gaussian Process
      A stochastic process where any finite collection of random variables has a
      multivariate normal distribution. Often assumed in random vibration theory.

   Rayleigh Distribution
      A probability distribution often used to model the magnitude of a vector whose
      components are normally distributed. Relevant for modeling peaks in random processes.

   Rice Distribution
      A probability distribution that generalizes the Rayleigh distribution to include
      a non-zero mean. Sometimes used in advanced peak factor calculations.

   Stationary Process
      A stochastic process whose statistical properties do not change over time.
      Ground motions are often modeled as approximately stationary over their
      strong motion duration.

Acronyms and Abbreviations
--------------------------

.. glossary::

   ASCE
      American Society of Civil Engineers

   BSSA
      Bulletin of the Seismological Society of America

   CLI
      Command-Line Interface

   FAS
      Fourier Amplitude Spectrum

   GMPE
      Ground Motion Prediction Equation

   HVSR
      Horizontal-to-Vertical Spectral Ratio

   NGAWest2
      Next Generation Attenuation West 2 project

   PGA
      Peak Ground Acceleration

   PGV
      Peak Ground Velocity

   PSA
      Pseudo-Spectral Acceleration (often used synonymously with Spectral Acceleration)

   PSHA
      Probabilistic Seismic Hazard Analysis

   RVT
      Random Vibration Theory

   SA
      Spectral Acceleration

   SDOF
      Single-Degree-of-Freedom

   USGS
      United States Geological Survey

Contributing to the Glossary
----------------------------

If you encounter terms in pyRVT documentation that are not defined here,
or if you think additional terms should be included:

1. Create an issue on `GitHub <https://github.com/arkottke/pyrvt/issues>`_
2. Suggest the term and provide a definition
3. Consider contributing a pull request to add the term

The glossary aims to be accessible to users with varying backgrounds in
earthquake engineering and random vibration theory.
