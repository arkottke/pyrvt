pyRVT Documentation
===================

A Python library and command-line application for using random vibration
theory to transform between acceleration Fourier amplitude spectrum and
acceleration response spectrum.

*pyRVT* provides random vibration theory (RVT) models for use in earthquake
ground motion models. It provides multiple peak factor models in a common
framework such that they can be compared and tested. Additionally, it provides
an interface to define RVT based ground motion models through specification of
the Fourier amplitude spectrum, acceleration response spectrum, or calculated by
seismological models.

.. grid:: 2

   .. grid-item-card:: Getting Started
      :link: getting-started/index
      :link-type: doc

      New to pyRVT? Start here for installation and basic usage.

   .. grid-item-card:: User Guide
      :link: user-guide/index
      :link-type: doc

      Learn how to use pyRVT for your projects with tutorials and examples.

   .. grid-item-card:: API Reference
      :link: api/index
      :link-type: doc

      Complete API documentation for all modules and functions.

   .. grid-item-card:: Developer Guide
      :link: developer-guide/index
      :link-type: doc

      Contributing guidelines and development information.

Quick Start
-----------

Install pyRVT with pip:

.. code-block:: bash

   $ python -m pip install pyrvt

Basic usage:

.. code-block:: python

   import pyrvt
   import numpy as np
   
   # Define frequency and Fourier amplitude spectrum
   freqs = np.logspace(-1, 2, 100)
   fourier_amps = np.ones_like(freqs)
   duration = 10.0
   
   # Create motion and calculate response spectrum
   motion = pyrvt.motions.RvtMotion(freqs, fourier_amps, duration)
   osc_freqs = np.logspace(-1, 1.5, 50)
   resp_spec = motion.calc_osc_accels(osc_freqs, damping=0.05)

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   :hidden:

   getting-started/index
   getting-started/installation
   getting-started/quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   :hidden:

   user-guide/index
   user-guide/tutorials
   user-guide/examples
   user-guide/cli

.. toctree::
   :maxdepth: 2
   :caption: Reference
   :hidden:

   api/index
   references
   glossary

.. toctree::
   :maxdepth: 2
   :caption: Development
   :hidden:

   developer-guide/index
   developer-guide/contributing
   developer-guide/changelog
   developer-guide/release-notes

.. toctree::
   :maxdepth: 1
   :caption: About
   :hidden:

   about/authors
   about/license
   about/citation

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
