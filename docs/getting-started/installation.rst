Installation
============

pyRVT supports Python 3.10 and later versions.

Quick Installation
------------------

The easiest way to install pyRVT is using pip:

.. code-block:: bash

   $ python -m pip install pyrvt

This will install pyRVT and all required dependencies.

Development Installation
------------------------

If you want to contribute to pyRVT or use the latest development version:

1. **Clone the repository**:

   .. code-block:: bash

      $ git clone https://github.com/arkottke/pyrvt.git
      $ cd pyrvt

2. **Install development dependencies**:

   pyRVT uses `uv <https://docs.astral.sh/uv/>`_ for dependency management. First install uv:

   .. code-block:: bash

      # On macOS and Linux
      $ curl -LsSf https://astral.sh/uv/install.sh | sh

      # Or with pip
      $ pip install uv

3. **Install pyRVT in development mode**:

   .. code-block:: bash

      $ ./scripts.sh install

   This will install pyRVT with all development dependencies (testing, documentation, etc.).

4. **Verify the installation**:

   .. code-block:: bash

      $ python -c "import pyrvt; print(pyrvt.__version__)"

Verify Installation
-------------------

To verify that pyRVT is installed correctly, you can run a simple test:

.. code-block:: python

   import pyrvt
   import numpy as np
   
   # Create a simple motion
   freqs = np.logspace(-1, 2, 100)
   fourier_amps = np.ones_like(freqs)
   duration = 10.0
   
   motion = pyrvt.motions.RvtMotion(freqs, fourier_amps, duration)
   print("pyRVT installed successfully!")
   print(f"Version: {pyrvt.__version__}")

If this runs without error, pyRVT is installed correctly!

Dependencies
------------

pyRVT requires the following packages:

**Core dependencies:**
   - numpy
   - numba  
   - scipy
   - pyexcel, pyexcel-io, pyexcel-xlsx (for Excel file support)

**Optional dependencies:**
   - matplotlib (for plotting)
   - pandas (for data manipulation)
   - jupyter (for notebook examples)

All core dependencies are automatically installed when you install pyRVT.

Troubleshooting
---------------

**ImportError: No module named 'pyrvt'**
   Make sure you've activated the correct Python environment where pyRVT was installed.

**ModuleNotFoundError: No module named 'numba'**
   Try updating pip and reinstalling: ``pip install --upgrade pip && pip install --force-reinstall pyrvt``

**Issues with Excel file support**
   The pyexcel dependencies are required for reading Excel files. They should be installed automatically.

**Performance issues**
   pyRVT uses Numba for just-in-time compilation. The first run of certain functions may be slower
   as Numba compiles the code, but subsequent runs will be much faster.

Getting Help
------------

If you encounter issues:

1. Check the :doc:`../user-guide/index` for common usage patterns
2. Look at the :doc:`../api/index` for detailed function documentation  
3. Search or create an issue on `GitHub <https://github.com/arkottke/pyrvt/issues>`_
