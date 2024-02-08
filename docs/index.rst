pyRVT
=====

.. toctree::
    :hidden:
    :maxdepth: 2

    user-guide
    advanced-usage
    reference/index
    contributing
    history
    zreferences

pyrvt is a Python library and command-line application for using random vibration theory
to transform between acceleration Fourier amplitude spectrum and acceleration
response spectrum.

pyrvt includes implementations of the following peak factor equations:

- :cite:t:`cartwright56`
- :cite:t:`davenport64`
- :cite:t:`vanmarcke75`
- :cite:t:`boore84`
- :cite:t:`igusa85`
- :cite:t:`toro87`
- :cite:t:`liu99`
- :cite:t:`boore12`
- :cite:t:`boore15`
- :cite:t:`wang2018`

Installing
----------

pyrvt can be installed with `pip <https://pip.pypaio>`_

.. code-block:: bash

   $ python -m pip install pyrvt

Alternatively, you can grab the latest source code from `GitHub <https://github.com/arkottke/pyrvt>`_:

.. code-block:: bash

  $ git clone https://github.com/arkottke/pyrvt.git
  $ cd pyrvt
  $ pip install .

Usage
-----

The :doc:`user-guide` is the place to go to learn how to use the command-line interface and
accomplish common tasks. The more in-depth :doc:`advanced-usage` guide is the place to go for understanding using the Python library.

The :doc:`reference/index` documentation provides API-level documentation.

Citation
--------

Please cite the software using the Zenodo DOI: https://doi.org/10.5281/zenodo.3630729.

License
-------

pyrvt is made available under the MIT License. For more details, see `LICENSE.txt <https://github.com/arkottke/pyrvt/blob/master/LICENSE.txt>`_.
