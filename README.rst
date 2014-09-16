pyRVT
+++++

.. image:: https://travis-ci.org/arkottke/pyrvt.svg?branch=master   
    :target: https://travis-ci.org/arkottke/pyrvt

.. image:: https://coveralls.io/repos/arkottke/pyrvt/badge.png?branch=master 
    :target: https://coveralls.io/r/arkottke/pyrvt?branch=master 


Use random vibration theory to transform between Fourier amplitude spectrum and
response spectrum.

Dependencies
============

Prior to using pyRVT, Python and a number of packages need to be installed. In
addition to Python, the following packages need to be installed:

*Required:*

- numpy -- fast vector operations

- scipy -- indefinite integration

*Optional:*

-  matplotlib -- used for plotting

-  xlrd/xlwt -- reading and writing xls files

-  openpyxl -- reading xlsx files

Install Python dependencies is best accomplished with a package manager. On
Windows or OS-X, I recommend using Miniconda.

Minconda has installers for `Windows 32-bit`_, `Windows 64-bit`_, and `OS-X`_.
See the full list of installers `here` -- make sure you select Miniconda3 for
Python3.

.. _Windows 32-bit: http://repo.continuum.io/miniconda/Miniconda3-3.3.0-Windows-x86.exe
.. _Windows 64-bit: http://repo.continuum.io/miniconda/Miniconda3-3.3.0-Windows-x86_64.exe
.. _OS-X: http://repo.continuum.io/miniconda/Miniconda3-3.3.0-MacOSX-x86_64.sh

After the installer is finished, install the required dependencies by opening a
terminal. On Windows, this is best accomplished with ``Windows Key + r``, enter
"cmd". Next enter the following command:

::
 
  conda install --yes setuptools numpy scipy matplotlib nose openpyxl xlrd pip

On Windows, the text can copied and pasted if "Quick Edit" mode is enabled. To
enable this feature, right click on the icon in the upper left portion of the
window, and select "Properties", and then check the "Quick Edit Mode" check box
within the "Edit Options" group. Copy the text, and then paste it by click the
right mouse button.

Next, install or upgrade pyRVT:

::

  pip install --upgrade https://github.com/arkottke/pyrvt/archive/master.zip


Using pyRVT
===========

pyRVT is used by executing ``rvt_operator`` with a number of arguments. These
arguments can be found by running ``rvt_operator``, which will produce the
following output:

::
  
  C:\Users\arkottke\Documents\>rvt_operator --help
  usage: runner.py [-h] -i SRC [-o DST] [-d DAMPING] [-f] [-m METHOD]
                   {psa2fa,fa2psa}
  
  Compute response or Fourier amplitude spectra using RVT.
  
  positional arguments:
    {psa2fa,fa2psa}       Operation to be performed. [psa2fa] converts from
                          pseudo-spectral acceleration to Fourier amplitude.
                          [fa2psa] converts from Fourier amplitude to pseudo-
                          spectral acceleration.
  
  optional arguments:
    -h, --help            show this help message and exit
    -i SRC, --input SRC   Path containing the input file(s). Supported file
                          types are csv, xls, and xlsx -- provided the required
                          packages have been installed. A single file or glob
                          can be specified. An example of a glob would be
                          "input/*_sa.xls" for all files within directory
                          "input" ending in "_sa.xls".
    -o DST, --output DST  Path where the output files should be created. If this
                          directory does not exist it will be created. Default:
                          ./output
    -d DAMPING, --damping DAMPING
                          Oscillator damping in decimal. Default: 0.05.
    -f, --fixed-spacing   Fixed spacing of the oscillator period of 0.01 to 10
                          sec log-spaced with 100 points. Target SA values will
                          be interpolated if needed
    -m METHOD, --method METHOD
                          Specify the peak factor calculation method. Possible
                          options are: BJ84: Boore and Joyner (1984) BT12: Boore
                          and Thompson (2012) DK85: Der Kiureghian (1985) LP99:
                          Liu and Pezeshk (1999) [default] TM87: Toro and
                          McGuire (1987) V75: Vanmarcke (1975) If the BT12
                          method is used, then the magnitude, distance and
                          region must be provided. The LP99 is the default
                          method as it provides correction for oscillation
                          duration without the need for specifying magnitude,
                          distance, and region.

For example, to compute the Fourier amplitude spectra that were compatible with
target response spectrum the following command could be used: 
``rvt_operator psa2fa -i examples\example_targetSa.csv``

The required format for the events is best understood by looking at one of the
example event files. The name of the input file should include an underscore
(i.e., '_') during the creation of the output everything to the left of the
last underscore is repeated.

Frequently Asked Questions
==========================

1. How long does this take?

The depends on the frequency spacing and the number of events that are being
considered. If a large number of periods and/or events are specified, it might
take several minutes or longer.
