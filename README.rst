IRVT
++++

.. image:: https://travis-ci.org/arkottke/pyrvt.png?branch=master   :target: https://travis-ci.org/arkottke/pyrvt


Use random vibration theory to transform between Fourier amplitude spectrum and
response spectrum.

Dependencies
============

Prior to using IRVT, Python_ and a number of packages need to be installed. Two
important points: 

1. Make sure that you are downloading Python 2.7, and 
2. Make sure that downloads are consistent with respect to 32-bit or 64-bit. If
   you are unsure of the difference, I would suggest 32-bit.

In addition to Python, the following packages need to be installed:

*Required:*

- numpy -- fast vector operations

- scipy -- indefinite integration

*Optional:*

-  matplotlib -- used for plotting

-  xlrd/xlwt -- reading and writing xls files

-  openpyxl -- reading xlsx files

For a computer running Windows, I would recommend the 32-bit installation of
Python (python-2.7.3.msi_) and installing the required packages provided by
`Christoph Gohlke`_ (e.g., numpy_ and scipy_). The optional ``xlrd``, ``xlwt``,
and ``openpyxl`` packages are provides by `base package`_ and provide support
for reading xls and xlsx files, respectively. Again, make sure that you are
downloading the proper version. For example, the appropriate file name for
``scipy`` for a 32-bit installation is ``scipy-0.10.1.win32-py2.7.exe``.

.. _Python: http://python.org/download/
.. _python-2.7.3.msi: http://python.org/ftp/python/2.7.3/python-2.7.3.msi
.. _Christoph Gohlke: http://www.lfd.uci.edu/~gohlke/pythonlibs
.. _numpy: http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
.. _scipy: http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy
.. _base package: http://www.lfd.uci.edu/~gohlke/pythonlibs/#base

Using IRVT
==========

IRVT is used by executing ``rvt_operator.py`` with a number of arguments. These
arguments can be found by running ``rvt_operator.py``, which will produce the
following output:

::
  
  C:\Users\arkottke\Documents\python\irvt>rvt_operator.py --help
  usage: rvt_operator.py [-h] -i SRC [-o DEST] [-d DAMPING] [-f] {sa2fa,fa2sa}
  
  Compute reponse or Fourier amplitude spectra using RVT.
  
  positional arguments:
    {sa2fa,fa2sa}         Operation to be performed. [sa2fa] converts from
                          (psuedo)-spectral acceleration to Fourier amplitude.
                          [fa2sa] converts from Fourier amplitude to
                          (psuedo)-spectral acceleration.
  
  optional arguments:
    -h, --help            show this help message and exit
    -i SRC, --input SRC   Path containing the input file(s). Supported file
                          types are csv, xls, and xlsx -- provided the required
                          packages have been installed. A single file or glob
                          can be specified. An example of a glob would be
                          "input/*_sa.xls" for all files within directory
                          "input" ending in "_sa.xls".
    -o DEST, --output DEST
                          Path where the output files should be created. If this
                          directory does not exist it will be created. Default:
                          ./output
    -d DAMPING, --damping DAMPING
                          Oscillator damping in decimal. Default: 0.05.
    -f, --fixed-spacing   Fixed spacing of the oscillator period of 0.01 to 10
                          sec log-spaced with 100 points. Target SA values will
                          be interpolated if needed

For example, to compute the Fourier amplitude spectra that were compatible with
target response spectrum the following command could be used: 
``rvt_operator.py sa2fa -i examples\example_targetSa.csv``

The required format for the events is best understood by looking at one of the
example event files. The name of the input file should include an underscore
(i.e., '_') during the creation of the output everything to the right of the
last underscore is repeated.

Frequently Asked Questions
==========================

1. How long does this take?

The depends on the frequency spacing and the number of events that are being
considered. If a large number of periods and/or events are specified, it might
take several minutes or longer.
