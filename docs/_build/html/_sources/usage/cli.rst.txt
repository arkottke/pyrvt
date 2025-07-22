Command-line interface
======================

*pyRVT* can be used from the command line by executing ``pyrvt`` with a number of arguments. These
arguments can be found by running ``pyrvt``, which will produce the
following output:

.. code-block:: console

   $ pyrvt --help
   usage: pyrvt [-h] [--version] -i INPUT [-o OUTPUT] [-d DAMPING] [-f]
                [-m {BJ84,BT12,DK85,LP99,TM87,V75,WR18}]
                {psa2fa,fa2psa}

   Compute response or Fourier amplitude spectra using RVT.

   positional arguments:
     {psa2fa,fa2psa}       Operation to be performed: [psa2fa] converts
                           from pseudo-spectral acceleration to Fourier
                           amplitude, and [fa2psa] converts from Fourier
                           amplitude to pseudo-spectral acceleration.

   options:
     -h, --help            show this help message and exit
     --version             show program's version number and exit

     -i INPUT, --input INPUT
                           Path containing the input file(s). Supported
                           file types are csv, xls, and xlsx -- provided
                           the required packages have been installed. A
                           single file or glob can be specified. An
                           example of a glob would be "input/*_sa.xls"
                           for all files within directory "input" ending
                               in "_sa.xls".

     -o OUTPUT, --output OUTPUT
                           Path where the output files should be created.
                           If this directory does not exist it will be
                           created. Default: ./output

     -d DAMPING, --damping DAMPING
                           Oscillator damping in decimal. Default: 0.05.

     -f, --fixed-spacing   Fixed spacing of the oscillator period of 0.01
                           to 10 sec log-spaced with 100 points. Target
                           SA values will be interpolated if needed

     -m {BJ84,BT12,DK85,LP99,TM87,V75,WR18},
     --method {BJ84,BT12,DK85,LP99,TM87,V75,WR18}
                           Specify the peak factor calculation method.
                           Possible options are: [BJ84] Boore and Joyner
                           (1984), [BT12] Boore and Thompson (2012),

Usage Examples
--------------

Convert from response spectrum to Fourier amplitude spectrum:

.. code-block:: bash

   pyrvt psa2fa -i input_response_spectra.csv -o output/ -m V75

Convert from Fourier amplitude spectrum to response spectrum:

.. code-block:: bash

   pyrvt fa2psa -i input_fourier_spectra.csv -o output/ -d 0.05

Process multiple files with fixed spacing:

.. code-block:: bash

   pyrvt psa2fa -i "input/*_sa.xlsx" -o output/ -f -m WR18

Input File Format
-----------------

The input files should contain data in a tabular format with:

- First column: Frequency (Hz) or Period (sec) depending on conversion direction
- Subsequent columns: Amplitude values for different ground motions
- Headers are recommended but not required

Supported file formats include CSV, XLS, and XLSX (with appropriate dependencies installed).

Output
------

The output files will be saved in the specified output directory with the same format as the input files.
The converted values will maintain the same structure and naming convention with appropriate suffixes
to indicate the conversion performed.
