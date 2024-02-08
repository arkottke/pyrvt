# Command-line interface

`pyRVT` can be used from the command line by executing `pyrvt` with a number of arguments. These
arguments can be found by running `pyrvt`, which will produce the
following output:

```
❯ pyrvt --help
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
                        [DK85] Der Kiureghian (1985), [LP99] Liu and
                        Pezeshk (1999), [TM87] Toro and McGuire
                        (1987), [V75] Vanmarcke (1975), and [WR18]
                        Wang and Rathje (2018). If the BT12 method is
                        used, then the magnitude, distance and region
                        must be provided by the input files. If no
                        value is provided, then "V75" is used as the
                        default.
```

For example, to compute the Fourier amplitude spectra that were compatible with
target response spectrum the following command could be used:

```bash
$ pyrvt psa2fa -i examples\example_targetSa.csv
```

The required format for the events is best understood by looking at one of the
example event files. The name of the input file should include an underscore
(i.e., '\_') during the creation of the output everything to the left of the
last underscore is repeated.
