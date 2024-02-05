# Welcome to pyRVT Documentation

`pyRVT` provides random vibration theory (RVT) models for use in earthquake
ground motion models. It provides multiple peak factor models in a common
framework such that they can be compared and tested. Additionally, it provides
an interface to define RVT based ground motion models through specificiation of
the Fourier amplitude spectrum, acceleration response spectrum, or calculated by
a seismological models.

## Peak factor models

The peak calculators and associated functions are contained within
`pyrvt.peak_calculators`. The primary peak calculators are:

- [Wang and Rathje (2018)][pyrvt.peak_calculators.WangRathje2018][@wang18]
- [Boore and BooreThompson (2015)][pyrvt.peak_calculators.BooreThompson2015][@boore15]
- [Vanmarcke (1975)][pyrvt.peak_calculators.Vanmarcke1975][@vanmarcke75]

Additional peak calculators are also provided. See the peak calculator
[API][pyrvt.peak_calculators].

## Ground motion models

Motions for calculating the peak response can be specified by the following
mechanisms. The frequency content of the motion can be specified by providing
the Fourier amplitude spectrum and during in a
[RvtMotion][pyrvt.motions.RvtMotion] instance. Alternatively, a compatible
ground motion can be computed from an acceleration response spectrum and
duration in a [CompatibileRvtMotion][pyrvt.motions.CompatibleRvtMotion]
instance. Seismological models based on a $\omega^2$-model based on the
parameters selected by Campbell (2003)[@campbell03] in a
[SourceTheoryMotion][pyrvt.motions.SourceTheoryMotion] instance, or by the
optimized functional form provided by Stafford et al. (2022)[@stafford22] in a
[StaffordEtAl22Motion][pyrvt.motions.StaffordEtAl22Motion].

## Command-line interface

A command-line interface is provided to convert response spectra to Fourier
amplitude spectra, or vice versa. This inferface as used by Al Atik et al.
(2014)[@alatik14].
