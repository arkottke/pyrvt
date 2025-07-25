"""Tools for reading/writing of files and performing operations."""

from __future__ import annotations

import functools
import multiprocessing
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pyexcel

from pyrvt import motions
from pyrvt.peak_calculators import get_peak_calculator, get_region

PARAMETER_NAMES = [
    ("mag", "Magnitude"),
    ("dist", "Distance (km)"),
    ("vs30", "Vs30 (m/s)"),
    ("kappa", "Site Atten., Kappa0 (sec)"),  # Site Atten., κ₀
    ("duration", "Duration (sec)"),
    ("region", "Region"),
]


def get_fpaths(pat: str) -> list[Path]:
    """Iterate over file paths."""
    return Path(".").glob(pat) if "*" in pat else [Path(pat)]


def read_events(
    fpath: str | Path, response_type: str
) -> tuple[str, np.ndarray, list[dict]]:
    """Read data from the file an Excel work book.

    Parameters
    ----------
    fpath : str or `Path`
        Filename of the input file.
    response_type : str
        Type of response. Valid options are: 'psa' for psuedo-spectral
        acceleration, or 'fa' for Fourier amplitude.

    Returns
    -------
    ext : str
        Extension of input file
    reference : :class:`numpy.ndarray`
        Reference of the response. This is either period (sec) for
        response_type 'psa' or frequency (Hz) for response_type 'fa'
    events : List[dict]
        List of events read from the file. See ``Note`` in
        :func:`.calc_compatible_spectra` for more information on structure of
        the dictionaries.

    """
    assert response_type in ["psa", "fa"]
    fpath = Path(fpath)

    data = pyexcel.get_array(file_name=str(fpath))
    ext = fpath.suffix

    parameters = {key: data[i][1:] for i, (key, label) in enumerate(PARAMETER_NAMES)}

    event_row = len(parameters) + 1
    event_count = len(data[0]) - 1

    reference = np.array([row[0] for row in data[event_row:]])

    events = []
    for i in range(event_count):
        resps = np.array([row[i + 1] for row in data[event_row:]])
        # Extract the appropriate attributes
        e = {k: v[i] for k, v in parameters.items()}
        e[response_type] = resps

        if "region" in e:
            e["region"] = get_region(e["region"])

        events.append(e)

    return ext, reference, events


def write_events(
    fname: str | Path,
    reference: npt.ArrayLike,
    reference_label: str,
    response_type: str,
    response_label: str,
    events: list[dict],
):
    """Write the events to a file.

    Parameters
    ----------
    fname : str
        Save the events to this file. The directory is created if needed.
    reference : array_like
        Periods of the oscillator response shared across all events.
    reference_label : str
        Label of the reference (e.g., 'Frequency (Hz)').
    response_type : str
        Type of response. Valid options: `psa` for pseudo-spectral
        accleration, or `fa` for Fourier amplitude.
    response_label : str
        Label of the response type (e.g., 'Fourier Ampl. (g/sec)')
    events : List[dict]
        Events to write to file. See ``Note`` in
        :func:`.calc_compatible_spectra` for more information.

    Raises
    ------
    NotImplementedError:
        If extension is not supported

    """
    # Create the rows of output
    rows = []
    # Output the parameters
    for key, label in PARAMETER_NAMES:
        rows.append([label] + [e[key] for e in events])

    rows.append([reference_label] + len(events) * [response_label])

    # Output the response spectra
    for i in range(len(reference)):
        rows.append([reference[i]] + [e[response_type][i] for e in events])

    # Create the parent directory
    fpath = Path(fname)
    fpath.parent.mkdir(parents=True, exist_ok=True)

    # Write the file
    pyexcel.save_as(array=rows, dest_file_name=str(fpath))


def _calc_fa(target_freqs, damping, method, event):
    """Calculate the fourier amplitudes for an event.

    Note that this is intended as a helper function to be called by
    multiprocessing.Pool.
    """
    event_keys = ["mag", "dist", "region"]
    event_kwds = {key: event[key] for key in event_keys}
    crm = motions.CompatibleRvtMotion(
        target_freqs,
        event["psa"],
        duration=event["duration"],
        osc_damping=damping,
        event_kwds=event_kwds,
        peak_calculator=get_peak_calculator(method, event_kwds),
    )
    psa_calc = crm.calc_osc_accels(target_freqs, damping)
    return crm, psa_calc


def calc_compatible_spectra(
    method: str, periods: npt.ArrayLike, events: list[dict], damping: float = 0.05
) -> tuple[np.ndarray, list[dict]]:
    """Compute the response spectrum compatible motions.

    Parameters
    ----------
    method : str
        RVT peak factor method, see
        :func:`~.peak_calculators.get_peak_calculator`.
    periods : array_like
        Periods of the oscillator response shared across all events.
    events : List[dict]
        All events to consider. See ``Note``.
    damping : float, optional
        Fractional damping of the oscillator (decimal). Default value of 0.05
        for a damping ratio of 5%.

    Returns
    -------
    :class:`numpy.ndarray`
        Frequency of the computed Fourier amplitude spectra.

    dicts
        Modified events

    Note
    ----
    Each event dictionary should have the following keys:

    - **psa** : :class:`numpy.ndarray` -- pseudo-spectral accelerations. This
      is the target for the :class:`~.motions.CompatibleRvtMotion`.
    - **duration** : float, optional -- duration of the ground motion
    - **magnitude** : float, optional -- earthquake magnitude
    - **distance** : float, optional -- earthquake distance (km)
    - **region** : str, optional -- earthquake source region, see
      :func:`~.peak_calculators.get_region` If no duration is provided one
      is estimated from the magnitude, distance, and region.

    The `events` dictionary is modified by this function and adds the
    following keys:

    - **duration** : float -- duration of the ground motion if one was not
      specified
    - **fa** : :class:`numpy.ndarray` -- Fourier amplitude spectra in units of
      g*sec
    - **psa_calc** : :class:`numpy.ndarray` -- Pseudo-spectral acceleration
      calculated from `fa`. This will differ slightly from `psa_target`.

    """
    target_freqs = 1.0 / periods
    with multiprocessing.Pool() as pool:
        results = pool.map(
            functools.partial(_calc_fa, target_freqs, damping, method), events
        )

    # Copy values back into the dictionary
    modified = []
    for event, (crm, psa_calc) in zip(events, results):
        event = dict(event)
        if not event["duration"]:
            event["duration"] = crm.duration
        event["fa"] = crm.fourier_amps
        event["psa_calc"] = psa_calc

        modified.append(event)

    # Return the frequency from one of the computed motions.
    freqs = results[0][0].freqs
    return freqs, modified


def operation_psa2fa(
    src: str | Path,
    dst: str | Path,
    damping: float,
    method: str = "LP99",
    fixed_spacing: bool = True,
    verbose: bool = True,
):
    """Compute the acceleration response spectrum from a Fourier amplitude spectrum.

    Parameters
    ----------
    src : str
        Source for the pseudo-spectral accelerations (PSA). This can be a
        filename or pattern used in :func:`glob.glob`.
    dst : str
        Destination directory for the output PSA. The directory is created if
        it does not exist.
    damping : float
        Fractional damping of the oscillator ( decimal).
    method : str
        RVT peak factor method, see
        :func:`~.peak_calculators.get_peak_calculator`.
    fixed_spacing : bool, optional
        If `True`, then the periods are interpolated to 301 points equally
        space in log-space from 0.01 to 10.
    verbose : bool, optional
        Print status of calculation.

    """
    dst = Path(dst)
    dst.mkdir(parents=True, exist_ok=True)

    for fpath in get_fpaths(src):
        if verbose:
            print("Processing:", fpath)
        ext, periods, events = read_events(fpath, "psa")

        if fixed_spacing:
            # Interpolate the periods to a smaller range
            _periods = np.logspace(-2, 1, 301)

            for e in events:
                e["psa"] = np.exp(np.interp(_periods, periods, np.log(e["psa"])))

            periods = _periods

        # Compute the FA from the PSA
        freqs, events = calc_compatible_spectra(
            method, periods, events, damping=damping
        )

        basename = fpath.stem.rsplit("_", 1)[0]

        write_events(
            dst / (basename + "_sa" + ext),
            periods,
            "Period (s)",
            "psa_calc",
            "Sa (g)",
            events,
        )
        write_events(
            dst / (basename + "_fa" + ext),
            freqs,
            "Frequency (Hz)",
            "fa",
            "FA (g-s)",
            events,
        )


def _calc_psa(
    osc_freqs: npt.ArrayLike,
    damping: float,
    method: str,
    freqs: npt.ArrayLike,
    event: dict[str, float],
) -> np.ndarray:
    """Calculate the response spectra for an event.

    Note that this is intended as a helper function to be called by
    multiprocessing.Pool.
    """
    m = motions.RvtMotion(
        freqs=freqs,
        fourier_amps=event["fa"],
        duration=event["duration"],
        peak_calculator=get_peak_calculator(
            method,
            {
                "region": event["region"],
                "mag": event["mag"],
                "dist": event["dist"],
            },
        ),
    )
    return m.calc_osc_accels(osc_freqs, damping)


def operation_fa2psa(
    src: str | Path,
    dst: str | Path,
    damping: float,
    method: str = "LP99",
    fixed_spacing: bool = True,
    verbose: bool = True,
):
    """Compute the Fourier amplitude spectrum from a accel. response spectrum.

    Parameters
    ----------
    src : str
        Source for the Fourier amplitudes. This can be a filename or pattern
        used in :func:`glob.glob`.
    dst : str
        Destination directory for the output PSA. The directory is created if
        it does not exist.
    damping : float
        Fractional damping of the oscillator (decimal).
    method : str
        RVT peak factor method, see
        :func:`~.peak_calculators.get_peak_calculator`.
    fixed_spacing : bool, optional
        If `True`, then the periods are interpolated to 301 points equally
        space in log-space from 0.01 to 10.

    """
    if fixed_spacing:
        periods = np.logspace(-2, 1, 301)
        osc_freqs = 1.0 / periods

    dst = Path(dst)
    dst.mkdir(parents=True, exist_ok=True)

    for fpath in get_fpaths(src):
        if verbose:
            print("Processing:", fpath)
        ext, freqs, events = read_events(fpath, "fa")

        if not fixed_spacing:
            osc_freqs = freqs
            periods = 1.0 / osc_freqs

        with multiprocessing.Pool() as pool:
            psas = pool.map(
                functools.partial(_calc_psa, osc_freqs, damping, method, freqs), events
            )

        for event, psa in zip(events, psas):
            event["psa"] = psa

        basename = fpath.stem.rsplit("_", 1)[0]
        write_events(
            dst / (basename + "_sa" + ext),
            periods,
            "Period (s)",
            "psa",
            "Sa (g)",
            events,
        )


def read_at2(fname: str) -> tuple[float, np.ndarray]:
    """Read an AT2 formatted file."""

    with open(fname) as fp:
        # Skip the first three header lines
        for _ in range(3):
            next(fp)
        parts = next(fp).replace(",", "").split()
        time_step = float(parts[3].replace("SEC", ""))

        accels = np.array([float(part) for line in fp for part in line.split()])

    return time_step, accels
