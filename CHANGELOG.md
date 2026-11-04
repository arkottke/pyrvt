# Changelog

## 0.9.0 (unreleased)

- Changed: the documentation is reorganized around a user guide that describes how
  motions, peak calculators and tools fit together; the tutorials and quickstart
  code, which used keywords that do not exist, is corrected. The example notebooks,
  which used the removed `SourceTheoryMotion`, are replaced by two that are executed
  by the tests.
- Changed: `pyexcel`, `pyexcel-io` and `pyexcel-xlsx` are no longer required
  dependencies. They are only used to read and write spreadsheets and are in
  the new `cli` extra (`pip install pyRVT[cli]`); reading or writing a
  spreadsheet without them raises an `ImportError` that says so.
- Fixed: `RvtMotion.calc_arias_intensity` was half the correct value, and
  `calc_cav`, which is computed from it, was low by a factor of `2**0.496`.
- Fixed: the Boore and Thompson calculators raise `ValueError` outside the
  magnitude and distance of their coefficient tables instead of returning NaN,
  and `TM87` honours `use_nonstationarity_factor`.
- Added: `calc_osc_accels` warns when the frequency grid is too coarse to
  resolve the oscillator.
- Fixed: `read_at2` reads both NGA-West1 and NGA-West2 headers.
- Changed: adopted [SPEC 0](https://scientific-python.org/specs/spec-0000/) for
  version support -- Python minors from the last 36 months, core dependencies
  from the last 24 months. **Minimum Python is now 3.13** (was 3.10), and
  dependency lower bounds are declared explicitly: `numpy>=2.2`,
  `numba>=0.61`, `scipy>=1.15`. CI now runs 3.13 and 3.14. The dated
  drop schedule is in `pyproject.toml`.

  The old `>=3.10` floor was already unmet in practice: numpy, scipy, pandas,
  matplotlib and pint all require >=3.11, so a 3.10 install silently resolved
  to a much older stack (pandas 2.3 rather than 3.0, numpy 2.2 rather than
  2.5). That forced the lockfile to fork on `python_full_version < '3.11'` into 11 duplicated packages.

- Changed: dropped the `tomli ; python_version < '3.11'` test marker, now dead.

**Breaking changes**:

- Removed `SourceTheoryMotion` and `StaffordEtAl22Motion` — use
  `pygmm.fourier_spectrum.SourceTheoryModel` / `StaffordEtAl2022` instead.
  Pipe the result through `RvtMotion.from_fas(model)` to reproduce previous behaviour.
- Removed `calc_stress_drop` and `calc_geometric_spreading` module-level helpers.
- Removed `event_kwds` parameter from `CompatibleRvtMotion.__init__`; `duration` is
  now required. Use `CompatibleRvtMotion.from_response_spectrum(rs, duration)` or pass
  `duration` explicitly.

Added:

- `RvtMotion.from_fas(fas, peak_calculator=None)` — create a motion from any object
  with `.freqs`, `.fourier_amps`, `.duration` attributes.
- `CompatibleRvtMotion.from_response_spectrum(rs, duration, **kw)` — create a
  compatible motion from any object with `.periods`, `.spec_accels`, `.damping`.
- `RvtMotion.calc_pga()` and `.calc_pgv()`.

## 0.8.1 (2025-07-25)

- Fixed: Issue with naming in CLI tools.

## 0.8.0 (2025-07-25)

- Added Stafford et al. (2022) ground motion model.
- Added Seifried et al. (2025) peak factor model.
- Fix: Moved non-stationarity correction from Vanmarck to ToroMcGuire.
- Changed to use uv for package versioning.

## 0.7.0

- Added Boore & Thompson (2015)
- Added Wang & Rathje (2018)
- Re-factored peak factor calculators

## 0.6.3 (2017-10-10)

- Add scripts for profiling

## 0.6.2 (2017-05-07)

- Fixed PEP8 issues in docstrings.

## 0.6.0 (2017-05-05)

- Removed pyprind require and progress bars.

## 0.6.0 (2017-05-05)

- Improved performance with multiprocessing and numba.
- Changed CLI name from `rvt_operator` to `pyrvt`

## 0.5.8 (2016-11-16)

- Fixed: osc_freqs typo in tools.py
- Added: test cases for tools.py

## 0.5.7 (2016-07-14)

- Fixed version numbering.

## 0.5.6 (2016-07-14)

- Fixed manifest to include data directory.

## 0.5.5 (2016-07-11)

- Updated Travis build scripts.

## 0.5.4 (2016-07-09)

- Fixed building the documentation on RTFD

## 0.5.1 (2016-07-06)

- Fixed Travis deployment

## 0.5.0 (2016-07-06)

- Added ability to pass transfer function to PSA calculation.
- Removed extraneous files.
- Cleaned up documentation to use NumPy style.
- Fixed PEP8 issues.

## 0.1.0 (2016-03-04)

- First release on PyPI.
