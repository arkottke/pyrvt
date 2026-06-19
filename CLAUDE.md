# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run all tests
hatch run test:run

# Run a single test file or test
hatch run test:run -- tests/test_motions.py
hatch run test:run -- tests/test_motions.py::test_name

# Lint and format (ruff)
ruff check src/ tests/
ruff format src/ tests/

# Build docs locally
hatch run docs:build
hatch run docs:serve

# Install in dev mode (alternative to hatch envs)
pip install -e .
```

Pre-commit hooks run ruff, ruff-format, pyupgrade, and prettier on commit.

**Note**: `tests/test_peak_calculators.py` currently fails to collect because the pinned `pystrata` version references `pyrvt.motions.SourceTheoryMotion`, which was removed. Other test files run cleanly.

## Architecture

pyRVT converts between acceleration Fourier amplitude spectra (FAS) and pseudo-spectral acceleration (PSA) using random vibration theory (RVT). The standard unit system is **g-s** (acceleration in g, time in seconds; FAS units are g·sec).

### Core modules (`src/pyrvt/`)

**`motions.py`** — RVT motion classes:
- `RvtMotion` — base class. Holds `freqs` [Hz], `fourier_amps` [g·sec], and `duration` [sec]. Key methods: `calc_peak(transfer_func)`, `calc_osc_accels(osc_freqs, osc_damping)`, `calc_pga()`, `calc_pgv()`, `calc_attenuation(min_freq)`. Frequencies must be monotonically increasing (enforced in `__init__`).
- `CompatibleRvtMotion(RvtMotion)` — iteratively fits a FAS to match a target PSA. The `from_response_spectrum(rs, duration)` classmethod accepts any duck-typed object with `.periods`, `.spec_accels`, `.damping`.
- `RvtMotion.from_fas(fas)` — builds from any object exposing `.freqs`, `.fourier_amps`, `.duration` (e.g. a `pygmm` FAS model).

**`peak_calculators.py`** — Peak factor models:
- Abstract base `Calculator`; `__call__(duration, freqs, fourier_amps)` returns a tuple where `[0]` is the peak.
- Implemented calculators (referenced by `ABBREV`): `V75`, `D64`, `DK85`, `TM87`, `CLH56`, `BJ84`, `LP99`, `BT12`, `BT15`, `WR18`.
- `BT12`, `BT15`, `WR18` require `region`, `mag`, `dist` keyword arguments — passed via `calc_kwds`.
- Factory: `get_peak_calculator(method_abbrev_or_name, calc_kwds_dict)`.
- `DEFAULT_CALC = "V75"` is the module-level default used when no calculator is specified.
- `SquaredSpectrum` caches spectral moments (lazy-computed) to avoid redundant integration.
- Performance-critical integrands are compiled with `numba.cfunc` + `scipy.integrate.quad`.

**`tools.py`** — Batch file operations:
- `read_events(fpath, response_type)` / `write_events(...)` — read/write Excel or CSV spreadsheets with a fixed row layout: 6 parameter rows (magnitude, distance, Vs30, kappa, duration, region), then one header row, then the spectral data.
- `calc_compatible_spectra(method, periods, events, damping)` — parallelizes PSA→FAS via `multiprocessing.Pool`.
- `operation_psa2fa` / `operation_fa2psa` — top-level functions called by the CLI.

**`runner.py`** — CLI entrypoint (`pyrvt psa2fa|fa2psa -i input -o output -d damping -m method`).

**`_contracts.py`** — Frozen dataclasses `FourierSpectrum` and `ResponseSpectrum` mirroring `pygmm.contracts`. Kept as a local duplicate so pyrvt has no runtime dependency on pygmm.

### Key conventions

- Oscillator damping is always fractional (decimal), e.g. `0.05` for 5%.
- `calc_sdof_tf(freqs, osc_freq, osc_damping)` returns the SDOF transfer function that converts acceleration FAS → PSA when `|H| * FAS` is integrated.
- `sort_increasing(*arrays)` reorders arrays to be monotonically increasing; raises `NotImplementedError` if non-monotonic.
- `log_spaced_values(lower, upper, per_decade=512)` generates log-spaced frequency grids.
