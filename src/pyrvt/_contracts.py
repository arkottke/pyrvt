"""Consumer-side data contracts for pyrvt ↔ pygmm interop.

Canonical definitions live in ``pygmm.contracts``; this file is a private
~35-LOC duplicate so pyrvt has no runtime dependency on pygmm.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True)
class FourierSpectrum:
    """One-sided Fourier amplitude spectrum of acceleration.

    Field names match ``pygmm.contracts.FourierSpectrum``.
    """

    freqs: npt.NDArray[np.floating]
    fourier_amps: npt.NDArray[np.floating]
    duration: float


@dataclass(frozen=True)
class ResponseSpectrum:
    """Pseudo-spectral acceleration response spectrum.

    Field names match ``pygmm.contracts.ResponseSpectrum``.
    """

    periods: npt.NDArray[np.floating]
    spec_accels: npt.NDArray[np.floating]
    damping: float
    duration: Optional[float] = None
