Release Notes
=============

This page summarises the high-level changes in each major release.
For the complete, granular changelog see :doc:`changelog`.

Version 2.0.0 (unreleased)
---------------------------

**Breaking changes**

- Removed ``SourceTheoryMotion`` and ``StaffordEtAl22Motion``.
  Use ``pygmm.fourier_spectrum.SourceTheoryModel`` / ``StaffordEtAl2022`` to build a
  FAS object, then pass it to :meth:`~pyrvt.motions.RvtMotion.from_fas`.
- Removed module-level helpers ``calc_stress_drop`` and ``calc_geometric_spreading``.
- Removed ``event_kwds`` parameter from ``CompatibleRvtMotion.__init__``.
  ``duration`` is now a required argument.
  Use the new :meth:`~pyrvt.motions.CompatibleRvtMotion.from_response_spectrum`
  factory instead.

**New in 2.0.0**

- :meth:`~pyrvt.motions.RvtMotion.from_fas` — create a motion from any object
  with ``.freqs``, ``.fourier_amps``, ``.duration`` attributes.
- :meth:`~pyrvt.motions.CompatibleRvtMotion.from_response_spectrum` — create a
  compatible motion from any object with ``.periods``, ``.spec_accels``,
  ``.damping``.
- :meth:`~pyrvt.motions.RvtMotion.calc_pga` and
  :meth:`~pyrvt.motions.RvtMotion.calc_pgv`.

Version 0.8.x
-------------

- Added :class:`~pyrvt.peak_calculators.SeifriedEtAl2025` (``Sea25``) peak factor model.
- Added Stafford et al. (2022) Fourier amplitude spectrum model (via pygmm).
- Fixed non-stationarity correction (moved from Vanmarcke to ToroMcGuire).
- Migrated to uv for package management.

Version 0.7.x
-------------

- Added :class:`~pyrvt.peak_calculators.BooreThompson2015` (``BT15``).
- Added :class:`~pyrvt.peak_calculators.WangRathje2018` (``WR18``).
- Refactored peak factor calculator hierarchy.

Migration Guide
---------------

Upgrading from 1.x / 0.8.x to 2.0.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Replace removed motion classes:

.. code-block:: python

   # Before (1.x)
   from pyrvt.motions import SourceTheoryMotion
   motion = SourceTheoryMotion(mag=6.5, dist=20, region="wna")

   # After (2.0)
   import pygmm.fourier_spectrum as fs
   from pyrvt.motions import RvtMotion
   fas = fs.SourceTheoryModel(mag=6.5, dist=20, region="wna")
   motion = RvtMotion.from_fas(fas)

Replace ``CompatibleRvtMotion`` construction:

.. code-block:: python

   # Before (1.x)
   from pyrvt.motions import CompatibleRvtMotion
   motion = CompatibleRvtMotion(freqs, spec_accels, event_kwds={...})

   # After (2.0)
   motion = CompatibleRvtMotion.from_response_spectrum(rs_object, duration=20.0)
