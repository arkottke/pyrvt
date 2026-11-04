Peak calculators
================

A :class:`~pyrvt.peak_calculators.Calculator` computes the peak factor of a
motion. The calculators differ in the statistics that they assume for the
oscillator response, and in how they account for its duration.

Selecting a calculator
----------------------

A calculator is selected by an abbreviation, which is the initials of the authors
and the last two digits of the year. A motion accepts the abbreviation, and
:func:`~pyrvt.peak_calculators.get_peak_calculator` returns the calculator.

.. code-block:: python

   from pyrvt.peak_calculators import get_peak_calculator

   calc = get_peak_calculator("V75", None)

The second argument is a dictionary of the keywords of the calculator, which the
calculators that depend on the earthquake need: ``region``, ``mag`` and ``dist``.

.. code-block:: python

   calc = get_peak_calculator("BT15", {"region": "wna", "mag": 6.5, "dist": 20.0})

Available calculators
---------------------

.. list-table::
   :header-rows: 1
   :widths: 15 55 30

   * - Abbreviation
     - Calculator
     - Keywords
   * - ``BJ84``
     - :class:`~pyrvt.peak_calculators.BooreJoyner1984`
     -
   * - ``CLH56``
     - :class:`~pyrvt.peak_calculators.CartwrightLonguetHiggins1956`
     -
   * - ``D64``
     - :class:`~pyrvt.peak_calculators.Davenport1964`
     -
   * - ``DK85``
     - :class:`~pyrvt.peak_calculators.DerKiureghian1985`
     -
   * - ``LP99``
     - :class:`~pyrvt.peak_calculators.LiuPezeshk1999`
     -
   * - ``TM87``
     - :class:`~pyrvt.peak_calculators.ToroMcGuire1987`
     -
   * - ``V75``
     - :class:`~pyrvt.peak_calculators.Vanmarcke1975`
     -
   * - ``Sea25``
     - :class:`~pyrvt.peak_calculators.SeifriedEtAl2025`
     -
   * - ``BT12``
     - :class:`~pyrvt.peak_calculators.BooreThompson2012`
     - ``region``, ``mag``, ``dist``
   * - ``BT15``
     - :class:`~pyrvt.peak_calculators.BooreThompson2015`
     - ``region``, ``mag``, ``dist``
   * - ``WR18``
     - :class:`~pyrvt.peak_calculators.WangRathje2018`
     - ``region``, ``mag``, ``dist``

The same abbreviations are the choices of the ``-m`` option of the command line
(:doc:`cli`).

Examples
--------

- :doc:`../examples/example-01`: the response spectrum of one Fourier spectrum with
  each calculator.
- :doc:`../examples/example-02`: compatible Fourier spectra with each calculator.
