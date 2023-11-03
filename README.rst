ARC Echelle Spectroscopic Observation Pipeline (aesop)
------------------------------------------------------

.. image:: https://github.com/bmorris3/aesop/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/bmorris3/aesop/actions/workflows/ci.yml
   :alt: Testing status

.. image:: https://readthedocs.org/projects/arces/badge/?version=latest
    :target: https://arces.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

.. image:: https://joss.theoj.org/papers/10.21105/joss.00854/status.svg
    :target: https://doi.org/10.21105/joss.00854
    :alt: Paper

The ARC Echelle Spectroscopic Observation Pipeline, or ``aesop``, is a high resolution
spectroscopy software toolkit tailored for observations from the `ARC Echelle Spectrograph (ARCES)
<https://www.apo.nmsu.edu/arc35m/Instruments/ARCES/>`_ mounted on the
`ARC 3.5 m Telescope <https://www.apo.nmsu.edu/arc35m/>`_ at
`Apache Point Observatory <https://www.apo.nmsu.edu>`_. ``aesop`` picks up where the
traditional IRAF reduction scripts leave off, offering an open development,
object-oriented Pythonic analysis framework for echelle spectra.

Basic functionality of ``aesop`` includes: (1) blaze function normalization by polynomial
fits to observations of early-type stars, (2) an additional/alternative robust least-squares
normalization method, (3) radial velocity measurements (or offset removals) via
cross-correlation with model spectra, including barycentric radial velocity calculations,
(4) concatenation of multiple echelle orders into a simple 1D spectrum, and (5) approximate
flux calibration.


Installation
------------

You can install ``aesop`` from the source code by doing the following::

    git clone https://github.com/bmorris3/aesop.git
    cd aesop
    pip install .

For more information, `read the docs <https://arces.readthedocs.io/>`_.

License
-------

This project is Copyright (c) Brett Morris & Trevor Dorn-Wallenstein and licensed under
the terms of the MIT license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`__
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.
