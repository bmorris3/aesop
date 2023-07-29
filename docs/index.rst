
.. _aesop:

.. include:: references.txt

*****
aesop
*****

``aesop`` is a data reduction toolkit for echelle spectra from the `Astrophysical
Research Consortium (ARC) 3.5 m Telescope Echelle Spectrograph (ARCES)
<http://www.apo.nmsu.edu/arc35m/Instruments/ARCES/>`__ at Apache Point
Observatory (APO). You can view the `source code and submit issues via GitHub
<https://github.com/bmorris3/aesop>`_.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   install
   iraf
   getting_started
   api


What is ``aesop``?
==================

In the interest of reproducing established results, we choose not to
re-implement the canonical IRAF reduction routine
`ReduceARCES.cl <https://github.com/bmorris3/aesop/blob/master/data/ReduceARCES.cl>`_,
and to limit ``aesop``'s functionality to processing after the aperture
extraction and initial wavelength solution.

The current version of ``aesop`` is useful for the following tasks:

 - Managing and manipulating spectra with convenient data structures
   (`~aesop.EchelleSpectrum`, `~aesop.Spectrum1D`)

 - Removing radial velocities from echelle orders via cross-correlation with
   PHOENIX model spectra

 - Normalizing out the blaze-function in each echelle order given observations
   of a spectroscopic standard

 - Further normalizing each order for a flat continuum with robust least-squares

 - Concatenating the spectra in each order, and taking the mean between adjacent
   orders where they overlap, to produce a 1D spectrum from 3000-10000 Angstroms


Typical Workflow
================

In general, ``aesop`` users will follow this procedure:

 1. Install ``aesop``
 2. Run the ``ReduceARCES.cl`` IRAF script to reduce your echelle spectra
    (see :ref:`iraf`)
 3. Open the extracted 1D spectra for each spectral order with ``aesop`` (see
    :ref:`getting_started`)


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
