.. _getting_started:

.. include:: references.txt

***************
Getting started
***************

Contents
========

* :ref:`getting_started-example`
* :ref:`getting_started-norm`

.. _getting_started-example:

Opening example spectra
-----------------------

Example ARCES spectra are available online for you to work with. You can
download them via Python like this:

.. code-block:: python

    >>> from astropy.utils.data import download_file

    >>> target_url = 'http://staff.washington.edu/bmmorris/docs/KIC8462852.0065.wfrmcpc.fits'
    >>> spectroscopic_standard_url = 'http://staff.washington.edu/bmmorris/docs/BD28_4211.0034.wfrmcpc.fits'

    >>> target_path = download_file(target_url, show_progress=False) # doctest: +REMOTE_DATA
    >>> standard_path = download_file(spectroscopic_standard_url, show_progress=False) # doctest: +REMOTE_DATA

This will download temporary copies of a Kepler target, KIC 8462852, and
spectroscopic standard O star BD+28 4211. We first create an
`~aesop.EchelleSpectrum` object for each star:

.. code-block:: python

    >>> from aesop import EchelleSpectrum

    >>> target_spectrum = EchelleSpectrum.from_fits(target_path)
    >>> standard_spectrum = EchelleSpectrum.from_fits(standard_path)

You can check basic metadata for an `~aesop.EchelleSpectrum` object by printing
it:

.. code-block:: python

    >>> print(target_spectrum)
    <EchelleSpectrum: 107 orders, 3506.8-10612.5 Angstrom>

The `~aesop.EchelleSpectrum` object behaves a bit like a Python list -- it
supports indexing, where the index counts the order number, starting with index
zero for the order with the shortest wavelengths. Elements of the
`~aesop.EchelleSpectrum` are `~aesop.Spectrum1D` objects. Suppose you want to
make a quick plot of the 73rd order of the target's echelle spectrum:

.. code-block:: python

    order73 = target_spectrum[73]
    order73.plot()

.. plot::

    from astropy.utils.data import download_file

    target_url = 'http://staff.washington.edu/bmmorris/docs/KIC8462852.0065.wfrmcpc.fits'
    spectroscopic_standard_url = 'http://staff.washington.edu/bmmorris/docs/BD28_4211.0034.wfrmcpc.fits'

    target_path = download_file(target_url)
    standard_path = download_file(spectroscopic_standard_url)

    from aesop import EchelleSpectrum

    target_spectrum = EchelleSpectrum.from_fits(target_path)
    standard_spectrum = EchelleSpectrum.from_fits(standard_path)

    order73 = target_spectrum[73]
    order73.plot()

    import matplotlib.pyplot as plt
    plt.show()

You can see the H-alpha absorption in this O star at 6562 Angstroms.

.. _getting_started-norm:

Normalizing your spectra
------------------------

You can continuum-normalize your echelle spectra with the...
