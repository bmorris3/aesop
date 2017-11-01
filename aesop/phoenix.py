
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from itertools import product

import numpy as np
from astropy.utils.data import download_file
from astropy.io import fits
import astropy.units as u


__all__ = ['get_phoenix_model_spectrum']

phoenix_model_temps = np.array(
    [2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100,
     3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
     4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900,
     5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800,
     5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700,
     6800, 6900, 7000, 7200, 7400, 7600, 7800, 8000, 8200,
     8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000,
     10200, 10400, 10600, 10800, 11000, 11200, 11400, 11600, 11800,
     12000, 12500, 13000, 13500, 14000, 14500, 15000])

phoenix_model_metallicities = np.array([-4, -3, -2, -1.5, -1, -0.5, -0, 0.5,
                                        1.0])

phoenix_model_gravities = np.arange(0, 6.5, 0.5)


def get_solar_metallicity_url(T_eff, log_g):
    closest_grid_temperature = phoenix_model_temps[np.argmin(np.abs(phoenix_model_temps - T_eff))]

    url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
           'PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte{T_eff:05d}-{log_g:1.2f}-0.0.PHOENIX-'
           'ACES-AGSS-COND-2011-HiRes.fits').format(T_eff=closest_grid_temperature,
                                                    log_g=log_g)
    return url


def get_any_metallicity_url(T_eff, log_g, z):
    closest_grid_temperature = phoenix_model_temps[np.argmin(np.abs(phoenix_model_temps - T_eff))]

    if z > 0:
        z = "+{0:1.1f}".format(z)
    elif z == 0:
        z = "-{0:1.1f}".format(z)
    else:
        z = "{0:1.1f}".format(z)

    url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
           'PHOENIX-ACES-AGSS-COND-2011/Z{z}/lte{T_eff:05d}-'
           '{log_g:1.2f}{z}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
           ).format(T_eff=closest_grid_temperature, log_g=log_g, z=z)
    return url


def get_phoenix_model_wavelengths(cache=True):
    """
    Return the wavelength grid that the PHOENIX models were computed on,
    transformed into wavelength units in air (not vacuum).
    """
    wavelength_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/'
                      'HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
    wavelength_path = download_file(wavelength_url, cache=cache, timeout=30)
    wavelengths_vacuum = fits.getdata(wavelength_path)

    # Wavelengths are provided at vacuum wavelengths. For ground-based
    # observations convert this to wavelengths in air, as described in
    # Husser 2013, Eqns. 8-10:
    sigma_2 = (10**4 / wavelengths_vacuum)**2
    f = (1.0 + 0.05792105/(238.0185 - sigma_2) + 0.00167917 /
         (57.362 - sigma_2))
    wavelengths_air = wavelengths_vacuum / f
    return wavelengths_air


def get_phoenix_model_spectrum(T_eff, log_g=4.5, z=0, cache=True):
    """
    Download a PHOENIX model atmosphere spectrum for a star with given
    properties.

    Parameters
    ----------
    T_eff : float
        Effective temperature. The nearest grid-temperature will be selected.
    log_g : float
        This must be a log g included in the grid for the effective temperature
        nearest ``T_eff``.
    cache : bool
        Cache the result to the local astropy cache. Default is `True`.

    Returns
    -------
    spectrum : `~specutils.Spectrum1D`
        Model spectrum
    """
    url = get_any_metallicity_url(T_eff=T_eff, log_g=log_g, z=z)
    fluxes_path = download_file(url, cache=cache, timeout=30)
    fluxes = fits.getdata(fluxes_path)

    wavelengths_air = get_phoenix_model_wavelengths()

    mask_negative_wavelengths = wavelengths_air > 0

    from .spectra import Spectrum1D

    spectrum = Spectrum1D.from_array(wavelengths_air[mask_negative_wavelengths],
                                     fluxes[mask_negative_wavelengths],
                                     dispersion_unit=u.Angstrom)

    return spectrum


def download_phoenix_spectrum_grid(path, teff_min, teff_max, log_g_min,
                                   log_g_max, z_min, z_max):
    """
    Get a grid of PHOENIX model spectra.

    Parameters
    ----------
    teff_min : float
        Minimum effective temperature (inclusive)
    teff_max : float
        Maximum effective temperature (inclusive)
    log_g_min : float
        Minimum surface gravity (inclusive)
    log_g_max : float
        Maximum surface gravity (inclusive)
    z_min : float
        Minimum metallicity (inclusive)
    z_max : float
       Maximum metallicity (inclusive)

    Examples
    --------
    Currently working on using this function from the command line like this:

        >>> from aesop.phoenix import download_phoenix_spectrum_grid
        >>> download_phoenix_spectrum_grid('phoenix_grid.hdf5', 4000, 6000, 4, 5, -1, 1)
    """
    temps = ((phoenix_model_temps <= teff_max) &
             (phoenix_model_temps >= teff_min))
    gravs = ((phoenix_model_gravities <= log_g_max) &
             (phoenix_model_gravities >= log_g_min))
    metals = ((phoenix_model_metallicities <= z_max) &
              (phoenix_model_metallicities >= z_min))

    wavelengths = get_phoenix_model_wavelengths()

    import h5py
    archive = h5py.File(path, 'w')

    data_cube_shape = (len(wavelengths), np.count_nonzero(temps),
                       np.count_nonzero(gravs), np.count_nonzero(metals))
    dset = archive.create_dataset('spectra', shape=data_cube_shape,
                                  dtype=np.float64, compression='gzip')

    for i, t in enumerate(phoenix_model_temps[temps]):
        for j, g in enumerate(phoenix_model_gravities[gravs]):
            for k, z in enumerate(phoenix_model_metallicities[metals]):
                url = get_any_metallicity_url(t, g, z)
                tmp_path = download_file(url, cache=False, timeout=30)
                dset[:, i, j, k] = fits.getdata(tmp_path)
            archive.flush()

    archive.close()

