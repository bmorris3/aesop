
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import h5py
from astropy.utils.data import download_file
from astropy.io import fits
import astropy.units as u

__all__ = ['get_phoenix_model_spectrum', 'PHOENIXModelGrid']

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
                                   log_g_max, z_min, z_max, wavelength_min=3000,
                                   wavelength_max=10000):
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
    wavelength_min : float
        Minimum wavelength to save (angstrom)
    wavelength_max : float
        Maximum wavelength to save (angstrom)

    Examples
    --------
    Currently working on using this function from the command line like this:

        >>> from aesop.phoenix import download_phoenix_spectrum_grid
        >>> download_phoenix_spectrum_grid('phoenix_grid.hdf5', 4000, 6000, 4, 5, -1, 1)

    The HDF5 archive that gets saved is roughly 0.5 GB.
    """
    temps = ((phoenix_model_temps <= teff_max) &
             (phoenix_model_temps >= teff_min))
    gravs = ((phoenix_model_gravities <= log_g_max) &
             (phoenix_model_gravities >= log_g_min))
    metals = ((phoenix_model_metallicities <= z_max) &
              (phoenix_model_metallicities >= z_min))

    wavelengths = get_phoenix_model_wavelengths()

    wavelength_mask = ((wavelengths < wavelength_max) &
                       (wavelengths > wavelength_min))
    wavelengths_in_bounds = wavelengths[wavelength_mask]

    archive = h5py.File(path, 'w')

    data_cube_shape = (len(wavelengths_in_bounds), np.count_nonzero(temps),
                       np.count_nonzero(gravs), np.count_nonzero(metals))
    dset = archive.create_dataset('spectra', shape=data_cube_shape,
                                  dtype=np.float32, compression='gzip')

    dset.attrs['temperatures'] = phoenix_model_temps[temps]
    dset.attrs['gravities'] = phoenix_model_gravities[gravs]
    dset.attrs['metallicities'] = phoenix_model_metallicities[metals]
    #dset.attrs['wavelengths'] = wavelengths_in_bounds

    for i, t in enumerate(phoenix_model_temps[temps]):
        for j, g in enumerate(phoenix_model_gravities[gravs]):
            for k, z in enumerate(phoenix_model_metallicities[metals]):
                if np.all(dset[:, i, j, k] == 0):
                    url = get_any_metallicity_url(t, g, z)
                    tmp_path = download_file(url, cache=False, timeout=30)
                    spectrum = fits.getdata(tmp_path)[wavelength_mask]
                    dset[:, i, j, k] = spectrum
            archive.flush()

    archive.close()

default_archive_path = 'phoenix_grid.hdf5'
wavelength_min = 3000  # Angstrom
wavelength_max = 10000  # Angstrom


class PHOENIXModelGrid(object):
    def __init__(self, path=default_archive_path,
                 interp_wavelength_min=6562.8-5,
                 interp_wavelength_max=6562.8+5):
        self.path = path

        if not os.path.exists(path):
            raise ValueError('No such HDF5 archive {0}'.format(path))

        archive = h5py.File(path, 'r')
        dset = archive['spectra']
        self.temperatures = dset.attrs['temperatures']
        self.gravities = dset.attrs['gravities']
        self.metallicities = dset.attrs['metallicities']

        all_wavelengths = get_phoenix_model_wavelengths()
        wavelength_mask = ((all_wavelengths < wavelength_max) &
                           (all_wavelengths > wavelength_min))
        wavelengths_in_bounds = all_wavelengths[wavelength_mask]
        interp_bounds = ((wavelengths_in_bounds < interp_wavelength_max) &
                         (wavelengths_in_bounds > interp_wavelength_min))
        self.wavelengths = wavelengths_in_bounds[interp_bounds]


        points = (self.wavelengths, self.temperatures, self.gravities,
                  self.metallicities)
        values = dset[np.where(interp_bounds)[0], :, :, :][:]

        rgi = RegularGridInterpolator(points, values)
        self._rgi = rgi

    def interp(self, temperature, gravity, metallicity, wavelengths=None,
               method='linear'):
        if wavelengths is None:
            wavelengths = self.wavelengths
        xi = np.hstack([wavelengths[:, np.newaxis],
                        np.repeat([[temperature, gravity, metallicity]],
                                  len(wavelengths), axis=0)])
        return self._rgi(xi, method=method)

    @classmethod
    def download_phoenix_spectrum_grid(cls, teff_min, teff_max, log_g_min,
                                       log_g_max, z_min, z_max,
                                       wavelength_min=wavelength_min,
                                       wavelength_max=wavelength_max,
                                       path=default_archive_path):
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
        wavelength_min : float
            Minimum wavelength to save (angstrom)
        wavelength_max : float
            Maximum wavelength to save (angstrom)

        Notes
        -----
        The HDF5 archive that gets saved is roughly 0.5 GB.
        """
        temps = ((phoenix_model_temps <= teff_max) &
                 (phoenix_model_temps >= teff_min))
        gravs = ((phoenix_model_gravities <= log_g_max) &
                 (phoenix_model_gravities >= log_g_min))
        metals = ((phoenix_model_metallicities <= z_max) &
                  (phoenix_model_metallicities >= z_min))

        wavelengths = get_phoenix_model_wavelengths()

        wavelength_mask = ((wavelengths < wavelength_max) &
                           (wavelengths > wavelength_min))
        wavelengths_in_bounds = wavelengths[wavelength_mask]

        import h5py
        archive = h5py.File(path, 'w')

        data_cube_shape = (len(wavelengths_in_bounds), np.count_nonzero(temps),
                           np.count_nonzero(gravs), np.count_nonzero(metals))
        dset = archive.create_dataset('spectra', shape=data_cube_shape,
                                      dtype=np.float32, compression='gzip')

        dset.attrs['temperatures'] = phoenix_model_temps[temps]
        dset.attrs['gravities'] = phoenix_model_gravities[gravs]
        dset.attrs['metallicities'] = phoenix_model_metallicities[metals]

        for i, t in enumerate(phoenix_model_temps[temps]):
            for j, g in enumerate(phoenix_model_gravities[gravs]):
                for k, z in enumerate(phoenix_model_metallicities[metals]):
                    if np.all(dset[:, i, j, k] == 0):
                        url = get_any_metallicity_url(t, g, z)
                        tmp_path = download_file(url, cache=False, timeout=30)
                        spectrum = fits.getdata(tmp_path)[wavelength_mask]
                        dset[:, i, j, k] = spectrum
                archive.flush()

        archive.close()
        return cls(path=path)

    def spectrum(self, temperature, gravity, metallicity, wavelengths=None,
                 method='linear'):
        if wavelengths is None:
            wavelengths = self.wavelengths
        flux = self.interp(temperature, gravity, metallicity, method=method,
                           wavelengths=wavelengths)
        from .spectra import Spectrum1D
        return Spectrum1D(wavelengths if hasattr(wavelengths, 'unit') else
                          u.Quantity(wavelengths, u.Angstrom), flux)
