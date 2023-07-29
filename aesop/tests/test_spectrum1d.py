

import numpy as np
import astropy.units as u
import pytest
from astropy.units import UnitsError
from astropy.utils.data import download_file

from ..spectra import Spectrum1D, EchelleSpectrum


def test_constructor():
    wl = np.linspace(3000, 4000) * u.Angstrom
    f = 1 + 0.1 * np.random.randn()
    example = Spectrum1D(wavelength=wl, flux=f)

    # raise an error if there is no unit on the wavelength
    with pytest.raises(TypeError):
        example = Spectrum1D(wavelength=[0, 1, 2, 3], flux=f)

    # raise an error if the wavelength unit is not a unit of length
    with pytest.raises(UnitsError):
        example = Spectrum1D(wavelength=u.Quantity([0, 1, 2, 3], unit=u.day),
                             flux=f)


@pytest.mark.remote_data
def test_read_fits():
    url = ('https://stsci.box.com/shared/static/'
           'mu4fa1fmq1lw8boem12e2umyi99skbdl.fits')

    path = download_file(url, show_progress=False)

    echelle_spectrum = EchelleSpectrum.from_fits(path, format='iraf')

    assert hasattr(echelle_spectrum, 'header')
    assert hasattr(echelle_spectrum, 'time')
    assert hasattr(echelle_spectrum, 'name')
    assert str(echelle_spectrum) == ('<EchelleSpectrum: 107 orders, '
                                     '3506.8-10612.4 Angstrom>')

    # There should be more flux in the 40th order (redder) than 0th (bluer)
    assert echelle_spectrum[40].flux.mean() > echelle_spectrum[0].flux.mean()


def generate_target_standard_pairs():
    # realistic wavelength range, resolution for one order
    wl = np.linspace(6000, 6117, 1650) * u.Angstrom
    sig = 25 + 10 * (0.5 - np.random.rand())
    amp = 500 + 200 * (0.5 - np.random.rand())
    flux = amp * np.exp(-0.5 * (wl.value - wl.value.mean())**2 / sig**2)
    flux += 10 * np.random.randn(flux.shape[0])

    standard = Spectrum1D(wavelength=wl, flux=flux)

    for i in range(4):
        flux -= 100 * np.exp(-0.5 * (wl.value - wl.value.mean() -
                                     20 * (i - 1.5))**2 / 0.7**2)

    target = Spectrum1D(wavelength=wl, flux=flux)

    return target, standard


def test_continuum_norm():
    target_orders = []
    standard_orders = []
    for i in range(10):
        target, standard = generate_target_standard_pairs()
        target_orders.append(target)
        standard_orders.append(standard)

    target_spectrum = EchelleSpectrum(target_orders)
    standard_spectrum = EchelleSpectrum(standard_orders)

    # Make sure that the continuum normalization is always good to 2%
    polynomial_order = 8
    target_spectrum.continuum_normalize_from_standard(standard_spectrum,
                                                      polynomial_order)

    for order in target_spectrum.spectrum_list:
        assert abs(np.median(order.flux) - 1) < 0.02

    # Make sure spectral features are preserved in target spectra
    for order in target_spectrum.spectrum_list:
        assert np.mean(order.flux) < 1