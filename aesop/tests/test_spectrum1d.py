from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import numpy as np
import astropy.units as u
import pytest
from astropy.units import UnitsError

from ..spectra import Spectrum1D


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
