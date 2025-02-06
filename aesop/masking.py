import numpy as np
from scipy.optimize import fmin_l_bfgs_b
import matplotlib.pyplot as plt

__all__ = ['get_spectrum_mask']

def _gaussian(x, a, x0, sigma):
    return a * np.exp(-0.5 * (x - x0)**2 / sigma**2)


def _chi2(p, x, y):
    return np.sum((_gaussian(x, *p) - y) ** 2)


def get_spectrum_mask(spectrum, cutoff=1.5, plot=False):
    """
    Fit the raw, pre-normalized spectrum's blaze function
    with a Gaussian. Use a

    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        Spectrum to fit
    cutoff : float
        Mask channels greater than ``cuttoff``-sigma away from
        the peak flux.

    Returns
    -------
    mask : `~numpy.ndarray`
        Mask that excludes channels greater than ``cuttoff``-sigma away from
        the peak flux.
    """
    wavelength_ptp = np.ptp(spectrum.wavelength.value)
    initp = np.array([spectrum.flux.max().value,
                      spectrum.wavelength.mean().value,
                      wavelength_ptp/4])

    bestp = fmin_l_bfgs_b(_chi2, initp[:],
                          approx_grad=True,
                          args=(spectrum.wavelength.value,
                                spectrum.flux.value),
                          bounds=[(0, np.inf),
                                  (spectrum.wavelength.min().value,
                                   spectrum.wavelength.max().value),
                                  (wavelength_ptp/8,
                                   wavelength_ptp/2)])[0]

    best_a, best_x0, best_sigma = bestp

    mask = ((spectrum.wavelength.value > best_x0 - cutoff * best_sigma) &
            (spectrum.wavelength.value < best_x0 + cutoff * best_sigma))

    if plot:
        plt.figure()
        plt.plot(spectrum.wavelength, spectrum.flux, label='unmasked')
        plt.plot(spectrum.wavelength[mask], spectrum.flux[mask], label='masked')
        plt.legend()

    return np.logical_not(mask)
