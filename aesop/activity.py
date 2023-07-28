"""
Tools for measuring equivalent widths, S-indices.
"""

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.time import Time

from .catalog import query_catalog_for_object

__all__ = ['integrate_spectrum_trapz', 'true_h_centroid', 'true_k_centroid',
           'uncalibrated_s_index', 'StarProps', 'Measurement',
           'FitParameter']

true_h_centroid = 3968.4673 * u.Angstrom
true_k_centroid = 3933.6614 * u.Angstrom


def integrate_spectrum_trapz(spectrum, center_wavelength, width,
                             weighting=False, plot=False):
    """
    Integrate the area under a spectrum.

    Parameters
    ----------
    spectrum : `EchelleSpectrum`
        Spectrum to integrate under
    center_wavelength : `~astropy.units.Quantity`
        Center of region to integrate
    width : `~astropy.units.Quantity`
        Width about the center to integrate
    wavelength_offset : float
        Offset wavelengths by this amount, which is useful if a wavelength
        solution refinement as been made.
    weighting : bool
        Apply a triangular weighting function to the fluxes

    Returns
    -------
    integral : float
        Integral under the spectrum
    error : float
        Square-root of the sum of the fluxes within the bounds of the integral
    """
    wavelength = spectrum.wavelength
    wavelengths_increasing = wavelength[1] > wavelength[0]
    flux = spectrum.flux
    norm_const = spectrum.meta['normalization']

    if not wavelengths_increasing:
        wavelength = wavelength[::-1]
        flux = flux[::-1]

    # Assert fluxes are not negative
    flux.value[flux.value < 0] = 0.0

    if (not center_wavelength < wavelength.max() and
            not center_wavelength > wavelength.min()):
        raise ValueError("This spectral order does not contain"
                         "the center_wavelength given.")

    within_bounds = ((wavelength > center_wavelength - width/2) &
                     (wavelength < center_wavelength + width/2))

    if not weighting:
        integral = np.trapz(flux[within_bounds].value,
                            wavelength[within_bounds].value)

        lam = wavelength[within_bounds].value
        sigma_f = (np.sqrt(flux[within_bounds].value *
                           norm_const[within_bounds]) /
                   norm_const[within_bounds])
        sigma_f[np.isnan(sigma_f)] = np.nanmean(sigma_f)  # Fix negative fluxes
        a2_k = 0.25 * (lam[1:] - lam[:-1])**2 * (sigma_f[1:]**2 +
                                                 sigma_f[:-1]**2)
        error = np.sqrt(np.sum(a2_k))

    else:
        triangle_weights = triangle_weighting(wavelength[within_bounds],
                                              center_wavelength)
        integral = np.trapz(flux[within_bounds].value * triangle_weights,
                            wavelength[within_bounds].value)

        lam = wavelength[within_bounds].value
        sigma_f = (np.sqrt(flux[within_bounds].value *
                           norm_const[within_bounds]) * triangle_weights /
                   norm_const[within_bounds])
        sigma_f[np.isnan(sigma_f)] = np.nanmean(sigma_f)  # Fix negative fluxes
        a2_k = 0.25 * (lam[1:] - lam[:-1])**2 * (sigma_f[1:]**2 +
                                                 sigma_f[:-1]**2)
        error = np.sqrt(np.sum(a2_k))

    if plot:
        plt.figure()
        plt.plot(wavelength[spectrum.mask], flux[spectrum.mask])
        if weighting:
            triangle = triangle_weighting(wavelength[within_bounds],
                                          center_wavelength)
            plt.plot(wavelength[within_bounds], triangle, 'r', lw=2)
        plt.show()

    return Measurement(integral, err=error)


def triangle_weighting(x, x0, fwhm=1.09*u.Angstrom):
    """
    Compute the triangular weighting function used in CaII cores

    Parameters
    ----------
    x : `~astropy.units.Quantity`, array-like
        Wavelengths
    x0 : `~astropy.units.Quantity`
        Central wavelength of emission feature
    fwhm : `~astropy.units.Quantity`
        Full-width at half maximum of the weighting function
    """

    left_half = (x <= x0) & (x > x0 - fwhm)
    right_half = (x > x0) & (x < x0 + fwhm)

    weights = np.zeros_like(x.value)
    # float typecasting below resolves the unitless quantities
    weights[left_half] = ((x[left_half] - x0) / fwhm).value + 1
    weights[right_half] = ((x0 - x[right_half]) / fwhm).value + 1

    return weights


def uncalibrated_s_index(spectrum, plots=False):
    """
    Calculate the uncalibrated S-index from an Echelle spectrum.

    Parameters
    ----------
    spectrum : `EchelleSpectrum`
        Normalized target spectrum

    Returns
    -------
    s_ind : `SIndex`
        S-index. This value is intrinsic to the instrument you're using.
    """

    order_h = spectrum.get_order(89)
    order_k = spectrum.get_order(90)
    order_r = spectrum.get_order(91)
    order_v = spectrum.get_order(88)

    r_centroid = 3900 * u.Angstrom
    v_centroid = 4000 * u.Angstrom
    hk_fwhm = 1.09 * u.Angstrom
    hk_width = 2 * hk_fwhm
    rv_width = 20 * u.Angstrom

    h = integrate_spectrum_trapz(order_h, true_h_centroid, hk_width,
                                        weighting=True, plot=plots)
    k = integrate_spectrum_trapz(order_k, true_k_centroid, hk_width,
                                        weighting=True, plot=plots)
    r = integrate_spectrum_trapz(order_r, r_centroid, rv_width, plot=plots)
    v = integrate_spectrum_trapz(order_v, v_centroid, rv_width, plot=plots)

    s_ind = SIndex(h=h, k=k, r=r, v=v, time=spectrum.time)
    return s_ind


class SIndex(object):
    def __init__(self, h, k, r, v, k_factor=0.84, v_factor=1.0, time=None):
        """
        The pre-factors have been chosen to make the ``h`` and ``k`` values
        of the same order of magnitude; same for ``r`` and ``v``.

        Parameters
        -----------
        h : float
            CaII H feature emission flux
        k : float
            CaII K feature emission flux
        r : float
            Pseudo-continuum flux redward of CaII H&K
        v : float
            Pseudo-continuum flux blueward of CaII H&K
        k_factor : float
            Multiplicative factor for the K emission feature flux
            to make the H & K fluxes similar
        v_factor : float
            Multiplicative factor for the blue continuum region flux
            to make the r & v fluxes similar
        time : `~astropy.time.Time`
            Time this S-index measurement was taken.
        """
        self.r = r
        self.v = v
        self.h = h
        self.k = k

        self.k_factor = k_factor
        self.v_factor = v_factor

        self.time = time

    @property
    def uncalibrated(self):
        """
        Compute Eqn 2 of Isaacson+ 2010, for C1=1 and C2=0. This can be used
        to solve for C1 and C2.
        """
        uncalibrated_s_ind = ((self.h.value + self.k_factor * self.k.value) /
                              (self.r.value + self.v_factor * self.v.value))

        s_ind_err = (1 / (self.r.value + self.v_factor*self.v.value)**2 *
                     (self.h.err**2 + self.k_factor**2 * self.k.err**2) +
                     (self.h.value + self.k_factor*self.k.value)**2 /
                     (self.r.value + self.v_factor * self.v.value)**4 *
                     (self.r.err**2 + self.v_factor**2 * self.v.err**2)
                     )**0.5

        return Measurement(uncalibrated_s_ind, err=s_ind_err)

    def calibrated(self, c1, c2):
        """
        Calibrated S-index measurement (comparable with MWO S-indices).

        Uses the scaling constants as defined in Isaacson 2010+ (c1 and c2).

        Parameters
        ----------
        c1 : float
        c2 : float

        Returns
        -------
        Calibrated S-index.
        """
        return c1 * self.uncalibrated + c2

    @classmethod
    def from_dict(cls, dictionary):
        d = dictionary.copy()
        for key in dictionary:
            if key == 'time':
                d[key] = Time(float(d[key]), format='jd')
            elif key in ['h', 'k', 'r', 'v']:
                d[key] = Measurement.from_dict(d[key])
            else:
                d[key] = float(d[key])
        return cls(**d)

    def to_dict(self):
        d = dict()

        for attr in self.__dict__:
            value = getattr(self, attr)
            if isinstance(value, Measurement):
                value = value.__dict__
            elif isinstance(value, Time):
                value = str(value.jd)
            d[attr] = value
        return d


class StarProps(object):
    """
    S-index properties for a star
    """
    def __init__(self, name=None, s_apo=None, s_mwo=None, time=None):
        self.name = name
        self.s_apo = s_apo
        self._s_mwo = s_mwo
        self.time = time

    def get_s_mwo(self):
        obj = query_catalog_for_object(self.name)

        # Replace the uncertainty by twice the mean uncertainty
        # of the Duncan 1991 tables if no uncertainty is provided
        # (calculated by get_mean_mwo_error.py)
        error_Smean = obj['e_Smean'] if obj['e_Smean'] != 0 else 10 * 0.0205

        self._s_mwo = Measurement(obj['Smean'], err=error_Smean)

    @property
    def s_mwo(self):
        if self._s_mwo is None:
            self.get_s_mwo()
        return self._s_mwo

    @classmethod
    def from_dict(cls, dictionary):
        if dictionary['s_apo'] != 'None':
            if "value" in dictionary['s_apo']:
                s_apo = Measurement.from_dict(dictionary['s_apo'])
            else:
                s_apo = SIndex.from_dict(dictionary['s_apo'])
        else:
            s_apo = None

        if '_s_mwo' in dictionary and dictionary['_s_mwo'] != "None":
            s_mwo = Measurement.from_dict(dictionary['_s_mwo'])
        else:
            s_mwo = None

        if dictionary['time'] != 'None':
            dictionary['time'] = Time(float(dictionary['time']), format='jd')
        else:
            dictionary['time'] = None

        return cls(s_apo=s_apo, s_mwo=s_mwo, name=dictionary['name'],
                   time=dictionary['time'])


class Measurement(object):
    def __init__(self, value=None, err=None, time=None, default_err=1e10,
                 meta=None):

        if hasattr(value, '__len__'):
            self.value = np.asarray(value)
            self.err = np.asarray(err)

            self.err[self.err == 0] = default_err

        else:
            self.value = value
            if err == 0:
                self.err = default_err
            else:
                self.err = err

        if isinstance(time, Time):
            self.time = time
        elif hasattr(time, 'real'):  # if time is float
            self.time = Time(time, format='jd')
        elif isinstance(time, list):
            if hasattr(time[0], 'real'):
                self.time = Time(time, format='jd')
            else:
                self.time = Time(time)
        else:
            self.time = None

        self.meta = meta

    @classmethod
    def from_min_max(cls, min, max):
        mean = 0.5*(min + max)
        return cls(value=mean, err=max-mean)

    @classmethod
    def from_dict(cls, dictionary):
        kwargs = {key: float(dictionary[key])
                  if (dictionary[key] != "None" and
                      dictionary[key] is not None)
                  else None
                  for key in dictionary}
        return cls(**kwargs)

    def __repr__(self):
        return "<{0}: {1} +/- {2}>".format(self.__class__.__name__,
                                          self.value, self.err)

    def __getitem__(self, item):

        new_attrs = dict()

        for attr in ['value', 'err', 'time']:
            attr_value = getattr(self, attr)
            if attr_value is not None and hasattr(attr_value, '__len__'):
                new_attrs[attr] = attr_value[item]

        return Measurement(**new_attrs)

    def __len__(self):
        return len(self.value)

    def to_latex(self):
        return "${0:.3f} \pm {1}$".format(self.value, error_to_latex(self.err))


def error_to_latex(error):
    str_err = "{0:.2g}".format(error)
    if 'e' in str_err:
        str_err = "{0:.6f}".format(error)
    return str_err


class FitParameter(object):
    """
    Fit results for a fitting parameter, with asymmetrical errors.
    """
    def __init__(self, value, err_upper=None, err_lower=None, default_err=1e10):
        """

        Parameters
        ----------
        value : {`~astropy.units.Quantity`, `~numpy.ndarray`, float}
        err_upper : {`~astropy.units.Quantity`, `~numpy.ndarray`, float}
        err_lower : {`~astropy.units.Quantity`, `~numpy.ndarray`, float}
        default_err : float
        """
        if hasattr(value, '__len__'):
            value = np.asarray(value)
            err_upper = np.asarray(err_upper)
            err_lower = np.asarray(err_lower)

            self.value = value
            self.err_upper = err_upper
            self.err_lower = err_lower

        else:

            self.value = value
            if err_upper == 0 or err_lower == 0:
                self.err_upper = default_err
                self.err_lower = default_err
            else:
                self.err_upper = err_upper
                self.err_lower = err_lower

    @classmethod
    def from_text(cls, path):
        value, err_upper, err_lower = np.loadtxt(path, unpack=True)
        return cls(value, err_lower=err_lower, err_upper=err_upper)

    def to_text(self, path):
        np.savetxt(path, [self.value, self.err_upper, self.err_lower])

    def __repr__(self):
        return "<{0}: {1} +{2} -{3}>".format(self.__class__.__name__,
                                             self.value, self.err_upper,
                                             self.err_lower)
