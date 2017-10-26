from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .catalog import get_k2_epic_catalog

from astroquery.simbad import Simbad

__all__ = ['query_for_spectral_type', 'query_for_T_eff']

# Source: http://www.uni.edu/morgans/astro/course/Notes/section2/spectraltemps.html
effective_temperatures = {'A0': 9600,
 'A1': 9330,
 'A2': 9040,
 'A3': 8750,
 'A4': 8480,
 'A5': 8310,
 'A7': 7920,
 'B0': 29200,
 'B1': 23000,
 'B2': 21000,
 'B3': 17600,
 'B5': 15200,
 'B6': 14300,
 'B7': 13500,
 'B8': 12300,
 'B9': 11400,
 'F0': 7350,
 'F2': 7050,
 'F3': 6850,
 'F5': 6700,
 'F6': 6550,
 'F7': 6400,
 'F8': 6300,
 'G0': 6050,
 'G1': 5930,
 'G2': 5800,
 'G5': 5660,
 'G8': 5440,
 'K0': 5240,
 'K1': 5110,
 'K2': 4960,
 'K3': 4800,
 'K4': 4600,
 'K5': 4400,
 'K7': 4000,
 'L0': 2600,
 'L3': 2200,
 'L8': 1500,
 'M0': 3750,
 'M1': 3700,
 'M2': 3600,
 'M3': 3500,
 'M4': 3400,
 'M5': 3200,
 'M6': 3100,
 'M7': 2900,
 'M8': 2700,
 'O5': 54000,
 'O6': 45000,
 'O7': 43300,
 'O8': 40600,
 'O9': 37800,
 'T2': 1400,
 'T6': 1000,
 'T8': 800}


def query_for_spectral_type(identifier, only_first_two_characters=True,
                            default_sptype='G0'):
    """
    Search SIMBAD for the spectral type of a star.

    If no spectral type is found, the default return value is ``"G0"``.

    Parameters
    ----------
    identifier : str
        Name of target
    only_first_two_characters : bool
        Return only first two characters of spectral type?
    default_sptype : str
        Spectral type returned when none is found on SIMBAD

    Returns
    -------
    sptype : str
        Spectral type of the star.
    """
    customSimbad = Simbad()
    customSimbad.SIMBAD_URL = 'http://simbad.harvard.edu/simbad/sim-script'
    customSimbad.add_votable_fields('sptype')
    result = customSimbad.query_object(identifier)

    if len(result['SP_TYPE']) > 0:

        if only_first_two_characters:
            return result['SP_TYPE'][0][:2].strip().decode()
        else:
            return result['SP_TYPE'][0].decode()
    else:
        return default_sptype.decode()


def query_for_T_eff(identifier):
    """
    Get the approximate effective temperature of a star.

    Query SIMBAD for the spectral type of the target, convert
    spectral type to approximate effective temperature, in general.
    If the target is in the K2 EPIC, use the EPIC Teff.

    Parameters
    ----------
    identifier : str
        Name of target

    Returns
    -------
    T_eff : int
        Approximate effective temperature of the star.
    """
    if not identifier.startswith('EPIC'):
        sptype = query_for_spectral_type(identifier)
        while not sptype in effective_temperatures:
            letter, number = list(sptype)
            sptype = letter + str(int(number) - 1)

        T_eff = effective_temperatures[sptype]
    else:
        k2_epic_table = get_k2_epic_catalog()
        epic_number = int(identifier[4:]) # Remove the EPIC, make int
        T_eff = k2_epic_table.loc[epic_number]['Teff']

    return T_eff
