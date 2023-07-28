import os
import functools
import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.8f')
from glob import glob

import numpy as np

from .activity import Measurement, SIndex, StarProps

__all__ = ['glob_spectra_paths', 'stars_to_json', 'json_to_stars']


def glob_spectra_paths(data_dir, target_names):
    """
    Collect paths to spectrum FITS files.

    Parameters
    ----------
    data_dir : str or list
        Paths to the directories containing spectrum FITS files
    target_names : list
        String patterns that match the beginning of files with targets to
        collect.

    Returns
    -------
    spectra_paths : list
        List of paths to spectrum FITS files
    """
    if type(data_dir) != list:
        data_dir = [data_dir]

    all_spectra_paths = []

    for d_dir in data_dir:

        # Collect files for each target:
        spectra_paths_lists = [glob(os.path.join(d_dir,
                                                 '{0}*.wfrmcpc.fits'.format(name)))
                               for name in target_names]

        # Reduce to one list:
        spectra_paths = functools.reduce(list.__add__, spectra_paths_lists)

        all_spectra_paths.extend(spectra_paths)

    return all_spectra_paths


def combine_measurements(measurement_list):
    mean = np.mean([m.value for m in measurement_list])
    err = np.sqrt(np.sum(np.array([m.err for m in measurement_list])**2))
    return Measurement(value=mean, err=err, meta=len(measurement_list))


def floats_to_strings(d):
    dictionary = d.copy()
    for key in dictionary:
        dictionary[key] = str(dictionary[key])
    return dictionary


def stars_to_json(star_list, output_path='star_data.json'):
    """
    Save list of stellar properties to a JSON file.

    Parameters
    ----------
    star_list : list of `StarProps`
        Star properties to save to json
    output_path : str
        File path to output
    """
    stars_attrs = star_list[0].__dict__.keys()
    all_data = dict()

    for star in star_list:
        star_data = dict()

        for attr in stars_attrs:
            value = getattr(star, attr)

            if isinstance(value, Measurement):
                value = floats_to_strings(value.__dict__)
            elif isinstance(value, SIndex):
                value = value.to_dict()
            else:
                value = str(value)

            star_data[attr] = value

        all_data[star.name + '; ' + str(star.time.datetime)] = star_data

    with open(output_path, 'w') as w:
        json.dump(all_data, w, indent=4, sort_keys=True)


def json_to_stars(json_path):
    """
    Loads JSON archive into list of `StarProps` objects.

    Parameters
    ----------
    json_path : str
        Path to saved stellar properties

    Returns
    -------
    stars : list of `StarProps`
        List of stellar properties.
    """
    with open(json_path, 'r') as w:
        dictionary = json.load(w)

    stars = [StarProps.from_dict(dictionary[star]) for star in dictionary]
    return stars

