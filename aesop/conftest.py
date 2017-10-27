from astropy.tests.pytest_plugins import *

# also save a copy of the astropy hooks so we can use them below when
# overriding
from astropy.tests import pytest_plugins as astropy_pytest_plugins

import warnings

import os

# This is to figure out the affiliated package version, rather than
# using Astropy's
try:
    from .version import version
except ImportError:
    version = 'dev'

packagename = os.path.basename(os.path.dirname(__file__))
TESTED_VERSIONS[packagename] = version


# Comment out this line to avoid deprecation warnings being raised as
# exceptions
enable_deprecations_as_exceptions()

# Define list of packages for which to display version numbers in the test log
try:
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    del PYTEST_HEADER_MODULES['h5py']
except KeyError:
    pass


def pytest_configure(config):
    if hasattr(astropy_pytest_plugins, 'pytest_configure'):
        # sure ought to be true right now, but always possible it will change in
        # future versions of astropy
        astropy_pytest_plugins.pytest_configure(config)

    # activate image comparison tests only if the dependencies needed are installed:
    # matplotlib, nose, pytest-mpl
    try:
        import matplotlib
        import nose  # needed for the matplotlib testing tools
        HAS_MATPLOTLIB_AND_NOSE = True
    except ImportError:
        HAS_MATPLOTLIB_AND_NOSE = False

    if HAS_MATPLOTLIB_AND_NOSE and config.pluginmanager.hasplugin('mpl'):
        pass
        # TODO: turn image comparison tests back on once this issue is figured out:
        # https://github.com/astropy/astroplan/pull/104#issuecomment-137734007
        # config.option.mpl = True
        # config.option.mpl_baseline_path = 'astroplan/plots/tests/baseline_images'

