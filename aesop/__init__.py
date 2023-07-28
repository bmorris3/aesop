from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    __version__ = 'dev'
    pass  # pragma: no cover

del get_distribution, DistributionNotFound

from .phoenix import *
from .spectra import *
from .utils import *
from .masking import *
from .activity import *
from .spectral_type import *
