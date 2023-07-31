try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

from .phoenix import *
from .spectra import *
from .utils import *
from .masking import *
from .activity import *
from .spectral_type import *
