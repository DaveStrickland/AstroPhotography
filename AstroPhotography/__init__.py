"""
The AstroPhotography Python Package
------------------------------------

The AstroPhotography python package provides python classes and command line
applications for quick and convenient processing amateur astronomical images
from RAW digital camera format to common graphical image formats (e.g. PNG)
and astronomical FITS format, and for calibrating and combining series
of FITS images.
"""

# SPDX-License-Identifier: GPL-3.0-or-later

from .__version__ import __version__
from .__main__ import warn
from .core import *
from .api import *
from .util import *
