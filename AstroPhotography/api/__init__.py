"""
AstroPhotography api module
---------------------------

The ``api`` module contains functins and classes associated with camera
RAW format conversion and processing, e.g. through the ``dksraw`` command
line utility.
"""

from .split import main as split
from .grey import main as grey
from .rgb import main as rgb

__all__ = ["split", "grey", "rgb"]
