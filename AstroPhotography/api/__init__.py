"""
Application commands common to all interfaces.
"""

# This is just ugly, but it is how the cookiecutter package was set up.
from .split import main as split
from .grey import main as grey
from .rgb import main as rgb

__all__ = ["split", "grey", "rgb"]
