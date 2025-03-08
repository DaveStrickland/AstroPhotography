"""
File checking utilities.

Note that this module should NOT depend on any other AstroPhotography modules.
"""

from ..core.logger import logger
#import imageio
#import time
import os.path
#import yaml
from pathlib import Path
#from datetime import datetime, timezone
#from astropy.io import fits
from typing import Any
#import numpy as np
#import numpy.typing as npt

# AstroPhotography includes
from .. import __version__

def does_file_exist(filename: str, verbose: bool = False) -> bool:
    """
    Returns True if the file name or path exists, false otherwise

    Parameters
    ----------
    filename : str
        Name (optionally including path) of the file we want to check
        the existence of.
    verbose : bool, optional, default=False
        If True then writes to stdout.

    Returns
    -------
    exists : bool
        Returns True if the file path exists, False otherwise.
    """
    if verbose and not Path(filename).exists():
        print(f"Cannot find {filename}. Not a valid path or file.")
        return False
    else:
        print(f"Found {filename}.")
    return True


def _does_file_exist(filename: str, verbose: bool = False, a_logger: Any = None) -> bool:
    """
    Returns True if the file name or path exists, false otherwise

    Parameters
    ----------
    filename : str
        Name (optionally including path) of the file we want to check
        the existence of.
    verbose : bool, optional, default=False
        If True then writes to stdout.
    a_logger : Logging instance or None, optional, default=None
        Logger

    Returns
    -------
    exists : bool
        Returns True if the file path exists, False otherwise.
    """
    if verbose and not Path(filename).exists():
        if a_logger is not None:
            a_logger.error(f"Cannot find {filename}. Not a valid path or file.")
        return False
    else:
        if a_logger is not None:
            a_logger.debug(f"Found {filename}.")
    return True


def CapitalCase_to_snake_case(input_str: str) -> str:  # noqa: N802
    """
    Converts CapitalCase or camelCase strings to snake_case

    Based on: https://www.geeksforgeeks.org/python-program-to-convert-camel-case-string-to-snake-case/

    Limitations:

    - Messes up consecutive capital letters, e.g. ISO become i_s_o, or
      FocalLengthMM becomes focal_length_m_m

    :param input_str: Input CapitalCase or camelCase string
    """

    output_str = "".join(["_" + i.lower() if i.isupper() else i for i in input_str]).lstrip("_")

    return output_str


def determine_file_type(fname: str) -> str:
    """
    Determine if file is a common graphics format or FITS.
    """

    # Graphics file formats, delegated to imageio
    graphics_list = [".tif", ".tiff", ".jpg", ".jp2", "jpeg", ".png", ".gif"]

    # FITS, delegate to astropy
    fits_list = [".fits", ".ftz", ".fit", ".fits.gz"]

    # want file extension, e.g .png, .fits but also .fits.gz
    root, ext = os.path.splitext(fname.lower())
    if ".gz" in ext:
        root, ext2 = os.path.splitext(root)
        ext = ext2 + ext
    ext = ext.lower()

    if ext in fits_list:
        ftype = "fits"
    elif ext in graphics_list:
        ftype = "graphics"
    else:
        ftype = "unknown"

    return ftype
