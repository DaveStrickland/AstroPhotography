"""
Function based interface to FileWriter
"""

# SPDX-License-Identifier: GPL-3.0-or-later

from ..core.logger import logger

# import imageio
# import time
# import os.path
import yaml

# from pathlib import Path
# from datetime import datetime, timezone
from astropy.io import fits
from typing import Any
import numpy as np
import numpy.typing as npt

# AstroPhotography includes
from .. import __version__
from .check import _does_file_exist


def read_fits(
    image_filename: str,
    image_extension: str | int,
    a_logger: Any | None = None,
    remove_pedestal: str = "always",
) -> tuple[npt.NDArray, Any]:
    """
    Read a single extension's data and header from a FITS file

    Parameters
    ----------
    image_filename : str
        Name/path of the FITS file to be opended.
    image_extension : str or int
        FITS extension number or name from which to read the header
        and data.
    a_logger : logger instance or None, optional, default=None
        A logger instance, or None
    remove_pedestal : {'always', 'positive', 'never'}
        FITS ``PEDESTAL`` keyword handling. If the specified file/extension
        has the ``PEDESTAL`` keyword then ``always`` will remove the pedestal
        from the data using the value of the keyword, ``never`` will not
        modify the raw data array values irrespective of the value of the
        pedestal keyword, and ``positive`` will only remove the pedestal
        if doing so would leave a positive median pixel value. The latter
        option can be necessary in cases where other tools have modified
        the data but **not** removed the ``PEDSTAL`` keyword or set its
        value to zero.

    Returns
    -------
    ext_data : ndarray
        The data array associated with the specified file/extension
    ext_hdr : astropy.io.fits.Header
        The FITS header associated with the specified file/extension
    """

    if not _does_file_exist(image_filename, False):
        err_msg = ""
        raise RuntimeError(err_msg)

    if a_logger is not None:
        a_logger.info(f"Loading extension {image_extension} of FITS file {image_filename}")

    # open() parameters that can be important.
    # Default values used here.
    # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
    uint_handling = True
    image_scaling = False

    with fits.open(
        image_filename, uint=uint_handling, do_not_scale_image_data=image_scaling
    ) as hdu_list:
        ext_hdr = hdu_list[image_extension].header
        ext_data = hdu_list[image_extension].data

    ndim = ext_hdr["NAXIS"]
    cols = ext_hdr["NAXIS1"]
    rows = ext_hdr["NAXIS2"]
    bitpix = ext_hdr["BITPIX"]
    info_str = f"{ndim}-D BITPIX={bitpix} image with {cols} columns, {rows} rows"

    if ndim == 3:
        layers = ext_hdr["NAXIS3"]
        info_str += f", {layers} layers"

    if "BSCALE" in ext_hdr:
        bscale = ext_hdr["BSCALE"]
        info_str += f", BSCALE={bscale}"

    if "BZERO" in ext_hdr:
        bzero = ext_hdr["BZERO"]
        info_str += f", BZERO={bzero}"

    if a_logger is not None:
        a_logger.debug(info_str)

    if ndim == 3:
        err_msg = "Error, 3-D handling has not been implemented yet."
        if a_logger is not None:
            a_logger.error(err_msg)
        raise RuntimeError(err_msg)

    # Get data absolute limits.
    minval = np.amin(ext_data)
    maxval = np.amax(ext_data)
    medval = np.median(ext_data)
    if a_logger is not None:
        a_logger.debug(
            f"Raw data statistics are min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}"
        )

    # Is there a PEDESTAL value? MaximDL likes to add an offset, and
    # the PEDESTAL value is the value to ADD to the data to remove the
    # pedestal.
    if "PEDESTAL" in ext_hdr:
        if remove_pedestal in ["always", "positive"]:
            pedestal = float(ext_hdr["PEDESTAL"])
            condition_met = True
            if pedestal == 0:
                condition_met = False
            elif "positive" in remove_pedestal:
                # values would be on average negative if pedestal corrected for
                bg_negative = medval < pedestal
                if bg_negative:
                    condition_met = False
            if condition_met:
                if a_logger is not None:
                    a_logger.debug(f"Removing a PEDESTAL value of {pedestal} ADU.")
                ext_data += pedestal
                minval = np.amin(ext_data)
                maxval = np.amax(ext_data)
                medval = np.median(ext_data)
                if a_logger is not None:
                    a_logger.debug(
                        (
                            "After PEDESTAL removal, "
                            f"min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}"
                        )
                    )

    return ext_data, ext_hdr


def read_yaml_into_metadatadict(yamlfile: str) -> dict[str, tuple[Any, str]]:
    """
    Read a YaML file of FITS header metadata and converts that into a format
    that can easily be used with FITS files

    It is convenient to interact with FITS header records as a dictionary
    of keys (str) with assoiated tuples of value and comment,
    ``key_str: (value, comment_str)``. When output to YaML the tuples become
    lists, but astropy FITS I/O does not accept lists as the value of keyword
    keys. This function reads a YaML file (assumed to be in the correct format)
    and converts the lists back into tuples.

    For example, this metadata dictionary::

        meta_dict = {'OBJNAME':  ('M101 Supernova 2023ixf', 'Target'),
             'OBJECT':   ('M101 Supernova 2023ixf', 'Target'),
             'TELESCOP': ('iTelescope 24', 'Telescope name'),
             'INSTRUME': ('FLI-PL09000', 'Instrument name'),
             'XPIXSZ':   (12.0, 'Pixel Width in microns (after binning)'),
             'YPIXSZ':   (12.0, 'Pixel Height in microns (after binning)'),
             'EGAIN':    (1.53, 'Electronic gain in e-/ADU')}

    would appear as the following YaML::

        ---
        EGAIN:
        - 1.53
        - Electronic gain in e-/ADU
        INSTRUME:
        - FLI-PL09000
        - Instrument name
        OBJECT: &id001
        - M101 Supernova 2023ixf
        - Target
        OBJNAME: *id001
        TELESCOP:
        - iTelescope 24
        - Telescope name
        XPIXSZ:
        - 12.0
        - Pixel Width in microns (after binning)
        YPIXSZ:
        - 12.0
        - Pixel Height in microns (after binning)
        ...

    The normal ``yaml.safe_load`` converts that back into::

        # Reloaded data, back in python
        {'EGAIN': [1.53, 'Electronic gain in e-/ADU'],
            'INSTRUME': ['FLI-PL09000', 'Instrument name'],
            'OBJECT': ['M101 Supernova 2023ixf', 'Target'],
            'OBJNAME': ['M101 Supernova 2023ixf', 'Target'],
            'TELESCOP': ['iTelescope 24', 'Telescope name'],
            'XPIXSZ': [12.0, 'Pixel Width in microns (after binning)'],
            'YPIXSZ': [12.0, 'Pixel Height in microns (after binning)']
        }

    while this function aims to convert it back into the original format.

    Parameters
    ----------
    yamlfile : str
        Name of the yaml file to read

    Returns
    -------
    metadata_dict : dict[str, tuple[Any, str]]
        Metadata dictionary, of the form expected by (for example)
        :func:`ApProcess.process_all`
    """

    raw_dict = {}
    with open(yamlfile, "r") as f2:
        raw_dict = yaml.safe_load(f2)

    metadata_dict = {key: tuple(val) for key, val in raw_dict.items()}
    return metadata_dict
