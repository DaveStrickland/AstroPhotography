# -*- coding: utf-8 -*-
#
#  Contains the implementation of various Astrophotography name utilities
#
#  Copyright 2024 Dave Strickland <dave.strickland@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
# SPDX-License-Identifier: GPL-3.0-or-later

# 2024-04-07 dks : issue-002 Initial coding
# 2025-03-01 dks : Split off from ApUtil.py

##import logging
##import os.path
from typing import Any

##from matplotlib.patches import Ellipse
##import time
##from datetime import datetime

# AstroPhotography includes
##from .. import __version__


def namefn_calibrated_input(input_file, input_rootname, input_suffix, ap_filetype):
    """
    Given the name of a calibrated input file, generate the file name for
    and output directory for one of the subsequent output file types.

    For calibrated input files the the allowed output file types are:

    - ``input``: This is the calibrated input file itself.
    - ``srclist``: Star detection output FITS table file.
    - ``regfile``: ds9-format region file.
    - ``plotfile``: PNG plots of the input image with the detected star-like sources
      plotted as circles.
    - ``qualfile``: YaML source detection quality summary file.
    - ``fwhmplot``: PNG plot zooming in around a subset of the detected
      sources in the input image.
    - ``navfile``: A navigated and renamed copy of the input file, but now
      with a valid WCS solution.

    Parameters
    ----------
    input_file : str
        The name of the calibrated input file.
    input_rootname : str, optional, default=None
               The part of the input file name that are shared and that
               designate them being the calibrated files. If not specified
               it is assumed that the file is from iTelescope or was
               created by the AstroPhotography module itself, and will
               have an input_rootname of either 'Calibrated-iTelescope',
               'calibrated', or just 'cal' .
               The output file name from this function replace the rootname with
               a file-type specific prefix, as described in the class
               documentation and get_file_names function documentation.
    input_suffix : str, optional, default='.fits'
               String denoting the file type suffix of the input image file.
               For example, '.fits' or '.fits' or '.ftz' or 'fits.gz' or '.fits.bz2'
    ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile'}
        The output file type for which the name should be returned.
        This should be one of the stage names described above.

    Returns
    -------
    output_file : str
        File name that would correspond to the named output file type given
        the name of the input file. This file does not necessarily exist yet.
        The output file name does **not** include any subdirectory name.
        It is up to the caller to combine output_dir and output_file if
        necessary.
    output_dir : str
        Directory name in which output_file will be written, relative
        to the root directory established by the input files.

    Raises
    ------
    RuntimeError :
        If input_rootname is None but the input_file name does not match
        any of the expected iTelescope/AstroPhotography root names.

    Warnings
    --------
    This function assumes that the input file is not in a subdirectory,
    i.e. there is no path within in ``input_file`` string. In the longer
    term this should be rewritten using the pathlib module without such
    an assumption.
    """

    allowed_stages = ["input", "srclist", "regfile", "plotfile", "qualfile", "fwhmplot", "navfile"]
    if ap_filetype not in allowed_stages:
        err_msg = (
            f"Requested ap_filetype {ap_filetype} not one of the allowed values: {allowed_stages}"
        )
        raise RuntimeError(err_msg)

    # Determine the root name to use.
    default_roots = ["Calibrated-iTelescope", "calibrated", "cal"]
    file_root = input_rootname
    if file_root is None:
        for a_root in default_roots:
            if a_root in input_file:
                file_root = a_root
                break
        if file_root is None:
            err_msg = (
                f"Error in namefn_calibrated_input: input file {input_file}"
                f" matches none of the expected file roots: {default_roots}"
            )
            raise RuntimeError(err_msg)

    name_conv_dict = _get_name_conv_dict(file_root)

    conv = name_conv_dict[ap_filetype]
    output_file = input_file
    output_dir = conv["dir"]
    if conv["replace"] is not None:
        output_file = output_file.replace(conv["replace"], conv["with"])
    if conv["extension"] is not None:
        output_file = output_file.replace(input_suffix, conv["extension"])

    return output_file, output_dir


def namefn_getdir(ap_filetype: str) -> str:
    """
    Given a file type for star detection/astrometry, return the
    output directory such files would be placed in.

    For calibrated input files the the allowed output file types are:

    - ``input``: This is the calibrated input file itself.
    - ``srclist``: Star detection output FITS table file.
    - ``regfile``: ds9-format region file.
    - ``plotfile``: PNG plots of the input image with the detected star-like sources
      plotted as circles.
    - ``qualfile``: YaML source detection quality summary file.
    - ``fwhmplot``: PNG plot zooming in around a subset of the detected
      sources in the input image.
    - ``navfile``: A navigated and renamed copy of the input file, but now
      with a valid WCS solution.

    Parameters
    ----------
    ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile'}
        The AstroPhotography file type for which the directory should be returned.
        This should be one of the stage names described above.

    Returns
    -------
    output_dir : str
        Directory name in which output files will be written, relative
        to the root directory established by the input files.
    """

    allowed_stages = ["input", "srclist", "regfile", "plotfile", "qualfile", "fwhmplot", "navfile"]
    if ap_filetype not in allowed_stages:
        err_msg = (
            f"Requested ap_filetype {ap_filetype} not one of the allowed values: {allowed_stages}"
        )
        raise RuntimeError(err_msg)

    # True file root not needed
    file_root = "Calibrated-iTelescope"
    name_conv_dict = _get_name_conv_dict(file_root)
    conv = name_conv_dict[ap_filetype]
    output_dir = conv["dir"]
    return output_dir


def _get_name_conv_dict(file_root):
    """
    Utility function used by :func:`namefn_calibrated_input` that returns
    the file name and output directory dictionary given a file name root.

    Parameters
    ----------
    file_root : str
        The ``file_root`` is the unique string all input calibrated FITS files
        either start with, or include. For calibrated files produced by
        Astrophotography this might be the string `cal`, while iTelescope
        uses ``Calibrated`` or ``calibrated``.

    Returns
    -------
    name_conv_dict : dict
        Name conversion dictionary used to convert input calibrated file
        input names and paths into output file names.

    See Also
    --------
    :class:`ApProcess` : Batch processing of calibrated files
    """

    # Settings for output file names and output file paths
    # - This is difficult to fully automate without some assumptions about
    #   the input file names
    #   - I assume that all input file have the same file name prefix, e.g.
    #     Calibrated, or cal
    #   - I assume that all input files have the same file extension, e.g.
    #     fits or (sadness) fit
    # - I prefer to place output files of a certain type in a subdirectories.
    #   If you don't want all the diagnostic plots then it may just be simpler
    #   to not deal with subdirectories.
    # - If dir is not None then the various outputs files will be written to
    #   directories with the specified path relative to the **current** directory.
    #   The directory will be created if it not already present.
    name_conv_dict = {
        "input": {"replace": None, "with": None, "extension": None, "dir": "./"},
        "srclist": {
            "replace": file_root,
            "with": "srclist",
            "extension": ".fits",
            "dir": "./SourceLists/",
        },
        "regfile": {
            "replace": file_root,
            "with": "ds9",
            "extension": ".reg",
            "dir": "./SourceLists/",
        },
        "plotfile": {
            "replace": file_root,
            "with": "implot",
            "extension": ".png",
            "dir": "./SourceLists/",
        },
        "fwhmplot": {
            "replace": file_root,
            "with": "fwhmplot",
            "extension": ".png",
            "dir": "./SourceLists/",
        },
        "qualfile": {
            "replace": file_root,
            "with": "qual",
            "extension": ".yaml",
            "dir": "./MetaData/",
        },
        "navfile": {
            "replace": file_root,
            "with": "navigated",
            "extension": ".fits",
            "dir": "./NavigatedImages/",
        },
    }
    return name_conv_dict
