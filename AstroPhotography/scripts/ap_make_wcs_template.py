#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_make_wcs_template.py
#
#  Create an empty FITS file consisting of a WCS header defined by
#  user-specified values. The primary (only?) use of such a file is
#  as a template when resampling a real image with some given WCS to
#  another WCS.
#
#  Copyright 2025 Dave Strickland <dave.strickland@gmail.com>
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
#  2025-06-28 dks : Initial skeleton.

import argparse

##import sys
import logging
import AstroPhotography as ap
from typing import Any
from astropy.io import fits


def command_line_opts(argv):
    """Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(
        prog="ap_make_wcs_template",
        description=(
            "Generate a raster grey-scale or RGB bitmap image from an FITS image"
            ", applying user specified normalization and limits."
        ),
    )

    # defaults
    p_naxis1 = 2048  # Number of cols
    p_naxis2 = 2048  # Number of rows
    p_cdelt1 = -1.0  # East to the left, 1 arcsec (columns) pixels
    p_cdelt2 = 1.0  # North up, 1 arcsec (row) pixels
    p_crota = None  # Don't specify rotation
    p_no_cd_mtx = False  # By default, create a CD matrix instead of CDELT1, CDELT2

    parser.add_argument(
        "--output",
        type=str,
        required=True,
        metavar="OUTPUT_FITS_FILE",
        help="Name for the output FITS file to be created. Will be overwritten if it exists",
    )
    parser.add_argument(
        "--ra",
        type=float,
        required=True,
        metavar="RA_DEG",
        help="Right Ascension (J2000) of the center of the image in decimal degrees.",
    )
    parser.add_argument(
        "--dec",
        type=float,
        required=True,
        metavar="DEC_DEG",
        help="Declination (J2000) of the center of the image in decimal degrees.",
    )

    parser.add_argument(
        "--naxis1",
        default=p_naxis1,
        type=int,
        help=(f"The number of pixel columns in the image. (Default: {p_naxis1})."),
    )
    parser.add_argument(
        "--naxis2",
        default=p_naxis2,
        type=int,
        help=(f"The number of pixel rows in the image. (Default: {p_naxis1})."),
    )

    parser.add_argument(
        "--cdelt1",
        default=p_cdelt1,
        type=float,
        help=(f"Width of a single pixel column in arcseconds. (Default: {p_cdelt1})."),
    )
    parser.add_argument(
        "--cdelt2",
        default=p_cdelt2,
        type=float,
        help=(f"Height of a single pixel row in arcseconds. (Default: {p_cdelt2})."),
    )

    parser.add_argument(
        "--nocd_matrix",
        action="store_true",
        default=p_no_cd_mtx,
        help=(
            "If specified then a CD matrix will not be used, and the older"
            " method of using CDELT1, CDELT1, and CROTA2 will be used."
            " Using a CD matrix is the preferred standard."
            f" (Default: {p_no_cd_mtx})"
        ),
    )

    parser.add_argument(
        "--crota",
        default=p_crota,
        help=(
            "Rotation angle of the Y-axis away from North in decimal degrees."
            f" (Default: {p_crota})."
        ),
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help=(
            "If specified then a summary of the generated WCS will be written to"
            f" stdout. (Default: {False})"
        ),
    )

    args = parser.parse_args(argv)
    return args


def main(args=None):
    p_args = command_line_opts(args)
    if p_args.crota is not None:
        p_args.crota = float(p_args.crota)

    arcsec_to_deg = 1.0 / 3600.0
    naxis = (p_args.naxis1, p_args.naxis2)
    crval = (p_args.ra, p_args.dec)
    cdelt = (p_args.cdelt1 * arcsec_to_deg, p_args.cdelt2 * arcsec_to_deg)
    do_cd_mtx = not p_args.nocd_matrix

    ap.util.make_dothead_from_keywords(
        p_args.output,
        naxis,
        crval,
        cdelt,
        format="fits",
        verbose=p_args.verbose,
        do_cd_mtx=do_cd_mtx,
        crota=p_args.crota,
    )

    # Summary of WCS printed to STDOUT
    if p_args.verbose:
        w, _, _, _ = ap.util.load_wcs_from_file(p_args.output, extnum=0, verbose=False)
        _ = ap.util.summarize_wcs(w, verbose=True)
    return 0


if __name__ == "__main__":
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
