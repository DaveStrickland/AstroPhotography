#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_fits_rasterizer.py
#
#  Perform arithmetic operations on a FITS file, using either a constant
#  value applied to the whole primary array or the pixels values in
#  another FITS file of the same shape. This is similar to the Heatools
#  farith, or the cfitsio example imarith.
#
#  Copyright 2021 Dave Strickland <dave.strickland@gmail.com>
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
#  2021-05-02 dks : Initial skeleton.

import argparse
import sys
import logging
import AstroPhotography as ap


def command_line_opts(argv):
    """Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(
        prog="ap_fits_rasterizer",
        description=(
            "Generate a raster grey-scale or RGB bitmap image from an FITS image"
            ", applying user specified normalization and limits."
        ),
    )

    parser.add_argument(
        "-l", "--loglevel", default="INFO", help="Logging message level. Default: INFO"
    )

    # --------------------------------------------------------------------------
    # Commands
    subparsers = parser.add_subparsers(
        title="Commands", dest="command", required=True, help="Available subcommands. Choose one."
    )

    grey_parser = subparsers.add_parser(
        "grey",
        help=(
            "Produce a greyscale raster image from the data in the"
            " input FITS image using the specified backend to do the processing."
        ),
    )
    rgb_parser = subparsers.add_parser(
        "rgb",
        help=(
            "Produce a RGB raster image from the data in the three"
            " input FITS images using the specified backend to do the processing."
        ),
    )

    # --------------------------------------------------------------------------
    # grey command
    # Required
    grey_parser.add_argument(
        "input_fits",
        metavar="INPUT_FITS",
        help="Path/name of the input FITS image to perform arithmetic on.",
    )
    grey_parser.add_argument(
        "output_image",
        metavar="OUTPUT",
        help=(
            "Path/name of the output raster image. This will be over-written"
            " if it already exists. Note that some backends only produce certain"
            " files types, irrespective of the file extension you enter here."
        ),
    )

    # --------------------------------------------------------------------------
    # rgb command
    # Required
    rgb_parser.add_argument(
        "input_fits",
        metavar="INPUT_FITS",
        help=(
            "Path/names of the input FITS images."
            " The three input files must be specified with the FITS image to become"
            " the red channel first, then the image to become green, and finally"
            " the image to become the blue channel last."
        ),
    )
    rgb_parser.add_argument(
        "output_image",
        metavar="OUTPUT",
        help=(
            "Path/name of the output raster image. This will be over-written"
            " if it already exists. Note that some backends only produce certain"
            " files types, irrespective of the file extension you enter here."
        ),
    )

    args = parser.parse_args(argv)
    return args


def main(args=None):
    p_args = command_line_opts(args)
    if "grey" in p_args.command:
        print("grey")
    elif "rgb" in p_args.command:
        print("rgb")

    return 0


if __name__ == "__main__":
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
