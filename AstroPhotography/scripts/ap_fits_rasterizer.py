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

##import sys
import logging
import AstroPhotography as ap
from typing import Any


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

    # --------------------------------------------------------------------------
    # Commands
    def_backend = "stiff"
    subparsers = parser.add_subparsers(
        title="Commands",
        dest="command",
        required=True,
        help="Available subcommands. Choose one.",
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

    grey_parser.add_argument(
        "--backend",
        default=def_backend,
        type=str,
        help=(
            "Specify the backend framework used to generate the bitmap output. "
            f"(Default is {def_backend}). Note that the stiff backend can only produce"
            " TIFF files."
        ),
    )
    grey_parser.add_argument(
        "--extnum",
        default=0,
        help="Specify the HDU extension number or name (Default is 0).",
    )
    grey_min_group = grey_parser.add_mutually_exclusive_group()
    grey_min_group.add_argument(
        "--min_cut",
        type=float,
        default=None,
        help="The pixel value of the minimum cut level (Default is the image minimum).",
    )
    grey_max_group = grey_parser.add_mutually_exclusive_group()
    grey_max_group.add_argument(
        "--max_cut",
        type=float,
        default=None,
        help="The pixel value of the maximum cut level (Default is the image maximum).",
    )
    grey_min_group.add_argument(
        "--min_percent",
        type=float,
        default=None,
        help=(
            "The percentile value used to determine the minimum cut"
            " level (Default is None). For most astronomical images"
            " the sky background is around the 50th percentile."
        ),
    )
    grey_max_group.add_argument(
        "--max_percent",
        type=float,
        default=None,
        help=(
            "The percentile value used to determine the maximum cut"
            "  level (Default is None). For most astronomical images"
            "  a percentile level of 99.9%% or slightly lower is a good choice."
        ),
    )

    grey_parser.add_argument(
        "--negative",
        action="store_true",
        default=False,
        help=(
            "Invert the color table so that bright pixels are black and"
            " faint pixels are white, similar to photographic negatives."
        ),
    )

    grey_parser.add_argument(
        "--binning",
        type=int,
        default=1,
        help=(
            "Bin input pixels by this factor on each dimension before"
            " generating the output image. For example, with binning=2"
            " the output image will have half the number of rows and half"
            " the number of columns than the input image, and only a quarter as many pixels."
        ),
    )

    grey_parser.add_argument(
        "--ignore_null",
        action="store_true",
        default=False,
        help=(
            "Ignore null pixels when computing the pixel value statistics used"
            " by the --min_percent and --max_percent arguments. The null value"
            " to be ignored is specified using the --null_value argument."
        ),
    )

    p_nullval: float = 0
    grey_parser.add_argument(
        "--null_value",
        type=float,
        default=p_nullval,
        help=(
            f"Null pixel value that is ignored if --ignore_null is specified. Default: {p_nullval}"
        ),
    )

    grey_parser.add_argument(
        "-l", "--loglevel", default="INFO", help="Logging message level. Default: INFO"
    )

    # --------------------------------------------------------------------------
    # rgb command
    # Required
    rgb_parser.add_argument(
        "input_fits",
        metavar="INPUT_FITS",
        nargs=3,
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

    rgb_parser.add_argument(
        "--backend",
        default=def_backend,
        type=str,
        help=(
            "Specify the backend framework used to generate the bitmap output. "
            f"(Default is {def_backend}). Note that the stiff backend can only produce"
            " TIFF files."
        ),
    )
    rgb_parser.add_argument(
        "--extnum",
        default=0,
        help=(
            "Specify the HDU extension number or name (Default is 0)."
            " This should be the same for all the input files."
        ),
    )
    rgb_min_group = rgb_parser.add_mutually_exclusive_group()
    rgb_min_group.add_argument(
        "--min_cut",
        nargs="*",
        const=None,
        default=None,
        help=(
            "The pixel value (or three values) of the minimum cut level"
            " (Default is the image minimum)."
            " If one value is specified it is applied to all channels."
            " Alternatively three values may be specified in order of the"
            " red, green, and blue channels."
        ),
    )
    rgb_max_group = rgb_parser.add_mutually_exclusive_group()
    rgb_max_group.add_argument(
        "--max_cut",
        nargs="*",
        const=None,
        default=None,
        help=(
            "The pixel value (or three values) of the maximum cut level"
            " (Default is the image maximum)."
            " If one value is specified it is applied to all channels."
            " Alternatively three values may be specified in order of the"
            " red, green, and blue channels."
        ),
    )
    rgb_min_group.add_argument(
        "--min_percent",
        nargs="*",
        const=None,
        default=None,
        help=(
            "The percentile value (or three values) used to determine the minimum cut"
            " level (Default is None). For most astronomical images"
            " the sky background is around the 50th percentile."
            " If one value is specified it is applied to all channels."
            " Alternatively three values may be specified in order of the"
            " red, green, and blue channels."
        ),
    )
    rgb_max_group.add_argument(
        "--max_percent",
        nargs="*",
        const=None,
        default=None,
        help=(
            "The percentile value (or three values) used to determine the maximum cut"
            "  level (Default is None). For most astronomical images"
            "  a percentile level of 99.9%% or slightly lower is a good choice."
            " If one value is specified it is applied to all channels."
            " Alternatively three values may be specified in order of the"
            " red, green, and blue channels."
        ),
    )

    rgb_parser.add_argument(
        "--negative",
        action="store_true",
        default=False,
        help=(
            "Invert the color table so that bright pixels are dark and"
            " faint pixels are bright, similar to photographic negatives."
        ),
    )

    rgb_parser.add_argument(
        "--binning",
        type=int,
        default=1,
        help=(
            "Bin input pixels by this factor on each dimension before"
            " generating the output image. For example, with binning=2"
            " the output image will have half the number of rows and half"
            " the number of columns than the input image, and only a quarter as many pixels."
        ),
    )

    rgb_parser.add_argument(
        "--ignore_null",
        action="store_true",
        default=False,
        help=(
            "Ignore null pixels when computing the pixel value statistics used"
            " by the --min_percent and --max_percent arguments. The null value"
            " to be ignored is specified using the --null_value argument."
        ),
    )

    p_nullval = 0
    rgb_parser.add_argument(
        "--null_value",
        type=float,
        default=p_nullval,
        help=(
            f"Null pixel value that is ignored if --ignore_null is specified. Default: {p_nullval}"
        ),
    )
    rgb_parser.add_argument(
        "-l", "--loglevel", default="INFO", help="Logging message level. Default: INFO"
    )

    # --------------------------------------------------------------------------
    # any other commmands
    args = parser.parse_args(argv)
    if "rgb" in args.command:
        args.min_cut = parse_arg_list(args.min_cut)
        args.max_cut = parse_arg_list(args.max_cut)
        args.min_percent = parse_arg_list(args.min_percent)
        args.max_percent = parse_arg_list(args.max_percent)
    return args


def parse_arg_list(possible_arg_list: list[Any] | None):
    """
    Rationalize arguments that may be a single value, a list, or
    optionally a comma-separted string
    """

    # Need to rationalize min_cut, max_cut, min_percent and max_percent
    # because fits_to_rbg expects 3-element lists
    if possible_arg_list is None:
        clean_list = None
    else:
        clean_list = []
        numel = len(possible_arg_list)
        if numel == 1:
            if "," in possible_arg_list[0]:
                split_list = possible_arg_list[0].split(",")
                split_numel = len(split_list)
                if split_numel == 3:
                    for el in split_list:
                        clean_list.append(float(el))
                else:
                    err_msg = (
                        "Expecting a list with either 1 or 3 elements, or None."
                        f"Got {numel} elements from {possible_arg_list}"
                    )
                    raise RuntimeError(err_msg)
            else:
                # single element list with no comma
                val = float(possible_arg_list[0])
                clean_list = [val, val, val]
        elif numel == 3:
            clean_list = []
            for el in possible_arg_list:
                clean_list.append(float(el))
        else:
            err_msg = (
                "Expecting a list with either 1 or 3 elements, or None."
                f"Got {numel} elements from {possible_arg_list}"
            )
            raise RuntimeError(err_msg)

    return clean_list


def main(args=None):
    p_args = command_line_opts(args)
    rasterizer = ap.ApFitsRasterizer(p_args.loglevel)
    if "grey" in p_args.command:
        rasterizer.fits_to_greyscale(
            p_args.input_fits,
            p_args.output_image,
            p_args.backend,
            p_args.extnum,
            p_args.min_cut,
            p_args.max_cut,
            p_args.min_percent,
            p_args.max_percent,
            p_args.negative,
            p_args.binning,
            p_args.ignore_null,
            p_args.null_value,
        )
    elif "rgb" in p_args.command:
        rasterizer.fits_to_rgb(
            p_args.input_fits,
            p_args.output_image,
            p_args.backend,
            p_args.extnum,
            p_args.min_cut,
            p_args.max_cut,
            p_args.min_percent,
            p_args.max_percent,
            p_args.negative,
            p_args.binning,
            p_args.ignore_null,
            p_args.null_value,
        )
    return 0


if __name__ == "__main__":
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
