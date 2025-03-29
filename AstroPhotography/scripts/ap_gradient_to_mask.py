#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_gradient_to_mask.py
#
#  Perform gaussian gradient magnitude filtering on an image and
#  create a boolean mask where the gradient exceeds a configurable
#  threshold above the median image gradient. If a list of images is
#  provided the boolean masks of each are added to produce a single
#  mask.
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
#  2025-03-15 dks : Initial skeleton.

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
        prog="ap_gradient_to_mask",
        description=(
            "Perform gaussian gradient magnitude filtering on an image and"
            " create a boolean mask where the gradient exceeds a configurable"
            " threshold above the median image gradient. If a list of images is"
            " provided the boolean masks of each are added to produce a single"
            " mask."
        ),
    )

    parser.add_argument(
        "input_fits",
        metavar="INPUT_FITS",
        nargs="+",
        help=(
            "Path/names of the input FITS images to be gradient filtered."
            " This can be the name of one file, a comma-separated list"
            " of files, or a space separated list of files."
        ),
    )
    parser.add_argument(
        "--output",
        metavar="OUTPUT",
        help=("Path/name of the output mask. This will be over-written" " if it already exists."),
    )

    parser.add_argument(
        "--extnum",
        default=0,
        help=(
            "Specify the HDU extension number or name (Default is 0)."
            " This should be the same for all the input files."
        ),
    )

    p_fwhm = 1.5
    parser.add_argument(
        "--fwhm_pix",
        default=p_fwhm,
        type=float,
        help=(
            "This is the FWHM of the Guassian gradient filter in units of pixels"
            " that is applied to a median-normalized copy of the input image(s)."
            f" Default value: {p_fwhm:.2f} pixels."
        ),
    )

    p_sigma = 5.0
    parser.add_argument(
        "--sigma",
        default=p_sigma,
        type=float,
        help=(
            "The gradient map is thresholded such that pixels where the gradient"
            " is greater than sigma times the gradient standard deviation above"
            " the mean gradient value are defined as holes or the edges of holes."
            " (holes are in-filled and then have a binary dilation applied to them.)"
            f" Default value: {p_sigma:.2f}."
        ),
    )

    parser.add_argument(
        "-l", "--loglevel", default="INFO", help="Logging message level. Default: INFO"
    )

    # --------------------------------------------------------------------------
    # any other commmands
    args = parser.parse_args(argv)
    args.input_files = parse_arg_list(args.input_fits)  # handle comma-separated list
    return args


def parse_arg_list(possible_arg_list: list[str] | str) -> list[str]:
    """
    Rationalize arguments that may be a single value, a list, or
    optionally a comma-separted string
    """

    if isinstance(possible_arg_list, str):
        clean_list = [possible_arg_list]
    elif isinstance(possible_arg_list, list):
        clean_list = []
        numel = len(possible_arg_list)
        if numel == 1:
            if "," in possible_arg_list[0]:
                split_list = possible_arg_list[0].split(",")
                for el in split_list:
                    clean_list.append(str(el))
            else:
                # single element list with no comma
                clean_list = [possible_arg_list[0]]
        else:
            clean_list = []
            for el in possible_arg_list:
                clean_list.append(str(el))

    return clean_list


def main(args=None):
    p_args = command_line_opts(args)
    holefinder = ap.ApGradientToMask(p_args.loglevel)
    holefinder.mask_files(
        input_files=p_args.input_fits,
        outhole_mask=p_args.output,
        extnum=p_args.extnum,
        fwhm_pix=p_args.fwhm_pix,
        sigma=p_args.sigma,
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
