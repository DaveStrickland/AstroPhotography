#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ap_combine_darks.py
#
#  Generates a master dark or master bias file, given a directory 
#  that contains the raw dark or bias file to be combined. Use separate
#  directories for calcilbration files with different exposure time,
#  binning, CCD temperature, and from different telescopes.
#  
#  Copyright 2020 Dave Strickland <dave.strickland@gmail.com>
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
#  2020-10-01 dks : Initial coding begun. 

import argparse
import sys
import logging
from pathlib import Path
import math

import numpy as np
import matplotlib                # for rc
import matplotlib.pyplot as plt
from astropy.io import fits

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_combine_darks',
        description='Generates a master dark or master bias file from' + 
        ' all calibration FITS files in a given directory.')
        
    # Required
    parser.add_argument('rawcaldir',
        metavar='RAW_CAL_DIR',
        help='The directory in which the raw calibration' +
        ' files to be combined are to be found.')
    parser.add_argument('master_filename',
        metavar='MASTER_CAL_FILENAME',
        help='Output file name for master calibration file.')
        
    # Optional
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    return args
                
def main(args=None):
    p_args       = command_line_opts(args)
    p_root       = p_args.rawcaldir
    p_masterfile = p_args.master_filename
    p_loglevel   = p_args.loglevel
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
