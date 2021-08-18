#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_auto_badcol.py
#
#  Attempts to automatically detect statistically bad columns and rows
#  in a given FITS image, reporting the output in a format that can be
#  cut and pasted into a user badpixel YaML file.
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
#  2021-08-14 dks : Initial skeleton.

import argparse
import sys
import logging
import AstroPhotography as ap
import numpy as np

def command_line_opts(argv):
    """Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_auto_badcol',
        description=('Attempts to automatically detect statistically'
        ' bad columns and rows'
        ' in a given FITS image, reporting the output in a format that can be'
        ' cut and pasted into a user badpixel YaML file.'))
    
    # Required
    parser.add_argument('fitsimage',
        metavar='FITSIMAGE.FITS',
        help=('Path/name of input FITS image to look for bad columns or rows in.'
        ' The image data is assumed to be in the primary extension of the FITS file.'))
        
    # Optional
    p_sigma      = 5.0
    p_window_len = 11

    parser.add_argument('--sigma',
        metavar='NSIGMA',
        default=p_sigma,
        help=('Columns or rows are identified as being bad if the median'
        ' pixel value is more than NSIGMA standard deviations away from'
        ' the locally determined average.'
        f' Default: {p_sigma}'))
    parser.add_argument('--window',
        metavar='LENGTH',
        default=p_window_len,
        help=('Size of the moving average window function used to generate'
        ' the local estimate of the column or row value.'
        f' Default value: {p_window_len}'))
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    return args

def main(args=None):
    p_args       = command_line_opts(args)
    p_fitsimg    = p_args.fitsimage
    p_sigma      = p_args.sigma
    p_window     = p_args.window
    p_loglevel   = p_args.loglevel
    retcode      = 0
        
    # Create an instance of the ApAutoBadcols.
    auto_badcols = ap.ApAutoBadcols(p_loglevel)
    
    # Process based on a file
    badcols, badrows = auto_badcols.process_fits(p_fitsimg, 
        p_sigma, 
        p_window)
    # TODO output to STDOUT in yaml-like format.
    print(f'# Auto bad columns from {p_fitsimg}, sigma={p_sigma}, window_len={p_window}')
    if badcols is not None:
        if len(badcols) > 0:
            # Add one to get FITS-like indexing
            badcols = badcols + 1

            print('bad_columns:')
            for val in badcols:
                print(f'- {val:d}')
        else:
            print('bad_columns: {}') # show it is empty
    else:
        print('# No bad columns detected.')
    if badrows is not None:
        if len(badrows) > 0:
            # Add one to get FITS-like indexing
            badrows = badrows + 1

            print('bad_rows:')
            for val in badrows:
                print(f'- {val:d}')
        else:
            print('bad_rows: {}')   # show it is empty
    else:
        print('# No bad rows detected.')
            
    return retcode

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
