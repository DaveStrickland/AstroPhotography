#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_add_metadata.py
#
#  Adds metadata to the FITS header of a calibrated fits file.
#  Currently this is written to use the information in iTelescope
#  file names, plus information derived from that using astropy.
#  In the longer term other ways of inputting the data could be added.
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
#  2020-12-16 dks : Initial skeleton.

import argparse
import sys
import logging
import AstroPhotography as ap

def command_line_opts(argv):
    """Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_add_metadata',
        description=('Adds metadata to the FITS header of a calibrated fits file.'
        ' Currently this is written to use the information in iTelescope'
        ' file names, plus information derived from that using astropy.'
        ' In the longer term other ways of inputting the data could be added.'))
    
    # Required
    parser.add_argument('fitsimage',
        metavar='FITSIMAGE.FITS',
        help=('Path/name of the FITS image to added metadata to.'
        ' Note that the input file is modified.'))
        
    # Optional
    p_mode = 'iTelescope'

    parser.add_argument('--mode',
        metavar='MODE',
        default=p_mode,
        help=('Mode of operation. If "iTelescope" the telescope, target,'
        ', and observer will be parsed from the file name.'))
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    return args

def main(args=None):
    p_args       = command_line_opts(args)
    p_fitsimg    = p_args.fitsimage
    p_mode       = p_args.mode
    p_loglevel   = p_args.loglevel
    
    # Create an instance of the calibrator.
    meta_adder = ap.ApAddMetadata(p_loglevel)
    
    # Update the image.
    meta_adder.process(p_fitsimg, p_mode)
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
