#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_find_badpix.py
#
#  Generates a bad pixel mask given a master dark file.
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
#  2020-10-11 dks : Initial skeleton. 
#  2020-10-17 dks : Working version completed.
#  2020-12-31 dks : Move ApFindBadPixels out of script into core.

import argparse
import sys
import logging

import AstroPhotography as ap

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_find_badpix',
        description='Generates a bad pixel mask given a master dark or master bias FITS file.')
        
    # Required
    parser.add_argument('masterdark',
        metavar='IN_MASTER_DARK.FITS',
        help='Path/name of the input master dark/bias to use.')
    parser.add_argument('badpixfile',
        metavar='OUT_BADPIX.FITS',
        help='Path/name of the output badpix file to generate.')
        
        
    # Optional
    p_sigma = 4.0
    parser.add_argument('--sigma',
        metavar='NSIGMA',
        default=p_sigma,
        help=('Number of MAD standard deviations to use in sigma clipping.'
            ' After clipping pixels that are more than this number of sigma'
            f' from the median will be marked bad. Default: {p_sigma:.2f}'))
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    return args

                
def main(args=None):
    p_args      = command_line_opts(args)
    p_inmstrdrk = p_args.masterdark
    p_outbadpix = p_args.badpixfile
    p_sigma     = p_args.sigma
    p_loglevel  = p_args.loglevel
    
    mkbadpix = ap.ApFindBadPixels(p_inmstrdrk,
        p_sigma,
        p_loglevel)
    mkbadpix.write_mask(p_outbadpix)
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
