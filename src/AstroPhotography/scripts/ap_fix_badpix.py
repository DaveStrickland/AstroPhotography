#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_fix_badpix.py
#
#  Patches an image given a bad pixel mask, replacing bad pixels with
#  the median value of the surrounding good pixels.
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
#  2020-10-24 dks : Initial skeleton. 
#  2020-11-16 dks : Finalized along with ApFixBadPixels.

import argparse
import sys
import logging
import AstroPhotography as ap

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_fix_badpix',
        description=('Patches an image given a bad pixel mask,'
        ' replacing bad pixels with the median value of the'
        ' surrounding good pixels.'))
        
    # Required
    parser.add_argument('raw_image',
        metavar='INPUT_IMAGE.FITS',
        help='Path/name of the input image that has bad pixels to patch.')
    parser.add_argument('master_badpix',
        metavar='BADPIX.FITS',
        help=('Path/name of the master badpixel file.'
        ' This is an integer image with zero at good pixels, non-zero at'
        ' bad pixels, e.g. as generated by the ap_find_badpix.py script.'))
    parser.add_argument('fixed_image',
        metavar='OUTPUT_IMAGE.FITS',
        help='Path/name of the output patched image.')        
        
    # Optional
    p_delta = 2
    parser.add_argument('--deltapix',
        default=p_delta,
        type=int,
        help=('Half-width and half-height of the box surrounding a bad'
        ' pixel from which the median value of any good pixels will be'
        ' used to replace the value in the bad pixel. If deltapix=1 then'
        ' the 8 surrounding pixels will be used, if deltapix=2 then the'
        ' surrounding 24 pixels will be used. Values above 2 are not recommended.'
        f' Default: {p_delta} pixels.'))
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    return args

def main(args=None):
    p_args       = command_line_opts(args)
    p_inp_img    = p_args.raw_image
    p_inp_badpix = p_args.master_badpix
    p_out_img    = p_args.fixed_image
    
    p_deltapix   = p_args.deltapix
    p_loglevel   = p_args.loglevel
    
    # Create an instance of the object that can fix bad pixels.
    fixpix = ap.ApFixBadPixels(p_loglevel)
    
    # Tell it to fix data from files (rather than numpy arrays).
    fixpix.fix_files(p_inp_img, 
        p_inp_badpix,
        p_out_img,
        p_deltapix) 
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
