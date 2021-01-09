#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_fix_cosmic_rays.py
#
#  Attempt to clean cosmic rays from a well-sampled CCD image. 
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
#  2021-01-09 dks : Initial skeleton. 

import argparse
import sys
import logging

import AstroPhotography as ap

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_fix_cosmic_rays',
        description=('Clean cosmic rays from a CCD image using the'
        ' ccdproc implementation of the L.A. Cosmic algorithm.'
        ' The input image should have already undergone bad pixel'
        ' correction (e.g. using ap_fix_badpix.py) before attempting'
        ' cosmic ray correction.'))
        
    # Required
    parser.add_argument('input',
        metavar='INPUT.FITS',
        help=('Path/name of the input FITS image. Note that it is'
        ' recommended that bad pixel correction should have already'
        ' been performed on this data.'))
    parser.add_argument('output',
        metavar='OUTPUT.FITS',
        help='Path/name of the output cosmic-ray cleaned FITS image.')
        
    # Optional
    parser.add_argument('--crdiffim',
        metavar='CR_DIFFERENCE_IMG.FITS',
        default=None,
        help=('Name of the optional difference image, specifically the '
        'output CR cleaned image subtracted from the original image.'))
    parser.add_argument('--crmaskim',
        metavar='CR_MASK_IMG.FITS',
        default=None,
        help=('Name of the optional CR mask image, which is non-zero'
        'at the location of the identified cosmic rays.'))
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    return args

                
def main(args=None):
    p_args     = command_line_opts(args)
    p_inpfile  = p_args.input
    p_outfile  = p_args.output
    p_crdiff   = p_args.crdiffim
    p_crmask   = p_args.crmaskim
    p_loglevel = p_args.loglevel
    
    crfixer = ap.ApFixCosmicRays(p_loglevel)
    crfixer.process_file(p_inpfile, p_outfile)
    
    if p_crdiff is not None:
        crfixer.write_crdiff_img(p_crdiff)
    if p_crmask is not None:
        crfixer.write_crmask_img(p_crmask)
        
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
