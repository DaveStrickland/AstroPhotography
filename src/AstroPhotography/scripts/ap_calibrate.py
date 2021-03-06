#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_calibrate.py
#
#  Performs the following calibration steps on an input image:
#  - Bias subtraction
#  - Dark subtraction
#  - Optional flat fielding correction.
#  - Optional bad pixel removal
#  - Optional cosmic ray removal.
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
#  2020-12-05 dks : Initial skeleton, based on example_ccdproc.py and 
#                   ap_fix_badpix.py
#  2021-01-14 dks : Add cosmic ray removal to match ap_fix_cosmic_rays.py

import argparse
import sys
import logging
import AstroPhotography as ap

def command_line_opts(argv):
    """Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_calibrate',
        description=('Performs calibration of raw astronomical images'
        ' by applying bias and dark frame subtraction, along with'
        ' optional (but recommended) flat fielding, bad pixel correction,'
        ' and cosmic ray removal.'))
    
    # Required
    parser.add_argument('raw_image',
        metavar='INPUT_IMAGE.FITS',
        help='Path/name of the raw (uncalibrated) input image.')
    parser.add_argument('master_bias',
        metavar='MBIAS.FITS',
        help='Path/name of the master bias file.')
    parser.add_argument('master_dark',
        metavar='MDARK.FITS',
        help=('Path/name of the master dark file.'
        ' WARNING. It is assumed that the master dark has not had the'
        ' bias subtracted. File an issue if you have example where this'
        ' is not the case.'))
    parser.add_argument('calibrated_image',
        metavar='CALIBRATED_IMAGE.FITS',
        help='Path/name of the output calibrated image.')        
        
    # Optional
    p_delta = 2

    parser.add_argument('--master_flat',
        metavar='MFLAT.FITS',
        default=None,
        help=('Path/name of the master flat file.'
        ' If this is specified flat fielding will be applied to the final image.'))
    parser.add_argument('--master_badpix',
        metavar='BADPIX.FITS',
        default=None,
        help=('Path/name of the master badpixel file.'
        ' If supplied bad pixel correction will be applied.'
        ' This is an integer image with zero at good pixels, non-zero at'
        ' bad pixels, e.g. as generated by the ap_find_badpix.py script.'))
    parser.add_argument('--normflat',
        metavar='NORMALIZED_FLAT.FITS',
        default=None,
        help='Optional output of normalized flat used. Useful for diagnostic purposes.')
    parser.add_argument('--deltapix',
        default=p_delta,
        type=int,
        help=('Half-width and half-height of the box surrounding a bad'
        ' pixel from which the median value of any good pixels will be'
        ' used to replace the value in the bad pixel. If deltapix=1 then'
        ' the 8 surrounding pixels will be used, if deltapix=2 then the'
        ' surrounding 24 pixels will be used. Values above 2 are not recommended.'
        f' Default: {p_delta} pixels.'))
    parser.add_argument('--fixcosmic',
        default=False,
        action='store_true',
        help=('If specified, cosmic ray removal will be performed.'
        ' This removes many cosmic rays and remaining flickering bad'
        ' pixels, but does not work well on alpha particle CRs.'))
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')

    # To be added?
    # yaml file of additional observatory related keywords to be added
    # to fits headers
                
    args = parser.parse_args(argv)
    return args

def main(args=None):
    p_args       = command_line_opts(args)
    p_raw_img    = p_args.raw_image
    p_mbias      = p_args.master_bias
    p_mdark      = p_args.master_dark
    p_out_img    = p_args.calibrated_image
        
    p_mflat      = p_args.master_flat
    p_mbadpix    = p_args.master_badpix

    p_normflat   = p_args.normflat
    p_deltapix   = p_args.deltapix
    p_fixcosmic  = p_args.fixcosmic
    p_loglevel   = p_args.loglevel
    
    # Create an instance of the calibrator.
    calibrator = ap.ApCalibrate(p_mbias,
        p_mdark,
        p_mflat,
        p_mbadpix,
        p_loglevel)
    
    # Tell it to calibrate data from files (rather than numpy arrays).
    calibrator.calibrate(p_raw_img, 
        p_out_img,
        p_deltapix,
        p_normflat,
        p_fixcosmic) 
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
