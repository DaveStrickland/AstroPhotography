#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_astrometry.py
#  
#  Copyright 2020-2021 Dave Strickland <dave.strickland@gmail.com>
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
#  

# 2020-08-12 dks : Initial WIP
# 2020-08-29 dks : Update to keywords read from source list.
# 2020-12-20 dks : Fix read fits when BZERO undefined, increase timeout.
# 2020-12-28 dks : Disable SIP by default because swarp does not support it.
# 2021-01-19 dks : Moved ApAstrometry class over to core.
# 2022-02-08 dks : Make key optional, depending on whether already in config.

import argparse
import sys
import logging
from astroquery.astrometry_net import AstrometryNet

import AstroPhotography as ap
    
def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_astrometry',
        description='Generate an astrometric solution using astrometry.net')
        
    # Required
    parser.add_argument('inp_image',
        help='Input FITS image file without a vaid WCS solution.')
    parser.add_argument('source_list',
        help='Input FITS table of star positions found by ap_find_stars.')
    parser.add_argument('out_image',
        help='Ouput copy of the FITS image with WCS info (if found).')
        
    parser.add_argument('--key',
        default=None,
        type=str,
        metavar='ASTROMETRY_API_KEY',
        help=('Your personal Astrometry.net API key, if you have not already'
        ' added it to your ~/.astropy/config/astroquery.cfg config file.'))

    # Optional
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
    parser.add_argument('--image_extension',
        default=0,
        type=int,
        metavar='EXT_NUM',
        help='FITS extension number to read image data from. Default=0 (Primary)')
    parser.add_argument('--xy_extension',
        default='AP_XYPOS',
        metavar='EXT_NAME',
        help='FITS extension name for star X,Y position data. Default=AP_XYPOS')
    parser.add_argument('--use-sip',
        default=False,
        action='store_true',
        help=('Allow astrometry.net to fit SIP polynomial distortion terms.'
            ' This may be necessary for very large fields of view (>10 deg),'
            ' but SIP is not treated correctly by swarp (and possibly other'
            ' software).'))
    parser.add_argument('--user_scale',
        type=float,
        default=None,
        metavar="ARCSEC_PER_PIX",
        help=('Override the estimate plate scale in the source list file'
            ' and instead use a user spacified estimate of the plate scale.'
            ' The units are arcseconds/pixel.'))
    parser.add_argument('--scale_err_ratio',
        type=float,
        default=None,
        help=('The relative uncertainty in the estimated plate scale,'
            ' expressed as a ratio. This applies to either the default' 
            ' estimate from the source list, or a user-supplied plate'
            ' scale. For example, if the estimate plat scale is 2.0 arcsec/pix'
            ' and the scale_err_ratio=1.5, then the plate scale range that'
            ' will be search by Astrometry.net is 2/1.5 (=4/3) to 2*1.5 (=3)'
            ' arcseconds. If not specified ApAstrometry will use a value of 1.3.'
            ' Using a larger value can help in cases where astrometric'
            ' solutions fail, for example if incorrect telescope metadata'
            ' leads to inaccurate estimated plat scales.'))                
        
    args = parser.parse_args(argv)
    
    # Check that the Astrometry.ent key is present in the default
    # config if the user has not specified it here.
    if args.key is None:
        ast = AstrometryNet()
        if not ast.api_key:
            # Empty string evaluates False
            parser.error('Your astroquery config file does not contain an Astrometry.net key. You must correct that or supply one using --key.')
    
    return args    

def main(args=None):
    p_args        = command_line_opts(args)
    p_loglevel    = p_args.loglevel          # Loglevel
    p_inpimg      = p_args.inp_image         # Input fits image
    p_img_extnum  = p_args.image_extension   # Extension number for data in input fits
    p_srclist     = p_args.source_list       # Inout fits table containing srcs
    p_src_extname = p_args.xy_extension      # Extension number for data in input fits
    p_outimg      = p_args.out_image         # Name of output fits image
    p_astnetkey   = p_args.key               # Astrometry.net API key or None
    p_use_sip     = p_args.use_sip           # Allow SIP fit? 
    p_user_scale  = p_args.user_scale        # User specified pixel scale (arcsec)
    p_scale_err_ratio = p_args.scale_err_ratio # Scale error ratio

    ap_astrom = ap.ApAstrometry(p_inpimg, 
        p_img_extnum,
        p_srclist, 
        p_src_extname,
        p_outimg, 
        p_astnetkey,
        p_use_sip,
        p_user_scale,
        p_scale_err_ratio,
        p_loglevel)
    return ap_astrom.status()

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
