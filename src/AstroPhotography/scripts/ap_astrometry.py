#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ap_astrometry.py
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
#  

# 2020-08-12 dks : Initial WIP

import argparse
import sys
import logging
import numpy as np
import math
import os.path

from astropy.io import fits
from astropy.table import QTable, Table
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy import units as u

from astroquery.astrometry_net import AstrometryNet
    
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
        required=True,
        type=str,
        metavar='ASTROMETRY_API_KEY',
        help='Your personal Astrometry.net API key is required.')

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
                
    args = parser.parse_args(argv)
    return args

class ApAstrometry:
    """Uses astrometry.net and a list of valid stars to calculate an
       astrometric solution to a given FITS image, writing a copy of
       the image with valid WCS keywords.
    """
    
    def __init__(self, inp_img_fname, 
        inp_img_extnum,
        srclist_fname, 
        srclist_extname,
        out_img_fname, 
        astnet_key,
        loglevel):
        
        self._initialize_logger(loglevel)
        self._loglevel = loglevel
        
        # Read input image
        self._inp_fname = inp_img_fname
        self._idata, self._ihdr = self._read_fits_image(inp_img_fname, 
            inp_img_extnum)
            
        # Read star position data 
        self._src_fname = srclist_fname
        self._src_xytable, self._src_meta = self._read_srclist(srclist_fname,
            srclist_extname)
            
        # Check input image corresponds to the image the source list was
        # generated from.
        # TODO
        
        # Build additional astrometry settings based on source metadata
        # if possible. This should improve the speed of the solution.
        # TODO
        
        # Run query 
        # TODO
        self._out_fname  = out_img_fname
        self._astnet_key = astnet_key
        
        # If a WCS solution was found then create a copy of the input
        # image and add parts of the WCS header data to it.
        # TODO
        return
    
    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        self._logger = logging.getLogger(__name__)
        
        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: {}'.format(loglevel))
        self._logger.setLevel(numeric_level)
    
        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(numeric_level)
    
        # create formatter
        formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s | %(message)s')
    
        # add formatter to ch
        ch.setFormatter(formatter)
    
        # add ch to logger
        self._logger.addHandler(ch)
        return
        
    def _check_file_exists(self, filename):
        if not os.path.isfile(filename):
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
    
    def _read_fits_image(self, image_filename, image_extension):
        """Read a single extension's data and header from a FITS file
        """
            
        self._check_file_exists(image_filename)
        self._logger.info('Loading extension {} of FITS file {}'.format(image_extension, image_filename))
            
        # open() parameters that can be important.
        # Default values used here.
        # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
        uint_handling = True
        image_scaling = False
            
        with fits.open(image_filename, 
            uint=uint_handling, 
            do_not_scale_image_data=image_scaling) as hdu_list:
            ext_hdr  = hdu_list[image_extension].header
            ext_data = hdu_list[image_extension].data
            
        ndim     = ext_hdr['NAXIS']
        cols     = ext_hdr['NAXIS1']
        rows     = ext_hdr['NAXIS2']
        bitpix   = ext_hdr['BITPIX']
        bzero    = ext_hdr['BZERO']
        bscale   = ext_hdr['BSCALE']
        info_str = '{}-D BITPIX={} image with {} columns, {} rows, BSCALE={}, BZERO={}'.format(ndim, bitpix, cols, rows, bscale, bzero)
        
        if ndim == 3:
            layers = ext_hdr['NAXIS3']
            info_str = '{}-D BITPIX={} image with {} columns, {} rows, {} layers, BSCALE={}, BZERO={}'.format(ndim, bitpix, cols, rows, layers, bscale, bzero)
            
        self._logger.debug(info_str)
        if ndim == 3:
            self._logger.error('Error, 3-D handling has not been implemented yet.')
            sys.exit(1)
        return ext_data, ext_hdr
        
    def _read_srclist(self, srclist_fname, srclist_extname):
        """Read the XY positions from a source list generated by ap_find_stars
        """

        self._check_file_exists(srclist_fname)
         
        # Open the source list
        with fits.open(srclist_fname) as hdu_list:
            # Read keywords from the primary header
            srclist_hdr = hdu_list[0].header
            
            # Search for the srclist_extname and read the data.
            if srclist_extname in hdu_list:
                xy_table= Table.read(hdu_list[srclist_extname])
            else:
                err_msg = f'{filename} does not contain the expected {srclist_extname} extension.'
                self._logger.error(err_msg)
                raise RuntimeError(err_msg)

        num_srcs = len(xy_table)
        self._logger.debug(f'Read {num_srcs} source positions from {srclist_fname}')
        return xy_table, srclist_hdr
        
    def _query_astrometry_dot_net(self, astnetkey, cols, rows, xy_table):
        """Attempt get astrometric solution"""
    
        wcs = {}
        image_width  = cols
        image_height = rows
        
        ast = AstrometryNet()
        ast.api_key = p_astnetkey
        try_again = True
        submission_id = None
    
        while try_again:
            try:
                if not submission_id:
                    SELF._logger.debug('Submitting astrometry.net solve from source list with {} poisitions'.format(len(xy_table)))
                    wcs_header = ast.solve_from_source_list(xy_table['X'], 
                        xy_table['Y'],
                        image_width, 
                        image_height,
                        parity=2,
                        positional_error=10,
                        crpix_center=True,
                        publicly_visible='n',
                        submission_id=submission_id)
                else:
                    self._logger.debug('Monitoring astrometry.net submission {}'.format(submission_id))
                    wcs_header = ast.monitor_submission(submission_id,
                        solve_timeout=120)
            except TimeoutError as e:
                submission_id = e.args[1]
            else:
                # got a result, so terminate
                try_again = False
        
        if wcs_header:
            self._logger.info('Obtained a WCS header')
            wcs = wcs_header
        else:
            self._logger.error('Astrometry.net submission={} failed'.format(submission_id))
    
        return wcs
    
    def create_primary_header(kw_dict):
        """Creates a FITS primary header adding the keywords in the dict"""
        
        pri_hdr = fits.Header()
        for kw in kw_dict:
            pri_hdr[kw] = kw_dict[kw]
        pri_hdr['HISTORY'] = 'Created by ApFindStars'
        return pri_hdr
        
    def add_optional_keywords(hdr, kw_dict):
        """Add keywords to kw_dict from the FITS header hdr if present"""
        logger = logging.getLogger(__name__)
        
        # The following may exist.
        kw_comment_dict = {'EXPOSURE': '[seconds] Image exposure time',
            'DATE-OBS': 'Observation date and time',
            'OBJECT':   'Target object', 
            'TELESCOP': 'Telescope used',
            'RA':       'Requested right ascension', 
            'DEC':      'Requested declination', 
            'XPIXSZ':   '[micrometers] Stated X-axis pixel scale after binning',
            'YPIXSZ':   '[micrometers] Stated Y-axis pixel scale after binning',
            'FOCALLEN': '[mm] Stated telescope focal length', 
            'FILTER':   'Filter used', 
            'EGAIN':    '[e/ADU] Gain in electrons per ADU'}
        kw_missing_list = []
        for kw in kw_comment_dict:
            if kw in hdr:
                kw_dict[kw] = (hdr[kw], kw_comment_dict[kw])
            else:
                kw_missing_list.append(kw)
    
        logger.debug('FITS keywords found in image: {}'.format(kw_dict))
        logger.debug('FITS keywords missing from image: {}'.format(kw_missing_list))
        return
    

def main(args=None):
    p_args        = command_line_opts(args)
    p_loglevel    = p_args.loglevel          # Loglevel
    p_inpimg      = p_args.inp_image         # Input fits image
    p_img_extnum  = p_args.image_extension   # Extension number for data in input fits
    p_srclist     = p_args.source_list       # Inout fits table containing srcs
    p_src_extname = p_args.xy_extension      # Extension number for data in input fits
    p_outimg      = p_args.out_image         # Name of output fits image
    p_astnetkey   = p_args.key               # Astrometry.net API key or None

    ap_astrom = ApAstrometry(p_inpimg, 
        p_img_extnum,
        p_srclist, 
        p_src_extname,
        p_outimg, 
        p_astnetkey,
        p_loglevel)
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
