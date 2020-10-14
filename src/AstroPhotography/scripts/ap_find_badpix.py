#!/usr/bin/env python
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

import argparse
import sys
import logging
import os.path
import math

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

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

class ApFindBadPixels:
    """A class used to find bad pixels within dark or bias files based on
       deviation from an expected uniform mean or median value. The instance
       may be queried for properties of the input file, good or bad pixel
       numbers or count values, the bad pixel map can be extracted as
       a numpy array or written to a fits file.
    """
    
    def __init__(self,
        darkfile,
        sigma,
        loglevel):
        """Constructs an ApFindBadPixels object and performs preliminary
           processing on it.
        
        :param darkfile: Input dark or bias file to search for bad pixels.
        :param sigma: Number of standard deviations away from the clipped
          median a pixel must be (or more) to count as a bad pixel.
        :param loglevel: Logging level to use.
        """
    
        # Initialize logging
        self._loglevel = loglevel
        self._initialize_logger(self._loglevel)
        
        # Data file to read
        self._imfile   = darkfile
        self._imextnum = 0
                
        # Number of MAD std deviations for sigma clipping
        self._sigma = sigma

        # Process data.
        self._imdata, self._imhdr = self._read_fits(darkfile, 0)
        self._generate_sigmaclip_mask(self._imdata, self._sigma)
        return
        
    def _check_file_exists(self, filename):
        if not os.path.isfile(filename):
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return

    def _generate_sigmaclip_mask(self, data, sigma):
        """Creates a bad pixel mask based on sigma-clipped statistics
           of the input data array.
           
        This routine is most appropriate for images that expected to
        be relatively uniform, but with a small number of highly 
        discrepant values.
        
        The generated mask is True when pixels are BAD and False where
        pixels are GOOD. 
                   
        :param data: Input data array
        :param sigma: Number of sigma away from the median to consider
          a pixel bad.
        """
        
        self._logger.debug(f'Generating a bad pixel mask using sigma={sigma} clipping on the input image data values.')
        
        # Clip first input image
        mean, med, std   = sigma_clipped_stats(data, sigma=sigma)
        lothresh         = med - (sigma * std)
        hithresh         = med + (sigma * std)
        bad_lo           = data < lothresh
        bad_hi           = data > hithresh
        self._badpixmask = np.logical_or(bad_lo, bad_hi)
        self._logger.debug(f'After sigma clipping good pixels have values between {lothresh:.2f} and {hithresh:.2f} ADU.')
        return

    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        self._logger = logging.getLogger('ApFindBadPixels')
        
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
            
    def _read_fits(self, image_filename, image_extension):
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
        if 'BZERO' in ext_hdr:
            bzero    = ext_hdr['BZERO']
        else:
            bzero = 0
        if 'BSCALE' in ext_hdr:
            bscale   = ext_hdr['BSCALE']
        else:
            bscale = 1.0
        info_str = '{}-D BITPIX={} image with {} columns, {} rows, BSCALE={}, BZERO={}'.format(ndim, bitpix, cols, rows, bscale, bzero)
        
        if ndim == 3:
            layers = ext_hdr['NAXIS3']
            info_str = '{}-D BITPIX={} image with {} columns, {} rows, {} layers, BSCALE={}, BZERO={}'.format(ndim, bitpix, cols, rows, layers, bscale, bzero)
            
        self._logger.debug(info_str)
        if ndim == 3:
            self._loggererror('Error, 3-D handling has not been implemented yet.')
            sys.exit(1)
            
        # Get data absolute limits.
        minval = np.amin(ext_data)
        maxval = np.amax(ext_data)
        medval = np.median(ext_data)
        self._logger.debug(f'Raw data statistics are min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}')
        
        # Is there a PEDESTAL value? MaximDL likes to add an offset, and
        # the PEDESTAL value is the value to ADD to the data to remove the
        # pedestal.
        if 'PEDESTAL' in ext_hdr:
            pedestal = float( ext_hdr['PEDESTAL'] )
            if pedestal != 0:
                self._logger.debug(f'Removing a PEDESTAL value of {pedestal} ADU.')
                ext_data += pedestal
                minval = np.amin(ext_data)
                maxval = np.amax(ext_data)
                medval = np.median(ext_data)
                self._logger.debug(f'After PEDESTAL removal, min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}')
        
        return ext_data, ext_hdr
                
def main(args=None):
    p_args      = command_line_opts(args)
    p_inmstrdrk = p_args.masterdark
    p_outbadpix = p_args.badpixfile
    p_sigma     = p_args.sigma
    p_loglevel  = p_args.loglevel
    
    mkbadpix = ApFindBadPixels(p_inmstrdrk,
        p_sigma,
        p_loglevel)
    #mkbadpix.write_mask(p_outbadpix)
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
