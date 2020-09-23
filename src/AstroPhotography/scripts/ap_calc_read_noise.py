#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ap_calc_read_noise.py
#
#  Calculates an estimate of the detector read noise (electrons/pixel)
#  given two raw bias images and the electronic gain (electrons/ADU).
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
#  2020-09-16 dks : Initial coding. 
# 

import argparse
import sys
import logging
import os.path
import math

import numpy as np
from astropy.io import fits

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_calc_rread_noise',
        description='Calculates an estimate of the detector read noise' + 
        ' (e/pixel) given two raw bias images and the electronic gain' +
        ' (e/ADU).')
        
    # Required
    parser.add_argument('biasfile1',
        metavar='RAW_BIAS_1.FITS',
        help='Path/name of the first raw bias file to read.')
    parser.add_argument('biasfile2',
        metavar='RAW_BIAS_2.FITS',
        help='Path/name of the second raw bias file to read.')
        
    # Defaults.
    p_gain = 'EGAIN'
        
    # Optional
    parser.add_argument('--gain',
        default=p_gain,
        help=('FITS keyword listing the electronic gain (electrons/ADU)'
        ' within the first bias file, OR the numerical value of the'
        f' gain to use. Default: {p_gain}'))
    parser.add_argument('--noclip',
        dest='sigmaclip',
        action='store_false',
        default=True,
        help='If specified, extreme values of the bias pixel-to-pixel' +
        ' differences will NOT be removed using sigma clipping.')
    parser.add_argument('--histplot',
        default=None,
        help='If specified, the path/name to use for a histogram plot' +
        ' that this script will generate of the pixel-to-pixel' +
        ' differences between the two bias files.')
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    return args

class ApImageDifference:
    """Takes the difference between two images of the same shape, with
       optional internal sigma clipping-based masking or using externally
       supplied mask arrays. 
       
    The difference array, masks, and difference statistics can be 
    returned using accessors. An optional plot of the difference 
    histogram can be generated.
    """
    
    def __init__(self,
        loglevel):
    
        # Initialize logging
        self._loglevel = loglevel
        self._initialize_logger(self._loglevel)
        
        return
                
    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        self._logger = logging.getLogger('ApImageDifference')
        
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
        
class ApCalcReadNoise:
    """Calculates an estimate of detector read noise (e/pixel) given
       two raw bias files and the electronic gain (e/ADU).
       
    Given two bias frames, B1 and B2, and the difference between them
    B1-B2, the standard deviation of the difference image is related
    to the read noise of the detector by
    
    Read Noise = Gain * (sigma_(b1-b2)) / sqrt(2)
       
    See section 4.3 of "Handbook of CCD Astronomy", Howell, S.B.,
    2000 (Cambridge University Press, Cambridge UK)
    """
    
    def __init__(self,
        biasfile1,      # Name/path of first bias file.
        biasfile2,      # Name/path of second bias file.
        gain,           # Gain value or gain FITS keyword to use.
        loglevel):
        """The constructor stores the names of the files and gain
           keyword or value, and initializes the logger, but does not
           perform any file IO or perform calculations.
           
           The actual file IO and calculations are only performed when
           the user calls estimate_rn.
        """
    
        # Initialize logging
        self._loglevel  = loglevel
        self._biasfile1 = biasfile1
        self._biasfile2 = biasfile2
        self._gaininfo  = gain
        self._initialize_logger(self._loglevel)
        self._logger.debug( ( 'Initialized an ApCalcReadNoise instance with'
            f' biasfile1={biasfile1}, biasfile1={biasfile2}, gain={gain},'
            f' and loglevel={loglevel}' ) )
        return
                
    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        self._logger = logging.getLogger('ApImageDifference')
        
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
        bzero    = ext_hdr['BZERO']
        bscale   = ext_hdr['BSCALE']
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
        
        
    def estimate_rn(self, sigmaclip, histplot=None):
        """Estimate the read noise from the files, and optionally
           produce a histograme plot of the pixel-to-pixel ADU
           differences. 
        
        If sigmaclip is true pixels with outlier values in either
        bias image will be excluded from the pixel-to-pixel differencing.
        Note that outliers are identified in each bias file separately
        in terms of absolute intensity, and the resulting bad pixel maps
        are combined.
        """
        
        # Read the files files, assuming data is in primary.
        data1, hdr1 = self._read_fits(self._biasfile1, 0)
        data2, hdr2 = self._read_fits(self._biasfile2, 0)
        if data1.shape != data2.shape:
            err_msg = ('Error, data array shapes do not match:'
                f' First file={data1.shape},'
                f' second file={data2.shape}')
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)
        
        # Work out the gain value to use (e/ADU)
        self._gain = self._select_gain(hdr1, hdr2)
        self._logger.info(f'Adopted gain is {self._gain:.2f} electrons/ADU.')
        
        # Calculate the difference with or without sigma clipping and
        # histogram plotting, get the standard deviation (ADU)
        stddev = 0 # TODO
        
        # Calculate read noise estimate in e/pixel
        read_noise = self._gain * stddev / math.sqrt(2)
        
        return read_noise
        
    def _isfloat(self, value_str):
        """Returns True if value_str can be converted to a numeric float
           value.
        
        Will also return True for booleans and NaNs.
        """
        
        try:
            float(value_str)
            return True
        except ValueError:
            return False

    def _select_gain(self, hdr1, hdr2):
        """Determines what gain value to use, based on the two file headers
           and the gaininfo supplied at object construction.
           
        At construction a gain value (string representation of a floating
        point number) or the name of a FITS primary header keyword was
        supplied, and is stored as _gaininfo.
        - If a string representation of a floating point number was 
          supplied then we should use that.
        - If _gaininfo is not convertable to float, treat it as a FITS
          keyword.
          - If the keyword not present in both FITS headers then this
            is treated as an error.
          - If the value associated with the keyword is different beyond
            a numerical tolerance this is an error.
        
        If no gain information can be found then this function raises
        an exception.
        """
        
        if self._isfloat(self._gaininfo):
            # Use the value given instead of looking in the FITS files.
            return float(self._gaininfo)
        
        # Otherwise treat _gaininfo as a FITS header keyword
        
        return
        
                
def main(args=None):
    p_args      = command_line_opts(args)
    p_biasfile1 = p_args.biasfile1
    p_biasfile2 = p_args.biasfile2
    p_gain      = p_args.gain
    p_sigmaclip = p_args.sigmaclip
    p_histplot  = p_args.histplot
    p_loglevel  = p_args.loglevel
    
    read_noise_calculator = ApCalcReadNoise(p_biasfile1,
        p_biasfile2,
        p_gain,
        p_loglevel)
    
    rn1 = read_noise_calculator.estimate_rn(p_sigmaclip, p_histplot)
    print(f'Estimated read noise is {rn1:.2f} electrons/pixel.')
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
