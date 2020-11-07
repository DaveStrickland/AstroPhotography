"""Contains the implementation of the ApFixBadPixels class.
"""

import sys
import logging
import os.path
import math
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

class ApFixBadPixels:
    """A class used to fix pre-indentified bad pixels within an image
       by replacing them with the median value for surrounding good 
       pixels.
    """
    
    def __init__(self,
        loglevel):
        """Constructs an ApFixBadPixels object. No processing is performed.
        
        :param loglevel: Logging level to use.
        """
    
        self._name = 'ApFixBadPixels'
    
        # Initialize logging
        self._loglevel = loglevel
        self._initialize_logger(self._loglevel)
        
        self._logger.debug(f'{self._name} instance constructed.')
        return
        
    def _check_file_exists(self, filename):
        if not os.path.isfile(filename):
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return

    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        self._logger = logging.getLogger(self._name)
        
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
        
        # Used in cases where we get the same message twice or more
        # See https://stackoverflow.com/a/44426266
        self._logger.propagate = False
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

    def fix_files(self, inpdata_file, 
        badpixmask_file,
        outdata_file, 
        deltapix=1):
        """Fix bad pixels in the input FITS file based on the mask in an
           badpixel mask file, and write the updated file with a new 
           file name.
        
        :param inpdata_file: Input FITS data file affected by bad pixels.
        :param badpixmask_file: Input FITS file specifying the bad pixels.
          Bad pixels should have non-zero values, good pixels should have
          zero as the pixel value.
        :param outdata_file: Modified copy of the input data file where the
          bad pixels have had their data valus modified by the median
          of the surrounding good pixels.
        :param deltapix: Linear distance away from a bad pixel from which
          the median value of the good pixels will be drawn. If 1 then
          the median value of good pixels within the surrounding 8 pixels
          will be used. If 2 then the median of the good pixels within the
          surrounding 24 pixels will be used. Values above 2 are not
          recommended.
        """
        
        msg = (f'fix_files input data file={inpdata_file},'
            f' mask file={badpixmask_file},'
            f' output file={outdata_file}, deltapix={deltapix}')
        self._logger.info(msg)
        
        idata, ihdr     = self._read_fits(inpdata_file, 0)
        mskdata, mskhdr = self._read_fits(badpixmask_file, 0 )
        
        odata, odict = self.fix_bad_pixels(idata, mskdata, deltapix)
        return

    def fix_bad_pixels(self, data, badpixmask, deltapix=1):
        """Fix bad pixels in the input array based on the mask in an
           badpixel array, returning the modified data array, 
           and a 
           dictionary summarizing the number of bad pixels, pixels 
           corrected, and pixels than could not be corrected.
        
        :param data: Input ndarray representing the data affected by
          bad pixels.
        :param badpixmask: Input ndarray specifying the bad pixels.
          Bad pixels should have non-zero values, good pixels should have
          zero as the pixel value.
        :param deltapix: Linear distance away from a bad pixel from which
          the median value of the good pixels will be drawn. If 1 then
          the median value of good pixels within the surrounding 8 pixels
          will be used. If 2 then the median of the good pixels within the
          surrounding 24 pixels will be used. Values above 2 are not
          recommended.
        """
        
        # Info message about array shape and data type.
        self._logger.info(f'fix_bad_pixels: data has {data.shape[0]} rows x {data.shape[1]} columns, dtype={data.dtype}')
        self._logger.info(f'fix_bad_pixels: mask has {badpixmask.shape[0]} rows x {badpixmask.shape[1]} columns, dtype={badpixmask.dtype}')
        
        # Check that sizes match
        if (data.shape != badpixmask.shape):
            msg = (f'Error, the shape of the input data array ({data.shape})'
                f' does not match that of the bad pixel mask array ({badpixmask.shape}).')
            self._logger.error(msg)
            raise RunTimeError(msg)
        
        newdata     = data.copy()
        fixed_stats = {}
        return newdata, fixed_stats
