"""
Contains the implementation of the ApProcess class.
"""

# 2024-05-02 dks : Initial implementation.

import sys
import logging
from pathlib import Path
import math
import time
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits

# AstroPhotography includes    
from .. import __version__
from . import ApFixBadPixels
from .ApFixCosmicRays import ApFixCosmicRays as ApFixCosmicRays
# TODO: I have no idea why the ApFixBadPixels line works, but the same for
# ApFixCosmicRays did not, instead requiring the more verbose line shown.

class ApProcess:
    """
    Astronomical image processor that works on calibrated images to
    perform the following tasks:
    
    - 
    """
    
    
    def __init__(self, loglevel):
        """
        Initializes an ApProcess instance
           
        :param loglevel: Log level.
        """
        
        self._name     = 'ApProcess'
        self._version  = __version__
        self._loglevel           = loglevel
        self._initialize_logger(self._loglevel)
        
        return
        
    def _check_file_exists(self, filename):
        """
        Raises an exception if the name file does not exist.
        """
        
        if not Path(filename).exists():
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return

    def _img_stats(self, data, label, verbose):
        """
        Calculate and display some image statistics
        """
        
        minval  = np.nanmin(data)
        maxval  = np.nanmax(data)
        meanval = np.nanmean(data)
        
        # percentiles, 50th percentile is the median
        #         0    1    2    3   4   5   6   7   8   9   10
        ipctls = [0.1, 1.0, 5.0, 10, 25, 50, 75, 90, 95, 99, 99.9]
        opctls = np.nanpercentile(data, ipctls)
        medval = opctls[5]
        
        if verbose:
            self._logger.info(f'{label} data min={minval:.2f}, max={maxval:.2f}, mean={meanval:.2f}, median={medval:.2f} ADU.')
            self._logger.info(f'  90% of data between {opctls[2]:.2f} and {opctls[8]:.2f} ADU (5-95 precentiles)')
            self._logger.info(f'  98% of data between {opctls[1]:.2f} and {opctls[9]:.2f} ADU (1-99 precentiles)')
        return [minval, maxval, meanval, medval]

    def _initialize_logger(self, loglevel):
        """
        Initialize and return the logger
        """
        
        self._logger = logging.getLogger(self._name)
        
        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: {}'.format(loglevel))
        self._logger.setLevel(numeric_level)
        self._logger.propagate = False
            
        # check if handlers already present
        if not len(logger.handlers):
            # create console handler and set level to debug
            ch = logging.StreamHandler()
            ch.setLevel(numeric_level)
        
            # create formatter
            formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s | %(message)s')
        
            # add formatter to ch
            ch.setFormatter(formatter)
        
            # add ch to logger
            logger.addHandler(ch)
        return
        
    def _read_fits(self, image_filename, image_extension):
        """
        Read a single extension's data and header from a FITS file
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
        info_str = '{}-D BITPIX={} image with {} columns, {} rows'.format(ndim, bitpix, cols, rows)
        
        if ndim == 3:
            layers = ext_hdr['NAXIS3']
            info_str = '{}-D BITPIX={} image with {} columns, {} rows, {} layers'.format(ndim, bitpix, cols, rows, layers)

        if 'BSCALE' in ext_hdr:
            bscale   = ext_hdr['BSCALE']
            info_str += f', BSCALE={bscale}'
        
        if 'BZERO' in ext_hdr:
            bzero    = ext_hdr['BZERO']
            info_str += f', BZERO={bzero}'

            
        self._logger.debug(info_str)
        if ndim == 3:
            self._logger.error('Error, 3-D handling has not been implemented yet.')
            sys.exit(1)
            
        # Convert to 32-bit floating point if necessary
        if not np.issubdtype(ext_data.dtype, np.floating):
            orig_dtype = ext_data.dtype
            ext_data   = ext_data.astype(np.float32)
            self._logger.debug(f'  Converted data type from {orig_dtype} to float32')
            
        # Get data absolute limits.
        minval = np.nanmin(ext_data)
        maxval = np.nanmax(ext_data)
        medval = np.nanmedian(ext_data)
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

    def _remove_pedestal_kw(self, hdr):
        """
        Removes the PEDESTAL keyword from the input FITS header
        
        AstroPhotography always removes any artificial PEDESTAL applied
        to the data when reading a FITS file, so it is important to make
        sure that the FITS header keywords remain consistent.
        
        This function need only be applied when modified data is being
        written or rewritten to disk using a copy of an original FITS
        header.
        """
        
        if 'PEDESTAL' in hdr:
            self._logger.debug('Removing PEDESTAL keyword from FITS header.')
            del hdr['PEDESTAL']
        
        return 

    def _write_corrected_image(self, inpdata_file,
            ext_num,
            outdata_file,
            odata, 
            odict):
        """
        Writes the calibrated data to the specified output
        file, preserving all other items from the original input file.
           
        The output file differs from the original input data file in 
        having the bias, dark and flat corrected data in floating point
        format, and a new/updated set of FITS header keywords
        from the input dictionary.
        
        :param inpdata_file: Input FITS data file affected by bad pixels.
        :param ext_num: Extension number for data array and header.
        :param outdata_file: Bias/dark/flat corrected image.
        :param odata: Bias/dark/flat field corrected data array.
        :param odict: Additional FITS header keywords. The output file
          the original keywords from the input FITS image file, plus
          the keywords in this dictionary.
        """
        
        self._logger.debug(f'FITS header keywords added to output: {odict}')
        self._check_file_exists(inpdata_file)
            
        # open() parameters that can be important.
        # Default values used here.
        # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
        uint_handling = True
        image_scaling = False
            
        with fits.open(inpdata_file, 
            uint=uint_handling, 
            do_not_scale_image_data=image_scaling) as hdu_list:
            
            self._remove_pedestal_kw(hdu_list[ext_num].header)
            for kw in ['BSCALE', 'BZERO']:
               if kw in hdu_list[ext_num].header:
                   del hdu_list[ext_num].header[kw]
            
            # Modify data
            hdu_list[ext_num].data = odata
            
            # Modify header
            for kw, val in odict.items():
                hdu_list[ext_num].header[kw] = val
            
            tnow = datetime.now().isoformat(timespec='milliseconds')
            hdu_list[ext_num].header['HISTORY'] = f'Processed by {self._name} {self._version} at {tnow}'
            
            # Write to new file
            hdu_list.writeto(outdata_file, 
                output_verify='ignore',
                overwrite=True)
        
        self._logger.info(f'Wrote bias/dark/flat corrected file to {outdata_file}')
        return
        
