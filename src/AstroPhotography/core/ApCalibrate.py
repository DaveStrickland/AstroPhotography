"""Contains the implementation of the ApCalibrate class.
"""

#  2020-12-05 dks : Initial implementation.

import sys
import logging
from pathlib import Path
import math
import time
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits

from .. import __version__
from . import ApFixBadPixels

class ApCalibrate:
    """Astronomical CCD image calibrator that performs bias subtraction,
       dark subtraction, flat fielding and bad pixel correction.
       
       The bias, dark and flat fielding methods exactly reproduce
       ccdproc's treatment on a limited number of test images. To 
       potentially gain greater control over the flat fielding I've 
       chosen to use this implementation instead of delegating to
       ccdproc.
    """
    
    # Flat fielding method enumeration
    MEAN_FULL   = 0  # Mean value of entire flat field.
    MEDIAN_FULL = 1  # Median value of entire flat field. TBA
    
    def __init__(self, loglevel):
        """Initializes an ApCalibrate instance.
        """
        
        self._name     = 'ApCalibrate'
        self._version  = __version__
        self._loglevel = loglevel
        self._initialize_logger(self._loglevel)
        return
        
    def _check_file_exists(self, filename):
        """Raises an exception if the name file does not exist.
        """
        
        if not Path(filename).exists():
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return

    def _find_exptime_ratio(self, img_hdr, dark_hdr):
        """Calculate the ratio of the image to dark exposure times.
        
        This function searches for 'EXPOSURE' or 'EXPTIME' to denote
        the exposure times, and assumes the values are in units of
        seconds. 
        """
        
        img_exp  = None
        dark_exp = None
        for kw in ['EXPOSURE', 'EXPTIME']:
            if img_exp is None:
                if kw in img_hdr:
                    img_exp = float( img_hdr[kw] )
                    self._logger.debug(f'Image exposure time [seconds]: {img_exp:.2f}')
            if dark_exp is None:
                if kw in dark_hdr:
                    dark_exp = float( dark_hdr[kw] )
                    self._logger.debug(f'Dark exposure time [seconds]: {dark_exp:.2f}')

        # Throw error if exposure times could not be found.
        msg = None
        if (img_exp is None) and (dark_exp is None):
            msg = 'Could not determine exposure time for both image and dark.'
        elif (img_exp is None):
            msg = 'Could not determine exposure time for image (dark exposure found).'
        elif (dark_exp is None):
            msg = 'Could not determine exposure time for dark (img exposure found).'
        if msg is not None:
            self._logger.error(msg)
            raise RunTimeError(msg)

        # Calculate exposure ratio
        exp_ratio = img_exp / dark_exp
        
        self._logger.info(f'Image to dark exposure time ratio: {exp_ratio:.3f}')
        return exp_ratio

    def _generate_flat(self, flat_data, flat_method):
        """Generate a normalized flat field given the flat field image
           and a method
           
        Currently this only supports MEAN_FULL, the same method used by ccdproc,
        normalization by the mean value of the entire flat field image.
           
        :param flat_data: Input array-like flat field image.
        :param flat_method: Method to be used to normalize the input
          flat_data.
        """
        
        norm_factor = 1
        if flat_method == byhand.MEAN_FULL:
            self._logger.debug('Using mean value of input field image to normalize by.')
            norm_factor = np.nanmean(flat_data)
        else:
            msg = f'Error, flat field normalization method {flat_method} has not been implemented yet.'
            self._logger.error(msg)
            raise RunTimeError(msg)
            
        self._logger.info(f'Flat field normalization factor: {norm_factor:.2f}')
        out_flat = flat_data / norm_factor
        [meanval, medval, minval, maxval] = img_stats(out_flat, 'Normalized flat', True)
        return out_flat

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
        """Removes the PEDESTAL keyword from the input FITS header
        
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
        """Writes the calibrated data to the specified output
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
            hdu_list[ext_num].header['HISTORY'] = f'Generated by {self._name} {self._version} at {tnow}'
            
            # Write to new file
            hdu_list.writeto(outdata_file, 
                output_verify='ignore',
                overwrite=True)
        
        self._logger.info(f'Wrote bias/dark/flat corrected file to {outdata_file}')
        return
        
    def calibrate(self, raw_image,
        master_bias,
        master_dark,
        master_flat,
        byhand_image,
        byhand_flat):
        """Perform bias subtraction, dark subtraction, flat fielding 
           (optional), and bad pixel correction (optional).
           
        """
        
        # Convert to path objects
        raw_image   = Path(raw_image)
        master_bias = Path(master_bias)
        master_dark = Path(master_dark)
        master_flat = Path(master_flat)
        
        ext_num = 0
        raw_data,  raw_hdr  = self._read_fits(raw_image,   ext_num)
        bias_data, bias_hdr = self._read_fits(master_bias, ext_num)
        dark_data, dark_hdr = self._read_fits(master_dark, ext_num)
        flat_data, flat_hdr = self._read_fits(master_flat, ext_num)
        
        # Subtract bias from raw image and from dark.
        # Skipping the normal sanity checks I would do.
        img_sub_b  = raw_data - bias_data
        dark_sub_b = dark_data - bias_data
        
        # Scale bias-subtracted dark by ratio of exposure times.
        exp_ratio   = self._find_exptime_ratio(raw_hdr, dark_hdr)
        dark_scaled = exp_ratio * dark_sub_b
        img_sub_bd  = img_sub_b - dark_scaled 
        
        # Flat field correction
        flat_method = byhand.MEAN_FULL
        norm_flat = self._generate_flat(flat_data, flat_method)
        img_bdf   = np.where(norm_flat != 0,
            img_sub_bd / norm_flat,
            img_sub_bd)
        
        # updated keywords
        odict = {'BIASCORR': (True, 'True if bias subtracted.'),
            'BIASFILE': (master_bias.name, 'Master bias file used.'),
            'DARKCORR': (True, 'True if scaled dark subtracted.'),
            'DARKFILE': (master_dark.name, 'Master dark file used.'),
            'FLATCORR': (True, 'True if flat field applied.'),
            'FLATFILE': (master_flat.name, 'Master flat file used.'),
            'BUNIT':    ('adu', 'Pixel value units.')}
        self._write_corrected_image(raw_image, ext_num,
            byhand_image, img_bdf, odict)
            
        self._write_corrected_image(raw_image, ext_num,
            byhand_flat, norm_flat, {})
        
        return
