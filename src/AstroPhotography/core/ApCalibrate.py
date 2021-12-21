"""Contains the implementation of the ApCalibrate class.

    TODO: Allow calibration of images in memory, in addition to files.
"""

# 2020-12-05 dks : Initial implementation.
# 2020-12-08 dks : Working version.
# 2021-01-14 dks : Added calls to ApFixCosmicRays
# 2021-01-26 dks : Got the syntax for importing ApFixCosmicRays right.
# 2021-02-27 dks : Add dark_still_biased flag.

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
from .ApFixCosmicRays import ApFixCosmicRays as ApFixCosmicRays
# TODO: I have no idea why the ApFixBadPixels line works, but the same for
# ApFixCosmicRays did not, instead requiring the more verbose line shown.

class ApCalibrate:
    """Astronomical CCD image calibrator that performs bias subtraction,
       dark subtraction, flat fielding and bad pixel correction.
       
       The bias, dark and flat fielding methods exactly reproduce
       ccdproc's treatment on a limited number of test images. To 
       potentially gain greater control over the flat fielding I've 
       chosen to use this implementation instead of delegating to
       ccdproc.
       
       The calibrator object applies to a specific telescope/detector
       (temperature) combination via the master bias, master dark 
       and optional master bad pixel files. If a master flat is supplied
       then it is specific to a given filter too.
    """
    
    # Flat fielding method enumeration
    MEAN_FULL   = 0  # Mean value of entire flat field.
    MEDIAN_FULL = 1  # Median value of entire flat field. TBA
    
    def __init__(self, master_bias_file,
        master_dark_file,
        master_flat_file,
        master_badpix_file,
        loglevel,
        dark_still_biased=None):
        """Initializes an ApCalibrate instance with the master 
           calibration files.
           
        :param master_bias_file: Name/path of master bias.
        :param master_dark_file: Name/path of master dark.
        :param master_flat_file: Name/path of master flat (filter specific)
        :param master_badpix_file: Name/path of master bad pixel file.
        :param loglevel: Log level.
        :param dark_still_biased: True or False, or None.
            Use this flag to specify that the master dark file has NOT
            already had bias subtraction performed, and that ApCalibrate
            should bias-subtract the dark before scaling it by the
            exposure time ratio. By default the software assumes that
            the bias has ALREADY been subtracted from the master dark.
            iTelescope master darks with CALSTAT=M are still biased and
            require the use of this flag. Files with CALSTAT=BM have
            had bias subtraction and this flag can be omitted.
        """
        
        self._name     = 'ApCalibrate'
        self._version  = __version__
        self._master_bias_file   = master_bias_file
        self._master_dark_file   = master_dark_file
        self._master_flat_file   = master_flat_file
        self._master_badpix_file = master_badpix_file
        self._loglevel           = loglevel
        self._initialize_logger(self._loglevel)
        
        self._master_bias = Path(self._master_bias_file)
        self._master_dark = Path(self._master_dark_file)

        if dark_still_biased is not None:
            self._dark_still_biased = dark_still_biased
        else:
            self._dark_still_biased = False

        # Helpers.
        self._bpix  = None # No ApFixBadPixels by default.
        self._crfix = ApFixCosmicRays(self._loglevel)

        # Read bias and dark data
        # TODO check bias exposure time is zero.
        ext_num = 0
        self._bias_data, self._bias_hdr = self._read_fits(self._master_bias, ext_num)
        self._dark_data, self._dark_hdr = self._read_fits(self._master_dark, ext_num)

        # Read optional flat
        if self._master_flat_file is not None:
            self._master_flat = Path(self._master_flat_file)
            self._flat_method = ApCalibrate.MEAN_FULL
            
            self._logger.info(f'Reading master flat field {self._master_flat.name}')
            flat_data, flat_hdr = self._read_fits(self._master_flat, ext_num)
            self._norm_flat     = self._generate_flat(flat_data, self._flat_method)

        # Read optional bad pixel file
        if self._master_badpix_file is not None:
            self._bpix  = ApFixBadPixels(self._loglevel)
            self._master_bpix = Path(self._master_badpix_file)
            self._mskdata, self._mskhdr = self._read_fits(self._master_bpix, ext_num)
        else:
            self._bpix = None
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
        if flat_method == ApCalibrate.MEAN_FULL:
            self._logger.debug('Using mean value of input field image to normalize by.')
            norm_factor = np.nanmean(flat_data)
        else:
            msg = f'Error, flat field normalization method {flat_method} has not been implemented yet.'
            self._logger.error(msg)
            raise RunTimeError(msg)
            
        self._logger.info(f'Flat field normalization factor: {norm_factor:.2f}')
        out_flat = flat_data / norm_factor
        [meanval, medval, minval, maxval] = self._img_stats(out_flat, 'Normalized flat', True)
        return out_flat

    def _get_gain(self, hdr):
        """Try to find the gain (e/ADU) from the specified header
        
        If the gain cannot be found a value of 1.0 will be returned. 
        """
        
        # Determine gain
        gain = None
        for kw in ['GAIN', 'EGAIN']:
            if kw in hdr:
                gain = float( hdr[kw] )
                self._logger.debug(f'Read gain value of {gain:.3f} e/ADU from {kw} keyword.')
        if gain is None:
            gain = 1.0
            self._logger.warning(f'Could not find gain value in header. Assuming gain={gain:.3f} e/ADU.')
        
        return gain

    def _img_stats(self, data, label, verbose):
        """Calculate and display some image statistics"""
        
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
            hdu_list[ext_num].header['HISTORY'] = f'Processed by {self._name} {self._version} at {tnow}'
            
            # Write to new file
            hdu_list.writeto(outdata_file, 
                output_verify='ignore',
                overwrite=True)
        
        self._logger.info(f'Wrote bias/dark/flat corrected file to {outdata_file}')
        return
        
    def calibrate(self, raw_image,
        cal_image,
        delta_pix,
        norm_flat,
        fixcosmic):
        """Perform bias subtraction, dark subtraction, flat fielding 
           (optional), and bad pixel correction (optional) on the input
           image.
           
        :param raw_image: File name/path for raw light frame to be 
          processed.
        :param cal_image: File name/path for output calibrated image
          file that has the requested calibration stages applied.
          The file will be over written if it already exists.
        :param delta_pix: Linear distance away from a bad pixel from which
          the median value of the good pixels will be drawn. See
          ApFixBadPixels for additional explanation.
        :param norm_flat: If not None, then a FITS image with the 
          normalized field will be written to this file name/path.
        :param fixcosmic: Boolean that controls whether cosmic ray
          corrections using ApFixCosmicRays will be performed.
        """
        
        # Performance timer
        perf_time_start = time.perf_counter() # highest res timer, counts sleeps
        
        # Convert to path objects
        raw_image   = Path(raw_image)        
        ext_num = 0
        raw_data,  raw_hdr  = self._read_fits(raw_image,   ext_num)
        
        # Subtract bias from raw image and from dark if still biased.
        # Skipping the normal sanity checks I would do.
        img_sub_b  = raw_data - self._bias_data
        if self._dark_still_biased:
            self._logger.info('Subtracting bias from dark')
            dark_sub_b = self._dark_data - self._bias_data
        else:
            self._logger.debug('Dark assumed to already be bias-subtracted.')
            dark_sub_b = self._dark_data
        
        # Scale bias-subtracted dark by ratio of exposure times.
        # TODO check that dark exposure time is >= image exposure time.
        exp_ratio   = self._find_exptime_ratio(raw_hdr, self._dark_hdr)
        dark_scaled = exp_ratio * dark_sub_b
        img_sub_bd  = img_sub_b - dark_scaled 
        
        # Dictionary of keywords to add to output...
        odict = {'BIASCORR': (True, 'True if bias subtracted.'),
            'BIASFILE': (self._master_bias.name, 'Master bias file used.'),
            'DARKCORR': (True, 'True if scaled dark subtracted.'),
            'DARKFILE': (self._master_dark.name, 'Master dark file used.'),
            'BUNIT':    ('adu', 'Pixel value units.')}

        # Flat field correction
        if self._master_flat_file is not None:
            img_bdf   = np.where(self._norm_flat != 0,
                img_sub_bd / self._norm_flat,
                img_sub_bd)
            odict['FLATCORR'] = (True, 'True if flat field applied.')
            odict['FLATFILE'] = (self._master_flat.name, 'Master flat file used.')
            if norm_flat is not None:
                # TODO, set keywords to be used in normalized flat?
                self._logger.debug(f'Writing normalized flat field to {norm_flat}')
                self._write_corrected_image(raw_image, ext_num,
                    norm_flat, self._norm_flat, {})
        else:
            self._logger.info('No flat field correction applied.')
            img_bdf = img_sub_bd
                    
        # Optional bad pixel correction.
        if self._bpix is not None:
            img_bdf, bpix_odict = self._bpix.fix_bad_pixels(img_bdf, 
                self._mskdata, delta_pix)
                
            # Update output FITS header keyword dictionary.
            odict['BPIXFILE'] = (self._master_bpix.name,
                'Name of master bad pixel file used')
            for key, val in bpix_odict.items():
                if 'BPIX' in key:
                    odict[key] = val
        else:
            self._logger.info('No bad pixel correction applied.')

        # Optional cosmic ray removal.
        if fixcosmic:
            self._logger.info('Correcting cosmic rays...')
            gain = self._get_gain(raw_hdr)
            img_clean, cr_kw = self._crfix.process(img_bdf, gain)
            for kw, val in cr_kw.items():
                odict[kw] = val
            img_bdf = img_clean

        # Performance 
        perf_time_end = time.perf_counter() # highest res timer, counts sleeps
        run_time_secs = perf_time_end - perf_time_start

        # Write final image.
        self._logger.info(f'Writing calibrated image to {cal_image}')    
        self._write_corrected_image(raw_image, ext_num,
            cal_image, img_bdf, odict)

        self._logger.info(f'Calibrated {raw_image.name} in {run_time_secs:.3f} seconds.')                    
        return
