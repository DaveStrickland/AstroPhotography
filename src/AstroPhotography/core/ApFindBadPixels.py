"""Contains the implementation of the ApFindBadPixels class.
"""

# 2020-12-31 dks : Moved ApFindBadPixels into core from ap_find_badpix.py

import logging
import os.path
import math
from datetime import datetime, timezone
 
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from .. import __version__

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
    
        self._name = 'ApFindBadPixels'
    
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
        self._image_stats()
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
        
        npix = self._imdata.size
        self._logger.debug(f'Generating a bad pixel mask using sigma={sigma} clipping on the input image data values.')
        
        # Compute sigma clipped statistics
        mean, med, std   = sigma_clipped_stats(data, sigma=sigma)
        self._logger.debug(f'Sigma-clipped mean={mean:.2f}, median={med:.2f}, and madstddev={std:.2f} values (ADU).')
        
        lothresh         = med - (sigma * std)
        hithresh         = med + (sigma * std)
        self._logger.info(f'Good pixels have values between {lothresh:.2f} and {hithresh:.2f} ADU.')
        
        # Pixels below low threshold
        bad_lo           = data < lothresh
        nbad_lo          = np.sum(bad_lo)
        self._logger.debug(f'{nbad_lo} (out of {npix}) pixels below {lothresh:.2f} ADU.')
        
        # Pixels above high threshold
        bad_hi           = data > hithresh
        nbad_hi          = np.sum(bad_hi)
        self._logger.debug(f'{nbad_hi} (out of {npix}) pixels above {hithresh:.2f} ADU.')
        
        # Generate final mask
        self._badpixmask = np.logical_or(bad_lo, bad_hi)
        nbad             = np.sum(self._badpixmask)
        pct_bad          = 100 * (nbad/npix)
        msg              = f'Out of {npix} pixels, {nbad} are bad ({pct_bad:.4f}%).'
        self._logger.info(msg)
        return

    def _image_stats(self):
        """If logging level is DEBUG, print some statistics of the data.
        """
        
        if self._logger.getEffectiveLevel() == logging.DEBUG:
            minval = np.min(self._imdata)
            maxval = np.max(self._imdata)
            self._logger.debug(f'Data min={minval:.2f}, max={maxval:.2f} ADU.')
            
            # percentiles
            ipctls = [0.1, 1.0, 5.0, 10, 25, 50, 75, 90, 95, 99, 99.9]
            opctls = np.percentile(self._imdata, ipctls)
            for idx, pct in enumerate(ipctls):
                self._logger.debug(f'{pct:.2f}% of the data has value of {opctls[idx]:.2f} ADU or less.')
        
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
        
    def _update_header(self, hdu):
        """Updates the raw mask FITS primary header by adding select
           keywords from the input master dark/bias file.
           
        :param hdu: FITS hdu object to be modified.
        """
        
        self._logger.debug('Updating FITS primary HDU keywords.')
        
        # keyword dictionary to write to header
        kw_dict = {}
        
        # Copy select keywords: basically those that would identify the
        # telescope, instrument, and imaging mode used.
        copy_list = ['TELESCOP', 'INSTRUME', 'SET-TEMP', 'CCD-TEMP',
            'XPIXSZ', 'YPIXSZ', 'XBINNING', 'YBINNING',
            'XORGSUBF', 'YORGSUBF', 'SITELAT', 'SITELONG']
        
        imgtype          = 'BADPIX'
        tnow             = datetime.now().isoformat(timespec='milliseconds')
        creation_date    = datetime.now(timezone.utc)
        creation_datestr = creation_date.isoformat(timespec='seconds')

        kw_dict['IMAGETYP'] = (imgtype, 'Type of file')
        kw_dict['CREATOR']  = (self._name, 'Software that generated this file.')
        kw_dict['DATE']     = (creation_datestr, 'UTC creation time.')
        kw_dict['DATAFILE'] = (self._imfile, 'Data file used to identify bad pixels.')
    
        for kw in copy_list:
            if kw in self._imhdr:
                val = self._imhdr[kw]
                cmt = self._imhdr.comments[kw]
                kw_dict[kw] = (val, cmt)
            
        for kw in kw_dict:
            hdu.header[kw] = kw_dict[kw]
        hdu.header['HISTORY'] = f'Created by {self._name} {__version__} at {tnow}'
        return
    
    def write_mask(self, mask_file_name):
        """Write the bad pixel mask to a FITS file with the user 
           specified name/path.
           
        :param mask_file_name: File name/path for the output bad pixel
          mask. The file will be overwritten if it exists.
        """
        
        mask = self._badpixmask.astype('uint8').copy()
        hdu  = fits.PrimaryHDU(data=mask)
        self._update_header(hdu)
        
        hdu_list = fits.HDUList([hdu])
        hdu_list.writeto(mask_file_name, overwrite=True)
        self._logger.info(f'Wrote bad pixel mask to {mask_file_name}')
        return
    
