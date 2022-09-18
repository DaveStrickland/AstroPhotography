"""Contains the implementation of the ApFindBadPixels class.
"""

# 2020-12-31 dks : Moved ApFindBadPixels into core from ap_find_badpix.py
# 2020-01-09 dks : Added user-defined bad pixel processing. 

import logging
from pathlib import Path
import math
from datetime import datetime, timezone
 
import numpy as np
import yaml
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
    GOOD     = 0
    AUTO_BAD = 1
    USER_BAD = 2
    
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
        
        # Data file to read for automated bad pixel identification.
        self._imfile   = darkfile
        self._imextnum = 0
        
        # File for user-defined bad pixels.
        self._userfile = None
                
        # Number of pixels defined bad by algorithm and by user,
        # used for output file header metadata.
        self._nbad_auto = 0
        self._nbad_user = 0
                
        # Number of MAD std deviations for sigma clipping
        self._sigma = sigma

        # Process data.
        self._imdata, self._imhdr = self._read_fits(darkfile, 0)
        self._image_stats()
        self._generate_sigmaclip_mask(self._imdata, self._sigma)
        return
        
    def _add_bad_columns(self, bad_col_list):
        """Add bad columns to the _badpixmask, setting the mask values
           to ApFindBadPixels.USER_BAD
        """
        num_cols     = len(bad_col_list)
        num_user_bad = 0
        ncols        = self._badpixmask.shape[1]
        nrows        = self._badpixmask.shape[0]
        self._logger.info(f'Adding {num_cols} bad columns to {nrows} row, {ncols} column mask.')
        
        for col in bad_col_list:
            # Convert to python zero based index
            col1 = col - 1
            col2 = col
            if (col1 < 0) or (col1 >= ncols):
                # Skip as outside image range.
                msg = f'Warning, column {col} (1-based) outside image.'
                self._logger.warning(msg)
            else:
                self._logger.debug(f'Setting [:,{col1}:{col2}] to user-defined bad value.')
                self._badpixmask[:,col1:col2] += ApFindBadPixels.USER_BAD
                num_user_bad += nrows

        self._logger.debug(f'Set {num_user_bad} pixels to user-defined bad value.')
        return num_user_bad
        
    def _add_bad_rectangles(self, bad_rectangle_list):
        """Add bad rectangular regions to the _badpixmask, setting the mask values
           to ApFindBadPixels.USER_BAD
        """
        num_rect     = len(bad_rectangle_list)
        ncols        = self._badpixmask.shape[1]
        nrows        = self._badpixmask.shape[0]
        num_user_bad = 0
        self._logger.info(f'Adding {num_rect} bad rows to {nrows} row, {ncols} column mask.')
        
        for rect in bad_rectangle_list:
            if len(rect) != 4:
                msg = f'Error, expecting 4-element list, got {rect}. Skipping.'
                self._logger.warning(msg)
                continue
            
            # Convert to python zero based index. Here we convert from
            # inclusive 1-based FITS style to [) zero-based indices.
            row1 = rect[0] - 1
            row2 = rect[1]
            col1 = rect[2] - 1
            col2 = rect[3] 
            if (row1 < 0) or (row2 > nrows):
                # Skip as outside image range.
                msg = f'Warning, row range {row1}:{row2} (0-based) outside image.'
                self._logger.warning(msg)
            elif (col1 < 0) or (col2 > ncols):
                # Skip as outside image range.
                msg = f'Warning, column range {col1}:{col2} (0-based) outside image.'
                self._logger.warning(msg)
            else:
                self._logger.debug(f'Setting [{row1}:{row2},{col1}:{col2}] bad.')
                self._badpixmask[row1:row2,col1:col2] += ApFindBadPixels.USER_BAD
                num_user_bad += (row2-row1)*(col2-col1)
                
        self._logger.debug(f'Set {num_user_bad} pixels to user-defined bad value.')
        return num_user_bad
    
    def _add_bad_rows(self, bad_row_list):
        """Add bad rows to the _badpixmask, setting the mask values
           to ApFindBadPixels.USER_BAD
        """
        num_rows     = len(bad_row_list)
        ncols        = self._badpixmask.shape[1]
        nrows        = self._badpixmask.shape[0]
        num_user_bad = 0
        self._logger.info(f'Adding {num_rows} bad rows to {nrows} row, {ncols} column mask.')
        
        for row in bad_row_list:
            # Convert to python zero based index
            row1 = row - 1
            row2 = row
            if (row1 < 0) or (row1 >= nrows):
                # Skip as outside image range.
                msg = f'Warning, row {row} (1-based) outside image.'
                self._logger.warning(msg)
            else:
                self._logger.debug(f'Setting [{row1}:{row2},:] bad.')
                self._badpixmask[row1:row2,:] += ApFindBadPixels.USER_BAD
                num_user_bad += ncols
                
        self._logger.debug(f'Set {num_user_bad} pixels to user-defined bad value.')
        return num_user_bad
        
    def _check_file_exists(self, filename):
        """Checks the file exists and cleans up the path
        """
        
        fpath = Path(filename).expanduser() 
        if not fpath.exists():
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return fpath

    def _generate_sigmaclip_mask(self, data, sigma):
        """Creates a bad pixel mask based on sigma-clipped statistics
           of the input data array.
           
        This routine is most appropriate for images that expected to
        be relatively uniform, but with a small number of highly 
        discrepant values.
        
        The generated mask is of type uint8, with value 1 where pixels 
        are BAD and zero where pixels are GOOD. 
                   
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
        
        # Generate final mask in uint8 form.
        self._badpixmask = np.logical_or(bad_lo, bad_hi).astype('uint8')
        nbad             = np.sum(self._badpixmask)
        pct_bad          = 100 * (nbad/npix)
        msg              = f'Out of {npix} pixels, {nbad} are bad ({pct_bad:.4f}%).'
        self._logger.info(msg)
        
        # Update count of algorithmically-defined bad pixels.
        self._nbad_auto = nbad
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
        
        image_filename = self._check_file_exists(image_filename)
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
        
    def _read_user_badpix(self, user_badpix_file):
        """Read user-defined bad columns, row and rectangles from a YaML
           file.
        
        Returns a set of lists of the bad columns, bad rows, 
        and bad rectangles (the latter consisting of lists) if these
        are defined in the yaml file, or Nones otherwise.
        """
        
        badcols = None
        badrows = None
        badrect = None
        
        user_badpix_file = self._check_file_exists(user_badpix_file)
        with open(user_badpix_file) as bpfile:
            lines = bpfile.read()
            yobj  = yaml.safe_load(lines)
            
            if 'bad_columns' in yobj:
                badcols = yobj['bad_columns']
                self._logger.debug(f'There are {len(badcols)} bad columns in the user-defined badpixel file.')
            else:
                self._logger.debug('There were no bad columns in the user-defined badpixel file.')
                
            if 'bad_rows' in yobj:
                badrows = yobj['bad_rows']
                self._logger.debug(f'There are {len(badrows)} bad rows in the user-defined badpixel file.')
            else:
                self._logger.debug('There were no bad rows in the user-defined badpixel file.')
                
            if 'bad_rectangles' in yobj:
                badrect = yobj['bad_rectangles']
                self._logger.debug(f'There are {len(badrect)} bad rectangles in the user-defined badpixel file.')
            else:
                self._logger.debug('There were no bad rectangles in the user-defined badpixel file.')

        # We don't want to deal with zero length objects.
        if len(badcols) == 0:
            badcols = None
        if len(badrows) == 0:
            badrows = None
        if len(badrect) == 0:
            badrect = None
        
        return badcols, badrows, badrect
        
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
        if self._userfile is not None:
            kw_dict['USERFILE'] = (self._userfile.name, 'User-defined bad pixel file.')
        kw_dict['NBADAUTO'] = (self._nbad_auto, 'Number of algorithm-detected bad pixels.')
        kw_dict['NBADUSER'] = (self._nbad_user, 'Number of user-defined bad pixels.')
    
        for kw in copy_list:
            if kw in self._imhdr:
                val = self._imhdr[kw]
                cmt = self._imhdr.comments[kw]
                kw_dict[kw] = (val, cmt)
            
        for kw in kw_dict:
            hdu.header[kw] = kw_dict[kw]
        hdu.header['HISTORY'] = f'Processed by {self._name} {__version__} at {tnow}'
        return
    
    def add_user_badpix(self, user_badpix_file):
        """
        
        :param user_badpix_file: Path/name of a YaML containing user
          defined bad columns, bad rows, and/or bad rectangular regions.
        """
        
        user_badpix_file = Path(user_badpix_file).expanduser()
        
        self._logger.info(f'Processing user-defined bad pixels from {user_badpix_file}')
        badcols, badrows, badrect = self._read_user_badpix(user_badpix_file)
        self._userfile = user_badpix_file
        
        num_user_bad = 0
        if badcols is not None:
            num_user_bad += self._add_bad_columns(badcols)
        if badrows is not None:
            num_user_bad += self._add_bad_rows(badrows)
        if badrect is not None:
            num_user_bad += self._add_bad_rectangles(badrect)
            
        # Update count of user-defined bad pixels.
        self._nbad_user = num_user_bad
        self._logger.debug(f'Total number of user-defined bad pixels applied to mask: {num_user_bad}')
        return
    
    def get_mask(self):
        """Returns the bad pixel mask as a uint8 numpy array.
        """
        return self._badpixmask
    
    def write_mask(self, mask_file_name):
        """Write the bad pixel mask to a FITS file with the user 
           specified name/path.
           
        The output FITS data array contains pixels that can have the
        *sum* of the following numeric values. A given pixel may be
        flagged bad based on both the statistical deviation in the dark
        and the user, so such a pixel would have a value of 3.
        0: Good pixels.
        1: Algorithmically detected bad pixels.
        2: User-defined bad pixels.
        4: Reserved for future use.
        8: Reserved for future use.
        16: Reserved for future use.
        32: Reserved for future use.
        64: Reserved for future use.
        128: Reserved for future use.
           
        :param mask_file_name: File name/path for the output bad pixel
          mask. The file will be overwritten if it exists.
        """
        
        hdu  = fits.PrimaryHDU(data=self._badpixmask)
        self._update_header(hdu)
        
        hdu_list = fits.HDUList([hdu])
        hdu_list.writeto(mask_file_name, overwrite=True)
        self._logger.info(f'Wrote bad pixel mask to {mask_file_name}')
        return
    
