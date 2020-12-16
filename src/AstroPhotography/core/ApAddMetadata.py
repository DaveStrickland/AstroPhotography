"""Contains the implementation of the ApAddMetadata class.
"""

# 2020-12-16 dks : Initial implementation.

import sys
import logging
from pathlib import Path
import math
import time
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits

from .. import __version__

class ApAddMetadata:
    """Adds metadata to the FITS header of a calibrated fits file.
    
    Currently this is written to use the information in iTelescope
    file names, plus information derived from that using astropy.
    In the longer term other ways of inputting the data could be added.
    """
    
    def __init__(self,
        loglevel):
        """Initializes an ApAddMetadata instance.
        """
        
        self._name     = 'ApAddMetadata'
        self._version  = __version__
        self._loglevel           = loglevel
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
        """Read a single extension's data and header from a FITS file.
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
            
        return ext_data, ext_hdr

    def _write_corrected_image(self, fitsfile, kwdict):
        """Updates the header keywords of the specified file with 
           the computed values.
        """
        
        self._logger.debug(f'FITS header keywords to be added to output: {kwdict}')
        self._check_file_exists(inpdata_file)
            
        # open() parameters that can be important.
        # Default values used here.
        # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
        uint_handling = True
        image_scaling = False
        ext_num       = 0
            
        # Re-open the file in update mode.
        with fits.open(fitsfile, mode='update', 
            uint=uint_handling, 
            do_not_scale_image_data=image_scaling) as hdu_list:
                        
            # Modify header
            for kw, val in kwdict.items():
                hdu_list[ext_num].header[kw] = val
            
            tnow = datetime.now().isoformat(timespec='milliseconds')
            hdu_list[ext_num].header['HISTORY'] = f'Modified by {self._name} {self._version} at {tnow}'

            # Flush the changes back to disk.
            hdu_list.flush()
                    
        self._logger.info(f'Finished updating FITS header of {fitsfile}')
        return
    
    def process(self, fitsfile, mode):
        """Adds or updates FITS header keywords in the name FITS file.
        
        This function generates or populates the values of FITS header
        keywords that are required by ApFindStars, ApAstrometry and
        ApQualitySummarizer.
        
        This function supports the following modes:
        - 'iTelescope': The telescope, target, and observer are parsed 
          from the file name. In combination with the date/time of the
          observation this allows the target RA, Dec, site 
          lat/lon/elevation, and airmass to be determined.
        
        :param fitsfile: Path/name of the FITS file to be modified.
        :param mode: Mode. Must be one of the modes specified above.
        """
        
        if 'iTelescope' in mode:
            # TODO
            print('TODO')
        else:
            err_msg = f'Error, unexpected/unsupported mode {mode}.'
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)
        
        return
