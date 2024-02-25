"""Contains the implementation of the ApAutoBadcols class.
"""

# 2021-08-15 dks : Initial implementation.

import sys
import logging
from pathlib import Path
import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

# AstroPhotography includes    
from .. import __version__

class ApAutoBadcols:
    """
    Automatically detect the worst bad columns and bad rows in an 
    astronomical image.
    """
    
    def __init__(self,
        loglevel):
        """Initializes an ApAutoBadcols instance.
        """
        
        self._name     = 'ApAutoBadcols'
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
    
    def _sliding_stats_1d(self, idata, window_len):
        """
        Simple brute force 1-dimensional sliding window statistics
        
        :param data: 1-dimensional array
        :param window_len: Integer window length, must be positive odd number
          e.g. 11.
        """
        hw        = int((window_len - 1)/2)
        nvals     = idata.size
        mean_data = np.zeros(nvals)
        std_data  = np.zeros(nvals)
        
        for idx in range(nvals):
            min_idx    = max(0, idx - hw)          # inclusive
            max_idx    = min(nvals, idx + hw +1)   # exclusive
            local_data = idata[min_idx:max_idx]
            
            # Use sigma clipping, because even a single discrepant
            # value can through a normal standard deviation off.
            cmean, cmedian, cstd = sigma_clipped_stats(local_data)
            
            mean_data[idx] = cmean
            std_data[idx]  = cstd
        return mean_data, std_data
    
    def process_fits(self, fitsimg, nsigma=None, window_len=None):
        """
        Identify bad columns and rows in a numpy 2-dimensional array,
        returning a 1-d array of the (zero-based) bad column and row
        indices, or None if there are none.
        """
        ext_num = 0
        idata, ihdr = self._read_fits(fitsimg, ext_num)
        badcols, badrows = self.process(idata, nsigma, window_len)
        return badcols, badrows
    
    def process(self, data_array, nsigma=None, window_len=None):
        """
        Identify bad columns and rows in a numpy 2-dimensional array,
        returning a 1-d array of the (zero-based) bad column and row 
        indices, or None if there are none..
        """
        
        if nsigma is None:
            nsigma = 5.0        # Want to be sure these are really clearly bad.
        if window_len is None:
            window_len = 11
                    
        nrows = data_array.shape[0]
        ncols = data_array.shape[1]

        # Look for bad columns, median over row axis (collapse down)
        medn_over_cols = np.nanmedian(data_array, axis=0)
        badcols        = self._process(medn_over_cols, 0, nsigma, window_len)

        # Look for bad rows, median over column axis (collapse across)
        medn_over_rows = np.nanmedian(data_array, axis=1)
        badrows        = self._process(medn_over_rows, 1, nsigma, window_len)
        
        return badcols, badrows

    def _process(self, median_array, axis_used, nsigma, window_len):
        """
        Internal utility that performs the work of processing the 
        1-dimensional column or row inspection
        
        :param median_array: A 1-d numpy array of the medians. If the
          median over all rows (axis=0) is supplied for each column,
          then this will be used to find bad columns.
        :param axisname: The axis used when generating the median array.
          This is an integer that is either 0 or 1. If 0, then the median
          over all rows is performed per column, so this is used to find
          bad columns. If 1, the median over all columns was calculated,
          so this is used to identify bad rows.
        :param nsigma: Values greater than or equal to this number of 
          standard deviations away from local mean are considered bad.
        :param window_len: Total width of the sliding window used to
          assess the local mean value. This should be an odd number.
        """
        type_str  = 'column'
        short_str = 'col'
        if axis_used == 1:
            type_str  = 'row'
            short_str = 'row'
        elif axis_used > 1:
            raise ValueError(f'axis_used should be 0 or 1, not {axis_used}')
        
        nvals                 = median_array.size
        sldng_mean, sldng_std = self._sliding_stats_1d(median_array, window_len)
        nsigma_from_mean      = np.abs(median_array - sldng_mean) / sldng_std
        bad_mask              = nsigma_from_mean >= nsigma
        self._logger.info(f'Found {np.sum(bad_mask)} bad {type_str}s out of {nvals} {type_str}s.')
        
        # Create information for column by column (or row by row) debug level output.
        # This does the first nalways columns/rows irrespective of whether
        # they are good or bad, and then only the bad ones.
        nalways      = 40
        dbg_str_list = []
        hdr_str      = '{:>4s}, {:>10s}, {:>10s}, {:>10s}, {:>10s}, {:>6s}'.format(short_str, 
            'median', 'local_mean', 'local_std', 'nsigma', 'isbad?')
        dbg_str_list.append(f'Diagnostics for first {nalways} {type_str}s and all bad {type_str}:')
        dbg_str_list.append( hdr_str )
        for idx in range(nvals):
            if (idx < nalways) or (bad_mask[idx]):
                dbg_str = f'{idx:04d}, {median_array[idx]:10.2f}, {sldng_mean[idx]:10.2f}, {sldng_std[idx]:10.2f}, {nsigma_from_mean[idx]:10.2f}, {bad_mask[idx]}'
                dbg_str_list.append(dbg_str)
        self._logger.debug('\n'.join(dbg_str_list))
        
        # Convert mask to indices
        if np.sum(bad_mask) > 0:
            bad_indices = np.arange(nvals)[bad_mask]
        else:
            bad_indices = None
        
        return bad_indices
        
