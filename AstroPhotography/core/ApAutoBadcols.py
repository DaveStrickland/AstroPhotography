"""Contains the implementation of the ApAutoBadcols class.
"""

# 2021-08-15 dks : Initial implementation.
# 2024-04-11 dks : Issue-023 improve functionality

import sys
import logging
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

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
        """
        Initializes an ApAutoBadcols instance.
        """
        
        self._name     = 'ApAutoBadcols'
        self._version  = __version__
        self._loglevel = loglevel
        self._initialize_logger(self._loglevel)
        
        # Information about file or data processed, to be carried
        # forward to any output files.
        self._meta = None
        
        # Computed statistics
        self._badcols  = None
        self._badrows  = None
        self._colstats = None
        self._rowstats = None
        return
        
    def _check_file_exists(self, filename, throws=True):
        """
        Checks if the named file exists, returning True if it does.
        
        By default it raises an exception if the named file does not exist,
        but this can be disabled by setting throws=False
        """
        exists = Path(filename).exists()
        if (not exists) and throws:
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return exists

    def _get_metadata_from_hdr(self, filename, fitshdr):
        """
        Extract some file metadata from the file name and a FITS header
        """
        
        fpath = Path(filename)
        metadict = {'filename': fpath.name}
        for kw in ['BITPIX', 'NAXIS1', 'NAXIS2']:
            if kw in fitshdr:
                metadict[kw] = fitshdr[kw]
        return metadict

    def _initialize_logger(self, loglevel):
        """
        nitialize and return the logger
        """
        
        self._logger = logging.getLogger(self._name)
        
        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: {}'.format(loglevel))
        self._logger.setLevel(numeric_level)
    
        # check if handlers already present
        if not len(self._logger.handlers):
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
        
        self._meta = self._get_metadata_from_hdr(image_filename, ext_hdr)
        self._logger.debug(f'File metadata: {self._meta}')
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
            
            # In some cases where the data in the local mask is very similar
            # the sigma clipped standard deviation is returned as zero, 
            # which causes problems. In such case, return the normal
            # standard deviation within the masked area
            if ( cstd == 0 ):
                cstd = np.nanstd(local_data)
            
            mean_data[idx] = cmean
            std_data[idx]  = cstd
        return mean_data, std_data
    
    def process_fits(self, fitsimg, nsigma=5, window_len=11):
        """
        Identify bad columns and rows in a numpy 2-dimensional array,
        returning a 1-d array of the (zero-based) bad column and row
        indices, or None if there are none.
        """
        ext_num = 0
        idata, ihdr = self._read_fits(fitsimg, ext_num)
        badcols, badrows = self.process(idata, nsigma, window_len)
        return badcols, badrows
    
    def process(self, data_array, nsigma=5, window_len=11):
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
        
        # Store for possible later use
        if badcols is not None:
            self._badcols = badcols.copy()
        
        if badrows is not None:
            self._badrows = badrows.copy()
        
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
        
        # Create recarray, data not initialized
        values_arr = np.recarray((nvals,),
            dtype=[('idx', int),
                ('median', float),
                ('local_mean', float),
                ('local_std', float),
                ('nsigma', float),
                ('isbad', int)])
        idx_arr = np.arange(nvals, dtype=int)
        values_arr['idx'] = idx_arr
        values_arr['median'] = median_array
        values_arr['local_mean'] = sldng_mean
        values_arr['local_std'] = sldng_std
        values_arr['nsigma'] = nsigma_from_mean
        values_arr['isbad'] = bad_mask
        if 'column' in type_str:
            self._colstats = values_arr.copy()
        elif 'row' in type_str:
            self._rowstats = values_arr.copy()
        else:
            raise RuntimeError(f'Error, unexpected type_str ({type_str}) in _process')
        
        # Create information for column by column (or row by row) debug level output.
        # This does the first nalways columns/rows irrespective of whether
        # they are good or bad, and then only the bad ones.
        nalways      = 40
        dbg_str_list = []
        hdr_str      = '{:>4s}, {:>10s}, {:>10s}, {:>10s}, {:>10s}, {:>6s}'.format(short_str, 
            'median', 'local_mean', 'local_std', 'nsigma', 'isbad?')
        dbg_str_list.append(f'Diagnostics for first {nalways} {type_str}s and all bad {type_str}:')
        dbg_str_list.append( hdr_str )
        
        # Use primitave arrays instead of values_arr because of nested ''
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
        
        # store nsigma and window length used
        if self._meta is None:
            self._meta = {}
        self._meta['nsigma'] = nsigma
        self._meta['window_len'] = window_len
        return bad_indices
        
    def generate_stats_plot(self, plotstatfile):
        """
        Generate a two-panel plot of the column and row statistics
        """
        
        fig, ax_arr = plt.subplots(nrows=2, ncols=1)
        
        # Generate title string
        if 'filename' in self._meta:
            fname = self._meta['filename']
        else:
            fname = 'unknown file'
        nsigma = self._meta['nsigma']
        window_len = self._meta['window_len']
        title_str = f'Input file {fname}\nSigma threshold={nsigma:.2f}, sliding window_len={window_len}' 
        fig.suptitle(title_str, fontsize=7)
        
        # Plot bad columns
        idx = 0
        self._generate_stats_panel(ax_arr[idx], self._colstats,
            'Along-column statistics and bad columns',
            'column')
        idx = 1
        self._generate_stats_panel(ax_arr[idx], self._rowstats,
            'Along-row statistics and bad rows',
            'row')
                    
        fig.tight_layout()
        plt.savefig(plotstatfile,
            dpi=200,
            bbox_inches='tight')
        self._logger.info(f'Column/row statistics info plot written to {plotstatfile}')        
        return
        
    def _generate_stats_panel(self, ax, stats_arr, titlestr, coltype):
        """
        Generate one of the panels for the column/row statistics plot
        """

        if 'column' in coltype:
            xlabel_str = 'X-axis column index (pixels)'
        else:
            xlabel_str = 'Y-axis row index (pixels)'
        ylabel_str = f'Median along {coltype}'
        
        # Number of sigma around local sliding mean to plot envelope
        env_sigma = 3.0
        medn_sub_mean = stats_arr['median'] - stats_arr['local_mean']
        envelope_low = -env_sigma*stats_arr['local_std']
        envelope_hi  = env_sigma*stats_arr['local_std']
        
        mask           = stats_arr['isbad'] == 1
        numbad         = np.sum( stats_arr['isbad'] )
        masked_bad_idx = stats_arr['idx'][mask]
        masked_bad_median = stats_arr['median'][mask]
        masked_bad_med_sub_mean = medn_sub_mean[mask]
        
        # min/max median values
        minval  = 0.995 * np.nanmin( stats_arr['median'] )
        maxval  = 1.005 * np.nanmax( stats_arr['median'] )
        maxval2 = 1.005 * self._meta['nsigma'] * np.nanmax( stats_arr['local_std'] )
        minval2 = -maxval2

        ax.set_title(titlestr, fontsize=7, pad=2)
        ax.tick_params('both', labelsize=4, labelcolor='blue')
        ax.step(stats_arr['idx'], stats_arr['median'], color='blue', where='mid',
            label=ylabel_str, linewidth=0.5, alpha=0.7)
        ax.set_ylim( (minval, maxval) )
        ax.minorticks_on()
        ax.scatter(masked_bad_idx, masked_bad_median, s=1.5, c='red', marker='o', label=f'Bad {coltype}s ({numbad})')
        ax.set_xlabel(xlabel_str, fontsize=6)
        ax.set_ylabel(ylabel_str, fontsize=6, color='blue')

        plt.rc('legend', fontsize = 6)
        lgnd1 = ax.legend(loc='lower left')

        ax2 = ax.twinx()
        if 'column' in coltype:
            ylabel_str2 = 'Column median minus local mean' 
        else:
            ylabel_str2 = 'Row median minus local mean' 
        ax2.set_ylabel(ylabel_str2, fontsize=6, color='green')

        ##ax2.fill_between(stats_arr['idx'], envelope_low, envelope_hi, color='cyan', alpha=0.2)
        ax2.step(stats_arr['idx'], medn_sub_mean, where='mid', color='green', alpha=0.7, linewidth=0.75)
        ax2.scatter(masked_bad_idx, masked_bad_med_sub_mean, s=1.5, c='red', marker='o', label=f'Bad {coltype}s ({numbad})')

        ax2.set_ylim( (minval2, maxval2) )
        ax2.tick_params('both', labelsize=4, labelcolor='green')

        return
        
    def write_badcols_file(self, badcolfile, overwrite=False):
        """
        Write the bad columns and bad rows in a YaML-like format
        """
        
        exists = self._check_file_exists(badcolfile, False);
        if exists and not overwrite:
            self._logger.error('Bad column/row YaML file exists and overwrite=False')
            self._logger.error('  No output file will be written unless the file is first deleted or overwrite is specified.')
            return
            
        # TODO use an actual YaML writer?
        with open(badcolfile, 'w', encoding="utf-8") as f:
            # generate info string
            f.write( '---\n' )
            f.write( f'# Column statistics file generated by {__name__} version {__version__}\n' )
            if 'filename' in self._meta:
                fname = self._meta['filename']
            else:
                fname = 'unknown file'
            f.write( f'# Generated from {fname}\n' )
            nsigma = self._meta['nsigma']
            window_len = self._meta['window_len']
            f.write( f'# Processing parameters: badness sigma threshold={nsigma:.2f}, sliding window_len={window_len}\n' )

            if self._badcols is not None:
                if len(self._badcols) > 0:
                    # Add one to get FITS-like indexing
                    badcols = self._badcols + 1

                    f.write('bad_columns:\n')
                    for val in badcols:
                        f.write(f'- {val:d}\n')
                else:
                    f.write('bad_columns: {}\n') # show it is empty
            else:
                f.write('# No bad columns detected.\n')
            if self._badrows is not None:
                if len(self._badrows) > 0:
                    # Add one to get FITS-like indexing
                    badrows = self._badrows + 1

                    f.write('bad_rows:\n')
                    for val in badrows:
                        f.write(f'- {val:d}\n')
                else:
                    f.write('bad_rows: {}\n')   # show it is empty
            else:
                f.write('# No bad rows detected.\n')
                    
            f.write( '...\n' )
            self._logger.info(f'Wrote bad column/row YaML file to {badcolfile}' )
        return
        
    def write_stats(self, fcolstat, frowstat):
        """
        Write CSV files of the computed along-column and along-row statistics.
        
        :param fcolstat: Name for the output CSV of along-column statistics.
          Note that this will be overwritten if it already exists.
        :param frowstat: Name for the output CSV of along-row statistics.
          Note that this will be overwritten if it already exists.
        """
        
        # generate info string
        info_strings = [ f'# Column statistics file generated by {__name__} version {__version__}\n']
        if 'filename' in self._meta:
            fname = self._meta['filename']
        else:
            fname = 'unknown file'
        info_strings.append( f'# Generated from {fname}\n' )
        nsigma = self._meta['nsigma']
        window_len = self._meta['window_len']
        info_strings.append( f'# Processing parameters: badness sigma threshold={nsigma:.2f}, sliding window_len={window_len}\n')
        
        chdr_str = '{:3s},{:>10s},{:>10s},{:>10s},{:>10s},{:>5s}'.format('col', 'median', 'local_mean', 'local_std', 'nsigma', 'isbad')
        rhdr_str = chdr_str.replace('col', 'row')
        
        if fcolstat is not None:
            with open(fcolstat, 'w', encoding="utf-8") as f:
                f.writelines(info_strings)
                np.savetxt(f, self._colstats, header=chdr_str,
                    fmt=['%05d', '%10.2f', '%10.2f', '%10.2f', '%10.2f', '%5d'], delimiter=',')
                self._logger.debug(f'Wrote column statistics CSV data to {fcolstat}')

        if frowstat is not None:
            with open(frowstat, 'w', encoding="utf-8") as f:
                f.writelines(info_strings)
                np.savetxt(f, self._rowstats, header=rhdr_str,
                    fmt=['%05d', '%10.2f', '%10.2f', '%10.2f', '%10.2f', '%5d'], delimiter=',')
                self._logger.debug(f'Wrote row statistics CSV data to {frowstat}')
        
        return
