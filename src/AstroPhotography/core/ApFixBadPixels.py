"""Contains the implementation of the ApFixBadPixels class.
"""

#  2020-11-16 dks : Working version of ApFixBadPixels.

import sys
import logging
from pathlib import Path
import math
import time
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from .. import __version__

class ApFixBadPixels:
    """A class used to fix pre-indentified bad pixels within an image
       by replacing them with the median value for surrounding good 
       pixels.
    """
    
    # Class level constants
    MASK_GOOD = 0
    
    def __init__(self,
        loglevel):
        """Constructs an ApFixBadPixels object. No processing is performed.
        
        :param loglevel: Logging level to use.
        """
    
        self._name = 'ApFixBadPixels'
    
        # Initialize logging
        self._loglevel = loglevel
        self._initialize_logger(self._loglevel)
        
        # Minimum number of surrounding "good" pixels required to allow
        # patching of bad pixels by the median of the good pixels.
        # The smallest possible value is 1, but whether that gives 
        # scientifically valid output is not clear.
        self._min_valid = 4
        
        # Whether to replace unfixable pixel values with a artificial
        # fill value, or simply leave them unchanged (recommended)
        self._replace_unfixable       = False
        self._replace_unfixable_value = np.nan
        
        self._logger.debug(f'{self._name} instance constructed.')
        return
        
    def _check_file_exists(self, filename):
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
        """Writes the bad-pixel-corrected data to the specified output
           file, preserving all other items from the original input file.
           
        The output file differs from the original input data file in 
        having the bad-pixel corrected image, and additional header
        keywords:
        - BPIXCORR: Logical true denoting whether bad pixel correction
                    applied.
        - BPIXFILE: The name of the bad pixel file used, stripped of any
                    preceding path elements.
        - BPIX_MIN: Minimum number of surrounding good pixels required 
                    for a correction to be attempted. 
        - BPIXDPIX: Delta pixels, the distance around the bad pixel that
                    surrounding pixels are drawn from. 
        - BPIXNBAD: Number of bad pixels in the bad pixel file.
        - BPIXNFIX: Number of pixels successfully corrected.
        - BPIXNREM: Number of remaining bad pixels, the number that could
                    not be corrected using BPIXDPIX and BPIX_MIN
        
        :param inpdata_file: Input FITS data file affected by bad pixels.
        :param ext_num: Extension number for data array and header.
        :param outdata_file: Modified copy of the input data file where the
          bad pixels have had their data valus modified by the median
          of the surrounding good pixels.
        :param odata: Bad-pixel corrected data array.
        :param odict: Input dictionary generated by fix_bad_pixels (and
          fix_files), used to update output header as described above.
          Elements of the dictionary that have keywords containing
          'BPIX' are written to the output fits file, any other
          elements are not.
        """
        
        self._logger.debug(f'Bad pixel keywords added to output: {odict}')

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
            
            # Modify data
            hdu_list[ext_num].data = odata
            
            # Modify header
            for kw, val in odict.items():
                if 'BPIX' in kw:
                    hdu_list[ext_num].header[kw] = val
            
            tnow = datetime.now().isoformat(timespec='milliseconds')
            hdu_list[ext_num].header['HISTORY'] = f'Applied {self._name} {__version__} at {tnow}'
            
            # Write to new file
            hdu_list.writeto(outdata_file, 
                output_verify='ignore',
                overwrite=True)
        
        self._logger.info(f'Wrote bad pixel corrected file to {outdata_file}')
        return

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

        ext_num  = 0 # Really should have a better way of setting this.
        deltapix = int(deltapix)
        
        msg = (f'fix_files input data file={inpdata_file},'
            f' mask file={badpixmask_file},'
            f' output file={outdata_file}, deltapix={deltapix}')
        self._logger.info(msg)
        
        idata, ihdr     = self._read_fits(inpdata_file, ext_num)
        mskdata, mskhdr = self._read_fits(badpixmask_file, ext_num)
        
        odata, odict = self.fix_bad_pixels(idata, mskdata, deltapix)
        odict['BPIXFILE'] = (Path(badpixmask_file).name,
            'Name of master bad pixel file used')
        
        # Create copy of input image and write modified data to it
        # with an updated header
        self._write_corrected_image(inpdata_file, 
            ext_num,
            outdata_file,
            odata, 
            odict)
        return

    def fix_bad_pixels(self, data, badpixmask, deltapix=1):
        """Fix bad pixels in the input array based on the mask in an
           badpixel array, returning the modified data array, 
           and a dictionary summarizing the number of bad pixels, pixels 
           corrected, and pixels than could not be corrected.
           
        The output dictionary can be used to modify the FITS header of
        an output file, and consists of keyword: (value, comment) pairs.
        
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

        deltapix = int(deltapix)
        
        # Info message about array shape and data type.
        self._logger.info(f'fix_bad_pixels: data has {data.shape[0]} rows x {data.shape[1]} columns, dtype={data.dtype}')
        self._logger.info(f'fix_bad_pixels: mask has {badpixmask.shape[0]} rows x {badpixmask.shape[1]} columns, dtype={badpixmask.dtype}')
        
        # Check that sizes match
        if (data.shape != badpixmask.shape):
            msg = (f'Error, the shape of the input data array ({data.shape})'
                f' does not match that of the bad pixel mask array ({badpixmask.shape}).')
            self._logger.error(msg)
            raise RunTimeError(msg)
        
        if not np.issubdtype(data.dtype, np.floating):
            # Not a floating point datatype.
            msg = ( 'Pixel medians may suffer from casting truncation because'
                f' the input data is not a floating point datatype ({data.dtype}).')
            self._logger.warning(msg)

        
        newdata = data.copy()
        mask    = badpixmask != ApFixBadPixels.MASK_GOOD
        npix    = data.size
        nbad    = np.sum(mask)
        pctbad  = 100.0 * nbad / npix
        self._logger.debug(f'Percentage of pixels considered bad: {pctbad:.3f} ({nbad:d}/{npix:d})')
        fixed_stats = {'numpix': (npix, 'Total number of pixels in image'),
            'BPIXNBAD': (nbad,     'Total number of bad pixels in bad pixel file'),
            'pctbad':   (pctbad,   'Percentage of pixel defined bad'),
            'BPIX_MIN': (self._min_valid, 'Minimum number of good neighors needed'),
            'BPIXDPIX': (deltapix, 'Half height/width of collection region (pixels)')}
            
        # Mask to record which pixels were not fixed.
        newmask = mask.copy()

        # Formatting for diagnostic output
        nprint  = 40
        idx_fmt = '3d'
        pix_fmt = '4d'
        msg     = f'Diagnostics of the first {nprint} bad pixels follow:\n'
        hdr     =  '{}  {}  {}  {}  {}  {}  {}  {:>5s}  {:>10s}  \n'.format('idx',
            'row_',  'col_',  'rmin:rmax',  'cmin:cmax',
            'good',  'bad_', 'fix?',  'new_value')
        msg    += hdr
        line    = ('{:>' f'{idx_fmt}' '}  '                  # idx_
            '{:>' f'{pix_fmt}'  '}  '                        # row_
            '{:>' f'{pix_fmt}' '}  '                         # col_
            '{:>' f'{pix_fmt}' '}:{:<' f'{pix_fmt}' '}  '    # rmin:rmax
            '{:>' f'{pix_fmt}' '}:{:<' f'{pix_fmt}' '}  '    # cmin:cmax
            '{:>' f'{pix_fmt}' '}  '                         # good
            '{:>' f'{pix_fmt}' '}  '                         # bad_
            '{:>5s}'                                         # fix?
            '{:>10s}'                                        # new_value
            '\n'
            )
        ##print(line)

        # Iterate over pixels.
        # - Construct row and column index arrays, then select only the
        #   indices of the bad pixels and store those as 1-dimensional arrays.
        nrows = mask.shape[0]
        ncols = mask.shape[1]
        perf_time_start    = time.perf_counter() # highest res timer, counts sleeps
        row_idxs, col_idxs = np.mgrid[0:nrows, 0:ncols]
        bp_row_idxs        = row_idxs[mask].ravel()
        bp_col_idxs        = col_idxs[mask].ravel()
        for idx, pos in enumerate(zip(bp_row_idxs, bp_col_idxs)):
            ridx = pos[0]
            cidx = pos[1]
            rmin = max(0,     ridx - deltapix)     # included
            rmax = min(nrows, ridx + deltapix + 1) # excluded
            cmin = max(0,     cidx - deltapix)     # included
            cmax = min(ncols, cidx + deltapix + 1) # excluded
            
            # Get data and mask cut outs.
            # - Using data and mask, not newdata and newmask, means that
            #   previous pixel corrections aren't used when updating pixels.
            data_co  = data[rmin:rmax, cmin:cmax]
            mask_co  = mask[rmin:rmax, cmin:cmax]  # True where bad
            nbad_co  = np.sum(mask_co)
            ngood_co = mask_co.size - nbad_co
                        
            can_fix  = False
            if ngood_co >= self._min_valid:
                # Flip mask so good pixels are True and we can select them.
                good_mask = np.logical_not(mask_co)
                
                can_fix             = True
                newmask[ridx, cidx] = ApFixBadPixels.MASK_GOOD
                ##if idx < 3:
                ##    print(data_co)
                ##    print(good_mask)
                ##    print(data_co[good_mask])
                medval              = np.median( data_co[good_mask] )
                fix_str             = f'{medval:>10.2f}  '
                newdata[ridx, cidx] = medval
            else:
                fix_str = '{:>10s}  '.format('N/A')
            
            # Diagnostic string
            if idx < nprint:
                msg += line.format(idx, ridx, cidx, 
                    rmin, rmax, cmin, cmax,
                    ngood_co, nbad_co,
                    str(can_fix), fix_str)
                    
        perf_time_end   = time.perf_counter() # highest res timer, counts sleeps
        self._logger.debug(msg) 

        # Measure performance, in case we need to speed this up later.
        run_time   = perf_time_end - perf_time_start
        ms_per_pix = 1000 * run_time / nbad 
        msg = f'Processed {nbad} pixels in {run_time:.3f} s, {ms_per_pix:.2f} ms per bad pixel.'
        self._logger.info(msg)        

        # Work out how many pixels could not be fixed
        nnotfix = np.sum(newmask)
        fixed_stats['BPIXNREM'] = (nnotfix, 'Number of bad pixels not corrected')
        if nnotfix > 0:
            msg = (f'Could not fix {nnotfix} pixels as they had less'
                f' than {self._min_valid} good neighbors when deltapix={deltapix} pixels.')
            self._logger.warning(msg)
        nfixed = nbad - nnotfix

        any_fixed = False
        if nfixed > 0:
            any_fixed = True
        fixed_stats['BPIXCORR'] = (any_fixed, 'True if any bad pixels were corrected')

        fixed_stats['BPIXNFIX'] = (nfixed, 'Number of bad pixels corrected')

        return newdata, fixed_stats
