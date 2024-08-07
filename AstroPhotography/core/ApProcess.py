"""
Contains the implementation of the ApProcess class.
"""

# 2024-05-02 dks : Initial implementation.

import sys
import os
import logging
from pathlib import Path
import math
import time
from datetime import datetime, timezone

import numpy as np
from ccdproc import ImageFileCollection
from astropy.table import Table
from astropy import wcs
import astropy
from astropy.io import fits

# AstroPhotography includes    
from .. import __version__
from . import ApFixBadPixels
from .ApFixCosmicRays import ApFixCosmicRays as ApFixCosmicRays
# TODO: I have no idea why the ApFixBadPixels line works, but the same for
# ApFixCosmicRays did not, instead requiring the more verbose line shown.

class ApProcess:
    """
    Astronomical image processor that processes multiple calibrated images
    in one or more filters to obtain astrometric solutions, combine images
    to improve signal-to-noise, and optionally generate false-color composite
    images.
    
    Intro
    -----
    
    This class provides functions to perform the following batch processing:
    
    - Star Detection on a directory or specified set of calibrated 
      images, see find_stars().
    - Generate astrometric solutions and WCS headers for files that
      have had star detection performed on them.
      
    In addition to processing the files this class also has utility
    methods to 
    - return ccdproc.ImaegFileCollection instances corresponding
      to existing files that correspond to various processing stages
    - remove (a.k.a. "clean") any existing files that ApProcessor
      generated
    - return lists of file names/paths that it will process, or will generate,
      even if those files don't already exist. These functions can be
      used before processing to check that it will process the files you
      expect.
    
    Star Detection and Image Navigation
    -----------------------------------
    
    Combining multiple images of the same source in the same filter requires 
    that each image have a valid World Coordinate System solution, based on 
    an astrometric solution to the stars in the field of view. The 
    `AstroPhotography` package also refers to this process as image navigation. 

    The process consists of finding the pixel locations brightest N stars 
    in each image, then calling `Astrometry.net` to compute a WCS solution 
    based on those star locations, and creating a new "navigated" fits file 
    that combines the calibrated file and the WCS solution. By default 
    2-dimensional Gaussian profiles are also fitted to a representative 
    number of stars per image to assess the Full Width at Half Maximum 
    of the star images (in pixels), the results of which are shown graphically 
    in a plot and numerically in a YaML file describing the stars found and 
    fitted per image. At the end of the process a summary of the image (star) 
    statistics is generated in CSV format that also incorporates the plate 
    scale (arcseconds per pixel) found in the astrometric solution, which 
    allows a human to look for outliers such as images where the seeing 
    or tracking was especially bad.
    
    In terms of pseudo-code find_stars() performs the following operations:
    
    .. code-block:: python
    
        for cal-img in calibrated-images:
            ap_find_stars cal-img --> source-list source-plot fwhm-plot quality-report.yml ds9-sources.reg
            ap_astrometry cal-img source-list --> navigated-img
        ap_quality_summary (all quality-report.yml files) --> quality-summary.csv
    
    File Types and File Naming Conventions
    --------------------------------------

    TBA
    """
    
    def __init__(self, loglevel):
        """
        Initializes an ApProcess instance
           
        Parameters
        ----------
        loglevel: str
                  Standard logging level string, e.g. `INFO`
        """
        
        self._name     = 'ApProcess'
        self._version  = __version__
        self._loglevel           = loglevel
        self._initialize_logger(self._loglevel)
        
        return
        
    def _check_file_exists(self, filename):
        """
        Raises an exception if the name file does not exist.
        
        Parameters
        ----------
        filename: str
            Name of file to check the existence of.
                  
        Raises
        ------
        RuntimeError
            If the named file is not found.
        """
        
        if not Path(filename).exists():
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return

    def _img_stats(self, data, label, verbose):
        """
        Calculate and optionally display some image statistics, returning
        a list of the minimum, maximum, mean and median values
        
        If verbose=True then the computed statistics are also written
        to the log at INFO level along with the specified informative
        label.
        
        Parameters
        ----------
        data: ndarray
            Data array to compute statistics of
        label : str
            Descriptive label, only used if verbose=True
        verbose : bool
            If true the the computed statistics are also written
            to the log at INFO level
            
        Returns
        -------
        minval :float
            Minimum value in NaN-filtered input data
        maxval : float
            Maximum value in NaN-filtered input data
        meanval : float
            Mean value in NaN-filtered input data
        medval : float
            Median value in NaN-filtered input data
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
        Initialize the logger
        
        Parameters
        ----------
        loglevel: str
                  Standard logging level string, e.g. `INFO`
        """
        
        self._logger = logging.getLogger(self._name)
        
        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: {}'.format(loglevel))
        self._logger.setLevel(numeric_level)
        self._logger.propagate = False
            
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
        return
        
    def _read_fits(self, image_filename, image_extension):
        """
        Read a single extension's data and header from a FITS file
        
        Parameters
        ----------
        image_filename: str
            Name of input files file
        image_extension : int or str
            Image extension number or name. For example extension 0 is
            the primary extension, or `srclist` would be an extension
            named sourcelist.
                  
        Returns
        -------
        ext_data : ndarray
            Data stored in the requested extension of the file
        ext_hdr : fits.Header
            FITS header of the requested extension of the file
            
        Notes
        -----
        - Unsigned integer handling is performed.
        - PEDESTAL values are removed from the returned data
        - Image scaling is not performed. 
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
        
        Parameters
        -----------
        hdr : fits.Header
            The FITS header that will be modified in place
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
        
        Parameters
        ----------
        inpdata_file : 
            Input FITS data file affected by bad pixels.
        ext_num : 
            Extension number for data array and header.
        outdata_file : 
            Bias/dark/flat corrected image.
        odata : 
            Bias/dark/flat field corrected data array.
        odict : 
            Additional FITS header keywords. The output file
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
        
    def get_file_names(self, processing_stage,
        data_dir, include_pattern=None, exclude_pattern=None, 
        input_file_list=None):
        """
        Returns a list of the file names that would correspond to a certain
        processing stage given the other input parameters
        
        The files may or may not exist already.
        
        The allowed processing stages are:
        
        - 'inputs': These are input files that will have star detection performed on them.
        - 'srclists': Star detection output FITS table files for each input file.
        
        Parameters
        ----------
        processing_stage : {'inputs'}
            The processing stage for which names should be returned.
            This should be one of the stage names described above.
        data_dir : str or path
                   Path to base directory containing FITS files to process,
                   e.g. `./`.
        include_pattern : str, optional
                   Globbing pattern of files we want included, specified
                   relative to data_dir. If not specified all FITS files
                   will be included. This parameter is ignored if file_list
                   is not None.
        exclude_pattern : str, optional
                   Globbing pattern of files we want excluded, specified
                   relative to data_dir. If not specified no FITS files
                   will be excluded. This parameter is ignored if file_list
                   is not None.
        input_file_list : list of str, optional
                   An explicit list of files, paths relative to data_dir,
                   may be specified. If provided only those files in input_file_list
                   are looked for, and the include and exclude patterns are
                   ignored.
                   
        Returns:
        fname_list : list of str
            File names that correspond to the named processing stage.
        
        Raises
        ------
        RuntimeException :
            If the processing_stage is not one of the allowed stages defined above
        """

        allowed_stages = ['inputs']
        if processing_stage not in allowed_stages:
            err_msg = f'Requested processing_stage {processing_stage} not one of the allowed values: {allowed_stages}'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
            
        self._logger.debug(f'Current working directory: {os.getcwd()}')
        self._logger.debug(f'Attempting to find stars in FITS files within directory={data_dir}')
        if input_file_list is None:
            # Use patterns
            self._logger.debug(f'Using files that match include_pattern="{include_pattern}"')
            self._logger.debug(f'Excluding files that match exclude_pattern="{exclude_pattern}"')
            ifc_cal = ImageFileCollection(data_dir, glob_include=include_pattern, glob_exclude=exclude_pattern)
        else:
            # Use explicit file list
            self._logger.info(f'Using specified file list: {input_file_list}')
            ifc_cal = ImageFileCollection(data_dir, filenames=input_file_list)
            
        fname_list = []
        for file in ifc_cal.summary['file']:
            fname_list.append( file )
        return fname_list
        
    def find_stars(self, data_dir, include_pattern=None, exclude_pattern=None, 
        file_list=None):
        """
        Runs ApFindStars on a set of files
        
        Parameters
        ----------
        data_dir : str or path
                   Path to base directory containing FITS files to process,
                   e.g. `./`.
        include_pattern : str, optional
                   Globbing pattern of files we want included, specified
                   relative to data_dir. If not specified all FITS files
                   will be included. This parameter is ignored if file_list
                   is not None.
        exclude_pattern : str, optional
                   Globbing pattern of files we want excluded, specified
                   relative to data_dir. If not specified no FITS files
                   will be excluded. This parameter is ignored if file_list
                   is not None.
        file_list : list of str, optional
                   An explicit list of files, paths relative to data_dir,
                   may be specified. If provided only those files in file_list
                   are looked for, and the include and exclude patterns are
                   ignored.
        """

        # Generate an image file collection
        keys = ['naxis1', 'naxis2', 'imagetyp', 'object', 'filter', 'exposure']
        
        self._logger.debug(f'Current working directory: {os.getcwd()}')
        self._logger.info(f'Attempting to find stars in FITS files within directory={data_dir}')
        if file_list is None:
            # Use patterns
            self._logger.info(f'Using files that match include_pattern="{include_pattern}"')
            self._logger.info(f'Excluding files that match exclude_pattern="{exclude_pattern}"')
            ifc_cal = ImageFileCollection(data_dir, keywords=keys, glob_include=include_pattern, glob_exclude=exclude_pattern)
        else:
            # Use explicit file list
            self._logger.info(f'Using specified file list: {file_list}')
            ifc_cal = ImageFileCollection(data_dir, keywords=keys, filenames=file_list)
        print(ifc_cal.summary)
        
        return
        
    def _find_stars_wrapper(self, p_loglevel, p_fitsimg, p_fitstbl, 
        p_extnum=0, p_search_fwhm=3.0, p_search_nsigma=7.0,
        p_detector_bitdepth=16, p_sat_frac=0.8, p_max_sources=200,
        p_nosatmask=True, p_plotfile=None, p_quiet=False,
        p_fwhm_plot=None, p_qual_rprt=None, p_regfile=None):
        """
        Wrapper for finding stars in a single image using ApFindStars.
        """

        # Perform initial source detection using default parameters.
        find_stars = ap.ApFindStars(p_fitsimg, p_extnum, p_search_fwhm,
            p_search_nsigma, p_detector_bitdepth, 
            p_max_sources, p_nosatmask, p_sat_frac, p_loglevel,
            p_plotfile, p_quiet)
        
        # Measure 2-Gaussian FWHM for select stars, get average over x and y
        (p_new_fwhm, p_madstd_fwhm, p_npts) = find_stars.measure_fwhm(p_fwhm_plot, 'both')
        
        # Refine source detection
        find_stars.source_search(p_new_fwhm, p_search_nsigma)
        
        # Re-run photometry
        find_stars.aperture_photometry()
        
        # As the source searching and photometry was redone, we should redo
        # the plotting.
        if p_plotfile is not None:
            find_stars.plot_image(p_plotfile)
        
        # Write optional quality report
        if p_qual_rprt is not None:
            find_stars.write_quality_report(p_qual_rprt)

        # Write optional ds9 format region file
        if p_regfile is not None:
            find_stars.write_ds9_region_file(p_regfile)
        
        # Write final sourcelist with photometry.
        find_stars.write_source_list(p_fitstbl)

        return
