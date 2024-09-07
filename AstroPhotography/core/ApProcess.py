"""
Contains the implementation of the ApProcess class.
"""

# 2024-05-02 dks : Initial implementation.
# 2024-08-12 dks : Numpy format documentation

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
from astropy.table import unique
from astropy import wcs
import astropy
from astropy.io import fits

# AstroPhotography includes    
from .. import __version__
from . import ApFixBadPixels
from .ApFixCosmicRays import ApFixCosmicRays as ApFixCosmicRays
import AstroPhotography.util as util
from .ApFindStars import ApFindStars as ApFindStars
from .ApAstrometry import ApAstrometry as ApAstrometry
from .ApQualitySummarizer import ApQualitySummarizer as ApQualitySummarizer

# ApFixCosmicRays did not, instead requiring the more verbose line shown.

class ApProcess:
    """
    Astronomical image processor that processes multiple calibrated FITS images
    in one or more filters to obtain astrometric solutions, combine images
    to improve signal-to-noise, and optionally generate false-color composite
    images.
    
    If you have raw FITS images you should use :class:`ApCalibrate` to
    calibrate the images before using this class.
    
    Intro
    -----
    
    This class provides functions to perform the following batch processing:

    - Generate astrometric solutions and WCS headers for all calibrated images
      in a directory, or a specified set of images. See :func:`navigate_images`.    
      This requires astrometric solutions be found (:class:`ApAstrometry` is used),
      which in turn depends on detecting star-like point sources in
      the images (:class:`ApFindStars` is used).
    - Stack a series of images (that have WCS headers) to improve
      the signal to noise ratio, by resampling them onto a common
      WCS projection. See :func:`resample_images_to_match`.
      
    In addition to processing the files this class also has utility
    methods to 
    
    - Return lists of file names/paths that it will generate, has generated,
      or that will exist that correspond to a given set of inputs. See
      :func:`get_file_names` and :func:`get_directories`.
      These functions can be
      used before processing to check that it will process the files you
      expect.
    - Reapply bad pixel/column/row removal and/or adding image header
      metadata in cases where the input calibrated images are deficient.
      See TBA and TBA. These functions should be run before image
      navigation and resampling.
    
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
    
    In terms of pseudo-code navigate_images() performs the following operations:
    
    .. code-block:: python
    
        for cal-img in calibrated-images:
            ap_find_stars cal-img --> source-list source-plot fwhm-plot quality-report.yml ds9-sources.reg
            ap_astrometry cal-img source-list --> navigated-img
        ap_quality_summary (all quality-report.yml files) --> quality-summary.csv
    
    Note that processing is sequential and single threaded *at this time*.
    Even if this changes for source detection the astrometry part will
    remain sequential to avoid spamming the ``Astrometry.net`` servers.
    
    Image Resampling And Mosaicing (Stacking)
    -----------------------------------------
    
    TBA
    xxx
    
    .. _File Types:
    
    File Types and File Naming Conventions
    --------------------------------------

    Input Files
    ~~~~~~~~~~~

    Input files must be FITS format images. By default the file extension
    ``.fits`` is assumed and used also for output images, but the functions
    allow the possibility of non-standard input file extensions such as
    ``.fit``, ``.ftz``, ``.fits.gz``, and so on.
    
    Input files are assumed to be calibrated, with dark and bias subtraction
    and flat fielding already performed. Cosmic ray and bad pixel/colum/row
    detection and removal should ideally have been performed, but this
    class allows for additional bad pixel/colum/row correction to be run.
    
    By default it is assumed that the calibrated images either came from
    iTelescope, or were processed by the AstroPhotography package from
    raw fits files. iTelescope provided calibrated images often have incomplete
    bad pixel/column/row correction, which is why this class has an option
    for performing that step. (Note that this requires you have access to
    iTelescope master dark files.
    
    Output files
    ~~~~~~~~~~~~
    
    Output files correspond to one of the following ``ap_filetype`` types:
    
    * Detected star-like sources and fitted source parameters are written to
      ``srclist`` files, which are fits files with internal binary tables.
      These are generated by :func:`navigate_images`, which calls the necessary
      :class:`.ApFindStars` functions under the hood.
      These files are required as inputs in finding an astrometric solution.
    * ``regfile``: Region files in ds9's image (pixel) coordinate system used by 
      `SAOImage ds9 <https://sites.google.com/cfa.harvard.edu/saoimageds9>`_.
    * ``plotfile``: A png format plot of the entire input image (with asinh-scaling)
      and the detected sources. For information/diagnostic purposes only.
    * ``qualfile``: A YaML format summary of the star detection results,
      used as input to :class:`.ApQualitySummarizer`.
    * ``fwhmplot``: A png format plot of a subset of detected star-like sources,
      showing a small region around each source and the fitted 2-dimensional
      gaussian parameters. For information/diagnostic purposes only.
    * ``navfile``: Navigated images are copies of the input calibrated images
      with valid WCS coordinate systems added to them, based on the astrometric
      solution obtained using the ``srclist`` files. These are the output of running
      :class:`.ApAstrometry`. 
      
    The function :func:`get_file_names` can be used to return the file names
    associated with a given ``ap_filetype``. `ApProcessor` recognizes
    three states (``ap_filestate``) for a given output file name, each of which can be 
    queried for.
    
    1. ``conceptual``: File names that would result from processing a given
       set of input files with `ApProcessor`. Processing may or may not have
       occurred.
    2. ``processed``: Output files associated with inputs that *this instance*
       of `ApProcessor` has generated or would have generated. (Some
       processing stages will not regenerate an existing output file if
       it already exists, but the existing output file is counted as having
       being processed by this `ApProcess` instance.)
    3. ``existing``: Output files that exist on disk that correspond to the
       specified input file parameters. These file need not have been
       processed with this `ApProcess` instance. The :func:`set_file_names`
       function can be used to make the current `ApProcess` instance aware
       of these files and treat them as if it had generated them itself.     
      
    Note
    ----
    
    - Astrometric solutions currently use ``Astrometry.net``. You will
      need to create an `Astrometry.net account <https://nova.astrometry.net/>`_ 
      and obtain your personal `API key <https://nova.astrometry.net/api_help>`_.
    - Image resampling and mosaicing uses Astromatic ``swarp``, and three
      color compositing uses Astromatic ``stiff``. These external programs
      must be installed on the host computer. Prebuilt packages for these
      exist in most major Linux distributions, or they can be downloaded
      and compiled with only moderate difficulty.
    """
    
    def __init__(self, loglevel):
        """
        Initializes an ApProcess instance
           
        Parameters
        ----------
        loglevel: str
                  Standard logging level string, e.g. `INFO`
        """
        
        self._name     = 'ApProcess'        # str : class name
        self._version  = __version__        # str : class version
        self._loglevel = loglevel           # str : Logging level
        self._initialize_logger(self._loglevel) 
        
        # Processing state variables
        self._status_table = None           # Table : processing status
        self._input_ifc    = None           # ImageFileCollection : Input files that were actually used
        self._data_dir     = None           # pathlib.Path : to input ``data_dir``
        self._extnum       = 0              # Extension number or name for input data
        
        # Existing file of the specified types
        self._srclist_files  = []
        self._regfile_files  = []
        self._plotfile_files = []
        self._qualfile_files = []
        self._fwhmplot_files = [] 
        self._navfile_files = []
        
        self.qual_pref = 'qual'
        self.qual_suff = '.yaml'
        
        return
        
    def _check_file_exists(self, filename, throws=True):
        """
        Check if a file or path exists, by default raising an exception 
        if it does not exist, also returning a boolean True if the file exists 
        
        Parameters
        ----------
        filename: str
            Name of file to check the existence of.
        throws : bool, optional, default=True
            If throws is True then an exception is thrown if the file path
            does not exist. If throw is False then the function does not
            throw exceptions at all, and will return False to the caller
            if the file path does not exist.
        
                  
        Returns
        -------
        exists : bool
            Boolean flag that is True if the file exists, and false if
            it does not exist and throw is False.
                
        Raises
        ------
        RuntimeError
            If the named file is not found and throw is True
        """
        exists = Path(filename).exists()
        if not exists and throws:
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return exists

    def _img_stats(self, data, label, verbose=False):
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
        verbose : bool, optional, default=False
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
            the primary extension, or 'srclist' would be an extension
            named ''srclist''.
                  
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
        
    def get_file_names(self, ap_filetype, data_dir, 
        ap_filestate='conceptual', include_pattern=None, exclude_pattern=None, 
        input_file_list=None, input_rootname=None, input_suffix='.fits',
        name_and_dir=False):
        """
        Returns a list of the file names that would correspond to a certain
        AstroPhotography file type given the other input parameters
        
        The files may or may not exist already. By default `ap_filestate='conceptual'`
        is used,meaning that the returned file names need not neccessarily exist,
        i.e. they are the names that would be generated for a given set of
        inputs. To return only the file names that have been generated by
        this instance of `ApProcess` use `ap_filestate='processed'`, and to 
        return the file names that exist irrespective of what and when
        they were created use `ap_filestate='existing'`. 
        See `File Types`_ for more information.
        
        The allowed AstroPhotography file types are:
        
        - ``input``: These are input files that will have star detection performed on them.
        - ``srclist``: Star detection output FITS table files for each input file.
        - ``regfile``: ds9-format region files.
        - ``plotfile``: PNG plots of the input images with the detected star-like sources
          plotted as circles.
        - ``qualfile``: YaML source detection quality summary files.
        - ``fwhmplot``: PNG plots zooming in around a subset of the detected
          sources in each input image.
        - ``navfile``: Copies of the input files but with valid WCS solutions.
        
        Parameters
        ----------
        ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile'}
            The AstroPhotography file type for which names should be returned.
            This should be one of the stage names described above.
        data_dir : str or path
            Path to base directory containing FITS files to process,
            e.g. `./`.
        ap_filestate : {'conceptual', 'processed', 'existing'}, optional, default='conceptual'
            This switch controls whether only file names corresponding
            to already existing files are returned. By default this is
            ``'conceptual'``, i.e. the file names returned may or may not already
            exist, but that processing all the input files would result
            in the named files being generated. 
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
        input_rootname : str, optional, default=None
            The part of the input files names that are shared and that
            designate them being the calibrated files. If not specified
            it is assumed that the files are from iTelescope or were 
            created by the AstroPhotography module itself, and will
            have input_rootnames of either 'Calibrated-iTelescope', 
            'calibrated', or just 'cal' .
            The output files from this function replace the rootname with
            a file-type specific prefix, as described in the class 
            documentation and get_file_names function documentation.
        input_suffix : str, optional, default='.fits'
            String denoting the file type suffix of the input image file.
            For example, '.fits' or '.fits' or '.ftz' or 'fits.gz' or '.fits.bz2'
        name_and_dir : bool, optional, default=False
            If False then only the file name is returned. If true, then
            output directory relative to data_dir is return as well, 
            even if it is only ``./``.    
                   
        Returns
        -------
        output_fname_list : list of str
            File names that would correspond to the named AstroPhotography file type.
            
        See Also
        --------
        :func:`get_directories` : 
            Returns the directory path in which files for given processing
            stage are, or will be.
        
        Raises
        ------
        RuntimeException :
            If the ap_filetype is not one of the allowed stages defined above
            
        See Also
        --------
        `File Types`_ : File type and file state documentation
        """
            
        # Returning the files that this instance has already processed?
        must_exist = False
        if 'processed' in ap_filestate:
            output_fname_list = self._get_processed_file_names(ap_filetype, data_dir, name_and_dir=False)
            return output_fname_list
        else:
            if 'existing' in ap_filestate:
                must_exist = True
            
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
 
        if not name_and_dir:
            odir = self.get_directories(ap_filetype, data_dir, False)
            self._logger.debug(f'Note that output files will be a subdirectory {odir}')
                        
        fname_list = self._get_file_names(ap_filetype, ifc_cal, input_rootname, input_suffix, name_and_dir)
        if must_exist:
            # Check file exists before adding it to output file name list
            output_fname_list = []
            for fname_val in fname_list:
                if not name_and_dir:
                    check_fname_val = odir + fname_val  # The actual file path
                else:
                    check_fname_val = fname_val         # Already include the directory
                if self._check_file_exists(check_fname_val, False):
                    output_fname_list.append( fname_val )
        else:
            output_fname_list = fname_list
        return output_fname_list

    def _get_processed_file_names(self, ap_filetype, data_dir, name_and_dir=False):
        """
        Return the files of the specified file type that this instance is
        aware of, either by processing them itself or by calls to set_file_names.
        
        Parameters
        ----------
        ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile'}
            The AstroPhotography file type for which names should be returned.
        data_dir : str or path
            Path to base directory containing FITS files to process,
            e.g. `./`.
        name_and_dir : bool, optional, default=False
            If False then only the file name is returned. If true, then
            output directory relative to data_dir is return as well, 
            even if it is only ``./``.    
                
        Returns
        -------
        output_fname_list : list of str
            File names that would correspond to the named AstroPhotography file type.
                    
        Raises
        ------
        RuntimeException :
            If the ap_filetype is not one of the allowed stages defined above
        """
        
        allowed_stages = ['input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile']
        if ap_filetype not in allowed_stages:
            err_msg = f'Requested ap_filetype {ap_filetype} not one of the allowed values: {allowed_stages}'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        
        # Get file names as a lists from the member variables
        if ap_filetype in 'input':
            ftype_list = []
            for file in self._input_ifc.summary['file']:
                ftype_list.append( file )
        elif ap_filetype in 'srclist':
            ftype_list = self._srclist_files
        elif ap_filetype in 'regfile':
            ftype_list = self._regfile_files
        elif ap_filetype in 'plotfile':
            ftype_list = self._plotfile_files
        elif ap_filetype in 'fwhmplot':
            ftype_list = self._fwhmplot_files
        elif ap_filetype in 'qualfile':
            ftype_list = self._qualfile_files
        elif ap_filetype in 'navfile':
            ftype_list = self._navfile_files
        
        # The member variables store the file name and directory
        # with respect to data_dir, so if name_and_dir is False then
        # the directory name must be removed.
        if not name_and_dir:
            output_fname_list = []
            odir = self.get_directories(ap_filetype, data_dir, False)
            dir_str_len = len(odir)
            for file in ftype_list:
                dir_start_pos = file.find(odir)
                output_fname_list.append( file[dir_start_pos + dir_str_len:] )
        else:
            output_fname_list = ftype_list

        return output_fname_list


    def _get_file_names(self, ap_filetype, ifc_cal, 
        input_rootname=None, input_suffix='.fits', name_and_dir=False):
        """
        Returns a list of the file names that would correspond to a certain
        AstroPhotography file type given the existing set of input files.
        
        The files may or may not exist already. This corresponds to the
        ap_filestate='conceptual'.
        
        The allowed AstroPhotography file types are:
        
        - ``input``: These are input files that will have star detection performed on them.
        - ``srclist``: Star detection output FITS table files for each input file.
        - ``regfile``: ds9-format region files.
        - ``plotfile``: PNG plots of the input images with the detected star-like sources
          plotted as circles.
        - ``qualfile``: YaML source detection quality summary files.
        - ``fwhmplot``: PNG plots zooming in around a subset of the detected
          sources in each input image.
        - ``navfile``: Copies of the input files but with valid WCS solutions.
        
        Parameters
        ----------
        ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile'}
            The AstroPhotography file type for which names should be returned.
            This should be one of the stage names described above.
        ifc_cal : ccdproc.ImageFileCollection
            The ImageFileCollection corresponding to the `inputs`-type files.
        input_rootname : str, optional, default=None
            The part of the input files names that are shared and that
            designate them being the calibrated files. If not specified
            it is assumed that the files are from iTelescope or were 
            created by the AstroPhotography module itself, and will
            have input_rootnames of either 'Calibrated-iTelescope', 
            'calibrated', or just 'cal' .
            The output files from this function replace the rootname with
            a file-type specific prefix, as described in the class 
            documentation and get_file_names function documentation.
        input_suffix : str, optional, default='.fits'
            String denoting the file type suffix of the input image file.
            For example, '.fits' or '.fits' or '.ftz' or 'fits.gz' or '.fits.bz2'
        name_and_dir : bool, optional, default=False
            If False then only the file name is returned. If true, then
            output directory relative to data_dir is return as well, 
            even if it is only ``./``.    
                   
        Returns
        -------
        outputs_fname_list : list of str
            File names that would correspond to the named AstroPhotography file type.
            
        See Also
        --------
        :func:`get_directories` : 
            Returns the directory path in which files for given processing
            stage are, or will be.
        
        Raises
        ------
        RuntimeException :
            If the ap_filetype is not one of the allowed stages defined above
            
        See Also
        --------
        `File Types`_ : File type documentation
        """

        allowed_stages = ['input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile']
        if ap_filetype not in allowed_stages:
            err_msg = f'Requested ap_filetype {ap_filetype} not one of the allowed values: {allowed_stages}'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
                        
        outputs_fname_list = []
        for file in ifc_cal.summary['file']:
            ofile, output_dir = util.namefn_calibrated_input(file, input_rootname, input_suffix, ap_filetype)
            if name_and_dir:
                outputs_fname_list.append( output_dir + ofile )
            else:
                outputs_fname_list.append( ofile )
            
        return outputs_fname_list
        
    def get_directories(self, ap_filetype, data_dir, absolute=False):
        """
        Returns the path to the directory that will or does hold files
        corresponding to one of the AstroPhotography file types
        
        The directory may or may not exist already.
        
        The allowed file type are:
        
        - ``input``: These are input files that will have star detection performed on them.
        - ``srclist``: Star detection output FITS table files for each input file.
        - ``regfile``: ds9-format region files.
        - ``plotfile``: PNG plots of the input images with the detected star-like sources
          plotted as circles.
        - ``qualfile``: YaML source detection quality summary files.
        - ``fwhmplot``: PNG plots zooming in around a subset of the detected
          sources in each input image.
        - ``navfile``: Copies of the input files but with valid WCS solutions.
        
        Parameters
        ----------
        ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile'}
            The AstroPhotography file type for which names should be returned.
            This should be one of the stage names described above.
        data_dir : str or path
            Path to base directory containing FITS files to process,
            e.g. `./`.
        absolute : bool, optional, default=False
            If True then the absolute directory path is returned, instead
            of the path relative to ``data_dir``.
                  
        Returns
        -------
        output_dir : str 
            Directory name in which output_files will be written, relative
            to the root directory established by the input files.
        """
        
        allowed_stages = ['input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile']
        if ap_filetype not in allowed_stages:
            err_msg = f'Requested ap_filetype {ap_filetype} not one of the allowed values: {allowed_stages}'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
            
        # Relative path
        output_dir = util.namefn_getdir(ap_filetype)
        self._logger.debug(f'For ap_filetype {ap_filetype} the relative directory path is {output_dir}')
        if absolute:
            p = Path(data_dir) / Path(output_dir)
            output_dir = str(p.absolute())
            self._logger.debug(f'For ap_filetype {ap_filetype} the absolute directory path is {output_dir}')
        return output_dir
        
        return odir
        
    def _generate_all_output_names(self, inputfile, input_rootname=None, 
        input_suffix='.fits', name_and_dir=False, mkdir=False):
        """
        Given an input file name and a conversion dictionary, generate all the file names that might
        be used to run star finding and astrometry on that image.
        
        Parameters
        ----------
        inputfile : str
            Name of the input calibrated file being processed
        input_rootname : str, optional, default=None
            The part of the input files names that are shared and that
            designate them being the calibrated files. If not specified
            it is assumed that the files are from iTelescope or were 
            created by the AstroPhotography module itself, and will
            have input_rootnames of either 'Calibrated-iTelescope', 
            'calibrated', or just 'cal' .
            The output files from this function replace the rootname with
            a file-type specific prefix, as described in the class 
            documentation and get_file_names function documentation.
        input_suffix : str, optional, default='.fits'
            String denoting the file type suffix of the input image file.
            For example, '.fits' or '.fits' or '.ftz' or 'fits.gz' or '.fits.bz2'
        name_and_dir : bool, optional, default=False
            If False then only the file name is returned. If true, then
            output directory relative to data_dir is return as well, 
            even if it is only ``./``.    
        mkdir : bool, optional, default=False
            If True then the output directory will be generated if it does
            not already exist.
    
        Returns
        -------
        f_srclist : str
            File name or path to the FITS table of detected star parameters
        f_regfile : str
            File name or path to the ds9-format region file.
        f_plotfile : str
            File name or path to the PNG plot of the image with detected
            stars shown in circles.
        f_fwhmfile : str
            File name or path to the PNG cut-outs of selected stars and
            fitted parameters 
        f_qualfile : str
            File name or path to the YaML star detection summary file 
        f_navfile : str
            File name or path to the navigated copy of the input calibration
            file, now including a valid WCS header.
        """
        
        oname_type = ['srclist', 'regfile', 'plotfile', 'fwhmplot', 'qualfile', 'navfile']
        ofile_dict = {}
        for ap_filetype in oname_type:                
            ofile, output_dir = util.namefn_calibrated_input(inputfile, input_rootname, input_suffix, ap_filetype)

            if mkdir:
                if not os.access(output_dir, os.F_OK):
                    try:
                        os.mkdir(output_dir)
                    except:
                        self._logger.error(f'Error, mkdir failed trying to create {output_dir}')
                        raise
                    self._logger.info(f'Successfully created subdirectory {output_dir}')
                else:
                    self._logger.debug(f'Subdirectory {output_dir} already exists.')
                
            if name_and_dir:
                ofile = output_dir + ofile
            ofile_dict[ap_filetype] = ofile
        
        return ofile_dict['srclist'], ofile_dict['regfile'], ofile_dict['plotfile'], ofile_dict['fwhmplot'], ofile_dict['qualfile'], ofile_dict['navfile']
        
    def set_file_names(self, data_dir, include_pattern=None, exclude_pattern=None, 
        input_file_list=None, input_rootname=None,  input_suffix='.fits'):
        """
        Sets the input and output file names associated with `ApProcess`
        and finds all the existing files matching those prescription.
        
        This function can be used when picking up processing of a set of
        images that was partially completed at an earlier date using a
        different instance of `ApProcess`.
        
        The data_dir, include_pattern, exclude_pattern,
        and filelist parameters are used to specify which input files to
        process. The input_rootname and input_suffix parameters are *not*
        used to select inputs but instead are used to generate the 
        expected output file names based on the input file names.
        
        Once this has been done existing files matching those file
        names will be found and added to the internal member function
        file lists.
        
        Parameters
        ----------
        data_dir : str or path
            Path to base directory containing FITS files to process,
            e.g. `./`.
        include_pattern : str, optional
            Globbing pattern for input files we want included, specified
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
        input_rootname : str, optional, default=None
            The part of the input files names that are shared and that
            designate them being the calibrated files. If not specified
            it is assumed that the files are from iTelescope or were 
            created by the AstroPhotography module itself, and will
            have input_rootnames of either 'Calibrated-iTelescope', 
            'calibrated', or just 'cal' .
            The output files from this function replace the rootname with
            a file-type specific prefix, as described in the class 
            documentation and get_file_names function documentation.
        input_suffix : str, optional, default='.fits'
            String denoting the file type suffix of the input image file.
            For example, '.fits' or '.fits' or '.ftz' or 'fits.gz' or '.fits.bz2'
        """
        
        # Generate an image file collection
        keys = ['naxis1', 'naxis2', 'imagetyp', 'object', 'filter', 'exposure']

        self._logger.info(f'Finding input and output files starting in directory={data_dir}')        
        self._logger.debug(f'Current working directory: {os.getcwd()}')
        if input_file_list is None:
            # Use patterns
            self._logger.info(f'Using files that match include_pattern="{include_pattern}"')
            self._logger.info(f'Excluding files that match exclude_pattern="{exclude_pattern}"')
            ifc_cal = ImageFileCollection(data_dir, keywords=keys, glob_include=include_pattern, glob_exclude=exclude_pattern)
        else:
            # Use explicit file list
            self._logger.info(f'Using specified file list: {input_file_list}')
            ifc_cal = ImageFileCollection(data_dir, keywords=keys, filenames=input_file_list)
        self._input_ifc = ifc_cal
        self._data_dir  = Path(data_dir)
                
        # Generate status table, written to self._status_table
        num_inputs = len(ifc_cal.summary)
        self._make_status_table(num_inputs)
        
        self._logger.info(f'There are {num_inputs} input files matching the parameters given.')
        self._logger.debug(f'Input file collection:\n{ifc_cal}')
        
        # Now look for pre-existing files matching the expected names,
        # and add them to the file lists. Clear the existing lists first...
        self._logger.info('Looking for existing files matching name conventions...')
        name_and_dir = True
        oname_type = ['srclist', 'regfile', 'plotfile', 'fwhmplot', 'qualfile', 'navfile']
        for ftype in oname_type:
            # Clear old values
            if ftype in 'srclist':
                ftype_list = self._srclist_files
            if ftype in 'regfile':
                ftype_list = self._regfile_files
            if ftype in 'plotfile':
                ftype_list = self._plotfile_files
            if ftype in 'fwhmplot':
                ftype_list = self._fwhmplot_files
            if ftype in 'qualfile':
                ftype_list = self._qualfile_files
            if ftype in 'navfile':
                ftype_list = self._navfile_files
            ftype_list.clear()
            
            self._logger.debug(f'Looking for existing {ftype} files.')
            flist = self._get_file_names(ftype, self._input_ifc, input_rootname, input_suffix, name_and_dir)
            for fname_val in flist:
                if self._check_file_exists(fname_val, False):
                    if ftype in 'srclist':
                        self._srclist_files.append( fname_val )
                        continue
                    if ftype in 'regfile':
                        self._regfile_files.append( fname_val )
                        continue
                    if ftype in 'plotfile':
                        self._plotfile_files.append( fname_val )
                        continue
                    if ftype in 'fwhmplot':
                        self._fwhmplot_files.append( fname_val )
                        continue
                    if ftype in 'qualfile':
                        self._qualfile_files.append( fname_val )
                        continue 
                    if ftype in 'navfile':
                        self._navfile_files.append( fname_val )
                        continue        
            numfiles = len(ftype_list)
            self._logger.info(f'Found the following {numfiles} {ftype} files: {ftype_list}')
        return
        
    def navigate_images(self, data_dir, include_pattern=None, exclude_pattern=None, 
        input_file_list=None, input_rootname=None,  input_suffix='.fits',
        final_quality_file=None,
        clean_star_detection=False, clean_astrometry=False, stop_on_error=False,
        extnum=0, search_fwhm=3.0,
        search_nsigma=7.0, detector_bitdepth=16, 
        max_sources=200, nosatmask=True, sat_frac=0.8, quiet=True,
        srclist_extname='AP_XYPOS',
        astnet_key=None,
        use_sip=False,
        user_scale=None,
        scale_err_ratio=None):
        """
        Runs ApFindStars on a set of files
        
        Given a set of input files, star detection and astrometric solutions
        are computed for each file. The data_dir, include_pattern, exclude_pattern,
        and filelist parameters are used to specify which input files to
        process. The input_rootname and input_suffix parameters are *not*
        used to select inputs but instead are used to generate output file
        names based on the input file names.
        
        By default, if the expected output file for a given input file 
        already exists it is not regenerated and that stage of processing
        is skipped to save time. For star detection the key output file
        is the ``srclist`` output, and for astrometry to key output is
        the ``navfile``. The presence or absence of other output file
        types, e.g. ``fwhmplot`` files, is ignored.
        
        To force ``navigate_images`` to rerun star detection and/or 
        astrometry even in the presence of existing output files,
        the user can set either ``clean_star_detection`` and/or ``clean_astrometry``
        to ``True``.

        
        Parameters
        ----------
        data_dir : str or path
            Path to base directory containing FITS files to process,
            e.g. `./`.
        include_pattern : str, optional
            Globbing pattern for input files we want included, specified
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
        input_rootname : str, optional, default=None
            The part of the input files names that are shared and that
            designate them being the calibrated files. If not specified
            it is assumed that the files are from iTelescope or were 
            created by the AstroPhotography module itself, and will
            have input_rootnames of either 'Calibrated-iTelescope', 
            'calibrated', or just 'cal' .
            The output files from this function replace the rootname with
            a file-type specific prefix, as described in the class 
            documentation and get_file_names function documentation.
        input_suffix : str, optional, default='.fits'
            String denoting the file type suffix of the input image file.
            For example, '.fits' or '.fits' or '.ftz' or 'fits.gz' or '.fits.bz2'
        final_quality_file : str, optional, default=None
            If specified, this is the name for a quality summary file
            in CSV format 
            generated using :class:`ApQualitySummarizer`, which summarizes
            all YaML quality files found in the standard quality
            summary file directory (typically `./MetaData/`).
        clean_star_detection : bool, default=False
            If True then **existing** star detection outputs (``srclist``,
            ``regile``, ``plotfile``, ``qualfile``, and ``fwhmplot`` files)
            will be deleted and regenerated.
        clean_astrometry : bool, default=False
            If True then **existing** astrometry outputs (``navfile``files)
            will be deleted and regenerated.
        stop_on_error : bool, optional, default=False
            By default :func:`navigate_images` will continue attempting to
            process all the input files even if a given stage fails on one
            of the inputs. That failure will be recorded in the output
            ``status_table``. If you want it to instead stop processing
            when a failure is detected, set stop_on_error to True. 
                   
        Other Parameters
        ----------------
        
        extnum : int or str, optional, default=0
            Extension number or name for the extension holding the image data.
            Usually this is 0, for the ``PrimaryHDU``. (ApFindStars parameter)
        search_fwhm : float, optional, default=3.0
            Initial guess or estimate of the stellar PSF FWHM in pixels
            in this image (ApFindStars parameter)
        search_nsigma : float, optional, default=7.0
            Minimumn number of sigma above background for a detection
            (ApFindStars parameter)
        detector_bitdepth : int, optional, default=16
            Detector bit-depth, used in estimating which pixels are saturated.
            (16 for most CCDs, ?? for CMOS, ?? for camera)  (ApFindStars parameter)
        max_sources : int, optional, default=200
            Maximum number of sources to output or None. If not None
            then only the max_sources brightest sources will be output.
            Typically there is no advantage to having very large numbers
            of detected sources when performing astrometry, and may
            well slow it down. I normally use max_sources=200.
            (ApFindStars parameter)
        nosatmask : bool, optional, default=True
            If True, keep possibly saturated stars. (ApFindStars parameter)
        sat_frac : float, optional, default=0.8
            Fraction of full well at which we assume star saturated 
            (ApFindStars parameter)
        quiet : bool, optional, default=True
            If True this suppresses the runtime source list printing 
            to STDOUT (ApFindStars parameter)
        srclist_extname : str, optional, default='AP_XYPOS'
            FITS extension name for star X,Y position data. Default=``AP_XYPOS``
            (ApAstrometry parameter)
        astnet_key : str, optional, default=None
            Your personal Astrometry.net API key, if you have not already
            added it to your ~/.astropy/config/astroquery.cfg config file.
            (ApAstrometry parameter)
        use_sip : bool, optional, default=False
            Allow astrometry.net to fit SIP polynomial distortion terms.
            This may be necessary for very large fields of view (>10 deg),
            but SIP is not treated correctly by swarp (and possibly other
            software).
            (ApAstrometry parameter)
        user_scale : float, optional, default=None
            If specified, override the estimate plate scale in the source list file
            and instead use a user spacified estimate of the plate scale.
            The units are arcseconds/pixel.
            (ApAstrometry parameter)
        scale_err_ratio: float, optional, default=None
            The relative uncertainty in the estimated plate scale,
            expressed as a ratio. This applies to either the default 
            estimate from the source list, or a user-supplied plate
            scale. For example, if the estimate plat scale is 2.0 arcsec/pix
            and the scale_err_ratio=1.5, then the plate scale range that
            will be search by Astrometry.net is 2/1.5 (=4/3) to 2*1.5 (=3)
            arcseconds. If not specified ApAstrometry will use a value of 1.3.
            Using a larger value can help in cases where astrometric
            solutions fail, for example if incorrect telescope metadata
            leads to inaccurate estimated plate scales.
            (ApAstrometry parameter)
                   
        Returns
        -------
        status_table : astropy.table
            An astropy table containing the processing status of the images
            
        See Also
        --------
        get_file_names : Given input file names returns the file names associated
            with a given output file type.
        """

        # Generate an image file collection
        keys = ['naxis1', 'naxis2', 'imagetyp', 'object', 'filter', 'exposure']
        
        self._logger.debug(f'Current working directory: {os.getcwd()}')
        self._logger.info(f'Attempting to find stars in FITS files within directory={data_dir}')
        if input_file_list is None:
            # Use patterns
            self._logger.info(f'Using files that match include_pattern="{include_pattern}"')
            self._logger.info(f'Excluding files that match exclude_pattern="{exclude_pattern}"')
            ifc_cal = ImageFileCollection(data_dir, keywords=keys, glob_include=include_pattern, glob_exclude=exclude_pattern)
        else:
            # Use explicit file list
            self._logger.info(f'Using specified file list: {input_file_list}')
            ifc_cal = ImageFileCollection(data_dir, keywords=keys, filenames=input_file_list)
        self._input_ifc = ifc_cal
        self._data_dir  = Path(data_dir)
        self._extnum    = extnum
                
        # Generate status table, written to self._status_table
        num_inputs = len(ifc_cal.summary)
        self._make_status_table(num_inputs)
        
        self._logger.info(f'There are {num_inputs} input files matching the parameters given.')
        self._logger.debug(f'Input file collection:\n{ifc_cal}')
        
        # Iterate over the input files
        idx = 0
        name_and_dir = True     # We want _generate_all_output_names to include the sub-directory
        mkdir = True            # We want subdirectories made if not present
        proc_tstart = time.perf_counter()
        for hdu, fname in ifc_cal.hdus(return_fname=True):
            self._logger.debug(80*'-')
            self._logger.info(f'Processing input file {fname}')

            f_srclist, f_regfile, f_plotfile, f_fwhmplot, f_qualfile, f_navfile = self._generate_all_output_names(fname,
                 input_rootname, input_suffix, name_and_dir, mkdir)

            # Determine whether to perform star detection, based on p_clean and presence of
            dont_throw = False
            does_srclist_exist = self._check_file_exists(f_srclist, dont_throw)
            if clean_star_detection or not does_srclist_exist:    
                if clean_star_detection:
                    # Remove all find star outputs, not just the srclist
                    for file_to_remove in [f_srclist, f_regfile, f_plotfile, f_fwhmplot, f_qualfile]:
                        if self._check_file_exists(file_to_remove, dont_throw):
                            self._logger.debug(f'Removing existing {file_to_remove} as clean_star_detection={clean_star_detection}')
                            Path(file_to_remove).unlink()
                            
                self._logger.info(f'Finding stars, output source list: {f_srclist}')
                fs_tstart = time.perf_counter()
                status = 'pass'
                try:
                    self._find_stars_wrapper(fname, f_srclist, 
                        extnum, search_fwhm, search_nsigma,
                        detector_bitdepth, sat_frac, max_sources,
                        nosatmask, f_plotfile, quiet,
                        f_fwhmplot, f_qualfile, f_regfile)
                except:
                    self._logger.warning(f'Error, caught exception when processing {fname}')
                    status = 'fail'
                    if stop_on_error:
                        raise
                
                fs_tend     = time.perf_counter()
                fs_telapsed = fs_tend - fs_tstart # seconds
                fs_status   = status
            else:
                self._logger.info(f'Skipping star detection because {f_srclist} exists and clean_star_detection={clean_star_detection}')
                fs_status = 'skipped_exists'
                fs_telapsed = 0
                
            # Update lists of existing files of each type we've just generated
            # in source searching
            self._update_file_lists( {'srclist': f_srclist,
                'regfile': f_regfile,
                'plotfile': f_plotfile,
                'qualfile': f_qualfile,
                'fwhmplot': f_fwhmplot} )
            
            # Determine whether to perform astrometry, based on p_clean 
            # and presence of srclist and navfile
            does_navfile_exist =  self._check_file_exists(f_navfile, dont_throw)
            does_srclist_exist =  self._check_file_exists(f_srclist, dont_throw)
            if does_srclist_exist:
                if clean_astrometry or not does_navfile_exist: 
                    if clean_astrometry:
                        for file_to_remove in [f_navfile]:
                            if self._check_file_exists(file_to_remove, dont_throw):
                                self._logger.debug(f'Removing existing {file_to_remove} as clean_star_detection={clean_astrometry}')
                                Path(file_to_remove).unlink()
                    
                    ast_status = 'pass'
                    ast_tstart = time.perf_counter()
                    self._logger.info(f'Performing astrometry, output navigated image: {f_navfile}')
                    try:
            
                        ap_astrom = ApAstrometry(fname, 
                            f_srclist, 
                            f_navfile, 
                            inp_img_extnum=extnum,
                            srclist_extname=srclist_extname,
                            astnet_key=astnet_key,
                            use_sip=use_sip,
                            user_scale=user_scale,
                            scale_err_ratio=scale_err_ratio,
                            loglevel=self._loglevel)
                        p_status = ap_astrom.status()
                        self._logger.debug(f'ApAstrometry return status: {p_status}')
                        if p_status == ApAstrometry.NOMINAL:
                            ast_status = 'pass'
                        elif p_status == ApAstrometry.INPUT_ERROR:
                            ast_status = 'input_error'
                        else:
                            ast_status = 'fail'
                    except:
                        self._logger.warning(f'Error, caught exception when performing astrometry on {fname}')
                        ast_status = 'exception'
                        if stop_on_error:
                            raise
            
                    ast_tend     = time.perf_counter()
                    ast_telapsed = ast_tend - ast_tstart # seconds
                    ast_status   = status
                else:
                    self._logger.info(f'Skipping astrometry because {f_navfile} exists and clean_astrometry={clean_astrometry}')
                    ast_status = 'skipped_exists'
                    ast_telapsed = 0
                    
                self._update_file_lists( {'navfile': f_navfile} )
            else:
                # No source list
                self._logger.error(f'Skipping astrometry because {f_srclist} does not exist.')
                ast_status = 'skipped_no_srclist'
                ast_telapsed = 0

            # Fill in run info for this file
            result_tuple = (fname, fs_status, fs_telapsed, ast_status, ast_telapsed)
            self._status_table[idx] = result_tuple

            # update idx
            idx += 1
            
        proc_tend = time.perf_counter()
        proc_telapsed = proc_tend - proc_tstart
        self._logger.info(f'Finished processing {idx+1} files in {proc_telapsed:.3f} seconds.')
        
        # Generate overall quality file summary
        qual_file_dir = self.get_directories('qualfile', data_dir, absolute=False)
        self.create_quality_summary(qual_file_dir, final_quality_file)
        
        return self._status_table
        
    def create_quality_summary(self, qual_file_dir, quality_summary_file):
        """
        Runs :class:`ApQualitySummarizer` on any quality YaML files in the
        specified directory, writing a summary CSV to the main ``data_dir``.
        
        Notes
        -----
        :class:`ApQualitySummarizer` processes all quality files it finds in
        the specified directory, not just the ones generated by this instance
        of `ApProcess`.
        
        Parameters
        ----------
        qual_file_dir : str
            Directory path that contains the individual quality YaML files
            output by :func:`navigate_images`.
        quality_summary_file : str
            Name for the output quality summary CSV file.
        """
        
        walk_tree = False
        summarizer = ApQualitySummarizer(qual_file_dir, quality_summary_file,
            self._loglevel, walk_tree,
            self.qual_pref, self.qual_suff)
        return
        
    def _update_file_lists(self, added_file_dict):
        """
        Updates the internal lists of output files that were processed
        by this instance, or would have been processed but were skipped
        as they already exist.
        
        These file lists can be accessed by callers using :func:`get_file_names`
        by passing the function parameter ``ap_filestate='processed'``.
        
        Parameters
        ----------
        added_file_dict : dict(str, str)
            Dictionary of consisting of file type keys and file name
            values.
        """
        
        oname_type = ['srclist', 'regfile', 'plotfile', 'fwhmplot', 'qualfile', 'navfile']
        for key, val in added_file_dict.items():
            if key in oname_type:
                if self._check_file_exists(val, False):
                    if key in 'srclist':
                        self._srclist_files.append( val )
                        continue
                    if key in 'regfile':
                        self._regfile_files.append( val )
                        continue
                    if key in 'plotfile':
                        self._plotfile_files.append( val )
                        continue
                    if key in 'fwhmplot':
                        self._fwhmplot_files.append( val )
                        continue
                    if key in 'qualfile':
                        self._qualfile_files.append( val )
                        continue 
                    if key in 'navfile':
                        self._navfile_files.append( val )
                        continue
            else:
                self._logger.warning(f'Unexpected file type key ({key}) supplied to _update_file_lists')
        
        return
        
    def _find_stars_wrapper(self, a_fitsimg, a_fitstbl, 
        an_extnum=0, a_search_fwhm=3.0, a_search_nsigma=7.0,
        a_detector_bitdepth=16, a_sat_frac=0.8, a_max_sources=200,
        do_nosatmask=True, a_plotfile=None, do_quiet=False,
        a_fwhm_plot=None, a_qual_rprt=None, a_regfile=None):
        """
        Wrapper for finding stars in a single image using ApFindStars.
        
        Parameters
        ----------
        a_fitsimg : str 
            Name of the input FITS image to search for star-like sources.
        a_fitstbl : str
            Name of output FITS file containing detected source parameters.
            This is a FITS binary table file that is used along with the input
            image by :class:`ApAstrometry`.
        an_extnum : int or str, default=0 
            Extension number or name for the extension holding the image data. Usually this is 0, for the ``PrimaryHDU``.
        a_search_fwhm : float, default=3.0 
            Initial guess or estimate of the stellar PSF FWHM in pixels
            in this image
        a_search_nsigma : float, default=7.0
            Minimumn number of sigma above background for a detection
        a_detector_bitdepth : int, default=16 
            Detector bit-depth, used in estimating which pixels are 
            saturated. (16 for most CCDs, ?? for CMOS, ?? for camera)
        a_sat_frac : float, default=0.8 
            Fraction of full well at which we assume star saturated
        a_max_sources : int, default=200
            Maximum number of sources to output or None. If not None then 
            only the max_sources brightest sources will be output. 
            Typically there is no advantage to having very large numbers 
            of detected sources when performing astrometry, and may well 
            slow it down. I normally use max_sources=200.
        do_nosatmask : bool : default=True 
            If True, keep possibly saturated stars.
        a_plotfile : str, default=None 
            If not None then this is the name for an output PNG plot of the image with detected sources plotted as circles.
        do_quiet : bool, default=False
            If True this suppresses the runtime source list printing to STDOUT
        a_fwhm_plot : str, default=None 
            If not None then a PNG plot zooming in around a subset of the detected
            sources, and their best-fit parameters, is generated.
        a_qual_rprt : str, default=None 
            If not None then a YaML file summarizing the source detection
            outputs is generated.
        a_regfile : str, default=None
            If not None a ds9-format region file will be generated with
            the coordinates of the detected stars (in pixel coordinates).
        
        See Also
        --------
        :class:`ApFindStars`
        """

        # Perform initial source detection using default parameters.
        find_stars = ApFindStars(a_fitsimg, an_extnum, a_search_fwhm,
            a_search_nsigma, a_detector_bitdepth, 
            a_max_sources, do_nosatmask, a_sat_frac, self._loglevel,
            a_plotfile, do_quiet)
        
        # Measure 2-Gaussian FWHM for select stars, get average over x and y
        # We don't generate a plot for this step because we will rerun this
        # function later with updated source search results
        (a_new_fwhm, a_madstd_fwhm, a_npts) = find_stars.measure_fwhm(None, 'both')
        self._logger.info(f'Initial star detection and fitting: FWHM={a_new_fwhm:.3f} +/- {a_madstd_fwhm:.3f} pixels using {a_npts} stars.')
        
        # Refine source detection
        self._logger.debug(f'Updating source searching using initial FWHM={a_new_fwhm:.3f} +/- {a_madstd_fwhm:.3f} pixels using {a_npts} stars.')
        find_stars.source_search(a_new_fwhm, a_search_nsigma)
        
        # Re-run photometry
        find_stars.aperture_photometry()
        
        # Measure 2-Gaussian FWHM for select stars, get average over x and y
        (a_new_fwhm2, a_madstd_fwhm2, a_npts2) = find_stars.measure_fwhm(a_fwhm_plot, 'both')
        self._logger.debug(f'Final star detection and fitting: FWHM={a_new_fwhm2:.3f} +/- {a_madstd_fwhm2:.3f} pixels using {a_npts2} stars.')
        
        # As the source searching and photometry was redone, we should redo
        # the plotting.
        if a_plotfile is not None:
            find_stars.plot_image(a_plotfile)
        
        # Write optional quality report
        if a_qual_rprt is not None:
            find_stars.write_quality_report(a_qual_rprt)

        # Write optional ds9 format region file
        if a_regfile is not None:
            find_stars.write_ds9_region_file(a_regfile)
        
        # Write final sourcelist with photometry.
        find_stars.write_source_list(a_fitstbl)
        if self._check_file_exists(a_fitstbl, True):
            self._logger.info(f'Confirming output srclist {a_fitstbl} was created.')
        else:
            self._logger.error(f'Expected output srclist {a_fitstbl} not found.')
        return
        
    def resample_images_to_match(self, target_wcs_file, filter_list=None):
        """
        Resample and combine all navigated images to match the WCS defined
        in the specified file, separating outputs by FILTER.

        By default all filters found in the navigated images will be 
        processed, but a user-specified list of filters can be specified 
        instead.
        
        The function parameters ``resampled_dir``, ``resampled_prefix``,
        and ``resampled_suffix`` control the location and name of the 
        resampled stacked output images.
        The output resampled stacked images are placed in the 
        subdirectory ``resampled_dir`` with file names of the form
        ``<resampled_prefix>_<filter><resampled_suffix>``.
                
        Image Selection
        ~~~~~~~~~~~~~~~
                
        Image must have valid WCS headers in order to be resampled and
        stacked, either from external sources or generated by 
        `ApProcess.navigate_images`.
        
        By default this function will process ``navfile`` type files
        associated with earlier processing by this ApProcess instance,
        most likely by `ApProcess.navigate_images`.
        
        If you wish to resample pre-existing navigated images there are
        two ways to specify which files to process:
        
        * If the files were created by older runs of `ApProcess` and/or 
          follow the `ApProcess` file naming conventions you should first
          call :func:`set_file_names`. This will populate the internal
          file lists with all existing ``navfile``-type files matching
          the specified patterns. A subsequent call to `resample_images_to_match`
          will use those files.
        * If the navigated files do not follow the `ApProcess` file naming 
          conventions, for example they were generated by other software,
          then the simplest solution is to use the ``navigated_file_list``
          parameter to explicitly specify which files to use.
        
        See Also
        --------
        get_filter_list : 
            Returns a list of the filters found in the current image set
        set_file_names :
            Set the file name parameters for `ApProcess` and find all the
            existing files than match the expected file names.
        get_file_names :
            Returns the file names of different processing stages associated
            with `ApProcess`.
        
        Parameters
        ----------
        target_wcs_file : str
            File name, including path, to the FITS file than contains
            the WCS that all navigated images should be resampled to
            match.
        filter_list : list of str or None, optional, default=None
            If specified then only images with ``FILTER`` header keywords
            matching one of the input list strings will be resampled. For
            example, to only resample the H-alpha and luminance images
            the filter list would be ``['Lum', 'Ha']``
        navigated_file_list : list of str or None, optional
            An explicit list of navigated images, paths relative to data_dir,
            may be specified. If provided only those files in navigated_file_list
            are used.
            
        Returns
        -------
        info_table : astropy.table.Table
            Informational table of filters, number of navigated images,
            exposure time stats, and optionally any resampled outputs.
        """
        
        # Check file we want to match exists, throw if it does not
        self._logger.debug(f'Checking if WCS target file {target_wcs_file} exists.')
        self._check_file_exists(target_wcs_file, True)
        
        return
        
    def get_filter_list(self, ap_filetype='input'):
        """
        Returns a list of the unique ``FILTER`` keyword values found the
        specified file type collection.
        
        Parameters
        ----------
        ap_filetype : {'input', 'navfile'}, optional, default='input'
            The AstroPhotography file type for which ``FILTER`` keywords
            will be read and returned.
            This should be one of the stage names described above.
            
        Returns
        -------
        filter_list : list of str or None
            If the type of file have been specified **and** this 'ApProcess'
            instance has read or generated that type then a list of the
            unique ``FILTER`` header keywords is returned. If no inputs
            have been read or no files of that type have been generated, then
            None is returned.
        """
        
        filter_list = None
        ifc_filtertype = None
        keys = ['naxis1', 'naxis2', 'imagetyp', 'object', 'filter', 'exposure']
        
        if ap_filetype in 'input':
            if len(self._input_ifc.summary) == 0:
                err_msg = 'No input files have been specified.'
                self._logger.error(err_msg)
                return filter_list
            ifc_filtertype = self._input_ifc
        elif ap_filetype in 'navfile':
            if len(self._navfile_files) == 0:
                err_msg = f'No {ap_filetype} files have been generated by this instance of ApProcess.'
                self._logger.error(err_msg)
                return filter_list
            ifc_filtertype = ImageFileCollection(self._data_dir, keywords=keys, filenames=self._navfile_files)
        else:
            err_msg = f'Unrecognized file type {ap_filetype} specified.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        
        self._logger.debug(f'get_filter_list: there are {len(ifc_filtertype.summary)} files of type {ap_filetype} found under {self._data_dir}')

        # Get the unique filters, then build a table of filter, num_images, net_exposure_time
        uniq_filter_exposure = unique(ifc_filtertype.summary, keys=['filter'], keep='first', silent=True)
        filter_list = []
        for filter in uniq_filter_exposure['filter']:
            filter_list.append(filter)

        # Generate a valid summary table
        file_list = self._get_processed_file_names(ap_filetype, self._data_dir, name_and_dir=True)
        filter_image_summary = self.make_filter_image_summary(filter_list, file_list)
        table_str_list = filter_image_summary.pformat_all(max_lines=-1, max_width=-1)
        self._logger.debug(f'Unique filter and exposure time data:')
        for line in table_str_list:
            self._logger.debug(f'{line}')

        self._logger.info(f'There are {len(filter_list)} unique filters among {len(ifc_filtertype.summary)} files: {filter_list}')
        return filter_list
        
    def make_filter_image_summary(self, filter_list, img_file_list, resampled_image_info=False):
        """
        Generate a summary table of the total exposure time associated
        each filter in the images in the input ImageFileCollection, 
        and optionally the associated resampled images in those filters.
        
        Notes
        -----
        A user-provided list of files is required.
        There is a significant bug or limitation in ccdproc's ImageFileCollection handling
        of files in subdirectories that prevent it from accurately
        passing the full path, or iterating over files in subdirectories itself.
        So we have to do this the hard way.
        
        Parameters
        ----------
        filter_list : list of str
            List of ``FILTER`` names, e.g. ``['Lum', 'Ha', 'OIII']``.
        img_file_list : list of str
            A list of the file names, **including sub-directory names**,
            of the files for which we want a filter image summary.
            This list conatisn either all the ``input`` calibrated image
            file path and name or the same for the ``navfile`` images.
        resampled_image_info : bool, optional, default=False
            If True then add a column listing any combined resampled images
            (``ap_filetype='stacked'``) for those filters.
            
        Returns
        -------
        info_table : astropy.table.Table
            Informational table of filters, number of navigated images,
            exposure time stats, and optionally any resampled outputs.
        """
        
        # Columns are lists of
        # - Filter name
        # - Num (navigated) images
        # - Total exposure time (minutes)
        # - Min exposure time (minutes)
        # - Max exposure time (minutes)
        # - {optional) Resampled output image
        
        num_images = []
        net_exptime_min = []
        min_exptime_min = []
        max_exptime_min = []
                
        # TODO This loop is inefficient 
        for filter_name in filter_list:
            self._logger.debug(f'Finding files and exposure times for filter {filter_name}')
            sumt = 0
            numf = 0
            mint = None
            maxt = None
            files_in_filter = []
            for a_file in img_file_list:
                # TODO extension number matches input files, not necessarily
                # the navigated fies
                hdr = fits.getheader(a_file, extnum=self._extnum)
                
                # only process file if its FILTER matches expectation
                this_filter_name = hdr['FILTER']
                if this_filter_name == filter_name:
                    self._logger.debug(f'File {a_file} uses filter {this_filter_name}')
                    files_in_filter.append( a_file )
                    
                    exptime_s = util.get_exposure_time(hdr, False)
                    if exptime_s is None:
                        self._logger.warning(f'File {a_file} lacks exposure time information, and will be ignored.')
                        continue
                    else:
                        exptime_m = exptime_s / 60.0
                    sumt = sumt + exptime_m
                    numf = numf + 1
                    if mint is None:
                        mint = exptime_m
                    else:
                        mint = min(mint, exptime_m)
                    if maxt is None:
                        maxt = exptime_m
                    else:
                        maxt = max(maxt, exptime_m)
            num_images.append( numf )
            net_exptime_min.append( sumt )
            min_exptime_min.append( mint )
            max_exptime_min.append( maxt )
            self._logger.info(f'Found the following {len(files_in_filter)} files for filter {filter_name}: {files_in_filter}')        
            
        # Finally contruct the table
        col_names  = ['filter', 'num_navfiles', 'tot_exptime_m', 'min_exptime_m', 'max_exptime_m']
        col_dtypes = [object, np.int32, float, float, float]
        info_table = Table([filter_list, num_images, net_exptime_min, min_exptime_min, max_exptime_min],
            names=col_names,
            dtype=col_dtypes)
        
        if resampled_image_info:
            resampled_imgs  = []
            #temporary
            for idx, filter_name in enumerate(filter_list):
                resampled_imgs.append(f'TBA image {idx}')
            c = Column(data=resampled_imgs, name='stacked_image', dtype=object)
            info_table.add_column(c)
        
        return info_table
        
    def _make_status_table(self, num_input_ifc):
        """
        Create the processing stage status table with one row per input file
        
        Parameters
        ----------
        num_input_ifc : int
            The number of input files
        """

        # define a structure to hold the processing status of the images
        # Note: this is an astropy table, so either string length must be specified,
        # or you must use `object`

        arr_dtype = [('filename', object), ('find_stars_status', object), ('find_stars_time', float), ('astrometry_status', object), ('astrometry_time', float)]

        self._status_table = Table(data=np.empty(num_input_ifc, dtype=arr_dtype))
        self._status_table['find_stars_time'].info.format = '7.3f'
        self._status_table['astrometry_time'].info.format = '7.3f'
        return
