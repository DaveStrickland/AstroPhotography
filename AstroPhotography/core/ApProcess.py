"""
Contains the implementation of the ApProcess class.
"""

# 2024-05-02 dks : Initial implementation.
# 2024-08-12 dks : Numpy format documentation
# 2024-11-01 dks : Ruff/Mypy changes
# 2024-12-31 dks : Add ApFixHoles to preprocessing

from typing import Any
import sys
import os
import logging
from pathlib import Path
import time
from datetime import datetime
import subprocess
import shlex

import numpy as np
import numpy.typing as npt
from ccdproc import ImageFileCollection
from astropy.table import Table, Column
from astropy.table import unique
from astropy.io import fits
from astropy import wcs

# AstroPhotography includes
from .. import __version__
from .ApFixCosmicRays import ApFixCosmicRays as ApFixCosmicRays
import AstroPhotography.util as util
from .ApFindStars import ApFindStars as ApFindStars
from .ApAstrometry import ApAstrometry as ApAstrometry
from .ApQualitySummarizer import ApQualitySummarizer as ApQualitySummarizer
from .ApFixBadPixels import ApFixBadPixels as ApFixBadPixels
from .ApFixHoles import ApFixHoles as ApFixHoles
from .ApFitsRasterizer import ApFitsRasterizer as ApFitsRasterizer

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
    - Reapply one or more forms of pixel corrections (bad pixel correction,
      hole inpainting, and cosmic ray correction), and/or add image header
      metadata in cases where the input calibrated images are deficient.
      See :func:`preprocess_images`. These functions should be run before image
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

    In terms of pseudo-code :func:`navigate_images` performs the following operations:

    .. code-block:: python

        for cal-img in calibrated-images:
            ap_find_stars cal-img --> source-list source-plot fwhm-plot
                quality-report.yml ds9-sources.reg
            ap_astrometry cal-img source-list --> navigated-img
        ap_quality_summary (all quality-report.yml files) --> quality-summary.csv

    Note that processing is sequential and single threaded *at this time*.
    Even if this changes for source detection the astrometry part will
    remain sequential to avoid spamming the ``Astrometry.net`` servers.

    Image Resampling And Mosaicing (Stacking)
    -----------------------------------------

    Astromatic ``swarp`` is used to resample navigated images to a common
    pixel scale and pointing.

    By default all filters found in the navigated images will be
    processed, but a user-specified list of filters can be specified
    instead.

    Currently only a single method of selecting the common scale and pointing
    is provided, although more options will be added in future.

    - Resample and combine all navigated images to match the WCS defined
      in one specified navigated image, separating outputs by FILTER.

    Optional Preprocessing
    ----------------------

    An optional stage of preprocessing is allowed to modify the original
    calibrated input files before performing image navigation
    and astrometry, in cases where the input calibrated files are
    deficient in one or more regards:

    1. Added metadata to the FITS header to identify the observation,
       telescope, or instrument characteristics. This mode is controlled
       by the ``keyword_dict`` and ``replace_keywords`` parameters.
    2. Add an ``EXPOSURE`` keyowrd value based on the value of a less
       commonly used varient already present in the FITS headers
       (e.g. ``ONTIME``). This is controlled by the ``find_exposure_time``
       parameter.
    3. Perform additional bad pixel, bad row, and/or bad column
       correction based on a user-supplied bad pixel file and using
       the :class:`ApFixBadPixels` class
       This mode is controlled
       by the ``badpixelfile`` and ``deltapix`` parameters.
    4. Inpaint "holes" in the image with values matching the local pixel
       statistics based on a user-supplied hole mask and
       :class:`ApFixHoles`. This mode is controlled by the ``holemaskfile``
       parameter.
    5. Apply Cosmic Ray (CR) rejection using :class:`ApFixCosmicRays`.
       If you decide that additional bad pixel
       processing is necessary then it is likely you will also require
       additional CR rejection as well.

    Preprocessing does not modify the original input files. Instead it
    creates new files with names based on string replacement of the of
    original inputs files.

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

    Navigation Output files
    ~~~~~~~~~~~~~~~~~~~~~~~

    Output files from  star detection and astrometry processing correspond
    to one of the following ``ap_filetype`` types:

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

    The allowed AstroPhotography file types associated with image resampling
    and stacking are:

    - ``resampled``: Resampled stacked images corresponding to a single
        filter. These file names depend on additional user-specified choices
        for the file name prefix, suffix, and output directory.

    The allowed AstroPhotography file types associated with preprocessing
    are:

    - ``preproc_out``: The output files generated by running :func:`preprocess_images`
        on the ``input`` files. The file names are based on the ``input`` file
        names, but are distinguished from them by a user-defined textual replacement.

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
       The only ``ap_filestate`` that returns results from
       :func:`get_file_names` for both resampled and preprocessed output file types
       is ``processed``.
    3. ``existing``: Output files that exist on disk that correspond to the
       specified input file parameters. These file need not have been
       processed with this `ApProcess` instance. The :func:`set_file_names`
       function can be used to make the current `ApProcess` instance aware
       of these files and treat them as if it had generated them itself.

    Image Resampling and Stacking Output Files
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The stacked output images do not have a fixed naming convention. Instead
    the file names consist of a user supplied prefix, the filter name, and a
    user supplied suffix.

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

    Warnings
    --------

    - Only a few variants of the input file name specification methods have
      been tested with preprocessing. This may cause problems when using the
      :func:`process_all` function. Calling the processing functions separately
      (e.g. :func:`preprocess_images`, then :func:`navigate_images`, and so on)
      may be more robust if you experience problems with :func:`process_all`.
    """

    def __init__(self, loglevel: str) -> None:
        """
        Initializes an ApProcess instance

        Parameters
        ----------
        loglevel: str
            Standard logging level string, e.g. `INFO`
        """

        self._name: str = "ApProcess"  # str : class name
        self._version: str = __version__  # str : class version
        self._loglevel: str = loglevel  # str : Logging level
        self._initialize_logger(self._loglevel)

        # Processing state variables
        self._nav_status_table: Any = None  # Table : processing status
        self._input_ifc: Any = None  # ImageFileCollection : Input files that were actually used
        self._data_dir: str = None  # Path : to input ``data_dir``
        self._extnum: int | str = 0  # Extension number or name for input data

        # Existing file of the specified types
        self._srclist_files: list[str] = []
        self._regfile_files: list[str] = []
        self._plotfile_files: list[str] = []
        self._qualfile_files: list[str] = []
        self._fwhmplot_files: list[str] = []
        self._navfile_files: list[str] = []
        self._resampled_image_dict: dict[str, str] = {}  # generated by resample_images_to_match
        self._preprocessed_modified_files: list[str] = []
        self._filter_list: list[str] = []

        self.qual_pref: str = "qual"
        self.qual_suff: str = ".yaml"
        self._rasterizer = ApFitsRasterizer(self._loglevel)

        return

    def _check_dir_exists(self, output_dir: str, mkdir: bool = False, throws: bool = True) -> bool:
        """
        Check if a directory path exists, by default raising an exception
        if it does not exist, but optionally creating the directory if mkdir is True

        The functionreturns boolean True if the directory exists or if
        it was successfully created.

        Parameters
        ----------
        output_dir: str
            Name of directory to check the existence of.
        mkdir : bool, optional, default=False
            If True and the directory does not exist, then create the
            directory. In this case this function will throw an exception
            only if the directory creation fails.
        throws : bool, optional, default=True
            If throws is True then an exception is thrown if the directory
            path does not exist and ``mkdir=False``. If throw is False
            and ``mkdir=False`` then the function will return False to
            the caller if the directory path does not exist.


        Returns
        -------
        exists : bool
            Boolean flag that is True if the directory exists or was
            created, and False if it does not exist and throw is False.

        Raises
        ------
        RuntimeError
            If the named file is not found and throw is True
        """

        exists = os.access(output_dir, os.F_OK)
        if exists:
            self._logger.debug(f"Subdirectory {output_dir} already exists.")
        else:
            # Does not exist
            if mkdir:
                try:
                    os.mkdir(output_dir)
                except:
                    self._logger.error(f"Error, mkdir failed trying to create {output_dir}")
                    raise
                self._logger.info(f"Successfully created subdirectory {output_dir}")
            else:
                # Does not exist and don't try to create it
                if throws:
                    err_msg = f"Cannot find directory {output_dir}"
                    self._logger.error(err_msg)
                    raise RuntimeError(err_msg)
        return exists

    def _check_file_exists(self, filename: str, throws=True) -> bool:
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
            err_msg = f"Cannot find {filename}. Not a valid path or file."
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return exists

    def _img_stats(
        self, data: npt.NDArray, label: str, verbose=False
    ) -> tuple[float, float, float, float]:
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

        minval = np.nanmin(data)
        maxval = np.nanmax(data)
        meanval = np.nanmean(data)

        # percentiles, 50th percentile is the median
        #         0    1    2    3   4   5   6   7   8   9   10
        ipctls = [0.1, 1.0, 5.0, 10, 25, 50, 75, 90, 95, 99, 99.9]
        opctls = np.nanpercentile(data, ipctls)
        medval = opctls[5]

        if verbose:
            msg1: str = (
                f"{label} data min={minval:.2f}, max={maxval:.2f},"
                f" mean={meanval:.2f}, median={medval:.2f} ADU."
            )
            msg2: str = (
                f"  90% of data between {opctls[2]:.2f} and {opctls[8]:.2f} ADU (5-95 precentiles)"
            )
            msg3: str = (
                f"  98% of data between {opctls[1]:.2f} and {opctls[9]:.2f} ADU (1-99 precentiles)"
            )
            self._logger.info(msg1)
            self._logger.info(msg2)
            self._logger.info(msg3)
        return [minval, maxval, meanval, medval]

    def _initialize_logger(self, loglevel: str) -> None:
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
            raise ValueError("Invalid log level: {}".format(loglevel))
        self._logger.setLevel(numeric_level)
        self._logger.propagate = False

        # check if handlers already present
        if not len(self._logger.handlers):
            # create console handler and set level to debug
            ch = logging.StreamHandler()
            ch.setLevel(numeric_level)

            # create formatter
            formatter = logging.Formatter("%(asctime)s | %(name)s | %(levelname)s | %(message)s")

            # add formatter to ch
            ch.setFormatter(formatter)

            # add ch to logger
            self._logger.addHandler(ch)
        return

    def _read_fits(
        self, image_filename: str, image_extension: int | str
    ) -> tuple[npt.NDArray, Any]:
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
        self._logger.info(
            "Loading extension {} of FITS file {}".format(image_extension, image_filename)
        )

        # open() parameters that can be important.
        # Default values used here.
        # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
        uint_handling = True
        image_scaling = False

        with fits.open(
            image_filename, uint=uint_handling, do_not_scale_image_data=image_scaling
        ) as hdu_list:
            ext_hdr = hdu_list[image_extension].header
            ext_data = hdu_list[image_extension].data

        ndim = ext_hdr["NAXIS"]
        cols = ext_hdr["NAXIS1"]
        rows = ext_hdr["NAXIS2"]
        bitpix = ext_hdr["BITPIX"]
        info_str = "{}-D BITPIX={} image with {} columns, {} rows".format(ndim, bitpix, cols, rows)

        if ndim == 3:
            layers = ext_hdr["NAXIS3"]
            info_str = "{}-D BITPIX={} image with {} columns, {} rows, {} layers".format(
                ndim, bitpix, cols, rows, layers
            )

        if "BSCALE" in ext_hdr:
            bscale = ext_hdr["BSCALE"]
            info_str += f", BSCALE={bscale}"

        if "BZERO" in ext_hdr:
            bzero = ext_hdr["BZERO"]
            info_str += f", BZERO={bzero}"

        self._logger.debug(info_str)
        if ndim == 3:
            self._logger.error("Error, 3-D handling has not been implemented yet.")
            sys.exit(1)

        # Convert to 32-bit floating point if necessary
        if not np.issubdtype(ext_data.dtype, np.floating):
            orig_dtype = ext_data.dtype
            ext_data = ext_data.astype(np.float32)
            self._logger.debug(f"  Converted data type from {orig_dtype} to float32")

        # Get data absolute limits.
        minval = np.nanmin(ext_data)
        maxval = np.nanmax(ext_data)
        medval = np.nanmedian(ext_data)
        self._logger.debug(
            f"Raw data statistics are min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}"
        )

        # Is there a PEDESTAL value? MaximDL likes to add an offset, and
        # the PEDESTAL value is the value to ADD to the data to remove the
        # pedestal.
        if "PEDESTAL" in ext_hdr:
            pedestal = float(ext_hdr["PEDESTAL"])
            if pedestal != 0:
                self._logger.debug(f"Removing a PEDESTAL value of {pedestal} ADU.")
                ext_data += pedestal
                minval = np.amin(ext_data)
                maxval = np.amax(ext_data)
                medval = np.median(ext_data)
                self._logger.debug(
                    f"After PEDESTAL removal, min={minval:.2f}"
                    f", max={maxval:.2f}, median={medval:.2f}"
                )

        return ext_data, ext_hdr

    def _remove_pedestal_kw(self, hdr: Any) -> None:
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

        if "PEDESTAL" in hdr:
            self._logger.debug("Removing PEDESTAL keyword from FITS header.")
            del hdr["PEDESTAL"]

        return

    def _write_corrected_image(
        self,
        inpdata_file: str,
        ext_num: int | str,
        outdata_file: str,
        odata: npt.NDArray,
        odict: dict[str, Any],
    ) -> None:
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

        self._logger.debug(f"FITS header keywords added to output: {odict}")
        self._check_file_exists(inpdata_file)

        # open() parameters that can be important.
        # Default values used here.
        # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
        uint_handling = True
        image_scaling = False

        with fits.open(
            inpdata_file, uint=uint_handling, do_not_scale_image_data=image_scaling
        ) as hdu_list:
            self._remove_pedestal_kw(hdu_list[ext_num].header)
            for kw in ["BSCALE", "BZERO"]:
                if kw in hdu_list[ext_num].header:
                    del hdu_list[ext_num].header[kw]

            # Modify data
            hdu_list[ext_num].data = odata

            # Modify header
            for kw, val in odict.items():
                hdu_list[ext_num].header[kw] = val

            tnow = datetime.now().isoformat(timespec="milliseconds")
            hdu_list[ext_num].header["HISTORY"] = (
                f"Processed by {self._name} {self._version} at {tnow}"
            )

            # Write to new file
            hdu_list.writeto(outdata_file, output_verify="ignore", overwrite=True)

        self._logger.info(f"Wrote bias/dark/flat corrected file to {outdata_file}")
        return

    def get_file_names(
        self,
        ap_filetype: str,
        data_dir: str,
        ap_filestate: str = "conceptual",
        include_pattern: str | None = None,
        exclude_pattern: str | None = None,
        input_file_list: list[str] | None = None,
        input_rootname: str | None = None,
        input_suffix: str | None = ".fits",
        name_and_dir: bool = False,
    ) -> list[str]:
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

        The allowed AstroPhotography file types associated with image
        navigation are:

        - ``input``: These are input files that will have star detection performed on them.
        - ``srclist``: Star detection output FITS table files for each input file.
        - ``regfile``: ds9-format region files.
        - ``plotfile``: PNG plots of the input images with the detected star-like sources
          plotted as circles.
        - ``qualfile``: YaML source detection quality summary files.
        - ``fwhmplot``: PNG plots zooming in around a subset of the detected
          sources in each input image.
        - ``navfile``: Copies of the input files but with valid WCS solutions.

        The allowed AstroPhotography file types associated with image resampling
        and stacking are:

        - ``resampled``: Resampled stacked images corresponding to a single
          filter.

        The allowed AstroPhotography file types associated with preprocessing
        are:

        - ``preproc_out``: The output files generated by running :func:`preprocess_images`
          on the ``input`` files.

        Note
        ----

        For ``resampled`` files the names cannot easily be determined
        ahead of time because (a) the user controls the resampled file name prefix
        and suffix, and (b) the set of filters is not known ahead of performing
        image navigation. Similarly the modified input files generated by
        preprocessing have additional user-defined prefixes. To avoid over complicating this
        function the only ``ap_filestate`` that returns results from
        :func:`get_file_names` for both resampled and preprocessed output file types
        is ``processed``.

        Parameters
        ----------
        ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile', 'resampled', 'preproc_out'}
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

        Raises
        ------
        RuntimeException :
            If the ap_filetype is not one of the allowed stages defined above

        See Also
        --------
        :func:`get_directories` :
            Returns the directory path in which files for given processing
            stage are, or will be.
        `File Types`_ :
            File type and file state documentation
        """  # noqa: E501

        # Returning the files that this instance has already processed?
        must_exist = False
        if "processed" in ap_filestate:
            output_fname_list = self._get_processed_file_names(ap_filetype, data_dir, name_and_dir)
            return output_fname_list
        else:
            nonnav_stages = [
                "resampled",
                "preproc_out",
            ]
            if ap_filetype in nonnav_stages:
                # We only return non-empty data for ap_filestate=='processed'
                # for these file types, so...
                self._logger.warn(
                    (
                        "get_file_names only returns file names for "
                        f'filetype={ap_filetype} for processing stage "processed"'
                    )
                )

            if "existing" in ap_filestate:
                must_exist = True

        self._logger.debug(f"Current working directory: {os.getcwd()}")
        self._logger.debug(f"Attempting to find stars in FITS files within directory={data_dir}")
        if input_file_list is None:
            # Use patterns
            self._logger.debug(f'Using files that match include_pattern="{include_pattern}"')
            self._logger.debug(f'Excluding files that match exclude_pattern="{exclude_pattern}"')
            ifc_cal = ImageFileCollection(
                data_dir, glob_include=include_pattern, glob_exclude=exclude_pattern
            )
        else:
            # Use explicit file list
            self._logger.info(f"Using specified file list: {input_file_list}")
            ifc_cal = ImageFileCollection(data_dir, filenames=input_file_list)

        if not name_and_dir:
            odir = self.get_directories(ap_filetype, data_dir, False)
            self._logger.debug(f"Note that output files will be in subdirectory {odir}")

        fname_list = self._get_file_names(
            ap_filetype, ifc_cal, input_rootname, input_suffix, name_and_dir
        )
        if must_exist:
            # Check file exists before adding it to output file name list
            output_fname_list = []
            for fname_val in fname_list:
                if not name_and_dir:
                    check_fname_val = odir + fname_val  # The actual file path
                else:
                    check_fname_val = fname_val  # Already include the directory
                if self._check_file_exists(check_fname_val, False):
                    output_fname_list.append(fname_val)
        else:
            output_fname_list = fname_list
        return output_fname_list

    def _get_processed_file_names(
        self, ap_filetype: str, data_dir: str, name_and_dir: bool = False
    ) -> list[str]:
        """
        Return the files of the specified file type that this instance is
        aware of, either by processing them itself or by calls to set_file_names.

        Parameters
        ----------
        ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile', 'resampled', 'preproc_out'}
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
        """  # noqa: E501

        allowed_stages = [
            "input",
            "srclist",
            "regfile",
            "plotfile",
            "qualfile",
            "fwhmplot",
            "navfile",
            "resampled",
            "preproc_out",
        ]
        if ap_filetype not in allowed_stages:
            err_msg = (
                f"Requested ap_filetype {ap_filetype} is not"
                f" one of the allowed values: {allowed_stages}"
            )
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)

        # Get file names as a lists from the member variables
        if ap_filetype in "input":
            ftype_list = []
            for file in self._input_ifc.summary["file"]:
                ftype_list.append(file)
        elif ap_filetype in "srclist":
            ftype_list = self._srclist_files
        elif ap_filetype in "regfile":
            ftype_list = self._regfile_files
        elif ap_filetype in "plotfile":
            ftype_list = self._plotfile_files
        elif ap_filetype in "fwhmplot":
            ftype_list = self._fwhmplot_files
        elif ap_filetype in "qualfile":
            ftype_list = self._qualfile_files
        elif ap_filetype in "navfile":
            ftype_list = self._navfile_files
        elif ap_filetype in "resampled":
            ftype_list = [
                self._resampled_image_dict[key] for key in self._resampled_image_dict.keys()
            ]
        elif ap_filetype in "preproc_out":
            ftype_list = self._preprocessed_modified_files
        else:
            self._logger.error(f"Unexpected ap_filetype of {ap_filetype} specified.")

        # The member variables store the file name and directory
        # with respect to data_dir, so if name_and_dir is False then
        # the directory name must be removed.
        if not name_and_dir:
            output_fname_list = []
            odir = self.get_directories(ap_filetype, data_dir, False)
            dir_str_len = len(odir)
            for file in ftype_list:
                dir_start_pos = file.find(odir)
                output_fname_list.append(file[dir_start_pos + dir_str_len - 1 :])
        else:
            output_fname_list = ftype_list

        return output_fname_list

    def _get_file_names(
        self,
        ap_filetype: str,
        ifc_cal: Any,
        input_rootname: str | None = None,
        input_suffix: str | None = ".fits",
        name_and_dir: bool = False,
    ) -> list[str]:
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
        ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile', 'resampled', 'preproc_out'}
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
        """  # noqa: E501

        allowed_stages = [
            "input",
            "srclist",
            "regfile",
            "plotfile",
            "qualfile",
            "fwhmplot",
            "navfile",
            "resampled",
            "preproc_out",
        ]
        if ap_filetype not in allowed_stages:
            err_msg = (
                f"Requested ap_filetype {ap_filetype} is not"
                f" one of the allowed values: {allowed_stages}"
            )
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)

        outputs_fname_list = []
        for file in ifc_cal.summary["file"]:
            ofile, output_dir = util.namefn_calibrated_input(
                file, input_rootname, input_suffix, ap_filetype
            )
            if name_and_dir:
                outputs_fname_list.append(output_dir + ofile)
            else:
                outputs_fname_list.append(ofile)

        return outputs_fname_list

    def get_directories(
        self, ap_filetype: str, data_dir: str, absolute: bool = False, resampled_dir: str = r"./"
    ) -> str:
        """
        Returns the path to the directory that will or does hold files
        corresponding to one of the AstroPhotography file types

        The directory may or may not exist already.

        The allowed file type associated with image navigation are:

        - ``input``: These are input files that will have star detection performed on them.
        - ``srclist``: Star detection output FITS table files for each input file.
        - ``regfile``: ds9-format region files.
        - ``plotfile``: PNG plots of the input images with the detected star-like sources
          plotted as circles.
        - ``qualfile``: YaML source detection quality summary files.
        - ``fwhmplot``: PNG plots zooming in around a subset of the detected
          sources in each input image.
        - ``navfile``: Copies of the input files but with valid WCS solutions.

        The allowed AstroPhotography file types associated with image resampling
        and stacking are:

        - ``resampled``: Resampled stacked images corresponding to a single
          filter.

        The allowed AstroPhotography file types associated with preprocessing
        are:

        - ``preproc_out``: The output files generated by running :func:`preprocess_images`
          on the ``input`` files.

        Parameters
        ----------
        ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile', 'resampled', 'preproc_out'}
            The AstroPhotography file type for which names should be returned.
            This should be one of the stage names described above.
        data_dir : str or path
            Path to base directory containing FITS files to process,
            e.g. `./`.
        absolute : bool, optional, default=False
            If True then the absolute directory path is returned, instead
            of the path relative to ``data_dir``.
        resampled_dir : str, optional, default='.'
            Directory in which resampled output images will be written.
            If None then the current directory is used. If not None then
            the named directory will be created if it does not exist.

        Returns
        -------
        output_dir : str
            Directory name in which output_files will be written, relative
            to the root directory established by the input files.
        """  # noqa: E501

        allowed_stages = [
            "input",
            "srclist",
            "regfile",
            "plotfile",
            "qualfile",
            "fwhmplot",
            "navfile",
            "resampled",
            "preproc_out",
        ]
        nonnav_stages = [
            "resampled",
            "preproc_out",
        ]

        if ap_filetype not in allowed_stages:
            err_msg = (
                f"Requested ap_filetype {ap_filetype} is not"
                f" one of the allowed values: {allowed_stages}"
            )
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)

        # Relative path
        if ap_filetype in nonnav_stages:
            if ap_filetype == "resampled":
                output_dir = resampled_dir
            else:
                output_dir = util.namefn_getdir("input")
        else:
            output_dir = util.namefn_getdir(ap_filetype)
        self._logger.debug(
            f"For ap_filetype {ap_filetype} the relative directory path is {output_dir}"
        )
        if absolute:
            p = Path(data_dir) / Path(output_dir)
            output_dir = str(p.absolute())
            self._logger.debug(
                f"For ap_filetype {ap_filetype} the absolute directory path is {output_dir}"
            )
        return output_dir

    def _generate_all_output_names(
        self,
        inputfile: str,
        input_rootname: str | None = None,
        input_suffix: str = ".fits",
        name_and_dir: bool = False,
        mkdir: bool = False,
    ) -> tuple[str, str, str, str, str, str]:
        """
        Given an input file name and a conversion dictionary,
        generate all the file names that might
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

        oname_type = ["srclist", "regfile", "plotfile", "fwhmplot", "qualfile", "navfile"]
        ofile_dict = {}
        for ap_filetype in oname_type:
            ofile, output_dir = util.namefn_calibrated_input(
                inputfile, input_rootname, input_suffix, ap_filetype
            )

            dir_exists = self._check_dir_exists(output_dir, mkdir, False)
            if not dir_exists:
                # Directory doesn't exist
                self._logger.error(f"Directory {output_dir} does not exist and mkdir={mkdir}")

            if name_and_dir:
                ofile = output_dir + ofile
            ofile_dict[ap_filetype] = ofile

        srclist = ofile_dict["srclist"]
        regfile = ofile_dict["regfile"]
        plotfile = ofile_dict["plotfile"]
        fwhmplot = ofile_dict["fwhmplot"]
        qualfile = ofile_dict["qualfile"]
        navfile = ofile_dict["navfile"]

        self._logger.debug(
            (
                f"Output file names for input {inputfile},"
                f" input_rootname={input_rootname}, "
                f" input_suffix={input_suffix} are"
                f" srclist={srclist}, regfile={regfile}"
                f" plotfile={plotfile}, fwhmplot={fwhmplot}"
                f" qualfile={qualfile} and navfile={navfile}"
            )
        )

        return (srclist, regfile, plotfile, fwhmplot, qualfile, navfile)

    def set_file_names(
        self,
        data_dir: str,
        include_pattern: str | None = None,
        exclude_pattern: str | None = None,
        input_file_list: list[str] | None = None,
        input_rootname: str | None = None,
        input_suffix: str | None = ".fits",
    ) -> None:
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
        keys = ["naxis1", "naxis2", "imagetyp", "object", "filter", "exposure"]

        self._logger.info(f"Finding input and output files starting in directory={data_dir}")
        self._logger.debug(f"Current working directory: {os.getcwd()}")
        if input_file_list is None:
            # Use patterns
            self._logger.info(f'Using files that match include_pattern="{include_pattern}"')
            self._logger.info(f'Excluding files that match exclude_pattern="{exclude_pattern}"')
            ifc_cal = ImageFileCollection(
                data_dir, keywords=keys, glob_include=include_pattern, glob_exclude=exclude_pattern
            )
        else:
            # Use explicit file list
            self._logger.info(f"Using specified file list: {input_file_list}")
            ifc_cal = ImageFileCollection(data_dir, keywords=keys, filenames=input_file_list)
        self._input_ifc = ifc_cal
        self._data_dir = Path(data_dir)

        # Generate status table, written to self._nav_status_table
        num_inputs = len(ifc_cal.summary)
        self._make_nav_status_table(num_inputs)

        self._logger.info(f"There are {num_inputs} input files matching the parameters given.")
        self._logger.debug(f"Input file collection:\n{ifc_cal}")

        # Now look for pre-existing files matching the expected names,
        # and add them to the file lists. Clear the existing lists first...
        self._logger.info("Looking for existing files matching name conventions...")
        name_and_dir = True
        oname_type = ["srclist", "regfile", "plotfile", "fwhmplot", "qualfile", "navfile"]
        for ftype in oname_type:
            # Clear old values
            if ftype in "srclist":
                ftype_list = self._srclist_files
            if ftype in "regfile":
                ftype_list = self._regfile_files
            if ftype in "plotfile":
                ftype_list = self._plotfile_files
            if ftype in "fwhmplot":
                ftype_list = self._fwhmplot_files
            if ftype in "qualfile":
                ftype_list = self._qualfile_files
            if ftype in "navfile":
                ftype_list = self._navfile_files
            ftype_list.clear()

            self._logger.debug(f"Looking for existing {ftype} files.")
            flist = self._get_file_names(
                ftype, self._input_ifc, input_rootname, input_suffix, name_and_dir
            )
            for fname_val in flist:
                if self._check_file_exists(fname_val, False):
                    if ftype in "srclist":
                        self._srclist_files.append(fname_val)
                        continue
                    if ftype in "regfile":
                        self._regfile_files.append(fname_val)
                        continue
                    if ftype in "plotfile":
                        self._plotfile_files.append(fname_val)
                        continue
                    if ftype in "fwhmplot":
                        self._fwhmplot_files.append(fname_val)
                        continue
                    if ftype in "qualfile":
                        self._qualfile_files.append(fname_val)
                        continue
                    if ftype in "navfile":
                        self._navfile_files.append(fname_val)
                        continue
            numfiles = len(ftype_list)
            self._logger.info(f"Found the following {numfiles} {ftype} files: {ftype_list}")
        return

    def process_all(
        self,
        data_dir: str | None,
        include_pattern: str | None = None,
        exclude_pattern: str | None = None,
        input_file_list: list[str] | None = None,
        input_rootname: str | None = None,
        input_suffix: str | None = ".fits",
        final_quality_file: str | None = None,
        clean_star_detection: bool = False,
        clean_astrometry: bool = False,
        stop_on_error: bool = False,
        extnum: str | int = 0,
        search_fwhm: float = 3.0,
        search_nsigma: float = 7.0,
        detector_bitdepth: int = 16,
        max_sources: int = 200,
        nosatmask: bool = True,
        sat_frac: float = 0.8,
        quiet: bool = True,
        exclude_corner_pct: float | None = None,
        srclist_extname: str = "AP_XYPOS",
        astnet_key=None,
        use_sip: bool = False,
        user_scale: float | None = None,
        scale_err_ratio: float | None = None,
        target_wcs_file: str | None = None,
        resampled_file_prefix: str = "resampled_",
        resampled_file_suffix: str = "_resamp_weighted.fits",
        resampled_dir: str = r"./",
        filter_list: list[str] | None = None,
        resampled_summary_plot: str | None = None,
        composite_summary_plot: str | None = None,
        combine_type: str = "WEIGHTED",
        preprocess_replace: str | None = None,
        preprocess_with: str | None = None,
        find_exposure_time: bool = False,
        keyword_dict: dict[str, Any] | None = None,
        replace_keywords: bool = False,
        badpixelfile: str | None = None,
        deltapix: int = 2,
        holemaskfile: str | None = None,
        fix_cosmic_rays: bool = False,
        clean_preprocess: bool = True,
    ):
        """
        Performs all processing stages

        This function wraps :func:`navigate_images`, :func:`create_quality_summary`
        :func:`resample_images_to_match`, and optionally :func:`preprocess_images`.

        Parameters are listed in the same order as the functions listed above.
        Note that the preprocessing options are the final set of function
        parameters, but if activated preprocessing occurs before image
        navigation.

        Note
        ----

        **Image resampling:** The following parameters **must** specified
        in order to activate image resampling and stacking: ``target_wcs_file``,
        ``resampled_file_prefix``, ``resampled_file_suffix``, and
        ``resampled_dir``. The ``filter_list`` parameter is optional.

        **Image preprocessing:** The following parameters **must** specified in
        order to activate preprocessing: ``preprocess_replace`` and ``preprocess_with``
        along with at least one of the following parameters **not** being None:
        ``find_exposure_time``, ``keyword_dict``,
        ``badpixelfile``, ``holemaskfile``, or ``fix_cosmic_rays``.

        Parameters
        ----------
        data_dir : str or path
            Path to base directory containing FITS files to process,
            e.g. `./`. (navigate_image parameter, preprocess_images parameter)
        include_pattern : str, optional
            Globbing pattern for input files we want included, specified
            relative to data_dir. If not specified all FITS files
            will be included. This parameter is ignored if file_list
            is not None. (navigate_image parameter, preprocess_images parameter)
        exclude_pattern : str, optional
            Globbing pattern of files we want excluded, specified
            relative to data_dir. If not specified no FITS files
            will be excluded. This parameter is ignored if file_list
            is not None. (navigate_image parameter, preprocess_images parameter)
        input_file_list : list of str, optional
            An explicit list of files, paths relative to data_dir,
            may be specified. If provided only those files in input_file_list
            are looked for, and the include and exclude patterns are
            ignored. (navigate_image parameter, preprocess_images parameter)
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
            (navigate_image parameter)
        input_suffix : str, optional, default='.fits'
            String denoting the file type suffix of the input image file.
            For example, '.fits' or '.fits' or '.ftz' or 'fits.gz' or '.fits.bz2'
            (navigate_image parameter, preprocess_images parameter)
        final_quality_file : str, optional, default=None
            If specified, this is the name for a quality summary file
            in CSV format
            generated using :class:`ApQualitySummarizer`, which summarizes
            all YaML quality files found in the standard quality
            summary file directory (typically `./MetaData/`).
            (navigate_image parameter)
        clean_star_detection : bool, default=False
            If True then **existing** star detection outputs (``srclist``,
            ``regile``, ``plotfile``, ``qualfile``, and ``fwhmplot`` files)
            will be deleted and regenerated.
            (navigate_image parameter)
        clean_astrometry : bool, default=False
            If True then **existing** astrometry outputs (``navfile``files)
            will be deleted and regenerated.
            (navigate_image parameter)
        stop_on_error : bool, optional, default=False
            By default :func:`navigate_images` will continue attempting to
            process all the input files even if a given stage fails on one
            of the inputs. That failure will be recorded in the output
            ``nav_status_table``. If you want it to instead stop processing
            when a failure is detected, set stop_on_error to True.
            (navigate_image parameter)
        extnum : int or str, optional, default=0
            Extension number or name for the extension holding the image data.
            Usually this is 0, for the ``PrimaryHDU``.
            (navigate_image parameter, preprocess_images parameter
        search_fwhm : float, optional, default=3.0
            Initial guess or estimate of the stellar PSF FWHM in pixels
            in this image (navigate_image parameter)
        search_nsigma : float, optional, default=7.0
            Minimumn number of sigma above background for a detection
            (navigate_image parameter)
        detector_bitdepth : int, optional, default=16
            Detector bit-depth, used in estimating which pixels are saturated.
            (16 for most CCDs, ?? for CMOS, ?? for camera)  (navigate_image parameter)
        max_sources : int, optional, default=200
            Maximum number of sources to output or None. If not None
            then only the max_sources brightest sources will be output.
            Typically there is no advantage to having very large numbers
            of detected sources when performing astrometry, and may
            well slow it down. I normally use max_sources=200.
            (navigate_image parameter)
        nosatmask : bool, optional, default=True
            If True, keep possibly saturated stars. (navigate_image parameter)
        sat_frac : float, optional, default=0.8
            Fraction of full well at which we assume star saturated
            (navigate_image parameter)
        quiet : bool, optional, default=True
            If True this suppresses the runtime source list printing
            to STDOUT (navigate_image parameter)
        exclude_corner_pct : float or None, optional, default=None
            If a positive value in the range [0,100] is specified then semi-circular
            regions of radius (exclude_corner_pct/100)*max(img_rows, img_cols)
            around each corner will be excluded from source detection. This
            parameter should be used if the image suffers from poor flat-fielding
            or vignetting that creates spurious sources at the image corners.
            (navigate_image parameter)
        srclist_extname : str, optional, default='AP_XYPOS'
            FITS extension name for star X,Y position data. Default=``AP_XYPOS``
            (navigate_image parameter)
        astnet_key : str, optional, default=None
            Your personal Astrometry.net API key, if you have not already
            added it to your ~/.astropy/config/astroquery.cfg config file.
            (navigate_image parameter)
        use_sip : bool, optional, default=False
            Allow astrometry.net to fit SIP polynomial distortion terms.
            This may be necessary for very large fields of view (>10 deg),
            but SIP is not treated correctly by swarp (and possibly other
            software).
            (navigate_image parameter)
        user_scale : float, optional, default=None
            If specified, override the estimated plate  scale in the source list file
            and instead use a user spacified estimate of the plate scale.
            The units are arcseconds/pixel.
            (navigate_image parameter)
        scale_err_ratio: float, optional, default=None
            The relative uncertainty in the estimated plate scale,
            expressed as a ratio. This applies to either the default
            estimate from the source list, or a user-supplied plate
            scale. For example, if the estimated plate  scale is 2.0 arcsec/pix
            and the scale_err_ratio=1.5, then the plate scale range that
            will be search by Astrometry.net is 2/1.5 (=4/3) to 2*1.5 (=3)
            arcseconds. If not specified ApAstrometry will use a value of 1.3.
            Using a larger value can help in cases where astrometric
            solutions fail, for example if incorrect telescope metadata
            leads to inaccurate estimated plate scales.
            (navigate_image parameter)
        target_wcs_file : str, optional, default=None
            File name, including path, to the FITS file than contains
            the WCS that all navigated images should be resampled to
            match.
            If set to ``None`` then the first navigated file is used.
            (resample_images_to_match parameter)
        resampled_file_prefix : str. optional, default='resampled_'
            The part of the resampled image file name in front of the
            filter name. For example, for a ``Green`` filter output
            image name of ``M101_Supernova_2023ixf_Green_resamp_weighted.fits``
            then ``resampled_file_prefix='M101_Supernova_2023ixf_'``.
            (resample_images_to_match parameter)
        resampled_file_suffix : str, optional, default='_resamp_weighted.fits'
            The part of the resampled image file name after the filter
            name. **Note** that this should include the file suffix, which
            should include ``.fits`` (or further, e.g. ``.fits.gz``)
            For example, for a ``Green`` filter output
            image name of ``M101_Supernova_2023ixf_Green_resamp_weighted.fits``
            then ``resampled_file_suffix='_resamp_weighted.fits'``.
        resampled_dir : str, optional, default='.'
            Directory in which resampled output images will be written.
            If None then the current directory is used. If not None then
            the named directory will be created if it does not exist.
            (resample_images_to_match parameter)
        filter_list : list of str or None, optional, default=None
            If specified then only images with ``FILTER`` header keywords
            matching one of the input list strings will be resampled. For
            example, to only resample the H-alpha and luminance images
            the filter list would be ``['Lum', 'Ha']``
            (resample_images_to_match parameter)
        resampled_summary_plot: str or None, optional, default=None
            If not None then a PNG file with the specified name containing plots
            of each output resampled image will be generated. The plots are
            generated using :func:`util.load_three_images_and_plot`.
            Default intensity scaling and WCS plotting options will be used.
            (resample_images_to_match parameter)
        composite_summary_plot: str or None, optional, default=None
            If there are resampled images in all three Red, Green, and Blue
            filters and/or all three Ha, SII, and OIII filters and this
            argument is not None then a PNG file with summary plots of three
            color composities will be generated. If Red, Green, and Blue
            filters are present then an ``RGB`` three color composite plot
            will be generated. If the SII, H-alpha, and OIII filters
            are present then an ``SHO`` color-scheme composite plot
            will be generated. The plots are
            generated using :func:`util.plot_lupton_threecolor`.
            Default intensity scaling and WCS plotting options will be used.
            (resample_images_to_match parameter)
        combine_type : str, optional, default='WEIGHTED'
            One of the valid swarp ``COMBINE_TYPE`` values, i.e. one
            of ``AVERAGE, CHI2, MEDIAN, MIN, MAX, SUM, WEIGHTED``.
            I recommened including some information on the chosen type in your
            ``resampled_file_suffix``, e.g. in the example above the suffix was
            ``resampled_file_suffix='_resamp_weighted.fits'``.
            (resample_images_to_match parameter)
        preprocess_replace : str
            Substring common to all the input calibrated files that will
            be replaced with ``preprocess_with``. (preprocess_images parameter)
        preprocess_with : str
            String that will replace ``preprocess_replace`` in all
            input file names.
            (preprocess_images parameter)
        find_exposure_time : bool, optional, default=False
            If True then :func:`ApUtil.get_exposure_time` will be used to
            extract the exposure time from the input file headers and
            set the ``EXPOSURE`` keyword if it is not already present.
            (preprocess_images parameter)
        keyword_dict : dict, optional, default=None
            If not None, then supply a dictionary of
            ``keyword: (value, comment)`` entries that
            will be added to the processed file FITS headers.
            (preprocess_images parameter)
        replace_keywords : bool, optional, default=False
            If True then existing FITS header keys that match the keys
            in ``keyword_dict`` will be over-written with the new values.
            (preprocess_images parameter)
        badpixelfile : str, optional, default=None
            File path and name to the bad pixel file to apply. This file
            should conform to the format generated by ``ApFindBadPixels``
            and used by ``ApFixBadPixels``.
            (preprocess_images parameter)
        deltapix : int, optional, default=2
            Linear distance away from a bad pixel from which
            the median value of the good pixels will be drawn. If 1 then
            the median value of good pixels within the surrounding 8 pixels
            will be used. If 2 then the median of the good pixels within the
            surrounding 24 pixels will be used. Values above 2 are not
            recommended.
            (preprocess_images parameter)
        holemaskfile : str, optional, default=None
            File path and name to any hole mask file to apply. This file
            should conform to the format used by ``ApFixHoles``.
            (preprocess_images parameter)
        fix_cosmic_rays : bool, optional, default=False
            If True then perform Cosmic Ray rejection on the images.
            (preprocess_images parameter)
        clean_preprocess : bool, optional, default=True
            Overwrite any existing file with the same name as a preprocessed
            output file. If False, then the presence of the output file
            will cause preprocessing to be skipped for the assciated input
            file.
            (preprocess_images parameter)

        Returns
        -------
        nav_status_table : astropy.table
            An astropy table containing the processing status of the images
        resampled_images : dict of str, str
            A dictionary consiting of filter name (key), file path
            and name of the resample images (value) pairs, as
            generated by this function
        resampled_info_table : astropy.table.Table
            Informational table of filters, number of navigated images,
            exposure time stats, and optionally any resampled outputs.
        preprocessed_modified_files : list of str
            List of file names of the output modified files
        """

        procall_tstart = time.perf_counter()
        preprocessed_modified_files = None
        nav_status_table = None
        resampled_images = None
        resampled_info_table = None

        # TODO Add exception handling.

        # Check if we are going to do preprocessing...
        if (preprocess_replace is not None) and (preprocess_with is not None):
            # Check that at least one pre-processing option has been
            # selected, otherwise this is an error
            if (
                find_exposure_time
                or keyword_dict
                or badpixelfile
                or holemaskfile
                or fix_cosmic_rays
            ):
                preprocessed_modified_files = self.preprocess_images(
                    data_dir,
                    preprocess_replace,
                    preprocess_with,
                    include_pattern,
                    exclude_pattern,
                    input_file_list,
                    input_rootname,
                    input_suffix,
                    extnum,
                    find_exposure_time,
                    keyword_dict,
                    replace_keywords,
                    badpixelfile=badpixelfile,
                    deltapix=deltapix,
                    holemaskfile=holemaskfile,
                    fix_cosmic_rays=fix_cosmic_rays,
                    clean_preprocess=clean_preprocess,
                )
                self._logger.debug(
                    (
                        "Preprocessing has produced the following files: "
                        f"{preprocessed_modified_files}"
                    )
                )
            else:
                err_msg: str = (
                    "Preprocessing arguements preprocess_replace"
                    " and preprocess_with specified, but no preprocessing"
                    " options actually invoked. One or more of "
                    "find_exposure_time, keyword_dict, badpixelfile,"
                    " fix_cosmic_rays must be non-NULL."
                )
                self._logger.error(err_msg)
                raise RuntimeError(err_msg)

        # if we've preprocessed any files then we need to use a modified
        # call to navigate images, otherwise we use the default call.
        actual_incl_pattern: str | None = include_pattern
        actual_excl_pattern: str | None = exclude_pattern
        actual_input_rootname: str | None = input_rootname
        actual_input_file_list: list[str] | None = input_file_list
        if preprocessed_modified_files is not None:
            # Use the modified file list, and explicitly set input_rootname
            # Disable include and exclude patterns
            actual_input_file_list = preprocessed_modified_files
            actual_incl_pattern = None
            actual_excl_pattern = None
            if preprocessed_modified_files[0].startswith(preprocess_with):
                actual_input_rootname = preprocess_with

            msg: str = (
                "Navigating images using modified inputs: "
                f" include_pattern={actual_incl_pattern}"
                f" exclude_pattern={actual_excl_pattern}"
                f" input_rootname={actual_input_rootname}"
                f" input_file_list={actual_input_file_list}"
            )
            self._logger.info(msg)

        nav_status_table = self.navigate_images(
            data_dir,
            actual_incl_pattern,
            actual_excl_pattern,
            actual_input_file_list,
            actual_input_rootname,
            input_suffix,
            final_quality_file,
            clean_star_detection,
            clean_astrometry,
            stop_on_error,
            extnum,
            search_fwhm,
            search_nsigma,
            detector_bitdepth,
            max_sources,
            nosatmask,
            sat_frac,
            quiet,
            exclude_corner_pct,
            srclist_extname,
            astnet_key,
            use_sip,
            user_scale,
            scale_err_ratio,
        )

        # Resample and image mosaicing/stacking
        if target_wcs_file is None:
            target_wcs_file = self._navfile_files[0]
            self._logger.info(
                (
                    f"Adopting {target_wcs_file} as the target"
                    " WCS file for resample_images_to_match."
                )
            )
        else:
            self._logger.info(
                (
                    f"Using input {target_wcs_file} as the target"
                    " WCS file for resample_images_to_match."
                )
            )

        resampled_images, resampled_info_table = self.resample_images_to_match(
            target_wcs_file,
            resampled_file_prefix,
            resampled_file_suffix,
            resampled_dir,
            filter_list,
            resampled_summary_plot=resampled_summary_plot,
            composite_summary_plot=composite_summary_plot,
            combine_type=combine_type,
        )

        procall_tend = time.perf_counter()
        procall_telapsed = procall_tend - procall_tstart  # seconds
        self._logger.info((f"Finished. Processing took {procall_telapsed:.3f} seconds."))

        return (
            nav_status_table,
            resampled_images,
            resampled_info_table,
            preprocessed_modified_files,
        )

    def navigate_images(
        self,
        data_dir,
        include_pattern=None,
        exclude_pattern=None,
        input_file_list=None,
        input_rootname=None,
        input_suffix=".fits",
        final_quality_file=None,
        clean_star_detection=False,
        clean_astrometry=False,
        stop_on_error=False,
        extnum=0,
        search_fwhm=3.0,
        search_nsigma=7.0,
        detector_bitdepth=16,
        max_sources=200,
        nosatmask=True,
        sat_frac=0.8,
        quiet=True,
        exclude_corner_pct: float | None = None,
        srclist_extname="AP_XYPOS",
        astnet_key=None,
        use_sip=False,
        user_scale=None,
        scale_err_ratio=None,
    ):
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
            ``nav_status_table``. If you want it to instead stop processing
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
        exclude_corner_pct : float or None, optional, default=None
            If a positive value in the range [0,100] is specified then semi-circular
            regions of radius (exclude_corner_pct/100)*max(img_rows, img_cols)
            around each corner will be excluded from source detection. This
            parameter should be used if the image suffers from poor flat-fielding
            or vignetting that creates spurious sources at the image corners.
            (ApFindStars parameter)
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
            If specified, override the estimated plate  scale in the source list file
            and instead use a user spacified estimate of the plate scale.
            The units are arcseconds/pixel.
            (ApAstrometry parameter)
        scale_err_ratio: float, optional, default=None
            The relative uncertainty in the estimated plate scale,
            expressed as a ratio. This applies to either the default
            estimate from the source list, or a user-supplied plate
            scale. For example, if the estimated plate  scale is 2.0 arcsec/pix
            and the scale_err_ratio=1.5, then the plate scale range that
            will be search by Astrometry.net is 2/1.5 (=4/3) to 2*1.5 (=3)
            arcseconds. If not specified ApAstrometry will use a value of 1.3.
            Using a larger value can help in cases where astrometric
            solutions fail, for example if incorrect telescope metadata
            leads to inaccurate estimated plate scales.
            (ApAstrometry parameter)

        Returns
        -------
        nav_status_table : astropy.table
            An astropy table containing the processing status of the images

        See Also
        --------
        get_file_names : Given input file names returns the file names associated
            with a given output file type.
        """

        # Generate an image file collection
        keys = ["naxis1", "naxis2", "imagetyp", "object", "filter", "exposure"]

        self._logger.debug(f"Current working directory: {os.getcwd()}")
        self._logger.info(f"Attempting to find stars in FITS files within directory={data_dir}")
        if input_file_list is None:
            # Use patterns
            self._logger.debug(f'Using files that match include_pattern="{include_pattern}"')
            self._logger.debug(f'Excluding files that match exclude_pattern="{exclude_pattern}"')
            ifc_cal = ImageFileCollection(
                data_dir,
                keywords=keys,
                glob_include=include_pattern,
                glob_exclude=exclude_pattern,
                ext=extnum,
            )
        else:
            # Use explicit file list
            self._logger.debug(f"Using specified file list: {input_file_list}")
            ifc_cal = ImageFileCollection(
                data_dir, keywords=keys, filenames=input_file_list, ext=extnum
            )
        self._input_ifc = ifc_cal
        self._data_dir = Path(data_dir)
        self._extnum = extnum

        # Generate status table, written to self._nav_status_table
        num_inputs = len(ifc_cal.summary)
        self._make_nav_status_table(num_inputs)

        self._logger.info(f"There are {num_inputs} input files matching the parameters given.")
        self._logger.debug(f"Input file collection:\n{ifc_cal.summary}")

        # Iterate over the input files
        idx = 0
        iskipped_nav = 0  # Number of skipped navigations
        iskipped_ast = 0  # Number of skipped astrometic solutions
        name_and_dir = True  # We want _generate_all_output_names to include the sub-directory
        mkdir = True  # We want subdirectories made if not present
        proc_tstart = time.perf_counter()
        for hdu, fname in ifc_cal.hdus(return_fname=True):
            status = "pass"
            self._logger.debug(80 * "-")
            self._logger.info(f"Processing input file {fname}")

            f_srclist, f_regfile, f_plotfile, f_fwhmplot, f_qualfile, f_navfile = (
                self._generate_all_output_names(
                    fname, input_rootname, input_suffix, name_and_dir, mkdir
                )
            )

            # Determine whether to perform star detection, based on p_clean and presence of
            dont_throw = False
            does_srclist_exist = self._check_file_exists(f_srclist, dont_throw)
            if clean_star_detection or not does_srclist_exist:
                if clean_star_detection:
                    # Remove all find star outputs, not just the srclist
                    for file_to_remove in [
                        f_srclist,
                        f_regfile,
                        f_plotfile,
                        f_fwhmplot,
                        f_qualfile,
                    ]:
                        if self._check_file_exists(file_to_remove, dont_throw):
                            self._logger.debug(
                                (
                                    f"Removing existing {file_to_remove} because "
                                    f"clean_star_detection={clean_star_detection}"
                                )
                            )
                            Path(file_to_remove).unlink()

                self._logger.debug(f"Finding stars, output source list: {f_srclist}")
                fs_tstart = time.perf_counter()
                status = "pass"
                try:
                    self._find_stars_wrapper(
                        fname,
                        f_srclist,
                        extnum,
                        search_fwhm,
                        search_nsigma,
                        detector_bitdepth,
                        sat_frac,
                        max_sources,
                        nosatmask,
                        f_plotfile,
                        quiet,
                        exclude_corner_pct,
                        f_fwhmplot,
                        f_qualfile,
                        f_regfile,
                    )
                except Exception as e:
                    self._logger.warning(f"Error, caught exception when processing {fname}")
                    status = "fail"
                    if stop_on_error:
                        raise e

                fs_tend = time.perf_counter()
                fs_telapsed = fs_tend - fs_tstart  # seconds
                fs_status = status
            else:
                self._logger.debug(
                    (
                        f"Skipping star detection because {f_srclist} exists"
                        f" and clean_star_detection={clean_star_detection}"
                    )
                )
                fs_status = "skipped_exists"
                fs_telapsed = 0
                iskipped_nav += 1

            # Update lists of existing files of each type we've just generated
            # in source searching
            self._update_file_lists(
                {
                    "srclist": f_srclist,
                    "regfile": f_regfile,
                    "plotfile": f_plotfile,
                    "qualfile": f_qualfile,
                    "fwhmplot": f_fwhmplot,
                }
            )

            # Determine whether to perform astrometry, based on p_clean
            # and presence of srclist and navfile
            does_navfile_exist = self._check_file_exists(f_navfile, dont_throw)
            does_srclist_exist = self._check_file_exists(f_srclist, dont_throw)
            if does_srclist_exist:
                if clean_astrometry or not does_navfile_exist:
                    if clean_astrometry:
                        for file_to_remove in [f_navfile]:
                            if self._check_file_exists(file_to_remove, dont_throw):
                                msg = (
                                    f"Removing existing {file_to_remove}"
                                    f" as clean_star_detection={clean_astrometry}"
                                )
                                self._logger.debug(msg)
                                Path(file_to_remove).unlink()

                    ast_status = "pass"
                    ast_tstart = time.perf_counter()
                    self._logger.debug(
                        f"Performing astrometry, output navigated image: {f_navfile}"
                    )
                    try:
                        ap_astrom = ApAstrometry(
                            fname,
                            f_srclist,
                            f_navfile,
                            inp_img_extnum=extnum,
                            srclist_extname=srclist_extname,
                            astnet_key=astnet_key,
                            use_sip=use_sip,
                            user_scale=user_scale,
                            scale_err_ratio=scale_err_ratio,
                            loglevel=self._loglevel,
                        )
                        p_status = ap_astrom.status()
                        self._logger.debug(f"ApAstrometry return status: {p_status}")
                        if p_status == ApAstrometry.NOMINAL:
                            ast_status = "pass"
                        elif p_status == ApAstrometry.INPUT_ERROR:
                            ast_status = "input_error"
                        else:
                            ast_status = "fail"
                    except Exception as e:
                        self._logger.warning(
                            f"Error, caught exception when performing astrometry on {fname}"
                        )
                        ast_status = "exception"
                        if stop_on_error:
                            raise e

                    ast_tend = time.perf_counter()
                    ast_telapsed = ast_tend - ast_tstart  # seconds
                    ast_status = status
                else:
                    self._logger.debug(
                        (
                            f"Skipping astrometry because {f_navfile} exists"
                            f" and clean_astrometry={clean_astrometry}"
                        )
                    )
                    ast_status = "skipped_exists"
                    ast_telapsed = 0
                    iskipped_ast += 1

                self._update_file_lists({"navfile": f_navfile})
            else:
                # No source list
                self._logger.error(f"Skipping astrometry because {f_srclist} does not exist.")
                ast_status = "skipped_no_srclist"
                ast_telapsed = 0
                iskipped_ast += 1

            # Fill in run info for this file
            result_tuple = (fname, fs_status, fs_telapsed, ast_status, ast_telapsed)
            self._nav_status_table[idx] = result_tuple

            # update idx
            idx += 1

        # Generate overall quality file summary
        qual_file_dir = self.get_directories("qualfile", data_dir, absolute=False)
        self.create_quality_summary(qual_file_dir, final_quality_file)

        proc_tend = time.perf_counter()
        proc_telapsed = proc_tend - proc_tstart
        self._logger.info(f"Finished processing {idx} files in {proc_telapsed:.3f} seconds.")
        self._logger.debug(
            (
                f"Skipped {iskipped_nav}/{(idx)} navigations"
                f" and {iskipped_ast}/{(idx)} astrometric solutions."
            )
        )
        return self._nav_status_table

    def create_quality_summary(self, qual_file_dir: str, quality_summary_file: str) -> None:
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
        summarizer = ApQualitySummarizer(  # noqa: F841
            qual_file_dir,
            quality_summary_file,
            self._loglevel,
            walk_tree,
            self.qual_pref,
            self.qual_suff,
        )
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

        oname_type = ["srclist", "regfile", "plotfile", "fwhmplot", "qualfile", "navfile"]
        for key, val in added_file_dict.items():
            if key in oname_type:
                if self._check_file_exists(val, False):
                    if key in "srclist":
                        self._srclist_files.append(val)
                        continue
                    if key in "regfile":
                        self._regfile_files.append(val)
                        continue
                    if key in "plotfile":
                        self._plotfile_files.append(val)
                        continue
                    if key in "fwhmplot":
                        self._fwhmplot_files.append(val)
                        continue
                    if key in "qualfile":
                        self._qualfile_files.append(val)
                        continue
                    if key in "navfile":
                        self._navfile_files.append(val)
                        continue
            else:
                self._logger.warning(
                    f"Unexpected file type key ({key}) supplied to _update_file_lists"
                )

        return

    def _find_stars_wrapper(
        self,
        a_fitsimg,
        a_fitstbl,
        an_extnum=0,
        a_search_fwhm=3.0,
        a_search_nsigma=7.0,
        a_detector_bitdepth=16,
        a_sat_frac=0.8,
        a_max_sources=200,
        do_nosatmask=True,
        a_plotfile=None,
        do_quiet=False,
        exclude_corner_pct: float | None = None,
        a_fwhm_plot=None,
        a_qual_rprt=None,
        a_regfile=None,
    ):
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
            Extension number or name for the extension holding the image data.
            Usually this is 0, for the ``PrimaryHDU``.
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
            If not None then this is the name for an output PNG plot of the
            image with detected sources plotted as circles.
        do_quiet : bool, default=False
            If True this suppresses the runtime source list printing to STDOUT
        exclude_corner_pct : float or None, optional, default=None
            If a positive value in the range [0,100] is specified then semi-circular
            regions of radius (exclude_corner_pct/100)*max(img_rows, img_cols)
            around each corner will be excluded from source detection. This
            parameter should be used if the image suffers from poor flat-fielding
            or vignetting that creates spurious sources at the image corners.
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
        find_stars = ApFindStars(
            a_fitsimg,
            an_extnum,
            a_search_fwhm,
            a_search_nsigma,
            a_detector_bitdepth,
            a_max_sources,
            do_nosatmask,
            a_sat_frac,
            self._loglevel,
            a_plotfile,
            do_quiet,
            exclude_corner_pct=exclude_corner_pct,
        )

        # Measure 2-Gaussian FWHM for select stars, get average over x and y
        # We don't generate a plot for this step because we will rerun this
        # function later with updated source search results
        (a_new_fwhm, a_madstd_fwhm, a_npts) = find_stars.measure_fwhm(None, "both")
        self._logger.debug(
            (
                f"Initial star detection and fitting: FWHM={a_new_fwhm:.3f}"
                f" +/- {a_madstd_fwhm:.3f} pixels using {a_npts} stars."
            )
        )

        # If we could measure a FWHM then we should refine the source detection
        if a_npts > 0:
            # Refine source detection
            self._logger.debug(
                (
                    f"Updating source searching using initial FWHM={a_new_fwhm:.3f}"
                    f" +/- {a_madstd_fwhm:.3f} pixels using {a_npts} stars."
                )
            )
            find_stars.source_search(a_new_fwhm, a_search_nsigma)

            # Re-run photometry
            find_stars.aperture_photometry()

            # Measure 2-Gaussian FWHM for select stars, get average over x and y
            (a_new_fwhm2, a_madstd_fwhm2, a_npts2) = find_stars.measure_fwhm(a_fwhm_plot, "both")
            self._logger.debug(
                (
                    f"Final star detection and fitting: FWHM={a_new_fwhm2:.3f}"
                    f" +/- {a_madstd_fwhm2:.3f} pixels using {a_npts2} stars."
                )
            )

            # As the source searching and photometry was redone, we should redo
            # the plotting.
            if a_plotfile is not None:
                find_stars.plot_image(a_plotfile)
        else:
            err_msg = (
                "Failed to find measure any source extents."
                " Suggest reprocessing with different parameters"
            )
            raise RuntimeError(err_msg)

        # Write optional quality report
        if a_qual_rprt is not None:
            find_stars.write_quality_report(a_qual_rprt)

        # Write optional ds9 format region file
        if a_regfile is not None:
            find_stars.write_ds9_region_file(a_regfile)

        # Write final sourcelist with photometry.
        find_stars.write_source_list(a_fitstbl)
        if self._check_file_exists(a_fitstbl, True):
            self._logger.debug(f"Confirming output srclist {a_fitstbl} was created.")
        else:
            self._logger.error(f"Expected output srclist {a_fitstbl} not found.")
        return

    def resample_images_to_match(
        self,
        target_wcs_file: str,
        resampled_file_prefix: str,
        resampled_file_suffix: str,
        resampled_dir: str | None = None,
        filter_list: list[str] | None = None,
        resampled_summary_plot: str | None = None,
        composite_summary_plot: str | None = None,
        combine_type: str = "WEIGHTED",
    ) -> tuple[Any, Any]:
        """
        Resample and combine all navigated images to match the WCS defined
        in the specified file, separating outputs by FILTER.

        By default all filters found in the navigated images will be
        processed, but a user-specified list of filters can be specified
        instead.

        The function parameters ``resampled_dir``, ``resampled_file_prefix``,
        and ``resampled_file_suffix`` control the location and name of the
        resampled stacked output images.
        The output resampled stacked images are placed in the
        subdirectory ``resampled_dir`` with file names of the form
        ``<resampled_file_prefix><filter><resampled_file_suffix>``.

        Any existing files with the same name as the output will be deleted
        and regenerated. This behavior differs from :func:`navigate_images`
        because changing the optional parameters of ``resample_images_to_match``
        can change the output image drastically.

        In addition to the FITS-format resampled images, a rasterizer
        "quick look" version of each image will be generated using
        :class:`ApFitsRasterizer` and the ``swarp`` backend. The output
        rasterized images will have the same name as the resampled FITS
        output files bu with a ``.tiff`` extension.

        Parameters
        ----------
        target_wcs_file : str
            File name, including path, to the FITS file than contains
            the WCS that all navigated images should be resampled to
            match.
        resampled_file_prefix : str
            The part of the resampled image file name in front of the
            filter name. For example, for a ``Green`` filter output
            image name of ``M101_Supernova_2023ixf_Green_resamp_weighted.fits``
            then ``resampled_file_prefix='M101_Supernova_2023ixf_'``.
        resampled_file_suffix : str
            The part of the resampled image file name after the filter
            name. **Note** that this should include the file suffix, which
            should include ``.fits`` (or further, e.g. ``.fits.gz``)
            For example, for a ``Green`` filter output
            image name of ``M101_Supernova_2023ixf_Green_resamp_weighted.fits``
            then ``resampled_file_suffix='_resamp_weighted.fits'``.
        resampled_dir : str, optional, default='.'
            Directory in which resampled output images will be written.
            If None then the current directory is used. If not None then
            the named directory will be created if it does not exist.
        filter_list : list of str or None, optional, default=None
            If specified then only images with ``FILTER`` header keywords
            matching one of the input list strings will be resampled. For
            example, to only resample the H-alpha and luminance images
            the filter list would be ``['Lum', 'Ha']``
        resampled_summary_plot: str or None, optional, default=None
            If not None then a PNG file with the specified name containing plots
            of each output resampled image will be generated. The plots are
            generated using :func:`util.load_three_images_and_plot`.
            Default intensity scaling and WCS plotting options will be used.
        composite_summary_plot: str or None, optional, default=None
            If there are resampled images in all three Red, Green, and Blue
            filters and/or all three Ha, SII, and OIII filters and this
            argument is not None then a PNG file with summary plots of three
            color composities will be generated. If Red, Green, and Blue
            filters are present then an ``RGB`` three color composite plot
            will be generated. If the SII, H-alpha, and OIII filters
            are present then an ``SHO`` color-scheme composite plot
            will be generated. The plots are
            generated using :func:`util.plot_lupton_threecolor`.
            Default intensity scaling and WCS plotting options will be used.
        combine_type : str, optional, default='WEIGHTED'
            One of the valid swarp ``COMBINE_TYPE`` values, i.e. one
            of ``AVERAGE, CHI2, MEDIAN, MIN, MAX, SUM, WEIGHTED``.
            I recommened including some information on the chosen type in your
            ``resampled_file_suffix``, e.g. in the example above the suffix was
            ``resampled_file_suffix='_resamp_weighted.fits'``.

        Returns
        -------
        resampled_images : dict of str, str
            A dictionary consiting of filter name (key), file path
            and name of the resample images (value) pairs, as
            generated by this function
        resampled_info_table : astropy.table.Table
            Informational table of filters, number of navigated images,
            exposure time stats, and optionally any resampled outputs.

        See Also
        --------
        `swarp manual <https://raw.githubusercontent.com/astromatic/swarp/legacy_doc/prevdoc/swarp.pdf>`_ :
            PDf manual for Astromatic ``swarp``
        get_filter_list :
            Returns a list of the filters found in the current image set
        set_file_names :
            Set the file name parameters for `ApProcess` and find all the
            existing files than match the expected file names.
        get_file_names :
            Returns the file names of different processing stages associated
            with `ApProcess`.


        Notes
        -----

        Astromatic ``swarp`` must be installed. The rasterized TIFF-format
        quick look images also require that ``swarp`` is installed.

        Images must have valid WCS headers in order to be resampled and
        stacked, either from external sources or generated by
        `ApProcess.navigate_images`.

        By default this function will process ``navfile`` type files
        associated with earlier processing by this ApProcess instance,
        most likely by `ApProcess.navigate_images`.

        If you wish to resample pre-existing navigated images, this is
        possible if the files were created by older runs of `ApProcess` and/or
        follow the `ApProcess` file naming conventions.

        You should first
        call :func:`set_file_names`. This will populate the internal
        file lists with all existing ``navfile``-type files matching
        the specified patterns. A subsequent call to `resample_images_to_match`
        will use those files.

        Warnings
        --------

        If the navigated files do not follow the `ApProcess` file naming
        conventions, for example they were generated by other software,
        then there is *currently* no simple solution to process them
        with `ApProcess` without manipulating the file names by hand. A
        future version of `ApProcess` will provide a more generic
        `ApResample` class that can handle this case.

        Astromatic `swarp` has many options. This function currently
        only uses the ``COMBINE_TYPE=WEIGHTED`` and ``FSCALASTRO_TYPE=FIXED``
        methods, which appear to be the best "general" methods for images
        without significant distortion (as `swarp` cannot handle the SIP
        distortion types supported by `Astrometry.net`). This may be
        fine for visualization, but the statistical and photometric
        accuracy of the output has not been investigated enough to make
        any claims in that regard.
        """  # noqa: E501

        res_tstart = time.perf_counter()

        # Check file we want to match exists, throw if it does not
        self._logger.debug(f"Checking if WCS target file {target_wcs_file} exists.")
        self._check_file_exists(target_wcs_file, True)

        if len(self._navfile_files) == 0:
            err_msg1 = "There are no images of type navfile known to the ApProcess instance."
            err_msg2 = (
                "Use get_file_names('navfile', ap_filestate='processed') to see known files."
            )
            err_msg3 = (
                "Use set_file_names to set input files and pick up existing navigated images."
            )
            err_msg4 = "Or run navigate_images to generate navfile images from calibrated inputs."
            self._logger.error(err_msg1)
            self._logger.error(err_msg2)
            self._logger.error(err_msg3)
            self._logger.error(err_msg4)
            raise RuntimeError(err_msg1)

        if filter_list is None:
            filter_list = self.get_filter_list("navfile")
        self._filter_list = filter_list

        if len(filter_list) == 0:
            err_msg = "Empty list of unique filters for images of type navfile"
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)

        # swarp option info
        # - `CELESTIAL_TYPE`: This needs to be set to native to copy the first input image.
        # - `VERBOSE_TYPE`: By default `NORMAL`, options are `QUIET`, `LOG`, `NORMAL`, `FULL`.
        # - `COMBINE_TYPE`: Median is good for a first look but will increase
        #   the variance. If possible use `WEIGHTED`
        # - `PIXEL_SCALE_TYPE`: Various options, see `swarp`
        #   [manual](https://raw.githubusercontent.com/astromatic/swarp/legacy_doc/prevdoc/swarp.pdf).
        #   May be set to `MAX`
        # - - Note this is ignored if a `.head` file is supplied
        # - `CENTER_TYPE`: Instead of `MANUAL` we want `MOST` or `ALL`.
        #   The latter options also generate output images that are aligned
        #   North up, East to the left.
        # - - Hopefully this is ignored if we specify a `.head`?
        # - `CENTER` and `PIXEL_SCALE` wont be used if we change `PIXEL_SCALE_TYPE`
        #    and `CENTER_TYPE`
        # - `IMAGE_SIZE`: Use 0 for automatic sizing based on other parameter values
        # - `SUBTRACT_BACK`: Using `N` is appropriate if we have a lot of nebular
        #   emission, crowded fields, large galaxies, or have already performed
        #   background subtraction. You should really try with both 'Y' and 'N'
        #   to determine the effect of swarp's in-built background subtraction.
        # - `RESAMPLE_DIR`: Directory where temporary files are written. This must exist.
        # center_type : {'MANUAL', 'MOST', 'ALL', 'FIRST'}, optional, default='FIRST'
        # Controls how the output resampled image is centered given a set
        # of inputs that cover different parts of the sky.
        # See the swarp documentation for ``ALL, MOST, MANUAL``. The
        # ``FIRST`` option is special to this application, where
        # the output image footprint is exactly constrained to match
        # ``target_wcs_file``.
        # Note that use of ``MOST`` or ``ALL`` will force the output
        # image to be aligned North up, East left, irrespective of the
        # orientation in ``target_wcs_file``.

        # check the combine type
        allowed_combine_types = ["AVERAGE", "CHI2", "MEDIAN", "MIN", "MAX", "SUM", "WEIGHTED"]
        swarp_combine_type = combine_type.upper()  # = "MEDIAN"   # WEIGHTED MEDIAN etc
        if swarp_combine_type not in allowed_combine_types:
            err_msg = (
                f"Unexpected combine_type={combine_type} specified"
                f", should be one of {allowed_combine_types}"
            )
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)

        swarp_center_type = "FIRST"  # MANUAL, MOST, ALL, FIRST (MANUAL and ALL are swarp types)
        dothead_format = "fits"  # 'text' or 'fits
        swarp_verbose = False  # Echo swarp input string if True, echo .head for FIRST case
        resampled_images: dict[str, str] = {}
        raster_images: dict[str, str] = {}

        # Check target WCS file exists
        crot: float = 0
        if not self._check_file_exists(target_wcs_file, False):
            err_msg = "Error, target_wcs_file {target_wcs_file} not found."
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        else:
            # otherwise get some information on the file
            with fits.open(target_wcs_file) as hdulist:
                ipwcs = wcs.WCS(hdulist[self._extnum].header)
                dothead_str, cd1as, cd2as, xsiz1am, ysiz2am, crot = util.summarize_wcs(
                    ipwcs, False
                )
                self._logger.info(
                    f"Target WCS file is {xsiz1am:.2f} x {ysiz2am:.2f}"
                    f" arcmin (wxh) with {cd1as:.2f} x {cd2as:.2f} arcsec pixels"
                    f" at an angle of {crot:.2f} deg from North."
                )

        if resampled_dir is not None:
            mkdir = True
            dir_exists = self._check_dir_exists(resampled_dir, mkdir, False)
            if not dir_exists:
                err_msg = (
                    f"Warning, directory {resampled_dir} does not exist"
                    f" but you specified mkdir={mkdir}."
                )
                self._logger.warning(err_msg)

        for filter in filter_list:
            files_in_filter = self._get_files_matching(
                self._navfile_files, self._extnum, "FILTER", filter
            )
            self._logger.info(f"Processing {len(files_in_filter)} images for filter {filter}")

            # Names for output resampled image, net weights, and .head files
            ofilename = f"{resampled_file_prefix}{filter}{resampled_file_suffix}"
            if resampled_dir is not None:
                ofilename = resampled_dir + "/" + ofilename
            owgtsname = ofilename.replace(".fits", "_weights.fits")
            Path(ofilename).unlink(missing_ok=True)  # remove old file
            oheadname = ofilename.replace(".fits", ".head")
            Path(oheadname).unlink(missing_ok=True)  # remove old file

            filtered_netexp = 0
            filtered_numimg = 0
            inp_weights = []
            inp_files = []

            first_file_kw_dict = None
            for fname in files_in_filter:
                # Get a broad range of navigated file header keywords that
                # we want added to the eventual resampled mosaic
                if first_file_kw_dict is None:
                    first_file_kw_dict = self._read_navfile_keywords(fname, self._extnum)

                hdr = fits.getheader(fname, extnum=self._extnum)
                texp = util.get_exposure_time(hdr)
                if texp is None:
                    raise RuntimeError(f"Error, could not get exposure time from {fname}")
                    continue
                inp_files.append(fname)
                fscale = 1.0 / texp
                inp_weights.append(f"{fscale:.7f}")
                filtered_netexp += texp
                filtered_numimg += 1

            if len(inp_files) == 0:
                self._logger.warning(
                    f"No files with non-zero exposure times found for filter {filter}"
                )
                continue
            else:
                # If we could not read any keywords from the navigated file
                # we need to at least add the FILTER keyword to the output
                # resampled mosaic
                if first_file_kw_dict is None:
                    first_file_kw_dict = {"FILTER": (filter, "Filter used")}

                fs_tstart = time.perf_counter()
                self._logger.info(
                    (
                        f"Filter {filter} with {filtered_numimg} images"
                        f" totalling {filtered_netexp} seconds."
                    )
                )
                self._logger.debug(f"Input files:   {inp_files}")
                self._logger.debug(f"Input weights: {inp_weights}")

                file_str = " ".join(inp_files)
                fscale_str = ",".join(inp_weights)
                test_cmd1 = (
                    f"swarp {file_str} -FSCALASTRO_TYPE VARIABLE -FSCALE_DEFAULT {fscale_str}"
                )
                test_cmd2 = (
                    f"-VERBOSE_TYPE FULL -SUBTRACT_BACK N -COMBINE_TYPE {swarp_combine_type}"
                    " -GAIN_DEFAULT 1.0 -GAIN_KEYWORD EGAIN -RESAMPLING_TYPE LANCZOS3"
                    " -OVERSAMPLING 4 -PROJECTION_TYPE TAN"
                )
                if "MANUAL" in swarp_center_type:
                    Path(oheadname).unlink(missing_ok=True)
                    test_cmd3 = (
                        "-PIXELSCALE_TYPE MANUAL -PIXEL_SCALE 0.62"
                        " -CENTER_TYPE MANUAL -CENTER 210.8022671,54.3489500"
                        " -IMAGE_SIZE 4096,4096"
                    )
                elif "ALL" in swarp_center_type:
                    test_cmd3 = "-PIXELSCALE_TYPE MAX -CENTER_TYPE ALL"
                    Path(oheadname).unlink(missing_ok=True)
                elif "MOST" in swarp_center_type:
                    test_cmd3 = "-PIXELSCALE_TYPE MAX -CENTER_TYPE MOST"
                    Path(oheadname).unlink(missing_ok=True)
                elif "FIRST" in swarp_center_type:
                    # Need to generate a header based on the first file,
                    # with a name based on the output file name
                    # The pixelscale and centertype should be set from the .head file.
                    util.make_dothead_from_file(
                        target_wcs_file,
                        oheadname,
                        extnum=self._extnum,
                        format=dothead_format,
                        verbose=swarp_verbose,
                    )

                    # swarp does honor the imagesize if set in a FITS format
                    # .head file, but not in the ASCII format it describes.
                    test_cmd3 = ""
                else:
                    raise RuntimeError(
                        f"Error, center type {swarp_center_type} has not yet been implemented."
                    )

                test_cmd4 = (
                    f"-IMAGEOUT_NAME {ofilename} -WRITE_FILEINFO Y"
                    " -WRITE_XML N -DELETE_TMPFILES Y -RESAMPLE_DIR ./"
                    f" -WEIGHTOUT_NAME {owgtsname}"
                )
                test_cmd_list = (
                    shlex.split(test_cmd1)
                    + shlex.split(test_cmd2)
                    + shlex.split(test_cmd3)
                    + shlex.split(test_cmd4)
                )
                if swarp_verbose:
                    self._logger.info(
                        f"Command string to be passed to subprocess.run is {test_cmd_list}"
                    )

                # Run swarp
                fs_tstart = time.perf_counter()
                try:
                    result = subprocess.run(
                        test_cmd_list,
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                    )
                    fs_tend = time.perf_counter()
                    fs_telapsed = fs_tend - fs_tstart  # seconds
                    self._logger.info(
                        (
                            f"Success: swarp process took {fs_telapsed:.3f} seconds"
                            f", return code {result.returncode}"
                        )
                    )

                    # Update the output image header keywords
                    self._update_fits_header(ofilename, self._extnum, first_file_kw_dict)
                    resampled_images[filter] = ofilename

                    # Create a raster image
                    raster_file_name = None
                    try:
                        raster_file_name = ofilename.replace(".fits", ".tiff")
                        self._rasterizer.fits_to_greyscale(
                            ofilename,
                            raster_file_name,
                            backend="stiff",
                            extnum=self._extnum,
                            min_cut=None,
                            max_cut=None,
                            min_percent=50,
                            max_percent=99.9,
                            negative=False,
                            binning=1,
                        )
                        raster_images[ofilename] = raster_file_name
                    except RuntimeError as e:
                        err_msg = (
                            f"Failed to generate raster file from {ofilename}. Proceeding anyway."
                        )
                        self._logger.warning(err_msg)
                        self._logger.warning(f"Error message from ApFitsRasterized was: {e}")

                except subprocess.CalledProcessError as err:
                    fs_tend = time.perf_counter()
                    fs_telapsed = fs_tend - fs_tstart  # seconds
                    self._logger.error(
                        (
                            f"  Error, process took {fs_telapsed:.3f} seconds"
                            f", return code {err.returncode}"
                        )
                    )
                    self._logger.error(f"  Input args: {err.cmd}\n")
                    self._logger.error(f"  Stdout: {err.output}\n")
                    self._logger.error(f"  Stderr: {err.stderr}\n")
                    self._logger.error(f'  Command line equivalent command: {" ".join(err.cmd)}')

        # Generate a summary
        self._resampled_image_dict = resampled_images
        resampled_info_table = self.make_filter_image_summary(
            filter_list,
            self._navfile_files,
            resampled_image_info=True,
            raster_img_dict=raster_images,
        )
        res_tend = time.perf_counter()
        res_telapsed = res_tend - res_tstart  # seconds
        self._logger.info(
            (
                f"Finished. Generated {len(resampled_images)} resampled images"
                f" in {res_telapsed:.3f} seconds: {resampled_images}"
            )
        )

        # Optional plotting
        swap_axis = False  # assumes Y-axis roughly in line with north or south
        if (45.0 < crot <= 135.0) or (225.0 < crot <= 315.0):
            self._logger.debug(f"Setting swap_axis True because crot={crot:.2f} degrees.")
            swap_axis = True

        # determine reasonable angle tick spaceing based on image size
        angle_tick_spacing: float = 10  # arcmin, reasonable for iTelescope
        tick_scale_factor: float = 4.0
        imsizam: float = max(xsiz1am, ysiz2am)
        for tick_spacing in [
            0.083333333,
            0.25,
            1.0,
            2.0,
            3.0,
            5.0,
            10.0,
            15.0,
            30.0,
            60.0,
            90.0,
            120.0,
            180.0,
            240.0,
            300.0,
            360.0,
            450.0,
            600.0,
            900.0,
        ]:
            angle_tick_spacing = tick_spacing
            if imsizam < (tick_scale_factor * tick_spacing):
                break
        self._logger.info(
            f"Tick spacing chosen for image size of {imsizam:.2f} arcmin"
            f" is {tick_spacing:.2f} arcmin."
        )

        qval = 8
        stretch = 0.5
        if resampled_summary_plot is not None:
            imlist: list[str] = [resampled_images[key] for key in resampled_images.keys()]
            nimgs = len(imlist)

            self._logger.debug(
                f"Generating a summary plot of {nimgs} using "
                f"WCS info, tick spacing {angle_tick_spacing} arcmin,"
                f" and swap_radec_axis={swap_axis}."
            )
            util.load_imlist_and_plot(
                imlist,
                extnum=self._extnum,
                output=resampled_summary_plot,
                usewcs=True,
                vmin=None,
                vmax=None,
                xaxlim=None,
                yaxlim=None,
                verbose=True,
                angle_tick_spacing_am=angle_tick_spacing,
                swap_radec_axis=swap_axis,
            )
            self._logger.info(f"Resampled image summary plot generated: {resampled_summary_plot}")
        if composite_summary_plot is not None:
            self._logger.debug(
                "Generating a summary color composite plot using "
                f"WCS info, tick spacing {angle_tick_spacing} arcmin,"
                f" and swap_radec_axis={swap_axis}."
            )
            util.make_lupton_threecolor_plots(
                resampled_images,
                composite_summary_plot,
                extnum=self._extnum,
                usewcs=True,
                vmin=None,
                xaxlim=None,
                yaxlim=None,
                qval=qval,
                stretchval=stretch,
                verbose=True,
                angle_tick_spacing_am=angle_tick_spacing,
                swap_radec_axis=swap_axis,
            )
            self._logger.info(f"Summary color composite plot generated: {composite_summary_plot}")

        return resampled_images, resampled_info_table

    def get_filter_list(self, ap_filetype="input"):
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
        keys = ["naxis1", "naxis2", "imagetyp", "object", "filter", "exposure"]

        if ap_filetype in "input":
            if len(self._input_ifc.summary) == 0:
                err_msg = "No input files have been specified."
                self._logger.error(err_msg)
                return filter_list
            ifc_filtertype = self._input_ifc
        elif ap_filetype in "navfile":
            if len(self._navfile_files) == 0:
                err_msg = (
                    f"No {ap_filetype} files have been generated by this instance of ApProcess."
                )
                self._logger.error(err_msg)
                return filter_list
            ifc_filtertype = ImageFileCollection(
                self._data_dir, keywords=keys, filenames=self._navfile_files
            )
        else:
            err_msg = f"Unrecognized file type {ap_filetype} specified."
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)

        self._logger.debug(
            (
                f"get_filter_list: there are {len(ifc_filtertype.summary)} files"
                f" of type {ap_filetype} found under {self._data_dir}"
            )
        )

        # Get the unique filters, then build a table of filter, num_images, net_exposure_time
        uniq_filter_exposure = unique(
            ifc_filtertype.summary, keys=["filter"], keep="first", silent=True
        )
        filter_list = []
        for filter in uniq_filter_exposure["filter"]:
            filter_list.append(filter)

        # Generate a valid summary table
        file_list = self._get_processed_file_names(ap_filetype, self._data_dir, name_and_dir=True)
        filter_image_summary = self.make_filter_image_summary(filter_list, file_list)
        table_str_list = filter_image_summary.pformat_all(max_lines=-1, max_width=-1)
        self._logger.debug("Unique filter and exposure time data:")
        for line in table_str_list:
            self._logger.debug(f"{line}")

        self._logger.info(
            (
                f"There are {len(filter_list)} unique filters"
                f" among {len(ifc_filtertype.summary)} files: {filter_list}"
            )
        )
        return filter_list

    def make_filter_image_summary(
        self,
        filter_list: list[str],
        img_file_list: list[str],
        resampled_image_info: bool = False,
        raster_img_dict: dict[str, str] = None,
    ):
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
            This list contains either all the ``input`` calibrated image
            file path and name or the same for the ``navfile`` images.
        resampled_image_info : bool, optional, default=False
            If True then add a column listing any combined resampled images
            (``ap_filetype='stacked'``) for those filters.
        raster_img_dict : dict of str, str or None, optional, default=None
            If ``resampled_image_info`` is ``True`` the file name of
            any raster (tiff, png, etc.) version of a given resampled file
            can be passed in using the ``raster_img_dict``. This is a
            dictionary with keys equal to the resampled FITS image file
            names and values corresponding to any rasterized version.

        Returns
        -------
        resampled_info_table : astropy.table.Table
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
        # - {optional} Raster version of resampled FITS image

        num_images = []
        net_exptime_min = []
        min_exptime_min = []
        max_exptime_min = []
        keyword = "FILTER"

        # TODO This loop is inefficient
        for filter_name in filter_list:
            self._logger.debug(f"Finding files and exposure times for filter {filter_name}")
            sumt = 0
            numf = 0
            mint = None
            maxt = None
            files_in_filter = self._get_files_matching(
                img_file_list, self._extnum, keyword, filter_name
            )
            for a_file in files_in_filter:
                # TODO extension number matches input files, not necessarily
                # the navigated fies
                hdr = fits.getheader(a_file, extnum=self._extnum)

                exptime_s = util.get_exposure_time(hdr, False)
                if exptime_s is None:
                    self._logger.warning(
                        f"File {a_file} lacks exposure time information, and will be ignored."
                    )
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

            num_images.append(numf)
            net_exptime_min.append(sumt)
            min_exptime_min.append(mint)
            max_exptime_min.append(maxt)
            self._logger.info(
                (
                    f"Found the following {len(files_in_filter)} files"
                    f" for filter {filter_name}: {files_in_filter}"
                )
            )

        # Finally contruct the table
        col_names = ["filter", "num_navfiles", "tot_exptime_m", "min_exptime_m", "max_exptime_m"]
        col_dtypes = [object, np.int32, float, float, float]
        resampled_info_table = Table(
            [filter_list, num_images, net_exptime_min, min_exptime_min, max_exptime_min],
            names=col_names,
            dtype=col_dtypes,
        )

        if resampled_image_info:
            resampled_imgs = []
            for filter_name in filter_list:
                resampled_file = self._resampled_image_dict.get(filter_name, "None")
                resampled_imgs.append(resampled_file)

            c = Column(data=resampled_imgs, name="resampled_image", dtype=object)
            resampled_info_table.add_column(c)

            if raster_img_dict is not None:
                raster_img_list = []
                for fits_file in resampled_imgs:
                    raster_img = ""
                    if fits_file in raster_img_dict:
                        raster_img = raster_img_dict[fits_file]
                    raster_img_list.append(raster_img)

                c2 = Column(data=raster_img_list, name="raster_image", dtype=object)
                resampled_info_table.add_column(c2)

        return resampled_info_table

    def _update_fits_header(self, file_to_modify, extnum, kw_dict):
        """
        Update a FITS file header, adding the keywords in the dict

        Parameters
        ----------
        file_to_modify : str
            Name of the FITS files to modify the FITS header keyword values
            of.
        extnum : str or int
            Extension number or name for the HDU to modify the header of.
        kw_dict : dict
            Dictionary of FITS header keyword, (value, comment) pairs
            to add or modify in the named FITS file.
        """

        exists = self._check_file_exists(file_to_modify, False)
        if not exists:
            self._logger.error(
                f"Can not modify FITS header from {file_to_modify} as it cannot be found."
            )
            return

        with fits.open(file_to_modify, mode="update", do_not_scale_image_data=True) as hdulist:
            hdr = hdulist[extnum].header
            tnow = datetime.now().isoformat(timespec="milliseconds")
            for kw in kw_dict:
                hdr[kw] = kw_dict[kw]
            hdr["HISTORY"] = f"Created by {self._name} {__version__} at {tnow}"
            hdulist.flush()

        self._logger.info(f"Modified FITS header keywords in {file_to_modify}")
        return

    def _read_navfile_keywords(self, source_file, source_extnum):
        """
        Read optional keywords from the image file FITS header hdr,
        and add to the kw_dict dictionary if they are present.

        Parameters
        ----------
        source_file : str
            Name of the FITS files to get the FITS header keyword values
            from. In most cases the file used should be a navigated
            image file, a.k.a. a ``navfile``.
        source_extnum : str or int
            Extension number or name for the HDU to read the header from.

        Returns
        -------
        kw_dict : dict
            Dictionary of FITS header keyword, (value, comment) pairs
            read from the input FITS file.
        """

        kw_dict = {}

        exists = self._check_file_exists(source_file, False)
        if not exists:
            self._logger.warning(
                f"Can not read FITS header from {source_file} as it cannot be found."
            )
            return kw_dict

        # The following may exist. We add comment values for them.
        kw_comment_dict = {
            "EXPOSURE": "[seconds] Image exposure time",
            "DATE-OBS": "Observation date and time",
            "OBJECT": "Target object",
            "OBJNAME": "Target object",
            "TELESCOP": "Telescope used",
            "INSTRUME": "Detector used",
            "CCD-TEMP": "CCD temperature at start of exposure in C",
            "APTDIA": "[mm] Diameter of telescope aperture",
            "XPIXSZ": "[micrometers] X-axis pixel scale after binning",
            "YPIXSZ": "[micrometers] Y-axis pixel scale after binning",
            "FOCALLEN": "[mm] Stated telescope focal length",
            "FILTER": "Filter used",
            "EGAIN": "[e/ADU] Gain in electrons per ADU",
            "LAT-OBS": "[deg +N WGS84] Observatory Geodetic latitude",
            "LONG-OBS": "[deg +E WGS84] Observatory Geodetic longitude",
            "ALT-OBS": "[metres] Observatort altitude above mean sea level",
            "AIRMASS": "Airmass (multiple of zenithal airmass)",
        }

        with fits.open(source_file) as hdulist:
            self._logger.info(f"Reading FITS header keywords from {source_file}")
            hdr = hdulist[source_extnum].header
            kw_missing_list = []
            kw_found_list = []
            for kw in kw_comment_dict:
                if kw in hdr:
                    kw_dict[kw] = (hdr[kw], kw_comment_dict[kw])
                    kw_found_list.append(kw)
                else:
                    kw_missing_list.append(kw)

            self._logger.debug("FITS keywords found in image: {}".format(kw_found_list))
            self._logger.debug("FITS keywords missing from image: {}".format(kw_missing_list))
        return kw_dict

    def _get_files_matching(self, img_file_list, extnum, keyword, expected_value):
        """
        Return a list of files from the input list that have a FITS header
        keyword matching the expected value in the specified extension

        Parameters
        ----------
        img_file_list : list of str
            A list of the file names, **including sub-directory names**,
            of the files for which we want a filter image summary.
        extnum : int or str
            Extension number or name for the extension holding the image data.
            Usually this is 0, for the ``PrimaryHDU``. (ApFindStars parameter)
        keyword : str
            FITS header keyword to extract the values from and compare
            to the ``expected_value``.
        expected_value : str or int or float
            The keyword value that the headers of the returned files must
            match exactly.

        Returns
        -------
        matching_files : list of str
        """

        self._logger.info(
            (
                f"Searching {len(img_file_list)} files for headers in "
                f"extension={extnum} that have {keyword}={expected_value}"
            )
        )

        matching_files = []
        for a_file in img_file_list:
            hdr = fits.getheader(a_file, extnum)

            # only process file if its FILTER matches expectation
            value = hdr[keyword]
            if value == expected_value:
                self._logger.debug(f"File {a_file} has {keyword}={expected_value}")
                matching_files.append(a_file)

        self._logger.info(f"Found {len(matching_files)} files with {keyword}={expected_value}")
        return matching_files

    def _make_nav_status_table(self, num_input_ifc):
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

        arr_dtype = [
            ("filename", object),
            ("find_stars_status", object),
            ("find_stars_time", float),
            ("astrometry_status", object),
            ("astrometry_time", float),
        ]

        self._nav_status_table = Table(data=np.empty(num_input_ifc, dtype=arr_dtype))
        self._nav_status_table["find_stars_time"].info.format = "7.3f"
        self._nav_status_table["astrometry_time"].info.format = "7.3f"
        return

    def preprocess_images(
        self,
        data_dir: str,
        preprocess_replace: str,
        preprocess_with: str,
        include_pattern: str | None = None,
        exclude_pattern: str | None = None,
        input_file_list: list[str] | None = None,
        input_rootname: str | None = None,
        input_suffix: str = ".fits",
        extnum: str | int = 0,
        find_exposure_time: bool = False,
        keyword_dict: dict[str, Any] | None = None,
        replace_keywords: bool = False,
        badpixelfile: str | None = None,
        deltapix: int = 2,
        holemaskfile: str | None = None,
        fix_cosmic_rays: bool = False,
        clean_preprocess: bool = True,
    ):
        """
        Modify the calibrated input files before performing image navigation
        and astrometry.

        This function can be used to perform the following processing
        calibrated steps in cases where the input calibrated files are
        deficient in one or more regards:

        1. Added metadata to the FITS header to identify the observation,
           telescope, or instrument characteristics. This mode is controlled
           by the ``keyword_dict`` and ``replace_keywords`` parameters.
        2. Add an ``EXPOSURE`` keyowrd value based on the value of a less
           commonly used varient already present in the FITS headers
           (e.g. ``ONTIME``). This is controlled by the ``find_exposure_time``
           parameter.
        3. Perform additional bad pixel, bad row, and/or bad column
           correction based on a user-supplied bad pixel file and using
           the :class:`ApFixBadPixels` class
           This mode is controlled
           by the ``badpixelfile`` and ``deltapix`` parameters.
        4. Inpaint "holes" in the image with values matching the local pixel
           statistics based on a user-supplied hole mask and
           :class:`ApFixHoles`. This mode is controlled by the ``holemaskfile``
           parameter.
        5. Apply Cosmic Ray (CR) rejection using :class:`ApFixCosmicRays`.
           If you decide that additional bad pixel
           processing is necessary then it is likely you will also require
           additional CR rejection as well.

        Preprocessing does not modify the original input files. Instead it
        creates new files with names based on string replacement of the of
        original inputs files. The new file names are then either explicitly
        input as an ``input_file_list`` when running ``navigate_images``,
        *or* you can modify the ``include_pattern``. In both cases
        the ``input_rootname`` should be updated.

        The pseudo-code below shows a simplified example:

        >>> # Normal (no-preprocessing) Case
        >>> #------------------------------------------------------
        >>>
        >>> data_dir        = r'.'
        >>> include_pattern = "Calibrated-*-?.fits*"
        >>> exclude_pattern = None
        >>> processor       = ap.ApProcess(loglevel)
        >>>
        >>> status          = processor.navigate_images(data_dir,
        >>>     include_pattern, exclude_pattern, ...)
        >>>
        >>> # Preprocessing Case
        >>> # ------------------------------------------------------
        >>>
        >>> data_dir        = r'.'
        >>> include_pattern = "Calibrated-*-?.fits*"
        >>> exclude_pattern = None
        >>> processor       = ap.ApProcess(loglevel)
        >>>
        >>> preprocess_instr = 'Calibrated'
        >>> preprocess_outstr = 'recalibrated'
        >>> preprocessed_modified_files = processor.preprocess_images(data_dir,
        >>>     preprocess_instr, preprocess_outstr,
        >>>     include_pattern, exclude_pattern, ...)
        >>>
        >>> # Then we can either use the modified file list and reset the
        >>> # input_rootname, e.g.
        >>>
        >>> status = processor.navigate_images(data_dir, include_pattern=None,
        >>>     input_file_list=preprocessed_modified_files, input_rootname=preprocess_outstr, ...)
        >>>
        >>> # OR we can modify the input include_pattern, e.g.
        >>>
        >>> modified_pattern = include_pattern.replace(preprocess_instr, preprocess_outstr)
        >>> print(modified_pattern)     # produces "recalibrated-*-?.fits*"
        >>>
        >>> status = processor.navigate_images(data_dir,
        >>>     modified_pattern, exclude_pattern,
        >>>     input_file_list=None, input_rootname=preprocess_outstr, ...)
        >>>

        The python dictionary ``keyword_dict`` should consist of
        `key: (value, comment)` pairs, where `key` and `comment` are
        strings. The `key` keyword must conform to FITS file conventions,
        in particular the keyword cannot be longer than 8 letters, cannot
        start with a numeral, and cannot include whitespace.

        Parameters
        ----------
        data_dir : str or path
            Path to base directory containing FITS files to process,
            e.g. `./`.
        preprocess_replace : str
            Substring common to all the input calibrated files that will
            be replaced with ``preprocess_with``.
        preprocess_with : str
            String that will replace ``preprocess_replace`` in all
            input file names.
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
        extnum : int or str, optional, default=0
            Extension number or name for the extension holding the image data
            and the FITS header keywords.
            Usually this is 0, for the ``PrimaryHDU``.
        find_exposure_time : bool, optional, default=False
            If True then :func:`ApUtil.get_exposure_time` will be used to
            extract the exposure time from the input file headers and
            set the ``EXPOSURE`` keyword if it is not already present.
        keyword_dict : dict, optional, default=None
            If not None, then supply a dictionary of
            ``keyword: (value, comment)`` entries that
            will be added to the processed file FITS headers.
        replace_keywords : bool, optional, default=False
            If True then existing FITS header keys that match the keys
            in ``keyword_dict`` will be over-written with the new values.
        badpixelfile : str, optional, default=None
            File path and name to the bad pixel file to apply. This file
            should conform to the format generated by ``ApFindBadPixels``
            and used by ``ApFixBadPixels``.
        deltapix : int, optional, default=2
            Linear distance away from a bad pixel from which
            the median value of the good pixels will be drawn. If 1 then
            the median value of good pixels within the surrounding 8 pixels
            will be used. If 2 then the median of the good pixels within the
            surrounding 24 pixels will be used. Values above 2 are not
            recommended.
        holemaskfile : str, optional, default=None
            File path and name to any hole mask file to apply. This file
            should conform to the format used by ``ApFixHoles``.
            (preprocess_images parameter)
        fix_cosmic_rays : bool, optional, default=False
            If True then perform Cosmic Ray rejection on the images.
        clean_preprocess : bool, optional, default=True
            Overwrite any existing file with the same name as a preprocessed
            output file. If False, then the presence of the output file
            will cause preprocessing to be skipped for the assciated input
            file.

        Returns
        -------
        preprocessed_modified_files : list of str
            List of file names of the output modified files
        """

        preprocessed_modified_files = []
        # Generate an image file collection
        keys = ["naxis1", "naxis2", "imagetyp", "object", "filter", "exposure"]

        self._logger.debug(f"Current working directory: {os.getcwd()}")
        self._logger.info(f"Attempting to preprocess FITS files within directory={data_dir}")
        if input_file_list is None:
            # Use patterns
            self._logger.info(f'Using files that match include_pattern="{include_pattern}"')
            self._logger.info(f'Excluding files that match exclude_pattern="{exclude_pattern}"')
            ifc_cal = ImageFileCollection(
                data_dir,
                keywords=keys,
                glob_include=include_pattern,
                glob_exclude=exclude_pattern,
                ext=extnum,
            )
        else:
            # Use explicit file list
            self._logger.info(f"Using specified file list: {input_file_list}")
            ifc_cal = ImageFileCollection(
                data_dir, keywords=keys, filenames=input_file_list, ext=extnum
            )
        self._input_ifc = ifc_cal
        num_inputs = len(ifc_cal.summary)
        self._data_dir = Path(data_dir)
        self._extnum = extnum

        self._logger.info(f"There are {num_inputs} input files matching the parameters given.")
        self._logger.debug(f"Input file collection:\n{ifc_cal}")

        bad_pixel_fixer = None
        if badpixelfile is not None:
            bad_pixel_fixer = ApFixBadPixels(self._loglevel)
            msk_data, msk_hdr = self._read_fits(badpixelfile, extnum)

        hole_fixer = None
        if holemaskfile is not None:
            hole_fixer = ApFixHoles(self._loglevel)
            hole_msk_data, hole_msk_hdr = self._read_fits(holemaskfile, extnum)

        cr_fixer = None
        if fix_cosmic_rays:
            cr_fixer = ApFixCosmicRays(self._loglevel)

        # Iterate over the input files
        idx = 0
        proc_tstart = time.perf_counter()
        for hdu, fname in ifc_cal.hdus(return_fname=True):
            # NOTE hdu is defined by the ext value given to the ImageFileCollection ctor.
            # Generate output file name
            oname = fname.replace(preprocess_replace, preprocess_with)

            self._logger.debug(80 * "-")
            self._logger.info(f"Preprocessing input file {fname} into {oname}")
            if not clean_preprocess:
                ofile_exists = self._check_file_exists(oname, throws=False)
                if ofile_exists:
                    self._logger.debug(
                        f"Skipping preprocessing of {fname}"
                        f" because {oname} exists and clean_preprocess=False."
                    )
                    preprocessed_modified_files.append(oname)
                    continue

            inp_hdr = hdu.header
            inp_data = hdu.data

            if find_exposure_time:
                if "EXPOSURE" in inp_hdr:
                    expval = inp_hdr["EXPOSURE"]
                    self._logger.debug(f"EXPOSURE keyword already set to {expval}")
                else:
                    expval = util.get_exposure_time(inp_hdr)
                    if expval is None:
                        self._logger.warning(
                            (
                                f"Could not find an exposure related header value in {fname}."
                                " Not setting EXPOSURE"
                            )
                        )
                    else:
                        inp_hdr["EXPOSURE"] = (expval, "[s] Exposure time")
                        self._logger.debug(f"Set EXPOSURE keyword to {expval}")

            if keyword_dict is not None:
                self._logger.debug(f"Adding or updating the header keywords {keyword_dict.keys()}")
                for key, val in keyword_dict.items():
                    if (key not in inp_hdr) or replace_keywords:
                        inp_hdr[key] = val

            if badpixelfile is not None:
                out_data, out_dict = bad_pixel_fixer.fix_bad_pixels(inp_data, msk_data, deltapix)
                out_dict["BPIXFILE"] = (
                    Path(badpixelfile).name,
                    "Name of master bad pixel file used",
                )
                for key, val in out_dict.items():
                    inp_hdr[key] = val
                hdu.data = out_data
                hdu.header = inp_hdr

            if holemaskfile is not None:
                out_data, out_dict = hole_fixer.fix_holes(out_data, hole_msk_data)
                out_dict["HOLEFILE"] = (
                    Path(holemaskfile).name,
                    "Name of master hole mask file used",
                )
                for key, val in out_dict.items():
                    inp_hdr[key] = val
                hdu.data = out_data
                hdu.header = inp_hdr

            if fix_cosmic_rays:
                if "EGAIN" in inp_hdr:
                    gain = inp_hdr["EGAIN"]
                else:
                    gain = 1.0
                hdu.data, cr_kw_dict = cr_fixer.process(out_data, gain)
                for key, val in cr_kw_dict.items():
                    hdu.header[key] = val

            # Finally, save the file
            hdu.writeto(oname, overwrite=True)
            self._logger.debug(f"Wrote preprocessed file {oname}")
            preprocessed_modified_files.append(oname)
            idx = idx + 1

        proc_tend = time.perf_counter()
        proc_telapsed = proc_tend - proc_tstart
        self._logger.info(f"Finished preprocessing {idx} files in {proc_telapsed:.3f} seconds.")
        self._preprocessed_modified_files = preprocessed_modified_files

        return preprocessed_modified_files
