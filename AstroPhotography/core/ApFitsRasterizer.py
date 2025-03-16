"""
Contains the implementation of the ApFitsRasterizer class.
"""

# SPDX-License-Identifier: GPL-3.0-or-later

#  2025-02-08 dks : Initial work on ApFitsRasterizer.

# import sys
import logging
from pathlib import Path
from datetime import datetime  # , timezone
import subprocess
import shlex
import time
from typing import Any

import numpy as np
import numpy.typing as npt
from astropy.io import fits

# AstroPhotography includes
from .. import __version__
from ..util.read import read_fits


class ApFitsRasterizer:
    """ """

    def __init__(self, loglevel: str) -> None:
        """Constructs an ApFitsRasterizer object. No processing is performed.

        Parameters
        ----------
        loglevel : str
            Logging level to use.
        """

        self._name = "ApFitsRasterizer"

        # Initialize logging
        self._loglevel = loglevel
        self._logger = self._initialize_logger(self._loglevel)

        self._logger.debug(f"{self._name} instance constructed.")
        return

    def _check_file_exists(self, filename: str) -> None:
        """
        Checks a file exits, otherwise raises RuntimeError
        """
        if not Path(filename).exists():
            err_msg = f"Cannot find {filename}. Not a valid path or file."
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return

    def _initialize_logger(self, loglevel: str) -> Any:
        """
        Initialize and return the logger
        """

        logger = logging.getLogger("ApFitsRasterizer")

        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError("Invalid log level: {}".format(loglevel))
        logger.setLevel(numeric_level)
        logger.propagate = False

        # check if handlers already present
        if not len(logger.handlers):
            # create console handler and set level to debug
            ch = logging.StreamHandler()
            ch.setLevel(numeric_level)

            # create formatter
            formatter = logging.Formatter("%(asctime)s | %(name)s | %(levelname)s | %(message)s")

            # add formatter to ch
            ch.setFormatter(formatter)

            # add ch to logger
            logger.addHandler(ch)
        return logger

    def _remove_pedestal_kw(self, hdr: Any) -> None:
        """
        Removes the PEDESTAL keyword from the input FITS header

        AstroPhotography always removes any artificial PEDESTAL applied
        to the data when reading a FITS file, so it is important to make
        sure that the FITS header keywords remain consistent.

        This function need only be applied when modified data is being
        written or rewritten to disk using a copy of an original FITS
        header.
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
        ounits: str,
        ohistory_str: str,
    ) -> None:
        """
        Writes the modified data to the specified output
        file, preserving all other items from the original input file.

        The output file differs from the original input data file in
        having the arithmetically modified data, possibly modified
        BUNIT value, and a history entry summarising the arithmetical
        operation that was performed.

        :param inpdata_file: Input FITS data file affected by bad pixels.
        :param ext_num: Extension number for data array and header.
        :param outdata_file: Modified copy of the input data file where the
          bad pixels have had their data valus modified by the median
          of the surrounding good pixels.
        :param odata: Bad-pixel corrected data array.
        :param ounits: BUNIT keyword to be added to the output extension.
          This may be the same as the input value, None if the input did
          not have units, or a user-supplied string.
        :param ohistory_str: String summarising the arithmatic operation
          performed to generate the data, of the form input file name,
          operation, second file name or constant.
        """

        self._logger.debug(f"Output data BUNIT keyword value: {ounits}")

        self._check_file_exists(inpdata_file)

        # open() parameters that can be important.
        # Default values used here.
        # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
        uint_handling = True
        image_scaling = False

        # TODO check this works when over-writing an image.
        hdu_list = fits.open(
            inpdata_file, uint=uint_handling, do_not_scale_image_data=image_scaling
        )

        # We do NOT modify PEDESTAL values.
        ##self._remove_pedestal_kw(hdu_list[ext_num].header)

        # Modify data
        hdu_list[ext_num].data = odata

        # Modify header
        if ounits is not None:
            hdu_list[ext_num].header["BUNIT"] = (ounits, "Pixel value units")

        tnow = datetime.now().isoformat(timespec="milliseconds")
        hdu_list[ext_num].header["HISTORY"] = f"Applied {self._name} {__version__} at {tnow}"
        hdu_list[ext_num].header["HISTORY"] = ohistory_str

        # Write to new file
        hdu_list.writeto(outdata_file, output_verify="ignore", overwrite=True)
        hdu_list.close()

        self._logger.info(f"Wrote modified data to {outdata_file}")
        return

    def fits_to_greyscale(
        self,
        input_fits: str,
        output_image: str,
        backend: str = "stiff",
        extnum: str | int = 0,
        min_cut: float | None = None,
        max_cut: float | None = None,
        min_percent: float | None = None,
        max_percent: float | None = None,
        negative: bool = False,
        binning: int = 1,
    ) -> None:
        """
        Given an input FITS file with image data in a given extension, generate
        a rasterized image (e.g. PNG or TIFF file) using the specified backend,
        optionally limiting the range to certain values or percentiles and/or
        using an inverted (negative) color table.

        Parameters
        ----------
        input_fits : str
            Input FITS image
        output_image : str
            Output raster file name. Will be overwritten if already present.
        backend : str, optional, default='stiff'
            Backend used to generate a rasterized plot.
        extnum : str or int, optional, default=0
            Extension number or name containing the data in the ``input_fits``
        min_cut : float or None, optional, default=None
            If specified, the minimum pixel value displayed will be ``min_cut`` and
            not the actual data minimum. Specify only one of ``min_cut`` and
            ``min_percentile``, not both.
        max_cut : float or None, optional, default=None
            If specified, the maximum pixel value displayed will be ``min_cut`` and
            not the actual data minimum. Specify only one of ``max_cut`` and
            ``max_percentile``, not both.
        min_percent : float or None, optional, default=None
            If specified, the minimum pixel value displayed will correspond to
            the ``min_percent`` percentile of all the data in the image,
            not the actual data minimum. Specify only one of ``min_cut`` and
            ``min_percentile``, not both.
        max_percent : float or None, optional, default=None
            If specified, the maximum pixel value displayed will correspond to
            the ``mX_percent`` percentile of all the data in the image,
            not the actual data minimum. Specify only one of ``max_cut`` and
            ``max_percentile``, not both.
        negative : bool, optional, default=False
            If ``True``, then display the image similarly to a photographic negative
            where no light is white and intense light is black.
        binning : int, optional, default=1
            Bin input pixels by this factor on each dimension before generating
            the output image. For example, with ``binning=2`` the output image
            will have half the number of rows and half the number of columns
            than the input image, and only a quarter as many pixels.
        """

        self._check_file_exists(input_fits)
        imdata, imhdr = read_fits(
            image_filename=input_fits,
            image_extension=extnum,
            a_logger=self._logger,
            remove_pedestal="never",
        )

        # determine minimum and maximum values to use.
        stat_list, percentiles = self._img_stats(
            imdata, "Input image", verbose=True, input_percentiles=(min_percent, max_percent)
        )
        min_val = stat_list[0]
        max_val = stat_list[1]
        if min_cut is not None:
            min_val = min_cut
        elif min_percent is not None:
            min_val = percentiles[0]
        if max_cut is not None:
            max_val = max_cut
        elif max_percent is not None:
            max_val = percentiles[1]

        # run the selected backend
        allowed_backends = ["stiff"]
        if "stiff" in backend:
            # Stiff does the binning somewhat differently... it must average
            # if binning > 1:
            #     min_val *= float(binning * binning)
            #     max_val *= float(binning * binning)
            user = __name__
            description = input_fits
            verbose = True
            self._logger.info(
                f"Output raster data value limits are {min_val:.3f}"
                f" to {max_val:.3f} ADU per binned pixel"
            )
            if max_val <= min_val:
                self._logger.warning(
                    "Maximum output data value limits is less than minimum value."
                )

            self._stiff_greyscale_handler(
                input_fits,
                output_image,
                min_val,
                max_val,
                negative=negative,
                binning=binning,
                copyright=user,
                description=description,
                verbose=verbose,
            )
        else:
            err_msg = (
                f"Unsupported backend {backend} specified."
                f" Allowable backends are: {allowed_backends}"
            )
            raise RuntimeError(err_msg)

        self._logger.debug("File processing completed.")
        return

    def fits_to_rgb(
        self,
        input_fits_files: list[str],
        output_image: str,
        backend: str = "stiff",
        extnum: str | int = 0,
        min_cut: list[float] | None = None,
        max_cut: list[float] | None = None,
        min_percent: list[float] | None = None,
        max_percent: list[float] | None = None,
        negative: bool = False,
        binning: int = 1,
    ) -> None:
        """
        Given three input FITS files with image data in a given extension, generate
        a rasterized RGB image (e.g. PNG or TIFF file) using the specified backend,
        optionally limiting the range to certain values or percentiles and/or
        using an inverted (negative) color table.

        Parameters
        ----------
        input_fits_files : list of str
            List of three input FITS images, in red, green, blue channel order.
        output_image : str
            Output raster file name. Will be overwritten if already present.
        backend : str, optional, default='stiff'
            Backend used to generate a rasterized plot.
        extnum : str or int, optional, default=0
            Extension number or name containing the data in the ``input_fits``
        min_cut : list of float or None, optional, default=None
            If specified, the minimum pixel value for a given channel will be ``min_cut`` and
            not the actual data minimum. Specify only one of ``min_cut`` and
            ``min_percentile``, not both.
            If not None this must be a three element list of float values
            in red, green, and blue channel order.
        max_cut : list of float or None, optional, default=None
            If specified, the maximum pixel value for a given channel will be ``min_cut`` and
            not the actual data minimum. Specify only one of ``max_cut`` and
            ``max_percentile``, not both.
            If not None this must be a three element list of float values
            in red, green, and blue channel order.
        min_percent : list of float or None, optional, default=None
            If specified, the minimum pixel value for a given channel will correspond to
            the ``min_percent`` percentile of all the data in the image,
            not the actual data minimum. Specify only one of ``min_cut`` and
            ``min_percentile``, not both.
            If not None this must be a three element list of float values
            in red, green, and blue channel order.
        max_percent : list of float or None, optional, default=None
            If specified, the maximum pixel value for a given channel will correspond to
            the ``mX_percent`` percentile of all the data in the image,
            not the actual data minimum. Specify only one of ``max_cut`` and
            ``max_percentile``, not both.
            If not None this must be a three element list of float values
            in red, green, and blue channel order.
        negative : bool, optional, default=False
            If ``True``, then display the image similarly to a photographic negative
            where no light is white and intense light is black.
        binning : int, optional, default=1
            Bin input pixels by this factor on each dimension before generating
            the output image. For example, with ``binning=2`` the output image
            will have half the number of rows and half the number of columns
            than the input image, and only a quarter as many pixels.
        """

        if len(input_fits_files) != 3:
            err_msg = (
                "Expecting a list of 3 files names"
                f", got a list with {len(input_fits_files)} elements."
            )
            raise RuntimeError(err_msg)

        # Check all files are present before proceeding

        for input_fits in input_fits_files:
            self._check_file_exists(input_fits)

        # compute stats for all files
        min_val_list = []
        max_val_list = []
        for idx, input_fits in enumerate(input_fits_files):
            imdata, imhdr = read_fits(
                image_filename=input_fits,
                image_extension=extnum,
                a_logger=self._logger,
                remove_pedestal="never",
            )

            # determine minimum and maximum values to use.
            pct_min = None
            pct_max = None
            if min_percent is not None:
                pct_min = min_percent[idx]
            if max_percent is not None:
                pct_max = max_percent[idx]
            stat_list, percentiles = self._img_stats(
                imdata, "Input image", verbose=True, input_percentiles=(pct_min, pct_max)
            )
            min_val = stat_list[0]
            max_val = stat_list[1]
            if min_cut is not None:
                min_val = min_cut[idx]
            elif min_percent is not None:
                min_val = percentiles[0]
            if max_cut is not None:
                max_val = max_cut[idx]
            elif max_percent is not None:
                max_val = percentiles[1]

            self._logger.info(
                f"Output raster data value limits are {min_val:.3f}"
                f" to {max_val:.3f} ADU per binned pixel"
            )
            if max_val <= min_val:
                self._logger.warning(
                    "Maximum output data value limits is less than minimum value."
                )
            min_val_list.append(min_val)
            max_val_list.append(max_val)

        # run the selected backend
        allowed_backends = ["stiff"]
        if "stiff" in backend:
            # Stiff does the binning somewhat differently... it must average
            # if binning > 1:
            #     min_val *= float(binning * binning)
            #     max_val *= float(binning * binning)
            user = __name__
            description = input_fits
            verbose = True

            self._stiff_rgb_handler(
                input_fits_files,
                output_image,
                min_val_list,
                max_val_list,
                negative=negative,
                binning=binning,
                copyright=user,
                description=description,
                verbose=verbose,
            )
        else:
            err_msg = (
                f"Unsupported backend {backend} specified."
                f" Allowable backends are: {allowed_backends}"
            )
            raise RuntimeError(err_msg)

        self._logger.debug("File processing completed.")
        return

    def _stiff_greyscale_handler(
        self,
        input_fits: str,
        output_tiff: str,
        min_level: float,
        max_level: float,
        gamma_type: str = "POWER-LAW",
        gamma: float = 2.2,
        gamma_fac: float = 1.0,
        color_sat: float = 1.2,
        negative: bool = False,
        binning: int = 1,
        bits_per_channel: int = 8,
        copyright: str = __name__,
        description: str = "Astronominal image",
        copy_header: bool = True,
        verbose: bool = True,
    ) -> None:
        """
        Convert the input parameters into a stiff command line string

        Notes
        -----
        See the
        `stiff manual <https://raw.githubusercontent.com/astromatic/stiff/master/doc/stiff.pdf>`_

        Parameters
        ----------
        input_fits : str
            Input FITS image
        output_tiff : str
            Name for the output TIFF file. Will be overwritten if already present.
        min_level : float
            The minimum pixel value to be displayed.
        max_level : float
            The maximum pixel value to be displayed.
        gamma_type : str, optional, default="POWER-LAW"
            Gamma correction type
        gamma : float, optional, default=2.2
            Power law slope for the gamma correction if the type is "POWER-LAW".
        gamma_fac : float, optional, default=1.0
            Additional gamma correction factor. See the ``stiff`` manual.
        color_sat : float, optional, default = 1.2
            Color saturation factor. See the ``stiff`` manual.
        negative : bool, optional, default=False
            If ``True``, then display the image similarly to a photographic negative
            where no light is white and intense light is black.
        binning : int, optional, default=1
            Bin input pixels by this factor on each dimension before generating
            the output image. For example, with ``binning=2`` the output image
            will have half the number of rows and half the number of columns
            than the input image, and only a quarter as many pixels.
            bits_per_channel: int = 8,
        copyright: str, optional, default=__name__
            The copyright field of the EXIF header.
        description: str, optional, default = "Astronominal image"
            EXIF description of the image
        copy_header: bool, optional, default = True
            If True then the FITS header of the input FITS file will be
            written into the TIFF file EXIF header. This is primary for
            interoperability with astromatic software. See the ``stiff`` manual.
        verbose: bool, optional, default = True
            If True then addition diagnostic logging of the stiff subprocess
            will be generated.
        """

        verb_type = "QUIET"
        neg_type = "N"
        copy_type = "N"
        if verbose:
            verb_type = "FULL"
        if negative:
            neg_type = "Y"
        if copy_header:
            copy_type = "Y"

        test_cmd = (
            f"stiff {input_fits}"
            f" -OUTFILE_NAME {output_tiff}"
            f" -GAMMA_TYPE {gamma_type} -GAMMA {gamma} -GAMMA_FAC {gamma_fac}"
            f" -BINNING {binning} -NEGATIVE {neg_type} -COPY_HEADER {copy_type}"
            f" -COLOUR_SAT {color_sat} -MAX_TYPE MANUAL -MAX_LEVEL {max_level}"
            f" -MIN_TYPE MANUAL -MIN_LEVEL {min_level} -SATUR_LEVEL 1e6"
            f" -BITS_PER_CHANNEL {bits_per_channel} -VERBOSE_TYPE {verb_type}"
            f" -DESCRIPTION {description} -COPYRIGHT {copyright} -WRITE_XML N"
        )

        test_cmd_list = shlex.split(test_cmd)
        success: bool = self._run_stiff(test_cmd_list, verbose=verbose)
        if success:
            self._logger.info(f"Successfully generated rasterized image {output_tiff}")
        else:
            self._logger.error(f"Failed to generate rasterized image {output_tiff}")

        return

    def _stiff_rgb_handler(
        self,
        input_fits_list: list[str],
        output_tiff: str,
        min_level_list: list[float],
        max_level_list: list[float],
        gamma_type: str = "POWER-LAW",
        gamma: float = 2.2,
        gamma_fac: float = 1.0,
        color_sat: float = 1.2,
        negative: bool = False,
        binning: int = 1,
        bits_per_channel: int = 8,
        copyright: str = __name__,
        description: str = "Astronominal image",
        copy_header: bool = True,
        verbose: bool = True,
    ) -> None:
        """
        Convert the input parameters into a stiff command line string

        Notes
        -----
        See the
        `stiff manual <https://raw.githubusercontent.com/astromatic/stiff/master/doc/stiff.pdf>`_

        Parameters
        ----------
        input_fits_list : str
            Input FITS image
        output_tiff : str
            Name for the output TIFF file. Will be overwritten if already present.
        min_level_list : float
            The minimum pixel value to be displayed.
        max_level_list : float
            The maximum pixel value to be displayed.
        gamma_type : str, optional, default="POWER-LAW"
            Gamma correction type
        gamma : float, optional, default=2.2
            Power law slope for the gamma correction if the type is "POWER-LAW".
        gamma_fac : float, optional, default=1.0
            Additional gamma correction factor. See the ``stiff`` manual.
        color_sat : float, optional, default = 1.2
            Color saturation factor. See the ``stiff`` manual.
        negative : bool, optional, default=False
            If ``True``, then display the image similarly to a photographic negative
            where no light is white and intense light is black.
        binning : int, optional, default=1
            Bin input pixels by this factor on each dimension before generating
            the output image. For example, with ``binning=2`` the output image
            will have half the number of rows and half the number of columns
            than the input image, and only a quarter as many pixels.
            bits_per_channel: int = 8,
        copyright: str, optional, default=__name__
            The copyright field of the EXIF header.
        description: str, optional, default = "Astronominal image"
            EXIF description of the image
        copy_header: bool, optional, default = True
            If True then the FITS header of the input FITS file will be
            written into the TIFF file EXIF header. This is primary for
            interoperability with astromatic software. See the ``stiff`` manual.
        verbose: bool, optional, default = True
            If True then addition diagnostic logging of the stiff subprocess
            will be generated.
        """

        verb_type = "QUIET"
        neg_type = "N"
        copy_type = "N"
        if verbose:
            verb_type = "FULL"
        if negative:
            neg_type = "Y"
        if copy_header:
            copy_type = "Y"

        input_file_str = " ".join(input_fits_list)
        min_level_str = ",".join([f"{x}" for x in min_level_list])
        max_level_str = ",".join([f"{x}" for x in max_level_list])
        test_cmd = (
            f"stiff {input_file_str}"
            f" -OUTFILE_NAME {output_tiff}"
            f" -GAMMA_TYPE {gamma_type} -GAMMA {gamma} -GAMMA_FAC {gamma_fac}"
            f" -BINNING {binning} -NEGATIVE {neg_type} -COPY_HEADER {copy_type}"
            f" -COLOUR_SAT {color_sat} -MAX_TYPE MANUAL -MAX_LEVEL {max_level_str}"
            f" -MIN_TYPE MANUAL -MIN_LEVEL {min_level_str} -SATUR_LEVEL 1e6"
            f" -BITS_PER_CHANNEL {bits_per_channel} -VERBOSE_TYPE {verb_type}"
            f" -DESCRIPTION {description} -COPYRIGHT {copyright} -WRITE_XML N"
        )

        test_cmd_list = shlex.split(test_cmd)
        success: bool = self._run_stiff(test_cmd_list, verbose=verbose)
        if success:
            self._logger.info(f"Successfully generated rasterized image {output_tiff}")
        else:
            self._logger.error(f"Failed to generate rasterized image {output_tiff}")

        return

    def _run_stiff(self, cmd_list: list[str], verbose: bool = False) -> bool:
        """
        Run stiff as a subprocess with a command list

        Parameters
        ----------
        cmd_list : list of str
            The command line arguments to be supplied to the subprocess instance.
            This should be a ``shlex`` processed list.
        verbose : bool, optional, default=False
            If True the full command list will be echo to INFO level logging

        Returns
        -------
        success : bool
            Returns true if the subprocess command succeeded.
        """

        success = False
        if verbose:
            self._logger.info(f"Command string to be passed to subprocess.run is {cmd_list}")

        # Run swarp
        fs_tstart = time.perf_counter()
        try:
            result = subprocess.run(
                cmd_list,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            fs_tend = time.perf_counter()
            fs_telapsed = fs_tend - fs_tstart  # seconds
            self._logger.debug(
                (
                    f"Success: stiff process took {fs_telapsed:.3f} seconds"
                    f", return code {result.returncode}"
                )
            )
            success = True

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

        return success

    def _img_stats(
        self,
        data: npt.NDArray,
        label: str,
        verbose: bool = False,
        input_percentiles: tuple[Any, Any] = (None, None),
    ) -> tuple[list[float], tuple[Any, Any]]:
        """
        Calculate and optionally display some image statistics that may
        include two user-specified percentiles, returning both (a) a
        a list of the minimum, maximum, mean, standard deviation,
        and median values; and (b) a tuple of the values corresponding
        to the user-specified percentiles.

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
        input_percentiles: two element tuple of float or None, optional, default = (None, None)
            The percentiles corresponding to any of the non-None elements of
            the input tuple will be computed. Note that these are percentiles
            levels in the range ``[0,100]``, not quantiles in the range ``[0, 1]``.

        Returns
        -------
        minval :float
            Minimum value in NaN-filtered input data
        maxval : float
            Maximum value in NaN-filtered input data
        meanval : float
            Mean value in NaN-filtered input data
        stdval : float
            Standard deviation in NaN-filtered input data
        medval : float
            Median value in NaN-filtered input data
        percentile_values : two element tuple of float or None
            Data values corresponding the percentile of any of the
            non-None elements of the input tuple.
        """

        minval = np.nanmin(data)
        maxval = np.nanmax(data)
        meanval = np.nanmean(data)
        stdval = np.nanstd(data)

        # percentiles, 50th percentile is the median
        #         0    1    2    3   4   5   6   7   8   9   10
        ipctls = [0.1, 1.0, 5.0, 10, 25, 50, 75, 90, 95, 99, 99.9]
        opctls = np.nanpercentile(data, ipctls)
        medval = opctls[5]

        # User-specified percentiles
        first_val = None
        if input_percentiles[0] is not None:
            first_val = np.nanpercentile(data, input_percentiles[0])
        second_val = None
        if input_percentiles[1] is not None:
            second_val = np.nanpercentile(data, input_percentiles[1])
        percentile_values = (first_val, second_val)

        if verbose:
            msg1: str = (
                f"{label} data min={minval:.2f}, max={maxval:.2f},"
                f" mean={meanval:.2f} +/- {stdval:.2f}, median={medval:.2f} ADU."
            )
            msg2: str = (
                f"  90% of data between {opctls[2]:.2f} and {opctls[8]:.2f} ADU (5-95 precentiles)"
            )
            msg3: str = (
                f"  98% of data between {opctls[1]:.2f} and {opctls[9]:.2f} ADU (1-99 precentiles)"
            )
            self._logger.debug(msg1)
            self._logger.debug(msg2)
            self._logger.debug(msg3)
            msg4 = (
                f"Requested percentile values are {first_val} ({input_percentiles[0]})"
                f" and {second_val} ({input_percentiles[1]})"
            )
            self._logger.debug(msg4)
        return [minval, maxval, meanval, stdval, medval], percentile_values
