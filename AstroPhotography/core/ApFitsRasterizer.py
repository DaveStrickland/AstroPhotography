"""
Contains the implementation of the ApFitsRasterizer class.
"""

#  2025-02-08 dks : Initial work on ApFitsRasterizer.

# import sys
import logging
from pathlib import Path
from datetime import datetime  # , timezone
from typing import Any

import numpy as np
import numpy.typing as npt
from astropy.io import fits

# AstroPhotography includes
from .. import __version__
from ..util import read_fits


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
        backend: str,
        extnum: str | int,
        min_cut: float | None,
        max_cut: float | None,
        min_percent: float | None,
        max_percent: float | None,
        negative: bool,
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
        self._logger.info(
            f"Output raster data value limits are {min_val:.3f}" f" to {max_val:.3f} ADU per pixel"
        )
        if max_val <= min_val:
            self._logger.warning("Maximum output data value limits is less than minimum value.")

        self._logger.debug("File processing completed.")
        return

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
            print(msg1)
            print(msg2)
            print(msg3)
        return [minval, maxval, meanval, stdval, medval], percentile_values
