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

        self._logger.debug("File processing completed.")
        return
