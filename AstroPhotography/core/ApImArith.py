"""
Contains the implementation of the ApImArith class.
"""

# SPDX-License-Identifier: GPL-3.0-or-later

#  2021-05-02 dks : Initial work on ApImArith.
#  2025-01-26 dks : Ruff linting, switch to utility fits reader.

# import sys
import logging
from pathlib import Path

# import math
# import time
from datetime import datetime  # , timezone
from typing import Any

import numpy as np
import numpy.typing as npt
from astropy.io import fits
# from astropy.stats import sigma_clipped_stats

# AstroPhotography includes
from .. import __version__
from ..util import read_fits


class ApImArith:
    """
    A utility class to implement the basic functionality of the HEASOFT
    fimarith/fcarith tools.

    Notes
    -----
    Most AstroPhotography classes automatically apply any image ``PEDESTAL``
    specified in a FITS header. Unlike those classes, ApImArith **always ignores**
    ``PEDESTAL`` values because they can be incorrecr or misleading on mathematically
    manipulated data, or you may be needing to fix a file that mistakenly still has
    a ``PEDESTAL`` value in it.
    """

    def __init__(self, loglevel: str) -> None:
        """Constructs an ApImArith object. No processing is performed.

        Parameters
        ----------
        loglevel : str
            Logging level to use.
        """

        self._name = "ApImArith"
        self._allowed_ops = ["ADD", "SUB", "MUL", "DIV"]  #: Allowed math operators

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

        logger = logging.getLogger("ApImArith")

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

    def _is_float(self, value_str: Any) -> bool:
        """
        Check if a string can be converted to a floating point number,
        returning true if it can.
        """
        try:
            float(value_str)
            return True
        except ValueError:
            return False

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

    def _sanitize_operation(self, operation: str) -> str:
        """
        Check that the operation string is valid and clean it up if necessary

        :param operation: Input string defining the arithmatic operation
          to perform.
        """
        cleaned_operation = operation.strip().upper()
        if cleaned_operation not in self._allowed_ops:
            err_str = (
                f"Error, input operation {cleaned_operation} is not one"
                f" of the allowed operations: {self._allowed_ops}"
            )
            raise ValueError(err_str)

        return cleaned_operation

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

    def process_files(
        self, inp_img: str, operation: str, value: str | float, out_img: str, units: str
    ) -> None:
        """
        Process the input file with the given value (or image) using
        the specified mathematical operatioon, writing the result to
        the specified output file.

        :param inp_img: Input FITS file containing data. Only the data
          in the primary image extension will be read.
        :param operation: String specifying the mathematical operation
          to perform. Must be one of 'ADD', 'SUB', 'MUL', 'DIV'.
        :param value: A string that can be converted to a floating
          point value or the file name of another FITS file. If
          a single value is specified the same value is applied to
          all pixels in the input image primary image extension. If a
          FITS file is specified a pixel by pixel addition, subtraction,
          multiplication or division is applied using the values in the
          primary image extension of the second file. The array dimensions
          of the second file's primary image extension must match the
          input image.
        :param out_img: Output FITS file to write the result to. This
          will be a copy of the input image file with the primary
          extension header and data modified.
        :param units: An optional string specifying a new value for the
          BUNIT header keyword in the output FITS file.
        """

        operation = self._sanitize_operation(operation)

        remove_pedestal = "never"
        data1, hdr1 = read_fits(inp_img, 0, a_logger=self._logger, remove_pedestal=remove_pedestal)
        data2 = None
        second_data = None
        value_str = None

        if self._is_float(value):
            data2 = float(value)
            second_data = f"scalar {data2}"
            value_str = f"{data2}"
            self._logger.debug(f"Input second item can be cast to float: {data2:.3f}")
        else:
            # Check if a valid file name path.
            if not Path(value).exists():
                err_str = f"Error, {value} is not a scalar or a valid file path."
                raise ValueError(err_str)

            self._logger.debug(f"Input second item is an existing path: {value}")

            # Attempt to read as FITS.
            try:
                data2, hdr2 = read_fits(
                    value, 0, a_logger=self._logger, remove_pedestal=remove_pedestal
                )
                second_data = "array"
                # Check data is the same array dimension.
                if data1.shape != data2.shape:
                    err_str = (
                        f"Error, the dimension of the second data array does not match the first."
                        f" First image shape: {data1.shape}, second image shape: {data2.shape}"
                    )
                    raise RuntimeError(err_str)
            except:
                err_str = f"Error, {value} is not a valid FITS file."
                raise ValueError(err_str)

            value_str = Path(value).name

        # Perform operation.
        result = np.zeros(data1.shape, dtype=data1.dtype)
        if "ADD" in operation:
            np.add(data1, data2, out=result)
            self._logger.info(f"Added {second_data} to input image")
        elif "SUB" in operation:
            np.subtract(data1, data2, out=result)
            self._logger.info(f"Subtracted {second_data} from input image")
        elif "MUL" in operation:
            np.multiply(data1, data2, out=result)
            self._logger.info(f"Multiplied input image by {second_data}")
        elif "DIV" in operation:
            np.divide(data1, data2, out=result)
            self._logger.info(f"Divided input image by {second_data}")

        # String representing operation.
        ohistory_str = (
            f"{self._name} input1 operation input2 are:"
            f" {Path(inp_img).name} {operation} {value_str}"
        )

        # Write data
        self._write_corrected_image(inp_img, 0, out_img, result, units, ohistory_str)
        self._logger.debug("File processing completed.")
        return
