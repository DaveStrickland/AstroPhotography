"""
Contains the implementation of the ApGradientToMask class.
"""

# SPDX-License-Identifier: GPL-3.0-or-later

#  2025-03-15 dks : Initial work on ApGradientToMask.

# import sys
import logging
from pathlib import Path

##from datetime import datetime  # , timezone
##import time
from typing import Any

import numpy as np
import numpy.typing as npt
from astropy.io import fits
from scipy import ndimage

# AstroPhotography includes
from .. import __version__

##from ..util.read import read_fits
from ..util.datautils import img_stats


class ApGradientToMask:
    """
    Creates a binary mask of logically connected pixels associated with
    regions of large pixel value gradient in one or more input FITS
    images or arrays.

    The primary use case of this class within an astrophotography workflow
    is identifying regions in image flat fields affected by dust particles.
    These cause large gradients in the normalized flats that create artifacts
    in calibrated images that require correction.

    The algorithm is as follows:

    .. code-block:: python

    """

    def __init__(self, loglevel: str) -> None:
        """
        Constructs an ApGradientToMask object. No processing is performed.

        Parameters
        ----------
        loglevel : str
            Logging level to use.
        """

        self._name = "ApGradientToMask"

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

        logger = logging.getLogger("ApGradientToMask")

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

    def mask_files(
        self,
        input_files: list[str],
        outhole_mask: str,
        extnum: int | str = 0,
        fwhm_pix: float = 1.5,
        sigma: float = 5.0,
    ) -> None:
        """
        Fix holes by inpainting pixels with values matching the statistical
        properties of the local surroundings, based on a hole mask in an
        FITS file, and write the updated file with a new
        file name.

        See :func:`fix_holes` for the method used.

        Parameters
        ----------
        input_files : str
            A list of input FITS files that will be used to create a combined
            gradient-based mask. The user ensure that all images are of the
            same dimension and size.
        outhole_mask : str
            Output FITS file containing the computed hole mask.
            Will be over-written if it already exists.
        extnum : int or str, optional, default=0
            FITS extension number or name that contains the data in the
            input files.
        fwhm_pix : float, optional, default=1.5
            This is the FWHM of the Guassian gradient filter in units of pixels
            that is applied to a median-normalized copy of the input image(s).
        sigma : float, optional, default = 5.0
            The gradient map is thresholded such that pixels where the gradient
            is greater than sigma times the gradient standard deviation above
            the mean gradient value are defined as holes or the edges of holes.
            (holes are in-filled and then have a binary dilation applied to them.)
        """

        data_list: list[npt.NDArray] = []
        for filename in input_files:
            self._check_file_exists(filename)
            with fits.open(filename) as hdulist:
                data_list.append(hdulist[extnum].data.copy())
                self._logger.info(f"Read data from {filename}")

        hole_mask = self.make_gradient_mask(data_list, fwhm_pix=fwhm_pix, sigma=sigma)

        self._logger.info(f"Will write master hole mask to {outhole_mask}")
        hdu = fits.PrimaryHDU(data=hole_mask)
        hdulist = fits.HDUList([hdu])
        hdulist[0].data = hole_mask
        hdulist[0].header["IMAGETYP"] = ("HOLEMASK", "Good pixels = 0, holes = 1")
        hdulist[0].header["COMMENT"] = (
            "Hole mask from gaussian image gradient thresholding and binary dilation"
        )
        hdulist[0].header["COMMENT"] = f"Generated by {self._name} version {__version__}"
        hdulist.writeto(outhole_mask, output_verify="silentfix", overwrite=True, checksum=True)
        return

    def make_gradient_mask(
        self,
        data_list: list[npt.NDArray],
        fwhm_pix: float = 1.5,
        sigma: float = 5.0,
    ) -> npt.NDArray | None:
        """


        Parameters
        ----------
        data_list : list of ndarray
            A list of input ndarrays that will be processed.
        fwhm_pix : float, optional, default=1.5
            This is the FWHM of the Guassian gradient filter in units of pixels
            that is applied to a median-normalized copy of the input image(s).
        sigma : float, optional, default = 5.0
            The gradient map is thresholded such that pixels where the gradient
            is greater than sigma times the gradient standard deviation above
            the mean gradient value are defined as holes or the edges of holes.
            (holes are in-filled and then have a binary dilation applied to them.)

        Returns
        -------
        merged_mask : ndarray
            The merged mask, combined from the thresholded in-filled gradient
            maps created from each input data array.
        """

        merged_mask = None

        # check that all input ndarray shapes match
        num_inputs = len(data_list)
        if num_inputs == 0:
            err_msg = "Error, data_list must contain one or more NDArrays of the same shape."
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        elif num_inputs > 1:
            default_shape = data_list[0].shape
            for idx in range(1, num_inputs):
                if data_list[idx].shape != default_shape:
                    err_msg = (
                        f"Error, input data #{idx} (0-based) shape {data_list[idx].shape}"
                        " does not math first data shape {default_shape}"
                    )
                    self._logger.error(err_msg)
                    raise RuntimeError(err_msg)

        # loop through input data creating the masks
        mask_list: list[npt.NDArray] = []
        for fdata in data_list:
            a_mask = self._make_gradient_mask(fdata, fwhm_pix=fwhm_pix, sigma=sigma)
            mask_list.append(a_mask)

        # combine the masks
        for a_mask in mask_list:
            if merged_mask is None:
                merged_mask = a_mask.copy()
            else:
                merged_mask = np.logical_or(merged_mask, a_mask).astype(np.int16)

        assert merged_mask is not None
        nbad = np.sum(merged_mask)
        ntot = default_shape[0] * default_shape[1]
        pct = 100.0 * float(nbad) / float(ntot)
        self._logger.info(
            f"Number of bad pixels associated with holes: {nbad}/{ntot} ({pct:.2f}%)"
        )
        return merged_mask

    def _make_gradient_mask(
        self, imgdata: npt.NDArray, fwhm_pix: float, sigma: float
    ) -> npt.NDArray:
        """
        Create a mask consisting of filleds region that enclose where the Gaussian
        gradient of the input data is greater than sigma times the standard
        deviation of the gradient magnitude above the mean gradient magnitude.

        Parameters
        ----------
        imgdata : ndarray
            Input ndarray to be gradient filtered and masked.
        fwhm_pix : float, optional, default=1.5
            This is the FWHM of the Guassian gradient filter in units of pixels
            that is applied to a median-normalized copy of the input image(s).
        sigma : float, optional, default = 5.0
            The gradient map is thresholded such that pixels where the gradient
            is greater than sigma times the gradient standard deviation above
            the mean gradient value are defined as holes or the edges of holes.
            (holes are in-filled and then have a binary dilation applied to them.)

        Returns
        -------
        merged_mask : ndarray
            The merged mask, combined from the thresholded in-filled gradient
            maps created from each input data array.
        """

        # Normalize input data by median
        medval = np.median(imgdata)
        norm_flat_data = imgdata.astype(np.float64) / medval

        # Gaussian FWHM = 2.35482 sigma, so to smooth to FWHM pixels...
        sigma_pix = fwhm_pix / 2.35482
        flat_data_gradient = ndimage.gaussian_gradient_magnitude(norm_flat_data, sigma=sigma_pix)

        [minval, maxval, meanval, stdval, medval] = img_stats(
            flat_data_gradient, "flat_data gradient", False
        )
        thresh_gradient = meanval + sigma * stdval
        self._logger.debug(f"Mean normalized gaussian gradient = {meanval:.3f} +/- {stdval:.3f}")
        self._logger.debug(f"{sigma}-sigma gradient threshold = {thresh_gradient:.6f} ADU")

        # Threshold and fill
        grad_mask = (flat_data_gradient > thresh_gradient).astype(np.int16)
        filled = ndimage.binary_fill_holes(grad_mask).astype(np.int16)

        ##struct1 = ndimage.generate_binary_structure(2, 1)
        struct2 = ndimage.generate_binary_structure(2, 2)

        # TODO? To remove features smaller than the structuring element do an opening.
        # may be useful if we pick up bad pixels that we would otherwise correct
        # in a different way

        final_mask = ndimage.binary_dilation(filled, structure=struct2).astype(np.int16)
        self._logger.debug(f"Total number of pixels in mask: {np.sum(final_mask)}")
        return final_mask
