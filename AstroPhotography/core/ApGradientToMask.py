"""
Contains the implementation of the ApGradientToMask class.
"""

# SPDX-License-Identifier: GPL-3.0-or-later

#  2025-03-15 dks : Initial work on ApGradientToMask.

# import sys
import logging
from pathlib import Path
from datetime import datetime  # , timezone
import time
from typing import Any

import numpy as np
import numpy.typing as npt
from astropy.io import fits

# AstroPhotography includes
from .. import __version__
from ..util.read import read_fits
from ..util.datautils import regionprops_to_astropy_table


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

        self._logger.info(f"Will write master hold mask to {outhole_mask}")
        hdu = fits.PrimaryHDU(data=hole_mask)
        hdulist = fits.HDUList([hdu])
        hdulist[0].data = hole_mask
        hdulist[0].header["IMAGETYP"] = ("HOLEMASK", "Good pixels = 0, holes = 1")
        hdulist[0].header["COMMENT"] = (
            "Hole mask from gaussian image gradient thresholding and binary dilation"
        )
        hdulist.writeto(outhole_mask, output_verify="silentfix", overwrite=True, checksum=True)
        return

    def make_gradient_mask(
        self,
        data_list: list[npt.NDArray],
        fwhm_pix: float = 1.5,
        sigma: float = 5.0,
    ) -> npt.NDArray:
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

        # check that all input ndarray shapes match

        # loop through input data creating the masks

        # combine the masks

        # print some info?

        hfillmin = int(hfillmin)
        hfilltgt = float(hfilltgt)
        deltapix = 1  # number of pixels on each size that a bounding box is grown

        # Info message about array shape and data type.
        self._logger.info(
            f"fix_holes: data has {data.shape[0]} rows x {data.shape[1]} columns"
            f", dtype={data.dtype}"
        )
        self._logger.info(
            f"fix_holes: mask has {holemask.shape[0]} rows x {holemask.shape[1]} columns"
            f", dtype={holemask.dtype}"
        )

        # Check that sizes match
        if data.shape != holemask.shape:
            msg = (
                f"Error, the shape of the input data array ({data.shape})"
                f" does not match that of the bad pixel mask array ({holemask.shape})."
            )
            self._logger.error(msg)
            raise RuntimeError(msg)

        if not np.issubdtype(data.dtype, np.floating):
            # Not a floating point datatype.
            msg = (
                "Pixel medians may suffer from casting truncation because"
                f" the input data is not a floating point datatype ({data.dtype})."
            )
            self._logger.warning(msg)

        newdata = data.copy()
        mask = holemask != ApFixHoles.MASK_GOOD
        npix = data.size
        nbad = np.sum(mask)
        pctbad = 100.0 * nbad / npix

        # Identify unique holes using segmentation
        label_mask = skm.label(mask)
        properties = ["label", "area", "area_bbox", "bbox"]
        raw_holes_tbl = regionprops_to_astropy_table(
            label_mask, intensity_image=None, properties=properties
        )

        # Compute some of the statistics
        num_holes = len(raw_holes_tbl)
        self._logger.debug(f"Number of holes in mask: {num_holes}")
        self._logger.debug(
            f"Percentage of pixels considered bad: {pctbad:.3f} ({nbad:d}/{npix:d})"
        )
        fixed_stats = {
            "HOLENBAD": (num_holes, " Number of holes in the input mask"),
            "HOLENPIX": (nbad, " Number of bad pixels in the hole mask"),
            "HFILLMIN": (hfillmin, "Minimum number of good pixels sampled around each hole"),
            "HFILLTGT": (hfilltgt, "Target ratio of good pixels to hole pixels"),
        }

        # Formatting for diagnostic output
        nprint = 40
        hole_num_fmt = "8d"
        pix_fmt = "9d"
        val_fmt = "10.2f"
        if np.nanmean(data) < 1.0:
            val_fmt = "10.6f"
        msg = f"Diagnostics of the first {nprint} holes follow:\n"
        hdr = (
            "{0:8s}  {1:9s}  {2:9s}  {3:9s}  {4:9s}  {5:9s}"
            "  {6:10s}  {7:10s}  {8:9s}  {9:10s}  {10:10s} \n"
        ).format(
            "hole_num",
            "bbox_rmin",
            "bbox_rmax",
            "bbox_cmin",
            "bbox_cmax",
            "hole_npix",
            "hole_medn",
            "hole_std",
            "npix_good",
            "fix_medn",
            "fix_std",
        )
        msg += hdr
        line = (
            "{:<"
            f"{hole_num_fmt}"
            "}  "  # hole_num
            "{:>"
            f"{pix_fmt}"
            "}  "  # rmin
            "{:>"
            f"{pix_fmt}"
            "}  "  # rmax
            "{:>"
            f"{pix_fmt}"
            "}  "  # cmin
            "{:>"
            f"{pix_fmt}"
            "}  "  # cmax
            "{:>"
            f"{pix_fmt}"
            "}"  # hole_npix
            "{:>"
            f"{val_fmt}"
            "}  "  # hole_medn
            "{:>"
            f"{val_fmt}"
            "}  "  # hole_std
            "{:>"
            f"{pix_fmt}"
            "}"  # npix_good
            "{:>"
            f"{val_fmt}"
            "}  "  # fix_medn
            "{:>"
            f"{val_fmt}"
            "}  "  # fix_std
            "\n"
        )
        ##print(line)

        # Iterate over holes.
        # - Construct row and column index arrays, then select only the
        #   indices of the bad pixels and store those as 1-dimensional arrays.
        nrows = mask.shape[0]
        ncols = mask.shape[1]
        perf_time_start = time.perf_counter()  # highest res timer, counts sleeps

        ##row_idxs, col_idxs = np.mgrid[0:nrows, 0:ncols] # TODO cleanup?
        ##bp_row_idxs = row_idxs[mask].ravel() # TODO cleanup?
        ##bp_col_idxs = col_idxs[mask].ravel() # TODO cleanup?

        for idx in range(num_holes):
            label = raw_holes_tbl["label"][idx]
            rmin = raw_holes_tbl["bbox_rmin"][idx]
            rmax = raw_holes_tbl["bbox_rmax"][idx]
            cmin = raw_holes_tbl["bbox_cmin"][idx]
            cmax = raw_holes_tbl["bbox_cmax"][idx]

            # compute stats on number of hole pixels and background pixels
            mask_co = mask[rmin:rmax, cmin:cmax]  # True where bad
            nbad_co: int = np.sum(mask_co)
            init_nbad_co = nbad_co  # we can grow into other holes
            init_bbox_size: int = mask_co.size  # so we store some original stats
            bbox_size = init_bbox_size
            ngood_co: int = mask_co.size - nbad_co
            goodbad_ratio: float = float(ngood_co) / float(nbad_co)

            # iteratively adjust bounding box if it fails to meet the
            # hfillmin and hfilltgt constraints, but break out of loop if the
            # bounding box is not actually getting bigger.
            num_growth = 0
            while (ngood_co < hfillmin) or (goodbad_ratio < hfilltgt):
                last_bbox_size = bbox_size

                # grow by deltapix on each side
                rmin = max(0, rmin - deltapix)  # included
                rmax = min(nrows, rmax + deltapix + 1)  # excluded
                cmin = max(0, cmin - deltapix)  # included
                cmax = min(ncols, cmax + deltapix + 1)  # excluded

                mask_co = mask[rmin:rmax, cmin:cmax]  # True where bad
                nbad_co = np.sum(mask_co)
                bbox_size = mask_co.size
                ngood_co = mask_co.size - nbad_co
                goodbad_ratio = float(ngood_co) / float(nbad_co)
                num_growth += 1

                if bbox_size == last_bbox_size:
                    # bounding box has not grown
                    break

            info_str = (
                f"Hole {label}, {init_nbad_co} bad pix"
                f", initial bounding box {init_bbox_size} pix"
                f", final bounding box {bbox_size} pix"
                f" after {num_growth} rounds of growth"
                f", has {ngood_co} good pixels (good to bad pixel"
                f" ratio={goodbad_ratio:.2f})"
            )
            self._logger.debug(info_str)

            data_slice = data[rmin:rmax, cmin:cmax]
            mask_slice = mask[rmin:rmax, cmin:cmax]
            fixer = ApGaussianVariateFiller(self._loglevel, data_slice, mask_slice)
            fixed_slice = fixer.fill_with_median_stddev()
            newdata[rmin:rmax, cmin:cmax] = fixed_slice

            # Before and after stats
            hole_medval = np.median(data_slice[mask_slice])
            hole_stdval = np.std(data_slice[mask_slice])
            fixed_medval = np.median(fixed_slice[mask_slice])
            fixed_stdval = np.std(fixed_slice[mask_slice])

            # Diagnostic string
            if idx < nprint:
                msg += line.format(
                    label,
                    rmin,
                    rmax,
                    cmin,
                    cmax,
                    init_nbad_co,
                    hole_medval,
                    hole_stdval,
                    ngood_co,
                    fixed_medval,
                    fixed_stdval,
                )

        perf_time_end = time.perf_counter()  # highest res timer, counts sleeps
        if verbose:
            print(msg)

        # Measure performance, in case we need to speed this up later.
        run_time = perf_time_end - perf_time_start
        ms_per_pix = 1000 * run_time / num_holes
        msg = f"Processed {num_holes} pixels in {run_time:.3f} s, {ms_per_pix:.3f} ms per hole."
        self._logger.info(msg)

        fixed_stats["HOLECORR"] = (True, "True if any holes were corrected")
        fixed_stats["HOLENFIX"] = (nbad, "Number of bad pixels corrected")
        return newdata, fixed_stats
