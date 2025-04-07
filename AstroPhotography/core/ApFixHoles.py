"""Contains the implementation of the ApFixHoles class."""

# SPDX-License-Identifier: GPL-3.0-or-later

#  2024-12-26 dks : Original based on ApFixBadPixels.

import sys
import logging
from pathlib import Path
import time
from datetime import datetime
from typing import Any

import numpy as np
import numpy.typing as npt
import skimage.measure as skm
from skimage import restoration
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

# AstroPhotography includes
from .. import __version__
import AstroPhotography.util as util


class ApFixHoles:
    """
    A class used to fix pre-indentified holes or regions of bad data within an image
    represented by a mask with random values sharing the same median and standard
    deviation as pixels surrounding each hole.

    This class is intended for use with larger and/or more irregular regions of
    bad pixels than ApFixBadPixels, which operates on a pixel-by-pixel basis.
    Instead ApFixHoles uses image segmentation to address connected regions of
    bad pixels as a single region to be patched.

    The holes to be filled are specified by an input mask array (or FITS file)
    with integer or local values, and where non-zero (or True) denotes that the
    pixel is part of a hole that should be filled. Zero (or False) denotes that
    the pixel is good, and is not to be altered.
    """

    # Class level constants
    MASK_GOOD = 0
    MASK_BAD = 1

    def __init__(self, loglevel: str) -> None:
        """
        Constructs an ApFixHoles object. No processing is performed.

        Parameters
        ----------
        loglevel : str
            Logging level to use, e.g. 'INFO'
        """

        self._name = "ApFixHoles"

        # Initialize logging
        self._loglevel = loglevel
        self._initialize_logger(self._loglevel)
        self._logger.debug(f"{self._name} instance constructed.")
        return

    def _check_file_exists(self, filename: str) -> None:
        """
        Raises an exception if the specified file does not exist

        Parameters
        ----------
        filename : str
            File name to check the existence of

        Raises
        ------
        RuntimeError :
            If filename does not exist
        """

        if not Path(filename).exists():
            err_msg = f"Cannot find {filename}. Not a valid path or file."
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return

    def _initialize_logger(self, loglevel: str) -> None:
        """
        Initialize and return the logger

        Parameters
        ----------
        loglevel : str
            Logging level to use, e.g. 'INFO'
        """

        self._logger = logging.getLogger(self._name)

        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError("Invalid log level: {}".format(loglevel))
        self._logger.setLevel(numeric_level)

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

        # Used in cases where we get the same message twice or more
        # See https://stackoverflow.com/a/44426266
        self._logger.propagate = False
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
        if "BZERO" in ext_hdr:
            bzero = ext_hdr["BZERO"]
        else:
            bzero = 0
        if "BSCALE" in ext_hdr:
            bscale = ext_hdr["BSCALE"]
        else:
            bscale = 1.0
        info_str = (
            f"{ndim}-D BITPIX={bitpix} image with {cols} columns, {rows} rows,"
            f" BSCALE={bscale}, BZERO={bzero}"
        )

        if ndim == 3:
            layers = ext_hdr["NAXIS3"]
            info_str = (
                f"{ndim}-D BITPIX={bitpix} image with {cols} columns, {rows} rows, {layers} layers"
                f" BSCALE={bscale}, BZERO={bzero}"
            )

        self._logger.debug(info_str)
        if ndim == 3:
            self._logger.error("Error, 3-D handling has not been implemented yet.")
            sys.exit(1)

        # Get data absolute limits.
        minval = np.amin(ext_data)
        maxval = np.amax(ext_data)
        medval = np.median(ext_data)
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
                    f"After PEDESTAL removal, min={minval:.2f}, max={maxval:.2f},"
                    f" median={medval:.2f}"
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
        Writes the bad-pixel-corrected data to the specified output
        file, preserving all other items from the original input file.

        The output file differs from the original input data file in
        having the bad-pixel corrected image, and additional the header
        keywords specified in :func:`fix_holes`

        Parameters
        ----------
        inpdata_file : str
            Input FITS data file affected by bad pixels.
        ext_num : int or str
            Extension number or name for data array and header.
        outdata_file : str
            File name for the modified copy of the input data file where the
            bad pixels have had their data value inpainted.
        odata : ndarray
            Hole corrected data array.
        odict : dict
            Input dictionary generated by fix_holes (and
            fix_files), used to update output header as described above.
            Elements of the dictionary that have keywords containing
            'HOLE' or `HFILL` are written to the output fits file, any other
            elements are not.
        """

        self._logger.debug(f"Hole pixel keywords added to output: {odict}")

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

            # Modify data
            hdu_list[ext_num].data = odata

            # Modify header
            for kw, val in odict.items():
                if "HOLE" in kw or "HFILL" in kw:
                    hdu_list[ext_num].header[kw] = val

            tnow = datetime.now().isoformat(timespec="milliseconds")
            hdu_list[ext_num].header["HISTORY"] = f"Applied {self._name} {__version__} at {tnow}"

            # Write to new file
            hdu_list.writeto(outdata_file, output_verify="ignore", overwrite=True)

        self._logger.info(f"Wrote bad pixel corrected file to {outdata_file}")
        return

    def fix_files(
        self,
        inpdata_file: str,
        holemask_file: str,
        outdata_file: str,
        hfillmin: int = 20,
        hfilltgt: float = 2.0,
        verbose: bool = True,
    ) -> None:
        """
        Fix holes by inpainting pixels with values matching the statistical
        properties of the local surroundings, based on a hole mask in an
        FITS file, and write the updated file with a new
        file name.

        See :func:`fix_holes` for the method used.

        Parameters
        ----------
        inpdata_file : str
            Input FITS data file affected by bad pixels.
        holemask_file : str
            Input FITS file specifying the bad pixels.
            Bad pixels should have non-zero values, good pixels should have
            zero as the pixel value.
        outdata_file : str
            Modified copy of the input data file where the
            bad pixels have had their data valus modified by the median
            of the surrounding good pixels.
        hfillmin : int, optional, default=20
            This is the minimum number of good pixels around each hole to use
            when determining the statistical properties of the local good pixels.
            The background region around each hole will be adjusted in size to
            make sure there are at least ``hfillmin`` good (not masked) pixels.
        hfilltgt : floart, optional, default = 2.0
            This is target ratio of good pixels to hole (masked) pixels to aim
            for when determining the size of the background region around each
            hole. Each background region has at minimum this ratio of good
            to mask pixels, but can have more.
        verbose : bool, optional, default=True
            If True a table of information on each hole processed will be
            output to stdout.
        """

        ext_num = 0  # Really should have a better way of setting this.

        msg = (
            f"fix_files input data file={inpdata_file},"
            f" hole mask file={holemask_file},"
            f" output file={outdata_file}, hfillmin={hfillmin}, hfilltgt={hfilltgt:.2f}"
        )
        self._logger.info(msg)

        idata, ihdr = self._read_fits(inpdata_file, ext_num)
        mskdata, mskhdr = self._read_fits(holemask_file, ext_num)

        odata, odict = self.fix_holes(idata, mskdata, hfillmin, hfilltgt, verbose)
        odict["HOLEFILE"] = (Path(holemask_file).name, "Name of hole mask file used")

        # Create copy of input image and write modified data to it
        # with an updated header
        self._write_corrected_image(inpdata_file, ext_num, outdata_file, odata, odict)
        return

    def fix_holes(
        self,
        data: npt.NDArray,
        holemask: npt.NDArray,
        hfillmin: int = 20,
        hfilltgt: float = 2.0,
        verbose: bool = True,
    ) -> tuple[npt.NDArray, dict[str, Any]]:
        """
        Fix pixels associated with holes in the input array, based on the mask
        hole mask array, returning the modified data array,
        and a dictionary summarizing the number of affected holes, total number
        of affected pixels and so on.

        For each separate hole in the hold mask a surrounding retangular region
        is determined. A separate hole is a set of connected non-zero pixels in
        the hole mask that the ``scipy.ndimage`` segmentation functions identify.
        The statistical median and standard deviation of good
        pixels within each "background" region associated with each hold is
        determined. Then the pixels in the hole are randomly assigned a value
        based on a gaussian normal variate with a mean equal to the background
        median and standard deviation equal to the background standard deviation.

        The size of each background region is chosen so that the number of good
        (non hole mask) pixels within it is greater than or equal to hfillmin,
        and the number of good pixels divided by the number of hole pixels is
        greater than or equal to hfilltgt.

        The output dictionary can be used to modify the FITS header of
        an output file, and consists of keyword: (value, comment) pairs.

        - HOLECORR: Logical true denoting whether bad pixel correction
          applied.
        - HOLEFILE: The name of the hole mask file used, stripped of any
          preceding path elements. (This is only specified if ``fix_files``
          was used.)
        - HFILLMIN: Minimum number of good pixels to sample to get local median
          and standard deviation estimate for statistical inpainting
          of a hole. This parameter and HFILLTGT control the size
          of the "good" data background region used for each hole.
        - HFILLTGT: The minimum allowed ratio of good pixels to bad (hole) pixels
          to aim for when creating a background region around each hole.
          This parameter and HFILLMIN control the size
          of the "good" data background region used for each hole.
        - HOLENBAD: Number of separate holes (logically connected features)
        - HOLENPIX: Number of bad pixels in the hole mask file.
        - HOLENFIX: Number of pixels successfully corrected.

        Parameters
        ----------
        data : ndarray
            Input ndarray representing the data affected by holes.
        holemask : ndarray
            Input ndarray specifying the which pixels are considered holes.
            Hole pixels should have non-zero values, good pixels should have
            zero as the pixel value.
        hfillmin : int, optional, default=20
            This is the minimum number of good pixels around each hole to use
            when determining the statistical properties of the local good pixels.
            The background region around each hole will be adjusted in size to
            make sure there are at least ``hfillmin`` good (not masked) pixels.
        hfilltgt : float, optional, default = 2.0
            This is target ratio of good pixels to hole (masked) pixels to aim
            for when determining the size of the background region around each
            hole. Each background region has at minimum this ratio of good
            to mask pixels, but can have more.
        verbose : bool, optional, default=True
            If True a table of information on each hole processed will be
            output to stdout.

        Returns
        -------
        newdata : ndarray
            A modified copy of the input data array with the pixel values in the
            holes inpainted with values consistent with the local background
            statistics.
        fixed_stats : dict
            A dictionary of metadata associlated with the number of holes and
            pixels corrected.
        """

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
        raw_holes_tbl = util.regionprops_to_astropy_table(
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


class ApBiharmonicFiller:
    """
    Encapulsates the functionality to fill a masked region using scikit image's
    ``skimage.restoration.inpaint_biharmonic`` function.

    :func:`fill_biharmonic` in-paints the masked pixels in the image cutout using
    the values of surrounding "good" (non-masked) pixels. The masked region should
    be surrounded by at least 5 (?) good pixels on each side.

    Notes
    -----
    Biharmonic inpainting works on greyscale images with values in the range [0,1].
    This class handles scaling the input data to that range and restoring the
    original data range on output. **However** the quality of inpainting also
    depends on the expected data values in the region to be filled not being too
    close the minimum or maximum pixel values in the image.

    It is important to use a image slice that is local to the region to be
    inpainted, rather than attempting to inpaint all the holes on a large image.
    That does work on a normal computer graphics images, but will not work well
    on astronomical data.

    For example, if the inpainted region would be expected to have data values
    around 100 and the input image slice has a data range of [-110, 350] then
    Biharmonic inpainting will work well because the values to be inpainted are
    not near the extremes. If the input image slice had a few pixels with very
    large values, e.g. 40,000, then the inpainiting will not work well because
    the 40,000 will be scaled to 1 internally, while the expected region's value
    of 100 will be scaled to approximately 0.005 and may get inpainted with
    zeros that then scale back to the input image minimum.

    Warnings
    --------
    Work In Progress, not yet fully tested
    """

    def __init__(self, loglevel: str, data_cutout: npt.NDArray, mask_cutout: npt.NDArray) -> None:
        """

        Parameters
        ----------
        loglevel : str
            Standard logging level string such as 'debug' or 'info'
        data_cutout : ndarray
            A slice of a larger data-containing array that has a hole in those
            pixels where mask_cutout is non-zero or True. A slice is recommended
            as the statisticaly properties of the entire slice are computed.
        mask_cutout : ndarray
            A slice of a mask array (integer 0 or 1, or logical True or False)
            where 1/True represents pixels considered part of the hole to be filled
            and 0/False pixels are those considered to be good (whose statistical
            properties will be used to fill the hole).
        """

        self._name = "ApBiharmonicFiller"
        self._loglevel = loglevel
        self._initialize_logger(loglevel)
        self._inpdata = data_cutout
        self._inpmask = mask_cutout
        return

    def _initialize_logger(self, loglevel: str) -> None:
        """
        Initialize and return the logger

        Parameters
        ----------
        loglevel : str
            Logging level to use, e.g. 'INFO'
        """

        self._logger = logging.getLogger(self._name)

        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError("Invalid log level: {}".format(loglevel))
        self._logger.setLevel(numeric_level)

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

        # Used in cases where we get the same message twice or more
        # See https://stackoverflow.com/a/44426266
        self._logger.propagate = False
        return

    def fill_biharmonic(self) -> npt.NDArray:
        """
        Fill the hole using biharmonic equations, from the skimage.restoration module.

        Returns
        -------
        odata : ndarray
            Copy of input image where pixels within the mask have been inpainted
            using biharmonic equations
        """
        self._logger.warning("fill_biharmonic not yet fully implemented or tested")
        odata = self._inpaint_biharmonic_wrapper(self._inpdata, self._inpmask)
        return odata

    def _inpaint_biharmonic_wrapper(self, image: npt.NDArray, mask: npt.NDArray) -> npt.NDArray:
        """
        Wrapper around skimage.restoration.inpaint_biharmonic that handles the
        rescaling to and from [0,1] range it uses

        Parameters
        ----------
        image :  ndarray
            Input single-channel floating point image
        mask : ndarray
            Input logical or integer mask array where pixels with 0 or False values
            are considered good, while pixels with 1 or True are to be painted
            over.

        Returns
        -------
        odata : ndarray
            Copy of input image where pixels within the mask have been inpainted
            using biharmonic equations

        """

        # input image has range minval to maxval, so to scale to [0,1] range
        # we do
        # scaled = (image - minval)/(maxval-minval)
        #
        # to restore the fixed_scaled we do
        # fixed = (maxval-minval) * fixed_scaled + minval

        minval = np.min(image)
        maxval = np.max(image)
        delta = maxval - minval
        self._logger.debug(
            f"Input image has minimum={minval:f}, maximum={maxval:f}"
            f", range={delta:f} counts per pixel."
        )

        scaled = (image.astype(float) - minval) / delta
        chk_minval = np.min(scaled)
        chk_maxval = np.max(scaled)
        self._logger.debug(
            f"Rescaled image has minimum={chk_minval:f}, maximum={chk_maxval:f}"
            f", range={chk_maxval-chk_minval:f} counts per pixel."
        )

        fixed_scaled = restoration.inpaint_biharmonic(scaled, mask)

        fixed = (fixed_scaled * delta) + minval
        chk_minval = np.min(fixed)
        chk_maxval = np.max(fixed)
        self._logger.debug(
            f"Restored image has minimum={chk_minval:f}, maximum={chk_maxval:f}"
            f", range={chk_maxval-chk_minval:f} counts per pixel."
        )
        return fixed


class ApGaussianVariateFiller:
    """
    Encapulsates the functionality to fill a masked region with pixel values
    drawn from a gaussian distribution based on one of several statistical
    measurements of the surrounding "background" (i.e. not masked pixels)

    :func:`fill_with_median_stddev` generates random numbers with a mean equal
    to the median of the surrounding "background", and a standard deviation
    from the same background. This function does not attempt to remove outliers.

    :func:`fill_with_sigma_clipped_mean_stddev` generates random numbers with
    a mean and standard deviation based on the default astropy sigma-clipped
    stats implementation.
    """

    def __init__(self, loglevel: str, data_cutout: npt.NDArray, mask_cutout: npt.NDArray) -> None:
        """

        Parameters
        ----------
        loglevel : str
            Standard logging level string such as 'debug' or 'info'
        data_cutout : ndarray
            A slice of a larger data-containing array that has a hole in those
            pixels where mask_cutout is non-zero or True. A slice is recommended
            as the statisticaly properties of the entire slice are computed.
        mask_cutout : ndarray
            A slice of a mask array (integer 0 or 1, or logical True or False)
            where 1/True represents pixels considered part of the hole to be filled
            and 0/False pixels are those considered to be good (whose statistical
            properties will be used to fill the hole).
        """

        self._name = "ApGaussianVariateFiller"
        self._loglevel = loglevel
        self._initialize_logger(loglevel)
        self._inpdata = data_cutout
        self._inpmask = mask_cutout
        self._antimask = np.logical_not(mask_cutout)
        return

    def _initialize_logger(self, loglevel: str) -> None:
        """
        Initialize and return the logger

        Parameters
        ----------
        loglevel : str
            Logging level to use, e.g. 'INFO'
        """

        self._logger = logging.getLogger(self._name)

        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError("Invalid log level: {}".format(loglevel))
        self._logger.setLevel(numeric_level)

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

        # Used in cases where we get the same message twice or more
        # See https://stackoverflow.com/a/44426266
        self._logger.propagate = False
        return

    def fill_with_median_stddev(self) -> npt.NDArray:
        """
        Fill the hole with pixel values drawn from a gaussian distribution that
        has a mean equal to the median of the non-masked pixels and a standard
        devaiation equal to the standard deviation of the non-masked pixels.

        Returns
        -------
        odata : ndarray
            Copy of input image where pixels within the mask have been filled
            with random numbers consistent with the surrounding, non-masked,
            part of the image.
        """

        odata = self._inpdata.copy()
        hole_medval = np.median(self._inpdata[self._inpmask])
        hole_stdval = np.std(self._inpdata[self._inpmask])
        if hole_medval < 1.0:
            ffmt = ".6f"
        else:
            ffmt = ".2f"
        self._logger.debug(
            f"Input hole has median={hole_medval:{ffmt}}"
            f" +/- {hole_stdval:{ffmt}} counts per pixel."
        )

        # Calculate mean and std deviation from data not in the mask
        medval = np.median(self._inpdata[self._antimask])
        stdval = np.std(self._inpdata[self._antimask])
        if medval < 1.0:
            ffmt = ".6f"
        else:
            ffmt = ".2f"
        self._logger.debug(
            f"Surrounding good data has median={medval:{ffmt}}"
            f" +/- {stdval:{ffmt}} counts per pixel."
        )

        # generate a full input data sized array of random numbers
        rng = np.random.default_rng()
        randvals = rng.normal(medval, stdval, self._inpdata.size).reshape(self._inpdata.shape)

        odata[self._inpmask] = randvals[self._inpmask]
        hole_medval = np.median(self._inpdata[self._inpmask])
        hole_stdval = np.std(self._inpdata[self._inpmask])
        self._logger.debug(
            f"Filled hole now has median={hole_medval:{ffmt}}"
            f" +/- {hole_stdval:{ffmt}} counts per pixel."
        )
        return odata

    def fill_with_sigma_clipped_mean_stddev(self) -> npt.NDArray:
        """
        Fill the hole with pixel values drawn from a gaussian distribution that
        has a mean and standard deviation equal to the sigma-clipped statistics
        of the non-masked pixels.

        Returns
        -------
        odata : ndarray
            Copy of input image where pixels within the mask have been filled
            with random numbers consistent with the surrounding, non-masked,
            part of the image.
        """

        odata = self._inpdata.copy()
        hole_meanval = np.mean(self._inpdata[self._inpmask])
        hole_medval = np.median(self._inpdata[self._inpmask])
        hole_stdval = np.std(self._inpdata[self._inpmask])
        if hole_medval < 1.0:
            ffmt = ".6f"
        else:
            ffmt = ".2f"
        self._logger.debug(
            f"Input hole has median={hole_medval:{ffmt}}"
            f", mean={hole_meanval:{ffmt}}"
            f" +/- {hole_stdval:{ffmt}} counts per pixel."
        )

        # Calculate mean and std deviation from data not in the mask
        (meanval, medval, stdval) = sigma_clipped_stats(self._inpdata, self._inpmask.astype(bool))
        if medval < 1.0:
            ffmt = ".6f"
        else:
            ffmt = ".2f"
        self._logger.debug(
            f"Surrounding good data has sigma-clipped stats of median={medval:{ffmt}}"
            f", mean={meanval:{ffmt}}"
            f" +/- {stdval:{ffmt}} counts per pixel."
        )

        # generate a full input data sized array of random numbers
        rng = np.random.default_rng()
        randvals = rng.normal(meanval, stdval, self._inpdata.size).reshape(self._inpdata.shape)

        odata[self._inpmask] = randvals[self._inpmask]
        hole_meanval = np.mean(self._inpdata[self._inpmask])
        hole_medval = np.median(self._inpdata[self._inpmask])
        hole_stdval = np.std(self._inpdata[self._inpmask])
        self._logger.debug(
            f"Filled hole now has median={hole_medval:{ffmt}}"
            f", mean={hole_meanval:{ffmt}}"
            f" +/- {hole_stdval:{ffmt}} counts per pixel."
        )
        return odata
