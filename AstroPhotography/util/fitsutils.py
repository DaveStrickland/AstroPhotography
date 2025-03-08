# -*- coding: utf-8 -*-
#
#  Contains the implementation of various Astrophotography utilities
#
#  Copyright 2024 Dave Strickland <dave.strickland@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

# 2024-04-07 dks : issue-002 Initial coding
# 2025-03-01 dks : Split off from ApUtil.py

##import logging
##import os.path
from typing import Any
import numpy as np
import numpy.typing as npt


##from matplotlib.patches import Ellipse
import math
##import time
##from datetime import datetime

from astropy import wcs

##import astropy
from astropy.io import fits

##from astropy.table import QTable, Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units

# AstroPhotography includes
##from .. import __version__


def load_wcs_from_file(filename, extnum=0, verbose=False):
    """
    Load the astropy.wcs.WCS object from a given extension in a FITS file,
    returning it along with some summary statistics

    Modified WCS example based on the `astropy documentation <https://docs.astropy.org/en/stable/wcs/loading_from_fits.html>`_.

    Note
    ~~~~

    * Assumes a 2-dimensional image with angular coordinate axes

    Parameters
    ----------
    filename : str
        Name of the input FITS image.
    extnum : int or str, default=0
        Extension number or name for the extension holding the image data.
        Usually this is 0, for the ``PrimaryHDU``.
    verbose : bool, optional, default=False
        If True then print the axis plate scale, pixel area, number of
        pixels and image angular scale to stdout.

    Returns
    -------
    w : astropy.wcs.WCS
        WCS instance initialized from the FITS header
    pix_scales : array_like
        Floating point array of X and Y axis pixel scales in the default
        angular units of the image (in almost all cases both ``CUNIT1``
        and ``CUNIT2`` is degrees).
    pix_area : float
        Angular area of a pixel, in units of ``CUNIT1 * CUNIT2``.
    im_scales : The angular extent of the images along the X and Y axes,
        in the default angular units of the image.
    """
    w = None
    pix_scales = None  # in deg
    pix_area = None  # in deg^2
    im_scales = None  # in deg

    # Load the FITS hdulist using astropy.io.fits
    with fits.open(filename) as hdulist:
        # Parse the WCS keywords in the primary HDU
        w = wcs.WCS(hdulist[0].header)

        # Print out the "name" of the WCS, as defined in the FITS header
        if verbose:
            print(w.wcs.name)

            # Print out all of the settings that were parsed from the header
            w.wcs.print_contents()

        # Pixel scale at the CRPIX pixel location. Note, ideal, ignores distortions
        pix_scales = wcs.utils.proj_plane_pixel_scales(w)

        # Pixel area at the CRPIX pixel location. Again, ideal, ignores distortions
        pix_area = wcs.utils.proj_plane_pixel_area(w.celestial)

        if w.wcs.naxis != 2:
            print(
                f"WARNING, WCS from {filename} is {w.wcs.naxis}-dimensional, not 2-D as expected"
            )

        if len(pix_scales) != len(w.array_shape):
            raise RuntimeError(
                (
                    f"Shape of pixel scales ({len(pix_scales)})"
                    f" differs from shape of array ({len(w.array_shape)})"
                )
            )
        else:
            im_scales = pix_scales.copy()
            for idx in range(len(w.array_shape)):
                im_scales[idx] *= w.array_shape[idx]

        if verbose:
            print(f"Pixel angular scale ({w.wcs.cunit[0]}):    {pix_scales}")
            print(f"Pixel angular area ({w.wcs.cunit[0]}^2):   {pix_area}")
            print(f"NAXIS1={w.array_shape[0]}   NAXIS2={w.array_shape[1]}")
            print(f"Image angular scale ({w.wcs.cunit[0]}): {im_scales}")

    return w, pix_scales, pix_area, im_scales


def summarize_wcs(w: Any, verbose: bool = True) -> tuple[str, float, float, float, float, float]:
    """
    Given an astropy WCS object, summarise the properties of the WCS header.

    If verbose is True a full summary is printed to STDOUT. A smaller subset of
    key parameters is returned to the caller.

    This assumes that the images are two dimensional, and that angles are in degrees.

    Parameters
    ----------
    w : astropy.wcs.WCS
        A valid astropy WCS instance
    verbose : bool, optional, default=True
        If true, output information to STDOUT.

    Returns
    -------
    dothead_str : str
        Multiline string WCS summary in the text form that ``swarp`` is
        supposed to accept (but appears not to).
    cdelt1_as : float
        Equivalent X-axis pixel size in arcseconds.
    cdelt2_as : float
        Equivalent Y-axis pixel size in arcseconds.
    ximgsz_am : float
        Equivalent X-axis image size in arcminutes.
    yimgsz_am : float
        Equivalent Y-axis image size in arcminutes.
    crot : float
        Image rotation angle (degrees) equiavelnt to CROTA2, the angle between
        the Y-axis and true North.
    """

    ipwcs = w

    dothead_list = []

    # Numbr of dimensions
    val_str = f"{'NAXIS':8s} = {ipwcs.wcs.naxis}"
    dothead_list.append(val_str)

    keywords = ["NAXIS", "CTYPE", "CRVAL", "CRPIX"]
    values = [ipwcs.array_shape, ipwcs.wcs.ctype, ipwcs.wcs.crval, ipwcs.wcs.crpix]
    for keyword, value in zip(keywords, values):
        for idx in range(ipwcs.naxis):
            kw_str = f"{keyword}{1+idx}"
            if "CTYPE" in keyword:
                # Wrap strings in quotes
                val_str = f"{kw_str:8s} = '{value[idx]}'"
            else:
                val_str = f"{kw_str:8s} = {value[idx]}"
            dothead_list.append(val_str)

    naxis = ipwcs.wcs.naxis
    naxis1 = ipwcs.array_shape[0]
    naxis2 = ipwcs.array_shape[1]
    if naxis == 3:
        ##naxis3 = ipwcs.array_shape[2]
        print(f"Warning: summarize_wcs() written for 2-D images, not {naxis}-dimensional data")

    cdelt1 = None
    cdelt2 = None
    crot = None

    if hasattr(ipwcs.wcs, "cd"):
        for irow in range(ipwcs.naxis):
            for jcol in range(ipwcs.naxis):
                kw_str = f"CD{irow+1}_{jcol+1}"
                val_str = f"{kw_str:8s} = {ipwcs.wcs.cd[irow, jcol]}"
                dothead_list.append(val_str)
        cd = ipwcs.wcs.cd
        # From https://lweb.cfa.harvard.edu/~jzhao/SMA-FITS-CASA/docs/wcs88.pdf
        cd11 = cd[0, 0]
        cd12 = cd[0, 1]
        cd21 = cd[1, 0]
        cd22 = cd[1, 1]
        cdelt1_mag = math.sqrt(cd12 * cd12 + cd22 * cd22)
        cdelt2_mag = math.sqrt(cd12 * cd12 + cd22 * cd22)
        # the sign of cdelt1.cdelt2 = sign of (cd11*cd22 - cd12*cd21)
        # if the RHS is negative cdelt1 is negative by convention, so cdelt2 is always positive
        tmpa = cd11 * cd22 - cd12 * cd21
        cdelt1 = math.copysign(cdelt1_mag, tmpa)
        cdelt2 = cdelt2_mag
        sign = math.copysign(1, tmpa)
        crot = math.degrees(math.atan2((sign * cd12), cd22))

    elif hasattr(ipwcs.wcs, "pc"):
        for irow in range(ipwcs.naxis):
            for jcol in range(ipwcs.naxis):
                kw_str = f"PC{irow+1}_{jcol+1}"
                val_str = f"{kw_str:8s} = {ipwcs.wcs.pc[irow, jcol]}"
                dothead_list.append(val_str)
            kw_str = f"CDELT{1+irow}"
            val_str = f"{kw_str:8s} = {ipwcs.wcs.cdelt[irow]}"
            dothead_list.append(val_str)
            # Think that CD1_* = CDELT1 * PC1_* and CD2_* = CDELT2 * PC2_*
            pc = ipwcs.wcs.pc
            cdelt_vec = ipwcs.wcs.cdelt
            cd = pc.copy()
            for irow in range(ipwcs.naxis):
                cd[irow] = cdelt_vec[irow] * pc[irow]
            # then as above
            # From https://lweb.cfa.harvard.edu/~jzhao/SMA-FITS-CASA/docs/wcs88.pdf
            cd11 = cd[0, 0]
            cd12 = cd[0, 1]
            cd21 = cd[1, 0]
            cd22 = cd[1, 1]
            cdelt1_mag = math.sqrt(cd12 * cd12 + cd22 * cd22)
            cdelt2_mag = math.sqrt(cd12 * cd12 + cd22 * cd22)
            # the sign of cdelt1.cdelt2 = sign of (cd11*cd22 - cd12*cd21)
            # if the RHS is negative cdelt1 is negative by convention, so cdelt2 is always positive
            tmpa = cd11 * cd22 - cd12 * cd21
            cdelt1 = math.copysign(cdelt1_mag, tmpa)
            cdelt2 = cdelt2_mag
            sign = math.copysign(1, tmpa)
            crot = math.degrees(math.atan2((sign * cd12), cd22))

    else:
        # Assume we have a simple CDELT[12] case with CROTA
        for irow in range(ipwcs.naxis):
            kw_str = f"CDELT{1+irow}"
            val_str = f"{kw_str:8s} = {ipwcs.wcs.cdelt[irow]}"
            dothead_list.append(val_str)
            kw_str = f"CROTA{1+irow}"
            val_str = f"{kw_str:8s} = {ipwcs.wcs.crota[irow]}"
            dothead_list.append(val_str)
        cdelt1 = ipwcs.wcs.cdelt[0]
        cdelt2 = ipwcs.wcs.cdelt[1]
        crot = ipwcs.wcs.crota[1]  # CROTA2 is used, CROTA1 not used. Assume crota[0] == crota[1]

    # Convert to segagesimarl RA Dec. Assumes units are degrees, frame is ICRS
    skycrd = SkyCoord(
        ra=ipwcs.wcs.crval[0] * units.degree, dec=ipwcs.wcs.crval[1] * units.degree, frame="icrs"
    )
    skycrd_str = skycrd.to_string("hmsdms", precision=1, sep="::", pad=True)
    dothead_list.append(f"CRVAL1,2 RA, Dec as HMS, DMS     = {skycrd_str}")

    cdelt1_as = cdelt1 * 3600.0  # arcseconds
    cdelt2_as = cdelt2 * 3600.0  # arcseconds
    ximgsz_am = cdelt1 * naxis1 * 60  # arcminutes
    yimgsz_am = cdelt2 * naxis2 * 60  # arcminutes
    dothead_list.append(f"Pixel size equivalent CDELT1     = {cdelt1_as:.3f} arcseconds")
    dothead_list.append(f"Pixel size equivalent CDELT2     = {cdelt2_as:.3f} arcseconds")
    dothead_list.append(f"Image X-axis angular size        = {ximgsz_am:.3f} arcminutes")
    dothead_list.append(f"Pixel Y-axis angular size        = {yimgsz_am:.3f} arcminutes")
    dothead_list.append(f"Image rotation equivalent CROTA2 = {crot:.3f} degrees")

    dothead_list.append("END     ")
    dothead_str = "\n".join(dothead_list)
    if verbose:
        print(dothead_str)
    return (dothead_str, cdelt1_as, cdelt2_as, ximgsz_am, yimgsz_am, crot)


def make_dothead_from_keywords(
    output_swarp_dothead: str,
    naxis: list[int, int],
    crval: list[float, float],
    cdelt: list[float, float],
    crpix: list[float, float] | None = None,
    ctype: list[str, str] | None = None,
    cunit: list[str, str] | None = None,
    extnum: int | str = 0,
    format: str = "fits",
    verbose: bool = False,
) -> None:
    """
    Create the .head format file that swarp expects based on user-supplied values
    for the standard FITS NAXIS, CRVAL, CDELT, and optionally CRPIX keywords.

    When generating the ASCII format file this function uses the method shown
    in `WCS.printwcs <https://docs.astropy.org/en/stable/_modules/astropy/wcs/wcs.html#WCS.printwcs>`_.

    Note
    ~~~~

    * Input value lists are ordered followuing the FITS (Fortran) convention,
      with the X-axis value first and the Y-axis value second. For example,
      ``naxis=[3056, 2048]`` denotes an image with 2048 rows and 3056 columns.
    * Images are aligned North-up, East-left (if CDELT1 is positive). Currently
      rotations are not supported.
    * Assumes a 2-dimensional image with angular coordinate axes
    * ASCII format ``.head`` files do not appear to work with current
      versions of ``swarp``. Use the ``fits`` format instead.

    Parameters
    ----------
    output_swarp_dothead : str
        Name for the output ``.head`` file ``swarp`` will use.
    naxis : list of int
        A list containing the number of columns and rows in the final image,
        corresponding to the NAXIS1 and NAXIS2 FITS keywords.
    crval : list of int
        A list containing the Right Ascension and Declination of pixel coordinates
        corresponding to CRPIX1 and CRPIX2. Values should be in decimal degrees.
        Note that the ``crpix`` input is optional, and if not specified it is
        assumed that CRVAL1, CRVAL2 refer to the exact center of the image.
    cdelt : list of int
        A list containing the pixel size parameters CDELT1 and CDELT2.
        Values should be in decimal degrees.
    crpix : list of int, optional, default=None
        A list containing the pixel coordinates CRPIX1, CRPIX2, that correspond
        to the RA, Dec coordinates in ``crval``. This follows the FITS convention
        that the center of the first pixel is 1.0, and stretches from 0.5 to 1.5
        in pixel coordinates. If ``None`` then crpix will automatically be
        set to the center of the image, ``0.5 + 0.5*NAXIS1, 0.5 + 0.5*NAXIS2``.
    ctype : list of str, optional, default=None
        A list containing the ``CTYPE`` keywords for the X and Y axes
        respectively. If ``None`` then the default values of ``"RA---TAN", "DEC--TAN"``
        will be used.
    cunit : list of str, optional, default=None
        A list containing the ``CUNIT`` keywords for the X and Y axes
        respectively. If ``None`` then the default values of ``"deg", "deg"``
        will be used.
    extnum : int or str, optional, default=0
        Extension number or name for the extension holding the image data.
        Usually this is 0, for the ``PrimaryHDU``.
    format : {'text', 'fits'}
        If 'text' then an ASCII header will be created using the format
        specified in the Swarp manual. If 'fits' is specified, an empty FITS
        file consisting only of a primary HDU will be created.
    verbose : bool, optional, default=False
        If True then diagnostic information will be written to stdout.
    """

    # Parse the WCS keywords in the primary HDU
    ipwcs = wcs.WCS(naxis=2)
    ipwcs.array_shape = (naxis[1], naxis[0])  # Note python row, column order
    ipwcs.wcs.crval = crval
    ipwcs.wcs.cdelt = cdelt

    if ctype is None:
        ctype = ["RA---TAN", "DEC--TAN"]
    ipwcs.wcs.ctype = ctype

    if crpix is None:
        crpix = [0.5 + 0.5 * naxis[0], 0.5 + 0.5 * naxis[1]]
    ipwcs.wcs.crpix = crpix

    if cunit is None:
        cunit = ["deg", "deg"]
    ipwcs.wcs.cunit = cunit

    make_dothead_from_wcs(
        ipwcs,
        output_swarp_dothead=output_swarp_dothead,
        data=None,
        extnum=extnum,
        format=format,
        verbose=verbose,
    )
    return


def make_dothead_from_file(
    input_fits_with_wcs: str,
    output_swarp_dothead: str,
    extnum: int | str = 0,
    format: str = "fits",
    verbose: bool = False,
) -> None:
    """
    Create the .head format file that swarp expects based on the WCS header of an input files file.

    When generating the ASCII format file this function uses the method shown
    in `WCS.printwcs <https://docs.astropy.org/en/stable/_modules/astropy/wcs/wcs.html#WCS.printwcs>`_.

    Note
    ~~~~

    * Assumes a 2-dimensional image with angular coordinate axes
    * ASCII format ``.head`` files do not appear to work with current
      versions of ``swarp``. Use the ``fits`` format instead.

    Parameters
    ----------
    input_fits_with_wcs : str
        Name of existing FITS file with WCS in primary HDU that we want to emulate.
    output_swarp_dothead : str
        Name for the output ``.head`` file ``swarp`` will use.
    extnum : int or str, default=0
        Extension number or name for the extension holding the image data.
        Usually this is 0, for the ``PrimaryHDU``.
    format : {'text', 'fits'}
        If 'text' then an ASCII header will be created using the format
        specified in the Swarp manual. If 'fits' is specified, an empty FITS
        file consisting only of a primary HDU will be created.
    verbose : bool, optional, default=False
        If True then diagnostic information will be written to stdout.
    """

    with fits.open(input_fits_with_wcs) as hdulist:
        # Parse the WCS keywords in the primary HDU
        ipwcs = wcs.WCS(hdulist[0].header)
        make_dothead_from_wcs(
            ipwcs,
            output_swarp_dothead=output_swarp_dothead,
            data=None,
            extnum=extnum,
            format=format,
            verbose=verbose,
        )
    return


def make_dothead_from_wcs(
    ipwcs: Any,
    output_swarp_dothead: str,
    data: npt.NDArray = None,
    extnum: int | str = 0,
    format: str = "fits",
    verbose: bool = False,
) -> None:
    """
    Create the .head format file that swarp expects based on a WCS object.

    When generating the ASCII format file this function uses the method shown
    in `WCS.printwcs <https://docs.astropy.org/en/stable/_modules/astropy/wcs/wcs.html#WCS.printwcs>`_.

    Note
    ~~~~

    * Assumes a 2-dimensional image with angular coordinate axes
    * ASCII format ``.head`` files do not appear to work with current
      versions of ``swarp``. Use the ``fits`` format instead.

    Parameters
    ----------
    ipwcs : astropy.wcs
        Astropy WCS instance we want to emulate.
    output_swarp_dothead : str
        Name for the output ``.head`` file ``swarp`` will use.
    extnum : int or str, default=0
        Extension number or name for the extension holding the image data.
        Usually this is 0, for the ``PrimaryHDU``.
    format : {'text', 'fits'}
        If 'text' then an ASCII header will be created using the format
        specified in the Swarp manual. If 'fits' is specified, an empty FITS
        file consisting only of a primary HDU will be created.
    verbose : bool, optional, default=False
        If True then diagnostic information will be written to stdout.
    """

    if "text" in format:
        if verbose:
            print(
                "Generating an ASCII header using the format"
                " specified in the swarp documentation."
            )

        dothead_list = []
        keywords = ["NAXIS", "CTYPE", "CRVAL", "CRPIX"]
        values = [ipwcs.array_shape, ipwcs.wcs.ctype, ipwcs.wcs.crval, ipwcs.wcs.crpix]
        for keyword, value in zip(keywords, values):
            for idx in range(ipwcs.naxis):
                kw_str = f"{keyword}{1+idx}"
                if "CTYPE" in keyword:
                    # Wrap strings in quotes
                    val_str = f"{kw_str:8s} = '{value[idx]}'"
                else:
                    val_str = f"{kw_str:8s} = {value[idx]}"
                dothead_list.append(val_str)

        if hasattr(ipwcs.wcs, "pc"):
            for irow in range(ipwcs.naxis):
                for jcol in range(ipwcs.naxis):
                    kw_str = f"PC{irow+1}_{jcol+1}"
                    val_str = f"{kw_str:8s} = {ipwcs.wcs.pc[irow, jcol]}"
                    dothead_list.append(val_str)
                kw_str = f"CDELT{1+irow}"
                val_str = "f{kw_str:8s} = {pwcs.wcs.cdelt[irow]}"
                dothead_list.append(val_str)
        elif hasattr(ipwcs.wcs, "cd"):
            for irow in range(ipwcs.naxis):
                for jcol in range(ipwcs.naxis):
                    kw_str = f"CD{irow+1}_{jcol+1}"
                    val_str = f"{kw_str:8s} = {ipwcs.wcs.cd[irow, jcol]}"
                    dothead_list.append(val_str)

        dothead_list.append("END     ")
        dothead_str = "\n".join(dothead_list)

        if verbose:
            print("==== dothead output from make_dothead_from_wcs ====")
            print(dothead_str)
            print(f"==== about to write to {output_swarp_dothead} ====")

        with open(output_swarp_dothead, "w") as ofile:
            ofile.write(dothead_str)
            if verbose:
                print(f"Wrote Swarp ASCII-format .head file to {output_swarp_dothead}")
    elif "fits" in format:
        if verbose:
            print("Generating a FITS format header consisting of a PrimaryHDU only.")
        # based on https://docs.astropy.org/en/stable/wcs/example_create_imaging.html
        hdr = ipwcs.to_header()

        # NAXIS = x, NAXISx values set from data, not by manipulating header
        if data is None:
            data = np.zeros(ipwcs.array_shape, dtype=np.int16)

        hdu = fits.PrimaryHDU(header=hdr, data=data)
        hdu.writeto(output_swarp_dothead, overwrite=True, output_verify="ignore")
        if verbose:
            print(f"Wrote FITS header .head file to {output_swarp_dothead}")
    else:
        raise RuntimeError(
            f'Error, format ({format}) is not one of the allowed options: "text" "fits"'
        )
    return


def get_exposure_time(hdr, verbose=False):
    """
    Return the first of EXPTIME, EXPOSURE, ONTIME, or LIVETIME from a FITS header

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        Input FITS Header instance
    verbose : bool, optional, default=False
        If True then diagnostic information will be written to stdout.

    Returns
    -------
    exposure_time : float or None
        The exposure time in seconds, if found in the FITS header object.
        Otherwise None.

    See Also
    --------

    `HEASARC Commonly Used FITS Keywords <https://heasarc.gsfc.nasa.gov/docs/fcg/common_dict.html>`_
    """
    keywords = ["EXPTIME", "EXPOSURE", "ONTIME", "LIVETIME", "TELAPSE", "ELAPTIME"]
    exposure_time = None
    for key in keywords:
        if hdr.count(key) > 0:
            exposure_time = float(hdr[key])
            if verbose:
                print(f"Keyword {key} found, setting exposure_time to {exposure_time}")
            break
        else:
            if verbose:
                print(f"Keyword {key} not found, continuing search...")
    return exposure_time
