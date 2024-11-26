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

##import logging
##import os.path
from typing import Any
from pathlib import Path
import numpy as np

##import matplotlib  # for rc
import matplotlib.pyplot as plt

##from matplotlib.patches import Ellipse
import math
##import time
##from datetime import datetime

from astropy import wcs

##import astropy
from astropy.io import fits

##from astropy.table import QTable, Table, vstack
from astropy.coordinates import SkyCoord, Angle
from astropy import units
from astropy.visualization import (
    AsymmetricPercentileInterval,
    MinMaxInterval,
    ManualInterval,
    SqrtStretch,
    AsinhStretch,
    LinearStretch,
    ImageNormalize,
)
from astropy.visualization import make_lupton_rgb

# AstroPhotography includes
##from .. import __version__


def does_file_exist(filename: str, verbose: bool = False) -> bool:
    """
    Returns True if the file name or path exists, false otherwise

    Parameters
    ----------
    filename : str
        Name (optionally including path) of the file we want to check
        the existence of.
    verbose : bool, optional, default=False
        If True then writes to stdout.

    Returns
    -------
    exists : bool
        Returns True if the file path exists, False otherwise.
    """
    if verbose and not Path(filename).exists():
        print(f"Cannot find {filename}. Not a valid path or file.")
        return False
    else:
        print(f"Found {filename}.")
    return True


def load_image_and_plot(
    fname: str,
    extnum: int | str = 0,
    output: str | None = None,
    usewcs: bool = True,
    vmin: float = None,
    vmax: float = None,
    xaxlim: Any = None,
    yaxlim: Any = None,
    verbose: bool = True,
    angle_tick_spacing_am: float = 2.0,
    swap_radec_axis: bool = False,
):
    """
    Quick and dirty FITS image plot.

    By default `asinh` scaling will be used between a minimum and maximum data value. If vmin
    and/or vmnax are not specified then the 0.5th and 99.5th perciles will be used, as appropriate.
    A viridis color map is used by default, arag with a color bar.

    Note
    ----
    Axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis
    pixels.

    Astropy WCSAxes transforms do not reorient the image to be North up, East left, in the same
    way SAOImage ds9 does when the WCS is applied. It still plots the original data rows and
    column as a 2-D XY grid. Consequently, for images not already in NE alignment, the RA may
    change most rapidly along Y and Dec most rapidly along X, in contrast to normal expectation.
    To avoid highly confusing plots this function plots RA/Dec grids when `usewcs=True`. Line
    of constant RA are red, lines of constant Declination are blue. In cases where your data
    is aligned closer to 90 degrees or 270 degrees away from North up, East left, specifying
    `swap_radec_axis=True` will reduce confusion by showing RA tickmarks and values along the Y
    axes and Declination tickmarks and values along the X axis (contrary to the normal convention).

    If your `usewcs=True` plot appears without axis tick values and labels then it is likely you
    need to alter `swap_radec_axis` or `angle_tick_spacing_am`

    TODO: Find out how to specify them in RA and Dec.

    Parameters
    ----------

    fname : str
        Name of existing FITS file with WCS in HDU number extnum that we want to plot.
    extnum : int or str, optional, default=0
        Extension number (zero-based) or extension name (string) for data and WCS header
    output : str, 'interactive', or None, optional, default=None
        Output format for the plot. If 'interactive' is specified then the
        pyplot interactive plotting device pyplot.show() is invoked. This is a
        blocking operation. The 'interactive' mode may not work when running
        in jupyter.
        Use None when using a jupyter backend (e.g. ``%matplotlib inline``)
        as they automatically invoke a non-blocking show at the end of
        each cell.
        Specify a filename with recognizable image-format extension, e.g.
        ``my_output_plot.png`` to directly save the plot to a file.
    usewcs : bool, optional, default=True
        If True then plot RA/Dec axes using WCS information
    vmin : float or None, optional, default=None
        If not None then vmin is the minimum data value in the normalization interval.
    vmax : float or None, optional, default=None
        If not None then vmax is the maximum data value in the normalization interval.
    xaxlim : array_like or None, optional, default=None
        2-element list or tuple of the minimum to maximum x-axis coordinates to
        plot. If None then the default axis limits will be used. To see those limits run
        with verbose=True.
    yaxlim : array_like or None, optional, default=None
        2-element list or tuple of the minimum to maximum x-axis coordinates to
        plot. If None then the default axis limits will be used. To see those limits run
        with verbose=True.
    verbose : bool, optional, default=True
        If True then diagnostic information will be written to stdout.
    angle_tick_spacing_am : float, optional, default=2.0
        If usewcs is True, specify the spacing between RA/Dec axis
        tick values in units of arcminutes. Astropy won't label the RA/Dec axis without
        specifying a value for this.
    swap_radec_axis : bool, optional, default=False
        If True then plot RA tickmarks and values along the Y
        axes and Declination tickmarks and values along the X axis (contrary to the normal
        convention). This does not alter how the image data itself is plotted.
    """

    hdu = fits.open(fname)[extnum]
    w = None
    if usewcs:
        w = wcs.WCS(hdu.header)

    default_font_size = 7

    # Compute vmin and vmax if necessary
    # percentiles
    ipctls = [0.5, 99.5]
    opctls = np.nanpercentile(hdu.data, ipctls)
    ourmin = opctls[0]
    ourmax = opctls[1]
    if verbose:
        print(f"Input data 0.5th percentile = {ourmin:.4f}, 99.5th percentile = {ourmax:.4f}")

    if vmin is not None:
        ourmin = vmin
    if vmax is not None:
        ourmax = vmax
    if verbose:
        print(f"Applying asinh stretch between {ourmin:.4f} and {ourmax:.4f}")

    # Create an ImageNormalize object
    norm = ImageNormalize(
        hdu.data, interval=ManualInterval(ourmin, ourmax), stretch=AsinhStretch()
    )

    # Display the image
    fig = plt.figure()
    if usewcs:
        if verbose:
            print("Using WCS information for axis projection")
            print(w)
        ax = fig.add_subplot(1, 1, 1, projection=w)
    else:
        ax = fig.add_subplot(1, 1, 1)

    im = ax.imshow(hdu.data, origin="lower", norm=norm)
    cbar = fig.colorbar(
        im, ax=ax, extend="neither", spacing="proportional", orientation="vertical", shrink=0.85
    )
    cbar.set_label(r"Units TBA")
    cbar.ax.tick_params(labelsize=default_font_size)
    title_str = fname  # .replace("_", "\_")
    ax.set_title(f"{title_str}", fontsize=default_font_size)

    # Display default axis limits
    xlim_used = ax.get_xbound()
    ylim_used = ax.get_ybound()
    if verbose:
        print(f"Default X-axis limits: {xlim_used}")
        print(f"Default Y-axis limits: {ylim_used}")

    title_str = fname  # .replace("_", "\_")
    ax.set_title(f"{title_str}", fontsize=default_font_size)
    if usewcs:
        # Adapted rom @astrofrog at https://github.com/astropy/astropy/issues/13458#issuecomment-1242640539
        ra = ax.coords[0]
        dec = ax.coords[1]
        ra.set_major_formatter("hh:mm:ss.ss")  # RA in Hours, minutes, seconds,
        dec.set_major_formatter("dd:mm:ss.ss")
        ra.set_ticks(
            spacing=angle_tick_spacing_am * units.arcmin, color="red"
        )  # Ticks must be defined for axis tickvals to appear
        dec.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color="blue")
        ra.set_ticklabel(color="red", fontsize=default_font_size)
        dec.set_ticklabel(color="blue", fontsize=default_font_size)
        ra.grid(color="red", linestyle="--", alpha=0.6)
        dec.grid(color="blue", linestyle="--", alpha=0.6)
        ra.set_axislabel("Right Ascension (HMS)", fontsize=default_font_size, color="red")
        dec.set_axislabel("Declination (dms)", fontsize=default_font_size, color="blue")

        if swap_radec_axis:
            # Images where ra changes fastest on Y, not X
            dec.set_ticks_position("b")
            dec.set_ticklabel_position("b")
            dec.set_axislabel_position("b")
            ra.set_ticks_position("l")
            ra.set_ticklabel_position("l")
            ra.set_axislabel_position("l")

        # ax.set_xlabel('Right Ascension (deg)', fontsize=default_font_size)
        # ax.set_ylabel('Declination (deg)', fontsize=default_font_size)
    else:
        ax.set_xlabel("X-axis pixel number", fontsize=default_font_size)
        ax.set_ylabel("Y-axis pixel number", fontsize=default_font_size)

    # Modify the axis limits?
    if xaxlim is not None:
        loval = np.min(xaxlim)
        hival = np.max(xaxlim)
        ax.set_xbound(lower=loval, upper=hival)
        if verbose:
            print(f"Setting X-axis limits to [{loval}, {hival}]")
    if yaxlim is not None:
        loval = np.min(yaxlim)
        hival = np.max(yaxlim)
        ax.set_ybound(lower=loval, upper=hival)
        if verbose:
            print(f"Setting Y-axis limits to [{loval}, {hival}]")

    if output is not None:
        if "interactive" in output:
            print("Processing blocked while plotting window is open. Close it to procede.")
            fig.show()
        else:
            fig.savefig(output, dpi=200, bbox_inches="tight")
            if verbose:
                print(f"Wrote plot to {output}")
    return


def load_imlist_and_plot(
    imlist: list[str],
    extnum: int = 0,
    output: str | None = None,
    usewcs: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    xaxlim=None,
    yaxlim=None,
    verbose: bool = True,
    angle_tick_spacing_am: float = 2.0,
    swap_radec_axis: bool = False,
):
    """
    Quick and dirty multi-FITS image plotter.

    By default `asinh` scaling will be used between a minimum and maximum data value. If vmin
    and/or vmnax are not specified then the 0.5th and 99.5th perciles will be used, as appropriate.
    A viridis color map is used by default, along with a color bar.

    Note that axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis
    pixels. ( TODO: Find out how to specify axis limits in RA and Dec.)

    Note
    ----

    Axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis
    pixels.

    Astropy WCSAxes transforms do not reorient the image to be North up, East left, in the same
    way SAOImage ds9 does when the WCS is applied. It still plots the original data rows and
    column as a 2-D XY grid. Consequently, for images not already in NE alignment, the RA may
    change most rapidly along Y and Dec most rapidly along X, in contrast to normal expectation.
    To avoid highly confusing plots this function plots RA/Dec grids when `usewcs=True`. Line
    of constant RA are red, lines of constant Declination are blue. In cases where your data
    is aligned closer to 90 degrees or 270 degrees away from North up, East left, specifying
    `swap_radec_axis=True` will reduce confusion by showing RA tickmarks and values along the Y
    axes and Declination tickmarks and values along the X axis (contrary to the normal convention).

    If your `usewcs=True` plot appears without axis tick values and labels then it is likely you
    need to alter `swap_radec_axis` or `angle_tick_spacing_am`

    Parameters
    ----------
    imlist : list of str
        A list of names of existing FITS file we want plotted. All files must
        have the data in the same HDU number (or name) extnum. Files will
        be plotted in the order that they appear in the list.
    extnum : int or str, optional, default=0
        Extension number (zero-based) or extension name (string) for data and WCS header
    output : str, 'interactive', or None, optional, default=None
        Output format for the plot. If 'interactive' is specified then the
        pyplot interactive plotting device pyplot.show() is invoked. This is a
        blocking operation. The 'interactive' mode may not work when running
        in jupyter.
        Use None when using a jupyter backend (e.g. ``%matplotlib inline``)
        as they automatically invoke a non-blocking show at the end of
        each cell.
        Specify a filename with recognizable image-format extension, e.g.
        ``my_output_plot.png`` to directly save the plot to a file.
    usewcs : bool, optional, default=True
        If True then plot RA/Dec axes using WCS information
    vmin : float or None, optional, default=None
        If not None then vmin is the minimum data value in the normalization interval.
    vmax : float or None, optional, default=None
        If not None then vmax is the maximum data value in the normalization interval.
    xaxlim : array_like or None, optional, default=None
        2-element list or tuple of the minimum to maximum x-axis coordinates to
        plot. If None then the default axis limits will be used. To see those limits run
        with verbose=True.
    yaxlim : array_like or None, optional, default=None
        2-element list or tuple of the minimum to maximum x-axis coordinates to
        plot. If None then the default axis limits will be used. To see those limits run
        with verbose=True.
    verbose : bool, optional, default=True
        If True then diagnostic information will be written to stdout.
    angle_tick_spacing_am : float, optional, default=2.0
        If usewcs is True, specify the spacing between RA/Dec axis
        tick values in units of arcminutes. Astropy won't label the RA/Dec axis without
        specifying a value for this.
    swap_radec_axis : bool, optional, default=False
        If True then plot RA tickmarks and values along the Y
        axes and Declination tickmarks and values along the X axis (contrary to the normal
        convention). This does not alter how the image data itself is plotted.
    """

    # Determine subplot layout
    num_imgs: int = len(imlist)
    print(f"The are {num_imgs} FITS images in the input imlist")
    if num_imgs > 25:
        panel_rows: int = int(math.ceil(math.sqrt(num_imgs)))
        panel_cols: int = panel_rows
    else:
        # num imgs <= 1, 2, 4, 6, 8, 9,  12, 16, 20, 25
        nimlt_arr = [2, 3, 5, 7, 9, 10, 13, 17, 21, 26]
        nrows_arr = [1, 2, 2, 3, 4, 3, 4, 4, 5, 5]
        ncols_arr = [1, 1, 2, 2, 2, 3, 3, 4, 4, 5]
        for idx, thresh in enumerate(nimlt_arr):
            if num_imgs < thresh:
                panel_rows: int = nrows_arr[idx]
                panel_cols: int = ncols_arr[idx]
                break
    print(
        (
            f"  With {num_imgs} files the plot layout is "
            f"{panel_rows} rows x {panel_cols} columns."
        )
    )

    default_font_size = 7

    # Percentiles to use
    ipctls = [0.5, 99.5]

    # Display the image
    fig = plt.figure()

    for idx in range(num_imgs):
        fname = imlist[idx]
        print(f"Attempting to plot {fname}")
        with fits.open(fname) as hdulist:
            imhdu = hdulist[extnum]

            # Compute vmin and vmax if necessary.
            opctls = np.nanpercentile(imhdu.data, ipctls)
            ourmin = opctls[0]
            ourmax = opctls[1]
            if verbose:
                print(
                    (
                        f"  Input image data 0.5th percentile = {ourmin:.4f}"
                        f", 99.5th percentile = {ourmax:.4f}"
                    )
                )
            if vmin is not None:
                ourmin = vmin
            if vmax is not None:
                ourmax = vmax
            if verbose:
                print(f"  Applying asinh stretch between {ourmin:.4f} and {ourmax:.4f}")

            title_str = fname
            w = wcs.WCS(imhdu.header)
            if usewcs:
                ax = fig.add_subplot(panel_rows, panel_cols, idx + 1, projection=w)
            else:
                ax = fig.add_subplot(panel_rows, panel_cols, idx + 1)

            # Create an ImageNormalize object
            norm = ImageNormalize(
                imhdu.data, interval=ManualInterval(ourmin, ourmax), stretch=AsinhStretch()
            )

            im = ax.imshow(imhdu.data, origin="lower", norm=norm)
            cbar = fig.colorbar(
                im,
                ax=ax,
                extend="neither",
                spacing="proportional",
                orientation="vertical",
                shrink=0.85,
            )
            cbar.set_label(r"Units TBA")
            cbar.ax.tick_params(labelsize=default_font_size)

            # Display default axis limits
            xlim_used = ax.get_xbound()
            ylim_used = ax.get_ybound()
            if verbose:
                print(f"Default X-axis limits: {xlim_used}")
                print(f"Default Y-axis limits: {ylim_used}")

            ax.set_title(f"{title_str}", fontsize=default_font_size)
            if usewcs:
                # Adapted from @astrofrog at https://github.com/astropy/astropy/issues/13458#issuecomment-1242640539
                ra = ax.coords[0]
                dec = ax.coords[1]
                ra.set_major_formatter("hh:mm:ss.ss")  # RA in Hours, minutes, seconds,
                dec.set_major_formatter("dd:mm:ss.ss")
                ra.set_ticks(
                    spacing=angle_tick_spacing_am * units.arcmin, color="red"
                )  # Ticks must be defined for axis tickvals to appear
                dec.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color="blue")
                ra.set_ticklabel(color="red", fontsize=default_font_size)
                dec.set_ticklabel(color="blue", fontsize=default_font_size)
                ra.grid(color="red", linestyle="--", alpha=0.6)
                dec.grid(color="blue", linestyle="--", alpha=0.6)
                ra.set_axislabel("Right Ascension (HMS)", fontsize=default_font_size, color="red")
                dec.set_axislabel("Declination (dms)", fontsize=default_font_size, color="blue")

                if swap_radec_axis:
                    # Images where ra changes fastest on Y, not X
                    dec.set_ticks_position("b")
                    dec.set_ticklabel_position("b")
                    dec.set_axislabel_position("b")
                    ra.set_ticks_position("l")
                    ra.set_ticklabel_position("l")
                    ra.set_axislabel_position("l")

                # ax.set_xlabel('Right Ascension (deg)', fontsize=default_font_size)
                # ax.set_ylabel('Declination (deg)', fontsize=default_font_size)
            else:
                ax.set_xlabel("X-axis pixel number", fontsize=default_font_size)
                ax.set_ylabel("Y-axis pixel number", fontsize=default_font_size)

            # Modify the axis limits?
            if xaxlim is not None:
                loval = np.min(xaxlim)
                hival = np.max(xaxlim)
                ax.set_xbound(lower=loval, upper=hival)
                if verbose:
                    print(f"Setting X-axis limits to [{loval}, {hival}]")
            if yaxlim is not None:
                loval = np.min(yaxlim)
                hival = np.max(yaxlim)
                ax.set_ybound(lower=loval, upper=hival)
                if verbose:
                    print(f"Setting Y-axis limits to [{loval}, {hival}]")

    if output is not None:
        if "interactive" in output:
            print("Processing blocked while plotting window is open. Close it to procede.")
            fig.show()
        else:
            fig.savefig(output, dpi=200, bbox_inches="tight")
            if verbose:
                print(f"Wrote plot to {output}")
    return


def plot_lupton_threecolor(
    redfile: str,
    greenfile: str,
    bluefile: str,
    outpngfile: str | None,
    ax: Any = None,
    extnum: int | str = 0,
    usewcs: bool = True,
    vmin: float | None = None,
    xaxlim=None,
    yaxlim=None,
    qval: float = 10,
    stretchval: float = 0.5,
    verbose: bool = True,
    angle_tick_spacing_am: float = 2.0,
    swap_radec_axis: bool = False,
):
    """
    Quick and dirty 3-color composite using make_lupton_rgb

    - All images must share the same size, WCS orientation, and pixel size.
    - Axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis pixels.

    Note
    ----
    Axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis
    pixels. (TODO: Find out how to specify them in RA and Dec.)

    Astropy WCSAxes transforms do not reorient the image to be North up, East left, in the same
    way SAOImage ds9 does when the WCS is applied. It still plots the original data rows and
    column as a 2-D XY grid. Consequently, for images not already in NE alignment, the RA may
    change most rapidly along Y and Dec most rapidly along X, in contrast to normal expectation.
    To avoid highly confusing plots this function plots RA/Dec grids when `usewcs=True`. Line
    of constant RA are red, lines of constant Declination are blue. In cases where your data
    is aligned closer to 90 degrees or 270 degrees away from North up, East left, specifying
    `swap_radec_axis=True` will reduce confusion by showing RA tickmarks and values along the Y
    axes and Declination tickmarks and values along the X axis (contrary to the normal convention).

    If your `usewcs=True` plot appears without axis tick values and labels then it is likely you
    need to alter `swap_radec_axis` or `angle_tick_spacing_am`

    Parameters
    ----------
    redfile : str
        Name of existing FITS file with WCS in HDU number extnum that we want to
        appear as red in the RGB color composite.
    greenfile : str
        Name of existing FITS file with WCS in HDU number extnum that we want to
        appear as green in the RGB color composite.
    bluefile : str
        Name of existing FITS file with WCS in HDU number extnum that we want to
        appear as blue in the RGB color composite.
    outpngfile : str, 'interactive', or None, optional, default=None
        Output format for the plot. If 'interactive' is specified then the
        pyplot interactive plotting device pyplot.show() is invoked. This is a
        blocking operation. The 'interactive' mode may not work when running
        in jupyter.
        Use None when using a jupyter backend (e.g. ``%matplotlib inline``)
        as they automatically invoke a non-blocking show at the end of
        each cell.
        Specify a filename with recognizable image-format extension, e.g.
        ``my_output_plot.png`` to directly save the plot to a file.
    ax : matplotlib.Axes instance or None, optional, default=None
        This function supports being called from another function controlling
        a figure and set of Axes instances. If ax is not None then this function
        will not generated its own figure but instead use the specified Axes.
        Note that in such a case outpngfile should be None, **and** the
        caller must have created the Axes with the correct WCS project if
        WCS axes plotting is wanted.
    extnum : int or str, optional, default=0
        Extension number (zero-based) or extension name (string) for data and WCS header
    usewcs : bool, optional, default=True
        If True then plot RA/Dec axes using WCS information
    vmin : float or None, optional, default=None
        If not None then vmin is the minimum data value in the normalization interval.
    vmax : float or None, optional, default=None
        If not None then vmax is the maximum data value in the normalization interval.
    xaxlim : array_like or None, optional, default=None
        2-element list or tuple of the minimum to maximum x-axis coordinates to
        plot. If None then the default axis limits will be used. To see those limits run
        with verbose=True.
    yaxlim : array_like or None, optional, default=None
        2-element list or tuple of the minimum to maximum x-axis coordinates to
        plot. If None then the default axis limits will be used. To see those limits run
        with verbose=True.
    qval : float, optional, default=10
        Value of make_lupton_rgb Q parameter (Q in Lupton et al 2004)
    stretchval : float, optional, default=0.5
        Value of make_lupton_rgb stretch parameter (alpha in Lupton et al 2004)
    verbose : bool, optional, default=True
        If True then diagnostic information will be written to stdout.
    angle_tick_spacing_am : float, optional, default=2.0
        If usewcs is True, specify the spacing between RA/Dec axis
        tick values in units of arcminutes. Astropy won't label the RA/Dec axis without
        specifying a value for this.
    swap_radec_axis : bool, optional, default=False
        If True then plot RA tickmarks and values along the Y
        axes and Declination tickmarks and values along the X axis (contrary to the normal
        convention). This does not alter how the image data itself is plotted.

    Returns
    -------
    rgb_array : ndarray
        The 3-dimensional data array produced by astropy's make_lupton_rgb function.

    """

    hdulist_r = fits.open(redfile)
    hdulist_g = fits.open(greenfile)
    hdulist_b = fits.open(bluefile)

    hdur = hdulist_r[extnum]
    hdug = hdulist_g[extnum]
    hdub = hdulist_b[extnum]
    w = wcs.WCS(hdur.header)

    default_font_size = 7

    # Compute vmin and vmax if necessary. NOTE this only uses the red channel
    # percentiles
    ipctls = [0.5, 99.5]
    opctls = np.nanpercentile(hdur.data, ipctls)
    ourmin = opctls[0]
    ourmax = opctls[1]
    if verbose:
        print(
            (
                f"Input red channel data 0.5th percentile = {ourmin:.4f}"
                f", 99.5th percentile = {ourmax:.4f}"
            )
        )

    vmax = None
    if vmin is not None:
        ourmin = vmin
    if vmax is not None:
        ourmax = vmax
    if verbose:
        print(f"Using an image minimum of {ourmin} for all channels.")

    rgb_array = make_lupton_rgb(
        hdur.data,
        hdug.data,
        hdub.data,
        minimum=ourmin,
        stretch=stretchval,
        Q=qval,
        filename=outpngfile,
    )

    # Display the image
    panel_plotter = False
    if ax is None:
        # Use our own figure and axes
        fig = plt.figure()
        if usewcs:
            ax = fig.add_subplot(1, 1, 1, projection=w)
        else:
            ax = fig.add_subplot(1, 1, 1)
    else:
        # Plot into externally provided Axes, don't create a PNG ourselves
        panel_plotter = True
        fig = plt.gcf()

    im = ax.imshow(rgb_array, origin="lower")
    cbar = fig.colorbar(
        im, ax=ax, extend="neither", spacing="proportional", orientation="vertical", shrink=0.85
    )
    cbar.set_label(r"Units TBA", fontsize=default_font_size)
    cbar.ax.tick_params(labelsize=default_font_size)

    # Display default axis limits
    xlim_used = ax.get_xbound()
    ylim_used = ax.get_ybound()
    if verbose:
        print(f"Default X-axis limits: {xlim_used}")
        print(f"Default Y-axis limits: {ylim_used}")

    title_str = (
        f"RGB composite image stretch={stretchval}, Q={qval}, minimum={vmin}"
        f"\nRed = {redfile}"
        f"\nGreen = {greenfile}"
        f"\nBlue = {bluefile}\n"
    )
    ax.set_title(f"{title_str}", fontsize=default_font_size)
    if usewcs:
        # Adapted rom @astrofrog at https://github.com/astropy/astropy/issues/13458#issuecomment-1242640539
        ra = ax.coords[0]
        dec = ax.coords[1]
        ra.set_major_formatter("hh:mm:ss.ss")  # RA in Hours, minutes, seconds,
        dec.set_major_formatter("dd:mm:ss.ss")
        ra.set_ticks(
            spacing=angle_tick_spacing_am * units.arcmin, color="red"
        )  # Ticks must be defined for axis tickvals to appear
        dec.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color="blue")
        ra.set_ticklabel(color="red", fontsize=default_font_size)
        dec.set_ticklabel(color="blue", fontsize=default_font_size)
        ra.grid(color="red", linestyle="--", alpha=0.6)
        dec.grid(color="blue", linestyle="--", alpha=0.6)
        ra.set_axislabel("Right Ascension (HMS)", fontsize=default_font_size, color="red")
        dec.set_axislabel("Declination (dms)", fontsize=default_font_size, color="blue")

        if swap_radec_axis:
            # Images where ra changes fastest on Y, not X
            dec.set_ticks_position("b")
            dec.set_ticklabel_position("b")
            dec.set_axislabel_position("b")
            ra.set_ticks_position("l")
            ra.set_ticklabel_position("l")
            ra.set_axislabel_position("l")

        # ax.set_xlabel('Right Ascension (deg)', fontsize=default_font_size)
        # ax.set_ylabel('Declination (deg)', fontsize=default_font_size)
    else:
        ax.set_xlabel("X-axis pixel number", fontsize=default_font_size)
        ax.set_ylabel("Y-axis pixel number", fontsize=default_font_size)

    # Modify the axis limits?
    if xaxlim is not None:
        loval = np.min(xaxlim)
        hival = np.max(xaxlim)
        ax.set_xbound(lower=loval, upper=hival)
        if verbose:
            print(f"Setting X-axis limits to [{loval}, {hival}]")
    if yaxlim is not None:
        loval = np.min(yaxlim)
        hival = np.max(yaxlim)
        ax.set_ybound(lower=loval, upper=hival)
        if verbose:
            print(f"Setting Y-axis limits to [{loval}, {hival}]")

    output = outpngfile
    if not panel_plotter and output is not None:
        if "interactive" in output:
            print("Processing blocked while plotting window is open. Close it to procede.")
            fig.show()
        else:
            fig.savefig(output, dpi=200, bbox_inches="tight")
            if verbose:
                print(f"Wrote plot to {output}")

    hdulist_r.close()
    hdulist_g.close()
    hdulist_b.close()
    return rgb_array


def make_lupton_threecolor_plots(
    filter_imdict: dict[str, str],
    outpngfile: str,
    extnum: int | str = 0,
    usewcs: bool = True,
    vmin: float | None = None,
    xaxlim=None,
    yaxlim=None,
    qval: float = 10,
    stretchval: float = 0.5,
    verbose: bool = True,
    angle_tick_spacing_am: float = 2.0,
    swap_radec_axis: bool = False,
):
    """
    Create a single summary plot containing any RGB and/or SHO composite plots
    that can be generated from a set of images in multiple filters.

    Given a dictionary of ``dict{filter_name, image_name}`` this function will
    create a single panel color composite using :func:`plot_lupton_threecolor`
    for each of the following **complete** filter sets:

    - ``Red``, ``Green``, and ``Blue``: plotted in that order.
    - ``SII``, ``Ha``, and ``OIII``: plotted in the standard Hubble SHO order
      with SII as red, H-alpha as green, and OIII as blue.

    If images in the required filters for a given three color composite are
    not present then no panel will be generated.

    Parameters
    ----------
    filter_imdict : dict of str, str
        A dictionary consisting of filter name, image name pairs.
    outpngfile : str, 'interactive', or None, optional, default=None
        Output format for the plot. If 'interactive' is specified then the
        pyplot interactive plotting device pyplot.show() is invoked. This is a
        blocking operation. The 'interactive' mode may not work when running
        in jupyter.
        Use None when using a jupyter backend (e.g. ``%matplotlib inline``)
        as they automatically invoke a non-blocking show at the end of
        each cell.
        Specify a filename with recognizable image-format extension, e.g.
        ``my_output_plot.png`` to directly save the plot to a file.
    extnum : int or str, optional, default=0
        Extension number (zero-based) or extension name (string) for data and WCS header
    usewcs : bool, optional, default=True
        If True then plot RA/Dec axes using WCS information
    vmin : float or None, optional, default=None
        If not None then vmin is the minimum data value in the normalization interval.
    vmax : float or None, optional, default=None
        If not None then vmax is the maximum data value in the normalization interval.
    xaxlim : array_like or None, optional, default=None
        2-element list or tuple of the minimum to maximum x-axis coordinates to
        plot. If None then the default axis limits will be used. To see those limits run
        with verbose=True.
    yaxlim : array_like or None, optional, default=None
        2-element list or tuple of the minimum to maximum x-axis coordinates to
        plot. If None then the default axis limits will be used. To see those limits run
        with verbose=True.
    qval : float, optional, default=10
        Value of make_lupton_rgb Q parameter (Q in Lupton et al 2004)
    stretchval : float, optional, default=0.5
        Value of make_lupton_rgb stretch parameter (alpha in Lupton et al 2004)
    verbose : bool, optional, default=True
        If True then diagnostic information will be written to stdout.
    angle_tick_spacing_am : float, optional, default=2.0
        If usewcs is True, specify the spacing between RA/Dec axis
        tick values in units of arcminutes. Astropy won't label the RA/Dec axis without
        specifying a value for this.
    swap_radec_axis : bool, optional, default=False
        If True then plot RA tickmarks and values along the Y
        axes and Declination tickmarks and values along the X axis (contrary to the normal
        convention). This does not alter how the image data itself is plotted.
    """

    panel_cols: int = 1
    panel_rows: int = 0

    # Look for RGB or SHO
    rgb_satisfied: bool = False
    rgb_list: list[str] = []
    rgb_required: list[str] = ["Red", "Green", "Blue"]
    for filter_name in rgb_required:
        if filter_name in filter_imdict.keys():
            rgb_list.append(filter_imdict[filter_name])
    if len(rgb_list) == 3:
        rgb_satisfied = True
        panel_rows += 1
        print("There are enough files to generate a RGB threecolor composite.")
    else:
        print("NB: insufficient files to generate a RGB threecolor composite.")

    sho_satisfied: bool = False
    sho_list: list[str] = []
    sho_required: list[str] = ["SII", "Ha", "OIII"]
    for filter_name in sho_required:
        if filter_name in filter_imdict.keys():
            sho_list.append(filter_imdict[filter_name])
    if len(sho_list) == 3:
        sho_satisfied = True
        panel_rows += 1
        print("There are enough files to generate a SHO threecolor composite.")
    else:
        print("NB: insufficient files to generate a SHO threecolor composite.")

    if not (rgb_satisfied or sho_satisfied):
        print("Returning as cannot generated either RGB or SHO composities.")
        return

    # Use our own figure and axes
    fig = plt.figure()

    iplot: int = 0
    if rgb_satisfied:
        iplot += 1
        hdulist = fits.open(rgb_list[0])
        hdu1 = hdulist[extnum]

        if usewcs:
            w = wcs.WCS(hdu1.header)
            ax = fig.add_subplot(panel_rows, panel_cols, iplot, projection=w)
        else:
            ax = fig.add_subplot(panel_rows, panel_cols, iplot)

        _ = plot_lupton_threecolor(
            rgb_list[0],
            rgb_list[1],
            rgb_list[2],
            outpngfile=None,
            ax=ax,
            extnum=extnum,
            usewcs=usewcs,
            vmin=vmin,
            xaxlim=xaxlim,
            yaxlim=yaxlim,
            qval=qval,
            stretchval=stretchval,
            verbose=verbose,
            angle_tick_spacing_am=angle_tick_spacing_am,
            swap_radec_axis=swap_radec_axis,
        )
        hdulist.close()

    if sho_satisfied:
        iplot += 1
        hdulist = fits.open(rgb_list[0])
        hdu1 = hdulist[extnum]

        if usewcs:
            w = wcs.WCS(hdu1.header)
            ax = fig.add_subplot(panel_rows, panel_cols, iplot, projection=w)
        else:
            ax = fig.add_subplot(panel_rows, panel_cols, iplot)

        _ = plot_lupton_threecolor(
            sho_list[0],
            sho_list[1],
            sho_list[2],
            outpngfile=None,
            ax=ax,
            extnum=extnum,
            usewcs=usewcs,
            vmin=vmin,
            xaxlim=xaxlim,
            yaxlim=yaxlim,
            qval=qval,
            stretchval=stretchval,
            verbose=verbose,
            angle_tick_spacing_am=angle_tick_spacing_am,
            swap_radec_axis=swap_radec_axis,
        )
        hdulist.close()

    output = outpngfile
    if output is not None:
        if "interactive" in output:
            print("Processing blocked while plotting window is open. Close it to procede.")
            fig.show()
        else:
            fig.savefig(output, dpi=200, bbox_inches="tight")
            if verbose:
                print(f"Wrote plot to {output}")
    return


def namefn_calibrated_input(input_file, input_rootname, input_suffix, ap_filetype):
    """
    Given the name of a calibrated input file, generate the file name for
    and output directory for one of the subsequent output file types.

    For calibrated input files the the allowed output file types are:

    - ``input``: This is the calibrated input file itself.
    - ``srclist``: Star detection output FITS table file.
    - ``regfile``: ds9-format region file.
    - ``plotfile``: PNG plots of the input image with the detected star-like sources
      plotted as circles.
    - ``qualfile``: YaML source detection quality summary file.
    - ``fwhmplot``: PNG plot zooming in around a subset of the detected
      sources in the input image.
    - ``navfile``: A navigated and renamed copy of the input file, but now
      with a valid WCS solution.

    Parameters
    ----------
    input_file : str
        The name of the calibrated input file.
    input_rootname : str, optional, default=None
               The part of the input file name that are shared and that
               designate them being the calibrated files. If not specified
               it is assumed that the file is from iTelescope or was
               created by the AstroPhotography module itself, and will
               have an input_rootname of either 'Calibrated-iTelescope',
               'calibrated', or just 'cal' .
               The output file name from this function replace the rootname with
               a file-type specific prefix, as described in the class
               documentation and get_file_names function documentation.
    input_suffix : str, optional, default='.fits'
               String denoting the file type suffix of the input image file.
               For example, '.fits' or '.fits' or '.ftz' or 'fits.gz' or '.fits.bz2'
    ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile'}
        The output file type for which the name should be returned.
        This should be one of the stage names described above.

    Returns
    -------
    output_file : str
        File name that would correspond to the named output file type given
        the name of the input file. This file does not necessarily exist yet.
        The output file name does **not** include any subdirectory name.
        It is up to the caller to combine output_dir and output_file if
        necessary.
    output_dir : str
        Directory name in which output_file will be written, relative
        to the root directory established by the input files.

    Raises
    ------
    RuntimeError :
        If input_rootname is None but the input_file name does not match
        any of the expected iTelescope/AstroPhotography root names.

    Warnings
    --------
    This function assumes that the input file is not in a subdirectory,
    i.e. there is no path within in ``input_file`` string. In the longer
    term this should be rewritten using the pathlib module without such
    an assumption.
    """

    allowed_stages = ["input", "srclist", "regfile", "plotfile", "qualfile", "fwhmplot", "navfile"]
    if ap_filetype not in allowed_stages:
        err_msg = (
            f"Requested ap_filetype {ap_filetype} not one of the allowed values: {allowed_stages}"
        )
        raise RuntimeError(err_msg)

    # Determine the root name to use.
    default_roots = ["Calibrated-iTelescope", "calibrated", "cal"]
    file_root = input_rootname
    if file_root is None:
        for a_root in default_roots:
            if a_root in input_file:
                file_root = a_root
                break
        if file_root is None:
            err_msg = (
                f"Error in namefn_calibrated_input: input file {input_file}"
                f" matches none of the expected file roots: {default_roots}"
            )
            raise RuntimeError(err_msg)

    name_conv_dict = _get_name_conv_dict(file_root)

    conv = name_conv_dict[ap_filetype]
    output_file = input_file
    output_dir = conv["dir"]
    if conv["replace"] is not None:
        output_file = output_file.replace(conv["replace"], conv["with"])
    if conv["extension"] is not None:
        output_file = output_file.replace(input_suffix, conv["extension"])

    return output_file, output_dir


def namefn_getdir(ap_filetype: str) -> str:
    """
    Given a file type for star detection/astrometry, return the
    output directory such files would be placed in.

    For calibrated input files the the allowed output file types are:

    - ``input``: This is the calibrated input file itself.
    - ``srclist``: Star detection output FITS table file.
    - ``regfile``: ds9-format region file.
    - ``plotfile``: PNG plots of the input image with the detected star-like sources
      plotted as circles.
    - ``qualfile``: YaML source detection quality summary file.
    - ``fwhmplot``: PNG plot zooming in around a subset of the detected
      sources in the input image.
    - ``navfile``: A navigated and renamed copy of the input file, but now
      with a valid WCS solution.

    Parameters
    ----------
    ap_filetype : {'input', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile'}
        The AstroPhotography file type for which the directory should be returned.
        This should be one of the stage names described above.

    Returns
    -------
    output_dir : str
        Directory name in which output files will be written, relative
        to the root directory established by the input files.
    """

    allowed_stages = ["input", "srclist", "regfile", "plotfile", "qualfile", "fwhmplot", "navfile"]
    if ap_filetype not in allowed_stages:
        err_msg = (
            f"Requested ap_filetype {ap_filetype} not one of the allowed values: {allowed_stages}"
        )
        raise RuntimeError(err_msg)

    # True file root not needed
    file_root = "Calibrated-iTelescope"
    name_conv_dict = _get_name_conv_dict(file_root)
    conv = name_conv_dict[ap_filetype]
    output_dir = conv["dir"]
    return output_dir


def _get_name_conv_dict(file_root):
    """
    Utility function used by :func:`namefn_calibrated_input` that returns
    the file name and output directory dictionary given a file name root.

    Parameters
    ----------
    file_root : str
        The ``file_root`` is the unique string all input calibrated FITS files
        either start with, or include. For calibrated files produced by
        Astrophotography this might be the string `cal`, while iTelescope
        uses ``Calibrated`` or ``calibrated``.

    Returns
    -------
    name_conv_dict : dict
        Name conversion dictionary used to convert input calibrated file
        input names and paths into output file names.

    See Also
    --------
    :class:`ApProcess` : Batch processing of calibrated files
    """

    # Settings for output file names and output file paths
    # - This is difficult to fully automate without some assumptions about
    #   the input file names
    #   - I assume that all input file have the same file name prefix, e.g.
    #     Calibrated, or cal
    #   - I assume that all input files have the same file extension, e.g.
    #     fits or (sadness) fit
    # - I prefer to place output files of a certain type in a subdirectories.
    #   If you don't want all the diagnostic plots then it may just be simpler
    #   to not deal with subdirectories.
    # - If dir is not None then the various outputs files will be written to
    #   directories with the specified path relative to the **current** directory.
    #   The directory will be created if it not already present.
    name_conv_dict = {
        "input": {"replace": None, "with": None, "extension": None, "dir": "./"},
        "srclist": {
            "replace": file_root,
            "with": "srclist",
            "extension": ".fits",
            "dir": "./SourceLists/",
        },
        "regfile": {
            "replace": file_root,
            "with": "ds9",
            "extension": ".reg",
            "dir": "./SourceLists/",
        },
        "plotfile": {
            "replace": file_root,
            "with": "implot",
            "extension": ".png",
            "dir": "./SourceLists/",
        },
        "fwhmplot": {
            "replace": file_root,
            "with": "fwhmplot",
            "extension": ".png",
            "dir": "./SourceLists/",
        },
        "qualfile": {
            "replace": file_root,
            "with": "qual",
            "extension": ".yaml",
            "dir": "./MetaData/",
        },
        "navfile": {
            "replace": file_root,
            "with": "navigated",
            "extension": ".fits",
            "dir": "./NavigatedImages/",
        },
    }
    return name_conv_dict


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
        print(f"Warning: summarize_wcs() not written for 2-D images, {naxis}-dimensional data")

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
            val_str = "f{kw_str:8s} = {ipwcs.wcs.cdelt[irow]}"
            dothead_list.append(val_str)
            # Think that CD1_* = CDELT1 * PC1_* and CD2_* = CDELT2 * PC2_*
            pc = ipwcs.wcs.pc
            cdelt_vec = ipwcs.wcs.cdelt
            cd = pc.copy()
            for irow in range(ipwcs.naxis):
                cd[irow] = cdelt_vec[irow] * pc[irow]
            # then as above
        raise RuntimeError(
            "summarize_wcs needs to be updated to handle cases with a PC_ matrix and no CD_ matrix"
        )
    else:
        # Assume we have a simple CDELT[12] case with CROTA
        for irow in range(ipwcs.naxis):
            kw_str = f"CDELT{1+irow}"
            val_str = "f{kw_str:8s} = {ipwcs.wcs.cdelt[irow]}"
            dothead_list.append(val_str)
            kw_str = f"CROTA{1+irow}"
            val_str = "f{kw_str:8s} = {ipwcs.wcs.crota[irow]}"
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


def make_dothead_from_file(
    input_fits_with_wcs, output_swarp_dothead, extnum=0, format="fits", verbose=False
):
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

        if verbose:
            print(f"WCS created from primary header of {input_fits_with_wcs}")
            print(ipwcs)
            print(80 * "-")
            # print(ipwcs.wcs)
            # print(80*"-")

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
                print(f"--- dothead output from {input_fits_with_wcs} ---")
                print(dothead_str)
                print(f"--- about to write to {output_swarp_dothead} ---")

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
            olddata = hdulist[0].data
            data = 0 * olddata.astype(int)

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
