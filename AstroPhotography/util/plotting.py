# -*- coding: utf-8 -*-
#
#  Contains the implementation of various Astrophotography plotting utilities
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
# SPDX-License-Identifier: GPL-3.0-or-later

# 2024-04-07 dks : issue-002 Initial coding
# 2025-03-01 dks : Split off from ApUtil.py

from typing import Any
import numpy as np
import numpy.typing as npt

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


def load_image_and_plot(
    fname: str,
    extnum: int | str = 0,
    output: str | None = None,
    usewcs: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
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

    # Compute vmin and vmax if necessary from percentiles
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
        assert hasattr(ax, "coords")
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
            fig.savefig(output, dpi=100, bbox_inches="tight")
            if verbose:
                print(f"Wrote plot to {output}")
    plt.close(fig=fig)
    return


def load_imlist_and_plot(
    imlist: list[str],
    extnum: int | str = 0,
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
                panel_rows = nrows_arr[idx]
                panel_cols = ncols_arr[idx]
                break
    print(
        (f"  With {num_imgs} files the plot layout is {panel_rows} rows x {panel_cols} columns.")
    )

    default_font_size = 7
    cbar_font_size = 5

    # Percentiles to use
    ipctls = [0.5, 99.5]

    # Display the image
    fig = plt.figure(figsize=(7.5, 10.0))
    fig.subplots_adjust(hspace=0.4)

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

            title_str: str = fname.split("/")[-1]
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
            ax.tick_params(labelsize=default_font_size)
            cbar = fig.colorbar(
                im,
                ax=ax,
                extend="neither",
                spacing="proportional",
                orientation="vertical",
                shrink=0.85,
            )
            cbar.set_label(r"Units TBA", fontsize=cbar_font_size)
            cbar.ax.tick_params(labelsize=cbar_font_size)

            # Display default axis limits
            xlim_used = ax.get_xbound()
            ylim_used = ax.get_ybound()
            if verbose:
                print(f"Default X-axis limits: {xlim_used}")
                print(f"Default Y-axis limits: {ylim_used}")

            ax.set_title(f"{title_str}", fontsize=default_font_size)
            if usewcs:
                # Adapted from @astrofrog at https://github.com/astropy/astropy/issues/13458#issuecomment-1242640539
                assert hasattr(ax, "coords")
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
            fig.savefig(output, dpi=100, bbox_inches="tight")
            if verbose:
                print(f"Wrote plot to {output}")
    plt.close(fig=fig)
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
    add_cbar: bool = False,
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

    Warnings
    --------
    Currently this function does not handle different minimum and maximum levels
    for the different channels, which results in terrible quality plots if the
    images have different background levels or are in different units.

    This will be addressed in a later version of ``AstroPhotography`` when the
    required base version of ``astropy``will be bumped from 6.0 to 7.0.

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
    add_cbar : bool, optional, default=False
        Add a scale bar to the plot. For the astropy lupton plot all inputs
        values are scaled to the range ``[0, 255]`` so the scale bar is
        superfluous.

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

    default_font_size = 6
    cbar_font_size = 5

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

    vmax = None  # vmax is not an input function parameter
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

    ax.tick_params(labelsize=default_font_size)
    im = ax.imshow(rgb_array, origin="lower")
    if add_cbar:
        cbar = fig.colorbar(
            im,
            ax=ax,
            extend="neither",
            spacing="proportional",
            orientation="vertical",
            shrink=0.85,
        )
        cbar.set_label(r"Units TBA", fontsize=cbar_font_size)  # TODO get real data units
        cbar.ax.tick_params(labelsize=cbar_font_size)

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

    output = outpngfile

    if not panel_plotter and output is not None:
        if "interactive" in output:
            print("Processing blocked while plotting window is open. Close it to procede.")
            fig.show()
        else:
            fig.savefig(output, dpi=100, bbox_inches="tight")
            if verbose:
                print(f"Wrote plot to {output}")

    hdulist_r.close()
    hdulist_g.close()
    hdulist_b.close()
    return rgb_array


def plot_stiff_threecolor(
    redfile: str,
    greenfile: str,
    bluefile: str,
    outpngfile: str,
    rasterizer: Any,
    ax: Any | None = None,
    extnum: int | str = 0,
    min_cut: list[float] | None = None,
    max_cut: list[float] | None = None,
    min_percent: list[float] | None = None,
    max_percent: list[float] | None = None,
    negative: bool = False,
    binning: int = 1,
    xaxlim=None,
    yaxlim=None,
    verbose: bool = True,
):
    """
    Quick and dirty 3-color composite using ApFitsRasterizer and stiff

    - All images must share the same size, WCS orientation, and pixel size.
    - Axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis pixels.

    Note
    ----

    Unlike :func:`plot_lupton_threecolor` there is no mechanism to show WCS
    coordinates

    Axis limits (xaxlim and yaxlim) are defined in image x-axis and y-axis
    pixels.

    See Also
    --------
    :class:`ApFitsRasterizer`

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
    outpngfile : str
        Output file name for the 3-color composite plot generated by
        ``ApFitsRasterizer`` using ``stiff``.
        Specify a filename with recognizable image-format extension, e.g.
        ``my_output_plot.png`` to directly save the plot to a file.
    rasterizer : ApFitsRasterizer
        An existing ``ApFitsRasterizer`` instance to use.
    ax : matplotlib.Axes instance or None, optional, default=None
        This function supports being called from another function controlling
        a figure and set of Axes instances. If ax is not None then this function
        will not generated its own figure but instead use the specified Axes.
        Note that in such a case outpngfile should be None, **and** the
        caller must have created the Axes with the correct WCS project if
        WCS axes plotting is wanted.
    extnum : int or str, optional, default=0
        Extension number (zero-based) or extension name (string) for data and WCS header
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

    Returns
    -------
    rgb_array : ndarray
        The 3-dimensional data array produced by astropy's make_lupton_rgb function.

    """

    imlist = [redfile, greenfile, bluefile]
    rasterizer.fits_to_rgb(
        imlist,
        outpngfile,
        extnum=extnum,
        min_cut=min_cut,
        max_cut=max_cut,
        min_percent=min_percent,
        max_percent=max_percent,
        negative=negative,
        binning=binning,
    )

    default_font_size = 6
    cbar_font_size = 5

    rgb_array = plt.imread(outpngfile)

    # Display the image in a set of matplotlib Axes if necessary
    if ax is None:
        # Nothing to do
        pass
    else:
        # Also plot into externally provided matplotlib Axes, don't create a PNG ourselves
        fig = plt.gcf()

        ax.tick_params(labelsize=default_font_size)
        im = ax.imshow(rgb_array, origin="lower")

        # Display default axis limits
        xlim_used = ax.get_xbound()
        ylim_used = ax.get_ybound()
        if verbose:
            print(f"Default X-axis limits: {xlim_used}")
            print(f"Default Y-axis limits: {ylim_used}")

        title_str = (
            f"RGB composite image {outpngfile}"
            f"\nRed = {redfile}"
            f"\nGreen = {greenfile}"
            f"\nBlue = {bluefile}\n"
        )
        ax.set_title(f"{title_str}", fontsize=default_font_size)

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
    add_cbar: bool = False,
):
    """
    Create a single summary plot containing any RGB and/or SHO composite plots
    that can be generated from a set of images in multiple filters.

    Given a dictionary of ``dict{filter_name, image_name}`` this function will
    create a single panel color composite using :func:`plot_lupton_threecolor`
    for each of the following **complete** filter sets:

    - ``Red``, ``Green``, and ``Blue``: plotted in that order.
    - ``Ha``, ``Green``, and ``Blue``: plotted in that order.
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
    add_cbar : bool, optional, default=False
        Add a scale bar to the plot. For the astropy lupton plot all inputs
        values are scaled to the range ``[0, 255]`` so the scale bar is
        superfluous.
    """

    panel_cols: int = 1
    panel_rows: int = 0

    # Look for RGB, HGB, or SHO

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

    hgb_satisfied: bool = False
    hgb_list: list[str] = []
    hgb_required: list[str] = ["Ha", "Green", "Blue"]
    for filter_name in hgb_required:
        if filter_name in filter_imdict.keys():
            hgb_list.append(filter_imdict[filter_name])
    if len(hgb_list) == 3:
        hgb_satisfied = True
        panel_rows += 1
        print("There are enough files to generate a HaGB threecolor composite.")
    else:
        print("NB: insufficient files to generate a HaGB threecolor composite.")

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

    if not (rgb_satisfied or hgb_satisfied or sho_satisfied):
        print("Returning as cannot generated either RGB, HaGB, or SHO composities.")
        return

    # Use our own figure and axes
    fig = plt.figure(figsize=(7.5, 10.0))

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
            add_cbar=add_cbar,
        )
        hdulist.close()

    if hgb_satisfied:
        iplot += 1
        hdulist = fits.open(hgb_list[0])
        hdu1 = hdulist[extnum]

        if usewcs:
            w = wcs.WCS(hdu1.header)
            ax = fig.add_subplot(panel_rows, panel_cols, iplot, projection=w)
        else:
            ax = fig.add_subplot(panel_rows, panel_cols, iplot)

        _ = plot_lupton_threecolor(
            hgb_list[0],
            hgb_list[1],
            hgb_list[2],
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
            add_cbar=add_cbar,
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
            fig.savefig(output, dpi=100, bbox_inches="tight")
            if verbose:
                print(f"Wrote plot to {output}")
    plt.close(fig=fig)
    return


def make_stiff_threecolor_plots(
    filter_imdict: dict[str, str],
    outpngfile: str,
    rasterizer: Any,
    extnum: int | str = 0,
    min_cut: list[float] | None = None,
    max_cut: list[float] | None = None,
    min_percent: list[float] | None = None,
    max_percent: list[float] | None = None,
    negative: bool = False,
    binning: int = 1,
    xaxlim=None,
    yaxlim=None,
    verbose: bool = True,
):
    """
    Create a single summary plot containing any RGB and/or SHO composite plots
    that can be generated from a set of images in multiple filters.

    Given a dictionary of ``dict{filter_name, image_name}`` this function will
    create a single panel color composite using :func:`plot_lupton_threecolor`
    for each of the following **complete** filter sets:

    - ``Red``, ``Green``, and ``Blue``: plotted in that order.
    - ``Ha``, ``Green``, and ``Blue``: plotted in that order.
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
    rasterizer : ApFitsRasterizer, optional, default=None
        An existing ``ApFitsRasterizer`` instance to use. If None then
        a new instance will be created with ``INFO`` level logging.
    extnum : int or str, optional, default=0
        Extension number (zero-based) or extension name (string) for data and WCS header
    extnum : int or str, optional, default=0
        Extension number (zero-based) or extension name (string) for data and WCS header
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
    """

    panel_cols: int = 1
    panel_rows: int = 0

    # Look for RGB, HGB, or SHO

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

    hgb_satisfied: bool = False
    hgb_list: list[str] = []
    hgb_required: list[str] = ["Ha", "Green", "Blue"]
    for filter_name in hgb_required:
        if filter_name in filter_imdict.keys():
            hgb_list.append(filter_imdict[filter_name])
    if len(hgb_list) == 3:
        hgb_satisfied = True
        panel_rows += 1
        print("There are enough files to generate a HaGB threecolor composite.")
    else:
        print("NB: insufficient files to generate a HaGB threecolor composite.")

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

    if not (rgb_satisfied or hgb_satisfied or sho_satisfied):
        print("Returning as cannot generated either RGB, HaGB, or SHO composities.")
        return

    # Use our own figure and axes
    fig = plt.figure(figsize=(7.5, 10.0))

    iplot: int = 0
    if rgb_satisfied:
        iplot += 1
        ax = fig.add_subplot(panel_rows, panel_cols, iplot)

        _ = plot_stiff_threecolor(
            rgb_list[0],
            rgb_list[1],
            rgb_list[2],
            outpngfile="rgb.tiff",
            rasterizer=rasterizer,
            ax=ax,
            extnum=extnum,
            min_cut=min_cut,
            max_cut=max_cut,
            min_percent=min_percent,
            max_percent=max_percent,
            negative=negative,
            binning=binning,
            xaxlim=xaxlim,
            yaxlim=yaxlim,
            verbose=verbose,
        )

    if hgb_satisfied:
        iplot += 1
        ax = fig.add_subplot(panel_rows, panel_cols, iplot)

        _ = plot_stiff_threecolor(
            hgb_list[0],
            hgb_list[1],
            hgb_list[2],
            outpngfile="hgb.tiff",
            rasterizer=rasterizer,
            ax=ax,
            extnum=extnum,
            min_cut=min_cut,
            max_cut=max_cut,
            min_percent=min_percent,
            max_percent=max_percent,
            negative=negative,
            binning=binning,
            xaxlim=xaxlim,
            yaxlim=yaxlim,
            verbose=verbose,
        )

    if sho_satisfied:
        iplot += 1
        ax = fig.add_subplot(panel_rows, panel_cols, iplot)

        _ = plot_stiff_threecolor(
            sho_list[0],
            sho_list[1],
            sho_list[2],
            outpngfile="sho.tiff",
            rasterizer=rasterizer,
            ax=ax,
            extnum=extnum,
            min_cut=min_cut,
            max_cut=max_cut,
            min_percent=min_percent,
            max_percent=max_percent,
            negative=negative,
            binning=binning,
            xaxlim=xaxlim,
            yaxlim=yaxlim,
            verbose=verbose,
        )

    output = outpngfile
    if output is not None:
        if "interactive" in output:
            print("Processing blocked while plotting window is open. Close it to procede.")
            fig.show()
        else:
            fig.savefig(output, dpi=100, bbox_inches="tight")
            if verbose:
                print(f"Wrote plot to {output}")
    plt.close(fig=fig)
    return
