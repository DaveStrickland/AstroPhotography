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

import logging
import os.path
import pathlib
import numpy as np
import matplotlib                # for rc
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import math
import time
from datetime import datetime

from astropy import wcs
import astropy
from astropy.io import fits
from astropy.table import QTable, Table, vstack
from astropy.coordinates import SkyCoord, Angle
from astropy import units
from astropy.visualization import (AsymmetricPercentileInterval, 
                                   MinMaxInterval, 
                                   ManualInterval, 
                                   SqrtStretch, AsinhStretch, LinearStretch,
                                   ImageNormalize)
from astropy.visualization import make_lupton_rgb
        
# AstroPhotography includes    
from .. import __version__

def does_file_exist(filename, verbose=False):
    """
    Returns True if the file name or path exists, false otherwise
    """
    if verbose and not Path(filename).exists():
        print(f"Cannot find {filename}. Not a valid path or file.")
        return False
    else:
        print(f"Found {filename}.")
    return True

def load_image_and_plot(fname, extnum=0, output=None, 
                        usewcs=True, vmin=None, vmax=None, 
                        xaxlim=None, yaxlim=None, verbose=True,
                        angle_tick_spacing_am=2.0, swap_radec_axis=False):
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
        w   = wcs.WCS(hdu.header)

    # Compute vmin and vmax if necessary
    # percentiles
    ipctls = [0.5, 99.5]
    opctls = np.nanpercentile(hdu.data, ipctls)
    ourmin = opctls[0]
    ourmax = opctls[1]
    if verbose:
        print(f'Input data 0.5th percentile = {ourmin:.4f}, 99.5th percentile = {ourmax:.4f}')
    
    if vmin is not None:
        ourmin = vmin
    if vmax is not None:
        ourmax = vmax
    if verbose:
        print(f'Applying asinh stretch between {ourmin:.4f} and {ourmax:.4f}')
    
    # Create an ImageNormalize object
    norm = ImageNormalize(hdu.data, interval=ManualInterval(ourmin, ourmax),
                      stretch=AsinhStretch())

    # Display the image
    fig = plt.figure()
    if usewcs:
        if verbose:
            print('Using WCS information for axis projection')
            print(w)
        ax = fig.add_subplot(1, 1, 1, projection=w)
    else:
        ax = fig.add_subplot(1, 1, 1)
        
    im = ax.imshow(hdu.data, origin='lower', norm=norm)
    cbar = fig.colorbar(im, ax=ax, extend='neither', spacing='proportional',
                orientation='vertical', shrink=0.85)
    cbar.set_label(r"Units TBA")
    cbar.ax.tick_params(labelsize=8) 
    title_str = fname#.replace("_", "\_")
    ax.set_title(f'{title_str}', fontsize=8)

    # Display default axis limits
    xlim_used = ax.get_xbound()
    ylim_used = ax.get_ybound()
    if verbose:
        print(f'Default X-axis limits: {xlim_used}')
        print(f'Default Y-axis limits: {ylim_used}')
    
    title_str = fname#.replace("_", "\_")
    ax.set_title(f'{title_str}', fontsize=8)
    if usewcs:
        # Adapted rom @astrofrog at https://github.com/astropy/astropy/issues/13458#issuecomment-1242640539
        ra = ax.coords[0]
        dec = ax.coords[1]
        ra.set_major_formatter('hh:mm:ss.ss')   # RA in Hours, minutes, seconds,
        dec.set_major_formatter('dd:mm:ss.ss')
        ra.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color='red')      # Ticks must be defined for axis tickvals to appear
        dec.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color='blue')
        ra.set_ticklabel(color='red', fontsize=8)
        dec.set_ticklabel(color='blue', fontsize=8)
        ra.grid(color='red', linestyle='--', alpha=0.6)
        dec.grid(color='blue', linestyle='--', alpha=0.6)
        ra.set_axislabel('Right Ascension (HMS)', fontsize=8, color='red')
        dec.set_axislabel('Declination (dms)', fontsize=8, color='blue')

        if swap_radec_axis:
            # Images where ra changes fastest on Y, not X
            dec.set_ticks_position('b')
            dec.set_ticklabel_position('b')
            dec.set_axislabel_position('b')
            ra.set_ticks_position('l')
            ra.set_ticklabel_position('l')
            ra.set_axislabel_position('l')
    
        #ax.set_xlabel('Right Ascension (deg)', fontsize=8)
        #ax.set_ylabel('Declination (deg)', fontsize=8)
    else:
        ax.set_xlabel('X-axis pixel number', fontsize=8)
        ax.set_ylabel('Y-axis pixel number', fontsize=8)

    # Modify the axis limits?
    if xaxlim is not None:
        loval = np.min(xaxlim)
        hival = np.max(xaxlim)
        ax.set_xbound(lower=loval, upper=hival)
        if verbose:
            print(f'Setting X-axis limits to [{loval}, {hival}]')
    if yaxlim is not None:
        loval = np.min(yaxlim)
        hival = np.max(yaxlim)
        ax.set_ybound(lower=loval, upper=hival)
        if verbose:
            print(f'Setting Y-axis limits to [{loval}, {hival}]')
            
    if output is not None:
        if 'interactive' in output:
            print('Processing blocked while plotting window is open. Close it to procede.')
            fig.show()
        else:
            fig.savefig(output, dpi=200, bbox_inches='tight')
            if verbose:
                print(f'Wrote plot to {output}')
    return
    
# First, lets just plot each resampled image separately
def load_three_images_and_plot(redfile, greenfile, bluefile, extnum=0, 
                        output=None, usewcs=True, vmin=None, vmax=None, xaxlim=None, yaxlim=None, verbose=True,
                        angle_tick_spacing_am=2.0, swap_radec_axis=False):
    """
    Quick and dirty FITS image plot.

    By default `asinh` scaling will be used between a minimum and maximum data value. If vmin
    and/or vmnax are not specified then the 0.5th and 99.5th perciles will be used, as appropriate.
    A viridis color map is used by default, along with a color bar.

    Note that axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis
    pixels.

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
    redfile : str
        Name of existing FITS file with WCS in HDU number extnum that we want to
        appear as red in the RGB color composite.
    greenfile : str
        Name of existing FITS file with WCS in HDU number extnum that we want to
        appear as green in the RGB color composite.
    bluefile : str
        Name of existing FITS file with WCS in HDU number extnum that we want to
        appear as blue in the RGB color composite.
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
    
    hdur = fits.open(redfile)[extnum]
    hdug = fits.open(greenfile)[extnum]
    hdub = fits.open(bluefile)[extnum]

    # Compute vmin and vmax if necessary. NOTE this only uses the red channel
    # percentiles
    ipctls = [0.5, 99.5]
    opctls = np.nanpercentile(hdur.data, ipctls)
    ourmin = opctls[0]
    ourmax = opctls[1]
    if verbose:
        print(f'Input red channel data 0.5th percentile = {ourmin:.4f}, 99.5th percentile = {ourmax:.4f}')
    
    if vmin is not None:
        ourmin = vmin
    if vmax is not None:
        ourmax = vmax
    if verbose:
        print(f'Applying asinh stretch between {ourmin:.4f} and {ourmax:.4f}')
    
    # Display the image
    fig = plt.figure()

    for idx in range(3):
        if idx == 0:
            hdu       = hdur
            title_str = redfile
            w         = wcs.WCS(hdu.header)
            print(80*'-')
            print(f'WCS for {title_str}:')
            print(w)
        elif idx == 1:
            hdu       = hdug
            title_str = greenfile
            w         = wcs.WCS(hdu.header)
            print(80*'-')
            print(f'WCS for {title_str}:')
            print(w)
        elif idx == 2:
            hdu       = hdub
            title_str = bluefile
            w         = wcs.WCS(hdu.header)
            print(80*'-')
            print(f'WCS for {title_str}:')
            print(w)
        else:
            raise RuntimeError(f'Error, index {idx} should be in range 0..2 for 3-image plotting.')

        if usewcs:
            ax = fig.add_subplot(2, 2, idx+1, projection=w)
        else:
            ax = fig.add_subplot(2, 2, idx+1)

        # Create an ImageNormalize object
        norm = ImageNormalize(hdu.data, interval=ManualInterval(ourmin, ourmax),
                      stretch=AsinhStretch())
        
        im = ax.imshow(hdu.data, origin='lower', norm=norm)
        cbar = fig.colorbar(im, ax=ax, extend='neither', spacing='proportional',
                orientation='vertical', shrink=0.85)
        cbar.set_label(r"Units TBA")
        cbar.ax.tick_params(labelsize=8) 
    
        # Display default axis limits
        xlim_used = ax.get_xbound()
        ylim_used = ax.get_ybound()
        if verbose:
            print(f'Default X-axis limits: {xlim_used}')
            print(f'Default Y-axis limits: {ylim_used}')
        
        ax.set_title(f'{title_str}', fontsize=8)
        if usewcs:
            # Adapted rom @astrofrog at https://github.com/astropy/astropy/issues/13458#issuecomment-1242640539
            ra = ax.coords[0]
            dec = ax.coords[1]
            ra.set_major_formatter('hh:mm:ss.ss')   # RA in Hours, minutes, seconds,
            dec.set_major_formatter('dd:mm:ss.ss')
            ra.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color='red')      # Ticks must be defined for axis tickvals to appear
            dec.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color='blue')
            ra.set_ticklabel(color='red', fontsize=8)
            dec.set_ticklabel(color='blue', fontsize=8)
            ra.grid(color='red', linestyle='--', alpha=0.6)
            dec.grid(color='blue', linestyle='--', alpha=0.6)
            ra.set_axislabel('Right Ascension (HMS)', fontsize=8, color='red')
            dec.set_axislabel('Declination (dms)', fontsize=8, color='blue')

            if swap_radec_axis:
                # Images where ra changes fastest on Y, not X
                dec.set_ticks_position('b')
                dec.set_ticklabel_position('b')
                dec.set_axislabel_position('b')
                ra.set_ticks_position('l')
                ra.set_ticklabel_position('l')
                ra.set_axislabel_position('l')
        
            #ax.set_xlabel('Right Ascension (deg)', fontsize=8)
            #ax.set_ylabel('Declination (deg)', fontsize=8)
        else:
            ax.set_xlabel('X-axis pixel number', fontsize=8)
            ax.set_ylabel('Y-axis pixel number', fontsize=8)
    
        # Modify the axis limits?
        if xaxlim is not None:
            loval = np.min(xaxlim)
            hival = np.max(xaxlim)
            ax.set_xbound(lower=loval, upper=hival)
            if verbose:
                print(f'Setting X-axis limits to [{loval}, {hival}]')
        if yaxlim is not None:
            loval = np.min(yaxlim)
            hival = np.max(yaxlim)
            ax.set_ybound(lower=loval, upper=hival)
            if verbose:
                print(f'Setting Y-axis limits to [{loval}, {hival}]')
    
    if output is not None:
        if 'interactive' in output:
            print('Processing blocked while plotting window is open. Close it to procede.')
            fig.show()
        else:
            fig.savefig(output, dpi=200, bbox_inches='tight')
            if verbose:
                print(f'Wrote plot to {output}')
    return
    
def plot_lupton_threecolor(redfile, greenfile, bluefile, outpngfile, 
                        extnum=0, usewcs=True, vmin=None, xaxlim=None, yaxlim=None, qval=10, stretchval=0.5, verbose=True,
                        angle_tick_spacing_am=2.0, swap_radec_axis=False):
    """
    Quick and dirty 3-color composite using make_lupton_rgb

    - All images must share the same size, WCS orientation, and pixel size.
    - Axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis pixels.

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
        

    redfile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as red in the RGB color composite.
    greenfile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as green in the RGB color composite.
    bluefile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as blue in the RGB color composite.
    outpngfile: Name of PNG version of composite the generate.
    extnum: Extension number (zero-based) for data and WCS header
    usewcs: If True then plot using WCS information
    vmin: If not None then vmin is the minimum value the image, i.e. that will correspond to black
    xaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True. 
    yaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True.
    qval: Value of make_lupton_rgb Q parameter (Q in Lupton et al 2004)
    stretchval: Value of make_lupton_rgb stretch parameter (alpha in Lupton et al 2004)
    verbose: If True then diagnostic information will be written to stdout.
    """
    
    hdur = fits.open(redfile)[extnum]
    hdug = fits.open(greenfile)[extnum]
    hdub = fits.open(bluefile)[extnum]
    w    = wcs.WCS(hdur.header)

    # Compute vmin and vmax if necessary. NOTE this only uses the red channel
    # percentiles
    ipctls = [0.5, 99.5]
    opctls = np.nanpercentile(hdur.data, ipctls)
    ourmin = opctls[0]
    ourmax = opctls[1]
    if verbose:
        print(f'Input red channel data 0.5th percentile = {ourmin:.4f}, 99.5th percentile = {ourmax:.4f}')
    
    vmax = None
    if vmin is not None:
        ourmin = vmin
    if vmax is not None:
        ourmax = vmax
    if verbose:
        print(f'Using an image minimum of {ourmin} for all channels.')

    rgb_array = make_lupton_rgb(hdur.data, hdug.data, hdub.data,
                               minimum=ourmin, stretch=stretchval, Q=qval, filename=outpngfile)
    
    # Display the image
    fig = plt.figure()
    if usewcs:
        ax = fig.add_subplot(1, 1, 1, projection=w)
    else:
        ax = fig.add_subplot(1, 1, 1)
        
    im = ax.imshow(rgb_array, origin='lower')
    cbar = fig.colorbar(im, ax=ax, extend='neither', spacing='proportional',
                orientation='vertical', shrink=0.85)
    cbar.set_label(r"Units TBA")
    cbar.ax.tick_params(labelsize=8) 

    # Display default axis limits
    xlim_used = ax.get_xbound()
    ylim_used = ax.get_ybound()
    if verbose:
        print(f'Default X-axis limits: {xlim_used}')
        print(f'Default Y-axis limits: {ylim_used}')
    
    title_str = f'RGB composite image stretch={stretchval}, Q={qval}\nRed = {redfile}\nGreen = {greenfile}\nBlue = {bluefile}\n'
    ax.set_title(f'{title_str}', fontsize=8)
    if usewcs:
        # Adapted rom @astrofrog at https://github.com/astropy/astropy/issues/13458#issuecomment-1242640539
        ra = ax.coords[0]
        dec = ax.coords[1]
        ra.set_major_formatter('hh:mm:ss.ss')   # RA in Hours, minutes, seconds,
        dec.set_major_formatter('dd:mm:ss.ss')
        ra.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color='red')      # Ticks must be defined for axis tickvals to appear
        dec.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color='blue')
        ra.set_ticklabel(color='red', fontsize=8)
        dec.set_ticklabel(color='blue', fontsize=8)
        ra.grid(color='red', linestyle='--', alpha=0.6)
        dec.grid(color='blue', linestyle='--', alpha=0.6)
        ra.set_axislabel('Right Ascension ((HMS))', fontsize=8, color='red')
        dec.set_axislabel('Declination (dms)', fontsize=8, color='blue')

        if swap_radec_axis:
            # Images where ra changes fastest on Y, not X
            dec.set_ticks_position('b')
            dec.set_ticklabel_position('b')
            dec.set_axislabel_position('b')
            ra.set_ticks_position('l')
            ra.set_ticklabel_position('l')
            ra.set_axislabel_position('l')
    
        #ax.set_xlabel('Right Ascension (deg)', fontsize=8)
        #ax.set_ylabel('Declination (deg)', fontsize=8)
    else:
        ax.set_xlabel('X-axis pixel number', fontsize=8)
        ax.set_ylabel('Y-axis pixel number', fontsize=8)

    # Modify the axis limits?
    if xaxlim is not None:
        loval = np.min(xaxlim)
        hival = np.max(xaxlim)
        ax.set_xbound(lower=loval, upper=hival)
        if verbose:
            print(f'Setting X-axis limits to [{loval}, {hival}]')
    if yaxlim is not None:
        loval = np.min(yaxlim)
        hival = np.max(yaxlim)
        ax.set_ybound(lower=loval, upper=hival)
        if verbose:
            print(f'Setting Y-axis limits to [{loval}, {hival}]')
    
    output = outpngfile        
    if output is not None:
        if 'interactive' in output:
            print('Processing blocked while plotting window is open. Close it to procede.')
            fig.show()
        else:
            fig.savefig(output, dpi=200, bbox_inches='tight')
            if verbose:
                print(f'Wrote plot to {output}')
    return rgb_array
        
def namefn_calibrated_input(input_file, input_rootname, input_suffix, ap_filetype):
    """
    Given the name of a calibrated input file, generate the file name for
    and output directory for one of the subsequent output file types.
    
    For calibrated input files the the allowed output file types are:
    
    - ``inputs``: This is the calibrated input file itself.
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
    ap_filetype : {'inputs', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile'}
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
    This function assumes that the input files is not in a subdirectory,
    i.e. there is no path within in ``input_file`` string. In the longer
    term this should be rewritten using the pathlib module without such
    an assumption.
    """
    
    allowed_stages = ['inputs', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile']
    if ap_filetype not in allowed_stages:
        err_msg = f'Requested ap_filetype {ap_filetype} not one of the allowed values: {allowed_stages}'
        self._logger.error(err_msg)
        raise RuntimeError(err_msg)
    
    # Determine the root name to use.
    default_roots = ['Calibrated-iTelescope', 'calibrated', 'cal']
    file_root = input_rootname
    if file_root is None:
        for a_root in default_roots:
            if a_root in input_file:
                file_root = a_root
                break
        if file_root is None:
            err_msg = f'Error in namefn_calibrated_input: input file {input_file} matches none of the expected file roots: {default_roots}'
            raise RuntimeError(err_msg)
                
    name_conv_dict = _get_name_conv_dict(file_root)

    conv = name_conv_dict[ap_filetype]        
    output_file = input_file
    output_dir  = conv['dir']
    if conv['replace'] is not None:
            output_file = output_file.replace(conv['replace'], conv['with'])
    if conv['extension'] is not None:
            output_file = output_file.replace(input_suffix, conv['extension'])
    
    return output_file, output_dir

def namefn_getdir(ap_filetype):
    """
    Given a file type for star detection/astrometry, return the
    output directory such files would be placed in.
    
    For calibrated input files the the allowed output file types are:
    
    - ``inputs``: This is the calibrated input file itself.
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
    ap_filetype : {'inputs', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile'}
        The AstroPhotography file type for which the directory should be returned.
        This should be one of the stage names described above.
               
    Returns
    -------
    output_dir : str
        Directory name in which output files will be written, relative
        to the root directory established by the input files.
    """
    
    allowed_stages = ['inputs', 'srclist', 'regfile', 'plotfile', 'qualfile', 'fwhmplot', 'navfile']
    if ap_filetype not in allowed_stages:
        err_msg = f'Requested ap_filetype {ap_filetype} not one of the allowed values: {allowed_stages}'
        self._logger.error(err_msg)
        raise RuntimeError(err_msg)
    
    # True file root not needed
    file_root = 'Calibrated-iTelescope'
    name_conv_dict = _get_name_conv_dict(file_root)
    conv = name_conv_dict[ap_filetype]        
    output_dir  = conv['dir']
    return output_dir

def _get_name_conv_dict(file_root):
    """
    Utility function used by :func:`namefn_calibrated_input` that returns 
    the file name and output directory dictionary given a file name root.
    """
    
    # Settings for output file names and output file paths
    # - This is difficult to fully automate without some assumptions about the input file names
    #   - I assume that all input file have the same file name prefix, e.g. Calibrated, or cal
    #   - I assume that all input files have the same file extension, e.g. fits or (sadness) fit
    # - I prefer to place output files of a certain type in a subdirectories. If you don't want all
    #   the diagnostic plots then it may just be simpler to not deal with subdirectories.
    # - If dir is not None then the various outputs files will be written to directories with the specified
    #   path relative to the **current** directory. The directory will be created if it not already
    #   present.
    name_conv_dict = {'inputs':   {'replace': None,      'with': None,        'extension': None,    'dir': './'},
                      'srclist':  {'replace': file_root, 'with': 'srclist',   'extension': '.fits', 'dir': './SourceLists/'},
                      'regfile':  {'replace': file_root, 'with': 'ds9',       'extension': '.reg',  'dir': './SourceLists/'},
                      'plotfile': {'replace': file_root, 'with': 'implot',    'extension': '.png',  'dir': './SourceLists/'},
                      'fwhmplot': {'replace': file_root, 'with': 'fwhmplot',  'extension': '.png',  'dir': './SourceLists/'},
                      'qualfile': {'replace': file_root, 'with': 'qual',      'extension': '.yaml', 'dir': './MetaData/'},
                      'navfile':  {'replace': file_root, 'with': 'navigated', 'extension': '.fits', 'dir': './NavigatedImages/'}}
    return name_conv_dict
