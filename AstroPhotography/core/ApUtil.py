# -*- coding: utf-8 -*-
#
#  Contains the implementation of ApMeasureStars
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

#-----------------------------------------------------------------------
# Should all be moved to a subpackage
#-----------------------------------------------------------------------

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

from astropy.io import fits
from astropy.table import QTable, Table, vstack
from astropy.coordinates import SkyCoord, Angle
from astropy.visualization import (AsymmetricPercentileInterval, 
                                   MinMaxInterval, 
                                   ManualInterval, 
                                   SqrtStretch, AsinhStretch, LinearStretch,
                                   ImageNormalize)
from astropy.visualization import make_lupton_rgb
        
# AstroPhotography includes    
from .. import __version__

def does_file_exists(filename, verbose=False):
    """
    Returns True if the file name or path exists, false otherwise
    """
    if verbose and not Path(filename).exists():
        print(f"Cannot find {filename}. Not a valid path or file.")
        return False
    else:
        print(f"Found {filename}.")
    return True

def load_image_and_plot(fname, extnum=0, usewcs=True, vmin=None, vmax=None, xaxlim=None, yaxlim=None, verbose=True):
    """
    Quick and dirty FITS image plot.

    By default `asinh` scaling will be used between a minimum and maximum data value. If vmin
    and/or vmnax are not specified then the 0.5th and 99.5th perciles will be used, as appropriate.
    A viridis color map is used by default, along with a color bar.

    Note that axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis
    pixels.

    TODO: Find out how to specify them in RA and Dec.

    :param fname: Name of existing FITS file with WCS in HDU number extnum that we want to plot.
    :param extnum: Extension number (zero-based) for data and WCS header
    :param usewcs: If True then plot using WCS information
    :param vmin: If not None then vmin is the minimum value in the normalaztion interval.
    :param vmax: If not None then vmax is the maximum value in the normalaztion interval.
    :param xaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True. 
    :param yaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True.
    :param verbose: If True then diagnostic information will be written to stdout.
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
        ax = fig.add_subplot(1, 1, 1, projection=w)
    else:
        ax = fig.add_subplot(1, 1, 1)
        
    im = ax.imshow(hdu.data, origin='lower', norm=norm)
    fig.colorbar(im, ax=ax)

    # Display default axis limits
    xlim_used = ax.get_xbound()
    ylim_used = ax.get_ybound()
    if verbose:
        print(f'Default X-axis limits: {xlim_used}')
        print(f'Default Y-axis limits: {ylim_used}')
    
    title_str = fname#.replace("_", "\_")
    ax.set_title(f'{title_str}', fontsize=8)
    if usewcs:
        ax.set_xlabel('Right Ascension (deg)', fontsize=8)
        ax.set_ylabel('Declination (deg)', fontsize=8)
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
    return
    
# First, lets just plot each resampled image separately
def load_three_images_and_plot(redfile, greenfile, bluefile, extnum=0, usewcs=True, vmin=None, vmax=None, xaxlim=None, yaxlim=None, verbose=True):
    """
    Quick and dirty FITS image plot.

    By default `asinh` scaling will be used between a minimum and maximum data value. If vmin
    and/or vmnax are not specified then the 0.5th and 99.5th perciles will be used, as appropriate.
    A viridis color map is used by default, along with a color bar.

    Note that axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis
    pixels.

    TODO: Find out how to specify them in RA and Dec.

    :param redfile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as red in the RGB color composite.
    :param greenfile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as green in the RGB color composite.
    :param bluefile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as blue in the RGB color composite.
    :param extnum: Extension number (zero-based) for data and WCS header
    :param usewcs: If True then plot using WCS information
    :param vmin: If not None then vmin is the minimum value in the normalaztion interval.
    :param vmax: If not None then vmax is the maximum value in the normalaztion interval.
    :param xaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True. 
    :param yaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True.
    :param verbose: If True then diagnostic information will be written to stdout.
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
        fig.colorbar(im, ax=ax)

        # Display default axis limits
        xlim_used = ax.get_xbound()
        ylim_used = ax.get_ybound()
        if verbose:
            print(f'Default X-axis limits: {xlim_used}')
            print(f'Default Y-axis limits: {ylim_used}')
        
        ax.set_title(f'{title_str}', fontsize=8)
        if usewcs:
            ax.set_xlabel('Right Ascension (deg)', fontsize=8)
            ax.set_ylabel('Declination (deg)', fontsize=8)
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
    return
    
def plot_lupton_threecolor(redfile, greenfile, bluefile, outpngfile, extnum=0, usewcs=True, vmin=None, xaxlim=None, yaxlim=None, qval=10, stretchval=0.5, verbose=True):
    """
    Quick and dirty 3-color composite using make_lupton_rgb

    Note: 
    - All images must share the same size, WCS orientation, and pixel size.
    - Axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis pixels.

    TODO: Find out how to specify them in RA and Dec.

    :param redfile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as red in the RGB color composite.
    :param greenfile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as green in the RGB color composite.
    :param bluefile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as blue in the RGB color composite.
    :param outpngfile: Name of PNG version of composite the generate.
    :param extnum: Extension number (zero-based) for data and WCS header
    :param usewcs: If True then plot using WCS information
    :param vmin: If not None then vmin is the minimum value the image, i.e. that will correspond to black
    :param xaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True. 
    :param yaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True.
    :param qval: Value of make_lupton_rgb Q parameter (Q in Lupton et al 2004)
    :param stretchval: Value of make_lupton_rgb stretch parameter (alpha in Lupton et al 2004)
    :param verbose: If True then diagnostic information will be written to stdout.
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
    fig.colorbar(im, ax=ax)

    # Display default axis limits
    xlim_used = ax.get_xbound()
    ylim_used = ax.get_ybound()
    if verbose:
        print(f'Default X-axis limits: {xlim_used}')
        print(f'Default Y-axis limits: {ylim_used}')
    
    title_str = f'RGB composite image stretch={stretchval}, Q={qval}\nRed = {redfile}\nGreen = {greenfile}\nBlue = {bluefile}\n'
    ax.set_title(f'{title_str}', fontsize=8)
    if usewcs:
        ax.set_xlabel('Right Ascension (deg)', fontsize=8)
        ax.set_ylabel('Declination (deg)', fontsize=8)
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
    return rgb_array
    
