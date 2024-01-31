# -*- coding: utf-8 -*-
#
#  Contains the implementation of ApMeasureStars
#  
#  Copyright 2020-2021 Dave Strickland <dave.strickland@gmail.com>
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

# 2021-01-18 dks : Move ApMeasureStars into core in separate file.
# 2024-01-25 dks : Catch up to latest astropy/photutils changes

import logging
import os.path
import numpy as np
import matplotlib                # for rc
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import math
import time
import yaml
from datetime import datetime

from scipy import spatial

from astropy.io import fits
from astropy.table import QTable, Table, vstack
from astropy.coordinates import SkyCoord, Angle
from astropy.visualization import AsymmetricPercentileInterval, MinMaxInterval, ManualInterval
from astropy.visualization import SqrtStretch, AsinhStretch, LinearStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import sigma_clipped_stats
from astropy import units as u
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip, mad_std

from regions import PixCoord, CirclePixelRegion

from photutils import find_peaks, DAOStarFinder
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
        
class ApMeasureStars:
    """
    Measures stellar PSF size and symmetry in selected stars across
    and input image.
    
    The primary purpose of this class is to measure the size (full width
    at half maximum) of stars in the input image and sourcelist, and to
    attempt to detect if there is significant asymetry in the star images.
    This information can be used to:
    
    #. Measure the seeing (star FWHM) in the image.
    #. Refine source searching with ApFindStars based on the "true" FWHM,
    #. Identify images with poor seeing or tracking in comparison to
       other images of the same target from the same telescope. Such 
       images should be excluded from image stacking.
    
    This class uses astropy.modelling to fit 2-D gaussian plus a constant
    background model to cut-out regions around a selected subset of the
    stars specified in the input source list.
    
    A major current limitation is quantifying the errors in the fit (not
    currently done). Astropy modelling uses scipy under the hood, which
    has a very limited capability of expressing or investigating fit
    uncertainties.
    """
    
    # Threshold in number of sigma that fwhm_y must differ from
    # fwhm_x in order to be counted as NOT being circular.
    _circ_thresh_sigma = 3.0
    
    def __init__(self,
            img_data,
            srclist,
            init_fwhm,
            init_bglevel,
            full_srclist,
            fwhm_plot_file,
            fwhm_plot_title,
            loglevel,
            quiet):
        """
        ApMeasureStars constructor
        """
        
        self._img_data    = img_data        # Image data 2-D array
        self._init_fwhm   = init_fwhm       # Initial guess at FWHM in pixels
        self._init_bglvl  = init_bglevel    # Initial guess at BG level per pixel
        self._loglevel    = loglevel        # Logging level
        self._quiet       = quiet           # If True then table summary printing disabled.

        self._fit_table   = None            # Table to store fit results in.
        self._pixel_array = None            # 3-D array for star cut-outs
        self._fwhm_plot   = fwhm_plot_file  # None or PNG/JPG image file name.
        self._plot_title  = fwhm_plot_title # None or string to use in fwhm_plot_file
            
        # Set up logging
        self._logger = self._initialize_logger(self._loglevel)

        # Settings related to candidate selection.
        self._use_weights    = True
        self._num_per_reg    = 5           # Number of sources per region to fit
        self._skip_brightest = 0           # Skip the brightest N stars in each region
        self._logger.debug(f'Up to {self._num_per_reg} stars per sub-region will be fitted, excluding the brightest {self._skip_brightest} stars.')

        # Settings related to source fitting.
        self._use_weights    = True
        self._fit_for_bg     = True
        self._logger.debug(f'Count based weighting factors will used in fitting?: {self._use_weights}')
        self._logger.debug(f'Will the background level be fitted for?           : {self._fit_for_bg}')

        # POsitions for all possible sources in the image
        self._full_srcs = full_srclist

        # Slightly modify the input source list to exclude possibly
        # saturated stars, remove unneeded columns.
        self._init_srcs  = Table(srclist[srclist['psbl_sat']==False])         # Initial sourclist from ApFindStars
        self._init_srcs.remove_columns(['aperture_sum', 'psbl_sat', 'adu_per_sec'])
        self._full_srcs.remove_columns(['aperture_sum', 'psbl_sat', 'adu_per_sec'])
        if not self._quiet:
            print('Input table supplied to ApMeasureStars after saturated star filtering:')
            print(self._init_srcs)
            print('')
        self._logger.info(f'Size of input trimmed source list (filtered): {len(self._init_srcs)}')
        self._logger.info(f'Size of full source list used for neighbor removal: {len(self._full_srcs)}')

        self._cols       = img_data.shape[1]
        self._rows       = img_data.shape[0]
        
        # Set up variables related to fit box size and edge exclusion
        self._fit_box_initialization()
        
        # Select candidates and determine data extraction boxes.
        self._fit_table  = self._select_candidates()
        self._calculate_boxes()
        self._do_fitting()
        
        # Plot the fits
        if self._fwhm_plot is not None:
            self._plot_fits()
        else:
            self._logger.debug('Skipping plotting of star cutouts and best fits.')
        
        # Adjust the reported fitted positions back to image frame
        # from cutout fram.
        self._cutout_to_image_positions()
        return
        
    def _cutout_to_image_positions(self):
        """
        Convert the best fit positions from the cutout subimage frame
        to the full image pixel frame.
        """
        
        self._fit_table['xc_fit'] += self._fit_table['xmin']
        self._fit_table['yc_fit'] += self._fit_table['ymin']
        return
        
    def _calculate_boxes(self):
        """
        Calculate the pixel indices of the boxes from which data will
        be extracted for fitting each star.
        """
        
        half_width = self._box_width_pix / 2
        nrst_xpix  = np.rint( self._fit_table['xcenter'] ).astype(int)
        nrst_ypix  = np.rint( self._fit_table['ycenter'] ).astype(int)
        self._fit_table['xmin'] = nrst_xpix - half_width
        self._fit_table['xmax'] = nrst_xpix + half_width
        self._fit_table['ymin'] = nrst_ypix - half_width
        self._fit_table['ymax'] = nrst_ypix + half_width
        
        # Make the column printing reasonable
        # Note that astropy needs integers to be padded to conform with
        # the TDISP keyword format.
        col_fmt_dict = {'xmin': '%8d',
            'xmax': '%8d',
            'ymin': '%8d',
            'ymax': '%8d'}
        for col in col_fmt_dict:
            self._fit_table[col].info.format = col_fmt_dict[col]  
        
        if not self._quiet:
            print('Fit candidates with data extraction box indices:')
            print(self._fit_table)
            print('')
        return
        
    def _calculate_rchisq(self,
            data_array,         # cutout image pixel data
            x_grid,             # x-values used in model
            y_grid,             # y-values used in model
            stddev_arr,         # cutout image uncertainty values
            fitted_mod,         # model to evaluate
            num_fit_pars):      # number of free fit parameters.
        """
        """
        
        data_size  = data_array.size
        dof        = data_size - num_fit_pars
        data_fit   = fitted_mod(x_grid, y_grid)
        deviations = (data_array - data_fit) / stddev_arr
        chisq_arr  = np.square(deviations)
        chisq_sum  = np.sum(chisq_arr)
        rchisq     = chisq_sum / dof
        
        return rchisq
        
    def _do_fitting(self):
        """
        Performs 1-D and 2-G gaussian fits to the stars in the
        internal fit_table.
           
        A 2-dimensional Gaussian plus flat background model is fit to
        each cut-out image. The model is initialized with the global
        background level, and source centroid and peak value from the
        input photometry list. The Levenberg-Marquardt non-linear
        least squares fitting algorithm is used to fit the data.
        The data is assumed to have Gaussian errors that are the square
        root of the pixel counts.
           
        There are some limitations to this at present:
        
        - The reduced chi-squared values obtained are unrealistic. It is
          not clear whether the errors are underestimated and/or the
          fitting has trouble converging to the true solution.
        """
        
        # Setup
        num_stars         = len(self._fit_table)
        sigma_to_fwhm     = 2.35482
        fwhm_to_sigma     = 1.0 / sigma_to_fwhm
        max_iterations    = 500
        
        # Fit for x and y after initially just fitting for amplitude,
        # fwhm, and theta?
        is_pos_fitted_for = True
                
        # Extract data from main image, stores in _pixel_array.
        self._extract_cutouts()
        
        # Add extra columns to the output table for the fit results
        self._prepare_fit_columns()
        
        self._logger.info(f'Starting to fit Const2D+Gaussian2D model to {num_stars} star cutouts.')
        perf_time_start = time.perf_counter() # highest res timer, counts sleeps

        # X and Y pixel index grids needed by the fitter
        x_grid, y_grid = np.mgrid[0:self._box_width_pix, 0:self._box_width_pix]
                
        # We initialize the 2-D Gaussians based on the initial FWHM
        # but with a slight asymmetry sig_x/sig_y = 1.1, theta approx 30 degrees.
        def_axrat = 1.05
        sig_y     = self._init_fwhm * fwhm_to_sigma
        sig_x     = def_axrat * sig_y
        rotang    = 1.1 # radians, approx 60 degrees.
        xpos      = self._box_width_pix/2
        ypos      = self._box_width_pix/2
        ampl      = 1

        # The background model starts off fixed.
        bg_mod    = models.Const2D(amplitude=self._init_bglvl)
        bg_mod.amplitude.fixed = True
            
        # 2-D gaussian model for the star.
        star_mod  = models.Gaussian2D(amplitude=ampl,
            x_mean=xpos,
            y_mean=ypos,
            x_stddev=sig_x,
            y_stddev=sig_y,
            theta=rotang)
            
        # Constraint theta to -pi/2 to pi/2
        # NOTE Disabled as it results in the fitter failing very often
        ##half_pi               = 0.5 * math.pi
        ##star_mod.theta.bounds = (-half_pi, half_pi)
            
        # Start fits at the initial centroid positions
        star_mod.x_mean.fixed = True
        star_mod.y_mean.fixed = True
        comb_mod  = star_mod + bg_mod

        fitter = fitting.LevMarLSQFitter(calc_uncertainties=True)

        # Iterate over stars, first updating initial values, then fitting,
        # then storing the fit results.
        for idx in range(num_stars):
            # Get better estimate for initial parameters, in particular
            # the amplitude.
            # Note: you can set a Parameter object by value directly,
            # (e.g. "bob=3") but to access it you must access its value
            # attribute (e.g. "x = bob.value")
            starid = self._fit_table['id'][idx]
            xpos   = self._fit_table['xcenter'][idx] - self._fit_table['xmin'][idx]
            ypos   = self._fit_table['ycenter'][idx] - self._fit_table['ymin'][idx]
            ampl   = self._fit_table['peak_adu'][idx]
            bglvl  = self._fit_table['bgmed_per_pix'][idx]
            comb_mod[0].x_mean    = xpos
            comb_mod[0].y_mean    = ypos
            comb_mod[0].amplitude = ampl
            comb_mod[1].amplitude = bglvl
            current_num_free_pars = 4
            
            # Estimated standard devation of the data.
            var_arr = np.where(self._pixel_array[idx] > 0,
                    self._pixel_array[idx],
                    1)
            mean_variance = np.mean(var_arr[var_arr != 1])
            rms_stddev    = math.sqrt(mean_variance)
            std_arr = np.where(var_arr != 1,
                np.sqrt(var_arr),
                rms_stddev)
            # If used, the weight array should have values of 1/sigma
            if self._use_weights:
                weights_arr = 1.0 / std_arr
            else:
                weights_arr = None
            
            # Fit
            self._logger.debug(f'Fitting star #{idx} (id={starid}) at fixed xpos={xpos:.2f}, ypos={ypos:.2f}')
            
            fitted_mod, fit_status, fit_message, fit_ok, uncert_vals = self._do_single_fit(fitter,
                comb_mod,
                self._pixel_array[idx],
                weights_arr,
                x_grid,
                y_grid,
                max_iterations,
                idx)
             
            # If the ift is OK, and fit_for_bg is true, then fit for the BG level.
            if (fit_ok) and (self._fit_for_bg is True):
                fitted_mod[1].amplitude.fixed = False
                current_num_free_pars += 1
                fitted_mod, fit_status, fit_message, fit_ok, uncert_vals = self._do_single_fit(fitter,
                    fitted_mod,
                    self._pixel_array[idx],
                    weights_arr,
                    x_grid,
                    y_grid,
                    max_iterations,
                    idx)
             
            # If the fit did NOT fail we can relax the constraints on the
            # position of the gaussian.
            if fit_ok and is_pos_fitted_for:
                fitted_mod[0].x_mean.fixed = False
                fitted_mod[0].y_mean.fixed = False
                current_num_free_pars += 2
                fitted_mod, fit_status, fit_message, fit_ok, uncert_vals = self._do_single_fit(fitter,
                    fitted_mod,
                    self._pixel_array[idx],
                    weights_arr,
                    x_grid,
                    y_grid,
                    max_iterations,
                    idx)
            
            # Calculate reduced chi squared.
            rchisq = self._calculate_rchisq(self._pixel_array[idx],
                x_grid,
                y_grid,
                std_arr,
                fitted_mod,
                current_num_free_pars)
            
            # Extract fit parameters
            self._fit_table['ampl'][idx]   = fitted_mod[0].amplitude.value
            self._fit_table['xc_fit'][idx] = fitted_mod[0].x_mean.value
            self._fit_table['yc_fit'][idx] = fitted_mod[0].y_mean.value
            self._fit_table['fwhm_x'][idx] = sigma_to_fwhm * fitted_mod[0].x_stddev.value
            self._fit_table['fwhm_y'][idx] = sigma_to_fwhm * fitted_mod[0].y_stddev.value
            self._fit_table['theta'][idx]  = fitted_mod[0].theta.value
            self._fit_table['fit_ok'][idx] = fit_ok
            self._fit_table['rchisq'][idx] = rchisq
            
            # Get the errors.
            if fit_ok:
                # Following https://github.com/astropy/astropy/pull/10552
                self._fit_table['ampl_err'][idx]   = fitted_mod[0].amplitude.std
                self._fit_table['xc_err'][idx]     = fitted_mod[0].x_mean.std
                self._fit_table['yc_err'][idx]     = fitted_mod[0].y_mean.std
                self._fit_table['fwhm_x_err'][idx] = sigma_to_fwhm * fitted_mod[0].x_stddev.std
                self._fit_table['fwhm_y_err'][idx] = sigma_to_fwhm * fitted_mod[0].y_stddev.std
                self._fit_table['theta_err'][idx]  = fitted_mod[0].theta.std
            
            # Calculate axis ratio and circularity.
            # Note this calculation is not quite right because the two
            # parameters are NOT independent.
            fwhm_x    = self._fit_table['fwhm_x'][idx]
            fwhm_y    = self._fit_table['fwhm_y'][idx]
            axrat     = max(fwhm_x, fwhm_y) / min(fwhm_x, fwhm_y)
            axrat_err = 0
            
            if fit_ok:
                fwhm_xerr = self._fit_table['fwhm_x_err'][idx]
                fwhm_yerr = self._fit_table['fwhm_y_err'][idx]
                axrat_err = axrat * math.sqrt( (fwhm_xerr/fwhm_x)**2 +
                    (fwhm_yerr/fwhm_y)**2 )
                
                circular = self.is_circular(fwhm_x, fwhm_y, 
                    fwhm_xerr, fwhm_yerr)
                self._fit_table['circular'][idx] = circular

            self._fit_table['axrat'][idx]     = axrat
            self._fit_table['axrat_err'][idx] = axrat_err

        perf_time_end = time.perf_counter() # highest res timer, counts sleeps
        perf_time     = perf_time_end - perf_time_start  # seconds
        avg_time_ms   = 1000.0 * perf_time / num_stars   # milliseconds
        self._logger.debug(f'Fitting completed in {perf_time:.3f} s (average {avg_time_ms:.1f} ms/star)')
        if not self._quiet:
            print('Fit results:')
            print(self._fit_table)
            print('')
        return
        
    @classmethod
    def is_circular(cls, fwhm_x, fwhm_y, fwhm_xerr, fwhm_yerr):
        """
        Returns True if fwhm_y is within _circ_thresh_sigma
        standard deviations of fwhm_x, otherwise False
        """
        circular = True
        # Circularity is whether fwhm_y is within _circ_thresh_sigma
        # standard devaitions of fwhm_x
        d_fwhm = math.fabs(fwhm_y - fwhm_x)
        sigma  = d_fwhm / fwhm_yerr
        if sigma > ApMeasureStars._circ_thresh_sigma:
            circular = False
        return circular

    def _do_single_fit(self, 
        fitter, 
        input_mod,
        data_arr,
        weights_arr,
        x_grid,
        y_grid,
        max_iterations,
        index):
        """
        Utility to perform a single fit of a model to the data,
        checking for error conditions along the way.
        """
        
        fitted_mod = fitter(input_mod, 
            x=x_grid, 
            y=y_grid, 
            z=data_arr,
            weights=weights_arr,
            maxiter=max_iterations)
        
        # Check fit infomation to see if the fit completed successfully.
        fit_status  = fitter.fit_info['ierr']
        fit_message = fitter.fit_info['message']
        fit_ok      = False
        cov_arr     = fitter.fit_info['param_cov']
        is_array    = isinstance(cov_arr, np.ndarray)
                
        # diagnostic message
        ##print(f'Diagnostic: Fit status {fit_status} for star {index}')
        ##print(f'Diagnostic: Fitter info message was: [{fit_message}]')
        if (fit_status > 0) and (fit_status <= 4) and is_array:
            # These may be bogus without correctly weighting the data.
            # These are in the order of the model parameters, except that
            # fixed parameters are excluded. See https://github.com/astropy/astropy/issues/7202
            
            err_param = np.sqrt(np.diag(cov_arr))
            fit_ok    = True
        else:
            self._logger.warning(f'Non-nominal fit status {fit_status} for star {index}')
            self._logger.warning(f'Fitter info message was: [{fit_message}]')
            fit_ok    = False
            err_param = None
        return fitted_mod, fit_status, fit_message, fit_ok, err_param

    def _extract_cutouts(self):
        """
        Create 3-D data array to store N stars x MxM pixels
        """
        num_stars = len(self._fit_table)
        wdth      = self._box_width_pix
        numel     = num_stars * wdth * wdth
        tmp_arr = np.zeros(numel, 
            dtype=self._img_data.dtype).reshape(num_stars, wdth, wdth)
            
        for idx in range(num_stars):
            xmin         = int( self._fit_table['xmin'][idx] )
            xmax         = int( self._fit_table['xmax'][idx] )
            ymin         = int( self._fit_table['ymin'][idx] )
            ymax         = int( self._fit_table['ymax'][idx] )
            cutout       = self._img_data[ymin:ymax, xmin:xmax]
            tmp_arr[idx] = cutout.copy()
            
        maxval = np.amax(tmp_arr)
        minval = np.amin(tmp_arr)
        self._logger.debug(f'In {num_stars}x{wdth}x{wdth} pixel array min={minval:.2f}, max={maxval:.2f}')
        self._pixel_array = tmp_arr
        return

        
    def _fit_box_initialization(self):
        """
        Sets up variables related to the size of the region used in
        fitting and the edge exclusion used in candidate selection.
        """
        
        # We want the fit box to be at least 2x the initial estimated
        # FWHM, and also an even number of pixels. We also pad a bit 
        # in case the initial FWHM estimate is an underestimate.
        pad_frac            = 3.0
        min_width           = 12
        self._box_width_pix = 2 * int(pad_frac * self._init_fwhm)
        if self._box_width_pix < min_width:
            self._box_width_pix = min_width
        
        # Edge exclusion shouldbe >= half the box size, and also an even
        # number of pixels.
        self._edge_excl_pix = 2 * int(math.ceil(self._box_width_pix/4)) 
        self._logger.debug(f'Adopting a fit box width of {self._box_width_pix} pixels' + 
            f', edge exclusion of {self._edge_excl_pix} pixels.')
        return
        
    def _get_fit_artists(self, index):
        """
        Returns a list of matplotlib artists to be added to a plot
        based on the Gaussian fits
        """
        artist_list = []
        
        xc     = self._fit_table['xc_fit'][index]
        yc     = self._fit_table['yc_fit'][index]
        fwhm_x = self._fit_table['fwhm_x'][index]
        fwhm_y = self._fit_table['fwhm_y'][index]
        theta  = self._fit_table['theta'][index]
        ellipse = Ellipse(xy=(xc, yc),
                width=fwhm_x, 
                height=fwhm_y,
                angle=theta,
                edgecolor='r',
                facecolor=None,
                fill=False,
                linewidth=0.8)
        artist_list.append( ellipse )
        
        return artist_list
        
    def _get_subplot_fitinfo(self, index):
        """
        Extract a short summary of the fit results for the star at
        the given index in the fit_table, suitable for use in the
        subplots.
        """
        
        # Get info from fit
        fwhm_x      = self._fit_table['fwhm_x'][index]
        fwhm_y      = self._fit_table['fwhm_y'][index]
        fit_ok      = self._fit_table['fit_ok'][index]
        symmetric   = self._fit_table['circular'][index]
        if not fit_ok:
            symmetric = 'N/A'
        info_list   = [f'FWHM_X={fwhm_x:.1f}',
            f'FWHM_Y={fwhm_y:.1f}',
            f'FitOK={fit_ok}',
            f'Circular={symmetric}']
        fitinfo_str = '\n'.join(info_list)
        return fitinfo_str
        
    def _get_subplot_title(self, index):
        """
        Return a title string for a subplot showing a single fitted star.
        """
        xcen  = self._fit_table['xcenter'][index]
        ycen  = self._fit_table['ycenter'][index]
        idnum = self._fit_table['id'][index]
        
        xcen = int( round( xcen ) )
        ycen = int( round( ycen ) )
        title_str = f'Star {idnum:d} at x={xcen:d}, y={ycen:d}'
        return title_str
        
    def _initialize_logger(self, loglevel):
        """
        Initialize and return the logger
        """
        
        logger = logging.getLogger('ApMeasureStars')
        
        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: {}'.format(loglevel))
        logger.setLevel(numeric_level)
    
        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(numeric_level)
    
        # create formatter
        formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s | %(message)s')
    
        # add formatter to ch
        ch.setFormatter(formatter)
    
        # add ch to logger
        logger.addHandler(ch)
        return logger
        
    def _plot_fits(self):
        """
        Plot all the fits.
        """

        # Make axis box stand out as viridis can be dark.
        matplotlib.rc('axes', edgecolor='r')

        # Hand-tuning suggests a linear scale between the 
        # min max of these stars is the best. Global limits or
        # percentile limits tend to saturate out the stars.
        ##pct_interval = AsymmetricPercentileInterval(0.10, 99.9)
        valmin = self._init_bglvl
        valmax = np.amax(self._fit_table['peak_adu'])
        minmax_interval = ManualInterval(vmin=valmin,
            vmax=valmax) 
        
        norm_func = ImageNormalize(self._img_data, 
            interval=minmax_interval, 
            stretch=AsinhStretch())
        ##norm_func    = ImageNormalize(self._img_data, 
        ##    interval=pct_interval, 
        ##    stretch=SqrtStretch())

        # May be overly agressive for subplots
        ##norm_func   = ImageNormalize(self._img_data,
        ##    interval=pct_interval,
        ##    stretch=AsinhStretch())
            
        # Location of text label added internal to subplots, in
        # normalized 0:1 coordinate. This is where the bottom left
        # the text box will appear.
        text_x = 0.2
        text_y = 0.8
            
        region_list     = ['TL', 'TR', 'CN', 'BR', 'BL']
        region_row_dict = {'TL': 0,
            'TR': 1,
            'CN': 2,
            'BR': 3,
            'BL': 4}

        # Create a dictionary that has keys of the row index within
        # _fit_table (also the first index of the _pixel_array)
        # and the column index within the 2-D plotting array.
        #
        # This is necessary because there may be less that _num_per_reg
        # stars per region.
        #
        # This is inelegant but should work.
        num_stars     = len(self._fit_table)
        star_col_dict = {}         # To fill
        col_idx_dict  = {'CN': 0,  # working memory
            'TL': 0, 
            'TR': 0, 
            'BL': 0, 
            'BR': 0}
        for idx in range(num_stars):
            this_reg               = self._fit_table['region'][idx]
            this_col_idx           = col_idx_dict[this_reg]
            star_col_dict[idx]     = this_col_idx
            new_col_idx            = this_col_idx + 1
            col_idx_dict[this_reg] = new_col_idx
            
        # Create figure for plot
        fig_rows = len(region_list)
        fig_cols = self._num_per_reg

        # List of handles to the images
        plt_imgs = []
        
        fig, ax_arr = plt.subplots(nrows=fig_rows, 
            ncols=fig_cols,
            sharex=True, 
            sharey=True)
        
        (median_fwhm, madstd_fwhm, npts) = self.median_fwhm('both')
        
        title = 'Star PSF measurements using 2-D Gaussian fits'
        if self._plot_title is not None:
            title = self._plot_title
        title += f'\nMedian FWHM={median_fwhm:.2f} +/- {madstd_fwhm:.2f} (MAD stddev) pixels'
        fig.suptitle(title, fontsize=7)
        
        for idx in range(num_stars):
            this_reg = self._fit_table['region'][idx]
            row_idx  = region_row_dict[this_reg]
            col_idx  = star_col_dict[idx]

            # Plot the cutout image.
            cut_out = self._pixel_array[idx]
            plt_imgs.append( ax_arr[row_idx, col_idx].imshow(cut_out, 
                origin='lower', 
                norm=norm_func) )
                 
            # Titles, fit info etc
            title_str   = self._get_subplot_title(idx)
            fitinfo_str = self._get_subplot_fitinfo(idx) 
            ax_arr[row_idx, col_idx].set_title(title_str,
                fontsize=4,
                pad=3)
            ax_arr[row_idx, col_idx].text(text_x, 
                text_y, 
                fitinfo_str, fontsize=3, color='w')
            #ax_arr[row_idx, col_idx].set_xlabel('X-axis (pixels)', fontsize=6)
            ax_arr[row_idx, col_idx].set_ylabel(this_reg,
                fontsize=6)
                
            # Get artists related to best-fit parameters
            artist_list = self._get_fit_artists(idx)
            if len(artist_list) > 0:
                for art in artist_list:
                    ax_arr[row_idx, col_idx].add_artist(art)
        
        # Only show tick values on outer edge of grid.
        for ax in ax_arr.flat:
            ax.tick_params(axis='both', 
                labelsize=4, 
                direction='in', 
                color='r',
                length=3)
            ax.label_outer()
        
        plt.savefig(self._fwhm_plot,
            dpi=200,
            bbox_inches='tight')
        self._logger.info(f'Plotted star FWHM fit cut-outs to {self._fwhm_plot}')
        return

    def _prepare_fit_columns(self):
        """Utility function to add extra columns to the fit_table
           that will be used for the best-fit parameters.
        """
        
        # This adds a large number of columns, but we'll be storing the
        # data in a file so it doesn't matter.
        
        new_col_dict = {'xc_fit': 0.0,
            'xc_err':     0.0,
            'yc_fit':     0.0,
            'yc_err':     0.0,
            'ampl':       0.0,
            'ampl_err':   0.0,
            'fwhm_x':     0.0,
            'fwhm_x_err': 0.0,
            'fwhm_y':     0.0,
            'fwhm_y_err': 0.0,
            'theta':      0.0,
            'theta_err':  0.0,
            'axrat':      0.0,
            'axrat_err':  0.0,
            'circular':   True,
            'fit_ok':     True,
            'rchisq':     0.0}
        for newcol in new_col_dict:
            self._fit_table[newcol] = new_col_dict[newcol]
            if 'circular' in newcol:
                pass
            elif 'fit_ok' in newcol:
                pass
            elif ('fwhm' in newcol) or ('theta' in newcol):
                self._fit_table[newcol].info.format = '%.3f' 
            else:
                self._fit_table[newcol].info.format = '%.2f' 
        return
        
    def _select_candidates(self):
        """
        Select candidate stars in the center and four quadrants
        
        Fitting is computationally expensive and will not give good
        results on saturated stars, blended/crowded stars, and faint
        stars. Thus it makes sense to limit fitting to a subset of the
        available stars in the input source list. In order to measure
        whether the PSF changes across the image it is also useful to
        break the image up into a series of regions. In the center
        the PSF is expected to be the sharpest and most symmetric,
        but at the edges of the image asymmetries may be expected. Hence
        we select a set of 5 candidate stars in each of the regions
        shown schematically below.
        
        Note that in the case of a square image the area of the central
        region CN is very close to the remaining area of each of the outer 
        quadrants TL, TR, BL, BR (despite what the ASCII depiction below
        seems to show)
        
        For an image of height H, width W, H <= W, central region 
        diameter C=H then the area of the central circle is
        A_C = pi * H^2 / 16,
        and
        A_Q = H/4 * (W - (pi*H/16))::
        
            +-----------+-----------+
            [           |           ]
            [ TL        |        TR ]
            [          -+-          ]
            [         / | \         ]
            [        |CN|  |        ]
            +--------+--+--+--------+
            [        |  |  |        ]
            [ BL      \ | /      BR ]
            [          -+-          ]
            [           |           ]
            [           |           ]
            +-----------+-----------+
        
        Candidate stars should not be too close to the edge of the 
        detector. Stars with centroids within _edge_excl_pix of the 
        edge are not selected as candidates.
        
        Source Confusion:
        It is also important that candidate stars should not have another
        star withing the image cutout used for fitting. Source confusion
        is problematic in that the input estimated magnitudes will be
        incorrect and fitting will be compromized. A kdtree is used to
        find the nearest neighbor of each source in the trimmed _init_srcs
        list, based on the full set of sources. The trimmed list is further
        filtered to remove all stars having neighbors within a radius of
        _box_width.
        """
        
        # Algorithm:
        #
        # Identify nearest neighbors and trim the input source list of
        # star that have neighbors within a radius of _box_width.
        #
        # Iterate over each source i in the trimmed source list
        #   If not saturated, calculate dx = x_i - x_c, dy = y_i - y_c, r_i
        #   If r_i < H/2, add source id, mag to center list
        #   else if dx < 0, dy > 0 add to TL list, etc
        # The sort lists by magnitude
        # Pick the top N from each quadrant and place in a table
        # The table forms the basis for self._fit_table
        
        candidate_table = None
        
        self._trim_neighbors()
        
        num_srcs        = len(self._init_srcs)
        self._logger.info(f'Selecting candidate stars for fitting out of {num_srcs} stars in {self._rows} row x {self._cols} column image.')

        radius = float( min(self._cols, self._rows) ) / 4
        xcen   = float(self._cols) / 2
        ycen   = float(self._rows) / 2
        self._logger.debug(f'Center has radius {radius:.1f} pixels centered on x={xcen:.1f}, y={ycen:.1f}')
        
        self._init_srcs['dx'] = self._init_srcs['xcenter'] - xcen
        self._init_srcs['dy'] = self._init_srcs['ycenter'] - ycen
        xsq = np.square(self._init_srcs['dx'])
        ysq = np.square(self._init_srcs['dy'])
        self._init_srcs['radius'] = np.sqrt(xsq + ysq)
        self._init_srcs['region'] = 'XX'
            
        # Select region
        is_cn     = self._init_srcs['radius'] <= radius
        is_right  = self._init_srcs['dx'] >= 0
        is_left   = np.logical_not( is_right )
        is_top    = self._init_srcs['dy'] >= 0
        is_bottom = np.logical_not( is_top )
        
        is_tr     = np.logical_and(is_top,    is_right)
        is_tl     = np.logical_and(is_top,    is_left)
        is_br     = np.logical_and(is_bottom, is_right)
        is_bl     = np.logical_and(is_bottom, is_left)
        
        # This is technically superflous considering we have the flag
        # arrays, but its useful for human-readable analysis and we'll
        # use it later in the candidate list.
        self._init_srcs['region'][is_tr] = 'TR'
        self._init_srcs['region'][is_tl] = 'TL'
        self._init_srcs['region'][is_br] = 'BR'
        self._init_srcs['region'][is_bl] = 'BL'
        self._init_srcs['region'][is_cn] = 'CN'
        
        # Edge exclusion masks, taking into account python 0:n-1 indexing.
        xmin = self._edge_excl_pix - 1
        xmax = self._cols - self._edge_excl_pix - 1
        ymin = self._edge_excl_pix
        ymax = self._rows - self._edge_excl_pix - 1
        ok_left = self._init_srcs['xcenter'] >= xmin
        ok_rght = self._init_srcs['xcenter'] <= xmax
        ok_bot  = self._init_srcs['ycenter'] >= ymin
        ok_top  = self._init_srcs['ycenter'] <= ymax
        ok_mask = np.logical_and(ok_left, ok_rght)
        ok_mask = np.logical_and(ok_mask, ok_bot)
        ok_mask = np.logical_and(ok_mask, ok_top)
        srclist = self._init_srcs[ok_mask]
        num_ok  = len(srclist)
        self._logger.debug(f'There are {num_ok} stars more than {self._edge_excl_pix} pixels away from the edge of the detector.')
                            
        for reg in ['CN', 'TL', 'TR', 'BR', 'BL']:
            row_select = srclist['region'] == reg
            cand_table = srclist[row_select]
            cand_table.sort(['magnitude'], reverse=False) # want most negative first...
            num_cand   = len(cand_table)
            
            if num_cand == 0:
                self._logger.warning(f'There are no candidates in the {reg} region.')
                continue
            else:
                self._logger.debug(f'In region {reg} found {num_cand} candidates for fitting:')

            # Select the Nth to Mth brightest candidates in each region,
            # unless there are not enough stars to do so.
            n_start    = self._skip_brightest
            m_end      = n_start + self._num_per_reg
            while m_end > num_cand:
                n_start = max(0, n_start-1)
                m_end   = max(0, m_end-1)
                if m_end < num_cand:
                    self._logger.error('Logic error in _select_candidates?')
                    self._logger.error(f'Region={reg}, num_cand={num_cand}, n_start={n_start},' + 
                        f' m_end={m_end}, num_per_reg={self._num_per_reg}, skip_brightest={self._skip_brightest}')
                    continue
                
            cand_table = cand_table[n_start:m_end]
            if not self._quiet:
                print(cand_table)
                print('')
        
            # Now build the candidate table
            if candidate_table is None:
                candidate_table = cand_table
            else:
                candidate_table = vstack([candidate_table, cand_table])
            
        return candidate_table
        
    def _trim_neighbors(self):
        """
        Remove stars from _init_srcs that have a neighbor from 
        _full_srcs within a radius of _box_width pixels.
           
        This uses a kdtree to identify the nearest neighbors.
        
        Note: This algorithm can fail to identify some close neighbors
        if the input "full" sourclist has excluded saturated stars. 
        """
        
        rad       = self._box_width_pix
        init_size = len(self._init_srcs)
        self._logger.info(f'Preparing to trim the input source list of stars within neighbors within {rad} pixels.')
        
        x      = self._full_srcs['xcenter']
        y      = self._full_srcs['ycenter']
        xy_pts = np.column_stack((x,y))
        tree   = spatial.KDTree(xy_pts)
        self._logger.debug(f'Constructed a kdtree with {len(x)} data points.')
        
        # Add  a column to _init_srcs for the nn_dist
        self._init_srcs['nn_dist']  = 0.0
        self._init_srcs['nn_dist'].info.format = '%.2f' 
        
        for idx in range(init_size):
            # The nearest point to a point in the tree is itself, so we
            # are interested in k=2 (2 nearest points).
            # A euclidean distance metric is p=2
            xypos = [self._init_srcs['xcenter'][idx],
                self._init_srcs['ycenter'][idx] ]
            
            # d is a list with k elements, the distance to each of the 
            # k'th nearest neighbors.
            # i is a list with k elements, the indices of the neighbors
            # within the original data points.
            d, i = tree.query(xypos, k=2, p=2)
            
            nn_dist = d[1]
            self._init_srcs['nn_dist'][idx] = nn_dist

        print('Initial sources with nearest neighbor distance\n', 
            self._init_srcs)
            
        # Create trutch mask for nn_dist greater than exclusion radius
        mask            = self._init_srcs['nn_dist'] >= rad
        self._init_srcs = self._init_srcs[mask]
            
        final_size  = len(self._init_srcs)
        num_removed = init_size - final_size
        self._logger.info(f'Nearest neighbor filtering removed {num_removed} stars from consideration.')
        return
        
        
    def median_fwhm(self, direction):
        """
        Returns the sigma-clipped median fitted FWHM over both X 
        and Y in pixels over all stars that fitted successfully, 
        along with median absolute deviation (MAD) standard deviation.
           
        Note that the deviation return is standard deviation based on
        the MAD, not the MAD itself. See astropy.stats.mad_std
        
        Direction must be one of 'both', 'x', or  'y'
        """
        
        # Clip values that are more than this number of sigma from the
        # median.
        numsig = 3.0
        
        # Select only cases where the fitting appeared to work correctly.
        ok_x    = self._fit_table['fwhm_x'][self._fit_table['fit_ok']]
        ok_y    = self._fit_table['fwhm_y'][self._fit_table['fit_ok']]
        
        if 'both' in direction:        
            ok_fwhm = np.concatenate((ok_x, ok_y)) # Note inputs as tuple
        elif 'x' in direction:
            ok_fwhm = ok_x 
        elif 'y' in direction:
            ok_fwhm = ok_y
        
        clipped  = sigma_clip(ok_fwhm, sigma=numsig, masked=False)
        num_used = len(clipped)
        self._logger.debug(f'Estimating median FWHM (direction={direction}) using {num_used} FWHM measurements ({len(ok_fwhm)} OK fits before clipping).')
        
        median_fwhm = float( np.median(clipped) )
        madstd_fwhm = float( mad_std(clipped) )
        return (median_fwhm, madstd_fwhm, num_used)
        
    def results_table(self):
        """
        Return the fitting results as an astropy Table
        """
        return self._fit_table
