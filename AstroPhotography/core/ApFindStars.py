# -*- coding: utf-8 -*-
#
#  Contains the implementation of ApFindStars
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

# 2020-08-05 dks : Initial coding
# 2020-08-08 dks : Final script for writing XY source list and approx mags.
# 2020-08-09 dks : Added initial and temporary astrometry.net query.
# 2020-08-10 dks : Moved search params to CLI, flag possible saturation,
#                  use percentile interval in plotting/visualization.
# 2020-08-20 dks : Refactor into ApFindStars
# 2020-08-22 dks : Implemented plotfile bitmapped image with sources.
# 2020-08-23 dks : Add stellar FWHM fitting and plotting class.
# 2020-09-01 dks : Switch to yaml from json as json float formatting is bad.
# 2020-09-02 dks : Select candidates with no nearby stars using a kdtree.
# 2021-01-18 dks : Move ApFindStars into core, split ApMeasureStars into its
#                  own file.

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

from photutils import make_source_mask, find_peaks, DAOStarFinder
from photutils import CircularAperture, CircularAnnulus, aperture_photometry

# AstroPhotography includes    
from .ApMeasureStars import ApMeasureStars as ApMeasureStars
    
def yaml_float_representer(dumper, value):
    """Change default yaml float representation.
    
    From: https://stackoverflow.com/questions/33944299/how-to-round-numeric-output-from-yaml-dump-in-python
    """
    text = '{0:.6f}'.format(value)
    return dumper.represent_scalar(u'tag:yaml.org,2002:float', text)

class ApFindStars:
    """Find and characterize stars within a FITS image
    """

    # Class constants
    GOOD        = 0
    INPUT_ERROR = 1

    def __init__(self, fitsimg, extnum, search_fwhm,
        search_nsigma, detector_bitdepth, 
        max_sources, nosatmask, sat_frac, loglevel,
        plotfile, quiet):
        """Constructor for ApFindStars
        
        This reads the input image, and performs by default the following
        processing steps:
        - 
        
        The user may update the default source detection and photometry 
        by calling the public member functions of this class. Output,
        whether FITS table source list or optional quality report
        and ds9-format region file, must be initiated by the user using
        the public write_ member functions.
        """
        
        self._status        = ApFindStars.GOOD
        self._loglevel      = loglevel
        self._fitsimg       = fitsimg
        self._extnum        = extnum
        self._search_fwhm   = search_fwhm
        self._search_nsigma = search_nsigma        
        self._bitdepth      = detector_bitdepth
        self._max_sources   = max_sources
        self._nosatmask     = nosatmask
        self._sat_frac      = sat_frac
        self._plotfile      = plotfile
        self._max_adu       = math.pow(2, detector_bitdepth) - 1
        self._sat_thresh    = math.floor(sat_frac * self._max_adu)
        self._quiet         = quiet
        
        # Not calculated by default, only if a user calls measure_fwhm
        self._psf_table     = None
        self._fwhm_both     = None
        self._fwhm_x        = None
        self._fwhm_y        = None

        # Number of sources detected, and number of sources that had
        # photometry, and number of sources that had FWHM/PSFs fitted.
        self._nsrcs_detected  = 0
        self._nsrcs_photom    = 0
        self._nsrcs_fitted    = 0
        self._nsrcs_saturated = 0 

        # Hard-wired constants
        self._ap_fwhm_mult    = 2.0  # Aperture radius is this times search_fhwm

        # Set up logging
        self._logger = self._initialize_logger(self._loglevel)

        # Read the image
        self._data, self._hdr = self._read_fits(fitsimg, extnum)
    
        # Get initial background estimates using hard-wired constants.
        self._bg_mean, self._bg_median, self._bg_stddev = sigma_clipped_stats(self._data, 
            sigma=3.0)
        self._logger.debug('Sigma clipped image stats: mean={:.3f}, median={:.3f}, stddev={:.3f}'.format(self._bg_mean, self._bg_median, self._bg_stddev))  
    
        tmp_mask = make_source_mask(self._data, nsigma=2, npixels=5, dilate_size=11)
        self._bg_mean, self._bg_median, self._bg_stddev = sigma_clipped_stats(self._data, 
            sigma=3.0, 
            mask=tmp_mask)
        self._logger.debug('Source-masked image stats: mean={:.3f}, median={:.3f}, stddev={:.3f}'.format(self._bg_mean, self._bg_median, self._bg_stddev))  
        
        # Generate mask for DAOStarFinder. 0 is good, 1 is bad.
        self._mask = np.zeros(self._data.shape, dtype=bool)
        self._logger.debug('Checking for possibly saturated stars or regions.')
        saturated_locations = self._find_saturated(self._data, 
            self._sat_thresh, 
            self._search_fwhm)
        num_sat_candidates = 0
        if saturated_locations is not None:
            num_sat_candidates = len(saturated_locations)
        if not self._nosatmask:
            # Find possibly saturated stars and exclude them from source
            # detection and photometry.
            box_width = int(4 * self._search_fwhm)
            ncols     = self._data.shape[1]
            nrows     = self._data.shape[0]
    
            self._logger.debug(f'Excluding {num_sat_candidates} possibly saturated stars using mask boxes of half-width {box_width} pixels.')
            # Masking format is mask[rowmin:rowmax+1, colmin:colmax] = True
            # But as we're counting the center as 1 pixel the first included 
            # row is = rowcen - width + 1. The last included row is rowcen
            # + width - 1. Adding one to the last one because python has 
            # exclusive ranges gives us the following...
            for isat in range(num_sat_candidates):
                scol = saturated_locations['x_peak'][isat]
                colmin = max(0, scol - box_width + 1)
                colmax = min(ncols, scol + box_width)
                srow = saturated_locations['y_peak'][isat]
                rowmin = max(0, srow - box_width + 1)
                rowmax = min(nrows, srow + box_width)
                self._mask[rowmin:rowmax, colmin:colmax] = True
        else:
            # Leave possibly saturated stars in the image.
            self._logger.debug(f'Retaining {num_sat_candidates} possibly saturated stars in source searching and photometry.')
        self._nsrcs_saturated = num_sat_candidates
        
        
        # Search for stars using the supplied FWHM and threshold.
        self.source_search(self._search_fwhm, self._search_nsigma)
        
        # Use initial source list for aperture photometry, generate
        # self._phot_table
        self.aperture_photometry()
    
        if self._plotfile is not None:
            self.plot_image(self._plotfile)
        return
        
    def trim(self, max_srcs):
        """Reduce the size of the photometry and XY positions tables
           to at most max_srcs of the brightest sources.
           
        This function is useful in restricting the number of sources
        written to the output FITS source list. Supplying too many
        stars to astrometry.net (via ap_astrometry) can slow or break
        it.
        
        This function is called by aperture_photometry using the
        constructed value of max_sources unless it aperture_photometry
        was called with notrim=False.
        """
        
        # Trim and update apertures
        self._phot_table = self._trim_table(self._phot_table, 
            max_srcs)
        self._apertures, _, _  = self._make_apertures(self._phot_table['xcenter'], 
            self._phot_table['ycenter'])
        return
        
    def plot_image(self, plotfile):
        """Plot an asinh-stretched image with the current apertures to
           standard graphics bitmap-format file.
        """
        
        # Normally the range 0.5 percent to 99.5 percent clips off the
        # outliers.
        pct_interval = AsymmetricPercentileInterval(0.50, 99.5)

        #sqrt_norm    = ImageNormalize(self._data, 
        #    interval=pct_interval, 
        #    stretch=SqrtStretch())
        asinh_norm   = ImageNormalize(self._data,
            interval=pct_interval,
            stretch=AsinhStretch())
    
        fig, ax = plt.subplots()
        ax.tick_params(axis='both', labelsize=8)

        im = ax.imshow(self._data, 
            origin='lower', 
            norm=asinh_norm)
        self._apertures.plot(color='red', lw=1.5, alpha=0.5)
         
        # Clean up file name string to prevent _ becoming subscripts.
        #fname_str = self._fitsimg.replace('_', '\_') # Not necessary?
        fname_str = self._fitsimg
        
        # Subtitle
        num_stars = len(self._phot_table)
        info_str  = f'{num_stars} stars'
        if self._nosatmask is not None:
            info_str += ' excluding saturated stars'
        else:
            info_str += ' including saturated stars'
        if self._max_adu is not None:
            info_str = 'Brightest ' + info_str
            
        ax.set_title(f'{fname_str}\n{info_str}', fontsize=8)
        ax.set_xlabel('X-axis (pixels)', fontsize=8)
        ax.set_ylabel('Y-axis (pixels)', fontsize=8)
        #plt.show()
        plt.savefig(self._plotfile,
            dpi=200, quality=95, optimize=True,
            bbox_inches='tight')
        self._logger.info(f'Plotting asinh-stretched bitmap of image and sources to {self._plotfile}')
        return
        
    def _make_apertures(self, colpos, rowpos):
        """Create an apertures, bg annular, and annular mask object for 
           aperture photometry.
           
        See https://photutils.readthedocs.io/en/stable/aperture.html#sigma-clipped-median-within-a-circular-annulus
        """

        positions = np.transpose((colpos, rowpos))
        
        # Radius of circular source region radius
        ap_radius = math.ceil(self._ap_fwhm_mult * self._search_fwhm)
        
        # BG outer annulus radius. By making this 1.5 times ap_radius
        # the BG annulus has an area 1.25 times the inner circle, so
        # we should get decent sampling.
        outer_radius = math.ceil(1.5 * ap_radius) 
        
        self._logger.debug(f'Radius of circular aperture photometry is {ap_radius} pixels.')
        self._logger.debug(f'Local BG estimation in annulus outer radius {outer_radius} pixels, inner radius {ap_radius} pixels.')
        
        apertures         = CircularAperture(positions, r=ap_radius)
        bg_apertures      = CircularAnnulus(positions, r_in=ap_radius, 
            r_out=outer_radius)
        bg_aperture_masks = bg_apertures.to_mask(method='center')
        
        return apertures, bg_apertures, bg_aperture_masks

    def source_search(self, search_fwhm, search_nsigma):
        """Search for star-like objects with the given FWHM that have
           a statistical significance nsigma above the established 
           background noise level.
        """
        
        self._logger.debug(f'Running DAOStarFinder with FWHM={search_fwhm:.2f} pixels, threshold={search_nsigma} BG sigma.')
        daofind = DAOStarFinder(fwhm=search_fwhm, 
            threshold=search_nsigma*self._bg_stddev)  
        sources = daofind(self._data - self._bg_median, 
            mask=self._mask)  
            
        # Clean the formatting up a bit using a column format dictionary
        col_fmt_dict = {'xcentroid': '%.2f',
            'ycentroid': '%.2f',
            'id':   '%d',
            'npix': '%d',
            'sky':  '%d',
            'peak': '%.2f'}
        for col in sources.colnames:
            if col in col_fmt_dict:
                new_fmt = col_fmt_dict[col]
            else:
                new_fmt = '%.4f'
            sources[col].info.format = new_fmt  # for consistent table output
        num_srcs = len(sources)
        self._logger.info(f'Initial source list using FHWM={search_fwhm} pixels, threshold={search_nsigma} x BG stddev, found {num_srcs} sources.')
        
        # Check for saturated sources, as nosatmask may have been false
        # when initially search
        self._logger.debug(f'Sources with pixels exceding {self._sat_thresh} ADU will be flagged as possibly saturated.')
        sources['psbl_sat'] = sources['peak'] > self._sat_thresh
        n_sat = np.sum(sources['psbl_sat'])
        self._logger.debug(f'There are {n_sat} possibly saturated stars in the output source list.')

        if not self._quiet:
            print(sources)
            
        # Set official source list and number of sources detected
        self._sources = sources
        self._nsrcs_detected = len(sources)
        return

    def write_source_list(self, output_fits_table):
        """Write the current source position and photometry data to
        a FITS table, overwriting any existing file
        """
        
        # Build or rebuild keyword list useful for the FITS table as well as the
        # image quality report.
        self._kw_dict = self._build_keyword_dictionary(self._fitsimg,
            self._hdr,
            self._bg_median,
            self._bg_stddev)
        
        self._write_source_list(output_fits_table, 
            self._fitsimg,
            self._max_sources, 
            self._hdr, 
            self._kw_dict,
            self._phot_table,
            self._psf_table)
        return

    def aperture_photometry(self, notrim=None):
        """Perform aperature photometry using the existing sources
           and apertures.
           
        If called with notrim=True then the number of sources
        output will NOT be trimmed to the max_sources value
        supplied in the constructor. This allows the user to override
        trimming if they update the initial source searching and 
        photometry.
        """
        
        # Whether to trim to the constructor self._max_sources value.
        dont_trim = False
        if notrim is not None:
            dont_trim = notrim
        
        self._apertures, bg_apertures, bg_aperture_masks = self._make_apertures(self._sources['xcentroid'], 
            self._sources['ycentroid'])
        
        # Find the median BG level in the background annuli
        bkg_median = []
        for mask in bg_aperture_masks:
            annulus_data         = mask.multiply(self._data)
            annulus_data_1d      = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)
        
        # Perform photometry but get rid of quantities, as they complicate
        # life.
        phot_table = Table( aperture_photometry(self._data, 
            self._apertures) )
            
        # Expected background counts in source aperture
        bg_in_src_ap = bkg_median * self._apertures.area
        
        # Correct the current aperture sum for the expected BG
        phot_table['aperture_sum'] = phot_table['aperture_sum'] - bg_in_src_ap            
        phot_table['aperture_sum'].info.format = '%.4f'  # for consistent table output
        
        # Copy over peak ADU and possibly saturated flag from initial
        # sourcelist 
        phot_table['peak_adu'] = self._sources['peak']
        phot_table['psbl_sat'] = self._sources['psbl_sat']
        phot_table['bgmed_per_pix'] = bkg_median
        
        exposure = None
        for kw in ['EXPOSURE', 'EXPTIME']:
            if exposure is None:
                if kw in self._hdr:
                    exposure = float( self._hdr[kw] )
                    self._logger.debug(f'Image exposure time [seconds]: {exposure:.2f}')
        if exposure is None:
            self._logger.warning('EXPOSURE not found in FITS header. Assuming 1 second.')
            exposure = 1
            
        phot_table['adu_per_sec'] = phot_table['aperture_sum'] / exposure
        phot_table['magnitude']   = -2.5 * np.log10(phot_table['adu_per_sec'])
        
        phot_table['adu_per_sec'].info.format = '%.4f'
        phot_table['magnitude'].info.format   = '%.4f'
        phot_table['xcenter'].info.format     = '%.2f'
        phot_table['ycenter'].info.format     = '%.2f'
        phot_table['bgmed_per_pix'].info.format    = '%.2f'
    
        phot_table.sort(['adu_per_sec', 'xcenter', 'ycenter'], reverse=True)
        if not self._quiet:
            print(phot_table)

        # Set official photometry table
        self._phot_table = phot_table
        
        # Create a copy of the full set of sources, also used to
        # get source statistics
        self._full_srclist = phot_table.copy()
        self._create_photometry_statistics()
                
        # Trim to the constructed value of max_sources
        if not dont_trim:
            self.trim(self._max_sources)
            
        # Number of sources with photometry
        self._nsrcs_photom = len(self._phot_table)
        return phot_table

    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        logger = logging.getLogger(__name__)
        
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
        
    def measure_fwhm(self, fwhm_plot_file, direction=None):
        """Fit gaussians to a sub-selection of identified stars in the
           center of the image and four outer quadrants, returning the
           median 1-D full width half maximum in pixels.
           
        This function in intended to measure the true FWHM of the stars
        within the image, returning a value that can be used to perform
        a better source searching and aperture photometry than the 
        default FWHM supplied at object instantiation. It is also intended
        to populate image quality metrics that may be useful in excluding
        images that have poor seeing when building stacked images.
        
        Only a representative subset of stars is used, drawn from the
        central half of the image plus four outer quadrants. This 
        partitioning is used to quantify variations in the size and
        symmetry of the stellar PSF across the image.
        
        The function delegates work to an instance of ApMeasureStars.
        
        If fwhm_plot_file is not None then a plot of image cut-outs
        around each selected star, along with the 2-D and 1-D best fit
        Gaussians, is created.
        
        Although the fitting always performs fitting along X and Y 
        rotated by an angle of theta, the direction
        parameter can be used to control the returned values.
        If direction = None or 'both' reports the values over both axes.
        This is the default as it has the clearest meaning.
        If direction = 'x' then only the median, error and number of points
        over the X-axis are reported.
        If direction is 'y' then only the median, error and number of points
        over the Y-axis are reported.
        Note that theta is not controlled, so 'x' and 'y' do not mean the same
        thing as the row/column axis of the image, rendering the returned
        values somewhat ambiguous.
        """
        
        # Construct a meaningful title for the plot file based on the
        # name of the original image file
        fname = os.path.basename(self._fitsimg)
        plot_title = f'Star FWHM measurements for:\n{fname}'
        
        # Note: Cannot use the initial sources table as it lacks
        #   accurate brightness estimates.
        measure_stars   = ApMeasureStars(self._data,
            self._phot_table,
            self._search_fwhm,
            self._bg_median,
            self._full_srclist,
            fwhm_plot_file,
            plot_title,
            self._loglevel,
            self._quiet)
        self._psf_table    = measure_stars.results_table()
        self._nsrcs_fitted = len(self._psf_table)
        
        # Get global FWHM (median over both X and Y)
        self._fwhm_both = measure_stars.median_fwhm('both') 
        self._fwhm_x    = measure_stars.median_fwhm('x') 
        self._fwhm_y    = measure_stars.median_fwhm('y') 
        
        self._logger.info(f'Median FWHM (over x and y) is {self._fwhm_both[0]:.2f} +/- {self._fwhm_both[1]:.2f} pixels using {self._fwhm_both[2]} data points')
        
        if direction is not None:
            if direction == 'both':
                result = self._fwhm_both
            elif direction == 'x':
                result = self._fwhm_x
            elif direction == 'y':
                result = self._fwhm_y
            else:
                err_msg = (
                    f'Error, unexpected direction={direction} passed to measure_fwhm.'
                    f' Expecting one of "None", "both", "x" or "y".'
                    )
                self._logger.error(err_msg)
        else:
            direction = 'both'
            result    = self._fwhm_both
        return result
        
    def _check_file_exists(self, filename):
        if not os.path.isfile(filename):
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            self._status = ApFindStars.INPUT_ERROR
            raise RuntimeError(err_msg)
        return
    
    def _read_fits(self, image_filename, image_extension):
        """Read a single extension's data and header from a FITS file
        """
        
        self._check_file_exists(image_filename)
        self._logger.info('Loading extension {} of FITS file {}'.format(image_extension, image_filename))
            
        # open() parameters that can be important.
        # Default values used here.
        # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
        uint_handling = True
        image_scaling = False
            
        with fits.open(image_filename, 
            uint=uint_handling, 
            do_not_scale_image_data=image_scaling) as hdu_list:
            ext_hdr  = hdu_list[image_extension].header
            ext_data = hdu_list[image_extension].data
            
        ndim     = ext_hdr['NAXIS']
        cols     = ext_hdr['NAXIS1']
        rows     = ext_hdr['NAXIS2']
        bitpix   = ext_hdr['BITPIX']
        info_str = '{}-D BITPIX={} image with {} columns, {} rows'.format(ndim, bitpix, cols, rows)
        
        if ndim == 3:
            layers = ext_hdr['NAXIS3']
            info_str = '{}-D BITPIX={} image with {} columns, {} rows, {} layers'.format(ndim, bitpix, cols, rows, layers)

        if 'BSCALE' in ext_hdr:
            bscale   = ext_hdr['BSCALE']
            info_str += f', BSCALE={bscale}'
        
        if 'BZERO' in ext_hdr:
            bzero    = ext_hdr['BZERO']
            info_str += f', BZERO={bzero}'
            
        self._logger.debug(info_str)
        
        if ndim == 3:
            self._logger.error('Error, 3-D handling has not been implemented yet.')
            sys.exit(1)
            
        # Get data absolute limits.
        minval = np.amin(ext_data)
        maxval = np.amax(ext_data)
        medval = np.median(ext_data)
        self._logger.debug(f'Raw data statistics are min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}')
        
        # Is there a PEDESTAL value? MaximDL likes to add an offset, and
        # the PEDESTAL value is the value to ADD to the data to remove the
        # pedestal.
        if 'PEDESTAL' in ext_hdr:
            pedestal = float( ext_hdr['PEDESTAL'] )
            if pedestal != 0:
                self._logger.debug(f'Removing a PEDESTAL value of {pedestal} ADU.')
                ext_data += pedestal
                minval = np.amin(ext_data)
                maxval = np.amax(ext_data)
                medval = np.median(ext_data)
                self._logger.debug(f'After PEDESTAL removal, min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}')
        
        return ext_data, ext_hdr
    
    def _write_source_list(self, 
            p_sourcelist,   # Name of utput FIT table with source info
            p_fitsimg,      # Name of input FITS image that was processed
            p_max_sources,  # Max number of sources for output table
            hdr,            # FITS header of input image
            kw_dict,        # Dictionary of keyword values to add to sourcelist header
            src_table,      # Astropy table of source photometry
            psf_table):     # None or Astropy table of source PSF fitting.   
        """Write detected source information to a FITS table with info on the
           original data file
           
        This function also:
         - Prints a summary of values useful for astrometry.net at INFO level.
        """
        
        # Create an X and Y table for use with astrometry.net,
        # adding 1 to get FITS like pixel indices
        self._logger.debug('Converting python 0-based coordinates to FITS 1-based pixel coordinates for XY table.')
        x = src_table['xcenter'] + 1.0 * u.Unit('pix')
        y = src_table['ycenter'] + 1.0 * u.Unit('pix')
        xy_table = Table([x, y],
            names=('X', 'Y'))
    
        # Currently just write trimmed FITS-format XY pos and merged position/photomtry
        self._logger.info('Writing source list to FITS binary table {}'.format(p_sourcelist))
    
        pri_hdr = self._create_primary_header(kw_dict)
        pri_hdu = fits.PrimaryHDU(header=pri_hdr)
        
        tbl_hdu1 = fits.table_to_hdu(xy_table)
        tbl_hdu1.header['EXTNAME'] = 'AP_XYPOS'
        tbl_hdu1.header['COMMENT'] = 'Uses FITS 1-based pixel coordinate system.'
        
        tbl_hdu2 = fits.table_to_hdu(src_table)
        tbl_hdu2.header['EXTNAME'] = 'AP_L1MAG'
        tbl_hdu2.header['COMMENT'] = 'Aperature photometry using photutils within ApFindStars.'
        tbl_hdu2.header['COMMENT'] = 'Uses python 0-based pixel coordinate system.'
    
        # Always add these to the output file as we know they're present
        hdus_to_append = [pri_hdu, tbl_hdu1, tbl_hdu2]
        
        if psf_table is not None:
            tbl_hdu3 = fits.table_to_hdu(psf_table)
            tbl_hdu3.header['EXTNAME'] = 'AP_L1PSF'
            tbl_hdu3.header['COMMENT'] = 'PSF characterization using ApMeasureStars.'
            tbl_hdu3.header['COMMENT'] = 'Uses python 0-based pixel coordinate system.'
            hdus_to_append.append(tbl_hdu3)
    
        hdu_list = fits.HDUList(hdus_to_append)
        hdu_list.writeto(p_sourcelist, overwrite=True)
    
        return
    
    def _create_photometry_statistics(self):
        """Calculates and stores the brightness (adu_per_sec) of the
           brightest, faintest, and median source from the full 
           photometric data table. 
           
        These values are written to the quality report.
        """
        
        # As the list is sorted by brightness, the source with index (N/2)
        # is the median, the index 0 zource the brightest, and the N-1
        # source the faintest.
        num_srcs = len(self._full_srclist)
        idx_bright = 0
        idx_median = int(num_srcs/2)
        idx_faint  = num_srcs - 1
        
        adups_bright = float( self._full_srclist['adu_per_sec'][idx_bright] )
        adups_median = float( self._full_srclist['adu_per_sec'][idx_median] )
        adups_faint  = float( self._full_srclist['adu_per_sec'][idx_faint]  )

        self._logger.debug(f'Brightest source (idx={idx_bright}) generates {adups_bright:.2f} ADU/s.')
        self._logger.debug(f'Median source (idx={idx_median}) generates {adups_median:.2f} ADU/s.')
        self._logger.debug(f'Faintest source (idx={idx_faint}) generates {adups_faint:.2f} ADU/s.')
        
        # The format of the stored info is a tuple of tuples 
        # ((adups_bright, idx_bright), (adups_median, idx_median), (adups_faint, idx_faint))

        self._phot_stats = ((adups_bright, idx_bright), 
            (adups_median, idx_median), 
            (adups_faint, idx_faint))
        return
    
    def _create_primary_header(self, kw_dict):
        """Creates a FITS primary header adding the keywords in the dict
        """
        
        pri_hdr = fits.Header()
        tnow    = datetime.now().isoformat(timespec='milliseconds')
        for kw in kw_dict:
            pri_hdr[kw] = kw_dict[kw]
        pri_hdr['HISTORY'] = f'Created by ApFindStars at {tnow}'
        return pri_hdr
        
    def _read_optional_keywords(self, hdr, kw_dict):
        """Read optional keywords from the image file FITS header hdr,
           and add to the kw_dict dictionary if they are present.
        """
        
        # The following may exist. We add comment values for them.
        kw_comment_dict = {'EXPOSURE': '[seconds] Image exposure time',
            'DATE-OBS': 'Observation date and time',
            'OBJECT':   'Target object', 
            'OBJNAME':  'Target object', 
            'TELESCOP': 'Telescope used',
            'INSTRUME': 'Detector used',
            'CCD-TEMP': 'CCD temperature at start of exposure in C',
            'APTDIA':   '[mm] Diameter of telescope aperture',
            'RA':       'Target right ascension', 
            'DEC':      'Target declination', 
            'RA-OBJ':   'Target right ascension', 
            'DEC-OBJ':  'Target declination', 
            'XPIXSZ':   '[micrometers] X-axis pixel scale after binning',
            'YPIXSZ':   '[micrometers] Y-axis pixel scale after binning',
            'FOCALLEN': '[mm] Stated telescope focal length', 
            'FILTER':   'Filter used', 
            'EGAIN':    '[e/ADU] Gain in electrons per ADU',
            'LAT-OBS':  '[deg +N WGS84] Observatory Geodetic latitude',
            'LONG-OBS': '[deg +E WGS84] Observatory Geodetic longitude',
            'ALT-OBS':  '[metres] Observatort altitude above mean sea level',
            'AIRMASS':  'Airmass (multiple of zenithal airmass)'}
        kw_missing_list = []
        for kw in kw_comment_dict:
            if kw in hdr:
                kw_dict[kw] = (hdr[kw], kw_comment_dict[kw])
            else:
                kw_missing_list.append(kw)
    
        self._logger.debug('FITS keywords found in image: {}'.format(kw_dict))
        self._logger.debug('FITS keywords missing from image: {}'.format(kw_missing_list))
        return
    
    def _build_keyword_dictionary(self, img_name, # Name of image
            hdr,            # FITS header associated with image
            bg_median,      # Estimated median BG level [ADU]
            bg_stddev):     # Estimated BG level standard deviation [ADU]
        """Build a dictionary of useful keywords and values to be
           used in the output sourcelist and optional quality report.
            
        The dictionary is of the form key: (value, comment). Numpy types
        are converted to native float or int.
        """
        
        # Get keywords we'll want to write for info purposes into
        # the output FITS file, and added information about the original
        # data file, the background, and the software used.
        # TODO check for KW presence and handle...
        
        kw_dict = {'IMG_FILE': (img_name, 'Name of image file searched for stars')}
        
        # These must exist by definition.
        cols         = int(hdr['NAXIS1'])
        rows         = int(hdr['NAXIS2'])
        kw_dict['IMG_COLS'] = (cols, 'Number of columns in input image')
        kw_dict['IMG_ROWS'] = (rows, 'Number of rows in input image')
        self._logger.info('Image is {} cols x {} rows'.format(cols, rows))
                
        # Number of sources detected, that are in the final photometry
        # and source lists, and were optionally used in estimating the
        # stellar FWHM
        kw_dict['AP_NDET']   = (self._nsrcs_detected, 'Number of sources detected in the image.')
        kw_dict['AP_NPHOT']  = (self._nsrcs_photom,   'Number of sources final photometry.')
        kw_dict['AP_NFIT']   = (self._nsrcs_fitted,   'Number of sources used in FWHM fitting.')
        kw_dict['AP_NSIGM'] = (self._search_nsigma,  'Source searching threshold (sigma above background)')
                
        # Add keywords that may or may not exist
        # Each value of the dict is a (value, comment) tuple.
        self._read_optional_keywords(hdr, kw_dict)
            
        # Coordinates
        angl_ra  = None
        angl_dec = None
        if 'RA' in kw_dict:
            raw_ra       = kw_dict['RA'][0]
            angl_ra      = Angle(raw_ra, unit=u.hourangle)
        if 'DEC' in hdr:
            raw_dec      = kw_dict['DEC'][0]
            angl_dec     = Angle(raw_dec, unit=u.deg)
        if (angl_ra is not None) and (angl_dec is not None):
            raw_coord    = SkyCoord(ra=angl_ra, dec=angl_dec)
            self._logger.info('Approximate coordinates: ra={:.6f} hours, dec={:.6f} deg'.format(raw_coord.ra.hour, raw_coord.dec.degree))
            kw_dict['APRX_RA']  = (raw_coord.ra.degree, '[deg] Approximate image center RA')
            kw_dict['APRX_DEC'] = (raw_coord.dec.degree, '[deg] Approximate image center Dec')
        
        # Pixel and image size
        if (('FOCALLEN' in kw_dict) and ('XPIXSZ' in kw_dict) and ('YPIXSZ' in kw_dict)):
            focal_len_mm  = float(kw_dict['FOCALLEN'][0])                      # mm
            pixsiz_x_um   = float(kw_dict['XPIXSZ'][0])                        # micrometers
            pixsiz_x_rad  = (pixsiz_x_um*1.0e-6) / (focal_len_mm*1.0e-3) # radians
            pixsiz_x_deg  = math.degrees(pixsiz_x_rad)
            pixsiz_x_arcs = 3600.0 * pixsiz_x_deg
    
            pixsiz_y_um   = float(kw_dict['YPIXSZ'][0])                        # micrometers
            pixsiz_y_rad  = (pixsiz_y_um*1.0e-6) / (focal_len_mm*1.0e-3) # radians
            pixsiz_y_deg  = math.degrees(pixsiz_y_rad)
            pixsiz_y_arcs = 3600.0 * pixsiz_y_deg
    
            imgsiz_x_deg  = cols * pixsiz_x_deg
            imgsiz_y_deg  = rows * pixsiz_y_deg
            imgsiz_deg    = math.sqrt( (imgsiz_x_deg * imgsiz_x_deg) + (imgsiz_y_deg * imgsiz_y_deg) )
        
            self._logger.info('Approximate image field of view is {:.3f} degrees across.'.format(imgsiz_deg))
            self._logger.info('Approximate pixel sizes (arcseconds) are x={:.3f}, y={:.3f}'.format(pixsiz_x_arcs, pixsiz_y_arcs))
            
            kw_dict['APRX_FOV'] = (imgsiz_deg,    '[deg] Approximate diagonal size of image')
            kw_dict['APRX_XWD'] = (imgsiz_x_deg,  '[deg] Approximate X-axis width of image')
            kw_dict['APRX_YHG'] = (imgsiz_y_deg,  '[deg] Approximate Y-axis height of image')
            kw_dict['APRX_XPS'] = (pixsiz_x_arcs, '[arcseconds] Approximate X-axis plate scale')
            kw_dict['APRX_YPS'] = (pixsiz_y_arcs, '[arcseconds] Approximate Y-axis plate scale')
        
        # If the FWHM and estimated error have been measured, add them to
        # the info we write to the primary header.
        if self._fwhm_both is not None:
            kw_dict['AP_FWHM']  = (self._fwhm_both[0], '[pix] Median FWHM of fitted stars in image')
            kw_dict['AP_EFWHM'] = (self._fwhm_both[1], '[pix] MAD standard deviation of fitted FWHM') 
        
        # Background levels are useful downstream
        kw_dict['AP_BGMED'] = (float( self._bg_median ), '[ADU] Median source-masked background level')
        kw_dict['AP_BGSTD'] = (float( self._bg_stddev ), '[ADU] Std dev of source-masked background level')
        
        return kw_dict
    
    def _trim_table(self, src_table, max_sources):
        """Trims the table to contain at maximum max_sources if max_sources
           is not None
        """
        
        # Create a copy, filtering the brightest if necessary...
        nsrc = len(src_table)
        nuse = nsrc
        if max_sources is not None:
            nuse = min(nsrc, max_sources)
        out_table = src_table[0:nuse]
        self._logger.debug('Selecting {} sources out of {} to write to the output source list.'.format(nuse, nsrc))
    
        return out_table
    
    def _find_saturated(self, data, sat_thresh, search_fwhm):
        """Use threshold search to look for regions of saturated pixels
        """
        
        boxsize = int(4 * search_fwhm)
        self._logger.debug(f'Looking for possibly saturated regions above {sat_thresh} ADU separated by {boxsize} pixels.')    
        saturated_positions = find_peaks(data, 
            threshold=sat_thresh, 
            box_size=boxsize)
        print('Saturated position table:\n', saturated_positions)
        return saturated_positions
        
    def write_ds9_region_file(self, region_file):
        """Write a ds9-format region file (image coordinates) of
           the stars within the photometry table.
           
        """
        
        # TODO pick color based on Filter
        reg_color = 'red'
        pix_fmt   =  '.2f' # 0.01 pixels (for degrees use .6f)
        
        num_stars = len(self._phot_table)
        ap_radius = math.ceil(self._ap_fwhm_mult * self._search_fwhm)
        regions = []
        for x,y,id in zip(self._phot_table['xcenter'], 
            self._phot_table['ycenter'],
            self._phot_table['id']):
            # There is no need to convert from python 0-based to FITS
            # 1-based as the region class(es) do this for us.
            xpos      = x
            ypos      = y
            meta_dict = {'name': id}
            vis_dict  = {'color': reg_color}
            cpr       = CirclePixelRegion(center=PixCoord(xpos, ypos), 
                radius=ap_radius, 
                meta=meta_dict, 
                visual=vis_dict)
            regions.append( cpr )
        
        #reg_str = ds9_objects_to_string(regions, 
            # coordsys='image',
            # fmt=pix_fmt,
            # radunit='pix')
        reg_str = regions.serialize(format='ds9')
        with open(region_file, 'w', encoding='utf-8') as f_out:
            f_out.write(reg_str)
        regions.write
            
        self._logger.debug(f'Wrote ds9-format region file to {region_file}')
        return    
            
    def write_quality_report(self, quality_report_name):
        """Write a summary of the input image source and background
           properties to a YaML file
           
        The quality report is a dictionary of dictionaries, containing:
        - basic input image properties
        - background level statistics
        - detected source information
        - saturation related statistics
        - FWHM statistics if present
        """
        
        # Null val for missing data
        null_val = -999
        
        # Slightly better than the default float representation
        yaml.add_representer(float, yaml_float_representer)
                
        # Build or rebuild keyword list useful for the FITS table as well as the
        # image quality report.
        self._kw_dict = self._build_keyword_dictionary(self._fitsimg,
            self._hdr,
            self._bg_median,
            self._bg_stddev)
            
        # Note that:
        # - Data within the kw_dict is key: (value, comment)
        # - JSON stores tuples as lists, so you can't retrieve a tuple.
    
        im_info_dict  = {}
        bg_info_dict  = {}
        src_info_dict = {}
        sat_info_dict = {}
        psf_info_dict = {}
        
        # Image information definition, used to fill the im_info_dict.
        im_info_struct_dict = {'file':     {'kw': 'IMG_FILE', 'fmt': 's'},
            'ncols':                {'kw': 'IMG_COLS', 'fmt': 'd'},
            'nrows':                {'kw': 'IMG_ROWS', 'fmt': 'd'},
            'object':               {'kw': 'OBJECT',   'fmt': 's'},
            'telescope':            {'kw': 'TELESCOP', 'fmt': 's'},
            'filter':               {'kw': 'FILTER',   'fmt': 's'},
            'date-obs':             {'kw': 'DATE-OBS', 'fmt': 's'},
            'exposure':             {'kw': 'EXPOSURE', 'fmt': '.2f'},
            'ccd_temperature':      {'kw': 'CCD-TEMP', 'fmt': '.2f'},
            'electronic_gain':      {'kw': 'EGAIN',    'fmt': '.4f'},
            'airmass':              {'kw': 'AIRMASS',  'fmt': '.2f'},
            'approx_width_deg':     {'kw': 'APRX_XWD', 'fmt': '.3f'},
            'approx_height_deg':    {'kw': 'APRX_YHG', 'fmt': '.3f'},
            'approx_xpixsiz_arcs':  {'kw': 'APRX_XPS', 'fmt': '.3f'},
            'approx_ypixsiz_arcs':  {'kw': 'APRX_YPS', 'fmt': '.3f'}}
        
        for okey in im_info_struct_dict:
            fkw = im_info_struct_dict[okey]['kw']
            fmt = im_info_struct_dict[okey]['fmt']
            if fkw in self._kw_dict:
                val = self._kw_dict[fkw][0] 
                # Skip using the fmt as yaml_float_representer does a decent job. 
                im_info_dict[okey] = val
            else:
                self._logger.warning(f'Keyword {fkw} not found, not written to quality report.')
                
        # Background info
        bgmed = self._kw_dict['AP_BGMED'][0]
        bgstd = self._kw_dict['AP_BGSTD'][0]
        bg_info_dict['median'] = bgmed
        bg_info_dict['stddev'] = bgstd
        
        # Source information
        src_info_dict['num_detected']        = self._kw_dict['AP_NDET'][0]
        src_info_dict['num_with_photometry'] = self._kw_dict['AP_NPHOT'][0]
        src_info_dict['search_nsigma']       = self._kw_dict['AP_NSIGM'][0]
        src_info_dict['adups_brightest']     = self._phot_stats[0][0]
        src_info_dict['adups_median']        = self._phot_stats[1][0]
        src_info_dict['adups_faintest']      = self._phot_stats[2][0]
        
        
        # Saturation info
        num_sat_in_phot = int( np.sum(self._phot_table['psbl_sat']==True) )
        sat_info_dict['num_saturated_in_image']      = self._nsrcs_saturated
        sat_info_dict['num_saturated_in_photometry'] = num_sat_in_phot
                
        # PSF aka star fitted FWHM info in units of pixels and arcseconds
        psf_info_dict['num_fit'] = self._kw_dict['AP_NFIT'][0]
        
        # Check have approximate plate scale
        have_platescale = False
        if ('APRX_XPS' in self._kw_dict) and ('APRX_YPS' in self._kw_dict):
            have_platescale = True
        else:
            self._logger.warning('Skipping some quality reporting that requires an estimate platescale.')
        
        
        
        # Do conversions to estimated arcseconds
        if self._psf_table is not None:
            if have_platescale:
                aprx_xpixsize = self._kw_dict['APRX_XPS'][0]
                aprx_ypixsize = self._kw_dict['APRX_YPS'][0]
                avg_pixsize   = math.sqrt(0.5 * 
                    (aprx_xpixsize**2 + aprx_ypixsize**2) )
            else:
                aprx_xpixsize = null_val
                aprx_ypixsize = null_val
                avg_pixsize   = null_val
            
            # Is fwhm_x within 3 sigma of fwhm_y
            fwhm_x    = self._fwhm_x[0]
            fwhm_xerr = self._fwhm_x[1]
            fwhm_y    = self._fwhm_y[0]
            fwhm_yerr = self._fwhm_y[1]
            circular  = ApMeasureStars.is_circular(fwhm_x, fwhm_y,
                fwhm_xerr, fwhm_yerr)
            psf_info_dict['circular_psf'] = circular
            
            # Now create three dictionaries for both x&y, x and y
            for direction in ['both', 'x', 'y']:
                if 'both' in direction:
                    result_tuple = self._fwhm_both
                    dict_name    = 'fwhm_xandy'
                    if have_platescale:
                        pixsiz_arcs  = avg_pixsize
                    else:
                        pixsiz_arcs  = null_val
                elif 'x' in direction:
                    result_tuple = self._fwhm_x
                    dict_name    = 'fwhm_x'
                    if have_platescale:
                        pixsiz_arcs  = aprx_xpixsize
                    else:
                        pixsiz_arcs  = null_val
                elif 'y' in direction:
                    result_tuple = self._fwhm_y
                    dict_name    = 'fwhm_y'
                    if have_platescale:
                        pixsiz_arcs  = aprx_ypixsize
                    else:
                        pixsiz_arcs  = null_val
                
                psf_dir_dict = {}
                fwhm_pix      = result_tuple[0]
                fwhm_err      = result_tuple[1]
                fwhm_val_arcs = fwhm_pix * pixsiz_arcs
                fwhm_err_arcs = fwhm_err * pixsiz_arcs
                psf_dir_dict['fwhm_val_pix']  = fwhm_pix
                psf_dir_dict['fwhm_err_pix']  = fwhm_err
                psf_dir_dict['fwhm_val_arcs'] = fwhm_val_arcs
                psf_dir_dict['fwhm_err_arcs'] = fwhm_err_arcs
                psf_dir_dict['num_data_pts']  = result_tuple[2]
                psf_info_dict[dict_name] = psf_dir_dict
                
        qual_dict = {}
        qual_dict['image_info']      = im_info_dict
        qual_dict['background_info'] = bg_info_dict
        qual_dict['source_info']     = src_info_dict
        qual_dict['saturation_info'] = sat_info_dict
        qual_dict['psf_info']        = psf_info_dict
        
        with open(quality_report_name, 'w') as write_file:
            yaml.dump(qual_dict, write_file, indent=4, sort_keys=False)
        self._logger.info(f'Wrote image quality report to {quality_report_name}')
        return
