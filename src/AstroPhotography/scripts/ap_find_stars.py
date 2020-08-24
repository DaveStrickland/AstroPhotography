#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ap_find_stars.py
#  
#  Copyright 2020 Dave Strickland <dave.strickland@gmail.com>
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
# 2020-08-23 dks : Add stellar FWHM fitting and plotting.

import argparse
import sys
import logging
import os.path
import numpy as np
import matplotlib.pyplot as plt
import math

from astropy.io import fits
from astropy.table import QTable, Table, vstack
from astropy.coordinates import SkyCoord, Angle
from astropy.visualization import AsymmetricPercentileInterval
from astropy.visualization import SqrtStretch, AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import sigma_clipped_stats
from astropy import units as u

from photutils import make_source_mask, find_peaks, DAOStarFinder
from photutils import CircularAperture, aperture_photometry

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_find_stars',
        description='Find stars in a FITS image using photutils')
        
    # Required
    parser.add_argument('fits_image',
        help='Input FITS image file name to process.')
    parser.add_argument('source_list',
        help='File name for output FITS table of detected sources.')
        
    # Optional
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
    parser.add_argument('-e','--fits_extension',
        default=0,
        type=int,
        metavar='EXT_NUM',
        help='FITS extension number to read from. Default=0 (Primary)')
    parser.add_argument('-m', '--max_sources',
        default=None,
        type=int,
        metavar='NUM_SRCS',
        help='Limit output source list the brightest set of NUM_SRCS sources.' + 
        ' This also applies to the ds9 format region file, if any.'
        ' Default: Output all the sources.')
        
    def_search_fwhm   = 3.0
    def_search_nsigma = 7.0
    def_bitdepth      = 16
    def_sat_frac      = 0.80
    parser.add_argument('--search_fwhm',
        default=def_search_fwhm,
        type=float,
        help=f'Initial source search FWHM (pixels). Default: {def_search_fwhm}')
    parser.add_argument('--search_nsigma',
        default=def_search_nsigma,
        type=float,
        help='Initial source search threshold in numbers of sigma above the background (standard deviations).' + 
        f' Default: {def_search_nsigma}')
    parser.add_argument('--bitdepth',
        default=def_bitdepth,
        type=int,
        help=f'Detector bitdepth used in saturation calculation. Default: {def_bitdepth}')
    parser.add_argument('--sat_frac',
        default=def_sat_frac,
        type=float,
        help=f'Fraction of max ADU used in saturation calculation. Default: {def_sat_frac}')
    parser.add_argument('--retain_saturated', 
        action='store_true', 
        help='Do not exclude possibly saturated stars from the source list.' + 
        ' By default pixels within a box of width 8x search_fwhm centered on possibly' +
        ' saturated stars are excluded from source detection and fitting.')
    parser.add_argument('--plotfile',
        default=None,
        metavar='IMG_WITH_SRCS.PNG',
        help='If specified, a scaled bitmap of input image with the detected' +
        ' sources circled will be written to disk.')
    parser.add_argument('--quality_report',
        default=None,
        metavar='QUALITY_REPORT.TXT',
        help="If specified, an ASCII file summarizing the image's source" +
        " and background properties will be written to disk.")
    parser.add_argument('--fwhm_plot',
        default=None,
        metavar='FWHM_FITS.PNG',
        help='If specified, a bitmap plot of selected star cut-outs' +
        ' along with fitted gaussian model parameters such as FWHM' +
        ' will be written to disk.' + 
        ' (The fit results are output in the source list FITS table,' +
        ' and summarized in the qulity report.)')
    parser.add_argument('-d', '--ds9',
        default=None,
        metavar='ds9.reg',
        help='Name for optional ds9-format region file.')
    parser.add_argument('-q', '--quiet',
        action='store_true',
        default=False,
        help='Quiet mode suppresses the default mode that prints of a' +
        ' summary of the detected source lists to STDOUT while running.')
            
    args = parser.parse_args(argv)
    return args

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
        self._psf_table     = None

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
        num_sat_candidates  = len(saturated_locations)
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
        
        # Search for stars using the supplied FWHM and threshold.
        self._sources = self.source_search(self._search_fwhm, self._search_nsigma)
        self._apertures = self._make_apertures(self._sources['xcentroid'], 
            self._sources['ycentroid'])
        
        # Use initial source list for aperture photometry
        self._phot_table = self.aperature_photometry()
    
        # Trim and update apertures
        self._phot_table = self._trim_table(self._phot_table, 
            self._max_sources)
        self._apertures  = self._make_apertures(self._phot_table['xcenter'], 
            self._phot_table['ycenter'])
    
        if self._plotfile is not None:
            self._plot_image(self._plotfile)

        return
        
    def _plot_image(self, plotfile):
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
        """Create an apertures object for aperture photometry
        """

        # TODO better estimate of aperture to use
        positions = np.transpose((colpos, rowpos))
        ap_radius = math.ceil(1.5 * self._search_fwhm)
        self._logger.debug(f'Radius of circular aperture photometry is {ap_radius} pixels.')
        apertures = CircularAperture(positions, r=ap_radius)
        
        return apertures

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
        return sources

    def write_source_list(self, output_fits_table):
        """Write the current source position and photometry data to
        a FITS table, overwriting any existing file
        """
        self._write_source_list(output_fits_table, 
            self._fitsimg,
            self._max_sources, 
            self._hdr, 
            self._bg_median, 
            self._bg_stddev,
            self._phot_table)
        return

    def aperature_photometry(self):
        """Perform aperature photometry using the existing sources
           and apertures.
        """
        
        phot_table = aperture_photometry(self._data - self._bg_median, 
            self._apertures)
        phot_table['aperture_sum'].info.format = '%.4f'  # for consistent table output
        
        # Copy over peak ADU and possibly saturated flag from initial
        # sourcelist 
        phot_table['peak_adu'] = self._sources['peak']
        phot_table['psbl_sat'] = self._sources['psbl_sat']
        
        exposure = float(self._hdr['EXPOSURE'])
        phot_table['adu_per_sec'] = phot_table['aperture_sum'] / exposure
        phot_table['magnitude']   = -2.5 * np.log10(phot_table['adu_per_sec'])
        
        phot_table['adu_per_sec'].info.format = '%.4f'
        phot_table['magnitude'].info.format   = '%.4f'
        phot_table['xcenter'].info.format     = '%.2f'
        phot_table['ycenter'].info.format     = '%.2f'
    
        phot_table.sort(['adu_per_sec', 'xcenter', 'ycenter'], reverse=True)
        if not self._quiet:
            print(phot_table)
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
        
    def measure_fwhm(self, fwhm_plot_file):
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
        """
        
        # Note: Cannot use the initial sources table as it lacks
        #   accurate brightness estimates.
        measure_stars   = ApMeasureStars(self._data,
            self._phot_table,
            self._search_fwhm,
            self._bg_median,
            fwhm_plot_file,
            self._loglevel,
            self._quiet)
        self._psf_table = measure_stars.results_table()
        median_fwhm     = measure_stars.median_fwhm() 
        return median_fwhm
        
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
        bzero    = ext_hdr['BZERO']
        bscale   = ext_hdr['BSCALE']
        info_str = '{}-D BITPIX={} image with {} columns, {} rows, BSCALE={}, BZERO={}'.format(ndim, bitpix, cols, rows, bscale, bzero)
        
        if ndim == 3:
            layers = ext_hdr['NAXIS3']
            info_str = '{}-D BITPIX={} image with {} columns, {} rows, {} layers, BSCALE={}, BZERO={}'.format(ndim, bitpix, cols, rows, layers, bscale, bzero)
            
        self._logger.debug(info_str)
        if ndim == 3:
            self._loggererror('Error, 3-D handling has not been implemented yet.')
            sys.exit(1)
        return ext_data, ext_hdr
    
    def _write_source_list(self, p_sourcelist, p_fitsimg, 
            p_max_sources, hdr,
            bg_median, bg_stddev,
            src_table):
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
        xy_table = QTable([x, y],
            names=('X', 'Y'))
    
        # Get keywords we'll want to write for info purposes into
        # the output FITS file, and added information about the original
        # data file, the background, and the software used.
        # TODO check for KW presence and handle...
        
        kw_dict = {'IMG_FILE': (p_fitsimg, 'Name of image file searched for stars')}
        
        # These must exist by definition.
        cols         = int(hdr['NAXIS1'])
        rows         = int(hdr['NAXIS2'])
        kw_dict['IMG_COLS'] = (cols, 'Number of columns in input image')
        kw_dict['IMG_ROWS'] = (rows, 'Number of rows in input image')
        kw_dict['BGMEDIAN'] = (bg_median, '[ADU] Median BG level in sigma clipped masked image')
        kw_dict['BGSTDDEV'] = (bg_stddev, '[ADU] Std dev of BG level')
        self._logger.info('Image is {} cols x {} rows'.format(cols, rows))
        
        # Add keywords that may or may not exist
        # Each value of the dict is a (value, comment) tuple.
        self._add_optional_keywords(hdr, kw_dict)
            
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
            
            kw_dict['APRX_FOV'] = (imgsiz_deg, '[deg] Approximate diagonal size of image')
            kw_dict['APRX_XSZ'] = (pixsiz_x_arcs, '[arcseconds] Approximate X-axis plate scale')
            kw_dict['APRX_YSZ'] = (pixsiz_y_arcs, '[arcseconds] Approximate Y-axis plate scale')
        
        # Currently just write trimmed FITS-format XY pos and merged position/photomtry
        self._logger.info('Writing source list to FITS binary table {}'.format(p_sourcelist))
    
        pri_hdr = self._create_primary_header(kw_dict)
        pri_hdu = fits.PrimaryHDU(header=pri_hdr)
        
        tbl_hdu1 = fits.table_to_hdu(xy_table)
        tbl_hdu1.header['EXTNAME'] = 'AP_XYPOS'
        tbl_hdu1.header['COMMENT'] = 'Uses FITS 1-based pixel coordinate system.'
        
        tbl_hdu2 = fits.table_to_hdu(src_table)
        tbl_hdu2.header['EXTNAME'] = 'AP_L1MAG'
        tbl_hdu2.header['COMMENT'] = 'Uses python 0-based pixel coordinate system.'
    
        hdu_list = fits.HDUList([pri_hdu, tbl_hdu1, tbl_hdu2])
        hdu_list.writeto(p_sourcelist, overwrite=True)
    
        return
    
    def _create_primary_header(self, kw_dict):
        """Creates a FITS primary header adding the keywords in the dict
        """
        
        pri_hdr = fits.Header()
        for kw in kw_dict:
            pri_hdr[kw] = kw_dict[kw]
        pri_hdr['HISTORY'] = 'Created by ApFindStars'
        return pri_hdr
        
    def _add_optional_keywords(self, hdr, kw_dict):
        """Add keywords to kw_dict from the FITS header hdr if present
        """
        
        # The following may exist.
        kw_comment_dict = {'EXPOSURE': '[seconds] Image exposure time',
            'DATE-OBS': 'Observation date and time',
            'OBJECT':   'Target object', 
            'TELESCOP': 'Telescope used',
            'RA':       'Requested right ascension', 
            'DEC':      'Requested declination', 
            'XPIXSZ':   '[micrometers] X-axis pixel scale after binning',
            'YPIXSZ':   '[micrometers] Y-axis pixel scale after binning',
            'FOCALLEN': '[mm] Stated telescope focal length', 
            'FILTER':   'Filter used', 
            'EGAIN':    '[e/ADU] Gain in electrons per ADU'}
        kw_missing_list = []
        for kw in kw_comment_dict:
            if kw in hdr:
                kw_dict[kw] = (hdr[kw], kw_comment_dict[kw])
            else:
                kw_missing_list.append(kw)
    
        self._logger.debug('FITS keywords found in image: {}'.format(kw_dict))
        self._logger.debug('FITS keywords missing from image: {}'.format(kw_missing_list))
        return
    
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
        
class ApMeasureStars:
    """Measures stellar PSF size and symmetry in selected stars across
       and input image.
    """
    
    def __init__(self,
            img_data,
            srclist,
            init_fwhm,
            init_bglevel,
            fwhm_plot_file,
            loglevel,
            quiet):
        """ApMeasureStars constructor
        """
        
        self._img_data   = img_data        # Image data 2-D array
        self._init_fwhm  = init_fwhm       # Initial guess at FWHM in pixels
        self._init_bglvl = init_bglevel    # Initial guess at BG level per pixel
        self._plot_file  = fwhm_plot_file  # None or name of file to plot results to.
        self._loglevel   = loglevel        # Logging level
        self._quiet      = quiet           # If True then table summary printing disabled.

        self._fit_table  = None            # Table to store fit results in.
                
        # Set up logging
        self._logger = self._initialize_logger(self._loglevel)

        # Setting related to candidate selection.
        self._num_per_reg    = 5           # Number of sources per region to fit
        self._skip_brightest = 2           # Skip the brightest N stars in each region
        self._logger.debug(f'Up to {self:_num_per_reg} stars per sub-region will be fitted, excluding the brightest {self:_skip_brightest} stars.')

        # Slightly modify the input source list to exclude possibly
        # saturated stars, remove unneeded columns.
        self._init_srcs  = Table(srclist[srclist['psbl_sat']==False])         # Initial sourclist from ApFindStars
        #self._init_srcs.remove_columns(['sharpness', 'roundness1', 'roundness2', 'flux'])
        if not self._quiet:
            print('Input table supplied to ApMeasureStars after saturated star filtering:')
            print(self._init_srcs)
            print('')

        self._cols       = img_data.shape[1]
        self._rows       = img_data.shape[0]
        
        # Set up variables related to fit box size and edge exclusion
        self._fit_box_initialization()
        
        # Select candidates and determine data extraction boxes.
        self._fit_table  = self._select_candidates()
        self._calculate_boxes()
        
        return
        
    def _calculate_boxes(self):
        """Calculate the pixel indices of the boxes from which data will
           be extracted for fitting each star.
        """
        
        half_width = self._box_width_pix / 2
        nrst_xpix  = np.rint( self._fit_table['xcenter'] )
        nrst_ypix  = np.rint( self._fit_table['ycenter'] )
        self._fit_table['xmin'] = nrst_xpix - half_width
        self._fit_table['xmax'] = nrst_xpix + half_width
        self._fit_table['ymin'] = nrst_ypix - half_width
        self._fit_table['ymax'] = nrst_ypix + half_width
        
        # Make the column printing reasonable
        col_fmt_dict = {'xmin': '%d',
            'xmax': '%d',
            'ymin': '%d',
            'ymax': '%d'}
        for col in col_fmt_dict:
            self._fit_table[col].info.format = col_fmt_dict[col]  
        
        if not self._quiet:
            print('Fit candidates with data extraction box indices:')
            print(self._fit_table)
            print('')
        return
        
    def _fit_box_initialization(self):
        """Sets up variables related to the size of the region used in
           fitting and the edge exclusion used in candidate selection.
        """
        
        # We want the fit box to be at least 2x the FWHM, and an even 
        # number of pixels. We also pad a bit in case the initial FWHM
        # estimate is an underestimate.
        pad_frac            = 2.5
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
        
    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
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
        
    def _select_candidates(self):
        """Select candidate stars in the center and four quadrants
        
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
        A_Q = H/4 * (W - (pi*H/16))
        
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
        """
        
        # Algorithm:
        # Iterate over each source i in the input source list
        #   If not saturated, calculate dx = x_i - x_c, dy = y_i - y_c, r_i
        #   If r_i < H/2, add source id, mag to center list
        #   else if dx < 0, dy > 0 add to TL list, etc
        # The sort lists by magnitude
        # Pick the top N from each quadrant and place in a table
        # The table forms the basis for self._fit_table
        
        candidate_table = None
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
            cand_table.sort(['magnitude'], reverse=True)
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
        
    def median_fwhm(self):
        """Returns the median fitted 1-D FWHM in pixels over the entire image.
        """
        
        fwhm_med = None 
        return fwhm_med
        
    def results_table(self):
        """Return the fitting results as an astropy Table
        """
        return self._fit_table

def main(args=None):
    p_args      = command_line_opts(args)
    p_loglevel  = p_args.loglevel          # Loglevel
    p_fitsimg   = p_args.fits_image        # Input fits image
    p_extnum    = p_args.fits_extension    # Extension number for data in input fits
    p_fitstbl   = p_args.source_list       # Output fits table of srcs

    p_search_fwhm       = p_args.search_fwhm
    p_search_nsigma     = p_args.search_nsigma
    p_detector_bitdepth = p_args.bitdepth
    p_sat_frac          = p_args.sat_frac

    p_max_sources = p_args.max_sources       # Max number of sources to output or None
    p_nosatmask   = p_args.retain_saturated  # Keep possibly saturated stars.
    
    p_regfile     = p_args.ds9               # Name of ds9 region file or None
    p_plotfile    = p_args.plotfile          # None or name of output png image w sources.
    p_qual_rprt   = p_args.quality_report    # None or name of output ascii report
    p_quiet       = p_args.quiet             # Quiet mode suppresses STDOUT list printing

    p_fwhm_plot   = p_args.fwhm_plot         # PNG/JPG plot of PSF fits 

    # Perform initial source detection using default parameters.
    find_stars = ApFindStars(p_fitsimg, p_extnum, p_search_fwhm,
        p_search_nsigma, p_detector_bitdepth, 
        p_max_sources, p_nosatmask, p_sat_frac, p_loglevel,
        p_plotfile, p_quiet)
    
    # Measure 2-Gaussian FWHM for select stars
    p_new_fwhm = find_stars.measure_fwhm(p_fwhm_plot)
    
    # Refine source detection
    #TODO#find_stars.source_search(p_new_fwhm, p_search_nsigma)
    
    # Re-run photometry
    #TODO#find_stars.aperture_photmetry()
    
    # Write optional quality report
    #TODO#if p_qual_rprt is not None:
    #TODO#    find_stars.write_quality_report(p_qual_rprt)

    # Write optional ds9 format region file
    #TODO#if p_regfile is not None:
    #TODO#    find_stars.write_ds9_region_file(p_regfile)
    
    # Write final sourcelist with photometry.
    find_stars.write_source_list(p_fitstbl)
        
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
