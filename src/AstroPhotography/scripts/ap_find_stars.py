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

import argparse
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
import math

from astropy.io import fits
from astropy.table import QTable, Table
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.visualization import AsymmetricPercentileInterval
from astropy.visualization import SqrtStretch
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import sigma_clipped_stats
from astropy import units as u

from photutils import make_source_mask  
from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils import aperture_photometry

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
    parser.add_argument('-d', '--ds9',
        default=None,
        metavar='ds9.reg',
        help='Name for optional ds9-format region file.')
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

    # TODO move astrometry.et stuff to separate script
    parser.add_argument('-k', '--key',
        default=None,
        type=str,
        metavar='ASTROMETRY_API_KEY',
        help='Your personal Astrometry.net API key.' + 
        ' If supplied an astrometry.net solution will be attempted.')
            
    args = parser.parse_args(argv)
    return args

def initialize_logger(loglevel):
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

def read_fits(image_filename, image_extension):
    """Read a single extension's data and header from a FITS file
    """

    logger = logging.getLogger(__name__)
    
    # TODO check for file being present.
    
    logger.info('Loading extension {} of FITS file {}'.format(image_extension, image_filename))
        
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
        
    logger.debug(info_str)
    if ndim == 3:
        logger.error('Error, 3-D handling has not been implemented yet.')
        sys.exit(1)
    return ext_data, ext_hdr

def write_source_list(p_sourcelist, p_fitsimg, 
        p_max_sources, p_regfile, hdr,
        bg_median, bg_stddev,
        src_table,
        p_astnetkey):
    """Write detected source information to a FITS table with info on the
       original data file
       
    This function also:
     - Prints a summary of values useful for astrometry.net at INFO level.
     - Writes the ds9-format region file if p_regfile is not None.
    """
    
    logger = logging.getLogger(__name__)
    
    # Create an X and Y table for use with astrometry.net,
    # adding 1 to get FITS like pixel indices
    logger.debug('Converting python 0-based coordinates to FITS 1-based pixel coordinates for XY table.')
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
    logger.info('Image is {} cols x {} rows'.format(cols, rows))
    
    # Add keywords that may or may not exist
    # Each value of the dict is a (value, comment) tuple.
    add_optional_keywords(hdr, kw_dict)
        
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
        logger.info('Approximate coordinates: ra={:.6f} hours, dec={:.6f} deg'.format(raw_coord.ra.hour, raw_coord.dec.degree))
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
    
        logger.info('Approximate image field of view is {:.3f} degrees across.'.format(imgsiz_deg))
        logger.info('Approximate pixel sizes (arcseconds) are x={:.3f}, y={:.3f}'.format(pixsiz_x_arcs, pixsiz_y_arcs))
        
        kw_dict['APRX_FOV'] = (imgsiz_deg, '[deg] Approximate diagonal size of image')
        kw_dict['APRX_XSZ'] = (pixsiz_x_arcs, '[arcseconds] Approximate X-axis plate scale')
        kw_dict['APRX_YSZ'] = (pixsiz_y_arcs, '[arcseconds] Approximate Y-axis plate scale')
    
    # TODO move to another script
    if p_astnetkey is not None:
        wcs = query_astrometry_dot_net(p_astnetkey,
            cols, 
            rows, 
            Table(xy_table))
        if len(wcs) > 0:
            logger.debug('Non-zero length WCS header returned')
            print(wcs)
        else:
            logger.warn('Zero length WCS returned. Solution failed')
    
    
    # Currently just do Simplistic write with no added header
    logger.info('Writing source list to FITS binary table {}'.format(p_sourcelist))

    pri_hdr = create_primary_header(kw_dict)
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

def query_astrometry_dot_net(p_astnetkey, cols, rows, xy_table):
    """Attempt get astrometric solution"""

    logger = logging.getLogger(__name__)    
    wcs = {}
    image_width = 3073
    image_height = 2048
    
    from astroquery.astrometry_net import AstrometryNet

    ast = AstrometryNet()
    ast.api_key = p_astnetkey
    try_again = True
    submission_id = None

    while try_again:
        try:
            if not submission_id:
                logger.debug('Submitting astrometry.net solve from source list with {} poisitions'.format(len(xy_table)))
                wcs_header = ast.solve_from_source_list(xy_table['X'], 
                    xy_table['Y'],
                    image_width, 
                    image_height,
                    parity=2,
                    positional_error=10,
                    crpix_center=True,
                    publicly_visible='n',
                    submission_id=submission_id)
            else:
                logger.debug('Monitoring astrometry.net submission {}'.format(submission_id))
                wcs_header = ast.monitor_submission(submission_id,
                    solve_timeout=120)
        except TimeoutError as e:
            submission_id = e.args[1]
        else:
            # got a result, so terminate
            try_again = False
    
    if wcs_header:
        logger.info('Obtained a WCS header')
        wcs = wcs_header
    else:
        logger.warn('Astrometry.net submission={} failed'.format(submission_id))

    return wcs

def create_primary_header(kw_dict):
    """Creates a FITS primary header adding the keywords in the dict"""
    
    pri_hdr = fits.Header()
    for kw in kw_dict:
        pri_hdr[kw] = kw_dict[kw]
    pri_hdr['HISTORY'] = 'Created by ApFindStars'
    return pri_hdr
    
def add_optional_keywords(hdr, kw_dict):
    """Add keywords to kw_dict from the FITS header hdr if present"""
    logger = logging.getLogger(__name__)
    
    # The following may exist.
    kw_comment_dict = {'EXPOSURE': '[seconds] Image exposure time',
        'DATE-OBS': 'Observation date and time',
        'OBJECT':   'Target object', 
        'TELESCOP': 'Telescope used',
        'RA':       'Requested right ascension', 
        'DEC':      'Requested declination', 
        'XPIXSZ':   '[micrometers] Stated X-axis pixel scale after binning',
        'YPIXSZ':   '[micrometers] Stated Y-axis pixel scale after binning',
        'FOCALLEN': '[mm] Stated telescope focal length', 
        'FILTER':   'Filter used', 
        'EGAIN':    '[e/ADU] Gain in electrons per ADU'}
    kw_missing_list = []
    for kw in kw_comment_dict:
        if kw in hdr:
            kw_dict[kw] = (hdr[kw], kw_comment_dict[kw])
        else:
            kw_missing_list.append(kw)

    logger.debug('FITS keywords found in image: {}'.format(kw_dict))
    logger.debug('FITS keywords missing from image: {}'.format(kw_missing_list))
    return

def trim_table(src_table, max_sources):
    """Trims the table to contain at maximum max_sources if max_sources
       is not None
    """
    
    logger = logging.getLogger(__name__)
    
    # Create a copy, filtering the brightest if necessary...
    nsrc = len(src_table)
    nuse = nsrc
    if max_sources is not None:
        nuse = min(nsrc, max_sources)
    out_table = src_table[0:nuse]
    logger.debug('Selecting {} sources out of {} to write to the output source list.'.format(nuse, nsrc))

    return out_table

def main(args=None):
    p_args      = command_line_opts(args)
    p_loglevel  = p_args.loglevel          # Loglevel
    p_fitsimg   = p_args.fits_image        # Input fits image
    p_extnum    = p_args.fits_extension    # Extension number for data in input fits
    p_fitstbl   = p_args.source_list       # Output fits table of srcs
    p_regfile   = p_args.ds9               # Name of ds9 region file or None
    p_maxsrcs   = p_args.max_sources       # Max number of sources to output or None
    p_astnetkey = p_args.key               # Astrometry.net API key or None

    search_fwhm   = p_args.search_fwhm
    search_nsigma = p_args.search_nsigma
    detector_bitdepth = p_args.bitdepth
    sat_frac          = p_args.sat_frac
    
    logger    = initialize_logger(p_loglevel)   
    
    data, hdr = read_fits(p_fitsimg, p_extnum)
    
    pct_interval = AsymmetricPercentileInterval(0.50, 99.5)
    sqrt_norm    = ImageNormalize(data, 
        interval=pct_interval, 
        stretch=SqrtStretch())
    asinh_norm   = ImageNormalize(data,
        interval=pct_interval,
        stretch=AsinhStretch())
    
    # TODO should replace plt use with proper fig, ax etc...

    # Get initial background estimates...
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    logger.debug('Sigma clipped image stats: mean={:.3f}, median={:.3f}, stddev={:.3f}'.format(mean, median, std))  

    mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
    logger.debug('Source-masked image stats: mean={:.3f}, median={:.3f}, stddev={:.3f}'.format(mean, median, std))  
    
    # TODO better estimate of FWHM, sigma to use
    logger.debug(f'Running DAOStarFinder with FWHM={search_fwhm:.2f} pixels, threshold={search_nsigma} BG sigma.')
    daofind = DAOStarFinder(fwhm=search_fwhm, threshold=search_nsigma*std)  
    sources = daofind(data - median)  
    for col in sources.colnames:  
        sources[col].info.format = '%.8g'  # for consistent table output
    print('Initial source list using FHWM={} pixels, threshold={} x BG stddev'.format(search_fwhm, search_nsigma))
    print(sources)
    
    # TODO better estimate of aperture to use
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    ap_radius = math.ceil(1.5 * search_fwhm)
    logger.debug(f'Radius of circular aperture photometry is {ap_radius} pixels.')
    apertures = CircularAperture(positions, r=ap_radius)
    
    # Use initial source list for aperture photometry
    phot_table = aperture_photometry(data -median, apertures)
    phot_table['aperture_sum'].info.format = '%.4f'  # for consistent table output
    
    # Copy over peak flux from star finder output and use it to derive
    # possible saturation flag.
    max_adu           = math.pow(2, detector_bitdepth) - 1
    sat_thresh        = math.floor(sat_frac * max_adu)
    logger.debug(f'Sources with pixels exceding {sat_thresh} ADU will be flagged as possibly saturated.')
    phot_table['peak_adu'] = sources['peak']
    phot_table['peak_adu'].info.format = '%.2f'
    phot_table['psbl_sat'] = phot_table['peak_adu'] > sat_thresh
    n_sat = np.sum(phot_table['psbl_sat'])
    logger.debug(f'There are {n_sat} possibly saturated stars in this image.')
    
    exposure = float(hdr['EXPOSURE'])
    phot_table['adu_per_sec'] = phot_table['aperture_sum'] / exposure
    phot_table['magnitude']   = -2.5 * np.log10(phot_table['adu_per_sec'])
    
    phot_table['adu_per_sec'].info.format = '%.4f'
    phot_table['magnitude'].info.format   = '%.4f'
    phot_table['xcenter'].info.format     = '%.2f'
    phot_table['ycenter'].info.format     = '%.2f'

    phot_table.sort(['adu_per_sec', 'xcenter', 'ycenter'], reverse=True)
    print(phot_table)

    phot_table = trim_table(phot_table, p_maxsrcs)

    positions = np.transpose((phot_table['xcenter'], phot_table['ycenter']))
    apertures = CircularAperture(positions, r=4.)

    plt.imshow(data, origin='lower', norm=asinh_norm)
    apertures.plot(color='red', lw=1.5, alpha=0.5)
    
    plt.show()
    
    # Write source list or subset to a FITS table
    # with expected FITS-style coordinates...
    write_source_list(p_fitstbl, 
        p_fitsimg,
        p_maxsrcs, 
        p_regfile, 
        hdr, 
        median, 
        std,
        phot_table,
        p_astnetkey)
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
