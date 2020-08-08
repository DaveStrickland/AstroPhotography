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

import argparse
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
import math

from astropy.io import fits
from astropy.table import QTable
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
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

def write_source_list(cli_args, hdr,
        bg_median, bg_stddev,
        src_table):
    """Write detected source information to a FITS table with info on the
       original data file, and optionally only outputing the brightest stars.
       
    This function also:
     - Prints a summary of values useful for astrometry.net at INFO level.
     - Writes the ds9-format region file if cli_args.ds9 is not None.
    """
    
    logger = logging.getLogger(__name__)
    
    # Create a copy, filtering the brightest if necessary...
    nsrc = len(src_table)
    nuse = nsrc
    if cli_args.max_sources is not None:
        nuse = min(nsrc, cli_args.max_sources)
    out_table = src_table[0:nuse]
    logger.debug('Selecting {} sources out of {} to write to the output source list.'.format(nuse, nsrc))

    # Create an X and Y table for use with astrometry.net,
    # adding 1 to get FITS like pixel indices
    logger.debug('Converting python 0-based coordinates the FITS 1-based pixel coordinates.')
    x = out_table['xcenter'] + 1.0 * u.Unit('pix')
    y = out_table['ycenter'] + 1.0 * u.Unit('pix')
    xy_table = QTable([x, y],
        names=('X', 'Y'))

    # Get keywords we'll want to write for info purposes into
    # the output FITS file, and added information about the original
    # data file, the background, and the software used.
    # TODO check for KW presence and handle...
    
    kw_dict = {'IMG_FILE': cli_args.fits_image}
    
    # These must exist by definition.
    cols         = int(hdr['NAXIS1'])
    rows         = int(hdr['NAXIS2'])
    kw_dict['IMG_COLS'] = cols
    kw_dict['IMG_ROWS'] = rows
    logger.info('Image is {} cols x {} rows'.format(cols, rows))
    
    # Add keywords that may or may not exist
    add_optional_keywords(hdr, kw_dict)
        
    # Coordinates
    angl_ra  = None
    angl_dec = None
    if 'RA' in kw_dict:
        raw_ra       = kw_dict['RA']
        angl_ra      = Angle(raw_ra, unit=u.hourangle)
    if 'DEC' in hdr:
        raw_dec      = kw_dict['DEC']
        angl_dec     = Angle(raw_dec, unit=u.deg)
    if (angl_ra is not None) and (angl_dec is not None):
        raw_coord    = SkyCoord(ra=angl_ra, dec=angl_dec)
        logger.info('Approximate coordinates: ra={:.6f} hours, dec={:.6f} deg'.format(raw_coord.ra.hour, raw_coord.dec.degree))
        kw_dict['APRX_RA']  = raw_coord.ra.hour
        kw_dict['APRX_DEC'] = raw_coord.dec.degree
    
    # Pixel and image size
    if (('FOCALLEN' in kw_dict) and ('XPIXSX' in kw_dict) and ('YPIXSZ' in kw_dict)):
        focal_len_mm  = float(kw_dict['FOCALLEN'])                      # mm
        pixsiz_x_um   = float(kw_dict['XPIXSZ'])                        # micrometers
        pixsiz_x_rad  = (pixsiz_x_um*1.0e-6) / (focal_len_mm*1.0e-3) # radians
        pixsiz_x_deg  = math.degrees(pixsiz_x_rad)
        pixsiz_x_arcs = 3600.0 * pixsiz_x_deg

        pixsiz_y_um   = float(kw_dict['YPIXSZ'])                        # micrometers
        pixsiz_y_rad  = (pixsiz_y_um*1.0e-6) / (focal_len_mm*1.0e-3) # radians
        pixsiz_y_deg  = math.degrees(pixsiz_y_rad)
        pixsiz_y_arcs = 3600.0 * pixsiz_y_deg

        imgsiz_x_deg  = cols * pixsiz_x_deg
        imgsiz_y_deg  = rows * pixsiz_y_deg
        imgsiz_deg    = math.sqrt( (imgsiz_x_deg * imgsiz_x_deg) + (imgsiz_y_deg * imgsiz_y_deg) )
    
        logger.info('Image is approximately {:.3f} degrees across.'.format(imgsiz_deg))
        logger.info('Pixel size (arcseconds) x={:.3f}, y={:.3f}'.format(pixsiz_x_arcs, pixsiz_y_arcs))
        
        kw_dict['APRX_FOV'] = imgsiz_deg
        kw_dict['APRX_XSZ'] = pixsiz_x_arcs
        kw_dict['APRX_YSZ'] = pixsiz_y_arcs
    
    # Currently just do Simplistic write with no added header
    logger.info('Writing source list to FITS binary table {}'.format(cli_args.source_list))

    pri_hdr = create_primary_header(kw_dict)
    pri_hdu = fits.PrimaryHDU(header=pri_hdr)
    
    tbl_hdu1 = fits.table_to_hdu(xy_table)
    tbl_hdu1.header['EXTNAME'] = 'AP_XYPOS'
    tbl_hdu1.header['COMMENT'] = 'Uses FITS 1-based pixel coordinate system.'
    
    tbl_hdu2 = fits.table_to_hdu(out_table)
    tbl_hdu2.header['EXTNAME'] = 'AP_L1MAG'
    tbl_hdu2.header['COMMENT'] = 'Uses python 0-based pixel coordinate system.'

    hdu_list = fits.HDUList([pri_hdu, tbl_hdu1, tbl_hdu2])
    hdu_list.writeto(cli_args.source_list, overwrite=True)

    return
    
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
    kw_list = ['EXPOSURE', 'DATE-OBS', 'OBJECT', 'TELESCOP',
        'RA', 'DEC', 'XPIXSZ', 'YPIXSZ',
        'FOCALLEN', 'FILTER', 'EGAIN']
    kw_missing_list = []
    for kw in kw_list:
        if kw in hdr:
            kw_dict[kw] = hdr[kw]
        else:
            kw_missing_list.append(kw)

    logger.debug('FITS keywords found in image: {}'.format(kw_dict))
    logger.debug('FITS keywords missing from image: {}'.format(kw_missing_list))
    return

def main(args=None):
    p_args    = command_line_opts(args)
    logger    = initialize_logger(p_args.loglevel)    
    data, hdr = read_fits(p_args.fits_image, p_args.fits_extension)
    
    sqrt_norm  = ImageNormalize(stretch=SqrtStretch())
    asinh_norm = ImageNormalize(stretch=AsinhStretch())
    
    # TODO should replace plt use with proper fig, ax etc...

    # Get initial background estimates...
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    logger.debug('Sigma clipped image stats: mean={:.3f}, median={:.3f}, stddev={:.3f}'.format(mean, median, std))  

    mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
    logger.debug('Source-masked image stats: mean={:.3f}, median={:.3f}, stddev={:.3f}'.format(mean, median, std))  
    
    # TODO better estimate of FWHM, sigma to use
    search_fwhm   = 3.0
    search_nsigma = 7.0
    daofind = DAOStarFinder(fwhm=search_fwhm, threshold=search_nsigma*std)  
    sources = daofind(data - median)  
    for col in sources.colnames:  
        sources[col].info.format = '%.8g'  # for consistent table output
    print('Initial source list using FHWM={} pixels, threshold={} x BG stddev'.format(search_fwhm, search_nsigma))
    print(sources)
    
    # TODO better estimate of aperture to use
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=4.)

    plt.imshow(data, cmap='Greys', origin='lower', norm=asinh_norm)
    apertures.plot(color='red', lw=1.5, alpha=0.5)
    
    # Use initial source list for aperture photometry
    phot_table = aperture_photometry(data -median, apertures)
    phot_table['aperture_sum'].info.format = '%.4f'  # for consistent table output
    print(phot_table)
    
    exposure = float(hdr['EXPOSURE'])
    phot_table['adu_per_sec'] = phot_table['aperture_sum'] / exposure
    phot_table['magnitude']   = -2.5 * np.log10(phot_table['adu_per_sec'])
    
    phot_table['adu_per_sec'].info.format = '%.4f'
    phot_table['magnitude'].info.format   = '%.4f'
    phot_table['xcenter'].info.format     = '%.2f'
    phot_table['ycenter'].info.format     = '%.2f'

    phot_table.sort(['adu_per_sec', 'xcenter', 'ycenter'], reverse=True)
    print(phot_table)
    
    plt.show()
    
    # Write source list or subset to a FITS table
    # with expected FITS-style coordinates...
    write_source_list(p_args, hdr, 
        median, std,
        phot_table)
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
