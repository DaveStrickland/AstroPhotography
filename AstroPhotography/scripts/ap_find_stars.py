#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_find_stars.py
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
# 2021-01-18 dks : Move classes into core.

import argparse
import sys
import logging
import os.path
import time

import AstroPhotography as ap

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
        help="If specified, an YaML file summarizing the image's source" +
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
    
def main(args=None):
    perf_time_start = time.perf_counter() # highest res timer, counts sleeps
    prcs_time_start = time.process_time() # usr+system, excludes sleeps
    
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
    find_stars = ap.ApFindStars(p_fitsimg, p_extnum, p_search_fwhm,
        p_search_nsigma, p_detector_bitdepth, 
        p_max_sources, p_nosatmask, p_sat_frac, p_loglevel,
        p_plotfile, p_quiet)
    
    # Measure 2-Gaussian FWHM for select stars, get average over x and y
    (p_new_fwhm, p_madstd_fwhm, p_npts) = find_stars.measure_fwhm(p_fwhm_plot, 'both')
    
    # Refine source detection
    find_stars.source_search(p_new_fwhm, p_search_nsigma)
    
    # Re-run photometry
    find_stars.aperture_photometry()
    
    # As the source searching and photometry was redone, we should redo
    # the plotting.
    find_stars.plot_image(p_plotfile)
    
    # Write optional quality report
    if p_qual_rprt is not None:
        find_stars.write_quality_report(p_qual_rprt)

    # Write optional ds9 format region file
    if p_regfile is not None:
        find_stars.write_ds9_region_file(p_regfile)
    
    # Write final sourcelist with photometry.
    find_stars.write_source_list(p_fitstbl)
        
    # Assumed logger was set up elsewhere
    perf_time_end = time.perf_counter() # highest res timer, counts sleeps
    prcs_time_end = time.process_time() # usr+system, excludes sleeps
    perf_time     = perf_time_end - perf_time_start
    prcs_time     = prcs_time_end - prcs_time_start
    logging.getLogger(__name__).info(f'Elapsed wall clock runtime {perf_time:.3f} s (process time {prcs_time:.3f} s).')
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
cd 
