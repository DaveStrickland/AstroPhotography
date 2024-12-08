#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_combine_darks.py
#
#  Generates a master dark or master bias file, given a directory 
#  that contains the raw dark or bias file to be combined. Use separate
#  directories for calcilbration files with different exposure time,
#  binning, CCD temperature, and from different telescopes.
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
#  2020-10-01 dks : Initial coding begun. 
#  2020-10-11 dks : Working version (astropy still spams logger).
#  2024-10-01 dks : Moved ApMasterCal into core

import argparse
import AstroPhotography as ap
import logging

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_combine_darks',
        description='Generates a master dark or master bias file from' + 
        ' all calibration FITS files in a given directory.')
        
    # Required
    parser.add_argument('rawcaldir',
        metavar='RAW_CAL_DIR',
        help='The directory in which the raw calibration' +
        ' files to be combined are to be found.')
    parser.add_argument('master_filename',
        metavar='MASTER_CAL_FILENAME',
        help='Output file name for master calibration file.')
        
    # Default values
    p_temptol  = 0.5 # 0.5 degree C
    p_telescop = 'UNKNOWN'
    p_exclude  = 'master*'
        
    # Optional
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
    parser.add_argument('--exclude',
        dest='exclude_pattern',
        default=p_exclude,
        metavar='FILE_PATTERN',
        help=('A unix-style file pattern that may be used to exclude'
            ' files in the target directory from being processed by this'
            ' command. Usually this is used to exclude any master'
            ' calibration files from being processed.'
            f' Default: "{p_exclude}"'))
    parser.add_argument('--telescop',
        default=p_telescop,
        metavar='TELESCOPE_NAME',
        help=('If the input files TELESCOP keyword is missing or empty,' +
            ' write this string as TELESCOP in the output master calibration' +
            f' file, e.g. "iTelescope 5". Default: {p_telescop}'))
    parser.add_argument('--temptol',
        default=p_temptol,
        metavar='DEGREES_C',
        help=('Allowable temperature tolerance up to which the CCDTEMP'
            ' can vary from the SET-TEMP and still be considered as'
            f' representing that temperature. Default: {p_temptol} C.'))
    
    args = parser.parse_args(argv)
    return args
                             
def main(args=None):
    p_args       = command_line_opts(args)
    p_root       = p_args.rawcaldir
    p_masterfile = p_args.master_filename
    p_loglevel   = p_args.loglevel
    p_telescop   = p_args.telescop
    p_temptol    = p_args.temptol
    p_exclude    = p_args.exclude_pattern

    logger = logging.getLogger(__name__)

    try:
        mkcal = ap.ApMasterCal(p_root, 
            p_exclude,
            p_telescop, 
            p_temptol,
            p_loglevel)
        mkcal.make_master(p_masterfile)
    except RuntimeError as rte:
        logger.error(f'Shutting down due to exception raised by ApMasterCal: {rte}')
        return 1
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
