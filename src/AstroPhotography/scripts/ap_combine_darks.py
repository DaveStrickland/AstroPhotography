#!/usr/bin/env python
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

import argparse
import sys
import logging
from pathlib import Path
import math

import numpy as np
import matplotlib                # for rc
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import mad_std

import ccdproc as ccdp

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
        default=None,
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
             
class ApMasterCal:
    """Combines a series of darks, bias or flats into a master dark,
       master bias, or master flat.
    """
    
    def __init__(self, rootdir, 
        exclude_pattern,
        telescop, 
        temptol,
        loglevel):
        """Construct and fully initialize an instance of ApMasterCal
        """

        # Initialize logging
        self._loglevel = loglevel
        self._initialize_logger(loglevel)

        # Set keywords that will be used to sort and check files.
        self._summary_kw = ['file', 'date-obs',
            'telescop', 'imagetyp', 
            'filter', 'exptime', 
            'set-temp', 'ccd-temp', 
            'naxis1', 'naxis2']

        # Generate raw ImageFileCollection
        self._rootdir  = rootdir
        self._data_dir = Path(rootdir)        
        self._files = ccdp.ImageFileCollection(self._data_dir, 
            keywords=self._summary_kw,
            glob_exclude=exclude_pattern)
    
        self._logger.info(f'Found {len(self._files.summary)} files in {rootdir} matching pattern.')
        print(self._files.summary)


        for kw in self._summary_kw:
            if 'file' not in kw:
                uniq_val_list = self._files.values(kw, unique=True)
                msg = f'For keyword {kw} there are {len(uniq_val_list)} values: {uniq_val_list}'
                self._logger.debug(msg)                
    
    
        # Check and filter files
        self._telescop = telescop
        self._temptol  = temptol
    
        return
             
    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        self._logger = logging.getLogger('ApMasterCal')
        
        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: {}'.format(loglevel))
        self._logger.setLevel(numeric_level)
    
        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(numeric_level)
    
        # create formatter
        formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s | %(message)s')
    
        # add formatter to ch
        ch.setFormatter(formatter)
    
        # add ch to logger
        if not self._logger.handlers:
            self._logger.addHandler(ch)
        
        # Used in cases where we get the same message twice or more
        # See https://stackoverflow.com/a/44426266
        self._logger.propagate = False
        return
                
def main(args=None):
    p_args       = command_line_opts(args)
    p_root       = p_args.rawcaldir
    p_masterfile = p_args.master_filename
    p_loglevel   = p_args.loglevel
    p_telescop   = p_args.telescop
    p_temptol    = p_args.temptol
    p_exclude    = p_args.exclude_pattern

    mkcal = ApMasterCal(p_root, 
        p_exclude,
        p_telescop, 
        p_temptol,
        p_loglevel)
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
