#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ap_calc_read_noise.py
#
#  Calculates an estimate of the detector read noise (electrons/pixel)
#  given two raw bias images and the electronic gain (electrons/ADU).
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
#  2020-09-16 dks : Initial coding. 
# 

import argparse
import sys
import logging
from pathlib import Path

import numpy as np
from astropy.io import fits

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_calc_rread_noise',
        description='Calculates an estimate of the detector read noise' + 
        ' (e/pixel) given two raw bias images and the electronic gain' +
        ' (e/ADU).')
        
    # Required
    parser.add_argument('biasfile1',
        metavar='RAW_BIAS_1.FITS',
        help='Path/name of the first raw bias file to read.')
    parser.add_argument('biasfile2',
        metavar='RAW_BIAS_2.FITS',
        help='Path/name of the second raw bias file to read.')
        
    # Defaults.
    p_gain = 'EGAIN'
        
    # Optional
    parser.add_argument('--gain',
        default=p_gain,
        help=('FITS keyword listing the electronic gain (electrons/ADU)'
        ' within the first bias file, OR the numerical value of the'
        f' gain to use. Default: {p_gain}'))
    parser.add_argument('--noclip',
        dest='sigma_clip',
        action='store_false',
        default=True,
        help='If specified, extreme values of the bias pixel-to-pixel' +
        ' differences will NOT be removed using sigma clipping.')
    parser.add_argument('--histplot',
        default=None,
        help='If specified, the path/name to use for a histogram plot' +
        ' that this script will generate of the pixel-to-pixel' +
        ' differences between the two bias files.')
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    return args

class ApImageDifference:
    """TODO
    """
    
    def __init__(self,
        loglevel):
    
        # Initialize logging
        self._loglevel = loglevel
        self._initialize_logger(self._loglevel)
        
        return
                
    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        self._logger = logging.getLogger('ApImageDifference')
        
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
        self._logger.addHandler(ch)
        return
                
def main(args=None):
    p_args      = command_line_opts(args)
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
