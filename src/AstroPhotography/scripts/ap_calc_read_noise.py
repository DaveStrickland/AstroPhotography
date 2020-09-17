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
        dest='sigmaclip',
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
        
class ApCalcReadNoise:
    """Calculates an estimate of detector read noise (e/pixel) given
       two raw bias files and the electronic gain (e/ADU).
       
    Given two bias frames, B1 and B2, and the difference between them
    B1-B2, the standard deviation of the difference image is related
    to the read noise of the detector by
    
    Read Noise = Gain * (sigma_(b1-b2)) / sqrt(2)
       
    See section 4.3 of "Handbook of CCD Astronomy", Howell, S.B.,
    2000 (Cambridge University Press, Cambridge UK)
    """
    
    def __init__(self,
        biasfile1,
        biasfile2,
        gain,
        loglevel):
        """The constructor stores the names of the files and gain
           keyword or value, and initializes the logger, but does not
           perform any file IO or perform calculations.
           
           The actual file IO and calculations are only performed when
           the user calls estimate_rn.
        """
    
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
        
    def estimate_rn(self, sigmaclip, histplot=None):
        """Estimate the read noise from the files, and optionally
           produce a histograme plot of the pixel-to-pixel ADU
           differences. 
        
        If sigmaclip is true pixels with outlier values in either
        bias image will be excluded from the pixel-to-pixel differencing.
        Note that outliers are identified in each bias file separately
        in terms of absolute intensity, and the resulting bad pixel maps
        are combined.
        """
        
        
        
        return

        
                
def main(args=None):
    p_args      = command_line_opts(args)
    p_biasfile1 = p_args.biasfile1
    p_biasfile2 = p_args.biasfile2
    p_gain      = p_args.gain
    p_sigmaclip = p_args.sigmaclip
    p_histplot  = p_args.histplot
    p_loglevel  = p_args.loglevel
    
    read_noise_calculator = ApCalcReadNoise(p_biasfile1,
        p_biasfile2,
        p_gain,
        p_loglevel)
    
    rn1 = read_noise_calculator.estimate_rn(p_sigmaclip, p_histplot)
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
