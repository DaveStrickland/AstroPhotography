# -*- coding: utf-8 -*-
#
#  Contains the implementation of ApMeasureBackground
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

# 2020-05-29 dks : Initial coding of skeleton

import logging
import os.path
import numpy as np
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

from regions import PixCoord, CirclePixelRegion, ds9_objects_to_string

from photutils import make_source_mask, find_peaks, DAOStarFinder
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
    
class ApMeasureBackground:
    """
    Measures and outputs large scale non-uniform sky backgrounds left 
    by imperfect bias/dark/flat calibration.
    """

    def __init__(self, loglevel):
        """Constructor for ApMeasureBackground        
        """
        
        self._loglevel      = loglevel

        # Set up logging
        self._logger = self._initialize_logger(self._loglevel)

        return
        
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
        
