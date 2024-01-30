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

# 2021-05-29 dks : Initial coding of skeleton
# 2021-08-20 dks : Add configurable BG box parameters, update documentation
# 2024-01-25 dks : Catch up to latest astropy/photutils changes

import logging
import os.path
import numpy as np
import math
import time
import yaml
from datetime import datetime
from pathlib import Path

from scipy import spatial

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip, mad_std
from astropy.stats import SigmaClip

from regions import PixCoord, CirclePixelRegion

from photutils.segmentation import detect_threshold, detect_sources
from photutils import find_peaks, DAOStarFinder
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.background import Background2D, MedianBackground

from .. import __version__

class ApMeasureBackground:
    """
    Measures and outputs large scale non-uniform sky backgrounds left 
    by imperfect bias/dark/flat calibration.
    
    The input image is split into a configurable number of large-scale
    boxes, and the sky background is estimated in each box using the 
    photutils background estimator after identified stars have been 
    excluded. The values from the large-scale boxes are then used to
    generate a smooth background estimate over the entire image, 
    excluding boxes that are outliers.
    
    The number of rows and columns of background boxes is configurable
    and the user is encouraged to carefully choose the numbers based on
    the images to be processed. Experimentation may be required in the
    case of narrow-band imaging of nebulae to avoid the nebula 
    emission from being considered part of the background.
    """

    def __init__(self, loglevel):
        """Constructor for ApMeasureBackground        
        """
        
        self._name      = 'ApMeasureBackground'
        self._loglevel  = loglevel
        
        # Standard constants, may be modified by user when processing
        # Based on H-alpha image of the Cygnus Loop with strong 
        # differential BG gradient (Moon?)
        self._default_minwdth = 48
        self._default_minhght = 48
        self._default_nbgcols = 16
        self._default_nbgrows = 16
        self._default_filtsiz = 3       # issue-021 Astropy 6 wants odd sized filter
        self._default_pctlile = 25.0
        self._default_nsigma  = 3.0
        self._boxsize         = (self._default_minhght, self._default_minwdth)
        self._filtersize      = (self._default_filtsiz, self._default_filtsiz)
        self._exclude_pctile  = self._default_pctlile
        self._nsigma          = self._default_nsigma
        
        # Data from last processing
        self._imdata    = None
        self._imhdr     = None
        self._bgdata    = None

        # Set up logging
        self._logger = self._initialize_logger(self._loglevel)

        return
        
    def _check_file_exists(self, filename):
        """Checks the file exists and cleans up the path
        """
        
        fpath = Path(filename).expanduser() 
        if not fpath.exists():
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return fpath

        
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
    
    def _make_source_mask(self):
        """
        Creates a source mask using image smoothing, segmentation, and 
        source dilation.
        
        The background estimation depends on excluding genuine large
        pointlike sources such as big galaxies and bright stars
        that are NOT well detected/excluded by the methods used in
        ApFindStars. Luckily photutils provies another method, which
        is good for this purpose.
        """
        
        sigma_clip  = SigmaClip(sigma=3.0, maxiters=10)
        threshold   = detect_threshold(self._imdata, nsigma=2.0, sigma_clip=sigma_clip)
        segment_img = detect_sources(self._imdata, threshold, npixels=5)
        mask        = segment_img.make_source_mask(size=13)
        
        ##mask = make_source_mask(self._imdata, nsigma=2, npixels=5, 
        ##    filter_fwhm=2.0, filter_size=5, 
        ##    dilate_size=13)
        
        # Print some statistics about the mask.
        npix = mask.size
        npos = np.sum(mask)
        pct  = 100 * float(npos)/float(npix)
        self._logger.debug(f'Generated segmentation source mask covering {pct:.3f}% of the image ({npos}/{npix} pixels)')
        
        #import matplotlib.pyplot as plt
        #plt.imshow(mask, cmap='Greys', origin='lower')
        #plt.show()
        return mask
        
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
        if 'BZERO' in ext_hdr:
            bzero    = ext_hdr['BZERO']
        else:
            bzero = 0
        if 'BSCALE' in ext_hdr:
            bscale   = ext_hdr['BSCALE']
        else:
            bscale = 1.0
        info_str = '{}-D BITPIX={} image with {} columns, {} rows, BSCALE={}, BZERO={}'.format(ndim, bitpix, cols, rows, bscale, bzero)
        
        if ndim == 3:
            layers = ext_hdr['NAXIS3']
            info_str = '{}-D BITPIX={} image with {} columns, {} rows, {} layers, BSCALE={}, BZERO={}'.format(ndim, bitpix, cols, rows, layers, bscale, bzero)
            
        self._logger.debug(info_str)
        if ndim == 3:
            self._loggererror('Error, 3-D handling has not been implemented yet.')
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

    def _remove_pedestal_kw(self, hdr):
        """Removes the PEDESTAL keyword from the input FITS header
        
        AstroPhotography always removes any artificial PEDESTAL applied
        to the data when reading a FITS file, so it is important to make
        sure that the FITS header keywords remain consistent.
        
        This function need only be applied when modified data is being
        written or rewritten to disk using a copy of an original FITS
        header.
        """
        
        if 'PEDESTAL' in hdr:
            self._logger.debug('Removing PEDESTAL keyword from FITS header.')
            del hdr['PEDESTAL']
        
        return 
        
    def _set_bgbox_size(self, imrows, imcols, 
            nbg_rows, nbg_cols, min_bgheight, min_bgwidth):
        """
        Calculate the size of the box or cell used for background 
        estimation based on the image size, the class defaults, and any
        user-specified constraints.
        
        :param imrows: Number of rows of data in the input image.
        :param imcols: Number of columns of data in the input image.
        :param nbg_rows: The number of rows of BG boxes the input image
          will be split into to determine the local sky background.
          If not specified a default number will be used.
        :param nbg_cols: The number of columns of BG boxes the input image
          will be split into to determine the local sky background.
          If not specified a default number will be used.
        :param min_bgheight: Minimum height in pixels of a background cell.
          The height of a background region in an image of NROWS rows
          is max(MIN_BGHEIGHT, 2*int(NROWS/2*NUM_BGROWS) pixels.
          If not specified a default number will be used.
        :param min_bgwidth: Minimum width in pixels of a background cell.
          The width of a background region in an image of NCOLS rows
          is max(MIN_BGWIDTH, 2*int(NCOLS/2*NUM_BGCOLS) pixels.
          If not specified a default number will be used.
        """

        actual_nbg_rows = self._default_nbgrows
        if nbg_rows is not None:
            actual_nbg_rows = nbg_rows

        actual_nbg_cols = self._default_nbgcols
        if nbg_cols is not None:
            actual_nbg_cols = nbg_cols

        actual_min_bgheight = self._default_minhght
        if min_bgheight is not None:
            actual_min_bgheight = min_bgheight

        actual_min_bgwidth = self._default_minwdth
        if min_bgwidth is not None:
            actual_min_bgwidth = min_bgheight
            
        self._logger.debug(('Background box size constraints:'
            f' {actual_nbg_rows} rows x {actual_nbg_cols} columns of background boxes.'
            f' Minimum box size {actual_min_bgheight} pixels high x {actual_min_bgwidth} pixels wide.'))

        # Makes the box sizes a multiple of the quantum value (e.g. 2, 8 etc)
        quantum = 2
        box_height = max(actual_min_bgheight, 
            quantum*(1 + int(imrows/(quantum*actual_nbg_rows)) ) )
        box_width = max(actual_min_bgwidth, 
            quantum*(1 + int(imcols/(quantum*actual_nbg_cols)) ) )
            
        # Now make sure we don't end up with a set of small boxes along
        # the margins.
        nrows = actual_nbg_rows*box_height
        ncols = actual_nbg_cols*box_width
        if nrows < imrows:
            # Should only be missing out by < quantum * num_bg rows
            if (imrows - nrows) > quantum*actual_nbg_rows:
                self._logger.error(f'Error in _set_bgbox_size logic: missing {imrows-nrows} rows using box height of {box_height} rows.')
            else:
                box_height += quantum
        
        if ncols < imcols:
            # Should only be missing out by < quantum * num_bg_cols
            if (imcols - ncols) > quantum*actual_nbg_cols:
                self._logger.error(f'Error in _set_bgbox_size logic: missing {imcols-ncols} cols using box width of {box_width} cols.')
            else:
                box_width += quantum
        
        self._logger.info(f'Background box size is {box_height} pixels high x {box_width} pixels wide.')
        
        self._box_height = box_height
        self._box_width  = box_width
        self._boxsize    = (box_height, box_width)
        return
        
    def process_data(self, imdata, imhdr=None,
        nbg_rows=None, 
        nbg_cols=None,
        min_bgheight=None,
        min_bgwidth=None,
        bg_filter_width=None,
        bg_badbox_pctile=None,
        bg_sigmaclip=None):
        """
        Measure the 2-dimensional sky background in data from an input
        FITS image.
        
        By default sources such as stars and galaxies are automatically
        detected and excluded from the background model.
        An optional FITS file containing source list table can be used
        *instead* in cases where the simple automated source detection
        does not work well.
        
        :param imdata: 2-dimensional numpy array containing the image
          from which we want to estimate the background.
        :param imhdr: Optional FITS header associated with the image
          imdata, if any. This will be used to populate any output
          background image's header metadata.
        :param nbg_rows: The number of rows of BG boxes the input image
          will be split into to determine the local sky background.
          If not specified a default number will be used.
        :param nbg_cols: The number of columns of BG boxes the input image
          will be split into to determine the local sky background.
          If not specified a default number will be used.
        :param min_bgheight: Minimum height in pixels of a background cell.
          The height of a background region in an image of NROWS rows
          is max(MIN_BGHEIGHT, 2*int(NROWS/2*NUM_BGROWS) pixels.
          If not specified a default number will be used.
        :param min_bgwidth: Minimum width in pixels of a background cell.
          The width of a background region in an image of NCOLS rows
          is max(MIN_BGWIDTH, 2*int(NCOLS/2*NUM_BGCOLS) pixels.
          If not specified a default number will be used.
        :param bg_filter_width: Photutils Background2D filter size
          parameter, used to median filter the course background values.
        :param bg_bg_badbox_pctile: Photutils Background2D exclude_percentile
          parameter, used to mark the course background cells to bad if
          more than this percentage of their pixels are excluded or masked.
        :param bg_sigmaclip: Photutils Background2D sigma_clip parameter,
          used to sigma clipping on the pixel values within each course
          background box.
        """
        
        self._imdata = imdata
        self._imhdr  = imhdr
        
        srcmask = self._make_source_mask()
        
        if bg_filter_width is not None:
            self._filtersize = (bg_filter_width, bg_filter_width)
        
        if bg_badbox_pctile is not None:
            self._exclude_pctile = bg_badbox_pctile
        
        if bg_sigmaclip is not None:
            self._nsigma = bg_sigmaclip
        
        sigma_clipper = SigmaClip(sigma=self._nsigma)
        bkg_estimator = MedianBackground()
        
        self._set_bgbox_size(imdata.shape[0], imdata.shape[1], 
            nbg_rows, nbg_cols, min_bgheight, min_bgwidth)
        
        self._logger.debug(('Background estimator settings used: estimator=MedianBackground'
            f', clipper=SigmaClip with sigma={self._nsigma}'
            f', boxsize={self._boxsize}, filter_size={self._filtersize}'
            f', exclude_percentile={self._exclude_pctile}'))
        
        bkg = Background2D(self._imdata, 
            self._boxsize,
            filter_size=self._filtersize,
            mask=srcmask,
            exclude_percentile=self._exclude_pctile,
            sigma_clip=sigma_clipper, 
            bkg_estimator=bkg_estimator)

        self._bgdata       = bkg.background
        self._bgmedian     = bkg.background_median
        self._bgmedian_rms = bkg.background_rms_median
        self._logger.info(f'Estimated median background level: {self._bgmedian:.3f}+/-{self._bgmedian_rms:.3f}')
        self._logger.debug(f'Background mesh shape: {bkg.mesh_nmasked.shape}')        
        return
        
    def process_files(self, input_fits, 
        srclist_fits=None,
        nbg_rows=None, 
        nbg_cols=None,
        min_bgheight=None,
        min_bgwidth=None,
        bg_filter_width=None,
        bg_badbox_pctile=None,
        bg_sigmaclip=None):
        """
        Measure the 2-dimensional sky background in data from an input
        FITS image.
        
        By default sources such as stars and galaxies are automatically
        detected and excluded from the background model.
        An optional FITS file containing source list table can be used
        *instead* in cases where the simple automated source detection
        does not work well.
        
        :param input_fits: Input FITS image containing the data.
        :param srclist_fits: Input FITS table containing a list of 
          detected sources, as generated by ap_find_stars.py.
        :param nbg_rows: The number of rows of BG boxes the input image
          will be split into to determine the local sky background.
          If not specified a default number will be used.
        :param nbg_cols: The number of columns of BG boxes the input image
          will be split into to determine the local sky background.
          If not specified a default number will be used.
        :param min_bgheight: Minimum height in pixels of a background cell.
          The height of a background region in an image of NROWS rows
          is max(MIN_BGHEIGHT, 2*int(NROWS/2*NUM_BGROWS) pixels.
          If not specified a default number will be used.
        :param min_bgwidth: Minimum width in pixels of a background cell.
          The width of a background region in an image of NCOLS rows
          is max(MIN_BGWIDTH, 2*int(NCOLS/2*NUM_BGCOLS) pixels.
          If not specified a default number will be used.
        :param bg_filter_width: Photutils Background2D filter size
          parameter, used to median filter the course background values.
        :param bg_bg_badbox_pctile: Photutils Background2D exclude_percentile
          parameter, used to mark the course background cells to bad if
          more than this percentage of their pixels are excluded or masked.
        :param bg_sigmaclip: Photutils Background2D sigma_clip parameter,
          used to sigma clipping on the pixel values within each course
          background box.
        """
        
        data, hdr = self._read_fits(input_fits, 0)
        
        # Generate source mask
        if srclist_fits is not None:
            self._logger.warning('Source list processing not yet implemented.')
            
        # Process the data.
        self.process_data(data, hdr, 
            nbg_rows, nbg_cols, min_bgheight, min_bgwidth,
            bg_filter_width, bg_badbox_pctile, bg_sigmaclip)
        return
        
    def get_bgimage(self):
        """
        Returns the background image calculated by process_data().
        """
        return self._bgdata
        
    def write_bgimage(self, output_bgfits):
        """
        Write the final 2-dimensional background map to a FITS file. 
        
        The header of the original input image will be used to populate
        the metadata of this file.
        
        :param output_bgfits: Name/path of the FITS file that the 2-D
          background image will be written to. This will be over-written
          if it already exists.
        """
        
        if self._bgdata is None:
            raise RuntimeError(f'Error, you can not write a background image before generating one using process_data or process_files.')
        
        # create hdulist
        ext_num  = 0
        hdu      = fits.PrimaryHDU(self._bgdata, header=self._imhdr)
        hdu_list = fits.HDUList([hdu])
                    
        self._remove_pedestal_kw(hdu_list[ext_num].header)
        
        tnow = datetime.now().isoformat(timespec='milliseconds')
        hdu_list[ext_num].header['HISTORY']  = f'Applied {self._name} {__version__} at {tnow}'
        hdu_list[ext_num].header['IMAGETYP'] = 'Background Sky'          
                    
        # Write to new file
        hdu_list.writeto(output_bgfits, 
            output_verify='ignore',
            overwrite=True)
        hdu_list.close()
    
        self._logger.info(f'Wrote estimated background data to {output_bgfits}')
        
        return
        
