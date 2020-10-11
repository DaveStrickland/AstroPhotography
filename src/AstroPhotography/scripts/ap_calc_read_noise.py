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
#  2020-09-27 dks : Added plotting.

import argparse
import sys
import logging
import os.path
import math

import numpy as np
import matplotlib                # for rc
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_calc_read_noise',
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
    """Takes the difference between two images of the same shape, with
       optional internal sigma clipping-based masking or using externally
       supplied mask arrays. 
       
    The difference array, masks, and difference statistics can be 
    returned using accessors. An optional plot of the difference 
    histogram can be generated.
    """
    
    def __init__(self,
        imdata1, 
        imdata2, 
        sigmaclip,
        loglevel, 
        mask1=None, 
        mask2=None):
        """Take the difference of two images and calculate its statistical
           properties, after optionally masking the input images in one
           of two ways.
        
        :param imdata1:
        :param imdata2:
        :param sigmaclip:
        :param loglevel:
        :param mask1:
        :param mask2:
        """
    
        # Initialize logging
        self._loglevel = loglevel
        self._initialize_logger(self._loglevel)
        
        # Check input images and perform differencing. Note that we
        # must convert to a type than won't have underflow/overflow issues.
        self._check_input_images(imdata1, imdata2)
        self._diff_img = imdata1.astype(float) - imdata2.astype(float)
        
        # Generate the mask
        self._generate_mask(imdata1, imdata2, sigmaclip, mask1, mask2)
        return
        
    def _check_input_images(self, imdata1, imdata2):
        """Check that the input images match each other in shape
           and data type
        """
        
        if imdata1.shape != imdata2.shape:
            err_msg = ('Error, data array shapes do not match:'
                f' First file={imdata1.shape},'
                f' second file={imdata2.shape}')
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)
        
        if imdata1.dtype != imdata2.dtype:
            err_msg = ('Error, data types do not match:'
                f' First file={imdata1.dtype},'
                f' second file={imdata2.dtype}')
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)

        return
        
    def _combine_existing_bad_masks(self, mask1, mask2):
        """Combine one or more existing bad pixel masks into a good 
           pixel mask
           
        Note that the output good pixel mask is True where the pixels are
        good.
        
        :param mask1: None or array like with non-zero values denoting
          bad pixels.
        :param mask2: None or array like with non-zero values denoting
          bad pixels.
        """
        
        self._logger.debug('Using existing image masks.')
        if (mask1 is None) and (mask2 is None):
            err_msg = 'Error, no input masks were passed to _combine_existing_bad_masks'
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)
        elif (mask1 is not None) and (mask2 is not None):
            # Check the shapes are the same.
            if mask1.shape != mask2.shape:
                err_msg = f'Error, mask shapes do not match. mask1 is {mask1.shape}, while mask2 is {mask2.shape}'
                self._logger.error(err_msg)
                raise RunTimeError(err_msg)
            good1 = mask1 == 0
            good2 = mask2 == 0
        else:
            # Only one mask suuplied.
            # Generate an all-good mask to replace the missing mask.
            if mask1 is not None:
                good1 = mask1 == 0
                good2 = np.ones(mask1.shape, dtype=bool)
            else:                                   # mask2 is not None
                good2 = mask2 == 0
                good1 = np.ones(mask2.shape, dtype=bool)
            
        self._good_pixel_mask = np.logical_and(good1, good2)
        return

    def _generate_mask(self, data1, data2, sigmaclip, mask1, mask2):
        """Generate a mask to be used to select pixels when calculating
           the image difference statistics.
           
        There are two exclusive modes of masking, with the first being
        more appropriate where each image is largely uniform with a few
        outliers (e.g. a bias image or dark):
        - If sigmaclip is True then outliers in each input image (not
          the difference) are detected by sigma clipping. The final mask
          excludes outliers in either input image.
        - If sigmaclip is None or False and either of mask1 or mask2
          is not None, then non-zero pixels in the input mask are 
          assumed to denote bad pixels. The two masks are combined to
          generate a final mask.
        - If sigmaclip is None or False, and mask1 and mask2 are None,
          then all pixels in the difference image are assumed to be good.

        Note that:
        - The input masks should be non-zero where the data is
          to be excluded, i.e. zero is good data.
        - Specifying sigmaclip as True will override and thus ignore any
          input masks.
        
        Masking (excluding pixels from statistics) can be disabled by
        setting all inputs to None.
        
        :param data1: First input image.
        :param data2: Second input image.
        :param sigmaclip: None or boolean.
        :param mask1: None, or array like with non-zero values denoting
          a bad pixel. If None then the first input image is assumed
          to be all good. 
        :param mask2: None, or array like with non-zero values denoting
          a bad pixel. If None then the second input image is assumed
          to be all good. 
        """
        
        if sigmaclip is None:
            sigmaclip = False
            
        if sigmaclip is True:
            self._generate_sigmaclip_mask(data1, data2)
        else:
            if (mask1 is not None) or (mask2 is not None):
                # Check that the shapes match.
                ## TODO
                self._combine_existing_bad_masks(mask1, mask2)
            else:
                self._good_pixel_mask = np.ones(data1.shape, dtype=bool)
            
        # Report mask statistics. Note the mask we use is True where
        # pixels are good and False where they are bad.
        self._numpix  = self._good_pixel_mask.size
        self._numgood = np.sum(self._good_pixel_mask)
        numbad        = self._numpix - self._numgood
        self._logger.debug(f'Final good pixel mask has {self._numgood} good pixels out of {self._numpix} pixels ({numbad} bad).')
        
        return
        
    def _generate_sigmaclip_mask(self, data1, data2):
        """Creates a good pixel mask based on sigma-clipped statistics
           of each input image data array.
           
        This routine is most appropriate for images that expected to
        be relaitively uniform, but with a small number of highly 
        discrepant values.
        
        The generated mask is True when pixels are good and False where
        pixels are bad. If a pixel is bad in either input image it is
        considered bad in the combined mask.
           
        :param data1: First input image
        :param data2: Second input image
        """
        
        sigma             = 3.0
        self._logger.debug(f'Generating a good pixel mask using sigma={sigma} clipping on input image data values.')
        
        # Clip first input image
        mean1, med1, std1 = sigma_clipped_stats(data1, sigma=sigma)
        lo1               = med1 - (sigma * std1)
        hi1               = med1 + (sigma * std1)
        good1_lo          = data1 >= lo1
        good1_hi          = data1 <= hi1
        good1             = np.logical_and(good1_lo, good1_hi)
        
        # Clip second input image
        mean2, med2, std2 = sigma_clipped_stats(data2, sigma=sigma)
        lo2               = med2 - (sigma * std2)
        hi2               = med2 + (sigma * std2)
        good2_lo          = data2 >= lo2
        good2_hi          = data2 <= hi2
        good2             = np.logical_and(good2_lo, good2_hi)
        
        self._logger.debug(f'Good pixels in the first image have pixel values between {lo1:.2f} and {hi1:.2f} ADU.')
        self._logger.debug(f'Good pixels in the second image have pixel values between {lo2:.2f} and {hi2:.2f} ADU.')
        
        self._good_pixel_mask = np.logical_and(good1, good2)
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

    def data(self):
        """Return the image difference array
        
        Note that the returned array always has a floating point data type.
        """
        
        return self._diff_img

    def good_pixel_mask(self):
        """Return the good pixel mask
        
        Note that the returned boolean array is True where pixels are good.
        """
        
        return self._good_pixel_mask
        
    def stddev(self):
        """Return the standard deviation ofthe masked image difference.
        """
                
        stddev = np.std(self._diff_img[self._good_pixel_mask])
        return stddev
    
    def min(self):
        """Return the minimum value in the masked image difference.
        """
        
        minval = np.amin(self._diff_img[self._good_pixel_mask])
        return minval
        
    def max(self):
        """Return the maximum value in the masked image difference.
        """
        
        maxval = np.amax(self._diff_img[self._good_pixel_mask])
        return maxval

    def mean(self):
        """Return the mean value in the masked image difference.
        """
        
        meanval = np.mean(self._diff_img[self._good_pixel_mask])
        return meanval
        
    def median(self):
        """Return the median value in the masked image difference.
        """
        
        medval = np.median(self._diff_img[self._good_pixel_mask])
        return medval

    def numpix(self):
        """Returns the number of good pixels and the total number of 
           pixels.
        """
        return self._numgood, self._numpix
        
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
        biasfile1,      # Name/path of first bias file.
        biasfile2,      # Name/path of second bias file.
        gain,           # Gain value or gain FITS keyword to use.
        loglevel):
        """The constructor stores the names of the files and gain
           keyword or value, and initializes the logger, but does not
           perform any file IO or perform calculations.
           
           The actual file IO and calculations are only performed when
           the user calls estimate_rn.
        """
    
        # Initialize logging
        self._loglevel  = loglevel
        self._biasfile1 = biasfile1
        self._biasfile2 = biasfile2
        self._gaininfo  = gain
        self._initialize_logger(self._loglevel)
        self._logger.debug( ( 'Initialized an ApCalcReadNoise instance with'
            f' biasfile1={biasfile1}, biasfile1={biasfile2}, gain={gain},'
            f' and loglevel={loglevel}' ) )
        return
                
    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        self._logger = logging.getLogger('ApCalcReadNoise')
        
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
        
    def _check_file_exists(self, filename):
        if not os.path.isfile(filename):
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
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
        
        # Read the files files, assuming data is in primary.
        data1, hdr1 = self._read_fits(self._biasfile1, 0)
        data2, hdr2 = self._read_fits(self._biasfile2, 0)
        if data1.shape != data2.shape:
            err_msg = ('Error, data array shapes do not match:'
                f' First file={data1.shape},'
                f' second file={data2.shape}')
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)
        
        # Work out the gain value to use (e/ADU)
        self._gain = self._select_gain(hdr1, hdr2)
        self._logger.info(f'Adopted gain is {self._gain:.2f} electrons/ADU.')
        
        # Calculate the difference with or without sigma clipping and
        # histogram plotting, get the standard deviation (ADU)
        im_diff               = ApImageDifference(data1, 
            data2, 
            sigmaclip, 
            self._loglevel,
            mask1=None, 
            mask2=None)
        stddev                = im_diff.stddev()
        mindiff               = im_diff.min()
        maxdiff               = im_diff.max()
        npix_good, npix_total = im_diff.numpix()
        pct_bad               = 100 * (npix_total - npix_good)/npix_total
        self._logger.info(f'Standard deviation={stddev:.2f} ADU using {npix_good}/{npix_total} pixels ({pct_bad:.3f} % bad).')
        
        # Generate plot
        if histplot is not None:
            self._plot_difference_histogram(im_diff, histplot)
        
        # Calculate read noise estimate in e/pixel
        read_noise = self._gain * stddev / math.sqrt(2)
        self._logger.info(f'Estimated resd noise is {read_noise:.2f} e/ADU')
        
        return read_noise
        
    def _isfloat(self, value_str):
        """Returns True if value_str can be converted to a numeric float
           value.
        
        Will also return True for booleans and NaNs.
        """
        
        try:
            float(value_str)
            return True
        except ValueError:
            return False

    def _plot_difference_histogram(self, im_diff, histplot):
        """Generate a figure showing a histogram of the image differences.
        """
        
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
        ax1.tick_params(axis='both', labelsize=6)
        ax2.tick_params(axis='both', labelsize=6)
        title_fsize  = 7
        axis_fsize   = 6
        legend_fsize = 5
        title_str = ('Histogram of pixel-to-pixel differences:\n' 
            f'  First image: {self._biasfile1}\n' 
            f'  Second image: {self._biasfile2}')
        
        # Get min and max values to use for histogram.
        minval  = im_diff.min()
        maxval  = im_diff.max()
        stddev  = im_diff.stddev()
        meanval = im_diff.mean()
        medval  = im_diff.median()
        bins    = np.arange(minval, maxval+1) - 0.5
        
        # Get the data and apply the mask.
        rawdata  = im_diff.data()
        rawmask  = im_diff.good_pixel_mask()
        gooddata = rawdata[rawmask].copy()
        
        for ax in [ax1, ax2]:
            plt.sca(ax)
            nperbin, outbins, patches = plt.hist(gooddata.ravel(), 
                bins=bins, 
                density=False, 
                weights=None,
                label='Masked difference',
                alpha=0.6)
            nperbin, outbins, patches = plt.hist(rawdata.ravel(), 
                bins=bins, 
                density=False, 
                weights=None,
                label='Raw image difference',
                alpha=0.3)

            if ax == ax2:
                ax.set_yscale('log')
                
            # Get plotted axis limits for plot.
            ymin, ymax = ax.get_ylim()
            ax.vlines(meanval, ymin, ymax, 
                linestyles='dashed', linewidth=0.8, label='Masked data mean')
            ax.vlines(medval, ymin, ymax, 
                linestyles='dotted', linewidth=0.8, label='Masked data median')
                
            ax.set_xlabel('Pixel difference (ADU)', fontsize=axis_fsize)
            ax.set_ylabel('Number of pixels', fontsize=axis_fsize)
            ax.legend(fontsize=legend_fsize, loc='center right')
        
        
        fig.suptitle(title_str, fontsize=title_fsize)
        fig.savefig(histplot, dpi=200, bbox_inches='tight')
        self._logger.info(f'Wrote histogram plot to {histplot}')
        
        return

    def _select_gain(self, hdr1, hdr2):
        """Determines what gain value to use, based on the two file headers
           and the gaininfo supplied at object construction.
           
        At construction a gain value (string representation of a floating
        point number) or the name of a FITS primary header keyword was
        supplied, and is stored as _gaininfo.
        - If a string representation of a floating point number was 
          supplied then we should use that.
        - If _gaininfo is not convertable to float, treat it as a FITS
          keyword.
          - If the keyword is NOT present in both FITS headers then this
            is treated as an error.
          - If the value associated with the keyword is different beyond
            a numerical tolerance this is an error.
        
        If no gain information can be found then this function raises
        an exception.
        """
        
        if self._isfloat(self._gaininfo):
            # Use the value given instead of looking in the FITS files.
            return float(self._gaininfo)
        
        # Otherwise treat _gaininfo as a FITS header keyword
        gain1 = None
        if self._gaininfo in hdr1:
            gain1 = float( hdr1[self._gaininfo] )
        gain2 = None
        if self._gaininfo in hdr2:
            gain2 = float( hdr2[self._gaininfo] )
        
        # If we can't find either or both this is a fatal error.
        if (gain1 is None) or (gain2 is None):
            err_msg = f'Error, {self._gaininfo} gain keyword not found in'
            if (gain1 is None) and (gain2 is None):
                err_msg += ' both FITS files.'
            elif gain1 is None:
                err_msg += ' the first FITS file.'
            else:                                       # gain2 is None
                err_msg += ' the second FITS file.'
            
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)

        # Now check the gains are the same to within a tolerance.
        # Gains can differ by this much and be OK (e/ADU)
        tolerance = 0.001
        if math.fabs(gain1-gain2) > tolerance:
            err_msg = f'Error, gains differ by more than {tolerance:.3f} e/ADU, where gain1={gain1:.3f}, gain2={gain2:.3f}.'
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)
        
        return gain1
        
                
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
    print(f'Estimated read noise is {rn1:.2f} electrons/pixel.')
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
