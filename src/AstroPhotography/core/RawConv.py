""" Implements the RawConv, which encapsulates raw data manipulations.
"""

from .logger import logger
import rawpy
import numpy as np
import os.path

class RawConv:
    """ Higher level conversions and processing of RAW image files to graphics formats 
        and/or FITS image format.
    """
    
    def __init__(self, rawfile):
        """RawConv constructor loads the raw file and performs initial processing.
        """
        
        self._supported_colors = ['RGBG']
        self._load(rawfile)
        return
        
    def __del__(self):
        """ RawConv class destructor
        """
        
        if self._rawpy is not None:
            self._rawpy.close()
        return

    def _unsupported_colors(self):
        """ Raises an exception if the image color description is not 
            supported.
        """
        
        err_str = 'Image colors ({}) not in supported list of {}.'.format(self._color_desc, 
            self._supported_colors)
        raise NotImplementedError(err_str)        
        return

    def _load(self, rawfile):
        self._rawfile    = rawfile
        self._rawpy      = None
        if os.path.isfile(rawfile):
            logger.debug('RawConv loading {}'.format(rawfile))
            self._rawpy      = rawpy.imread(rawfile)
        else:
            err_msg = 'RawConv cannot load {}. Not a valid path or file.'.format(rawfile)
            logger.error(err_msg)
            raise RuntimeError(err_msg)
        
        # Common image characteristics
        self._nrows        = self._rawpy.raw_image_visible.shape[0]
        self._ncols        = self._rawpy.raw_image_visible.shape[1]
        self._color_desc = str(self._rawpy.color_desc)
        if self._color_desc not in self._supported_colors:
            self._unsupported_colors
        # Know that for RGBG the color indices are 0, 1, 2, 3
        self.R  = 0
        self.G1 = 1
        self.B  = 2
        self.G2 = 3
        self._color_map  = self._rawpy.raw_colors_visible.copy()
        self._mask_r     = self._color_map == self.R
        self._mask_g1    = self._color_map == self.G1
        self._mask_b     = self._color_map == self.B
        self._mask_g2    = self._color_map == self.G2
        
        self._rawim_r    = np.where(self._mask_r,  self._rawpy.raw_image_visible, 0)
        self._rawim_g1   = np.where(self._mask_g1, self._rawpy.raw_image_visible, 0)
        self._rawim_b    = np.where(self._mask_b,  self._rawpy.raw_image_visible, 0)
        self._rawim_g2   = np.where(self._mask_g2, self._rawpy.raw_image_visible, 0)

        self._black_levels = self._rawpy.black_level_per_channel
        self._black_subtracted = False
        self._default_whitebalances()
        return

    def _default_whitebalances(self):
        """Convert default camera and daylight whitebalances into
           a more consistent form
           
        For the Canon Digital Rebel XTi, the format of the camera
        whitebalance and the daylight whitebalance format are
        inconsistent:
        - camera:   1997.0, 1080.0, 2333.0, 1080.0
        - daylight: 2.4238, 0.9213, 1.1510, 0.0

        This function converts them into a uniform set of R G B G
        factors where the lowest value is 1.0.
        """
        self._wb_camera   = []
        self._wb_daylight = []

        tmp_cam = self._rawpy.camera_whitebalance
        tmp_min = min(tmp_cam)
        for idx, val in enumerate(tmp_cam):
            self._wb_camera.append( val/tmp_min )
        
        # If last item is zero and corresponds to G, then use other G.
        tmp_day = self._rawpy.daylight_whitebalance
        if tmp_day[3] == 0 and self._color_desc[3] == 'G':
            tmp_day[3] = tmp_day[1]
        tmp_min = min(tmp_day)
        for idx, val in enumerate(tmp_day):
            self._wb_daylight.append( val/tmp_min )
        
        return

    def _subtract_black_levels(self):
        """Subtract black levels from sub-bands.
        """
        if self._black_subtracted:
            # Already black subtracted
            return
            
        r_bg  = int( self._black_levels[self.R]  )
        g1_bg = int( self._black_levels[self.G1] )
        b_bg  = int( self._black_levels[self.B]  )
        g2_bg = int( self._black_levels[self.G2] )
        print('r_bg={} g1_bg={} b_bg={} g2_bg={}'.format(r_bg,
            g1_bg, b_bg, g2_bg), type(r_bg))

        self._rawim_r  = self._safe_subtract(self._rawim_r,
            self._mask_r, 
            r_bg)
        self._rawim_g1 = self._safe_subtract(self._rawim_g1,
            self._mask_g1, 
            g1_bg)
        self._rawim_b  = self._safe_subtract(self._rawim_b,
            self._mask_b, 
            b_bg)
        self._rawim_g2 = self._safe_subtract(self._rawim_g2,
            self._mask_g2, 
            g2_bg)
        self._black_subtracted = True

    def _safe_subtract(self, data_arr, data_mask, val_to_subtract):
        """Safe subtraction of one array or scalar from another array
           avoiding unsigned integer wrap-around.
        
        It is possible for a value within a data array to be less
        than the value to subtract, so we need to carefully reset those
        pixels to the black level to avoid them wrapping around to
        a very large number.
        
        :param data_array: 2-D uint array to subract values from. Only
          those pixels within data_mask will be modified. 
        :param data_mask:  2-D boolean mask specifying where data_array
          is valid.
        :param val_to_subtract: 2-D uint array of values to subtract from
          data_array, or a scalar that will be broadcast to match the shape
          of data_array.
        :returns: 2-D uint array.
        """
        odd_mask = data_arr < val_to_subtract
        num_odd  = np.sum(odd_mask)
        pct_odd  = 100 * num_odd / data_arr.size
        logger.warn('Input array has {} pixels ({:.2f}%) less than array we will subract from it.'.format(num_odd, pct_odd))
        if num_odd > 0:
            # Reset those pixels to the value we're going to subtract.
            tmp_raw = np.where(odd_mask, val_to_subtract, data_arr)
        else:
            # Just use normal rawim
            tmp_raw = data_arr.copy()
        
        # np.subtract will leave original output array values (if out 
        # supplied and not None) where the where mask is False.
        # Elsewhere, the pixels will be the first input - second input.
        return np.subtract(tmp_raw, 
            val_to_subtract, 
            where=data_mask).copy()

    def get_whitebalance(self, wb_method):
        """Return the whitebalance multiplies for a given white-balance
           calculation method.
        """
        
        allowed_methods = ['daylight', 'camera']
        if wb_method not in allowed_methods:
            err_str = 'Unexpected white balance method {} not ones of the allowed method: {}'.format(wb_method, allowed_methods)
            logger.error(err_str)
            raise RuntimeError(err_str)
        
        wb_list = [1, 1, 1, 1]
        if wb_method == 'daylight':   
            wb_list = self._wb_daylight
        elif wb_method == 'camera':
            wb_list = self._wb_camera
            
        logger.debug('White balance values adopted: {}'.format(wb_list))
        return wb_list

    def grey(self, luminance_method='direct', subtract_black=False, wb_list=None, 
        verbose=False):
        """Create a luminance image using the supplied white-balance
           values, with or without black-level subtraction.
        
        For each band (R, G etc) the only non-zero pixels will be that 
        band, and pixels associated with other bands will be set to zero.
        
        :param luminance_method: Luminance method to use. At present only the
           following method(s) are supported:
           - direct: Each channel is multiplied by its white-balance factor
             and the added to the final image. There is no interpolation.
        :param subtract_black: If true the camera black levels will be 
          subtracted from the channel data.
        :param wb_list: List of white-balance multipliers to apply to
          each channel.
        :param verbose: TBA
        """
        
        allowed_methods=['direct']
        if luminance_method not in allowed_methods:
            err_str = 'Unexpected luminance calculate method supplied to RawConv.grey: {}. Allowed methods are: {}'.format(luminance_method, allowed_methods)
            logger.error(err_str)
        
        # Should check that wb_list is the correct size.
        
        
        if subtract_black:
            self._subtract_black_levels()
           
        r_im  = np.where(self._mask_r,  self._rawim_r,  0)
        g1_im = np.where(self._mask_g1, self._rawim_g1, 0)
        b_im  = np.where(self._mask_b,  self._rawim_b,  0)
        g2_im = np.where(self._mask_g2, self._rawim_g2, 0)
        
        # Perform calculations as double
        grey_im = np.zeros(r_im.shape[0:2], dtype=np.float64)
        grey_im += wb_list[self.R]  * r_im
        grey_im += wb_list[self.G1] * g1_im
        grey_im += wb_list[self.B]  * b_im
        grey_im += wb_list[self.G2] * g2_im
        
        # Should really do some sanity calculations here.
        
        return grey_im.astype(np.uint16)

    def split(self, subtract_black=False, verbose=False):
        """Exports the raw, unprocessed, bayer RGBG as four separate 
        uint16 numpy arrays.
        
        For each band (R, G etc) the only non-zero pixels will be that 
        band, and pixels associated with other bands will be set to zero.
        
        :param subtract_black: If true the camera black levels will be 
          subtracted from the channel data.
        :param verbose: TBA 
        """
        
        print('black_levels', self._black_levels)
        print('camera_wb', self._wb_camera, type(self._wb_camera))
        print('daylight_wb', self._wb_daylight, type(self._wb_daylight))
        
        if subtract_black:
            self._subtract_black_levels()
           
        r_im  = np.where(self._mask_r,  self._rawim_r,  0)
        g1_im = np.where(self._mask_g1, self._rawim_g1, 0)
        b_im  = np.where(self._mask_b,  self._rawim_b,  0)
        g2_im = np.where(self._mask_g2, self._rawim_g2, 0)
        
        # --> Should be moved to custom whitebalance
        # minr=400
        # maxr=500
        # minc=2800
        # maxc=2900
        # print('raw_r=', r_im[minr:maxr,minc:maxc], r_im.dtype)
        # print('raw_g1=', g1_im[minr:maxr,minc:maxc], g1_im.dtype)
        
        # --> Should be moved to RGB
        #r_im3d  = rawim.postprocess(output_color=rawpy.ColorSpace.raw,
        #    gamma=(1, 1), 
        #    no_auto_bright=True, 
        #    no_auto_scale=True,
        #    dcb_enhance=False,
        #    user_wb=[1, 0, 0, 0], 
        #    output_bps=16, 
        #    user_flip=0)
            
        # --> Should be moved to stats...
        # print('R  min, max: ', np.nanmin(r_im[self._mask_r]),   np.nanmax(r_im[self._mask_r]))
        # print('G1 min, max: ', np.nanmin(g1_im[self._mask_g1]), np.nanmax(g1_im[self._mask_g1]))
        # print('B  min, max: ', np.nanmin(b_im[self._mask_b]),   np.nanmax(b_im[self._mask_b]))
        # print('G2 min, max: ', np.nanmin(g2_im[self._mask_g2]), np.nanmax(g2_im[self._mask_g2]))
        return r_im, g1_im, b_im, g2_im
