""" Implements the RawConv class
"""

from .logger import logger
import rawpy
import numpy as np

class RawConv:
    """ Higher level conversions and processing of RAW image files to graphics formats 
        and/or FITS image format.
    """
    
    def __init__(self, rawfile):
        self._supported_colors = ['RGBG']
        self._load(rawfile)
        return
        
    def __del__(self):
        """ RawConv class destructor
        """
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
        print('RawConv on {}'.format(rawfile))
        self._rawfile      = rawfile
        self._rawpy        = rawpy.imread(rawfile)
        
        # Common image characteristics
        self._nrows        = self._rawpy.raw_image_visible.shape[0]
        self._ncols        = self._rawpy.raw_image_visible.shape[1]
        self._color_desc   = str(self._rawpy.color_desc)
        if self._color_desc not in self._supported_colors:
            self._unsupported_colors
        # Know that for RGBG the color indices are 0, 1, 2, 3
        self.R  = 0
        self.G1 = 1
        self.B  = 2
        self.G2 = 3
        self._color_map    = self._rawpy.raw_colors_visible.copy()
        self._mask_r       = self._color_map == self.R
        self._mask_g1      = self._color_map == self.G1
        self._mask_b       = self._color_map == self.B
        self._mask_g2      = self._color_map == self.G2

        self._black_levels = self._rawpy.black_level_per_channel
        self.default_whitebalances()
        return

    def default_whitebalances(self):
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

    def split(self, subtract_black=False, verbose=False):
        """ Exports the raw, unprocessed, bayer RGBG as four separate 
            uint16 numpy arrays.
        
        For each band (R, G etc) the only non-zero pixels will be that 
        band, and pixels associated with other bands will be set to zero.
        
        :param subtract_black: If true the camera black levels will be 
          subtracted from the channel data.
        :param verbose: TBA 
        """
        
        rawim      = self._rawpy
        print('black_levels', self._black_levels)
        print('camera_wb', self._wb_camera, type(self._wb_camera))
        print('daylight_wb', self._wb_daylight, type(self._wb_daylight))
        
        tmp_raw = rawim.raw_image_visible.copy()
        
        if subtract_black:
            r_bg  = int( self._black_levels[self.R]  )
            g1_bg = int( self._black_levels[self.G1] )
            b_bg  = int( self._black_levels[self.B]  )
            g2_bg = int( self._black_levels[self.G2] )
            print('r_bg={} g1_bg={} b_bg={} g2_bg={}'.format(r_bg,
                g1_bg, b_bg, g2_bg), type(r_bg))
        else:
            r_bg = g1_bg = b_bg = g2_bg = int(0)

        # Note python indexing is inclusive to exclusive, so the 2 is excluded.
        r_im = np.zeros(rawim.raw_image_visible.shape[0:2], 
            dtype=np.uint16)
        # np.subtract will leave original out array values if it 
        # supplied and not None, and where the where mask is False
        np.subtract(tmp_raw, r_bg, 
            out=r_im, where=self._mask_r)
        
        g1_im = np.zeros(rawim.raw_image_visible.shape[0:2], 
            dtype=np.uint16)
        np.subtract(tmp_raw, g1_bg, 
            out=g1_im, where=self._mask_g1)
            
        b_im = np.zeros(rawim.raw_image_visible.shape[0:2], 
            dtype=np.uint16)
        np.subtract(tmp_raw, b_bg, 
            out=b_im, where=self._mask_b)
                
        g2_im = np.zeros(rawim.raw_image_visible.shape[0:2], 
            dtype=np.uint16)
        np.subtract(tmp_raw, g2_bg, 
            out=g2_im, where=self._mask_g2)
            
        # where is good for selecting one or the other, not for operations.
        # r_im  = np.where(self._color_map == R,  tmp_raw_r,  0)
        # g1_im = np.where(self._color_map == G1, tmp_raw_g1, 0)
        # b_im  = np.where(self._color_map == B,  tmp_raw_b,  0)
        # g2_im = np.where(self._color_map == G2, tmp_raw_g2, 0)
        
        print('raw_r=', r_im, r_im.dtype)
        print('raw_g1=', g1_im, g1_im.dtype)
        
        #r_im3d  = rawim.postprocess(output_color=rawpy.ColorSpace.raw,
        #    gamma=(1, 1), 
        #    no_auto_bright=True, 
        #    no_auto_scale=True,
        #    dcb_enhance=False,
        #    user_wb=[1, 0, 0, 0], 
        #    output_bps=16, 
        #    user_flip=0)
            
        # Note python indexing is inclusive to exclusive, so the 2 is excluded.
        #r_im = np.zeros(r_im3d.shape[0:2], dtype=np.uint16)
            
        #print(r_im3d.shape, r_im3d.dtype)
        print(r_im.shape, r_im.dtype)
        print('R  min, max: ', np.nanmin(r_im),  np.nanmax(r_im))
        print('G1 min, max: ', np.nanmin(g1_im), np.nanmax(g1_im))
        print('B  min, max: ', np.nanmin(b_im),  np.nanmax(b_im))
        print('G2 min, max: ', np.nanmin(g2_im), np.nanmax(g2_im))
        return r_im, g1_im, b_im, g2_im
