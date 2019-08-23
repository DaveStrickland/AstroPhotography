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
        print('RawConv on {}'.format(rawfile))
        self._rawfile = rawfile
        self._rawpy = rawpy.imread(rawfile)
        return
        
    def __del__(self):
        """ RawConv class destructor
        """
        self._rawpy.close()
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
        color_desc = str(rawim.color_desc)
        black_levels = rawim.black_level_per_channel
        print('black_levels', black_levels)
        wb_list    = rawim.camera_whitebalance
        print('white_levels', wb_list)
        
        if 'RGBG' not in color_desc:
            raise NotiImplementedError('RawConv.split() not implemented for {} images.'.format(color_desc))
        
        num_rows = rawim.raw_image_visible.shape[0]
        num_cols = rawim.raw_image_visible.shape[1]
        
        # Know that for RGBG the color indices are 0, 1, 2, 3
        tmp_raw = rawim.raw_image_visible.copy()
        tmp_col = rawim.raw_colors_visible.copy()
        tmp_pat = rawim.raw_pattern.copy()
        
        #r_im = np.zeros([num_rows, num_cols], dtype=np.uint16)
        r_im  = np.where(tmp_col == 0, tmp_raw, 0)
        g1_im = np.where(tmp_col == 1, tmp_raw, 0)
        b_im  = np.where(tmp_col == 2, tmp_raw, 0)
        g2_im = np.where(tmp_col == 3, tmp_raw, 0)
        
        print('raw=', tmp_raw, tmp_raw.dtype)
        print('colors=', tmp_col, tmp_col.dtype)
        print('pattern=', tmp_pat, tmp_pat.dtype)
        
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
        print(np.nanmin(r_im), np.nanmax(r_im))
        print(np.nanmin(g1_im), np.nanmax(g1_im))
        print(np.nanmin(b_im), np.nanmax(b_im))
        print(np.nanmin(g2_im), np.nanmax(g2_im))
        return r_im, g1_im, b_im, g2_im
