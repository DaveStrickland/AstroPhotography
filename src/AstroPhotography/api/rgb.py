"""
Implements the rgb command using RawConv
"""
import time
from ..core.logger import logger
from ..core.RawConv import *
from ..core.file_writer import file_writer

def main(rawfile, output, method, keepblack, whitebalance, renormalize):
    """ Execute the rgb command.
    
    :param rawfile: RAW input file to process.
    :param output: Name of the output file to write the 3-channel RGB image to. This
      may be a graphics file format or the astronomical FITS format. The correct
      writer will be used depending on the user-supplied file extension.
    :param method: Method used to construct the 3-channel image from 
      Bayer channels. See RawConv.rgb() for details on allowable methods.
    :param whitebalance: Whitebalance method to use when converting R, G and B
      channels.
    :param keepblack: If true the camera black levels (similar to a CCD bias?)
      will NOT be subtracted from the channel data.
    :param renormalize: If true the output image will be linearly renormalized
      to fill the full dynamic range of the output image. This is useful
      for quick inspection of images, but not recommended if the images are
      to be used for image combination or amateur scientific purposes. 
    """
    logger.info(f'Executing rgb command on {rawfile:s}')
    t_start = time.perf_counter()
    
    # Read and process the file.
    rawconv = RawConv(rawfile)
    
    subblack = not keepblack
    stats    = True
    
    rgb_im, exif_dict = rawconv.rgb(luminance_method=method, 
        subtract_black=subblack, 
        wb_method=whitebalance,
        print_stats=stats,
        renorm=renormalize)
    file_writer(output, rgb_im, exif_dict)
    
    t_end = time.perf_counter()
    t_run = t_end - t_start
    logger.info(f'Completed rgb command in {t_run:.2f} seconds.')
    return 
