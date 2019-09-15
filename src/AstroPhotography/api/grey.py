""" Implements the grey command.

"""
from ..core.logger import logger
from ..core.RawConv import *
import imageio

def main(rawfile, output, black, whitebalance):
    """ Execute the grey command.
    
    :param rawfile: RAW input file to process.
    :param output: Name of the output file to write the luminance image to. This
      may be a graphics file format or the astronomical FITS format. The correct
      writer will be used depending on the user-supplied file extension.
    :param whitebalance: Whitebalance method to use when converting R, G and B
      channels to monochrome
    :param black: If true the camera black levels will be subtracted
      from the channel data.
    """
    logger.info("Executing grey command on {:s}".format(rawfile))
    
    # Read and process the file.
    rawconv = RawConv(rawfile)
    wb_vals = rawconv.get_whitebalance(whitebalance)
    grey_im = rawconv.grey(subtract_black=black, wb_list=wb_vals)
    file_writer(grey_im, output)
    return 
