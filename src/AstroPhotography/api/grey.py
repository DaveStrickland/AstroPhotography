""" Implements the grey command.

"""
from ..core.logger import logger
from ..core.RawConv import *
from ..core.file_writer import file_writer

def main(rawfile, output, method, black, whitebalance):
    """ Execute the grey command.
    
    :param rawfile: RAW input file to process.
    :param output: Name of the output file to write the luminance image to. This
      may be a graphics file format or the astronomical FITS format. The correct
      writer will be used depending on the user-supplied file extension.
    :param method: Luminance method used to construct monochrome image from 
      Bayer channels. See RawConv.grey() for details on allowable methods.
    :param whitebalance: Whitebalance method to use when converting R, G and B
      channels to monochrome.
    :param black: If true the camera black levels will be subtracted
      from the channel data.
    """
    logger.info("Executing grey command on {:s}".format(rawfile))
    
    # Read and process the file.
    rawconv = RawConv(rawfile)
    grey_im = rawconv.grey(luminance_method=method, 
        subtract_black=black, 
        wb_method=whitebalance)
    file_writer(grey_im, output)
    return 
