""" Implements the split command.

"""
from ..core.logger import logger
from ..core.RawConv import RawConv
from ..core.file_writer import file_writer

def main(rawfile, output, black, extension):
    """ Execute the split command.
    
    :param rawfile: RAW input file to process.
    :param output: File prefix for output images. Output image 
      file names will consist of the prefix followed by '_r.', 
      '_g1.', '_b.' and '_g2.', followed by the file name
      extension.
    :param black: If true the camera black levels will be subtracted
      from the channel data.
    :param extension: File type and extension to use, e.g. 'png', 'jpg'
      'jp2', or 'fits'
    """
    logger.info("Executing split command on {:s}".format(rawfile))
    
    # Read and process the file.
    rawconv = RawConv(rawfile)
    r_im, g1_im, b_im, g2_im = rawconv.split(subtract_black=black)
    
    r_fname  = output + '_r.'  + extension
    g1_fname = output + '_g1.' + extension
    b_fname  = output + '_b.'  + extension
    g2_fname = output + '_g2.' + extension
    
    file_writer(r_im,  r_fname)
    file_writer(g1_im, g1_fname)
    file_writer(b_im,  b_fname)
    file_writer(g2_im, g2_fname)
    return 
