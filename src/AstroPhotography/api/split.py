""" Implements the split command.

"""
from ..core.logger import logger
from ..core.RawConv import *
import imageio

def main(rawfile, output):
    """ Execute the split command.
    
    :param rawfile: RAW input file to process
    :param output: File prefix for output images. Output image file names
      will consist of the prefix followed by '_r.png', '_g1.png', '_b.png'
      and '_g2.png'
    """
    logger.info("Executing split command on {:s}".format(rawfile))
    
    # Read and process the file.
    rawconv = RawConv(rawfile)
    r_im, g1_im, b_im, g2_im = rawconv.split()
    
    r_fname  = output + '_r.png'
    g1_fname = output + '_g1.png'
    b_fname  = output + '_b.png'
    g2_fname = output + '_g2.png'
    
    logger.debug("Wrote R image to {:s}".format(r_fname))
    imageio.imsave(r_fname, r_im)
    
    logger.debug("Wrote G1 image to {:s}".format(g1_fname))
    imageio.imsave(g1_fname, g1_im)
    
    logger.debug("Wrote B image to {:s}".format(b_fname))
    imageio.imsave(b_fname, b_im)
    
    logger.debug("Wrote G2 image to {:s}".format(g1_fname))
    imageio.imsave(g2_fname, g2_im)

    return 
