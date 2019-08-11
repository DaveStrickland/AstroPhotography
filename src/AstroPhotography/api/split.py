""" Implements the split command.

"""
from ..core.logger import logger
from ..core.RawConv import *
import imageio

def main(rawfile, output):
    """ Execute the split command.
    
    :param name: name to use in greeting
    """
    logger.debug("Executing split command on {:s}".format(rawfile))
    
    # Read and process the file.
    rawconv = RawConv(rawfile)
    r_im, g1_im, b_im, g2_im = rawconv.split()
    
    rawname = get_rawname(rawfile, output)
    print('Output rawname is {}'.format(rawname))
    r_fname  = rawname + '_r.png'
    g1_fname = rawname + '_g1.png'
    b_fname  = rawname + '_b.png'
    g2_fname = rawname + '_g2.png'
    imageio.imsave(r_fname, r_im)
    imageio.imsave(g1_fname, g1_im)
    imageio.imsave(b_fname, b_im)
    imageio.imsave(g2_fname, g2_im)

    return 

def get_rawname(rawfile, output):
    """ Return the root name for the output R, G, B, G frames based on
        either the name of the RAW file itself or a user-supplied root name.
    """
    rawname = 'bob'
    return rawname
