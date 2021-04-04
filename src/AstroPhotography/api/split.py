"""
Implements the split command calling RawConv.
"""
import time
from ..core.logger import logger
from ..core.RawConv import RawConv
from ..core.file_writer import file_writer

def main(rawfile, output, keepblack, extension):
    """ Execute the split command.
    
    :param rawfile: RAW input file to process.
    :param output: File prefix for output images. Output image 
      file names will consist of the prefix followed by '_r.', 
      '_g1.', '_b.' and '_g2.', followed by the file name
      extension.
    :param keepblack: If true the camera black levels (similar to a CCD bias?)
    :param extension: File type and extension to use, e.g. 'png', 'jpg'
      'jp2', or 'fits'
    """
    logger.info(f'Executing split command on {rawfile:s}')
    t_start = time.perf_counter()
    subblack = not keepblack
    
    # Read and process the file.
    rawconv = RawConv(rawfile)
    r_im, g1_im, b_im, g2_im, exif_dict = rawconv.split(subtract_black=subblack)
    
    r_fname  = output + '_r.'  + extension
    g1_fname = output + '_g1.' + extension
    b_fname  = output + '_b.'  + extension
    g2_fname = output + '_g2.' + extension
    
    file_writer(r_fname,  r_im,  exif_dict)
    file_writer(g1_fname, g1_im, exif_dict)
    file_writer(b_fname,  b_im,  exif_dict)
    file_writer(g2_fname, g2_im, exif_dict)
    
    t_end = time.perf_counter()
    t_run = t_end - t_start
    logger.info(f'Completed split command in {t_run:.2f} seconds.')
    return 
