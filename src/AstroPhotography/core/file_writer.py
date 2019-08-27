"""Function based interface to FileWriter
"""

from .logger import logger
#import .FileWriter
import imageio
import time
import os.path

def file_writer(data_array, out_file):
    """
    """
    nrows = data_array.shape[0]
    ncols = data_array.shape[1]
    dtype = data_array.dtype
    
    logger.info('Writing {}x{} {} image to {:s}'.format(nrows,
        ncols,
        dtype,
        out_file))
        
    # TODO make sure write is succesful.
    t_start = time.perf_counter()
    imageio.imsave(out_file, data_array)
    t_end = time.perf_counter()
    
    io_size = os.path.getsize(out_file) / 1e6 # size in MB
    io_time = t_end - t_start                 # seconds
    io_rate = io_size / io_time               # MB/s
    logger.debug('Wrote {:.3f} MB file successfully in {:.3f} seconds ({:.3f} MB/s)'.format(io_size,
        io_time,
        io_rate))
    return
    
