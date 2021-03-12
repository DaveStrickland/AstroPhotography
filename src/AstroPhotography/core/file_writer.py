"""
Function based interface to FileWriter
"""

from .logger import logger
import imageio
import time
import os.path
from astropy.io import fits

def file_writer(out_file, data_array, exif_dict):
    """
    Abstracts away details of writing an image to a graphics
    file format or astronomical FITS file.

    :param out_file: File path/name for output file.       
    :param data_array: 2-D (monochrome) or 3-D (rgb) ndarray of
      data values. Usually but not always of dtype.uint16
    :param exif_dict: Dictionary of EXIF tags and values in the format
      produced by ExifRead. This will be processed to provide a format
      suitable for the exif package or astropy FITS to write.
    """
    
    ftype = determine_file_type(out_file)
    
    # Image properties.
    ndim  = data_array.ndim
    nrows = data_array.shape[0]
    ncols = data_array.shape[1]
    nlayers = 1
    if ndim == 3:
        nlayers = data_array.shape[2]
    elif ndim == 1 or ndim > 3:
        err_str = 'Error: file_writer cannot handle {}-dimension arrays.'.format(ndim)
        logger.error(err_str)
        return # maybe throw?
        
    dtype = data_array.dtype
    
    if ndim == 2:
        info_str = 'Writing {}x{} {} image to {:s}'.format(nrows,
            ncols,
            dtype,
            out_file)
    elif ndim == 3:
        info_str = 'Writing {}x{}x{} {} image to {:s}'.format(nrows,
            ncols,
            ndim,
            dtype,
            out_file)
    logger.info(info_str)
        
    # TODO make sure write is succesful.
    t_start = time.perf_counter()
    if ftype == 'graphics':
        # Special handling of format?
        our_format=None
        if 'jp2' in out_file:
            our_format='JPEG2000-PIL'
        imageio.imsave(uri=out_file, 
            im=data_array,
            format=our_format)
    elif ftype == 'fits':
        hdu_list = fits.PrimaryHDU(data_array)
        hdu_list.writeto(out_file, 
            output_verify='warn',
            overwrite=True)
    else:
        raise RuntimeError('Could not determine type for output file {}'.format(out_file))
    t_end = time.perf_counter()
    
    io_size = os.path.getsize(out_file) / 1e6 # size in MB
    io_time = t_end - t_start                 # seconds
    io_rate = io_size / io_time               # MB/s
    logger.debug('Wrote {:.3f} MB file successfully in {:.3f} seconds ({:.3f} MB/s)'.format(io_size,
        io_time,
        io_rate))
        
    # This should be moved to a function with error/exception checking.
    if ftype == 'graphics':
        logger.warning('Writing EXIF metadata not yet supported.')
    return
        
def CapitalCase_to_snake_case(input_str):
    """
    Converts CapitalCase or camelCase strings to snake_case
    
    Based on: https://www.geeksforgeeks.org/python-program-to-convert-camel-case-string-to-snake-case/
    
    Limitations:
    - Messes up consecutive capital letters, e.g. ISO become i_s_o, or
      FocalLengthMM becomes focal_length_m_m
    
    :param input_str: Input CapitalCase or camelCase string
    """
    
    output_str = ''.join(['_'+i.lower() if i.isupper() 
               else i for i in input_str]).lstrip('_')
               
    return output_str
    
def determine_file_type(fname):
    """
    Determine if file is a common graphics format or FITS.
    """
    
    # Graphics file formats, delegated to imageio
    graphics_list = ['.tif', '.tiff', '.jpg', '.jp2', 'jpeg', '.png', '.gif']
    
    # FITS, delegate to astropy
    fits_list = ['.fits', '.ftz', '.fit', '.fits.gz'] 
    
    # want file extension, e.g .png, .fits but also .fits.gz
    root, ext = os.path.splitext(fname)
    if '.gz' in ext:
        root, ext2 = os.path.splitext(root)
        ext = ext2 + ext
    ext = ext.lower()
        
    if ext in fits_list:
        ftype = 'fits'
    elif ext in graphics_list:
        ftype = 'graphics'
    else:    
        ftype = 'unknown'
    
    return ftype
