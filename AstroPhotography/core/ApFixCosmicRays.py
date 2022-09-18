"""Contains the implementation of the ApFixCosmicRays class.
"""

# 2021-01-09 dks : Initial start on coding.

import logging
from pathlib import Path
import math
from datetime import datetime, timezone
import time         # For performance counter
 
import numpy as np
import ccdproc
from astropy.io import fits

from .. import __version__

class ApFixCosmicRays:
    """
    Clean cosmic rays from a CCD image, either as a numpy array, or
    from a FITS file, using the ccdproc implementation of the L.A. 
    Cosmic algorithm.
    
    The input image should have already undergone bad pixel
    correction (e.g. using ap_fix_badpix.py) before attempting
    cosmic ray correction, so as to remove most of the 1-pixel
    artifacts. Per bad pixel this algorithm is approximately 250
    times slower than ApFixBadPixels, and only operates over a fixed
    number of iterations, making prior removal of as many bad pixels
    imperative.
    
    The L.A. Cosmic algorithm does a good job of removing most small
    CR artifacts as well as any remaining flickering bad pixel
    artifacts. However it does not do a good job on large transient
    artifacts on the scale of the stellar PSF or larger. These are
    presumably (?) alpha particle hits.
    """
    GOOD     = 0
    AUTO_BAD = 1
    USER_BAD = 2
    
    def __init__(self, loglevel):
        """Constructs an ApFixCosmicRays object.
        
           :param loglevel: Logging level to use.
        """
    
        self._name = 'ApFixCosmicRays'
    
        # Initialize logging
        self._loglevel = loglevel
        self._initialize_logger(self._loglevel)
        
        # Related to FITS file processing
        self._imhdr  = None
        self._imdata = None
        return
                
    def _check_file_exists(self, filename):
        """Checks the file exists and cleans up the path
        """
        
        fpath = Path(filename).expanduser() 
        if not fpath.exists():
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
        return fpath

    def _generate_sigmaclip_mask(self, data, sigma):
        """Creates a bad pixel mask based on sigma-clipped statistics
           of the input data array.
           
        This routine is most appropriate for images that expected to
        be relatively uniform, but with a small number of highly 
        discrepant values.
        
        The generated mask is of type uint8, with value 1 where pixels 
        are BAD and zero where pixels are GOOD. 
                   
        :param data: Input data array
        :param sigma: Number of sigma away from the median to consider
          a pixel bad.
        """
        
        npix = self._imdata.size
        self._logger.debug(f'Generating a bad pixel mask using sigma={sigma} clipping on the input image data values.')
        
        # Compute sigma clipped statistics
        mean, med, std   = sigma_clipped_stats(data, sigma=sigma)
        self._logger.debug(f'Sigma-clipped mean={mean:.2f}, median={med:.2f}, and madstddev={std:.2f} values (ADU).')
        
        lothresh         = med - (sigma * std)
        hithresh         = med + (sigma * std)
        self._logger.info(f'Good pixels have values between {lothresh:.2f} and {hithresh:.2f} ADU.')
        
        # Pixels below low threshold
        bad_lo           = data < lothresh
        nbad_lo          = np.sum(bad_lo)
        self._logger.debug(f'{nbad_lo} (out of {npix}) pixels below {lothresh:.2f} ADU.')
        
        # Pixels above high threshold
        bad_hi           = data > hithresh
        nbad_hi          = np.sum(bad_hi)
        self._logger.debug(f'{nbad_hi} (out of {npix}) pixels above {hithresh:.2f} ADU.')
        
        # Generate final mask in uint8 form.
        self._badpixmask = np.logical_or(bad_lo, bad_hi).astype('uint8')
        nbad             = np.sum(self._badpixmask)
        pct_bad          = 100 * (nbad/npix)
        msg              = f'Out of {npix} pixels, {nbad} are bad ({pct_bad:.4f}%).'
        self._logger.info(msg)
        
        # Update count of algorithmically-defined bad pixels.
        self._nbad_auto = nbad
        return

    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        self._logger = logging.getLogger(self._name)
        
        # Check that the input log level is legal
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: {}'.format(loglevel))
        self._logger.setLevel(numeric_level)
    
        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(numeric_level)
    
        # create formatter
        formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s | %(message)s')
    
        # add formatter to ch
        ch.setFormatter(formatter)
    
        # add ch to logger
        self._logger.addHandler(ch)
        return
            
    def _read_fits(self, image_filename, image_extension):
        """Read a single extension's data and header from a FITS file
        """
        
        image_filename = self._check_file_exists(image_filename)
        self._logger.info('Loading extension {} of FITS file {}'.format(image_extension, image_filename))
            
        # open() parameters that can be important.
        # Default values used here.
        # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
        uint_handling = True
        image_scaling = False
            
        with fits.open(image_filename, 
            uint=uint_handling, 
            do_not_scale_image_data=image_scaling) as hdu_list:
            ext_hdr  = hdu_list[image_extension].header
            ext_data = hdu_list[image_extension].data
            
        ndim     = ext_hdr['NAXIS']
        cols     = ext_hdr['NAXIS1']
        rows     = ext_hdr['NAXIS2']
        bitpix   = ext_hdr['BITPIX']
        info_str = '{}-D BITPIX={} image with {} columns, {} rows'.format(ndim, bitpix, cols, rows)
        
        if ndim == 3:
            layers = ext_hdr['NAXIS3']
            info_str = '{}-D BITPIX={} image with {} columns, {} rows, {} layers'.format(ndim, bitpix, cols, rows, layers)

        if 'BSCALE' in ext_hdr:
            bscale   = ext_hdr['BSCALE']
            info_str += f', BSCALE={bscale}'
        
        if 'BZERO' in ext_hdr:
            bzero    = ext_hdr['BZERO']
            info_str += f', BZERO={bzero}'
            
        self._logger.debug(info_str)
        
        if ndim == 3:
            self._logger.error('Error, 3-D handling has not been implemented yet.')
            sys.exit(1)
            
        # Get data absolute limits.
        minval = np.amin(ext_data)
        maxval = np.amax(ext_data)
        medval = np.median(ext_data)
        self._logger.debug(f'Raw data statistics are min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}')
        
        # Is there a PEDESTAL value? MaximDL likes to add an offset, and
        # the PEDESTAL value is the value to ADD to the data to remove the
        # pedestal.
        if 'PEDESTAL' in ext_hdr:
            pedestal = float( ext_hdr['PEDESTAL'] )
            if pedestal != 0:
                self._logger.debug(f'Removing a PEDESTAL value of {pedestal} ADU.')
                ext_data += pedestal
                minval = np.amin(ext_data)
                maxval = np.amax(ext_data)
                medval = np.median(ext_data)
                self._logger.debug(f'After PEDESTAL removal, min={minval:.2f}, max={maxval:.2f}, median={medval:.2f}')
        
        return ext_data, ext_hdr
                
    def _update_header(self, hdu):
        """Updates the raw mask FITS primary header by adding select
           keywords from the input master dark/bias file.
           
        :param hdu: FITS hdu object to be modified.
        """
        
        self._logger.debug('Updating FITS primary HDU keywords.')
        
        # keyword dictionary to write to header
        kw_dict = self._cr_kw
                
        tnow             = datetime.now().isoformat(timespec='milliseconds')
        creation_date    = datetime.now(timezone.utc)
        creation_datestr = creation_date.isoformat(timespec='seconds')

        kw_dict['CREATOR']  = (self._name, 'Software that generated this file.')
        kw_dict['DATE']     = (creation_datestr, 'UTC creation time.')
    
        for kw in kw_dict:
            hdu.header[kw] = kw_dict[kw]
        hdu.header['HISTORY'] = f'Processed by {self._name} {__version__} at {tnow}'
        return

    def get_crdiff(self):
        """Returns the difference between the original image and the
           comsic-ray cleaned image array.
        """
        return self._crdiff
        
    def get_crmask(self):
        """Returns the cosmic ray mask as a uint8 numpy array.
        """
        return self._crmask
    
    def process(self, inpdata, gain):
        """Apply L.A. Cosmic ray rejection algorithm to the data in the
           input numpy array.
          
        Assumes that the data units are ADU. The output cleaned image is
        also returned in ADU.
        
        :param inpdata: Input numpy data array. 
        :param gain: Gain, in electrons per ADU.
        
        Returns the cleaned image, in the same units as the input,
        and a dictionary of FITS keyword (value, comment) pairs that
        can be used to update a FITS header.
        """

        perf_time_start = time.perf_counter() # highest res timer, counts sleeps

        # Store inputs.
        self._imdata = inpdata
        self._gain   = gain

        # Process
        # TODO allow inputs. Currently set to be closer to iTelescope t05
        # values for t05.
        la_niter     = 6             # number of iterations
        la_readnoise = 12.0          # electrons.
        la_verbose   = False
        la_fwhm      = 3.5           # pixels
        la_satlevel  = gain * 65535  # electrons.
        la_gainapply = True
        la_sigclip   = 4.5           # default 4.5
        la_fsmode    = 'convolve'    # default 'median', but convolve works better.
        
        settings_dict = {'gain': self._gain, 
            'sigclip':      la_sigclip,
            'readnoise':    la_readnoise,
            'psffwhm':      la_fwhm, 
            'verbose':      la_verbose, 
            'gain_apply':   la_gainapply,
            'satlevel':     la_satlevel,
            'niter':        la_niter,
            'fsmode':       la_fsmode}
            
        msg = f'Running cosmicray_lacosmic with the following settings: {settings_dict}'
        self._logger.debug(msg)
        
        # Run the algorithm.
        crimg, crmask = ccdproc.cosmicray_lacosmic(inpdata,
            **settings_dict)

        # Convert back from electrons to ADU
        self._cleandata = crimg / gain

        # Convert mask to uint8
        self._crmask = crmask.astype('uint8')
        
        # TODO stats of CRs. 
        # Currently just count pixels set bad.
        numbad = np.sum(self._crmask)
        self._logger.info(f'{numbad} pixels in image identified as affected by cosmic rays.')

        # Set some metadata that a user can retrieve.
        kw_dict = {'CR_CLEAN': (True, 'Has cosmic ray removal been performed?'),
            'CR_NPIX':  (numbad, 'Number of pixels modified by lacosmic.')}

        # Compute difference image and store
        self._crdiff = self._imdata - self._cleandata      
        self._cr_kw  = kw_dict  
        
        # Performance metric
        perf_time_end = time.perf_counter() # highest res timer, counts sleeps
        run_time_secs = perf_time_end - perf_time_start
        ms_per_pix    = 1000 * run_time_secs / float(numbad)
        self._logger.debug(f'Finished CR processing in {run_time_secs:.3f} s, {ms_per_pix:.3f} ms per CR pixel.')
        return self._cleandata, kw_dict
    
    def process_file(self, inpfile, outfile):
        """
        Apply the L.A. Cosmic CR identification and removal algorithm to
        the FITS image in inpfile, writing the corrected image to outfile.
        
        The input image should have a header keyword specifying the gain
        in e/ADU, either GAIN or EGAIN. If not present then a value of 1.0
        will be assumed. The output image is returned in the same units
        as the original image.
        
        Note: assumes the input image units are ADU.
        
        :param inpfile: Input FITS image.
        :param outfile: Output FITS image, will be over-written if present.
        """
        
        self._imdata, self._imhdr = self._read_fits(inpfile, 0)
        
        # Determine gain
        gain = None
        for kw in ['GAIN', 'EGAIN']:
            if kw in self._imhdr:
                gain = float( self._imhdr[kw] )
                self._logger.debug(f'Read gain value of {gain:.3f} e/ADU from {kw} keyword.')
        if gain is None:
            gain = 1.0
            self._logger.warning(f'Could not find gain value in header. Assuming gain={gain:.3f} e/ADU.')
            
        # Perform CR rejection. (Here we can ignore the returned data as we
        # already store it.)
        self.process(self._imdata, gain)

        # Generate output image, based off header of input image.
        hdu  = fits.PrimaryHDU(data=self._cleandata,
            header=self._imhdr)
        self._update_header(hdu)
        
        hdu_list = fits.HDUList([hdu])
        hdu_list.writeto(outfile, overwrite=True)
        self._logger.info(f'Wrote cosmic ray cleaned image to {outfile}')
        return

    def write_crmask_img(self, mask_file_name):
        """Write the cosmic ray mask to a FITS file with the user 
           specified name/path.
                      
        :param mask_file_name: File name/path for the output cosmic ray
          pixel mask. The file will be overwritten if it exists.
        """
        
        hdu  = fits.PrimaryHDU(data=self._crmask)
        #self._update_header(hdu)
        
        hdu_list = fits.HDUList([hdu])
        hdu_list.writeto(mask_file_name, overwrite=True)
        self._logger.info(f'Wrote cosmic ray pixel mask to {mask_file_name}')
        return
    
    def write_crdiff_img(self, diff_file_name):
        """Write the cosmic ray difference image to a FITS file with the user 
           specified name/path.
                      
        The difference image is the cosmic ray cleaned data subtracted
        from the original image. This will be non-zero only at the location
        of the pixels the algorithm identified as cosmic rays.
                      
        :param diff_file_name: File name/path for the output cosmic ray
          difference image. The file will be overwritten if it exists.
        """
        
        hdu = fits.PrimaryHDU(data=self._crdiff)
        #self._update_header(hdu)
        
        hdu_list = fits.HDUList([hdu])
        hdu_list.writeto(diff_file_name, overwrite=True)
        self._logger.info(f'Wrote cosmic ray difference image to {diff_file_name}')
        return
