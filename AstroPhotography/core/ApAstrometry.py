# -*- coding: utf-8 -*-
#
#  ApAstrometry.py
#  
#  Copyright 2020-2021 Dave Strickland <dave.strickland@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

# 2020-08-12 dks : Initial WIP
# 2020-08-29 dks : Update to keywords read from source list.
# 2020-12-20 dks : Fix read fits when BZERO undefined, increase timeout.
# 2020-12-28 dks : Disable SIP by default because swarp does not support it.
# 2021-01-19 dks : Move ApAstrometry into core.
#                  Modify how scale estimates are used.

import sys
import logging
import numpy as np
import math
import os.path
import warnings
import json

from astropy.io import fits
from astropy import wcs
from astropy.table import QTable, Table
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy import units as u
from astropy import log as aplog
aplog.setLevel('ERROR')
from astroquery.astrometry_net import AstrometryNet
aplog.setLevel('INFO')
import astroquery.exceptions
    
class ApAstrometry:
    """
    Uses astrometry.net and a list of valid stars to calculate an
    astrometric solution to a given FITS image, writing a copy of
    the image with valid WCS keywords.
    
    If the astrometric solution is obtained the photometry table
    in the input sourcelist will also be updated with source RA
    and Dec values.
    """
    
    NOMINAL     = 0 #: Nominal or WCS solution obtained
    INPUT_ERROR = 1 #: Non-nominal input, missing file, etc
    NO_SOLUTION = 2 #: Astrometry.net did not find a solution
    
    def __init__(self, inp_img_fname, 
        inp_img_extnum,
        srclist_fname, 
        srclist_extname,
        out_img_fname, 
        astnet_key,
        use_sip,
        user_scale=None,
        scale_err_ratio=None,
        loglevel="INFO"):
        
        self._initialize_logger(loglevel)
        self._loglevel = loglevel
        self._status   = ApAstrometry.NOMINAL
        self._use_sip  = use_sip
        
        # User-define pixel scale (arcsec/pixel)? If None then
        # use pixel scale estimate from sourcelist
        self._user_scale = user_scale
        
        # Pixel scale error ratio to use. The default is equivalent
        # to a 30% scale error.
        if scale_err_ratio is None:
            self._scale_err_ratio = 1.3
        else:
            self._scale_err_ratio = scale_err_ratio

        # Output image cannot be the same as the input to avoid 
        # over-writing it.
        if (inp_img_fname == out_img_fname):
            err_msg = f'Output image cannot be the same as the input image. Choose a different name or path.'
            self._logger.error(err_msg)
            self._status = ApAstrometry.INPUT_ERROR
            return
        
        # Read input image
        self._inp_fname  = inp_img_fname
        self._inp_extnum = inp_img_extnum
        self._idata, self._ihdr = self._read_fits_image(inp_img_fname, 
            inp_img_extnum)
        self._cols = self._idata.shape[1]
        self._rows = self._idata.shape[0]
            
        # Read star position data 
        self._src_fname = srclist_fname
        self._src_xytable, self._src_meta = self._read_srclist(srclist_fname,
            srclist_extname)
            
        # Check input image corresponds to the image the source list was
        # generated from.
        self._sanity_check(self._inp_fname, self._src_meta)
        
        # Build additional astrometry settings based on source metadata
        # if possible. This should improve the speed of the solution.
        self._astnet_hint_dict = self._generate_hints(self._src_meta)
        
        # Run query 
        self._astnet_key = astnet_key
        self._wcs        = self._query_astrometry_dot_net(self._astnet_key,
            self._cols,
            self._rows,
            self._src_xytable,
            self._astnet_hint_dict)
        
        # If a WCS solution was found then create a copy of the input
        # image and add parts of the WCS header data to it.
        self._out_fname  = out_img_fname
        if len(self._wcs) > 0:
            self._write_fits_image(self._wcs)        # Creates _out_fname
            self._update_sourcelist(self._src_fname, 
                self._out_fname)
        else:
            self._logger.warn(f'No output as no WCS found. Output image {out_img_fname} NOT created.')
            self._status = ApAstrometry.NO_SOLUTION

        return
    
    def _initialize_logger(self, loglevel):
        """
        Initialize and return the logger
        """
        
        self._logger = logging.getLogger(__name__)
        
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
        
    def _check_file_exists(self, filename):
        if not os.path.isfile(filename):
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            self._status = ApAstrometry.INPUT_ERROR
            
    def _generate_hints(self, srclist_hdr):
        """
        Generate a dictionary of additional hints for astrometry.net
        
        This relies on keywords written to the primary header of the
        source list by ap_find_stars.py. If the expected keywords are not
        present then the defaults should still work, but the solution may
        take longer.
        
        This function returns a dictionary corresponding to valid
        astroquery astrometry net parameter/value pairs. See
        https://astroquery.readthedocs.io/en/latest/astrometry_net/astrometry_net.html
        """
        
        hints_dict   = {}
        aprx_ra      = None
        aprx_dec     = None
        aprx_fov     = None
        aprx_xpixsiz = None
        aprx_ypixsiz = None

        if 'RA-OBJ' in srclist_hdr:
            aprx_ra      = float( srclist_hdr['RA-OBJ'] )
        elif 'APRX_RA' in srclist_hdr:
            aprx_ra      = float( srclist_hdr['APRX_RA'] )

        if 'DEC-OBJ' in srclist_hdr:
            aprx_dec     = float( srclist_hdr['DEC-OBJ'] )
        elif 'APRX_DEC' in srclist_hdr:
            aprx_dec     = float( srclist_hdr['APRX_DEC'] )

        if self._user_scale is None:
            if 'APRX_FOV' in srclist_hdr:
                aprx_fov     = float( srclist_hdr['APRX_FOV'] )

            if 'APRX_XPS' in srclist_hdr:
                aprx_xpixsiz = float( srclist_hdr['APRX_XPS'] )

            if 'APRX_YPS' in srclist_hdr:
                aprx_ypixsiz = float( srclist_hdr['APRX_YPS'] )
        else:
            # Estimate image size from user-supplied scale and number
            # or rows and cols.
            if 'IMG_COLS' in srclist_hdr:
                img_cols = int( srclist_hdr['IMG_COLS'] )
            else:
                # Reasonable guess
                img_cols = 4096
            if 'IMG_ROWS' in srclist_hdr:
                img_rows = int( srclist_hdr['IMG_ROWS'] )
            else:
                # Reasonable guess
                img_cols = 4096
            xsiz = img_cols * self._user_scale / 3600.0 # deg
            ysiz = img_rows * self._user_scale / 3600.0 # deg
            aprx_fov = math.sqrt(xsiz*xsiz + ysiz*ysiz)
            aprx_xpixsiz = self._user_scale
            aprx_ypixsiz = self._user_scale
                            
        # center_ra, center_dec, radius: all must be supplied
        if (aprx_ra is not None) and (aprx_dec is not None):
            hints_dict['center_ra']  = aprx_ra
            hints_dict['center_dec'] = aprx_dec
            if aprx_fov is None:
                # Make a guess. Assuming this is an iTelescope image
                # the FOV cannot be more than approx 4 degrees.
                aprx_fov = 4.0
            # Place an upper bound on the radius assuming the RA, Dec
            # could be at the edge of the image, and the FOV could be 
            # wrong by a factor 2.
            radius = math.ceil(aprx_fov * 1.5 * self._scale_err_ratio)
            hints_dict['radius'] = radius
            self._logger.debug(f'Estimated center_ra={aprx_ra:.3f}, center_dec={aprx_dec:.3f}, radius={radius:.3f} deg.')
        else:
            self._logger.warning(f'Could not estimate center_ra, center_dec, radius.')
            
        # scale_units, scale_type, scale_est, scale_err
        if (aprx_xpixsiz is not None) and (aprx_ypixsiz is not None):
            scale_units               = 'arcsecperpix'
            hints_dict['scale_units'] = scale_units
            mean_pix_size = math.sqrt( (aprx_xpixsiz*aprx_xpixsiz +
                aprx_ypixsiz*aprx_ypixsiz) / 2 )

            # Giving a range works better than a value with a fixed 
            # +/-percent_error, so we use the 'ul' mode and not the
            # 'ev' mode.
            lolim_pix_size = mean_pix_size / self._scale_err_ratio
            uplim_pix_size = mean_pix_size * self._scale_err_ratio

            hints_dict['scale_type']  = 'ul'
            hints_dict['scale_lower'] = lolim_pix_size
            hints_dict['scale_upper'] = uplim_pix_size
            self._logger.debug(f'Estimated plate scale={mean_pix_size:.3f}, range={lolim_pix_size:.3f} to {uplim_pix_size:.3f} {scale_units}')
        else:
            self._logger.warning(f'Could not generate scale hints.')
        
        
        self._logger.debug(f'Astrometry.net hints dictionary: {hints_dict}')
        return hints_dict
    
    def _read_fits_image(self, image_filename, image_extension):
        """Read a single extension's data and header from a FITS file
        """
            
        self._check_file_exists(image_filename)
        self._logger.info('Loading extension {} of FITS file {}'.format(image_extension, image_filename))
            
        # open() parameters that can be important.
        # Default values used here.
        # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
        uint_handling = True
        image_scaling = False
            
        with fits.open(image_filename,
            mode='readonly',
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
            
        return ext_data, ext_hdr
        
    def _read_srclist(self, srclist_fname, srclist_extname):
        """
        Read the XY positions from a source list generated by ap_find_stars
        """

        self._check_file_exists(srclist_fname)
         
        # Open the source list
        with fits.open(srclist_fname) as hdu_list:
            # Read keywords from the primary header
            srclist_hdr = hdu_list[0].header
            
            # Search for the srclist_extname and read the data.
            if srclist_extname in hdu_list:
                xy_table= Table.read(hdu_list[srclist_extname])
            else:
                err_msg = f'{filename} does not contain the expected {srclist_extname} extension.'
                self._logger.error(err_msg)
                raise RuntimeError(err_msg)

        num_srcs = len(xy_table)
        self._logger.info(f'Read {num_srcs} source positions from {srclist_fname}')
        return xy_table, srclist_hdr
        
    def _query_astrometry_dot_net(self, astnetkey, cols, rows, xy_table, hints_dict):
        """
        Attempt to get an astrometric solution using astrometry.net
        
        Returns a FITS header if successful, an empty dictionary if not.
        """

        # Try to stop astroquery astrometry warnings. This won't work
        # yet, because the astrometry.net package is using logging instead
        # of warnings.
        warnings.filterwarnings("ignore", module='astroquery.astrometry_net')
        warnings.filterwarnings("ignore", module='astroquery.astrometry_net.core')
        warnings.filterwarnings("ignore", module='astroquery.astrometry_net.AstrometryNet')
    
        wcs           = {}
        image_width   = cols
        image_height  = rows
        aplog.setLevel('ERROR')
        ast           = AstrometryNet()
        aplog.setLevel('INFO')
        if astnetkey is not None:
            # Use user-specified key.
            ast.api_key   = astnetkey
        
        # Check that we have a key now    
        if not ast.api_key:
            # Empty string evaluates False
            err_msg = 'Your astroquery config file does not contain an Astrometry.net key and you did not supply one.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
            
        try_again     = True
        submission_id = None
        pos_err_pix   = 10               # TODO get better estimate from srclist?
        timeout       = 180
    
        sip_order = 0
        if self._use_sip:
            sip_order = 2
            self._logger.debug(f'Allowing fitting of SIP distortion polynomial of order {sip_order}')
            self._logger.warning('Some downstream software, e.g. swarp, may not handle SIP correctly.')
    
        # TODO: Robustness...
        # This block is a modified version of the example online.
        # - Both the original example and this version seem to fail to try
        #   again.
        # - If Astrometry.net is not responding we can different exceptions
        #   than just the JSONDecodeError currently handled.
        while try_again:
            try:
                if submission_id is None:
                    self._logger.debug('Submitting astrometry.net solve from source list with {} positions'.format(len(xy_table)))
                    wcs_header = ast.solve_from_source_list(xy_table['X'], 
                        xy_table['Y'],
                        image_width, 
                        image_height,
                        solve_timeout=timeout,
                        parity=2,
                        positional_error=pos_err_pix,
                        crpix_center=True,
                        publicly_visible='n',
                        tweak_order=sip_order,
                        submission_id=submission_id,
                        **hints_dict)
                else:
                    self._logger.debug('Monitoring astrometry.net submission {}'.format(submission_id))
                    try_again  = False
                    wcs_header = ast.monitor_submission(submission_id,
                        solve_timeout=timeout)
            except json.JSONDecodeError as e:
                err_msg = 'Caught JSONDecodeError. Usually this means Astrometry.net is down or login failed.'
                self._logger.error(err_msg)
                raise RuntimeError(err_msg)
            except astroquery.exceptions.TimeoutError as e:
                if (submission_id is not None) and (try_again):
                    self._logger.warning(f'Astronomy solve (id={submission_id}) from source list timed out after {timeout} seconds.')
                    submission_id = e.args[1]
            else:
                # got a result or failed twice, so terminate
                try_again = False
        
        if wcs_header:
            self._logger.info('Obtained a WCS header')
            wcs = wcs_header
        else:
            self._logger.error('Astrometry.net submission={} failed'.format(submission_id))
    
        return wcs
    
    def _sanity_check(self, image_filename, srclist_metadata):
        """
        Warn user if the sourcelist was not produced using the input image.
        """
        
        key = 'IMG_FILE'
        if key in srclist_metadata:
            src_img_file = srclist_metadata[key]
            if image_filename in src_img_file:
                msg = 'The source list was generated from this input image. Proceding.'
                self._logger.debug(msg)
                return
            else:
                msg = f'Source list NOT generated from the input image, but instead from {src_img_file}'
        else:
            msg = f'Source list file lacks {key} keyword. Cannot assess if it was produced from {image_filename}'
        self._logger.warn(msg)
        self._logger.warn('Proceding, but WCS may be incorrect.')
        return
    
    def _update_sourcelist(self, srclist, wcs_image):
        """Update the photometric table with positions from the WCS
           in an image.
        
        This function uses the WCS coordinates in the header of wcs_image
        to update the photometry table in the srclist file.
        """
        
        self._logger.info(f'Updating photometry in {srclist} with WCS soluton.')
        
        # Hardwired extension name
        ext_name = 'AP_L1MAG'
        
        with fits.open(wcs_image, mode='readonly') as im_hdulist:
            # Parse the WCS keywords in the primary HDU
            w = wcs.WCS(im_hdulist[self._inp_extnum].header)
            
            with fits.open(srclist, mode='update') as src_hdulist:
                # Search for the srclist_extname and read the data.
                if ext_name in src_hdulist:
                    phot_table = Table.read( src_hdulist[ext_name] )
                    
                    ra_dec_tupl = w.all_pix2world(phot_table['xcenter'], 
                        phot_table['ycenter'], 
                        0, # origin, zero=based 
                        ra_dec_order=True)
                        
                    phot_table['ra']  = ra_dec_tupl[0]
                    phot_table['dec'] = ra_dec_tupl[1]
                    
                    # Flush changes to 
                    src_hdulist[ext_name] = fits.table_to_hdu(phot_table)
                    src_hdulist.flush()
                else:
                    err_msg = f'{srclist} does not contain the expected {ext_name} extension.'
                    self._logger.error(err_msg)
                    raise RuntimeError(err_msg)
    

        return
    
    def _write_fits_image(self, wcs):
        """
        Create a copy of the input file updated with the WCS solution
        """

        # That should not be copied over from the WCS into the output
        # as they overwrite existing values that are important or they
        # are overly long....
        
        bad_keys = ['SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 
            'HISTORY', 'DATE', 'COMMENT']

        with fits.open(self._inp_fname, mode='readonly') as hdu_list:
            ext_hdr  = hdu_list[self._inp_extnum].header
 
            for key in wcs:
                if key not in bad_keys:
                    val_cmt      = (wcs[key], wcs.comments[key])
                    ext_hdr[key] = val_cmt
                    self._logger.debug(f'Setting {key} to {val_cmt}')
            
            hdu_list[self._inp_extnum].header = ext_hdr
            hdu_list.writeto(self._out_fname, overwrite=True)
        self._logger.info(f'Wrote updated file with WCS to {self._out_fname}')
        return
    
    # Public member functions
    def status(self):
        """Return status"""
        return self._status
    
    def create_primary_header(kw_dict):
        """Creates a FITS primary header adding the keywords in the dict"""
        
        pri_hdr = fits.Header()
        for kw in kw_dict:
            pri_hdr[kw] = kw_dict[kw]
        pri_hdr['HISTORY'] = 'Created by ApFindStars'
        return pri_hdr
        
    def add_optional_keywords(hdr, kw_dict):
        """Add keywords to kw_dict from the FITS header hdr if present"""
        logger = logging.getLogger(__name__)
        
        # The following may exist.
        kw_comment_dict = {'EXPOSURE': '[seconds] Image exposure time',
            'DATE-OBS': 'Observation date and time',
            'OBJECT':   'Target object', 
            'TELESCOP': 'Telescope used',
            'RA':       'Requested right ascension', 
            'DEC':      'Requested declination', 
            'XPIXSZ':   '[micrometers] Stated X-axis pixel scale after binning',
            'YPIXSZ':   '[micrometers] Stated Y-axis pixel scale after binning',
            'FOCALLEN': '[mm] Stated telescope focal length', 
            'FILTER':   'Filter used', 
            'EGAIN':    '[e/ADU] Gain in electrons per ADU'}
        kw_missing_list = []
        for kw in kw_comment_dict:
            if kw in hdr:
                kw_dict[kw] = (hdr[kw], kw_comment_dict[kw])
            else:
                kw_missing_list.append(kw)
    
        logger.debug('FITS keywords found in image: {}'.format(kw_dict))
        logger.debug('FITS keywords missing from image: {}'.format(kw_missing_list))
        return
    
