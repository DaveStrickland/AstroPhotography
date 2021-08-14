"""Contains the implementation of the ApAddMetadata class.
"""

# 2020-12-16 dks : Initial implementation.
# 2020-12-20 dks : Working version.
# 2021-08-14 dks : Add capabilty to remove Telescopius mosaic suffixes
#                  and also accept a user-supplied target name for name
#                  resolution

import sys
import logging
from pathlib import Path
import math
import time
from datetime import datetime, timezone
import re

import numpy as np

from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation

from astroplan import Observer
from astroplan import FixedTarget

from .. import __version__

class ApAddMetadata:
    """
    Adds metadata to the FITS header of calibrated fits files.
    
    Currently this is written to use the information in iTelescope
    file names, plus information derived from that using astropy.
    Successful processing depends on successful target name resolution. 
    
    In cases where it is known that the target name in the FITS file
    name will fail CDS/Simbad name resolution the user may supply the
    target name as a string.
    
    You can check target name resolution at http://simbad.u-strasbg.fr/simbad/sim-fid
    
    Mosaic pointings, particularly those planned using Telescopius, often
    append the sub-pointing or field identifier as a predictable suffix 
    that can cause name resolution to fail. To avoid having to specify 
    the target name by hand this class will attempt to automatically 
    detect such cases and remove the suffix.
    
    For example the first field of a set of mosaic pointing of the 
    Cygnus Loop might have a file name including `CygnusLoop_x1_y1.`
    The `x* y*` suffix will be stripped so that only `CygnusLoop` 
    is sent to CDS/Simbad for name resolution.
    """
    
    def __init__(self,
        loglevel):
        """Initializes an ApAddMetadata instance.
        """
        
        self._name     = 'ApAddMetadata'
        self._version  = __version__
        self._loglevel           = loglevel
        self._initialize_logger(self._loglevel)
        return
        
    def _check_file_exists(self, filename):
        """Raises an exception if the name file does not exist.
        """
        
        if not Path(filename).exists():
            err_msg = f'Cannot find {filename}. Not a valid path or file.'
            self._logger.error(err_msg)
            raise RuntimeError(err_msg)
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
        
        # Used in cases where we get the same message twice or more
        # See https://stackoverflow.com/a/44426266
        self._logger.propagate = False
        return
        
    def _get_itelescope_site(self, telescope):
        """Return an astroplan Observer instance for a given iTelescope
           telescope.
        """
        
        nm_loc = {'mpc': 'H06',
            'tz':   'UTC Minus 7:00. Daylight savings time is observed.',
            'lat':  '+32d54m11.91s',
            'lon':  '-105d31m43.32s',
            'elev': 2222 * u.m} 
        es_loc = {'mpc':    'I89',
            'tz':   'UTC +1:00 Madrid Daylight Savings Time is Observed',
            'lat':  '+38d09m56s', 
            'lon':  '-2d19m37s',  
            'elev': 1607 * u.m}
        au_loc = {'mpc':    'Q62',
            'tz':   'UTC +10:00 New South Wales, Australia Daylight savings time is observed.',
            'lat':  '-31d16m24s', 
            'lon':  '149d04m11s',
            'elev': 1118 * u.m}
        ca_loc = {'mpc':    'U69',
            'tz':   'UTC Minus 8:00. Daylight savings time is observed.',
            'lat':  '37d04m13s', 
            'lon':  '-119d24m47s',
            'elev': 1403 * u.m}
        
        # Dictory of iTelescope to site. Note, used lower case
        tel_site_dict = {'t02': 'mayhill',
            't02': 'mayhill',
            't05': 'mayhill',
            't11': 'mayhill',
            't14': 'mayhill',
            't20': 'mayhill',
            't21': 'mayhill',
            't68': 'mayhill',
            't24': 'auberry',
            't08': 'sidingspring',
            't09': 'sidingspring',
            't12': 'sidingspring',
            't17': 'sidingspring',
            't30': 'sidingspring',
            't31': 'sidingspring',
            't32': 'sidingspring',
            't33': 'sidingspring',
            't07': 'nerpio',
            't16': 'nerpio',
            't18': 'nerpio'}
        
        mayhill_nm = EarthLocation.from_geodetic(nm_loc['lon'],
            nm_loc['lat'], nm_loc['elev'])
        nerpio_es = EarthLocation.from_geodetic(es_loc['lon'],
            es_loc['lat'], es_loc['elev'])
        sidingspring_au = EarthLocation.from_geodetic(au_loc['lon'],
            au_loc['lat'], au_loc['elev'])
        auberry_ca = EarthLocation.from_geodetic(ca_loc['lon'],
            ca_loc['lat'], ca_loc['elev'])
        
        mayhill_obs = Observer(name='iTelescope New Mexico',
               location=mayhill_nm,
               description="iTelescope at Mayhill, NM")
        nerpio_obs = Observer(name='iTelescope Astrocamp',
               location=nerpio_es,
               description="iTelescope at Nerpio, Spain")
        sidingspring_obs = Observer(name='iTelescope Siding Spring',
               location=sidingspring_au,
               description="iTelescope at Siding Spring, Australia")
        auberry_obs = Observer(name='iTelescope Sierra Remote',
               location=auberry_ca,
               description="iTelescope at Auberry, CA")
        
        err_msg = None
        if telescope.lower() in tel_site_dict:
            loc = tel_site_dict[telescope.lower()]
            if 'mayhill' in loc:
                observatory = mayhill_obs
            elif 'auberry' in loc:
                observatory = auberry_obs
            elif 'nerpio' in loc:
                observatory = nerpio_obs
            elif 'sidingspring' in loc:
                observatory = sidingspring_obs
            else:
                err_msg = f'Error, unexpected location {loc} found for telescope {telescope.lower()}.'
        else:
            err_msg = f'Error, telescope {telescope.lower()} not in list of iTelescope telescopes.'
        
        if err_msg is not None:
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)
        
        return observatory
        
        
    def _parse_itelescope_filename(self, filename):
        """
        Extract the telescope name, observer, and target from an
        iTelescope-generated FITS file name.
           
        The input file name is expected to be of similar form to the
        raw iTelescope files, e.g. 
        `raw-T05-davestrickland-NGC_6888-20200716-231744-Ha-BIN1-E-180-001.fit`.
        In particular dashes (`-`) are used as a field separator, and there
        is only a single field before the telescope string. In this
        particular case this function will return telescope=T05,
        observer=davestrickland, and target=NGC_6888.
        
        :param filename: File name string to process. 
        """
        
        split_list = filename.split('-')
        num_fields = len(split_list)
        if num_fields > 3:
            telescope  = split_list[1]
            observer   = split_list[2]
            # Replace underscores with spaces.
            target     = split_list[3].replace('_', ' ')
            
            # Check if target contains a Telescopius-like mosaic suffix
            # of the form '* x\d+ y\d+', in which case we want to strip
            # the trailing x and y parts.
            mosaic_re = re.compile(' x\d+ y\d+')
            mosaic_check = mosaic_re.search(target)
            if mosaic_check:
                # Strip off the target part.
                self._logger.debug(f'Mosaic suffix pattern found in file-based target name [{target}]')
                target = target[0:mosaic_check.start()]
        else:
            err_msg = (f'Error, splitting the file name {filename}'
            f' only generated {num_fields} fields, expecting > 3 fields.'
            f' Split file name: {split_list}')
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)
        
        return telescope, observer, target
        
    def _read_fits(self, image_filename, image_extension):
        """Read a single extension's data and header from a FITS file.
        """
        
        self._check_file_exists(image_filename)
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
            
        return ext_data, ext_hdr

    def _write_corrected_header(self, fitsfile, kwdict):
        """Updates the header keywords of the specified file with 
           the computed values.
        """
        
        self._logger.debug(f'FITS header keywords to be added to output: {kwdict}')
        self._check_file_exists(fitsfile)
            
        # open() parameters that can be important.
        # Default values used here.
        # See https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open
        uint_handling = True
        image_scaling = False
        ext_num       = 0
            
        # Re-open the file in update mode.
        with fits.open(fitsfile, mode='update', 
            uint=uint_handling, 
            do_not_scale_image_data=image_scaling) as hdu_list:
                        
            # Modify header
            for kw, val in kwdict.items():
                hdu_list[ext_num].header[kw] = val
            
            tnow = datetime.now().isoformat(timespec='milliseconds')
            hdu_list[ext_num].header['HISTORY'] = f'Processed by {self._name} {self._version} at {tnow}'

            # Flush the changes back to disk.
            hdu_list.flush()
                    
        self._logger.info(f'Finished updating FITS header of {fitsfile}')
        return
    
    def process(self, fitsfile, mode, user_target_name=None):
        """
        Adds or updates FITS header keywords in the specified FITS file.
        
        This function generates or populates the values of FITS header
        keywords that are required by ApFindStars, ApAstrometry and
        ApQualitySummarizer.
        
        This function supports the following modes:
        - 'iTelescope': The telescope, target, and observer are parsed 
          from the file name. In combination with the date/time of the
          observation this allows the target RA, Dec, site 
          lat/lon/elevation, and airmass to be determined.
        
        :param fitsfile: Path/name of the FITS file to be modified.
        :param mode: Mode. Must be one of the modes specified above.
        :param user_target_name: If specified this string will be used
          for target name and coordinate resolution instead of the being
          parsed from the FITS file name.
        """
        
        self._logger.info(f'Adding metadata to {fitsfile}, mode={mode}.')
        
        # FITS keyword (value, comment) dictionary
        kwdict = {}
        
        telescope_str = None
        target_str    = None
        observer_str  = None
        site          = None
        target        = None
        if 'iTelescope' in mode:
            telescope_str, observer_str, target_str = self._parse_itelescope_filename(fitsfile)
            
            if user_target_name is not None:
                self._logger.debug(f'User specified target name [{user_target_name}] will be used instead of filename derived target name [{target_str}].')
                target_str = user_target_name
                
            self._logger.info(f'Telescope={telescope_str}, observer={observer_str}, and target={target_str}.')
            site          = self._get_itelescope_site(telescope_str)
            target        = FixedTarget.from_name(target_str)
            telescope_str = 'iTelescope ' + telescope_str
        else:
            err_msg = f'Error, unexpected/unsupported mode {mode}.'
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)
        
        # Observer-related keywords
        kwdict['OBSERVER'] = (observer_str, 'Name of observer')
        
        # observatory related keywords:
        lat      = site.location.lat.deg
        lon      = site.location.lon.deg
        elev     = site.location.height.value
        site_str = site.name
        kwdict['OBSERVAT'] = (site_str, 'Observatory.')
        kwdict['LAT-OBS']  = (lat, '[deg] Latitude of observatory.')
        kwdict['LON-OBS']  = (lon, '[deg] Latitude of observatory.')
        kwdict['ALT-OBS']  = (elev, '[m] Height of observatory.')
        kwdict['TELESCOP'] = (telescope_str, 'Name of telescope used.')
        
        # Target related keywords
        kwdict['OBJECT']  = (target_str, 'Target of observation')
        kwdict['OBJNAME'] = kwdict['OBJECT']
        kwdict['RA-OBJ']  = (target.ra.deg, '[deg] Right Ascension of target')
        kwdict['DEC-OBJ'] = (target.dec.deg, '[deg] Declination of target')
        
        # Get FITS header and date of observation
        ext_num = 0
        fdata, fhdr = self._read_fits(fitsfile, ext_num)
        date_obs = None
        if 'DATE-OBS' in fhdr:
            date_obs = Time(fhdr['DATE-OBS'])
            self._logger.debug(f'Date of observation start: {date_obs}')
        else:
            err_msg = f'Error, {fitsfile} header did not contain DATE-OBS keyword.'
            self._logger.error(err_msg)
            raise RunTimeError(err_msg)

        # Airmass
        airmass = site.altaz(date_obs, target).secz.value
        kwdict['AIRMASS'] = (airmass, 'Airmass at start of observation')

        # Update original file.
        self._write_corrected_header(fitsfile, kwdict)
        self._logger.debug(f'Finished processing {fitsfile}')
        return
