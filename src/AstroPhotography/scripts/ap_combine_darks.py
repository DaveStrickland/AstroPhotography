#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ap_combine_darks.py
#
#  Generates a master dark or master bias file, given a directory 
#  that contains the raw dark or bias file to be combined. Use separate
#  directories for calcilbration files with different exposure time,
#  binning, CCD temperature, and from different telescopes.
#  
#  Copyright 2020 Dave Strickland <dave.strickland@gmail.com>
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
#  2020-10-01 dks : Initial coding begun. 

import argparse
import sys
import logging
from pathlib import Path
import math

import numpy as np
import matplotlib                # for rc
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import mad_std

import ccdproc as ccdp

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_combine_darks',
        description='Generates a master dark or master bias file from' + 
        ' all calibration FITS files in a given directory.')
        
    # Required
    parser.add_argument('rawcaldir',
        metavar='RAW_CAL_DIR',
        help='The directory in which the raw calibration' +
        ' files to be combined are to be found.')
    parser.add_argument('master_filename',
        metavar='MASTER_CAL_FILENAME',
        help='Output file name for master calibration file.')
        
    # Default values
    p_temptol  = 0.5 # 0.5 degree C
    p_telescop = 'UNKNOWN'
    p_exclude  = 'master*'
        
    # Optional
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
    parser.add_argument('--exclude',
        dest='exclude_pattern',
        default=p_exclude,
        metavar='FILE_PATTERN',
        help=('A unix-style file pattern that may be used to exclude'
            ' files in the target directory from being processed by this'
            ' command. Usually this is used to exclude any master'
            ' calibration files from being processed.'
            f' Default: "{p_exclude}"'))
    parser.add_argument('--telescop',
        default=None,
        metavar='TELESCOPE_NAME',
        help=('If the input files TELESCOP keyword is missing or empty,' +
            ' write this string as TELESCOP in the output master calibration' +
            f' file, e.g. "iTelescope 5". Default: {p_telescop}'))
    parser.add_argument('--temptol',
        default=p_temptol,
        metavar='DEGREES_C',
        help=('Allowable temperature tolerance up to which the CCDTEMP'
            ' can vary from the SET-TEMP and still be considered as'
            f' representing that temperature. Default: {p_temptol} C.'))
    
    args = parser.parse_args(argv)
    return args
             
class ApMasterCal:
    """Combines a series of darks, bias or flats into a master dark,
       master bias, or master flat.
       
       This class is currently limited in that it cannot handle the 
       presence of different types of file (e.g. a directory containing
       both darks and biases, or darks of different exposure time or
       binning). If the input directory does contain different types of
       FITS file this class will log it and throw a RunTimeError
       exception.
    """
    
    def __init__(self, rootdir, 
        exclude_pattern,
        telescop, 
        temptol,
        loglevel):
        """Construct and fully initialize an instance of ApMasterCal
        """

        # Initialize logging
        self._loglevel = loglevel
        self._initialize_logger(loglevel)

        self._rootdir   = rootdir
        self._telescop = telescop
        self._temptol  = temptol

        # Set keywords that will be used to sort and check files.
        self._summary_kw = ['file', 'date-obs',
            'telescop', 'imagetyp', 
            'filter', 'exptime', 
            'set-temp', 'ccd-temp', 
            'naxis1', 'naxis2']

        # Generate initial ImageFileCollection
        self._data_dir  = Path(self._rootdir)
        self._file_list = None
        self._files     = self._create_file_collection(self._data_dir,
            exclude_pattern, self._file_list)
            
        # Check files, reread file collection
        list_all        = True
        self._file_list = self._check_files(list_all)
        self._files     = self._create_file_collection(self._data_dir,
            exclude_pattern, self._file_list)
            
    
        self._logger.debug('Finished')
        return

    def _check_files(self, list_all=None):
        """Check the file collection for type, exposure, size, and 
           temperature consistency, generating a final file_list for
           inclusion.
        
        This function checks the files in the existing image file collection
        and returns a list of file names to use in a second call to
        _create_file_collection, unless files of a different type, size
        or exposure are detected.
        
        These checks are of two types, one that will result in program
        termination if files do not match, and one that will result in
        the file being excluding from the file list used in the second
        call to _create_file_collection.
        
        - File type, size, exposure differences: All files must have the
          same exposure time, size (naxis1 and naxis2), and imagetyp. If
          not an exception will be thrown, as we can't programmatically
          decide what the user wants when confronted with (for example)
          a directory of five bias frames and five differently size darks. 
          If some files have `set-temp` but others don't this
          also counts as a fatal inconsistency.
        - CCD temperature: Files not outside the temperature tolerance 
          will be excluded from final file list. If there is a unique 
          `set-temp` then each `ccd-temp` is compared to it, otherwise
          the median `ccd-temp` is computed and each file is compared to
          that.
          
        :param list_all: If True then the unique values of all keywords
          within the _summary_kw list will be logged at DEBUG level
          for diagnostic purposes.
        """
        
        if list_all is not None:
            diag = list_all
        else:
            diag = False
        
        good_file_list = []
        raw_file_list = self._files.values('file')
        
        # Get unique values of all keywords, also logging them for 
        # diagnostic purposes if diag is True
        uniq_dict = {}
        for kw in self._summary_kw:
            uniq_val_list = self._files.values(kw, unique=True)
            uniq_dict[kw] = uniq_val_list
            if diag and ('file' not in kw):
                msg = f'For keyword {kw} there are {len(uniq_val_list)} values: {uniq_val_list}'
                self._logger.debug(msg)                

        # Check the number of unique values of:
        # - imgtype
        # - size naxis1 and naxis2
        # - exposure
        # - set-temp
        # If the number of unique values is > 1 then this is an error
        for kw in ['imagetyp', 'naxis1', 'naxis2', 'exptime', 'set-temp']:
            nuniq = len(uniq_dict[kw])
            if nuniq > 1:
                msg = f'Error, there are {nuniq} unique values of {kw} in the files being processed: {uniq_dict[kw]}'
                self._logger.error(msg)
                raise RuntimeError(msg)
        
        # Check set-temp has a value, if so, use it.
        set_temp_is_set = False
        set_temperature = None
        val = uniq_dict['set-temp'][0]  # Know list has a single value
        if isinstance(val, str): 
            # Is it an empty string?
            if not bool(val.strip()):
                msg = 'No numeric value found for SET-TEMP. Will use median of CCD-TEMP instead.'
                self._logger.warning(msg)
        else:
            try:
                set_temperature = float(val)
                set_temp_is_set = True
                self._logger.debug(f'Using SET-TEMP value of {set_temperature} degrees C.')
            except Exception as e:
                msg = f'Error, could not convert SET-TEMP value of {val} to a float. Will use median of CCD-TEMP instead.'
                self._logger.error(msg)
                self._logger.error(f'The exception that was caught is: {e}')
        
        # Check ccd-temp, convert values to float.
        # - If there is one unique value and its an empty string this
        #   is not ideal but we can assume the data is from a device
        #   such as a camera, and proceed.
        # - Otherwise we'll check each file's temperature individually using
        #   the summary table.
        uniq_temps = []
        nuniq = len(uniq_dict['ccd-temp'])
        if nuniq == 1:
            val = uniq_dict['ccd-temp'][0]
            # An empty string evaluates to False
            if not bool(val.strip()):
                msg = ('No files contain CCD-TEMP metadata. Continuing'
                    ' assuming all files obtained at the same temperature.')
                self._logger.warning(msg)
                return raw_file_list
            
        # Do we need to calculate set_temperature from ccd-temp?
        if not set_temp_is_set:
            set_temperature = np.median(self._files.summary['ccd-temp'])
            self._logger.debug(f'Using median of CCD-TEMP value: {set_temperature} degrees C.')

        # Temperature limits:
        temp_min = set_temperature - self._temptol
        temp_max = set_temperature + self._temptol
        self._logger.info(f'Selecting only files with CCD-TEMP between {temp_min:.2f} and {temp_max:.2f} degrees C.')
            
        # Number of files:
        nfiles = len(self._files.summary)
        good_file_list = []
        for idx in range(nfiles):
            fname = self._files.summary['file'][idx]
            temp  = self._files.summary['ccd-temp'][idx]
            if (temp >= temp_min) and (temp <= temp_max):
                good_file_list.append(fname)
            else:
                self.logger.warning(f'Excluding {fname} as CCD-TEMP={temp:2.f} outside allowed range.')
        
        self._logger.info(f'Updated file list contains {len(good_file_list)} files ({len(raw_file_list)} before filtering).')
        return good_file_list

    def _create_file_collection(self, data_dir, 
        exclude_pattern,
        file_list):
        """Create a FITS ImageFileCollection for a directory, optionally
           including only the files named in fname_list.
           
        :param data_dir: Path to the directory containing the files.
        :param exclude_pattern: Pattern used to exclude certain files.
        :param file_list: None, or list of file names to read instead
          of reading all files in the data_dir (excluding files that
          match the exclude_pattern).
        """
        
        if file_list is not None:
            msg = (f'Looking for FITS files in {data_dir},'
                f' including only files in the list: {file_list}')
        else:
            msg = (f'Looking for FITS files in {data_dir},'
                f' excluding files matching the pattern "{exclude_pattern}"')
        self._logger.info(msg)

        file_collection = ccdp.ImageFileCollection(data_dir, 
            keywords=self._summary_kw,
            glob_exclude=exclude_pattern,
            filenames=file_list)

        self._logger.info(f'Found {len(file_collection.summary)} FITS files matching the constraints.')


        return file_collection
             
    def _initialize_logger(self, loglevel):
        """Initialize and return the logger
        """
        
        self._logger = logging.getLogger('ApMasterCal')
        
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
        if not self._logger.handlers:
            self._logger.addHandler(ch)
        
        # Used in cases where we get the same message twice or more
        # See https://stackoverflow.com/a/44426266
        self._logger.propagate = False
        return
                
def main(args=None):
    p_args       = command_line_opts(args)
    p_root       = p_args.rawcaldir
    p_masterfile = p_args.master_filename
    p_loglevel   = p_args.loglevel
    p_telescop   = p_args.telescop
    p_temptol    = p_args.temptol
    p_exclude    = p_args.exclude_pattern

    logger = logging.getLogger(__name__)

    try:
        mkcal = ApMasterCal(p_root, 
            p_exclude,
            p_telescop, 
            p_temptol,
            p_loglevel)
        ##mkcal.make_master(p_masterfile)
    except RuntimeError as rte:
        logger.error(f'Shutting down due to exception raised by ApMasterCal: {rte}')
        return 1
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
