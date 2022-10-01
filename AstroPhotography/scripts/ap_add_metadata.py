#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_add_metadata.py
#
#  Adds metadata to the FITS header of a calibrated fits file.
#  Currently this is written to use the information in iTelescope
#  file names, plus information derived from that using astropy.
#  In the longer term other ways of inputting the data could be added.
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
#  2020-12-16 dks : Initial skeleton.
#  2021-08-12 dks : Updated documentation regarding name resolution.
#                   Added --target option to match ApAddMetadata.py
#  2022-09-30 dks : Added yamlkeyval mode and yamlfile.

import argparse
import sys
import logging
import AstroPhotography as ap
import astropy.coordinates.name_resolve

def command_line_opts(argv):
    """Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_add_metadata',
        description=('Adds metadata to the FITS header of a calibrated fits file.'
        ' Currently two methods of adding metadata exist.'
        ' The iTelescope mode extracts information from iTelescope FITS'
        ' file names, plus information derived from that using astropy.'
        ' In the longer term other ways of inputting the data could be added.'
        ' The target name must resolve correctly using the Object search from'
        ' http://simbad.u-strasbg.fr/simbad/sim-fid or processing will fail.'
        ' This will not work on iTelescope premium images. In that and more'
        ' general cases use the yamlkeyval mode along with a yaml file of'
        ' key/value pairs.'))
    
    # Required
    parser.add_argument('fitsimage',
        metavar='FITSIMAGE.FITS',
        help=('Path/name of the FITS image to added metadata to.'
        ' Note that the input file is modified.'))
        
    # Optional
    p_mode = 'iTelescope'

    parser.add_argument('--mode',
        metavar='MODE',
        default=p_mode,
        help=('Mode of operation. If "iTelescope" the telescope, target,'
        ' and observer will be parsed from the file name.'
        ' If the target name parsed from the iTelescope files contains'
        ' mosaic-like values, e.g. CygnusLoop x1 y1, then the x* and y*'
        ' parts will be dropped prior to attempted name resolution.'
        ' If "yamlkeyval" then "--yamlfile" must also be specified.'))
    parser.add_argument('--target',
        metavar='TARGET NAME STRING',
        default=None,
        help=('A target name to be used for name resolution in cases'
        ' where the file name based target name fails or will fail.'
        ' This is a string so use quotation marks to enclose names'
        ' with spaces. If specified this name is used irrespective of'
        ' whether the file name based target name would work or not.'))
    parser.add_argument('--yamlfile',
        metavar='YAMLFILE',
        default=None,
        help=('If "mode==yamlkeyval" then this parameter must be the'
        ' name/path of a yaml file containing `key: value` pairs,'
        ' order one per line. The keys will be converted to upper'
        ' case and used as new FITS header keywords with the specified'
        ' value. Note that these can overwrite existing FITS header'
        ' keywords.'))
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    
    if ('yamlkeyval' in args.mode) and (args.yamlfile is None):
        parser.error('--yamlfile must be specified if "mode=iTelescope"')
    return args

def main(args=None):
    p_args       = command_line_opts(args)
    p_fitsimg    = p_args.fitsimage
    p_mode       = p_args.mode
    p_target     = p_args.target
    p_yamlfile   = p_args.yamlfile
    p_loglevel   = p_args.loglevel
    retcode      = 0
    
    # Create an instance of the calibrator.
    meta_adder = ap.ApAddMetadata(p_loglevel)
    
    # Update the image.
    try:
        meta_adder.process(p_fitsimg, p_mode, p_target, p_yamlfile)
    except astropy.coordinates.name_resolve.NameResolveError as err:
        retcode = 1
        # TODO: Have logger pick up log format used in rest of AstroPhotography?
        logging.getLogger(__name__).error(f"Name resolution failed: {err}")
    
    return retcode

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
