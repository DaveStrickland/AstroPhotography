#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_imarith.py
#
#  Perform arithmetic operations on a FITS file, using either a constant
#  value applied to the whole primary array or the pixels values in 
#  another FITS file of the same shape. This is similar to the Heatools
#  farith, or the cfitsio example imarith.
#  
#  Copyright 2021 Dave Strickland <dave.strickland@gmail.com>
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
#  2021-05-02 dks : Initial skeleton. 

import argparse
import sys
import logging
import AstroPhotography as ap

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_imarith',
        description=('Perform arithmetic operations on a FITS file,'
            ' using either a constant value applied to the whole primary'
            ' array or the pixels values in another FITS file of the same'
            ' shape. This is similar to the Heatools farith, or the'
            ' cfitsio example imarith.'))
        
    allowed_ops = ['ADD', 'SUB', 'MUL', 'DIV']
        
    # Required
    parser.add_argument('input_image',
        metavar='INPUT_IMAGE.FITS',
        help='Path/name of the input image to perform arithmetic on.')
    parser.add_argument('operation',
        metavar='OPERATION',
        help=('The mathematical operation to perform, represented as'
            ' a three letter uppercase string. Allowed values are:'
            f' {allowed_ops}'))
    parser.add_argument('value',
        metavar='VALUE_OR_IMAGE',
        help=('A floating point value OR the path to another fits image'
            ' of the same shape as the input image. If an image is used'
            ' its data type will be converted to double while performing'
            ' any calculations.'))
    parser.add_argument('output_image',
        metavar='OUTPUT_IMAGE.FITS',
        help=('Path/name of the output image.'
            ' The output image will retain the structure and keyords of'
            ' the input image. If --units is specified the output image'
            ' BUNIT record will be updated. The operation this script'
            ' performed is recorded as a HISTORY header record.'
            ' The data type of the output image will be the same as'
            ' that of input image.'))        
        
    # Optional
    parser.add_argument('--units',
        default=None,
        help=('String specifying the output image units, e.g. "adu/s".'
        ' The BUNIT fits header keyword will be set to this value if'
        ' specified. If not specified then the value of the BUNIT'
        ' keyword in the original input image (if any), will be used.'))
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    return args

def main(args=None):
    p_args       = command_line_opts(args)
    p_inp_img    = p_args.input_image
    p_operation  = p_args.operation
    p_value      = p_args.value
    p_out_img    = p_args.output_image
    
    p_units      = p_args.units
    p_loglevel   = p_args.loglevel
    
    # Create an instance of the object that can fix bad pixels.
    imarith = ap.ApImArith(p_loglevel)
    
    # Tell it to fix data from files (rather than numpy arrays).
    imarith.process_files(p_inp_img,
        p_operation,
        p_value,
        p_out_img,
        p_units) 
    
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
