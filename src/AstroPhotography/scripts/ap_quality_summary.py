#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ap_quality_summary.py
#
#  Calls ApQualitySummarizer to reads the quality files produced 
#  by ApFindStars/ap_find_stars.py and writes a summary CSV.
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
#  2020-09-04 dks : Initial coding.
#  2021-01-19 dks : Move ApQualitySummarizer into core

import argparse
import sys
import logging
from pathlib import Path
import yaml

import numpy as np
from astropy.table import Table
import AstroPhotography as ap

def command_line_opts(argv):
    """ Parse command line arguments.

    :param argv: argument list to parse
    """
    parser = argparse.ArgumentParser(prog='ap_quality_summary',
        description='Generate a CSV summary of all quality yaml files' +
        ' within a specified directory or directory tree.')
        
    # Required
    parser.add_argument('qualfile_dir',
        help='Path of directory under which all quality files are to be found.' +
            ' By default only files within that specific directory will be used.')
    parser.add_argument('summary_file',
        metavar='OUTPUT.CSV',
        help='Name of output summary CSV to generate.')
        
    # Defaults.
    p_qual_pref = 'qual'
    p_qual_suff = '.yaml'
        
    # Optional
    parser.add_argument('--prefix',
        default=p_qual_pref,
        help=f'File prefix for quality files. Default: {p_qual_pref}')
    parser.add_argument('--suffix',
        default=p_qual_suff,
        help=f'File suffix for quality files. Default: {p_qual_suff}')
    parser.add_argument('--walk_tree',
        default=False,
        action='store_true',
        help=' If specified all directories under qualfile_dir' + 
            ' will be traversed looking for quality files.')
    parser.add_argument('-l', '--loglevel', 
        default='INFO',
        help='Logging message level. Default: INFO')
                
    args = parser.parse_args(argv)
    return args
        
def main(args=None):
    p_args      = command_line_opts(args)
    p_qualdir   = p_args.qualfile_dir
    p_sumfile   = p_args.summary_file
    p_loglevel  = p_args.loglevel 
    p_walk_tree = p_args.walk_tree
    p_qual_pref = p_args.prefix
    p_qual_suff = p_args.suffix
    
    summarizer = ap.ApQualitySummarizer(p_qualdir, p_sumfile,
        p_loglevel, p_walk_tree,
        p_qual_pref, p_qual_suff)
    return 0

if __name__ == '__main__':
    try:
        status = main()
    except:
        logging.getLogger(__name__).critical("Shutting down due to fatal error")
        raise  # print stack trace
    else:
        raise SystemExit(status)
