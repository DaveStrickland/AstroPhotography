#!/usr/bin/env python
# coding: utf-8

# # iTelescope premium image set processing, the hard way
# 
# **Author: Dave Strickland**
# 
# **Version: 0.5.2-beta1**
# 
# This notebook illustrates how to process a premium image set from iTelescope using the AstroPhotography package. The python environment used corresponds to a miniconda emvironment `ap-env.yml`.
# 
# The example dataset used is the iTelescope Plan-20 premium dataset of the M101 supernova [SN 2023 ixf](https://en.wikipedia.org/wiki/SN_2023ixf). Note that these files are not provided with this package. However the notebook should work with your own files if you specify a valid fits-file containing directory at the prompt below.
# 
# ## Notebook environment setup

# In[1]:


import os
import pathlib
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time
import math
import subprocess
import shlex

from ccdproc import ImageFileCollection
from astropy.table import Table
from astropy import wcs
import astropy
from astropy.io import fits

import AstroPhotography as ap


# In[2]:


def print_module_version(mod):
    """
    Convenience function to pretty-print Python module mod version.

    Arguments:
        mod {module} -- Imported Python module

    Returns:
        str -- module nam
e and version
    """
    print(f"Using module {mod.__name__:30s}  version: {mod.__version__}")
    return


# In[3]:


# Get Version information
print(f'Python version: {sys.version}')
print_module_version(np)
print_module_version(matplotlib)
print_module_version(astropy)
print_module_version(ap)


# In[4]:


# Enable inline plotting for graphics
# %matplotlib inline
# Set default figure size to be larger
# this may only work in matplotlib 2.0+!
from IPython.core.interactiveshell import InteractiveShell
matplotlib.rcParams['figure.figsize'] = [10.0, 6.0]
# Enable multiple outputs from jupyter cells
InteractiveShell.ast_node_interactivity = "all"


# # Doing it by hand
# 
# ## What files do we have to work with?
# 
# Let us change the working directory from the location of this notebook to a directory holding some iTelescope premium image set files. In my case this is the Plan-20 M101 dataset.
# 
# **Note:** The following processing is fairly low-level and requires a lot of human interaction. In a later section we'll process the same set of files in a more automated way.

# In[5]:


print('Enter the path to the directory containing the premium image set')
wdir = input('Image set path:').strip()
try:
    os.chdir(wdir)
    print('Switched directory to ' + os.getcwd())
except:
    print(f'Error, os.chdir threw an exception changing to {wdir}')
    print('Check that the path you supplied is a valid filesystem path.')
    raise


# We now want a pattern that can be used to select the FITS files we want to process, while ignoring any other files or directories we have in that directory (e.g. processed files, configurations files, and so on).
# 
# Listing the fits files from the shell I get:
# ```bash
# ls ../Plan-20-M101SN/*fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-180s-Lum-1.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-180s-Lum-2.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-180s-Lum-3.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-300s-Blue-1.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-300s-Blue-2.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-300s-Blue-3.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-300s-Green-1.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-300s-Green-2.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-300s-Green-3.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-300s-Red-1.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-300s-Red-2.fits
# ../Plan-20-M101SN/Calibrated-iTelescope-M101_Supernova_2023ixf-300s-Red-3.fits
# ```
# 
# So a valid pattern that picks up the original files is `Calibrated-iTelescope*.fits`

# In[6]:


default_pattern = 'Calibrated-iTelescope*.fits'
print(f'Default file pattern for input files: {default_pattern}')
msg = 'Enter new input file pattern (or return to accept default pattern)'
file_pattern = input(msg).strip() or default_pattern
print(f'Using "{file_pattern}" as the input file pattern.')


# In[7]:


p_dir = pathlib.Path(r'./')
flist = list(sorted(p_dir.glob(file_pattern)))
print(f'{len(flist)} files match pattern "{file_pattern}" in current directory.')
for fpath in flist:
    print(f'  {fpath.name}')


# ### Do these have full headers?
# 
# Lets have a look at the FITS header for the first file in this list.

# In[8]:


hdr = None
with fits.open(flist[0], mode='readonly') as hdulist:
    print(f'Opened {flist[0].name}')
    hdulist.verify('fix')
    print(hdulist.info())
    hdr = hdulist[0].header
    print(f'Printing the FITS header, length={len(hdr)}')
    hdr


# In[9]:


# What keywords do we have (non-comment, non-history)
sorted(list( hdr.keys() ))


# Lets look at the keyword differences in more detail
# 
# ```python
# def fkeywords(fname):
#     with fits.open(fname) as hdulist:
#         hdr=hdulist[0].header
#         print(f'FITS header for primary HDU of {fname}')
#         print( sorted( list( hdr.keys() ) ) )
#     return
# 
# f1='Calibrated-iTelescope-M101_Supernova_2023ixf-180s-Lum-1.fits'
# ```
# 
# On the M101 files we've been looking at:
# ```python
# fkeywords(f1)
# FITS header for primary HDU of Calibrated-iTelescope-M101_Supernova_2023ixf-180s-Lum-1.fits
# ['APTDIA', 'BITPIX', 'DEC', 'EXPTIME', 'EXTEND', 'FILTER', 'FOCALLEN', 'NAXIS', 'NAXIS1', 'NAXIS2', 'OBJCTDEC', 'OBJCTRA', 'OBSTARGET', 'OWNER', 'RA', 'SIMPLE', 'SRCLINK']
# 
# # iTelescope calibration system
# In [17]: f2_itel='calibrated-T05-davestrickland-M82-20210219-002600-Red-BIN1-W-300-001.fit'
# 
# # iTelescope raw file
# In [18]: f2_raw='raw-T05-davestrickland-M82-20210219-002600-Red-BIN1-W-300-001.fit'
# 
# # Calibrated file produced by AstroPhotography from the raw file above...
# In [19]: f2_apcal='cal-T05-davestrickland-M82-20210219-002600-Red-BIN1-W-300-001.fits'
# 
# In [20]: fkeywords(f2_raw)
# FITS header for primary HDU of raw-T05-davestrickland-M82-20210219-002600-Red-BIN1-W-300-001.fit
# ['APTAREA', 'APTDIA', 'BITPIX', 'BSCALE', 'BZERO', 'CBLACK', 'CCD-TEMP', 'CSTRETCH', 'CWHITE', 'DATE-OBS', 'EGAIN', 'EXPOSURE', 'EXPTIME', 'FILTER', 'FLIPSTAT', 'FOCALLEN', 'IMAGETYP', 'INSTRUME', 'JD', 'JD-HELIO', 'NAXIS', 'NAXIS1', 'NAXIS2', 'NOTES', 'OBJECT', 'OBSERVER', 'PEDESTAL', 'READOUTM', 'ROWORDER', 'SBSTDVER', 'SET-TEMP', 'SIMPLE', 'SITELAT', 'SITELONG', 'SWCREATE', 'SWOWNER', 'SWSERIAL', 'TELESCOP', 'TRAKTIME', 'XBINNING', 'XORGSUBF', 'XPIXSZ', 'YBINNING', 'YORGSUBF', 'YPIXSZ']
# 
# In [21]: fkeywords(f2_itel)
# FITS header for primary HDU of calibrated-T05-davestrickland-M82-20210219-002600-Red-BIN1-W-300-001.fit
# ['AIRMASS', 'ALT-OBS', 'APTAREA', 'APTDIA', 'BITPIX', 'BSCALE', 'BZERO', 'CALSTAT', 'CBLACK', 'CCD-TEMP', 'CSTRETCH', 'CWHITE', 'DATE', 'DATE-HP', 'DATE-OBS', 'DEC', 'EGAIN', 'EXPOSURE', 'EXPTIME', 'FILTER', 'FLIPSTAT', 'FOCALLEN', 'HISTORY', 'HISTORY', 'HISTORY', 'HISTORY', 'HISTORY', 'HISTORY', 'HISTORY', 'IMAGETYP', 'INSTRUME', 'JD', 'JD-HELIO', 'LAT-OBS', 'LONG-OBS', 'NAXIS', 'NAXIS1', 'NAXIS2', 'NOTES', 'OBJCTDEC', 'OBJCTRA', 'OBJECT', 'OBSERVAT', 'OBSERVER', 'PEDESTAL', 'PIERSIDE', 'RA', 'RADECSYS', 'READOUTM', 'ROWORDER', 'RollAngle', 'SBSTDVER', 'SET-TEMP', 'SIMPLE', 'SITELAT', 'SITELONG', 'ST', 'SWCREATE', 'SWMODIFY', 'SWOWNER', 'SWSERIAL', 'TELESCOP', 'TIME-OBS', 'TIMESYS', 'TRAKTIME', 'USERNAME', 'UT', 'XBINNING', 'XORGSUBF', 'XPIXSZ', 'YBINNING', 'YORGSUBF', 'YPIXSZ', 'iTelescope', 'iTelescopePlateScaleH', 'iTelescopePlateScaleV']
# 
# In [22]: fkeywords(f2_apcal)
# FITS header for primary HDU of cal-T05-davestrickland-M82-20210219-002600-Red-BIN1-W-300-001.fits
# ['AIRMASS', 'ALT-OBS', 'APTAREA', 'APTDIA', 'BIASCORR', 'BIASFILE', 'BITPIX', 'BPIXCORR', 'BPIXDPIX', 'BPIXFILE', 'BPIXNBAD', 'BPIXNFIX', 'BPIXNREM', 'BPIX_MIN', 'BUNIT', 'CBLACK', 'CCD-TEMP', 'CR_CLEAN', 'CR_NPIX', 'CSTRETCH', 'CWHITE', 'DARKCORR', 'DARKFILE', 'DATE-OBS', 'DEC-OBJ', 'EGAIN', 'EXPOSURE', 'EXPTIME', 'FILTER', 'FLATCORR', 'FLATFILE', 'FLIPSTAT', 'FOCALLEN', 'HISTORY', 'HISTORY', 'HISTORY', 'HISTORY', 'HISTORY', 'HISTORY', 'IMAGETYP', 'INSTRUME', 'JD', 'JD-HELIO', 'LAT-OBS', 'LON-OBS', 'NAXIS', 'NAXIS1', 'NAXIS2', 'NOTES', 'OBJECT', 'OBJNAME', 'OBSERVAT', 'OBSERVER', 'RA-OBJ', 'READOUTM', 'ROWORDER', 'SBSTDVER', 'SET-TEMP', 'SIMPLE', 'SITELAT', 'SITELONG', 'SWCREATE', 'SWOWNER', 'SWSERIAL', 'TELESCOP', 'TRAKTIME', 'XBINNING', 'XORGSUBF', 'XPIXSZ', 'YBINNING', 'YORGSUBF', 'YPIXSZ']
# 
# ```

# After some rearranging by hand we get the following table of header keywords. The premium file is clearly missing a large number of header keywords, although this doesn't tell us which ones are needed and which ones are superfluous.
# 
# |iTelescope raw|iTelescope calibrated|AP calibrated|M101 premium|
# |--------------|---------------------|-------------|------------|
# |              |AIRMASS              |AIRMASS      |            |
# |              |ALT-OBS              |ALT-OBS      |            |
# |APTAREA       |APTAREA              |APTAREA      |            |
# |APTDIA        |APTDIA               |APTDIA       |APTDIA      |
# |              |                     |BIASCORR     |            |
# |              |                     |BIASFILE     |            |
# |BITPIX        |BITPIX               |BITPIX       |BITPIX      |
# |              |                     |BPIXCORR     |            |
# |              |                     |BPIXDPIX     |            |
# |              |                     |BPIXFILE     |            |
# |              |                     |BPIXNBAD     |            |
# |              |                     |BPIXNFIX     |            |
# |              |                     |BPIXNREM     |            |
# |              |                     |BPIX_MIN     |            |
# |BSCALE        |BSCALE               |             |            |
# |BZERO         |BZERO                |             |            |
# |              |                     |BUNIT        |            |
# |              |CALSTAT              |             |            |
# |CBLACK        |CBLACK               |CBLACK       |            |
# |CCD-TEMP      |CCD-TEMP             |CCD-TEMP     |            |
# |              |                     |CR_CLEAN     |            |
# |              |                     |CR_NPIX      |            |
# |CSTRETCH      |CSTRETCH             |CSTRETCH     |            |
# |CWHITE        |CWHITE               |CWHITE       |            |
# |              |                     |DARKCORR     |            |
# |              |                     |DARKFILE     |            |
# |              |DATE                 |             |            |
# |              |DATE-HP              |             |            |
# |DATE-OBS      |DATE-OBS             |DATE-OBS     |            |
# |              |DEC                  |             |DEC         |
# |              |                     |DEC-OBJ      |            |
# |EGAIN         |EGAIN                |EGAIN        |            |
# |EXPOSURE      |EXPOSURE             |EXPOSURE     |            |
# |EXPTIME       |EXPTIME              |EXPTIME      |EXPTIME     |
# |              |                     |             |EXTEND      |
# |FILTER        |FILTER               |FILTER       |FILTER      |
# |              |                     |FLATCORR     |            |
# |              |                     |FLATFILE     |            |
# |FLIPSTAT      |FLIPSTAT             |FLIPSTAT     |            |
# |FOCALLEN      |FOCALLEN             |FOCALLEN     |FOCALLEN    |
# |IMAGETYP      |IMAGETYP             |IMAGETYP     |            |
# |INSTRUME      |INSTRUME             |INSTRUME     |            |
# |JD            |JD                   |JD           |            |
# |JD-HELIO      |JD-HELIO             |JD-HELIO     |            |
# |              |LAT-OBS              |LAT-OBS      |            |
# |              |LONG-OBS             |LON-OBS      |            |
# |NAXIS         |NAXIS                |NAXIS        |NAXIS       |
# |NAXIS1        |NAXIS1               |NAXIS1       |NAXIS1      |
# |NAXIS2        |NAXIS2               |NAXIS2       |NAXIS2      |
# |NOTES         |NOTES                |NOTES        |            |
# |              |OBJCTDEC             |             |OBJCTDEC    |
# |              |OBJCTRA              |             |OBJCTRA     |
# |OBJECT        |OBJECT               |OBJECT       |            |
# |              |                     |OBJNAME      |            |
# |              |OBSERVAT             |OBSERVAT     |            |
# |OBSERVER      |OBSERVER             |OBSERVER     |            |
# |              |                     |             |OBSTARGET   |
# |              |                     |             |OWNER       |
# |PEDESTAL      |PEDESTAL             |             |            |
# |              |PIERSIDE             |             |            |
# |              |RA                   |             |RA          |
# |              |RADECSYS             |             |            |
# |              |                     |RA-OBJ       |            |
# |READOUTM      |READOUTM             |READOUTM     |            |
# |ROWORDER      |ROWORDER             |ROWORDER     |            |
# |              |RollAngle            |             |            |
# |SBSTDVER      |SBSTDVER             |SBSTDVER     |            |
# |SET-TEMP      |SET-TEMP             |SET-TEMP     |            |
# |SIMPLE        |SIMPLE               |SIMPLE       |SIMPLE      |
# |SITELAT       |SITELAT              |SITELAT      |            |
# |SITELONG      |SITELONG             |SITELONG     |            |
# |              |                     |             |SRCLINK     |
# |              |ST                   |             |            |
# |SWCREATE      |SWCREATE             |SWCREATE     |            |
# |SWOWNER       |SWMODIFY             |SWOWNER      |            |
# |SWSERIAL      |SWOWNER              |SWSERIAL     |            |
# |              |SWSERIAL             |             |            |
# |TELESCOP      |TELESCOP             |TELESCOP     |            |
# |              |TIME-OBS             |             |            |
# |              |TIMESYS              |             |            |
# |TRAKTIME      |TRAKTIME             |TRAKTIME     |            |
# |              |USERNAME             |             |            |
# |              |UT                   |             |            |
# |XBINNING      |XBINNING             |XBINNING     |            |
# |XORGSUBF      |XORGSUBF             |XORGSUBF     |            |
# |XPIXSZ        |XPIXSZ               |XPIXSZ       |            |
# |YBINNING      |YBINNING             |YBINNING     |            |
# |YORGSUBF      |YORGSUBF             |YORGSUBF     |            |
# |YPIXSZ        |YPIXSZ               |YPIXSZ       |            |
# |              |iTelescope           |             |            |
# |              |iTelescopePlateScaleH|             |            |
# |              |iTelescopePlateScaleV|             |            |
# 

# ## Astrometry (aka Navigation)
# 
# Combining multiple images of the same source in the same filter requires that each image have a valid World Coordinate System solution, based on an astrometric solution to the stars in the field of view. The `AstroPhotography` package also refers to this process as image navigation. 
# 
# The process consists of finding the pixel locations brightest N stars in each image, then calling `Astrometry.net` to compute a WCS solution based on those star locations, and creating a new "navigated" fits file that combines the calibrated file and the WCS solution. By default 2-dimensional Gaussian profiles are also fitted to a representative number of stars per image to assess the Full Width at Half Maximum of the star images (in pixels), the results of which are shown graphically in a plot and numerically in a YaML file describing the stars found and fitted per image. At the end of the process a summary of the image (star) statistics is generated in CSV format that also incorporates the plate scale (arcseconds per pixel) found in the astrometric solution, which allows a human to look for outliers such as images where the seeing or tracking was especially bad.
# 
# In terms of pseudo-code
# ```python
# for cal-img in calibrated-images:
#     ap_find_stars cal-img --> source-list source-plot fwhm-plot quality-report.yml ds9-sources.reg
#     ap_astrometry cal-img source-list --> navigated-img
# ap_quality_summary (all quality-report.yml files) --> quality-summary.csv
# ```
# Using the (current but deprecated) bash script `navigate_all.sh` the output for a simple calibrated file is as follows:
# ```bash
# Processing ./20210306/cal-T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.fits at Sun Feb 11 04:17:23 PM EST 2024
#   Cleaning all files...
# removed 'SourceLists/srclist_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.fits'
# removed 'SourceLists/implot_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.png'
# removed 'SourceLists/ds9_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.reg'
# removed 'MetaData/qual_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.yaml'
# removed 'MetaData/srclog_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.log'
# removed 'NavigatedImages/nav_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.fits'
#   Performing star detection on ./20210306/cal-T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.fits
#     About the run the following command python3 /home/dks/git/AstroPhotography/AstroPhotography/scripts/ap_find_stars.py -l DEBUG -m 200 ./20210306/cal-T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.fits SourceLists/srclist_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.fits --retain_saturated --plotfile=SourceLists/implot_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.png --fwhm_plot=SourceLists/fwhmplot_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.png --quality_report=MetaData/qual_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.yaml --ds9=SourceLists/ds9_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.reg
#     Star detection completed successfully at Sun Feb 11 04:17:28 PM EST 2024
#       Log file written to MetaData/srclog_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.log
#   Performing astrometry on ./20210306/cal-T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.fits
#     About the run the following command:  python3 /home/dks/git/AstroPhotography/AstroPhotography/scripts/ap_astrometry.py -l DEBUG ./20210306/cal-T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.fits SourceLists/srclist_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.fits NavigatedImages/nav_T05-davestrickland-M81-20210306-235058-B-BIN1-W-300-001.fits
#     Astrometry completed successfully at Sun Feb 11 04:17:38 PM EST 2024
# ```

# ### Navigating a single image
# 
# Let's try detecting the stars in a single premium image, the positions of which can then be used with Astrometry.net to get an astrometric solution.

# In[10]:


# parameters used with source detection
p_loglevel  = 'DEBUG'                                                     # Loglevel
p_fitsimg   = str(flist[0])                                               # Input fits image
p_extnum    = 0                                                           # Extension number for data in input fits
p_fitstbl   = p_fitsimg.replace('Calibrated-iTelescope', 'srclist')       # Output fits table of srcs
p_search_fwhm       = 3.0                                                 # initial estimate at star fwhm in pixels
p_search_nsigma     = 7.0                                                 # min detection sigma above background
p_detector_bitdepth = 16                                                  # detector bitdepth (16 for most CCDs, ?? for CMOS, ?? for camera)
p_sat_frac          = 0.8                                                 # Fraction of full well at which we assume star saturated
p_max_sources = 200                                                       # Max number of sources to output or None
p_nosatmask   = True                                                      # Keep possibly saturated stars.
p_regfile     = p_fitsimg.replace('Calibrated-iTelescope', 'ds9').replace('fits', 'reg')       # Name of ds9 region file or None
p_plotfile    = p_fitsimg.replace('Calibrated-iTelescope', 'implot').replace('fits', 'png')    # None or name of output png image w sources.
p_qual_rprt   = p_fitsimg.replace('Calibrated-iTelescope', 'qual').replace('fits', 'yaml')     # None or name of output ascii report
p_quiet       = False                                                                          # Quiet mode suppresses STDOUT list printing
p_fwhm_plot   = p_fitsimg.replace('Calibrated-iTelescope', 'fwhmplot').replace('fits', 'png')  # PNG/JPG plot of PSF fits

# parameters associated with astrometric solution, if not already specified above.
# Astrometry.Net API key will be dealt with elsewhere
p_inpimg          = p_fitsimg
p_srclist         = p_fitstbl
p_src_extname     = 'AP_XYPOS'
p_outimg          = p_fitsimg.replace('Calibrated-iTelescope', 'navigated')
p_use_sip         = False
p_user_scale      = None
p_scale_err_ratio = None

# Regenerate files even if they exist, set p_clean = True
p_clean = False


# In[11]:


def does_file_exists(filename, verbose=False):
    """
    Returns True if the file name or path exists, false otherwise
    """
    if verbose and not Path(filename).exists():
        print(f"Cannot find {filename}. Not a valid path or file.")
        return False
    else:
        print(f"Found {filename}.")
    return True


# In[12]:


def find_stars_wrapper(p_loglevel, p_fitsimg, p_fitstbl, 
    p_extnum=0, p_search_fwhm=3.0, p_search_nsigma=7.0,
    p_detector_bitdepth=16, p_sat_frac=0.8, p_max_sources=200,
    p_nosatmask=True, p_plotfile=None, p_quiet=False,
    p_fwhm_plot=None, p_qual_rprt=None, p_regfile=None):
    """
    Wrapper for finding stars using ApFindStars within a python session,
    instead of using the command line ap_find_stars.py script.
    """

    # Perform initial source detection using default parameters.
    find_stars = ap.ApFindStars(p_fitsimg, p_extnum, p_search_fwhm,
        p_search_nsigma, p_detector_bitdepth, 
        p_max_sources, p_nosatmask, p_sat_frac, p_loglevel,
        p_plotfile, p_quiet)
    
    # Measure 2-Gaussian FWHM for select stars, get average over x and y
    (p_new_fwhm, p_madstd_fwhm, p_npts) = find_stars.measure_fwhm(p_fwhm_plot, 'both')
    
    # Refine source detection
    find_stars.source_search(p_new_fwhm, p_search_nsigma)
    
    # Re-run photometry
    find_stars.aperture_photometry()
    
    # As the source searching and photometry was redone, we should redo
    # the plotting.
    if p_plotfile is not None:
        find_stars.plot_image(p_plotfile)
    
    # Write optional quality report
    if p_qual_rprt is not None:
        find_stars.write_quality_report(p_qual_rprt)

    # Write optional ds9 format region file
    if p_regfile is not None:
        find_stars.write_ds9_region_file(p_regfile)
    
    # Write final sourcelist with photometry.
    find_stars.write_source_list(p_fitstbl)

    return


# In[13]:


find_stars_wrapper(p_loglevel, p_fitsimg, p_fitstbl, 
    p_extnum, p_search_fwhm, p_search_nsigma,
    p_detector_bitdepth, p_sat_frac, p_max_sources,
    p_nosatmask, p_plotfile, p_quiet,
    p_fwhm_plot, p_qual_rprt, p_regfile)


# Rather surprisingly this ran without error, despite the limited set of header keywords. The logging did note that "FITS keywords missing from image: ['EXPOSURE', 'DATE-OBS', 'OBJECT', 'OBJNAME', 'TELESCOP', 'INSTRUME', 'CCD-TEMP', 'RA-OBJ', 'DEC-OBJ', 'XPIXSZ', 'YPIXSZ', 'EGAIN', 'LAT-OBS', 'LONG-OBS', 'ALT-OBS', 'AIRMASS']"
# 
# #### Missing Keywords
# 
# Note that a warning is emitted stating several keywords, in particular `APRX_XWD APRX_YHG APRX_XPS APRX_YPS`, will not be written to the quality report. These represent the approximate angular width and height of the image in decimal degrees, and the approximate plate scale in X and Y. These could not be computed because some of the required telescope metadata was not present in the input image header. In this case, the CCD detector XPIXSZ and YPIXSZ were missing, so the approximate plate scale could not be computed. Without that information the approximate field of view (angular width & height) of the image could not be computed.
# 
# #### Is That Large Background Real?
# 
# Another item to note is a message early in the output: `Source-masked image stats: mean=489.139, median=488.849, stddev=25.656`. The average signal per pixel after excluding the region around the initially detected point sources is 489 +/- 26 per pixel. The units are not listed in the header, so they could be ADU, or could have been converted into electrons. One concern is that the numerical values iTelescope images often include an artificial positive `OFFSET` or `PEDESTAL`. This is listed in the raw FITS files, and if not removed acts as a large increase in sky brightness. It is not clear whether this is present in these premium images. 
# 
# We're not going to try to solve that question, or subtract the background, in this notebook.

# #### Astrometry
# 
# To attempt astrometry using `Astrometry.net` we need an **API key**.  If you have set up the astroquery config file with one already (see `doc/iTelescope_processing.md`) you can hope that it is picked up, or you can eneter one in the following prompt.

# In[14]:


msg = 'Enter your Astrometry.net API key, or hit return to use one preconfigured with astroquery: '
p_astnetkey = input(msg).strip() or None
if p_astnetkey is None:
    print("Warning: Assuming that 'api_key' is specified in your ~/.astropy/config/astroquery.cfg file.")
    print("  If this is not the case processing will fail.")


# In[15]:


ap_astrom = ap.ApAstrometry(p_inpimg, 
    p_srclist,
    p_outimg, 
    inp_img_extnum=p_extnum,
    srclist_extname=p_src_extname,
    astnet_key=p_astnetkey,
    use_sip=p_use_sip,
    user_scale=p_user_scale,
    scale_err_ratio=p_scale_err_ratio,
    loglevel=p_loglevel)
p_status = ap_astrom.status()
print(f'ApAstrometry return status: {p_status}')
if p_status == ap.ApAstrometry.NOMINAL:
    print('  Astrometric solution succeeded.')
elif p_status == ap.ApAstrometry.INPUT_ERROR:
    print('  Error, incorrect or missing input.')
elif p_status == ap.ApAstrometry.NO_SOLUTION:
    print('  Error, Astrometry.net failed to find a solution and/or timed out.')
else:
    print('  Error, unexpected error code returned by Astrometry.net')


# Astrometry.net can plate solve an image (or a set of star coordinates drawn from an image, as in this case) without any other information. This can be time consuming, as it has to try solutions of a wider variety of scales and pointings (image covering the whole sky down to images covering a few arcminutes, centered anyway on the sky). The solution can fail to return an answer, if the process run time exceeds limits set by the Astrometry.net owners.
# 
# Providing it a hint at the approximate pointing and/or approximate angular size of the image can considerably speed up plate solution. The hint does not need to be particularly accurate, so even a rough guess is better than nothing. Note that in this example ApAstrometry picked up on the `RA`, `DEC` in the image header, but defaulted to a guess of 8 degrees for the image size. The value of 8 degrees is a reasonable median between camera images covering tens of degrees and amateur telescopes covering 30 arcminutes to a few degrees. 
# 
# After the solution we can see that the image field of view is roughly 0.5 degrees on a side, so providing a hint of 8 degrees did not compromise the solution. You can also see that the approximate image center `(210.80166666666665, 54.349444444444444)`, does differ from the astrometric solution `(210.754804653, 54.4174059435)` by approximately 0.05 degrees (3 arcminutes) in RA and Dec.

# #### A brief aside about the CDi_j matrix
# 
# The WCS returned by Astrometry.net (mostly) follows the official FITS standard. In particular a `CD` matrix is used instead of the older `CDELT1`, `CDELT2`, `CROTA2`. The older keywords are more informative when you want an idea of the traditional plate scale (arcseconds per pixel). The conversion back from a simple `CDi_j` matrix into CDELT and CROTA is given in https://lweb.cfa.harvard.edu/~jzhao/SMA-FITS-CASA/docs/wcs88.pdf. Note that the conversion assumes an idealized case of square pixels at the tangent point.
# 
# If you are only interested in the magnitude `CDELT1`, `CDELT2`, and not their signs or the rotation angle `CROTA2`, then the conversion is as simple as:
# ```
# |CDELT1| = sqrt( CD1_1^2 + CD2_1^2 )    # Note: not CD1_1^2 and CD1_2^2
# 
# |CDELT2| = sqrt( CD1_2^2 + CD2_2^2 )    # Note: not CD2_1^2 and CD2_2^2
# ``` 

# # Processing all the images
# 
# We'll now run though processing all the images, but still doing most of the coding work by hand.

# In[16]:


p_dir = pathlib.Path(r'./')
flist = list(sorted(p_dir.glob(file_pattern)))
print(f'{len(flist)} files match pattern "{file_pattern}" in current directory.')
for fpath in flist:
    print(f'  {fpath.name}')


# In[17]:


keys = ['naxis1', 'naxis2', 'imagetyp', 'object', 'filter', 'exposure']
ifc_cal = ImageFileCollection(p_dir, keywords=keys, glob_include=file_pattern) # no glob_excludes at this time...
print(ifc_cal.summary)


# **Note:** here we see one of the issues with the lack of a strict FITS standard. The common fits keyword `EXPOSURE` is not defined in these headers, and in fact the premium images use `EXPTIME`. In fact, both `EXPOSURE` and `EXPTIME` are valid.
# 
# For a list of commonly used (in professional settings) FITS keywords see https://heasarc.gsfc.nasa.gov/docs/fcg/common_dict.html. Note that the [official standard](https://fits.gsfc.nasa.gov/fits_standard.html) is defined elsewhere. This [HEASARC webpage has a subset of the standard keywords](https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html) but note that it does not include the recommendations from the FITS papers.
# 
# This is not a huge issue for getting an Astrometric solution, where each image is processed separately. However, when later combining and mosaicing images we are likely to want to select on the exposure time of the observation in order to maximize the final signal to noise ratio.

# In[18]:


# define a structure to hold the processing status of the images
# Note: this is an astropy table, so it has that annoying string length must be specified
# ir you must use `object` limitation

arr_size  = len(ifc_cal.summary)
arr_dtype = [('filename', object), ('find_stars_status', object), ('find_stars_time', float), ('astrometry_status', object), ('astrometry_time', float)]

processing_results = Table(data=np.empty(arr_size, dtype=arr_dtype))
processing_results['find_stars_time'].info.format = '7.3f'
processing_results['astrometry_time'].info.format = '7.3f'

#processing_results = Table(names=('filename', 'find_stars_status', 'find_stars_time', 'astrometry_status', 'astrometry_time'), dtype=(str, str, float, str, float))
print(processing_results)


# In[19]:


# Assorted processing settings# parameters used with source detection
p_loglevel  = 'DEBUG'                                                     # Loglevel
p_fitsimg   = str(flist[0])                                               # Input fits image
p_extnum    = 0                                                           # Extension number for data in input fits
p_search_fwhm       = 3.0                                                 # initial estimate at star fwhm in pixels
p_search_nsigma     = 7.0                                                 # min detection sigma above background
p_detector_bitdepth = 16                                                  # detector bitdepth (16 for most CCDs, ?? for CMOS, ?? for camera)
p_sat_frac          = 0.8                                                 # Fraction of full well at which we assume star saturated
p_max_sources = 200                                                       # Max number of sources to output or None
p_nosatmask   = True                                                      # Keep possibly saturated stars.
p_quiet       = False                                                                          # Quiet mode suppresses STDOUT list printing

# parameters associated with astrometric solution, if not already specified above.
# Astrometry.Net API key will be dealt with elsewhere
p_inpimg          = p_fitsimg
p_srclist         = p_fitstbl
p_src_extname     = 'AP_XYPOS'
p_use_sip         = False
p_user_scale      = None
p_scale_err_ratio = None


# In[20]:


# Settings for output file names and output file paths
# - This is difficult to fully automate without some assumptions about the input file names
#   - I assume that all input file have the same file name prefix, e.g. Calibrated, or cal
#   - I assume that all input files have the same file extension, e.g. fits or (sadness) fit
# - I prefer to place output files of a certain type in a subdirectories. If you don't want all
#   the diagnostic plots then it may just be simpler to not deal with subdirectories.
# - If dir is not None then the various outputs files will be written to directories with the specified
#   path relative to the **current** directory. The directory will be created if it not already
#   present.

name_conv_dict = {'srclist':  {'replace': 'Calibrated-iTelescope', 'with': 'srclist',   'extension': None,   'dir': 'SourceLists'},
                  'regfile':  {'replace': 'Calibrated-iTelescope', 'with': 'ds9',       'extension': 'reg',  'dir': 'SourceLists'},
                  'plotfile': {'replace': 'Calibrated-iTelescope', 'with': 'implot',    'extension': 'png',  'dir': 'SourceLists'},
                  'fwhmfile': {'replace': 'Calibrated-iTelescope', 'with': 'fwhmplot',  'extension': 'png',  'dir': 'SourceLists'},
                  'qualfile': {'replace': 'Calibrated-iTelescope', 'with': 'qual',      'extension': 'yaml', 'dir': 'MetaData'},
                  'navfile':  {'replace': 'Calibrated-iTelescope', 'with': 'navigated', 'extension': None,   'dir': 'NavigatedImages'}}

# Settings for reprocessing. 
# If p_clean == True then existing output files will be erased
# If p_clean == False then if a file with the expected sourcelist output file name exists
#   (for find stars) or navigated image output name exists (for astrometry) then that stage
#   of processing is assumed completed and will not be done again. This is a considerable
#   time saver when rerunning the pipeline to fix one or two files that failed previously.
p_clean = False


# In[21]:


def mk_output_file_path(inpfile, prefix_replace_this, prefix_with_this, new_file_extension=None, sub_dir=None, mkdir=False, verbose=True):
    """
    Given an input file name, modify the file prefix, optionally modifying the file
    extension and adding a subdirectory to the front of the returned path.

    :param inpfile: String containing the input file name we want to use as a 
      basis for further file names, e.g. `Calibrated-iTelescope-M101_Supernova_2023ixf-180s-Lum-1.fits`
    :param prefix_replace_this: A part of the input file name that will be replaced with prefix_with_this,
      for example `Calibrated-iTelescope`. It is up to to the user to ensure that this substring
      is present in the input file name.
    :param prefix_with_this: What to replace prefix_replace_this with. For example `srclist`.
    :param new_file_extension: If not None then replace the `fits` file extension with this.
      Do not include the `.`.
    :param sub_dir: 
    :param mkdir: 
    :param verbose:

    Returns a string representing a new file name, for use by other programs. 
    """

    output_file_path = None
    if verbose:
        print(f'Input file name:  {inpfile}')
    
    if sub_dir is not None:
        if not os.access(sub_dir, os.F_OK):
            if mkdir:
                try:
                    os.mkdir(sub_dir)
                except:
                    print(f'Error, mkdir failed trying to create {sub_dir}')
                    raise
                if verbose:
                    print(f'Successfully created subdirectory {sub_dir}')
            else:
                raise RuntimeError(f'Subdirectory {sub_dir} cannot be accessed, but mkdir=False')
        else:
            if verbose:
                 print(f'Subdirectory {sub_dir} already exists.')

    tmp_str = inpfile.replace(prefix_replace_this, prefix_with_this)
    if new_file_extension is not None:
        tmp_str = tmp_str.replace('fits', new_file_extension)
    if sub_dir is not None:
        output_file_path = f'{sub_dir}/{tmp_str}'
    else:
        output_file_path = tmp_str

    if verbose:
        print(f'Output file name: {output_file_path}')
    if inpfile == output_file_path:
        print(f'Warning: input file name ({inpfile}) matches output file name ({output_file_path}).')
                
    return output_file_path

def mk_all_file_names(inpfile, name_conv_dict, mkdir=True, verbose=False):
    """
    Given an input file name and a conversion dictionary, generate all the file names that might
    be used to run star finding and astrometry on that image.

    Returns f_srclist, f_regfile, f_plotfile, f_fwhmfile, f_qualfile, f_navfile
    """
    oname = ['srclist', 'regfile', 'plotfile', 'fwhmfile', 'qualfile', 'navfile']
    ofile = {}
    for ftype in name_conv_dict.keys():
        ofile[ftype] = mk_output_file_path(inpfile,
            name_conv_dict[ftype]['replace'],
            name_conv_dict[ftype]['with'],
            name_conv_dict[ftype]['extension'],
            name_conv_dict[ftype]['dir'],
            mkdir,
            verbose)
    
    return ofile['srclist'], ofile['regfile'], ofile['plotfile'], ofile['fwhmfile'], ofile['qualfile'], ofile['navfile']


# In[22]:


# test the function
mkdir   = True
verbose = True
for ftype in name_conv_dict.keys():
    print(f'For file type {ftype} ---------------------------------------------------')
    ofile   = mk_output_file_path(ifc_cal.summary['file'][0],
                    name_conv_dict[ftype]['replace'],
                    name_conv_dict[ftype]['with'],
                    name_conv_dict[ftype]['extension'],
                    name_conv_dict[ftype]['dir'],
                    mkdir,
                    verbose)

print('\nAll the file names:')
f_srclist, f_regfile, f_plotfile, f_fwhmplot, f_qualfile, f_navfile = mk_all_file_names(ifc_cal.summary['file'][0],
                                                                                       name_conv_dict, mkdir, False)
print(f'  srclist={f_srclist}\n  regfile={f_regfile}\n  plotfile={f_plotfile}\n  fwhmplot={f_fwhmplot}\n  qualfile={f_qualfile}\n  navfile={f_navfile}')


# In[23]:


def check_file_exists(filename):
    """Return true if the specified file path exists"""
    if not os.path.isfile(filename):
        return False
    return True


# In[24]:


# iterate over files in the image file collection
#
# TODO work out how to direct output or log to a file...
p_clean     = False
p_loglevel  = 'INFO'                                            
p_quiet     = True
import time
idx = 0
proc_tstart = time.perf_counter()
for hdu, fname in ifc_cal.hdus(return_fname=True):
    print(80*'-')
    print(f'Processing input file {fname}')

    f_srclist, f_regfile, f_plotfile, f_fwhmplot, f_qualfile, f_navfile = mk_all_file_names(fname,
        name_conv_dict, mkdir, False)

    # Determine whether to perform star detection, based on p_clean and presence of
    does_srclist_exist = check_file_exists(f_srclist)
    if p_clean or not does_srclist_exist:    
        print(f'\n  Finding stars, output source list: {f_srclist}')
        fs_tstart = time.perf_counter()
        status = 'pass'
        try:
            find_stars_wrapper(p_loglevel, fname, f_srclist, 
                p_extnum, p_search_fwhm, p_search_nsigma,
                p_detector_bitdepth, p_sat_frac, p_max_sources,
                p_nosatmask, f_plotfile, p_quiet,
                f_fwhmplot, f_qualfile, f_regfile)
        except:
            print(f'  Error, caught exception when processing {fname}')
            status = 'fail'
        
        fs_tend     = time.perf_counter()
        fs_telapsed = fs_tend - fs_tstart # seconds
        fs_status   = status
    else:
        print(f'\n  Skipping star detection because {f_srclist} exists and clean={p_clean}')
        fs_status = 'skipped_exists'
        fs_telapsed = 0

    # Determine whether to perform star detection, based on p_clean and presence of
    does_navfile_exist = check_file_exists(f_navfile)
    if does_navfile_exist:
        if p_clean or not does_navfile_exist: 
            ast_status = 'pass'
            ast_tstart = time.perf_counter()
            print(f'\n  Performing astrometry, output navigated image: {f_navfile}')
            try:
    
                ap_astrom = ap.ApAstrometry(fname, 
                    f_srclist,
                    f_navfile, 
                    inp_img_extnum=p_extnum,
                    srclist_extname=p_src_extname,
                    astnet_key=p_astnetkey,
                    use_sip=p_use_sip,
                    user_scale=p_user_scale,
                    scale_err_ratio=p_scale_err_ratio,
                    loglevel=p_loglevel)
                p_status = ap_astrom.status()
                print(f'  ApAstrometry return status: {ast_status}')
                if p_status == ap.ApAstrometry.NOMINAL:
                    ast_status = 'pass'
                elif p_status == ap.ApAstrometry.INPUT_ERROR:
                    ast_status = 'input_error'
                else:
                    ast_status = 'fail'
            except:
                print(f'  Error, caught exception when processing {fname}')
                ast_status = 'exception'
    
            ast_tend     = time.perf_counter()
            ast_telapsed = ast_tend - ast_tstart # seconds
            ast_status   = status
        else:
            print(f'\n  Skipping astrometry because {f_navfile} exists and clean={p_clean}')
            ast_status = 'skipped_exists'
            ast_telapsed = 0
    else:
        # No source list
        print(f'\n  Skipping astrometry because {f_srclist} does not exist.')
        ast_status = 'skipped_no_srclist'
        ast_telapsed = 0

    # Fill in run info for this file
    result_tuple = (fname, fs_status, fs_telapsed, ast_status, ast_telapsed)
    processing_results[idx] = result_tuple

    # update idx
    idx += 1
    
proc_tend = time.perf_counter()
proc_telapsed = proc_tend - proc_tstart
print(f'Finished processing {idx+1} files in {proc_telapsed:.3f} seconds.')


# In[25]:


#print(processing_results)
processing_results.pprint(max_lines=-1, max_width=-1)


# ## Image Resampling and Mosaicing
# 
# Assuming that all the input images have been successfully navigated, i.e. have valid WCS solutions, we can move on to combining the images.
# 
# It is worth checking that the navigated images do have good solutions using SAOImage `ds9`. For example:
# ```bash
# ds9 -align yes -asinh -zoom to fit -cmap viridis -blink interval 2 -blink yes -tile no -single \
#     -scale mode 99.5 NavigatedImages/navigated*.fits
# ```
# 
# That should have loaded the images into separate frames, but only shown one at a time with the program blinking through them. However it seems to always tile them with multiple inputs. You need to perform the following actions in the ds9 GUI.
# 
# - resize the ds9 window to make it larger by dragging on the window edge
# - Click through the menu items `Frame`, then `Single Frame`
# - Click through the menu items `Frame`, then `Match`, then `Frame`, then `WCS`
# - Click through the menu items `Frame`, then `Blink Frames`
# 
# The same stars/features in each separate image should appear to be in the same place on the screen, even as the brightness and noise levels change between frames. 

# In[26]:


# Make sure we're still in the original directory, then move into the directory with the navigated images
os.chdir(wdir) 
ftype   = 'navfile'
fdir    =  name_conv_dict[ftype]['dir']
fprefix =  name_conv_dict[ftype]['with']
fext    =  name_conv_dict[ftype]['extension']
if fext is None:
    fext='.fits'
try:
    os.chdir(fdir)
    print(f'Successfully changed directory to {os.getcwd()}')
except:
    print('Error, failed to change directory to expected navfile directory {fdir}')
    print(f'  Current directory: {os.getcwd()}')


# ### Caveats
# 
# **Note:** The following code assumes that the navigated images in the current directory correspond to **a single target**.
# 
# All images in the directory will be resampled to a single pointing and platescale. By default this will be selected based on the first navigated image in the image file collection.
# 
# We would not necessarily want to resample images of the same `filter` but differing `exposure` together as their noise characteristics might be quite distinct, and the resampling method used might not be able to correctly weight them. Similarly cases where the seeing is significantly poorer than average, as determined by the fits performed by `ApMeasureStars` and reported in the quality files, should be excluded from being combined.
# 
# Such ideals run afoul of cases where the fits files don't have those keywords or are using different keywords. For now I don't have a work around for such cases and for automatically incorporating the seeing/FWHM measurements.

# In[27]:


nav_file_pattern=f'{fprefix}*{fext}'
p_dir = pathlib.Path(r'./')
keys = ['naxis1', 'naxis2', 'imagetyp', 'object', 'filter', 'exposure']
ifc_nav = ImageFileCollection(p_dir, keywords=keys, glob_include=nav_file_pattern) # no glob_excludes at this time...
print(ifc_nav.summary)


# In[28]:


from astropy.table import unique
uniq_filter_exposure = unique(ifc_nav.summary, keys=['filter', 'exposure'], keep='first', silent=True)
uniq_filter_exposure.pprint(max_lines=-1, max_width=-1)


# In[29]:


filter_list = []
for filter in uniq_filter_exposure['filter']:
    filter_list.append(filter)
print(f'There are {len(filter_list)} unique filters among {len(ifc_nav.summary)} files: {filter_list}')


# ### What Pointing and Plate Scale Should We Resample To?
# 
# There are any number of possible choices we can make regarding the pointing and plate scale that all the images will be resampled to. A few possible choices are shown:
# 
# 1. Resample all images to match the WCS of another image, for example the first image in the list of navigated images, irrespective of how uniform the platescale is or whether the image is North-up, East-left aligned.
# 2. Require the user to explicitly specify the RA and Dec of the image center, the X (RA) and Y (Dec) `CDELT*` values, and the number of output image rows and columns.
# 3. Require the user to specify a single plate scale (arcsec/pixel) and a RA/Dec or target identifiers, and calculate the values needed for #2.
# 
# There are advantages and disadvantages to each of these approaches. Method #1 may be superior for a final pretty display image where alignment is unimportant, while #2 and #3 are likely better for scientific use of the combined images. Method #2 requires user knowledge of how to specify a WCS and the likely platescale of the images, and is really only good for experts.
# 
# For now AstroPhotography will concentrate on supporting methods #1 and #3, starting with method #1.
# 
# **Note:** In image resampling the question of whether the input is well sampled is important. We do not have the space to cover it here, except to note that in general it is better to have pixels in the final (resampled) image be square and as large or indeed slightly larger than the pixels in the input images.
# 
# **Note:** For technical reasons most resampling methods assume the input data pixels represent surface brightness. If the input image has a large field of view it is likely distorted, in particular with images from DLSR or Mirrorless cameras, and the sky area covered by each pixel is not uniform. `Astrometry.net` allows fitting of `SIP` distortion coefficients. However, correctly calculation surface brightness in such cases is difficult (this package does not currently support it, and I don't know of any that do). Furthermore many of the resampling methods available don't correctly handle the non-uniform spatial aspect of SIP distortions either.
# 
# #### Other methods: resample_all.sh shell script
# 
# The bash script `resample_all.sh` does perform resampling using Astromatic's `swarp` resampler. Note that the script needs to be hard-wired for each set of images by specifying the RA, Dec of the image center, the plate scale, and the number of pixels along each axis.
# 
# ```bash
# cd /old_lnx/home/dks/Downloads/iTelescopeScratch/Plan-20-M101SN
# ~/git/AstroPhotography/AstroPhotography/scripts/resample_all.sh M101sn median clean verbose
# ... [lots of output omitted]...
# ----------------------------------------------------------------------
# Processing Lum filter images:
#   Using existing temporary file directory TmpSwarp
#   Found 3 navigated fits files.
#   File navigated-M101_Supernova_2023ixf-180s-Lum-3.fits EXPOSURE  180.000 FSCALE 0.005556 
#   File navigated-M101_Supernova_2023ixf-180s-Lum-1.fits EXPOSURE  180.000 FSCALE 0.005556 
#   File navigated-M101_Supernova_2023ixf-180s-Lum-2.fits EXPOSURE  180.000 FSCALE 0.005556 
#   Flux scalings: 0.005556,0.005556,0.005556
#   Total exposure in Lum filter: 9 minutes.
# 
#   Beginning swarp for filter Lum at Wed Feb 28 09:25:14 PM EST 2024
# Running the following command:  swarp .././NavigatedImages/navigated-M101_Supernova_2023ixf-180s-Lum-3.fits .././NavigatedImages/navigated-M101_Supernova_2023ixf-180s-Lum-1.fits .././NavigatedImages/navigated-M101_Supernova_2023ixf-180s-Lum-2.fits -SUBTRACT_BACK N -COMBINE_TYPE MEDIAN -FSCALASTRO_TYPE VARIABLE -VERBOSE_TYPE FULL -NOPENFILES_MAX p_nopenmax -GAIN_KEYWORD EGAIN -FSCALE_DEFAULT 0.005556,0.005556,0.005556 -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE 0.62 -CENTER_TYPE MANUAL -CENTER 210.8022671,54.3489500 -RESAMPLING_TYPE LANCZOS3 -OVERSAMPLING 4 -IMAGE_SIZE 4096,4096 -PROJECTION_TYPE TAN -IMAGEOUT_NAME M101_Supernova_2023ixf_Lum_4096x4096_resamp_median.fits -WRITE_FILEINFO Y -WRITE_XML N -DELETE_TMPFILES Y -RESAMPLE_DIR TmpSwarp -WEIGHTOUT_NAME M101_Supernova_2023ixf_Lum_4096x4096_weights_median.fits
#   Generated resampled co-added image M101_Supernova_2023ixf_Lum_4096x4096_resamp_median.fits
#   Generated co-added weights image M101_Supernova_2023ixf_Lum_4096x4096_weights_median.fits
#   Swarp for filter Lum completed at Wed Feb 28 09:25:19 PM EST 2024, took 5.099 seconds.
# 
# ----------------------------------------------------------------------
# Processing Luminance filter images:
#   Using existing temporary file directory TmpSwarp
#   Found 0 navigated fits files.
#   Moving on to next filter...
# --------------------------------------------------------------------------
# Run summary:
# Filter      Resampled Output                                                   Ninput   TexpMin  SwrpTSec  Status
# Red         M101_Supernova_2023ixf_Red_4096x4096_resamp_median.fits                 3     15.00      5.20       0
# Green       M101_Supernova_2023ixf_Green_4096x4096_resamp_median.fits               3     15.00      5.10       0
# Blue        M101_Supernova_2023ixf_Blue_4096x4096_resamp_median.fits                3     15.00      5.09       0
# Ha          None                                                                    0      0.00      0.00       0
# OIII        None                                                                    0      0.00      0.00       0
# SII         None                                                                    0      0.00      0.00       0
# B           None                                                                    0      0.00      0.00       0
# V           None                                                                    0      0.00      0.00       0
# I           None                                                                    0      0.00      0.00       0
# Clear       None                                                                    0      0.00      0.00       0
# Lum         M101_Supernova_2023ixf_Lum_4096x4096_resamp_median.fits                 3      9.00      5.10       0
# Luminance   None                                                                    0      0.00      0.00       0
# Master log file: /old_lnx/home/dks/Downloads/iTelescopeScratch/Plan-20-M101SN/tmp_resamplecd_all.log
# 
# --------------------------------------------------------------------------
#   Temporary files occupy 1 MB in /old_lnx/home/dks/Downloads/iTelescopeScratch/Plan-20-M101SN/Resampled/TmpSwarp
#   Commands were logged to /old_lnx/home/dks/Downloads/iTelescopeScratch/Plan-20-M101SN/tmp_resamplecd_all.log
#   Resampled images were written to ./Resampled
# resample_all.sh finished at Sun Feb 25 08:33:13 AM EST 2024, run time 23.638 seconds.
# ```

# #### Resampling To Match Another Image
# 
# This is the simplest conceptual method, but does allow some flexibility in that the other image need not be one of the navigated images. For example we could resample to match the WCS of a Chandra X-ray image, which would be ideal for scientific multi-waveband analysis.
# 
# First let's specify the file we want to match...

# In[30]:


f_default_match = ifc_nav.summary['file'][0]
msg = f'Name of file to match WCS of in resampled images (hit return for default {f_default_match})'
f_to_match = input(msg).strip() or f_default_match
print(f'WCS of resampled files will match {f_to_match}')


# In[31]:


from astropy import wcs
import math
def load_wcs_from_file(filename, verbose=False):
    """
    Modified WCS example based on astropy documentation

    https://docs.astropy.org/en/stable/wcs/loading_from_fits.html

    Limitations: Assumes 2-dimensional image
    """
    w          = None
    pix_scales = None # in deg
    pix_area   = None # in deg^2
    im_scales  = None # in deg
    
    # Load the FITS hdulist using astropy.io.fits
    with fits.open(filename) as hdulist:

        # Parse the WCS keywords in the primary HDU
        w = wcs.WCS(hdulist[0].header)
    
        # Print out the "name" of the WCS, as defined in the FITS header
        if verbose:
            print(w.wcs.name)
    
            # Print out all of the settings that were parsed from the header
            w.wcs.print_contents()
    
        # Pixel scale at the CRPIX pixel location. Note, ideal, ignores distortions
        pix_scales = wcs.utils.proj_plane_pixel_scales(w)
    
        # Pixel area at the CRPIX pixel location. Again, ideal, ignores distortions
        pix_area = wcs.utils.proj_plane_pixel_area(w.celestial)

        if w.wcs.naxis != 2:
            print(f'WARNING, WCS from {filename} is {w.wcs.naxis}-dimensional, not 2-D as expected')
        
        if len(pix_scales) != len(w.array_shape):
            raise RuntimeError(f'Shape of pixel scales ({len(pix_scales)}) differs from shape of array ({len(w.array_shape)})')
        else:
            im_scales = pix_scales.copy()
            for idx in range(len(w.array_shape)):
                im_scales[idx] *= w.array_shape[idx]
    
        if verbose:
            print(f'Pixel angular scale ({w.wcs.cunit[0]}):    {pix_scales}')
            print(f'Pixel angular area ({w.wcs.cunit[0]}^2):   {pix_area}')
            print(f'NAXIS1={w.array_shape[0]}   NAXIS2={w.array_shape[1]}')
            print(f'Detector angular scale ({w.wcs.cunit[0]}): {im_scales}')
        
    
    return w, pix_scales, pix_area, im_scales


# In[32]:


def summarize_wcs(w):
    """
    Given an astropy WCS object, summarise the properties of the WCS header,
    printing to STDOUT.

    This assumes that the images are two dimensional, and that angles are in degrees.

    :param w: astropy.wcs.WCS instance
    """

    ipwcs = w
    
    dothead_list = []

    # Numbr of dimensions
    val_str = f"{'NAXIS':8s} = {ipwcs.wcs.naxis}"
    dothead_list.append( val_str )
    
    keywords = ["NAXIS", "CTYPE", "CRVAL", "CRPIX"]
    values = [ipwcs.array_shape, ipwcs.wcs.ctype, ipwcs.wcs.crval, ipwcs.wcs.crpix]
    for keyword, value in zip(keywords, values):
        for idx in range(ipwcs.naxis):
            kw_str  = f"{keyword}{1+idx}"
            if 'CTYPE' in keyword:
                # Wrap strings in quotes
                val_str = f"{kw_str:8s} = '{value[idx]}'"
            else:
                val_str = f"{kw_str:8s} = {value[idx]}"
            dothead_list.append( val_str )

    naxis  = ipwcs.wcs.naxis
    naxis1 = ipwcs.array_shape[0]
    naxis2 = ipwcs.array_shape[1]
    if naxis == 3:
        naxis3 = ipwcs.array_shape[2]
        print(f'Warning: summarize_wcs() not written for 2-D images, {naxis}-dimensional data')

    cdelt1 = None
    cdelt2 = None
    crot   = None
    
    if hasattr(ipwcs.wcs, "cd"):
        for irow in range(ipwcs.naxis):
            for jcol in range(ipwcs.naxis):
                kw_str = f'CD{irow+1}_{jcol+1}'
                val_str = f'{kw_str:8s} = {ipwcs.wcs.cd[irow, jcol]}'
                dothead_list.append( val_str )
        cd = ipwcs.wcs.cd
        # From https://lweb.cfa.harvard.edu/~jzhao/SMA-FITS-CASA/docs/wcs88.pdf
        cd11 = cd[0,0]
        cd12 = cd[0,1]
        cd21 = cd[1,0]
        cd22 = cd[1,1]
        cdelt1_mag = math.sqrt( cd12*cd12 + cd22*cd22 )
        cdelt2_mag = math.sqrt( cd12*cd12 + cd22*cd22 )
        # the sign of cdelt1.cdelt2 = sign of (cd11*cd22 - cd12*cd21)
        # if the RHS is negative cdelt1 is negative by convention, so cdelt2 is always positive
        tmpa =  cd11*cd22 - cd12*cd21
        cdelt1 = math.copysign(cdelt1_mag, tmpa)
        cdelt2 = cdelt2_mag
        sign   = math.copysign(1, tmpa)
        crot   = math.degrees( math.atan2( (sign*cd12), cd22) )
            
    elif hasattr(ipwcs.wcs, "pc"):
        for irow in range(ipwcs.naxis):
            for jcol in range(ipwcs.naxis):
                kw_str = f'PC{irow+1}_{jcol+1}'
                val_str = f'{kw_str:8s} = {ipwcs.wcs.pc[irow, jcol]}'
                dothead_list.append( val_str )
            kw_str = f'CDELT{1+irow}'
            val_str = 'f{kw_str:8s} = {ipwcs.wcs.cdelt[irow]}'
            dothead_list.append( val_str )
            # Think that CD1_* = CDELT1 * PC1_* and CD2_* = CDELT2 * PC2_*
            pc = ipwcs.wcs.pc
            cdelt_vec = ipwcs.wcs.cdelt
            cd = pc.copy()
            for irow in range(ipwcs.naxis):
                cd[irow] = cdelt_vec[irow] * pc[irow]
            # then as above
        raise RuntimeError('summarize_wcs needs to be updated to handle cases with a PC_ matrix and no CD_ matrix')
    else:
        # Assume we have a simple CDELT[12] case with CROTA
        for irow in range(ipwcs.naxis):
            kw_str = f'CDELT{1+irow}'
            val_str = 'f{kw_str:8s} = {ipwcs.wcs.cdelt[irow]}'
            dothead_list.append( val_str )
            kw_str = f'CROTA{1+irow}'
            val_str = 'f{kw_str:8s} = {ipwcs.wcs.crota[irow]}'
            dothead_list.append( val_str )
        cdelt1 = ipwcs.wcs.cdelt[0]
        cdelt2 = ipwcs.wcs.cdelt[1]
        crot   = ipwcs.wcs.crota[1] # CROTA2 is used, CROTA1 not used. Assume crota[0] == crota[1]

    cdelt1_as = cdelt1 * 3600.0      # arcseconds
    cdelt2_as = cdelt2 * 3600.0      # arcseconds
    ximgsz_am = cdelt1 * naxis1 * 60 # arcminutes
    yimgsz_am = cdelt2 * naxis2 * 60 # arcminutes
    dothead_list.append( f'Pixel size equivalent CDELT1     = {cdelt1_as:.3f} arcseconds' )
    dothead_list.append( f'Pixel size equivalent CDELT2     = {cdelt2_as:.3f} arcseconds' )
    dothead_list.append( f'Image X-axis angular size        = {ximgsz_am:.3f} arcminutes' )
    dothead_list.append( f'Pixel Y-axis angular size        = {yimgsz_am:.3f} arcminutes' )
    dothead_list.append( f'Image rotation equivalent CROTA2 = {crot:.3f} degrees' )
    
    dothead_list.append('END     ')
    dothead_str = '\n'.join(dothead_list)
    print(dothead_str)
    return


# In[33]:


# Test load_wcs_from_file
w, pix_scales, pix_area, im_scales = load_wcs_from_file(f_to_match, False)
pix_scales_as = pix_scales * 3600.0
root_area_as = math.sqrt(pix_area) * 3600.0
im_scales_am = im_scales * 60.0
print(f'Pixel scales (arcsec):                 {pix_scales_as}')
print(f'Image angular extent (arcmin):         {im_scales_am}')
print(f'Square pixel equivalent size (arcsec): {root_area_as}')

summarize_wcs(w)


# In[34]:


def make_dothead_from_file(input_fits_with_wcs, output_swarp_dothead, format='fits', verbose=False):
    """
    Create the .head format file that swarp expected based on the WCS header of an input files file.
    
    Internally this uses the method from https://docs.astropy.org/en/stable/_modules/astropy/wcs/wcs.html#WCS.printwcs
    
    Limitations: Assumes the WCS information is in the primary header  

    :param input_fits_with_wcs: Name of existing FITS file with WCS in primary HDU that we want to emulate.
    :param output_swarp_dothead:
    :param format: Either 'text' or 'fits'. If 'text' then an ASCII header will be created using the format
      specified in the Swarp manual. If 'fits' is specified, an empty FITS file consisting only of a primary
      HDU will be created.
    :param verbose: If True then diagnostic information will be written to stdout.
    """ 
    
    with fits.open(input_fits_with_wcs) as hdulist:
        # Parse the WCS keywords in the primary HDU
        ipwcs = wcs.WCS(hdulist[0].header)

        if verbose:
            print(f'WCS created from primary header of {input_fits_with_wcs}')
            print(ipwcs)
            print(80*"-")
            #print(ipwcs.wcs)
            #print(80*"-")

        if 'text' in format:
            if verbose:
                print('Generating an ASCII header following the format specified in the swarp documentation.')
        
            dothead_list = []
            keywords = ["NAXIS", "CTYPE", "CRVAL", "CRPIX"]
            values = [ipwcs.array_shape, ipwcs.wcs.ctype, ipwcs.wcs.crval, ipwcs.wcs.crpix]
            for keyword, value in zip(keywords, values):
                for idx in range(ipwcs.naxis):
                    kw_str  = f"{keyword}{1+idx}"
                    if 'CTYPE' in keyword:
                        # Wrap strings in quotes
                        val_str = f"{kw_str:8s} = '{value[idx]}'"
                    else:
                        val_str = f"{kw_str:8s} = {value[idx]}"
                    dothead_list.append( val_str )
    
            if hasattr(ipwcs.wcs, "pc"):
                for irow in range(ipwcs.naxis):
                    for jcol in range(ipwcs.naxis):
                        kw_str = f'PC{irow+1}_{jcol+1}'
                        val_str = f'{kw_str:8s} = {ipwcs.wcs.pc[irow, jcol]}'
                        dothead_list.append( val_str )
                    kw_str = f'CDELT{1+irow}'
                    val_str = 'f{kw_str:8s} = {pwcs.wcs.cdelt[irow]}'
                    dothead_list.append( val_str )
            elif hasattr(ipwcs.wcs, "cd"):
                for irow in range(ipwcs.naxis):
                    for jcol in range(ipwcs.naxis):
                        kw_str = f'CD{irow+1}_{jcol+1}'
                        val_str = f'{kw_str:8s} = {ipwcs.wcs.cd[irow, jcol]}'
                        dothead_list.append( val_str )
            
            dothead_list.append('END     ')
            dothead_str = '\n'.join(dothead_list)
    
            if verbose:
                print(f'--- dothead output from {input_fits_with_wcs} ---')
                print(dothead_str)
                print(f'--- about to write to {output_swarp_dothead} ---')
            
            with open(output_swarp_dothead, 'w') as ofile:
                ofile.write(dothead_str)
                if verbose:
                    print(f'Wrote Swarp ASCII-format .head file to {output_swarp_dothead}')
        elif 'fits' in format:
            if verbose:
                print('Generating a FITS format header consisting of a PrimaryHDU only.')
            # based on https://docs.astropy.org/en/stable/wcs/example_create_imaging.html
            hdr = ipwcs.to_header()

            # NAXIS = x, NAXISx values set from data, not by manipulating header 
            olddata = hdulist[0].data
            data = 0 * olddata.astype(int)

            hdu = fits.PrimaryHDU(header=hdr, data=data)
            hdu.writeto(output_swarp_dothead, overwrite=True, output_verify='ignore')
            if verbose:
                print(f'Wrote FITS header .head file to {output_swarp_dothead}')
        else:
            raise RuntimeError(f'Error, format ({format}) is not one of the allowed options: "text" "fits"')    
    return


# In[35]:


# Test make_dothead_from_file
ofilename = 'M101_Supernova_2023ixf_Lum_4096x4096_resamp_weighted.fits'
oheadname = ofilename.replace('.fits','.head')
make_dothead_from_file(f_default_match, oheadname, 'fits', True)


# ##### Resampling Using Python
# 
# We'll try to recreate the bash shell script resampling of the luminance images, although with the same header as an input image instead of a user-defined pointing. As a reminder the shell command the script used earlier was:
# ```
# swarp .././NavigatedImages/navigated-M101_Supernova_2023ixf-180s-Lum-3.fits .././NavigatedImages/navigated-M101_Supernova_2023ixf-180s-Lum-1.fits .././NavigatedImages/navigated-M101_Supernova_2023ixf-180s-Lum-2.fits -SUBTRACT_BACK N -COMBINE_TYPE MEDIAN -FSCALASTRO_TYPE VARIABLE -VERBOSE_TYPE FULL -NOPENFILES_MAX p_nopenmax -GAIN_KEYWORD EGAIN -FSCALE_DEFAULT 0.005556,0.005556,0.005556 -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE 0.62 -CENTER_TYPE MANUAL -CENTER 210.8022671,54.3489500 -RESAMPLING_TYPE LANCZOS3 -OVERSAMPLING 4 -IMAGE_SIZE 4096,4096 -PROJECTION_TYPE TAN -IMAGEOUT_NAME M101_Supernova_2023ixf_Lum_4096x4096_resamp_median.fits -WRITE_FILEINFO Y -WRITE_XML N -DELETE_TMPFILES Y -RESAMPLE_DIR TmpSwarp -WEIGHTOUT_NAME M101_Supernova_2023ixf_Lum_4096x4096_weights_median.fits
# ```
# 
# The `swarp` command line parameters we'll want to change are:
# 
# - `CELESTIAL_TYPE`: This needs to be set to native to copy the first input image.
# - `VERBOSE_TYPE`: By default `NORMAL`, options are `QUIET`, `LOG`, `NORMAL`, `FULL`.
# - `COMBINE_TYPE`: Median is good for a first look but will increase the variance. If possible use `WEIGHTED`
# - `PIXEL_SCALE_TYPE`: Various options, see `swarp` [manual](https://raw.githubusercontent.com/astromatic/swarp/legacy_doc/prevdoc/swarp.pdf). Maybe set to `MAX`
# - - Note this is ignored if a `.head` file is supplied
# - `CENTER_TYPE`: Instead of `MANUAL` we want `MOST` or `ALL`. The latter options also generate output images that are aligned North up, East to the left.
# - - Hopefully this is ignored if we specify a `.head`?
# - `CENTER` and `PIXEL_SCALE` wont be used if we change `PIXEL_SCALE_TYPE` and `CENTER_TYPE`
# - `IMAGE_SIZE`: Use 0 for automatic sizing based on other parameter values
# - `SUBTRACT_BACK`: Using `N` is appropriate if we have a lot of nebular emission, crowded fields, large galaxies, or have already performed background subtraction. You should really try with both 'Y' and 'N' to determine the effect of swarp's in-built background subtraction.
# - `RESAMPLE_DIR`: Directory where temporary files are written. This must exist.
# 
# In the case of the M101 SN premium images `SUBTRACT_BACK` should be `N`. The premium dataset is not a great test of the different center type options because the three input images very closely overlap.

# In[36]:


# ***You may want to disable this, as it uses assumed file names and directory paths*** 
# It also assumes a Unix like operating system!!!!
do_simple_test = False
do_better_test = True
do_better_test_center_type = 'FIRST' # MANUAL, MOST, ALL, FIRST (MANUAL and ALL are swarp types)
dothead_format = 'fits' # 'text' or 'fits

ofilename = 'M101_Supernova_2023ixf_Lum_4096x4096_resamp_weighted.fits'
pathlib.Path(ofilename).unlink(missing_ok=True)  # remove old file 
oheadname = ofilename.replace('.fits','.head')

if do_simple_test:
    # Make sure to delete any .head file
    pathlib.Path(oheadname).unlink(missing_ok=True)
    
    # test command, long form, based of the previous example from the resample_all.sh script
    test_cmd = 'swarp navigated-M101_Supernova_2023ixf-180s-Lum-1.fits navigated-M101_Supernova_2023ixf-180s-Lum-2.fits navigated-M101_Supernova_2023ixf-180s-Lum-3.fits -VERBOSE_TYPE FULL -SUBTRACT_BACK N -COMBINE_TYPE WEIGHTED -FSCALASTRO_TYPE VARIABLE -VERBOSE_TYPE FULL -GAIN_DEFAULT 1.0  -GAIN_KEYWORD EGAIN -FSCALE_DEFAULT 0.005556,0.005556,0.005556 -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE 0.62 -CENTER_TYPE MANUAL -CENTER 210.8022671,54.3489500 -RESAMPLING_TYPE LANCZOS3 -OVERSAMPLING 4 -IMAGE_SIZE 4096,4096 -PROJECTION_TYPE TAN -IMAGEOUT_NAME M101_Supernova_2023ixf_Lum_4096x4096_resamp_weighted.fits -WRITE_FILEINFO Y -WRITE_XML N -DELETE_TMPFILES Y -RESAMPLE_DIR ./ -WEIGHTOUT_NAME M101_Supernova_2023ixf_Lum_4096x4096_weights_weighted.fits'
    test_cmd_list = shlex.split(test_cmd)
    print(f'Command string to be passed to subprocess.run is {test_cmd_list}')
    fs_tstart     = time.perf_counter()
    result = subprocess.run(test_cmd_list, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    fs_tend     = time.perf_counter()
    fs_telapsed = fs_tend - fs_tstart # seconds
    print(f'  Process took {fs_telapsed:.3f} seconds, return code {result.returncode}')
    print(f'  Input args: {result.args}')
    print(f'  Stdout: {result.stdout}')
    print(f'  Stderr: {result.stderr}')

if do_better_test:
    # test command, split into related sub-sections
    
    test_cmd1 = 'swarp navigated-M101_Supernova_2023ixf-180s-Lum-1.fits navigated-M101_Supernova_2023ixf-180s-Lum-2.fits navigated-M101_Supernova_2023ixf-180s-Lum-3.fits -FSCALASTRO_TYPE VARIABLE -FSCALE_DEFAULT 0.005556,0.005556,0.005556' 
    test_cmd2 = '-VERBOSE_TYPE FULL -SUBTRACT_BACK N -COMBINE_TYPE WEIGHTED -GAIN_DEFAULT 1.0 -GAIN_KEYWORD EGAIN -RESAMPLING_TYPE LANCZOS3 -OVERSAMPLING 4 -PROJECTION_TYPE TAN'
    if 'MANUAL' in do_better_test_center_type:
        pathlib.Path(oheadname).unlink(missing_ok=True)
        test_cmd3 = '-PIXELSCALE_TYPE MANUAL -PIXEL_SCALE 0.62 -CENTER_TYPE MANUAL -CENTER 210.8022671,54.3489500 -IMAGE_SIZE 4096,4096'
    elif 'ALL' in do_better_test_center_type:
        test_cmd3 = '-PIXELSCALE_TYPE MAX -CENTER_TYPE ALL'
        pathlib.Path(oheadname).unlink(missing_ok=True)
    elif 'MOST' in do_better_test_center_type:
        test_cmd3 = '-PIXELSCALE_TYPE MAX -CENTER_TYPE MOST'
        pathlib.Path(oheadname).unlink(missing_ok=True)
    elif 'FIRST'  in do_better_test_center_type:
        # Need to generate a header based on the first file, with a name based on the output file name
        # The pixelscale and centertype should be set from the .head file...
        make_dothead_from_file(f_default_match, oheadname, dothead_format, True)
        test_cmd3 = ''
    else:
        raise RuntimeError(f'Error, center type {do_better_test_center_type} has not yet been implemented.')
        
    test_cmd4 = f'-IMAGEOUT_NAME {ofilename} -WRITE_FILEINFO Y -WRITE_XML N -DELETE_TMPFILES Y -RESAMPLE_DIR ./ -WEIGHTOUT_NAME M101_Supernova_2023ixf_Lum_4096x4096_weights_weighted.fits'
    test_cmd_list = shlex.split(test_cmd1) + shlex.split(test_cmd2) + shlex.split(test_cmd3) + shlex.split(test_cmd4)
    print(f'Command string to be passed to subprocess.run is {test_cmd_list}')
    fs_tstart     = time.perf_counter()
    try:
        result = subprocess.run(test_cmd_list, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        fs_tend     = time.perf_counter()
        fs_telapsed = fs_tend - fs_tstart # seconds
        print(f'  Success: process took {fs_telapsed:.3f} seconds, return code {result.returncode}')
        #print(f'  Input args: {result.args}')
        #print(f'  Stdout: {result.stdout}')
        #print(f'  Stderr: {result.stderr}')

    except subprocess.CalledProcessError as err:
        fs_tend     = time.perf_counter()
        fs_telapsed = fs_tend - fs_tstart # seconds
        print(f'  Error, process took {fs_telapsed:.3f} seconds, return code {err.returncode}')
        print(f'  Input args: {err.cmd}\n')
        print(f'  Stdout: {err.output}\n')
        print(f'  Stderr: {err.stderr}\n')
        print(f'  Command line equivalent command: {" ".join(err.cmd)}')
    


# ##### Issues with swarp's ASCII `.head` file format
# 
# Information retention:  Using swarp's `.head` format with the `END     ` does not seem to work. The following `.head` should work
# ```
# CTYPE1   = RA---TAN
# CTYPE2   = DEC--TAN
# CRVAL1   = 210.754804653
# CRVAL2   = 54.4174059435
# CRPIX1   = 1528.5
# CRPIX2   = 1528.5
# CD1_1    = -4.64681630389e-06
# CD1_2    = 0.000175223532872
# CD2_1    = -0.000175223532872
# CD2_2    = -4.64681630389e-06
# END     
# ```
# but we get the following error:
# ```
# ----- SWarp 2.38.0 started on 2024-03-13 at 09:04:47 with 16 threads
# 
# Examining input data ...
# Looking for navigated-M101_Supernova_2023ixf-180s-Lum-1.fits ...
# Looking for navigated-M101_Supernova_2023ixf-180s-Lum-2.fits ...
# Looking for navigated-M101_Supernova_2023ixf-180s-Lum-3.fits ...
# Creating NEW output image ...
# 
# > *FATAL ERROR*: Unknown FITS type in fitswrite()
# ```
# This is emitted by the `fitswrite` keyword writing code in `src/fits/fitsutil.c`. 
# 
# ##### Trying with an actual FITS header
# 
# Added an option to `make_dothead_from_file` to generate FITS format header. **This works**. Note that swarp can eat inwards from the edge...

# ##### Resampling All The Filters
# 
# Now that we have the method for resampling to match a given image we'll move onto to resampling all the bands (filters) we have.
# 
# The new challenge is to compute the appropriate scale factors and which header keywords need to be updated post-resampling (e.g. exposure related keywords, history, and so on). We'll also want some form of summary table.

# In[37]:


def get_exposure_time(hdr, verbose=False):
    """Return the first of EXPTIME, EXPOSURE, ONTIME, or LIVETIME from a FITS header"""
    keywords = ['EXPTIME', 'EXPOSURE', 'ONTIME', 'LIVETIME']
    exposure_time = None
    for key in keywords:
        if hdr.count(key) > 0:
            exposure_time = float( hdr[key] )
            if verbose:
                print(f'Keyword {key} found, setting exposure_time to {exposure_time}')
            break
        else:
            if verbose:
                print(f'Keyword {key} not found, continuing search...')
    return exposure_time


# In[38]:


do_better_test_center_type = "FIRST"  # MANUAL, MOST, ALL, FIRST (MANUAL and ALL are swarp types)
dothead_format = "fits"  # 'text' or 'fits
swarp_verbose  = False  # Echo swarp input string if True, echo .head for FIRST case
final_images   = []

for filter in filter_list:
    print(f"Processing images from filter  {filter}")

    # Names for output resampled image, net weights, and .head files
    ofilename = f"M101_Supernova_2023ixf_{filter}_resamp_weighted.fits"
    owgtsname = f"M101_Supernova_2023ixf_{filter}_weights_weighted.fits"
    pathlib.Path(ofilename).unlink(missing_ok=True)  # remove old file
    oheadname = ofilename.replace(".fits", ".head")
    pathlib.Path(oheadname).unlink(missing_ok=True)  # remove old file

    filtered_netexp = 0
    filtered_numimg = 0
    inp_weights = []
    inp_files = []

    for hdu, fname in ifc_nav.hdus(return_fname=True, filter=filter):
        inp_files.append(fname)
        texp = get_exposure_time(hdu.header)
        if texp is None:
            raise RuntimeError(f"Error, could not get exposure time from {fname}")
        fscale = 1.0 / texp
        inp_weights.append(f"{fscale:.7f}")
        filtered_netexp += texp
        filtered_numimg += 1

    if len(inp_files) == 0:
        print(f"  No files found for filter {filter}")
        continue
    else:
        fs_tstart = time.perf_counter()
        print(f"  {filter} {filtered_numimg} {filtered_netexp}")
        print(f"    inp_files:   {inp_files}")
        print(f"    inp_weights: {inp_weights}")

        file_str   = " ".join(inp_files)
        fscale_str = ",".join(inp_weights)
        test_cmd1  = f"swarp {file_str} -FSCALASTRO_TYPE VARIABLE -FSCALE_DEFAULT {fscale_str}"
        test_cmd2 = "-VERBOSE_TYPE FULL -SUBTRACT_BACK N -COMBINE_TYPE WEIGHTED -GAIN_DEFAULT 1.0 -GAIN_KEYWORD EGAIN -RESAMPLING_TYPE LANCZOS3 -OVERSAMPLING 4 -PROJECTION_TYPE TAN"
        if "MANUAL" in do_better_test_center_type:
            pathlib.Path(oheadname).unlink(missing_ok=True)
            test_cmd3 = "-PIXELSCALE_TYPE MANUAL -PIXEL_SCALE 0.62 -CENTER_TYPE MANUAL -CENTER 210.8022671,54.3489500 -IMAGE_SIZE 4096,4096"
        elif "ALL" in do_better_test_center_type:
            test_cmd3 = "-PIXELSCALE_TYPE MAX -CENTER_TYPE ALL"
            pathlib.Path(oheadname).unlink(missing_ok=True)
        elif "MOST" in do_better_test_center_type:
            test_cmd3 = "-PIXELSCALE_TYPE MAX -CENTER_TYPE MOST"
            pathlib.Path(oheadname).unlink(missing_ok=True)
        elif "FIRST" in do_better_test_center_type:
            # Need to generate a header based on the first file, with a name based on the output file name
            # The pixelscale and centertype should be set from the .head file...
            make_dothead_from_file(f_default_match, oheadname, dothead_format, swarp_verbose)

            # swarp does honor the imagesize if set in a FITS format .head file, but not if test format.
            test_cmd3 = ""
        else:
            raise RuntimeError(f"Error, center type {do_better_test_center_type} has not yet been implemented.")

        test_cmd4 = f"-IMAGEOUT_NAME {ofilename} -WRITE_FILEINFO Y -WRITE_XML N -DELETE_TMPFILES Y -RESAMPLE_DIR ./ -WEIGHTOUT_NAME {owgtsname}"
        test_cmd_list = (shlex.split(test_cmd1)
            + shlex.split(test_cmd2)
            + shlex.split(test_cmd3)
            + shlex.split(test_cmd4))
        if swarp_verbose:
            print(f"Command string to be passed to subprocess.run is {test_cmd_list}")

        # Run swarp
        fs_tstart = time.perf_counter()
        try:
            result = subprocess.run(
                test_cmd_list,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            fs_tend = time.perf_counter()
            fs_telapsed = fs_tend - fs_tstart  # seconds
            print(f"  Success: process took {fs_telapsed:.3f} seconds, return code {result.returncode}")
            final_images.append( ofilename )

        except subprocess.CalledProcessError as err:
            fs_tend = time.perf_counter()
            fs_telapsed = fs_tend - fs_tstart  # seconds
            print(f"  Error, process took {fs_telapsed:.3f} seconds, return code {err.returncode}")
            print(f"  Input args: {err.cmd}\n")
            print(f"  Stdout: {err.output}\n")
            print(f"  Stderr: {err.stderr}\n")
            print(f'  Command line equivalent command: {" ".join(err.cmd)}')

# Generate a summary
print(f'Finished. Generated {len(final_images)} resampled images:')
show_full_wcs = True
for img in final_images:
    extnum = 0
    hdu    = fits.open(img)[extnum]
    w      = wcs.WCS(hdu.header)
    naxis1 = hdu.header['NAXIS1']
    naxis2 = hdu.header['NAXIS2']
    print(f'  {img}  naxis1={naxis1}  naxis2={naxis1}')
    if show_full_wcs:
        print(w)
        print(80*'-')


# ### Inspecting The Resampled Images
# 
# We can load one or more of the resampled images and plot them... Here we'll just pick the first one in the `final_images` list.

# In[39]:


from astropy.visualization import (ManualInterval, AsinhStretch,
                                   ImageNormalize)
from astropy import units
def load_image_and_plot(fname, extnum=0, usewcs=True, vmin=None, vmax=None, xaxlim=None, yaxlim=None, verbose=True,
                       angle_tick_spacing_am=2.0, swap_radec_axis=False):
    """
    Quick and dirty FITS image plot.

    By default `asinh` scaling will be used between a minimum and maximum data value. If vmin
    and/or vmnax are not specified then the 0.5th and 99.5th perciles will be used, as appropriate.
    A viridis color map is used by default, arag with a color bar.

    Note
    ----
    Axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis
    pixels.

    Astropy WCSAxes transforms do not reorient the image to be North up, East left, in the same
    way SAOImage ds9 does when the WCS is applied. It still plots the original data rows and
    column as a 2-D XY grid. Consequently, for images not already in NE alignment, the RA may
    change most rapidly along Y and Dec most rapidly along X, in contrast to normal expectation.
    To avoid highly confusing plots this function plots RA/Dec grids when `usewcs=True`. Line 
    of constant RA are red, lines of constant Declination are blue. In cases where your data
    is aligned closer to 90 degrees or 270 degrees away from North up, East left, specifying
    `swap_radec_axis=True` will reduce confusion by showing RA tickmarks and values along the Y
    axes and Declination tickmarks and values along the X axis (contrary to the normal convention).

    If your `usewcs=True` plot appears without axis tick values and labels then it is likely you
    need to alter `swap_radec_axis` or `angle_tick_spacing_am`

    TODO: Find out how to specify them in RA and Dec.

    Parameters
    ----------
    :param fname: Name of existing FITS file with WCS in HDU number extnum that we want to plot.
    :param extnum: Extension number (zero-based) for data and WCS header
    :param usewcs: If True then plot using WCS information
    :param vmin: If not None then vmin is the minimum value in the normalaztion interval.
    :param vmax: If not None then vmax is the maximum value in the normalaztion interval.
    :param xaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True. 
    :param yaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True.
    :param verbose: If True then diagnostic information will be written to stdout.
    :param angle_tick_spacing_am: If usewcs is True, specify the spacing between RA/Dec axis
        tick values in units of arcminutes. Astropy won't label the RA/Dec axis without 
        specifying a value for this.
    :param swap_radec_axis: If True then plot RA tickmarks and values along the Y
        axes and Declination tickmarks and values along the X axis (contrary to the normal 
        convention). This does not alter how the image data itself is plotted.
    """
    
    hdu = fits.open(fname)[extnum]
    w   = wcs.WCS(hdu.header)

    # Compute vmin and vmax if necessary
    # percentiles
    ipctls = [0.5, 99.5]
    opctls = np.nanpercentile(hdu.data, ipctls)
    ourmin = opctls[0]
    ourmax = opctls[1]
    if verbose:
        print(f'Input data 0.5th percentile = {ourmin:.4f}, 99.5th percentile = {ourmax:.4f}')
    
    if vmin is not None:
        ourmin = vmin
    if vmax is not None:
        ourmax = vmax
    if verbose:
        print(f'Applying asinh stretch between {ourmin:.4f} and {ourmax:.4f}')
    
    # Create an ImageNormalize object
    norm = ImageNormalize(hdu.data, interval=ManualInterval(ourmin, ourmax),
                      stretch=AsinhStretch())

    # Display the image
    fig = plt.figure()
    if usewcs:
        if verbose:
            print('Using WCS information for axis projection')
            print(w)
        ax = fig.add_subplot(1, 1, 1, projection=w)
    else:
        ax = fig.add_subplot(1, 1, 1)
        
    im = ax.imshow(hdu.data, origin='lower', norm=norm)
    cbar = fig.colorbar(im, ax=ax, extend='neither', spacing='proportional',
                orientation='vertical', shrink=0.85)
    cbar.set_label(r"Units TBA")
    cbar.ax.tick_params(labelsize=8) 
    title_str = fname#.replace("_", "\_")
    ax.set_title(f'{title_str}', fontsize=8)

    # Display default axis limits
    xlim_used = ax.get_xbound()
    ylim_used = ax.get_ybound()
    if verbose:
        print(f'Default X-axis limits: {xlim_used}')
        print(f'Default Y-axis limits: {ylim_used}')
    
    if usewcs:
        # Adapted rom @astrofrog at https://github.com/astropy/astropy/issues/13458#issuecomment-1242640539
        ra = ax.coords[0]
        dec = ax.coords[1]
        ra.set_major_formatter('dd:mm:ss.ss')
        dec.set_major_formatter('dd:mm:ss.ss')
        ra.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color='red')      # Ticks must be defined for axis tickvals to appear
        dec.set_ticks(spacing=angle_tick_spacing_am * units.arcmin, color='blue')
        ra.set_ticklabel(color='red', fontsize=8)
        dec.set_ticklabel(color='blue', fontsize=8)
        ra.grid(color='red', linestyle='--', alpha=0.6)
        dec.grid(color='blue', linestyle='--', alpha=0.6)
        ra.set_axislabel('Right Ascension (dms)', fontsize=8, color='red')
        dec.set_axislabel('Declination (dms)', fontsize=8, color='blue')

        if swap_radec_axis:
            # Images where ra changes fastest on Y, not X
            dec.set_ticks_position('b')
            dec.set_ticklabel_position('b')
            dec.set_axislabel_position('b')
            ra.set_ticks_position('l')
            ra.set_ticklabel_position('l')
            ra.set_axislabel_position('l')
    
        #ax.set_xlabel('Right Ascension (deg)', fontsize=8)
        #ax.set_ylabel('Declination (deg)', fontsize=8)
    else:
        ax.set_xlabel('X-axis pixel number', fontsize=8)
        ax.set_ylabel('Y-axis pixel number', fontsize=8)

    # Modify the axis limits?
    if xaxlim is not None:
        loval = np.min(xaxlim)
        hival = np.max(xaxlim)
        ax.set_xbound(lower=loval, upper=hival)
        if verbose:
            print(f'Setting X-axis limits to [{loval}, {hival}]')
    if yaxlim is not None:
        loval = np.min(yaxlim)
        hival = np.max(yaxlim)
        ax.set_ybound(lower=loval, upper=hival)
        if verbose:
            print(f'Setting Y-axis limits to [{loval}, {hival}]')
    return


# In[40]:


# 2024-08-23 Something weird with matplotlib plots 
# Basically they don't rotate the image to N up E left the way ds9 does,
# so the coord labelling is odd
# see https://github.com/astropy/astropy/issues/13458


# ---
# **WARNING**
# 
# Astropy's WCS axis plotting projection code does *NOT* work the way you might expect based on other FITS viewers.
# 
# The image itself is not projected to align with the Right Ascension/Declination coordinate system *North up, East left),
# e.g. the way `ds9`'s  Align mode works.
# Instead, the axis tickmark labeling and any over-plotted grid is aware of the actual image WCS, and can correctly label
# and mark the RA and Dec values. This isn't apparent in most cases because many images are provided with the data already
# in North up and East left mode.
# 
# See:
# 
# - Astropy issue discussing this: https://github.com/astropy/astropy/issues/13458 This also discussed getting the axis labels and tickmarks right
# - Astropy request for it to do the right thing: https://github.com/astropy/astropy/issues/8423
# - Work-around using montage to resample the data to NE-aligned: https://stackoverflow.com/questions/51952436/plotting-sdss-images-with-python/51967697#51967697
# 
# ---
# 

# In[41]:


print(f'Final images: {final_images}')
output_plot = 'delme.png'
ap.util.load_image_and_plot(final_images[0], 0, output_plot, True, angle_tick_spacing_am=10, swap_radec_axis=True)


# **Axis limits:** Note that the axis limits are in terms of pixel coordinates, not astronomical RA or Dec.
# 
# **For now, `xaxlim` and `yaxlim` passed to `load_image_and_plot` must be in pixels.**

# In[42]:


# We can plot without the WCS axes, just raw pixels..
output_plot = None
ap.util.load_image_and_plot(final_images[0], 0, output=output_plot, usewcs=False)


# So what do we see here?
# 
# **Problem 1:** Despite specifying a WCS header with a non-zero rotation it appears that `swarp` will always force A North-up, East to the left alignment. Looking at the header of any of the resampled files we can see
# ```bash
# SIMPLE  =                    T / This is a FITS file                            
# BITPIX  =                  -32 /                                                
# NAXIS   =                    2 /                                                
# NAXIS1  =                 3326 / NUMBER OF ELEMENTS ALONG THIS AXIS             
# NAXIS2  =                 3301 / NUMBER OF ELEMENTS ALONG THIS AXIS             
# EXTEND  =                    T / This file may contain FITS extensions          
# EQUINOX =        2000.00000000 / Mean equinox                                   
# RADESYS = 'ICRS    '           / Astrometric system                             
# CTYPE1  = 'RA---TAN'           / WCS projection type for this axis              
# CUNIT1  = 'deg     '           / Axis unit                                      
# CRVAL1  =   2.107644916490E+02 / World coordinate on this axis                  
# CRPIX1  =   1.663500000000E+03 / Reference pixel on this axis                   
# CD1_1   =  -1.752763491822E-04 / Linear projection matrix                       
# CD1_2   =   0.000000000000E+00 / Linear projection matrix                       
# CTYPE2  = 'DEC--TAN'           / WCS projection type for this axis              
# CUNIT2  = 'deg     '           / Axis unit                                      
# CRVAL2  =   5.441565675623E+01 / World coordinate on this axis                  
# CRPIX2  =   1.651000000000E+03 / Reference pixel on this axis                   
# CD2_1   =   0.000000000000E+00 / Linear projection matrix                       
# CD2_2   =   1.752763491822E-04 / Linear projection matrix  
# ```
# which differs from the header we specified in the `.head` file, e.g. in `M101_Supernova_2023ixf_Green_resamp_weighted.head`
# ```bash
#  fitsheader  M101_Supernova_2023ixf_Green_resamp_weighted.head 
# # HDU 0 in M101_Supernova_2023ixf_Green_resamp_weighted.head:
# SIMPLE  =                    T / conforms to FITS standard                      
# BITPIX  =                   64 / array data type                                
# NAXIS   =                    2 / number of array dimensions                     
# NAXIS1  =                 3056                                                  
# NAXIS2  =                 3056                                                  
# WCSAXES =                    2 / Number of coordinate axes                      
# CRPIX1  =               1528.5 / Pixel coordinate of reference point            
# CRPIX2  =               1528.5 / Pixel coordinate of reference point            
# PC1_1   =   -4.64681630389E-06 / Coordinate transformation matrix element       
# PC1_2   =    0.000175223532872 / Coordinate transformation matrix element       
# PC2_1   =   -0.000175223532872 / Coordinate transformation matrix element       
# PC2_2   =   -4.64681630389E-06 / Coordinate transformation matrix element       
# CDELT1  =                  1.0 / [deg] Coordinate increment at reference point  
# CDELT2  =                  1.0 / [deg] Coordinate increment at reference point  
# CUNIT1  = 'deg'                / Units of coordinate increment and value        
# CUNIT2  = 'deg'                / Units of coordinate increment and value        
# CTYPE1  = 'RA---TAN'           / Right ascension, gnomonic projection           
# CTYPE2  = 'DEC--TAN'           / Declination, gnomonic projection               
# CRVAL1  =        210.754804653 / [deg] Coordinate value at reference point      
# CRVAL2  =        54.4174059435 / [deg] Coordinate value at reference point      
# LONPOLE =                180.0 / [deg] Native longitude of celestial pole       
# LATPOLE =        54.4174059435 / [deg] Native latitude of celestial pole        
# MJDREF  =                  0.0 / [d] MJD of fiducial time                       
# RADESYS = 'FK5'                / Equatorial coordinate system                   
# EQUINOX =               2000.0 / [yr] Equinox of equatorial coordinates         
# ```
# 
# **Problem 2:** The sky background in the image is really bright (approx 0.6 per pixel) compared to the maximum data values from the stars and galaxy (about 0.78 per pixel at the 99.5th percentile). We can plot using a hand-chosen minimum, e.g.

# In[43]:


ap.util.load_image_and_plot(final_images[0], 0, output_plot, True, vmin=0.6, vmax=None,  angle_tick_spacing_am=10, swap_radec_axis=True)


# ...at which point we can start to see signs of the third problem.
# 
# Zooming in on a smaller region we see a repeated pattern of two bright points, separated by approximately the maximum offset between the raw frames used to make the composite image. The points are also more compact than real stars in the image.

# In[44]:


ap.util.load_image_and_plot(final_images[0], 0, output_plot, True, vmin=0.6, vmax=1.0, xaxlim=[1100, 1500], yaxlim=[350, 750], verbose=True, angle_tick_spacing_am=2, swap_radec_axis=True)


# In[45]:


# Zoom in even further
ap.util.load_image_and_plot(final_images[0], 0, output_plot, True, vmin=0.6, vmax=1.0, xaxlim=[1240, 1440], yaxlim=[500, 700], verbose=True,  angle_tick_spacing_am=2, swap_radec_axis=True)


# So we have uncorrected hot pixels!
# 
# For now we'll ignore that, and move on.

# ## Three color composite images
# 
# There are two ways we can generate 3-color composite images from the navigated resampled and combined images. 
# 
# One option is to to Astropy's native [make_lupton_rgb](https://docs.astropy.org/en/stable/visualization/rgb.html#astropy-visualization-rgb) function. This is based on [Lupton et al 2004](https://ui.adsabs.harvard.edu/abs/2004PASP..116..133L/abstract). The advantages to this algorithm are it uses the asinh stretch function, arguably the most versatile stretch and the color is unique and not brightness dependent. The downside to this is the lack of exact control over the scaling options, in particular per-channel minimum and maximum values, and also the lack of a 16-bit color depth option. It also has an builtin luminance function L = (R + G + B)/3, which is neither human eye equivalent nor adjustable. In practice Q is typically around 10 for good-loking images, and stretch is typically less than 1 (e.g. 0.5, but 0.02 in some cases).
# 
# The second option is to install Astromatic [stiff](https://vizier.cfa.harvard.edu/vizier/catstd/man/stiff.pdf) and call it via the subprocess module in the same way we used `swarp`. This gives us more options.

# In[46]:


# First, lets just plot each resampled image separately
def plot_threecolor(redfile, greenfile, bluefile, extnum=0, usewcs=True, vmin=None, vmax=None, xaxlim=None, yaxlim=None, verbose=True):
    """
    Quick and dirty FITS image plot.

    By default `asinh` scaling will be used between a minimum and maximum data value. If vmin
    and/or vmnax are not specified then the 0.5th and 99.5th perciles will be used, as appropriate.
    A viridis color map is used by default, along with a color bar.

    Note that axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis
    pixels.

    TODO: Find out how to specify them in RA and Dec.

    :param redfile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as red in the RGB color composite.
    :param greenfile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as green in the RGB color composite.
    :param bluefile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as blue in the RGB color composite.
    :param extnum: Extension number (zero-based) for data and WCS header
    :param usewcs: If True then plot using WCS information
    :param vmin: If not None then vmin is the minimum value in the normalaztion interval.
    :param vmax: If not None then vmax is the maximum value in the normalaztion interval.
    :param xaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True. 
    :param yaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True.
    :param verbose: If True then diagnostic information will be written to stdout.
    """
    
    hdur = fits.open(redfile)[extnum]
    hdug = fits.open(greenfile)[extnum]
    hdub = fits.open(bluefile)[extnum]

    # Compute vmin and vmax if necessary. NOTE this only uses the red channel
    # percentiles
    ipctls = [0.5, 99.5]
    opctls = np.nanpercentile(hdur.data, ipctls)
    ourmin = opctls[0]
    ourmax = opctls[1]
    if verbose:
        print(f'Input red channel data 0.5th percentile = {ourmin:.4f}, 99.5th percentile = {ourmax:.4f}')
    
    if vmin is not None:
        ourmin = vmin
    if vmax is not None:
        ourmax = vmax
    if verbose:
        print(f'Applying asinh stretch between {ourmin:.4f} and {ourmax:.4f}')
    
    # Display the image
    fig = plt.figure()

    for idx in range(3):
        if idx == 0:
            hdu       = hdur
            title_str = redfile
            w         = wcs.WCS(hdu.header)
            print(80*'-')
            print(f'WCS for {title_str}:')
            print(w)
        elif idx == 1:
            hdu       = hdug
            title_str = greenfile
            w         = wcs.WCS(hdu.header)
            print(80*'-')
            print(f'WCS for {title_str}:')
            print(w)
        elif idx == 2:
            hdu       = hdub
            title_str = bluefile
            w         = wcs.WCS(hdu.header)
            print(80*'-')
            print(f'WCS for {title_str}:')
            print(w)
        else:
            raise RuntimeError(f'Error, index {idx} should be in range 0..2 for 3-image plotting.')

        if usewcs:
            ax = fig.add_subplot(2, 2, idx+1, projection=w)
        else:
            ax = fig.add_subplot(2, 2, idx+1)

        # Create an ImageNormalize object
        norm = ImageNormalize(hdu.data, interval=ManualInterval(ourmin, ourmax),
                      stretch=AsinhStretch())
        
        im = ax.imshow(hdu.data, origin='lower', norm=norm)
        fig.colorbar(im, ax=ax)

        # Display default axis limits
        xlim_used = ax.get_xbound()
        ylim_used = ax.get_ybound()
        if verbose:
            print(f'Default X-axis limits: {xlim_used}')
            print(f'Default Y-axis limits: {ylim_used}')
        
        ax.set_title(f'{title_str}', fontsize=8)
        if usewcs:
            ax.set_xlabel('Right Ascension (deg)', fontsize=8)
            ax.set_ylabel('Declination (deg)', fontsize=8)
        else:
            ax.set_xlabel('X-axis pixel number', fontsize=8)
            ax.set_ylabel('Y-axis pixel number', fontsize=8)
    
        # Modify the axis limits?
        if xaxlim is not None:
            loval = np.min(xaxlim)
            hival = np.max(xaxlim)
            ax.set_xbound(lower=loval, upper=hival)
            if verbose:
                print(f'Setting X-axis limits to [{loval}, {hival}]')
        if yaxlim is not None:
            loval = np.min(yaxlim)
            hival = np.max(yaxlim)
            ax.set_ybound(lower=loval, upper=hival)
            if verbose:
                print(f'Setting Y-axis limits to [{loval}, {hival}]')
    return


# In[47]:


bluef  = 'M101_Supernova_2023ixf_Blue_resamp_weighted.fits'
greenf = 'M101_Supernova_2023ixf_Green_resamp_weighted.fits'
lumf   = 'M101_Supernova_2023ixf_Lum_resamp_weighted.fits'
redf   = 'M101_Supernova_2023ixf_Red_resamp_weighted.fits'
valmin  = None
valmax  = None
ap.util.load_imlist_and_plot([redf, greenf, bluef], extnum=0, output=output_plot, usewcs=True, vmin=0.6, vmax=1.6, xaxlim=None, yaxlim=None, verbose=True, angle_tick_spacing_am=10, swap_radec_axis=True)


# ### Astropy make_lupton_rgb

# In[48]:


from astropy.visualization import make_lupton_rgb
def plot_lupton_threecolor(redfile, greenfile, bluefile, outpngfile, extnum=0, usewcs=True, vmin=None, xaxlim=None, yaxlim=None, qval=10, stretchval=0.5, verbose=True):
    """
    Quick and dirty 3-color composite using make_lupton_rgb

    Note: 
    - All images must share the same size, WCS orientation, and pixel size.
    - Axis limits (xaxlim and yaxlim) are currently defined in image x-axis and y-axis pixels.

    TODO: Find out how to specify them in RA and Dec.

    :param redfile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as red in the RGB color composite.
    :param greenfile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as green in the RGB color composite.
    :param bluefile: Name of existing FITS file with WCS in HDU number extnum that we want to
      appear as blue in the RGB color composite.
    :param outpngfile: Name of PNG version of composite the generate.
    :param extnum: Extension number (zero-based) for data and WCS header
    :param usewcs: If True then plot using WCS information
    :param vmin: If not None then vmin is the minimum value the image, i.e. that will correspond to black
    :param xaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True. 
    :param yaxlim: 2-element list or tuple of the minimum to maximum x-axis coordinates to
      plot. If None then the default axis limits will be used. To see those limits run
      with verbose=True.
    :param qval: Value of make_lupton_rgb Q parameter (Q in Lupton et al 2004)
    :param stretchval: Value of make_lupton_rgb stretch parameter (alpha in Lupton et al 2004)
    :param verbose: If True then diagnostic information will be written to stdout.
    """
    
    hdur = fits.open(redfile)[extnum]
    hdug = fits.open(greenfile)[extnum]
    hdub = fits.open(bluefile)[extnum]
    w    = wcs.WCS(hdur.header)

    # Compute vmin and vmax if necessary. NOTE this only uses the red channel
    # percentiles
    ipctls = [0.5, 99.5]
    opctls = np.nanpercentile(hdur.data, ipctls)
    ourmin = opctls[0]
    ourmax = opctls[1]
    if verbose:
        print(f'Input red channel data 0.5th percentile = {ourmin:.4f}, 99.5th percentile = {ourmax:.4f}')
    
    vmax = None
    if vmin is not None:
        ourmin = vmin
    if vmax is not None:
        ourmax = vmax
    if verbose:
        print(f'Using an image minimum of {ourmin} for all channels.')

    rgb_array = make_lupton_rgb(hdur.data, hdug.data, hdub.data,
                               minimum=ourmin, stretch=stretchval, Q=qval, filename=outpngfile)
    
    # Display the image
    fig = plt.figure()
    if usewcs:
        ax = fig.add_subplot(1, 1, 1, projection=w)
    else:
        ax = fig.add_subplot(1, 1, 1)
        
    im = ax.imshow(rgb_array, origin='lower')
    fig.colorbar(im, ax=ax)

    # Display default axis limits
    xlim_used = ax.get_xbound()
    ylim_used = ax.get_ybound()
    if verbose:
        print(f'Default X-axis limits: {xlim_used}')
        print(f'Default Y-axis limits: {ylim_used}')
    
    title_str = f'RGB composite image stretch={stretchval}, Q={qval}\nRed = {redfile}\nGreen = {greenfile}\nBlue = {bluefile}\n'
    ax.set_title(f'{title_str}', fontsize=8)
    if usewcs:
        ax.set_xlabel('Right Ascension (deg)', fontsize=8)
        ax.set_ylabel('Declination (deg)', fontsize=8)
    else:
        ax.set_xlabel('X-axis pixel number', fontsize=8)
        ax.set_ylabel('Y-axis pixel number', fontsize=8)

    # Modify the axis limits?
    if xaxlim is not None:
        loval = np.min(xaxlim)
        hival = np.max(xaxlim)
        ax.set_xbound(lower=loval, upper=hival)
        if verbose:
            print(f'Setting X-axis limits to [{loval}, {hival}]')
    if yaxlim is not None:
        loval = np.min(yaxlim)
        hival = np.max(yaxlim)
        ax.set_ybound(lower=loval, upper=hival)
        if verbose:
            print(f'Setting Y-axis limits to [{loval}, {hival}]')
    return rgb_array


# In[49]:


bluef  = 'M101_Supernova_2023ixf_Blue_resamp_weighted.fits'
greenf = 'M101_Supernova_2023ixf_Green_resamp_weighted.fits'
lumf   = 'M101_Supernova_2023ixf_Lum_resamp_weighted.fits'
redf   = 'M101_Supernova_2023ixf_Red_resamp_weighted.fits'
qval    = 8
stretch = 0.5
valmin  = None
outf    = 'M101_lupton_RedGreenBlue.png'
rgb = ap.util.plot_lupton_threecolor(redf, greenf, bluef, outf, extnum=0, usewcs=True, vmin=valmin, xaxlim=None, yaxlim=None, qval=qval, stretchval=stretch, verbose=True, angle_tick_spacing_am=10, swap_radec_axis=True)


# **Analysis:** That looks reasonably decent from a professional astronomer's point of view, but isn't quite good enough for the more display-oriented astrophotographer's point of view: 
#  - The cast if the background is too green, given the higher background in the green channel compared to the other images.
#  - The galaxy is too red.
#  - The lack of control over the image maximum value means the galaxy is too faint (reducing Q from the default would help). The feature of the Lupton method not saturating the stars is counterproductive in this case.
#  - My python wrapper has a bug with the axis tickvals and labels, and the color bar is unnecessary.
# 
# Of course we can wrap the method in some more pre-processing to fix some of those issues, but showing the default is informative.

# ### Astromatic STIFF
# 
# If you have stiff installed then we can use the python `subprocess` module to run it. Invoking the commands by hand on the command line this would look as follows:
# ```bash
# stiff M101_Supernova_2023ixf_Red_resamp_weighted.fits M101_Supernova_2023ixf_Green_resamp_weighted.fits M101_Supernova_2023ixf_Blue_resamp_weighted.fits \
#     -OUTFILE_NAME M101_Supernova_2023ixf_RedGreenBlue_resamp_weighted_gf10_cs10_b8.tiff -GAMMA_TYPE POWER-LAW -GAMMA 2.2 -GAMMA_FAC 1.0 \
#     -COLOUR_SAT 1.0 -MAX_TYPE =QUANTILE,QUANTILE,QUANTILE -MAX_LEVEL 0.999,0.999,0.999 -MIN_TYPE =QUANTILE,QUANTILE,QUANTILE -MIN_LEVEL 0.60,0.60,0.60 \
#     -BITS_PER_CHANNEL 8 -VERBOSE_TYPE FULL -DESCRIPTION M101_Supernova_2023ixf -COPYRIGHT dks -WRITE_XML N
# 
# > WARNING: stiff.conf not found, using internal defaults
# 
# ----- STIFF 2.7.1 started on 2024-04-04 at 14:33:23 with 16 threads
# 
# ----- Inputs:
# M101_Supernova_2023ixf_Red_resamp_weighted.fits: "no ident"  3056x3056   32 bits (floats)
# Background level: 0.634998    Min level: 0.618508    Max level: 3.07458   
# M101_Supernova_2023ixf_Green_resamp_weighted.fits: "no ident"  3056x3056   32 bits (floats)
# Background level: 0.700548    Min level: 0.683503    Max level: 2.88235   
# M101_Supernova_2023ixf_Blue_resamp_weighted.fits: "no ident"  3056x3056   32 bits (floats)
# Background level: 0.607307    Min level: 0.590825    Max level: 2.50276   
# 
# ----- Output:
# M101_Supernova_2023ixf_RedGreenBlue_resamp_weighted_gf10_cs10_b8.tiff:    3056x3056        3x8  bits (integers) gamma: x1.00  compression: LZW 
# 
#       Generated M101_Supernova_2023ixf_RedGreenBlue_resamp_weighted_gf10_cs10_b8.tiff
#       Stiff took 0.713 seconds.
# 
# ```

# In[50]:


# where are we?
cwd = pathlib.Path('.')
print(f'Current directory: {cwd.resolve()}')
tiff_file = cwd / 'M101_Supernova_2023ixf_RedGreenBlue_resamp_weighted_gf10_cs10_b8.tiff'
if tiff_file.exists():
    fig, (ax1, ax2) = plt.subplots(1,2)
    print(f'Loading {tiff_file.name}')
    tiff_img = plt.imread(tiff_file)
    print(f'Plotting full image {tiff_file}')
    ax1.imshow(tiff_img)

    # Zoom in
    print(f'Plotting zoomed-in part of {tiff_file}')
    ax2.imshow(tiff_img[600:1000,1600:2000])    
else:
    print(f'File not found: {tiff_file.name}')
    print('  Skipping display of STIFF generated file...')


# The full image looks a little better, although note that the setting for the image minimum (black-level) are quite different from the make_lupton_rgb case. Zooming in again shows the bad pixels.
# 
# Anyway, on to running STIFF from python...

# In[51]:


stiff_verbose  = True  # Echo stiff input string if True

# Hard-coded, same as with the Lupton case, sorry...
bluef  = 'M101_Supernova_2023ixf_Blue_resamp_weighted.fits'
greenf = 'M101_Supernova_2023ixf_Green_resamp_weighted.fits'
redf   = 'M101_Supernova_2023ixf_Red_resamp_weighted.fits'
description = 'M101_Supernova_2023ixf'
user   = os.getenv("USER")

# Name for output TIFF image
ofilename = f"M101_Supernova_2023ixf_RedGreenBlue_resamp_weighted_gf10_cs10_b8.tiff"

# Quantiles Red channel green channel blue channel
# Note these are fraction in the rnage [0:1], not percentages
quantiles_lo = [0.600, 0.600, 0.600]
quantiles_hi = [0.999, 0.999, 0.999]

# See STIFF manual at https://vizier.cfa.harvard.edu/vizier/catstd/man/stiff.pdf
gamma = 2.2
gamma_factor = 1.0
color_sat = 1.0
bitsperpix = 8

fs_tstart = time.perf_counter()
test_cmd1  = f"stiff {redf} {greenf} {bluef} -OUTFILE_NAME {ofilename}"
test_cmd2  = f"-GAMMA_TYPE POWER-LAW -GAMMA {gamma} -GAMMA_FAC {gamma_factor} -COLOUR_SAT {color_sat}"
test_cmd3  = f"-MAX_TYPE =QUANTILE,QUANTILE,QUANTILE -MAX_LEVEL {quantiles_hi[0]:.4f},{quantiles_hi[1]:.4f},{quantiles_hi[2]:.4f}"
test_cmd4  = f"-MIN_TYPE =QUANTILE,QUANTILE,QUANTILE -MIN_LEVEL {quantiles_lo[0]:.4f},{quantiles_lo[1]:.4f},{quantiles_lo[2]:.4f}"
test_cmd5  = f"-BITS_PER_CHANNEL {bitsperpix} -VERBOSE_TYPE FULL -DESCRIPTION {description} -COPYRIGHT {user} -WRITE_XML N"
test_cmd_list = (shlex.split(test_cmd1)
    + shlex.split(test_cmd2)
    + shlex.split(test_cmd3)
    + shlex.split(test_cmd4)
    + shlex.split(test_cmd5))
if stiff_verbose:
    print(f"Command string to be passed to subprocess.run is {test_cmd_list}")

# Run stiff
fs_tstart = time.perf_counter()
try:
    result = subprocess.run(
        test_cmd_list,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    fs_tend = time.perf_counter()
    fs_telapsed = fs_tend - fs_tstart  # seconds
    print(f"  Success: process took {fs_telapsed:.3f} seconds, return code {result.returncode}")
    final_images.append( ofilename )

except subprocess.CalledProcessError as err:
    fs_tend = time.perf_counter()
    fs_telapsed = fs_tend - fs_tstart  # seconds
    print(f"  Error, process took {fs_telapsed:.3f} seconds, return code {err.returncode}")
    print(f"  Input args: {err.cmd}\n")
    print(f"  Stdout: {err.output}\n")
    print(f"  Stderr: {err.stderr}\n")
    print(f'  Command line equivalent command: {" ".join(err.cmd)}')

# Generate a summary
tiff_file = cwd / ofilename
if tiff_file.exists():
    fig, (ax1, ax2) = plt.subplots(1,2)
    print(f'Loading {tiff_file.name}')
    tiff_img = plt.imread(tiff_file)
    print(f'Plotting full image {tiff_file}')
    ax1.imshow(tiff_img)

    # Zoom in
    print(f'Plotting zoomed-in part of {tiff_file}')
    ax2.imshow(tiff_img[600:1000,1600:2000])    
else:
    print(f'Error, file not found: {tiff_file.name}')
    print('  STIFF generated file missing')


# #### Other methods: composite_all.sh shell script
# 
# The bash script `composite_all.sh` uses Astromatic's `stiff` to combine FITS images from multiple filters into RGB tiff files. It iterates through multiple combinations of output image bit-depth, color saturation and gamma factor. The user must then choose the best version to work on further using an image manipulation program such as the `gimp`. 
# 
# ```bash
# cd /old_lnx/home/dks/Downloads/iTelescopeScratch/Plan-20-M101SN/Resampled/
# ~/git/AstroPhotography/AstroPhotography/scripts/composite_all.sh M101_Supernova_2023ixf 4096x4096_resamp_median.fits rgb
# /home/dks/git/AstroPhotography/AstroPhotography/scripts/composite_all.sh started at Sun Feb 25 09:05:32 AM EST 2024
# Running from /old_lnx/home/dks/Downloads/iTelescopeScratch/Plan-20-M101SN/Resampled
# Logging file:     tmp_composite_all.log
# File prefix:      M101_Supernova_2023ixf
# Fill suffix:      4096x4096_resamp_median.fits
# Color selections: rgb
# Using stiff from /usr/local/bin/stiff
# -----------------------------------------------------------------
# Processing 3-color combination rgb
#   Filters in red, green, blue order: Red Green Blue
#   Checking that input FITS files exist:
#     Found M101_Supernova_2023ixf_Red_4096x4096_resamp_median.fits
#     Found M101_Supernova_2023ixf_Green_4096x4096_resamp_median.fits
#     Found M101_Supernova_2023ixf_Blue_4096x4096_resamp_median.fits
#     Bits-per-pixel=8, gamma-fac=1.0, color_sat=1.0
#     Running stiff, generating M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf10_cs10_b8.tiff at Sun Feb 25 09:05:32 AM EST 2024
#       Generated M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf10_cs10_b8.tiff
#       Stiff took 0.797 seconds.
# 
#     Bits-per-pixel=8, gamma-fac=1.0, color_sat=1.5
#     Running stiff, generating M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf10_cs15_b8.tiff at Sun Feb 25 09:05:33 AM EST 2024
#       Generated M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf10_cs15_b8.tiff
#       Stiff took 0.731 seconds.
# 
#     Bits-per-pixel=8, gamma-fac=1.0, color_sat=2.0
#     Running stiff, generating M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf10_cs20_b8.tiff at Sun Feb 25 09:05:33 AM EST 2024
#       Generated M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf10_cs20_b8.tiff
#       Stiff took 0.729 seconds.
# 
#     Bits-per-pixel=8, gamma-fac=1.2, color_sat=1.0
#     Running stiff, generating M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf12_cs10_b8.tiff at Sun Feb 25 09:05:34 AM EST 2024
#       Generated M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf12_cs10_b8.tiff
#       Stiff took 0.734 seconds.
# 
#     Bits-per-pixel=8, gamma-fac=1.2, color_sat=1.5
#     Running stiff, generating M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf12_cs15_b8.tiff at Sun Feb 25 09:05:35 AM EST 2024
#       Generated M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf12_cs15_b8.tiff
#       Stiff took 0.716 seconds.
# 
#     Bits-per-pixel=8, gamma-fac=1.2, color_sat=2.0
#     Running stiff, generating M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf12_cs20_b8.tiff at Sun Feb 25 09:05:36 AM EST 2024
#       Generated M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf12_cs20_b8.tiff
#       Stiff took 0.721 seconds.
# 
#     Bits-per-pixel=8, gamma-fac=1.4, color_sat=1.0
#     Running stiff, generating M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf14_cs10_b8.tiff at Sun Feb 25 09:05:36 AM EST 2024
#       Generated M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf14_cs10_b8.tiff
#       Stiff took 0.732 seconds.
# 
#     Bits-per-pixel=8, gamma-fac=1.4, color_sat=1.5
#     Running stiff, generating M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf14_cs15_b8.tiff at Sun Feb 25 09:05:37 AM EST 2024
#       Generated M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf14_cs15_b8.tiff
#       Stiff took 0.714 seconds.
# 
#     Bits-per-pixel=8, gamma-fac=1.4, color_sat=2.0
#     Running stiff, generating M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf14_cs20_b8.tiff at Sun Feb 25 09:05:38 AM EST 2024
#       Generated M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf14_cs20_b8.tiff
#       Stiff took 0.726 seconds.
# 
# 
# -----------------------------------------------------------------
# Run summary:
# Img  Resampled Output                                                  Time  Status
# 0    M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf10_cs10_b8.tiff     0.797       0
# 1    M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf10_cs15_b8.tiff     0.731       0
# 2    M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf10_cs20_b8.tiff     0.729       0
# 3    M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf12_cs10_b8.tiff     0.734       0
# 4    M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf12_cs15_b8.tiff     0.716       0
# 5    M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf12_cs20_b8.tiff     0.721       0
# 6    M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf14_cs10_b8.tiff     0.732       0
# 7    M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf14_cs15_b8.tiff     0.714       0
# 8    M101_Supernova_2023ixf_RedGreenBlue_4096x4096_resamp_median_gf14_cs20_b8.tiff     0.726       0
# Master log file: tmp_composite_all.log
# /home/dks/git/AstroPhotography/AstroPhotography/scripts/composite_all.sh finished at Sun Feb 25 09:05:39 AM EST 2024, run time 6.716 seconds.
# 
# ```

# # Summary and Conclusions
# 
# We've processed a (small) set of iTelescope Premium pre-calibrated images using a lot of hand-written python wrapping around individual `AstroPhotography` package functions: star detection, astrometric solutions, stacking, and processing color composite images from the stacked images in each band. The results are OK, but still leave a lot to be desired.
# 
# - This particular set of premium images, although calibrated, still include a lot of hot/bad pixels.
#    - Those would be fixed if we were processing from the raw files using our own calibration functions.
#    - We'll need a way of fixing those in any eventual calibrated images pipeline, along with recognizing whether bad pixel correction has been already applied at calibration.
# - The "sky" background in the images is high. This might be a pedestal, or just the ambient conditions.
#    - Again, we'll need a way of fixing that in any eventual calibrated images pipeline, along with recognizing whether sky background correction (and/or pedestal removal) has been already applied at calibration.

# ## Versions and Changes
# 
# | Version | Date | Description |
# |:--------|------|-------------|
# | 0.5.2-alpha | 2024-04-26 | Partial version includes by-hand processing, but not first version of pipeline script |
# | 0.6.0        | 2024-11-26 | Final version worked under issue-002 |

# In[ ]:




