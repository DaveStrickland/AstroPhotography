# AstroPhotography Release Notes

Release notes for versions starting at v0.2.2 and later.

# Version 0.4

Version 0.4 returns to issues related to processing of FITS images,
in particular from iTelescope.

## v0.4.1 2021-07-24

- Updated the requirements.txt file with long term aim of supporting python39.
- Added ap_measure_background.py and ApMeasureBackground.py to resolve
  issue #8. These measure and output an image of the large scale 
  background in an image that can be used to correct for imperfect 
  bias/dark/flat calibration or moon illumination. 
- Big changes in the ApFindStars interface to support it being used by
  external classes instead of doing all I/O on files.

## v0.4.0 2021-05-22

- Added ap_imarith.py and ApImArith to supply a native way of performing
  simple arithmetical operations on one or two FITS files.
- Added example `bash` driver scripts that can process multiple images
  from the observation of a single target. They're somewhat clunky to 
  set up, so they're only provided as an example.
  
Limitations:
- ap_imarith.py and ApImArith assume the data is always the first 
  extension (Primary HDU) of the FITS file. This does not allow for 
  multi-extension files (MEF).

# Version 0.3

Version 0.3 is an interim release that attempts to finally provide
the basic "required" functionality to the RAW file conversion 
tool `dksraw`, before attention switches back to FITS file processing.  

## v0.3.0 2021-04-04

- FITS files can now be written by the `dksraw grey`, `dksraw rgb`,
  and `dksraw split` commands. This expects you to specify file names 
  ending in in the full `.fits` (not `.fit`).
- There is a new command `dksraw rgb`, that generates 3-channel RGB
  images from RAW files. If the output file is FITS format the red,
  green and blue channels are written to separate 2-D image extensions
  within the same file. They can be loaded as a single RGB image in
  `ds9` using the command line, e.g. 
  `ds9 -rgb -red test_rgb.fits[0] -green test_rgb.fits[1] -blue test_rgb.fits[2]`
- The `dksraw grey` and new `dksraw rgb` commands now have a (default)
  `linear` method for Bayer demasking and looks much better than the 
  old `direct` method.
- The `dksraw grey` and `dksraw rgb` commands have a `--renormalize` 
  option suitable for generating quick-look images that fill the dynamic 
  range of the output file.
- RAW file EXIF metadata is now read within `RawConv`, using the python 
  package `ExifRead`. This is less capable than a command line tool like 
  `exiv2`. Nevertheless a limited set of this metadata is written to 
  FITS files (if used as output).  
- Minor CLI improvements: 
  - Renamed `--warn` to `--loglevel` for consistency with the ap_ scripts,
    and changed how it was specified in argparse to avoid having to specify
    it before the subcommand.
  - Made "black" level subtraction the default in `dksraw`, so now you have
    to specify `--keepblack` if you want the black levels to be retained.
- The old `TODO.md` file has been removed as it was very out of date. See
  the github Issue for a clearer view of development priorities.

Limitations:

- EXIF metadata writing to graphics files formats (tiff, jpg, png etc)
  has not been successfully implemented, as there is a disconnect 
  between python packages that can read metadata from RAW file and those
  that claim to write valid EXIF metadata to graphics file formats
  that I have not had time to solve.
- Sphinx documentation still does not work.
- We are still using `imageio` for graphics-file format output, which is
  **very slow for formats other than TIFF**. 

# Version 0.2

## v0.2.4 2021-02-27

The emphasis of this release improvements in the calibration process,
in particular bad pixel removal and the newly added cosmic ray removal.

- Added cosmic ray removal using ccdproc's lacosmic algorithm in 
  ApFixCosmicRays, invoked either by ap_fix_cosmic_rays.py or by
  ApCalibrate/ap_calibrate.
- Disabled the default tweak_order=2 SIP polynomial fit in astrometry.net
- Added option to apply user-defined bad pixels in ap_find_badpix.py and
  ApFindBadPixels.
- Fix issue #5, allowing darks that either have bias subtraction or
  have not (--dark_still_biased).

## v0.2.3 2020-12-22

The emphasis of this release is on the capability to process FITS images
from amatuer astronomical observations, in particular iTelescope.net.

The ap_ scripts now allow iTelescope.net data to be calibrated from raw
fits frames, through bias/dark/flat correction, star detection and WCS
solutions.

Outstanding Issues:
- The ap_ scripts are not installed as console scripts, and some still
  don't have their Ap implementation classes available externally.
- sphinx documentation is not working.
- RAW conversion still doesn't output FITS.
- The Ap classes and ap_ scripts were tested by analysis, not with 
  automated unit tests.

## v0.2.2 2020-07-12

Continued work on conversion of camera RAW files, including conversion
to grey-scale 16-bit formats, splitting by Bayer channels and white 
balance related functions. Output is standard imaging formats, e.g. png.

Outsanding Issues:
- Conversions to FITS format not yet implemented, although this should
  be trivial.
- Python imageio is simple and easy to work with, but doesn't seem to
  properly support metadata.
- Unit tests are clunky and depend on external "truth" calculations 
  made in octave/Matlab
