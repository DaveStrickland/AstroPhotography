# AstroPhotography Release Notes

Release notes for versions starting at v0.2.2 and later.

# 0.2

## Development version

The emphasis of this release improvements in the calibration process,
in particular bad pixel removal and the newly added cosmic ray removal.

- Added cosmic ray removal using ccdproc's lacosmic algorithm in 
  ApFixCosmicRays, invoked either by ap_fix_cosmic_rays.py or by
  ApCalibrate/ap_calibrate.
- Disabled the default tweak_order=2 SIP polynomial fit in astrometry.net
- Added option to apply user-defined bad pixels in ap_find_badpix.py and
  ApFindBadPixels.

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
