# AstroPhotography Release Notes

Release notes for versions starting at v0.2.2 and later.

**Note:** Standard sematic versioning API stability is not guaranteed
until version 1.0 or later.

## Version 0.6

Version 0.6 continues work on FITS file processing. Work on photographic
RAW file conversion will resume at some point in the future.

The objectives of version 0.6 are:

- Working `sphinx` documentation!
- Semi-tutorial jupyyer notebooks.
- Easier and more pythonic configuration for processing a set of images.
  The example `bash` shell scripts provided in the `scripts` directory
  are neither flexable nor maintainable.
- The ability to process the metadata-less "Premium" datasets now 
  available from iTelescope.

### Version 0.6.2 WIP

> [!CAUTION]
> Updated conda/python environment. 
> The updated code appears to still work with the older `python 3.10` and `astropy 6.0`
> configuration within the limited coverage of my batch processing tests. If you
> experience problems file an issue or create a new conda environment with the current
> `ap-env.yml` and `requirements.txt` files. The older configuration can be found in
> the `etc/` directory.

- #37 Add `ap_gradient_to_mask.py` and `ApGradientToMask`. 
- #30 Updated conda environment and pip `requirements.txt` to `python 3.12` and `astropy 7.0`.
- #30 Various minor code changes forced by newer versions of `ccdproc` and `numpy`.
- Updated `iTelescope_processing.md`
- More `ruff` linting and `mypy` type-hinting edits.

### Version 0.6.1 2025-03-07

- More Ruff linting and type hinting
- Reorganize `util` into more submodules to avoid circular dependencies. 
- #27 Added `ApFitsRasterizer` and CLI `ap_fits_rasterizer`, with only
  the `stiff` backend currently.
- #27  Added `plot_stiff_threecolor` and `make_stiff_threecolor_plots`
  to `ApUtil`.
- #27 `ApProcessor` will generate TIFF format quick look image of resampled output.
- #27 `ApProcessor` can skip regenerating already preprocessed images.
- #27 Fix bug introduced into `ApMasterCal` when linting.
- #27 Add option to `ApProcessor` to allow different swarp `COMBINE_TYPE` when 
  resampling, in particular allow `MEDIAN` instead of `WEIGHTED`. `MEDIAN` strongly
  reduces residual hard-to-correct CR noise but increases background noise 
  levels by 60%, and requires input images are similar in exposure time 
  (and maybe pixel scale and coverage). 
- #32 Improve summary plot generation for multi-degree wide images.
- #32 Better handing of cases where no candidates can be fitted in
  `ApMeasureStars` and `ApFindStars`, with helpful error messaging.
- #32 Include nearest neighbor distance estimate in fit box size calculation.
- Added `patching_holes_from_bad_flats.ipynb`
- Split make_dothead_from into file, keywords, wcs options
- Fix summarize_wcs when there is a `PC?_?` matrix.
- #31 Added ApFixHoles to fill in larger scale user-defined holes.
- #28 Added option to exclude source detection at the corners of images.

### Version 0.6.0 2024-11-26

- Python wrapper for processing multiple calibrated image sets, e.g.
  iTelescope premium image sets. (ApProcessor, ap_processor)
- Walk throughs of processing calibrated images and generating
  bad pixel maps using jupyter notebooks.
- API breaking package re-organization with creation of `util` subpackage.
- API breaking changes to other classes, e.g. ApAstrometry.
- Various minor fixes to ApAstrometry, ApFindBadPixels, etc.
- Continuing work on improving the Sphinx documentation.
- In process switch to Ruff for formatting/linting.

## Version 0.5

### Version 0.5.1 2024-01-30

- Improve documentation now that we can run sphinx on the module.
- Resolved #17, added `yamlkeyval` mode to ApAddMetadata / ap_add_metadata
  that works on iTelescope premium images. 
- Issue #21, catch up to latest astropy (>=6.0) and affiliated package changes.
  Fixes to `ApFindStars`, `ApMeasureStars`, `ApFixCosmicRays`, and
  `ApMeasureBackground`.
- Add AP_CALIB_DIR to `calibrate_all.sh` script, document in `doc/iTelescope_processing.md`

### Version 0.5.0 2022-09-18

Version 0.5 continued work on FITS file processing. Work on photographic
RAW file conversion will resume at some point in the future.

The objectives of version 0.5 were:

- Catch up to latest astropy (>=6.0) and affiliated package changes.
- Working `sphinx` documentation!
- *Package layout reorganization* to move past some limitations imposed
  by the original `cookiecutter` template, that are no longer helpful.

Changes:

- The module itself has moved from `./src/AstroPhotography` to `./AstroPhotography`.
- `sphinx`-generated documenation, including the module docstrings, now
  works. Change directory into the `doc` directory and run `make html`,
  then view the file in `_build/html/index.html`.

Limitations:

- Python docstrings may be incorrect in some places, but that will be 
  addressed at a later date.
- The default font size used in the sphinx `sizzle` theme is too large.

## Version 0.4

Version 0.4 returns to issues related to processing of FITS images,
in particular from iTelescope.

### v0.4.2 2022-09-11

- WIP on describing how to prepare for data processing and new calibration
  files in iTelescope_processing.md.
- WIP on capability to process telescopes other than T05 in 
  `calibrate_all.sh` (although note this script is a temporary clunky
  solution for batch processing).
- Added more options to `ApMeasureBackground` and `ap_measure_background.py`
  that allow user-driven finetuning of the background estimations, based
  on experience with T20 observations of the Cygnus Loop.
- Added `ApAutoBadcols` and `ap_auto_badcols` to automatically detect
  the majority of bad columns and rows when **first** generating a
  user-defined bad pixel calibration file. (Issue #14)
- Improved handling of failed name resolution in `ap_add_metadata.py`
- Added capability to specify target name on `ap_add_metadata.py` in
  cases where file name parsing will not work.
- `ApAddMetadata` will now attempt to remove the mosaic subfield 
  identifiers generated by Telescopius, e.g. `CygnusLoop_x1_y1` is no
  longer treated as `CygnusLoop x1 y1` but as `CygnusLoop` instead.
- Added `ap_fix_itelescope_dirs.sh` and `ap_rename_files_with_spaces.sh`

### v0.4.1 2021-07-24

- Updated the requirements.txt file with long term aim of supporting python39.
- Added ap_measure_background.py and ApMeasureBackground.py to resolve
  issue #8. These measure and output an image of the large scale 
  background in an image that can be used to correct for imperfect 
  bias/dark/flat calibration or moon illumination. 
- Big changes in the ApFindStars interface to support it being used by
  external classes instead of doing all I/O on files.

### v0.4.0 2021-05-22

- Added ap_imarith.py and ApImArith to supply a native way of performing
  simple arithmetical operations on one or two FITS files.
- Added example `bash` driver scripts that can process multiple images
  from the observation of a single target. They're somewhat clunky to 
  set up, so they're only provided as an example.
  
Limitations:
- ap_imarith.py and ApImArith assume the data is always the first 
  extension (Primary HDU) of the FITS file. This does not allow for 
  multi-extension files (MEF).

## Version 0.3

Version 0.3 is an interim release that attempts to finally provide
the basic "required" functionality to the RAW file conversion 
tool `dksraw`, before attention switches back to FITS file processing.  

### v0.3.0 2021-04-04

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

## Version 0.2

### v0.2.4 2021-02-27

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

### v0.2.3 2020-12-22

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

### v0.2.2 2020-07-12

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
