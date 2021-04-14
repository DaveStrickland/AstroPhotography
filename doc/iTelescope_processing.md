# Processing images from iTelescope.net

## Processing Stages

The following table outlines conceptual pipeline processing stages
for astronomical FITS imagery, primarily those obtained using iTelescope.

| Stage                  | Stage Name          | Description                     | Python Class | Package Used      | Command Line Tool?     |
| ---------------------- | ------------------- | ------------------------------- | ------------ | ----------------- | ---------------------- |
| Calibration Generation | Build Master Cal    | Build master calibration files  | Yes          | astropy.ccdproc   | ap_combine_darks.py    |
| Calibration Generation | Calc Read Noise     | Calculate detector read noise   | Yes          |                   | ap_calc_read_noise.py  |
| Calibration Generation | Find Bad Pixels     | Generate a bad pixel map        | Yes          |                   | ap_find_badpix.py      |
| Calibration            | Apply Master Cal    | Apply bias/dark/flat correction | Yes          |                   | ap_calibrate.py        |
| Calibration            | Add Metadata        | Add metadata to calibrated FITS | Yes          | astroplan         | ap_add_metadata.py     |
| Calibration            | Fix Bad Pixels      | Correct bad pixels              | Yes          |                   | ap_fix_badpix.py       |
| Calibration            | Cosmic Ray Reject   | Find and correct cosmic rays    | Yes          | astropy.ccdproc   | Not yet                |
| Image Processing       | Find Stars And PSF  | Find stars, measure PSF         | Yes          | astropy.photutils | ap_find_stars.py       |
| Image Processing       | Astrometry          | Astrometric solution per image  | Yes          | astrometry.net    | ap_astrometry.py       |
| Image Processing       | Quality Filter      | Identify "poor quality" images  | Yes          |                   | ap_find_stars.py       |
|                        |                     | and summarize set of images     | Yes          |                   | ap_quality_summary.py  |
| Image Processing       | Resample            | Resample an image to a WCS      | Not yet      | astromatic swarp* | (resample_all.sh**)    |
| Image Combination      | Continuum Subtract  | Continuum scaling/subtraction   | Not yet      |                   | Not yet                |
| Image Combination      | RGB Color Composite | Make 3-color composites         | Not yet      | astromatic stiff* | (composite_all.sh**)   | 
| Process All Images     | Calibrate All       | Calibrate all raw images        | Not yet      |                   | (calibrate_all.sh**)   |
| Process All Images     | Navigate All        | Astrometry/WCS on all images    | Not yet      |                   | (navigate_all.sh**)    |
| Process All Images     | Resample All        | Resample an image to a WCS      | Not yet      | astromatic swarp* | (resample_all.sh**)    |

Notes:

 - (*) Temporarily using an external non-python-based tool. 
 - (**) A temporary bash implementation.

## Pipeline Processing

