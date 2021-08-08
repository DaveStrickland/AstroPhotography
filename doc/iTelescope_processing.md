# Processing images from iTelescope.net

## Processing Stages

The following table outlines conceptual pipeline processing stages
for astronomical FITS imagery, primarily those obtained using iTelescope.

| Stage                  | Stage Name            | Description                     | Python Class | Package Used      | Command Line Tool?       |
| ---------------------- | --------------------- | ------------------------------- | ------------ | ----------------- | ------------------------ |
| Calibration Generation | Build Master Cal      | Build master calibration files  | Yes          | astropy.ccdproc   | ap_combine_darks.py      |
| Calibration Generation | Calc Read Noise       | Calculate detector read noise   | Yes          |                   | ap_calc_read_noise.py    |
| Calibration Generation | Find Bad Pixels       | Generate a bad pixel map        | Yes          |                   | ap_find_badpix.py        |
| Calibration            | Apply Master Cal      | Apply bias/dark/flat correction | Yes          |                   | ap_calibrate.py          |
| Calibration            | Add Metadata          | Add metadata to calibrated FITS | Yes          | astroplan         | ap_add_metadata.py       |
| Calibration            | Fix Bad Pixels        | Correct bad pixels              | Yes          |                   | ap_fix_badpix.py         |
| Calibration            | Cosmic Ray Reject     | Find and correct cosmic rays    | Yes          | astropy.ccdproc   | Not yet                  |
| Image Processing       | Find Stars And PSF    | Find stars, measure PSF         | Yes          | astropy.photutils | ap_find_stars.py         |
| Image Processing       | Astrometry            | Astrometric solution per image  | Yes          | astrometry.net    | ap_astrometry.py         |
| Image Processing       | Quality Filter        | Identify "poor quality" images  | Yes          |                   | ap_find_stars.py         |
|                        |                       | and summarize set of images     | Yes          |                   | ap_quality_summary.py    |
| Image Processing       | Image Arithmetic      | Arithmetic on 1 or 2 images     | Yes          |                   | ap_imarith.py            |
| Image Processing       | Background Estimation | Estimate sky background         | Yes          |                   | ap_measure_background.py |
| Image Processing       | Resample              | Resample an image to a WCS      | Not yet      | astromatic swarp* | (resample_all.sh**)      |
| Image Combination      | Continuum Subtract    | Continuum scaling/subtraction   | Not yet      |                   | Not yet                  |
| Image Combination      | RGB Color Composite   | Make 3-color composites         | Not yet      | astromatic stiff* | (composite_all.sh**)     | 
| Process All Images     | Calibrate All         | Calibrate all raw images        | Not yet      |                   | (calibrate_all.sh**)     |
| Process All Images     | Navigate All          | Astrometry/WCS on all images    | Not yet      |                   | (navigate_all.sh**)      |
| Process All Images     | Resample All          | Resample an image to a WCS      | Not yet      | astromatic swarp* | (resample_all.sh**)      |

Notes:

 - (*) Temporarily using an external non-python-based tool. 
 - (**) A temporary bash implementation.

## Data Preparation

After downloading the new data (or calibration data) from the iTelescope
FTP site you should run `ap_fix_itelescope_dirs.sh` to correct the
directory permissions and remove any spaces from the directory names.

### Observation Data

To unpack the data, for example for a new `T05` observation of `MyTarget`
on `yyyymmdd`, unpack the zip files while also removing problematic space
from file names (the following commands assume the `bash` shell):

```bash
cd T05/MyTarget/yyyymmdd/
bash $PATH_TO_AP/src/AstroPhotography/scripts/ap_rename_files_with_spaces.sh
for file in *.zip; do unzip $file; echo ""; done
bash $PATH_TO_AP/src/AstroPhotography/scripts/ap_rename_files_with_spaces.sh
rm *.zip
```

### Calibration Data Preparation

A few addition steps must be taken when downloading new calibration
data from the iTelescope website:

- Fix directory permissions using `ap_fix_itelescope_dirs.sh`
- Remove spaces from calibration file names using `ap_rename_files_with_spaces.sh`
- Generate an initial master bad pixel file from the master dark file
  using `ap_find_badpixel.py`.
- After processing real data with the initial master bad pixel file
  you will likely find additional bad pixels and/or columns that arise
  at dates or exposures different from those used for the calibration
  files. Careful inspection of the processed observation data (e.g. with
  `ds9`) can be used generate a user-defined bad pixel file (in yaml 
  format) which can then be used with `ap_find_badpixel.py` to generate
  an updated master bad pixel file. (The observation data can then be 
  reprocessed to remove the user-identified bad-pixels and columns.)

## Pipeline Processing

