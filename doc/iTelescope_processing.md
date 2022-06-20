# Processing images from iTelescope.net

## Processing Stages

The following table outlines conceptual pipeline processing stages
for astronomical FITS imagery, primarily those obtained using iTelescope.

| Stage                  | Stage Name            | Description                     | Python Class | Package Used      | Command Line Tool?       |
| ---------------------- | --------------------- | ------------------------------- | ------------ | ----------------- | ------------------------ |
| Calibration Generation | Build Master Cal      | Build master calibration files  | Yes          | astropy.ccdproc   | ap_combine_darks.py      |
| Calibration Generation | Calc Read Noise       | Calculate detector read noise   | Yes          |                   | ap_calc_read_noise.py    |
| Calibration Generation | Find Bad Pixels       | Generate a bad pixel map        | Yes          |                   | ap_find_badpix.py        |
| Calibration Generation | Find Bad Columns      | Find bad columns for badpix.yml | Yes          |                   | ap_auto_badcols.py       |
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
 
### Expected FITS Header Keywords 

The following FITS header keywords are required.

| Ap Class         | FITS Keywords       | Notes                            |
| ---------------- | ------------------- | -------------------------------- |
| ApAddMetadata    | (none)              |                                  |
| ApAstrometry     | IMG_FILE (a)        | Fom ApFindStars source list file |
| ApAutoBadcols    | (none)              |                                  |
| ApCalibrate      | EXPOSURE or EXPTIME |                                  |
|                  | GAIN or EGAIN (o)   |                                  |
| ApFindBadPixels  | (none)              |                                  |
| ApFindStars      | EXPOSURE            |                                  |

Notes:

- (a) Required, but generates by another stage of AstroPhotography processing
- (o) Optional, but not ideal if not present.


## Software Preparation

The `requirements.txt` file should handle all software requirements except
`Astromatic` software. However use of `Astrometry.net` requires a key
and that key is best stored in your local `astroquery` configuration file.

### Astrometry.net Key Configuration

Note that this is the `astroquery` config file, not the `astropy` config
file. Edit `~/.astropy/config/astroquery.cfg` and set `api_key` in 
the `[astrometry_net]` section.



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
- Determine whether the master dark files have had the bias values
  subtracted or not.
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

#### Generating An Initial Master Bad Pixel File

The following example demonstates how to generate an initial set of bad
pixel files. In this case for iTelescope T20 with the 2020-04 master files.

```bash
cd calibration-library/T20/Masters/Darks/2020-04/
python3 ~/git/AstroPhotography/src/AstroPhotography/scripts/ap_find_badpix.py -l DEBUG \
    Master_Dark_1_4008x2672_Bin1x1_Temp-15C_ExpTime900s.fit \
    Master_Badpix_1_4008x2672_Bin1x1_Temp-15C_ExpTime900s.fit 
python3 ~/git/AstroPhotography/src/AstroPhotography/scripts/ap_find_badpix.py -l DEBUG \
    Master_Dark_2_2004x1336_Bin2x2_Temp-15C_ExpTime900s.fit \
    Master_Badpix_2_2004x1336_Bin2x2_Temp-15C_ExpTime900s.fit 
```

Comparison of the resulting bad pixel files to the darks themselves show
that the default bad pixel detection works well on the isolated high data
value bad pixels, but misses some potentially problematic bad columns.

#### User-defined Bad Pixel File

These bad columns are most easily found by calibrating some real data 
with the initial badpix file, and inspecting the calibrated images for
artifacts.

The `ap_auto_badcol.py` script, which uses `ApAutoBadcols`, can be used
to find the majority of the most obvious bad columns and rows at thise 
stage. Its output can then be copy and pasted into a user-defined 
badpixel file (see template in `etc/user_badpixels.yml`). Rerunning
`ap_find_badpix.py` using both the master dark and the user-defined
bad pixels generates an updated master badpix file, which can then be
used to recalibrate your images. For example:

```bash
# copy template yml file to t20_user_badpixels_1_4008x2672_2020.yml

# automatically detect the worse bad columns in an initially calibrated
# image
ap_auto_badcol.py 20210708/cal-T20-davestrickland-CygnusLoop_x1_y1-20210708-220604-Ha-BIN1-E-300-001.fits

# copy output into t20_user_badpixels_1_4008x2672_2020.yml

# update the master bad pixel file
ap_find_badpix.py -l DEBUG \
    Master_Dark_1_4008x2672_Bin1x1_Temp-15C_ExpTime900s.fit \
    Master_Badpix_1_4008x2672_Bin1x1_Temp-15C_ExpTime900s.fit \
    --user_badpix t20_user_badpixels_1_4008x2672_2020.yml
```

You would then recalibrate your images using the updated cailbration
files. 

Inspecting those images will likely reveal other remaining bad columns
that are weaker, or only partially bad, that were missed by 
`ap_auto_badcols`. In my experience these are most easily visible in
long duration narrow-band images.

You should identify those by eye, and add them to the user bad pixel 
file. Then rerun `ap_find_badpix.py` and recalibrate your images,
and iterate this process until no obvious bad columns/rows/regions
are visible.

## Pipeline Processing

