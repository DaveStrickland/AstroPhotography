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
| Calibration            | Cosmic Ray Reject     | Find and correct cosmic rays    | Yes          | astropy.ccdproc   | ap_fix_cosmic_rays.py    |
| Calibration            | Identify Flat Defects | Find problem areas in flats     | No           |                   | Not yet                  |
| Calibration            | Fill In Holes         | Fill in user-defined "holes"    | Yes          |                   | Not yet                  |
| Image Processing       | Find Stars And PSF    | Find stars, measure PSF         | Yes          | astropy.photutils | ap_find_stars.py         |
| Image Processing       | Astrometry            | Astrometric solution per image  | Yes          | astrometry.net    | ap_astrometry.py         |
| Image Processing       | Quality Filter        | Identify "poor quality" images  | Yes          |                   | ap_find_stars.py         |
|                        |                       | and summarize set of images     | Yes          |                   | ap_quality_summary.py    |
| Image Processing       | Image Arithmetic      | Arithmetic on 1 or 2 images     | Yes          |                   | ap_imarith.py            |
| Image Processing       | Background Estimation | Estimate sky background         | Yes          |                   | ap_measure_background.py |
| Image Processing       | Resample              | Resample an image to a WCS      | [4]          | astromatic swarp* | Not yet                  |
| Image Processing       | Continuum Subtract    | Continuum scaling/subtraction   | Not yet      |                   | Not yet                  |
| Image Processing       | Rasterize a FIS image | Generate a nice TIFF or PNG     | Yes          | [1, 2]            | ap_fits_rasterizer.py    | 
| Image Processing       | RGB Color Composite   | Make 3-color composites         | Yes          | [1, 2]            | ap_fits_rasterizer.py    | 
| Process All Images     | Calibrate All         | Calibrate all raw images        | Not yet      |                   | (calibrate_all.sh [3])   |
| Process All Images     | Navigate All          | Astrometry/WCS on all images    | Yes          |                   | ap_fits_rasterizer.py    |
| Process All Images     | Resample All          | Resample an image to a WCS      | Yes          | astromatic swarp* | ap_fits_rasterizer.py    |

Notes:

 - (1) Temporarily using an external non-python-based tool.
 - (2) Alternate method available using `astropy`
 - (3) A temporary bash implementation.
 - (4) The pipeline processor resamples multiple files to a common WCS, but there
   is currently no tool to resample a single FITS image to another WCS.
 
### Expected FITS Header Keywords 

**Note:** This section is incomplete.

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

*Conda/miniconda/miniforge users:* (recommended) The provided `conda` environment
file (`ap-env*.yml`) should handle all software requirements except
`Astromatic` software. 


*Pip-only users:* The `requirements.txt` file should handle all software requirements except
`Astromatic` software. 

Use of `Astrometry.net` requires a key
and that key is best stored in your local `astroquery` configuration file.


### Astrometry.net Key Configuration

Note that this is the `astroquery` config file, not the `astropy` config
file. Edit `~/.astropy/config/astroquery.cfg` and set `api_key` in 
the `[astrometry_net]` section.

If this configuration file is not present on your system, use the [information
here (astroquery 0.4.6 documentation)](https://astroquery.readthedocs.io/en/stable/#default-configuration-file) 
to copy the [file here](https://astroquery.readthedocs.io/en/stable/configuration.html) into
a new file at the listed location.

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

Calibration data from the iTelescope website should be placed in a
directory of your choosing. The telescope and date specfic directory 
structure from iTelesope should be retained under this master calibration
directory.

To point the `bash` shell scripts to the calibration data set an
environment variable `AP_CAL_DIR` pointing to your calibration parent directory.
For example:
```bash
> cd /scratch/iTelescopeScratch/calibration-library/
> ls 
T05  T09  T14  T16  T20  T32  T33
> tree -d T05
T05
├── Masters
│   ├── Bias
│   │   ├── 2017-07
│   │   ├── 2017-11
│   │   ├── 2018-06-27
│   │   ├── 2018-06-29
│   │   ├── 2018-10-01
│   │   ├── 2018-10-18
│   │   ├── 2019-01
│   │   ├── 2019-03
│   │   ├── 2019-08
│   │   ├── 2020-03
│   │   └── 2021-02-14
│   ├── Darks
│   │   ├── 2017-07
│   │   ├── 2017-11
│   │   ├── 2018-06-27
│   │   ├── 2018-06-29
│   │   ├── 2018-10-01
│   │   ├── 2018-10-18
│   │   ├── 2019-01
│   │   ├── 2019-03
│   │   ├── 2019-08
│   │   ├── 2020-03
│   │   └── 2021-02-14
│   └── Flats
│       ├── 2018-01
│       ├── 2018-09
│       ├── 2019-01
│       ├── 2019-02
│       ├── 2019-04
│       ├── 2019-08
│       ├── 2020-03
│       └── 2021-02-14
└── Raws
    └── Darks
        └── 2020-03

38 directories
export AP_CAL_DIR=/scratch/iTelescopeScratch/calibration-library/
```

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
- Master flats can have high spatial frequency artifacts that can cause
  abnormal patches of brighter or fainter than expected pixels in
  calibrated images that are not present in your raw data. See
  `patching_holes_from_bad_flats.ipynb`.

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

#### Recent Changes To iTelescope Calibration Data

*NB*: As of mid-to-late 2024 iTelescope does not always (a) release master dark
and/or master bias files, and (b) does not place any master calibration
files it does generate to the `Masters` directory tree from the `Raw` directory tree.

This forces us to do much of the raw-to-master processing ourselves.

An example script to process the raw calibration files for `T33` is shown below.

```bash
#!/bin/bash
#-----------------------------------------------------------------------
#
# bash make_t33_masters.sh <t33_cal_base_dir>
#
# @brief Generate master calibration files for T33 2024_06 epoch.
#
# Script to process iTelescope T33 calibration files downloaded on
# 2024-12-01 in the./calibration_library/T33/Raw_*/2024_06/raw_*
# folders into master bias, mast dark, master flat (actually already
# present in the "raw" flats directories), and master bad pixel files.
#
# Pre-requisites:
# Conda or similar python virtual environment with AstroPhotography 
# v0.6.0 or higher installed and active.
#
# @history 2024-12-03 dks : Initial coding 
#
#-----------------------------------------------------------------------
#

p_usage="$0 <t33_cal_base_dir> <clean|noclean>"

# Directories relative to the base directory for the telescope
p_raw_flatdir="./Raw_Flats/2024_06/raw_flats/"
p_master_flatdir="./Masters/Flats/2024_06/"
raw_flat_list=(Master_Flat_Blue_1_Blue_4096x4096_Bin1x1_Temp-15C_ExpTime4s.fit.gz \
    Master_Flat_Green_1_Green_4096x4096_Bin1x1_Temp-15C_ExpTime5s.fit.gz \
    Master_Flat_Ha_1_Ha_4096x4096_Bin1x1_Temp-15C_ExpTime20s.fit.gz \
    Master_Flat_Luminance_1_Luminance_4096x4096_Bin1x1_Temp-15C_ExpTime3s.fit.gz \
    Master_Flat_OIII_1_OIII_4096x4096_Bin1x1_Temp-15C_ExpTime17s.fit.gz \
    Master_Flat_Red_1_Red_4096x4096_Bin1x1_Temp-15C_ExpTime4s.fit.gz \
    Master_Flat_SII_1_SII_4096x4096_Bin1x1_Temp-15C_ExpTime27s.fit.gz)

p_raw_biasdir="./Raw_Bias/2024_06/raw_bias/"
p_master_biasdir="./Masters/Bias/2024_06/"
p_bin1_master_bias=master_bias_t33_202406_4096x4096_bin1.fits
p_bin2_master_bias=master_bias_t33_202406_2048x2048_bin2.fits

p_raw_darkdir="./Raw_Darks/2024_06/raw_darks/"
p_master_darkdir="./Masters/Darks/2024_06/"
p_bin1_master_dark=master_dark_t33_202406_4096x4096_bin1.fits
p_bin2_master_dark=master_dark_t33_202406_2048x2048_bin2.fits

# bad pixel file
p_bin1_badpix=$(echo $p_bin1_master_dark | sed -e 's/dark/badpix/')
p_bin2_badpix=$(echo $p_bin2_master_dark | sed -e 's/dark/badpix/')

# user-defined bad pixels and auto-badcols output
p_bin1_user_badpix="user_badpix_t33_202406_4096x4096_bin1.yml"
p_bin2_user_badpix="user_badpix_t33_202406_2048x2048_bin2.yml"

# plot of bad cols
p_bin1_user_badpix_plot="user_badcolplot_t33_202406_4096x4096_bin1.png"
p_bin2_user_badpix_plot="user_badcolplot_t33_202406_2048x2048_bin2.png"

# Scripts
p_script_dir=~/git/AstroPhotography/AstroPhotography/scripts
p_make_master_script=$p_script_dir'/'ap_combine_darks.py
p_find_badpix_script=$p_script_dir'/'ap_find_badpix.py
p_auto_badcol_script=$p_script_dir'/'ap_auto_badcol.py

for script in $p_make_master_script; do
    if [ ! -e $script ]; then
        echo "Error, cannot find script $script"
        exit 4
    fi
done

make_op_dir(){
    if  [ ! -d "$1" ]; then
        echo "  Creating directory $1"
        mkdir -p "$1"
    else
        echo "  Directory $1 already exists."
    fi
}

do_flats(){
    # usage do_flats <base_dir> <raw_flat_dir> <master_flat_dir>
    local l_basedir="$1"
    local l_rawdir="$2"
    local l_masterdir="$3"
    cd $l_basedir
    make_op_dir "$l_masterdir"
    for file in ${raw_flat_list[@]}; do
        rsync -av "$l_rawdir"/$file "$l_masterdir"
    done
    echo "  Finished processing flats."
}

do_combine(){
    # usage do_combine <base_dir> <raw_bias_dir> <master_bias_dir> <bin1_master_bias> <bin2_master_bias>
    local l_basedir="$1"
    local l_rawdir="$2"
    local l_masterdir="$3"
    local l_bin1master="$4"
    local l_bin2master="$5"
    cd $l_basedir
    make_op_dir "$l_masterdir"
    
    if [ $p_do_clean -eq 1 ]; then
        for file in $l_masterdir/$l_bin1master  $l_masterdir/$l_bin2master; do 
            if [ -e $file ]; then
                echo "  Removing existing file $file"
                rm $file
            fi
        done
    fi
    
    # Process Bin2
    local l_masterfile=$l_masterdir/$l_bin2master
    if [ ! -e $l_masterfile ]; then
        $p_make_master_script -l DEBUG $l_rawdir $l_masterfile \
            --exclude="Apo*[Bb]in1.fit" --telescop="iTelescope 33"
        if [ $? -ne 0 ]; then
            echo "Error, ap_combine_darks failed"
            exit 10
        elif [ ! -e $l_masterfile ]; then
            echo "Error, failed to create $l_masterfile"
            exit 11
        else
            echo "  Successfully created $l_masterfile"
        fi
    else
        echo "  Skipping regeneration of $l_masterfile"
    fi

    # Process Bin1
    local l_masterfile=$l_masterdir/$l_bin1master
    if [ ! -e $l_masterfile ]; then
        $p_make_master_script -l DEBUG $l_rawdir $l_masterfile \
            --exclude="Apo*[Bb]in2.fit" --telescop="iTelescope 33"
        if [ $? -ne 0 ]; then
            echo "Error, ap_combine_darks failed"
            exit 10
        elif [ ! -e $l_masterfile ]; then
            echo "Error, failed to create $l_masterfile"
            exit 11
        else
            echo "  Successfully created $l_masterfile"
        fi
    else
        echo "  Skipping regeneration of $l_masterfile"
    fi
    
    echo "  Finished combining files into a master."
}

# Check number of command line arguments
if [ $# -ne 2 ]; then
    echo "Error, expecting 2 command line argument, got $#"
    echo "  usage: $p_usage"
    exit 1
fi

# Base directory containing T33-specific calibration files, both raw
# and masters.
p_basedir=$1
if [ ! -d $p_basedir ]; then
    echo "Error, input base directory for T33 calibration files $p_basedir does not exist or is not a directory."
    exit 2
else
    echo "Base directory containing T33 calibration files: $p_basedir"
fi

# Clean or don't clean. If clean then delete existing outputs.
p_do_clean=0
if [[ "$2" == "clean" ]]; then
    p_do_clean=1
    echo "Clean specified: existing outputs will be deleted and recreated."
else
    echo "Existing master files will be retained and not regenerated."
fi

p_tasks=("flats" "bias" "darks")
for task in ${p_tasks[@]}; do
    case $task in
        "flats")
            echo "Doing flats"
            do_flats $p_basedir $p_raw_flatdir $p_master_flatdir            
            ;;
        "darks")
            echo "Doing darks"
            do_combine $p_basedir $p_raw_darkdir $p_master_darkdir $p_bin1_master_dark $p_bin2_master_dark
            
            # initial badpixel computation without user bad pix/bad col/bad row etc
            pushd $p_master_darkdir
            $p_find_badpix_script $p_bin1_master_dark $p_bin1_badpix \
                --loglevel=DEBUG --sigma=4.0
            if [ $? -ne 0 ]; then
                echo "Error, ap_find_badpix.py failed on $p_bin1_master_dark"
                exit 8
            fi    
            
            $p_find_badpix_script $p_bin2_master_dark $p_bin2_badpix \
                --loglevel=DEBUG --sigma=4.0 
            if [ $? -ne 0 ]; then
                echo "Error, ap_find_badpix.py failed on $p_bin2_master_dark"
                exit 8
            fi 
            
            echo "About to run: $p_auto_badcol_script $p_bin1_master_dark --loglevel=DEBUG --user_badcol_file=$p_bin1_user_badpix --sigma=5.0 --plot_stats=$p_bin1_user_badpix_plot"
                            
            $p_auto_badcol_script $p_bin1_master_dark --loglevel=DEBUG \
                --user_badcol_file=$p_bin1_user_badpix --sigma=5.0 \
                --plot_stats=$p_bin1_user_badpix_plot
            if [ $? -ne 0 ]; then
                echo "Error, ap_auto_badcol.py failed on $p_bin1_master_dark"
                exit 16
            fi 
            
            echo "About to run: $p_auto_badcol_script $p_bin2_master_dark --loglevel=DEBUG --user_badcol_file=$p_bin2_user_badpix --sigma=5.0 --plot_stats=$p_bin2_user_badpix_plot"
            
            $p_auto_badcol_script $p_bin2_master_dark --loglevel=DEBUG \
                --user_badcol_file=$p_bin2_user_badpix --sigma=5.0 \
                --plot_stats=$p_bin2_user_badpix_plot
            if [ $? -ne 0 ]; then
                echo "Error, ap_auto_badcol.py failed on $p_bin2_master_dark"
                exit 16
            fi 
                      
            # Then add the badcols back into the badpix.
            echo "About to run: " $p_find_badpix_script $p_bin1_master_dark $p_bin1_badpix \
                --loglevel=DEBUG --sigma=4.0 --user_badpix=$p_bin1_user_badpix
            $p_find_badpix_script $p_bin1_master_dark $p_bin1_badpix \
                --loglevel=DEBUG --sigma=4.0 --user_badpix=$p_bin1_user_badpix
            if [ $? -ne 0 ]; then
                echo "Error, ap_find_badpix.py failed on $p_bin1_master_dark"
                exit 8
            fi    
            
            echo "About to run: "  $p_find_badpix_script $p_bin2_master_dark $p_bin2_badpix \
                --loglevel=DEBUG --sigma=4.0 --user_badpix=$p_bin2_user_badpix
            $p_find_badpix_script $p_bin2_master_dark $p_bin2_badpix \
                --loglevel=DEBUG --sigma=4.0  --user_badpix=$p_bin2_user_badpix
            if [ $? -ne 0 ]; then
                echo "Error, ap_find_badpix.py failed on $p_bin2_master_dark"
                exit 8
            fi 
                        
            popd
            ;;
        "bias")
            echo "Doing bias"
            do_combine $p_basedir $p_raw_biasdir $p_master_biasdir $p_bin1_master_bias $p_bin2_master_bias
            ;;
        *)
            echo "Error, unexpected task $task. Ignoring it and carrying on..."
            ;;
    esac
done

#-----------------------------------------------------------------------
#
#
exit 0

```

## Pipeline Processing

Pipeline processing currently described in the Jupyter notebook `itelescope_premium_the_easier_way.ipynb`
that can be found in the `AstroPhotography/notebooks/` sub-directory.

