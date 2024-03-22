#!/usr/bin/env bash
#-----------------------------------------------------------------------
# 
# resample_all.sh [target] [median|wgtavg|sum] [noclean|clean] [verbose|quiet]
# 
# Use swarp to resample a set of images to the same footprint and pixel
# scale.
#
# Assumes:
# That swarp in installed
# That astropy is installed and that fitsheader is called fitsheader
#
# Limitations:
# - swarp does not recognize SIP, and "eats in" from gaps in the data.
# - Script does not carry over useful header info into resampled image
#   or provide estimate of rough net exposure.
# - Bug in summary function, gets the wrong band if some bands were skipped.
# - clean option does work correctly.
#
# History:
# 2021-01-22 dks : Initial version begun.
# 2021-01-31 dks : Change file naming to reflect weighted avg or median
# 2021-02-09 dks : Add preliminary version of sum mode as well.
# 2021-02-28 dks : Added target
# 2021-08-27 dks : Added CygnusLoop, fix some directory creation bugs.
# 2024-01-30 dks : Comment out use of virtual environment.
# 2024-02-23 dks : Add M101sn, fix exposure time checking. Add verbose.
#
#-----------------------------------------------------------------------
# Initialization

# Timing function from dks_unixtime.sh
function dks_time(){
    # No arguments, return seconds.nanoseconds since unix epoch
    # 2 arguments, subtract the first time from the second time.
    if [[ $# -eq 0 ]]; then
        # Unix time in seconds and nanoseconds
        date +%s.%N
    elif [[ $# -eq 2 ]]; then
        # Time difference
        echo $1 $2 | awk '{printf"%.3f\n",$2-$1}'
    else
        echo "Error: Incorrect number of command line arguments."
        echo "  usage: $p_usage1"
        echo "  OR"
        echo "  usage: $p_usage2"
        exit 1
    fi
}

# Main script
p_usage="$0 [target] [median|wgtavg] [noclean|clean] [verbose|quiet]"

if [ $# -lt 4 ]; then
    echo "Error: Expecting at least 4 command line argument(s), got $#"
    echo "  usage: $p_usage"
    exit 1
fi

# Target
p_target=$1

# Median, weighted average, or simple sum co-addition?
if [[ "$2" == "median" ]]; then
    p_addmode=0
    p_add_str="median"
    p_combine_type="MEDIAN"
elif [[ "$2" == "wgtavg" ]]; then
    p_addmode=1
    p_add_str="wgtavg"
    p_combine_type="WEIGHTED"
elif [[ "$2" == "sum" ]]; then
    p_addmode=2
    p_add_str="sum"
    p_combine_type="SUM"
else
    echo "Error, unexpected co-addition mode supplied: $1"
    echo "  Expecting one of: median wgtavg sum"
    exit 2
fi
echo "Co-addition of images will use swarp COMBINE_TYPE $p_combine_type"

# Clean run, remove existing outputs before running again.
if [ -z $3 ]; then
    p_clean=0
else
    if [[ "$3" == "clean" ]]; then
        p_clean=1
    else
        p_clean=0
    fi
fi

# Verbose mode: echo swarp command to stdout
if [ -z $4 ]; then
    p_verbose=0
else
    if [[ "$4" == "verbose" ]]; then
        p_verbose=1
    else
        p_verbose=0
    fi
fi

# Log files for entire runs. Log file must have full path.
p_log=$(pwd)'/tmp_resamplecd_all.log'
if [ -e $p_log ]; then
    rm $p_log
fi

p_start=$(date)
p_start_time=$(dks_time)
echo "$0 started at $p_start" | tee -a $p_log
echo "Running from" $(pwd) | tee -a $p_log
echo "Logging to $p_log" | tee -a $p_log
##echo "Activating virtual environment." | tee -a $p_log
##source ~/venv/astro/bin/activate

# Name of swarp executable
p_swarp=swarp

# Name of astropy fitsheader. On macports can be fitsheader-3.7
p_fitsheader=fitsheader

# Check that the programs we need are present.
for p_exe in $p_swarp $p_fitsheader; do
    which $p_exe >& /dev/null
    if [ $? -ne 0 ]; then
        echo "Error, $p_exe not found in your path." | tee -a $p_log
        exit 1
    else
        echo "Using $p_exe from" $(which $p_exe) | tee -a $p_log
    fi
done

# General settings
p_gain_kw=EGAIN
p_oversampling=4
p_res_type=LANCZOS3
p_projection=TAN
p_fileinfo=Y
p_nopenmax=128
p_fscalastro_type=VARIABLE
p_delete_tmp=Y
p_tmpfile_dir=TmpSwarp


# Target specific info
if [[ $p_target == "n6888" ]]; then
    p_targ="n6888"
    p_ra=303.0272359
    p_dec=38.3549333
    p_pixscale=1.8  # 1920x1.8 = 0.96x0.54 deg
    p_image_size=1920,1080 # full image 2560,1920
elif [[  $p_target == "M82" ]] ; then
    p_targ="m82"
    p_ra=148.9684583
    p_dec=69.6797028
    p_pixscale=1.8
    p_image_size=2560,1920
elif [[  $p_target == "M81" ]] ; then
    p_targ="m81"
    p_ra=148.8882
    p_dec=69.065295
    p_pixscale=1.8
    p_image_size=2560,1920
elif [[  $p_target == "M82Group" ]] ; then
    p_targ="m82group"
    p_ra=148.9255
    p_dec=69.385417
    p_pixscale=1.8
    p_image_size=3000,3000
elif [[  $p_target == "CygnusLoop" ]] ; then
    p_targ="CygnusLoop"
    p_ra=312.75
    p_dec=30.67
    p_pixscale=3.6
    p_image_size=5000,5000
elif [[  $p_target == "M101sn" ]] ; then
    p_targ="M101_Supernova_2023ixf"
    p_ra=210.8022671
    p_dec=54.3489500
    p_pixscale=0.62
    p_image_size=4096,4096
else
    echo "Error, unexpected target $p_target" | tee -a $p_log
    echo "  This could be recoded to allow for general targets." | tee -a $p_log
    exit 3
fi

# Additional mode-specific parameter changes
if [ $p_addmode -eq 0 ]; then
    # median
    true
elif [ $p_addmode -eq 1 ]; then
    # weighted average
    true
elif [ $p_addmode -eq 2 ]; then
    # Sum. Technically this really only works if the output image has the same
    # pixel size as the input.
    p_fscalastro_type=NONE
    #true
else
    echo "Error, unexpected co-addition mode: p_addmode=$p_addmode" | tee -a $p_log
    exit 2
fi


# Filters to process:
p_filter_arr=( "Red" "Green" "Blue" "Ha" "OIII" "SII" "B" "V" "I" "Clear" "Lum" "Luminance" )
#p_filter_arr=( "Red" ) # Testing purposes only.

# Directory for input navigated images:
p_navdir=./NavigatedImages

# Directory for output resampled co-added images.
p_resdir=./Resampled

#-----------------------------------------------------------------------
# Main work area

# Check input and output directories exist. We can create the resampled
# dir if it exists but the navigated images directory must exist.
if [ ! -d $p_navdir ]; then
    echo "Error, can not find $p_navdir" | tee -a $p_log
    echo "  Current path: $cwd" | tee -a $p_log
    exit 2
else
    echo "Will read navigated images from $p_navdir" | tee -a $p_log
fi
if [ ! -d $p_resdir ]; then
    echo "Creating $p_resdir for resampled images." | tee -a $p_log
    mkdir $p_resdir
else
    echo "Reusing existing $p_resdir for resampled images." | tee -a $p_log
fi

# Array for file names and processing statuses.
p_filtname_arr=()
p_fname_arr=()
p_num_arr=()
p_time_arr=()
p_fstat_arr=()
p_texp_arr=()

# Original directory
p_odir=$(pwd)

# Begin main processing loop.
echo "---------------------------------------------------------------------"  | tee -a $p_log
echo "About to loop through the following filters: ${p_filter_arr[@]}" | tee -a $p_log
for filter in ${p_filter_arr[@]}; do
    echo "----------------------------------------------------------------------" | tee -a $p_log
    echo "Processing $filter filter images:" | tee -a $p_log 
    p_filtname_arr+=( $filter )
    
    # Change directories
    cd $p_resdir

    # Tempfile dir
    if [ ! -e $p_tmpfile_dir ]; then
        echo "  Temporary files will be written to new directory $p_tmpfile_dir" | tee -a $p_log
        mkdir -p $p_tmpfile_dir
    else
        echo "  Using existing temporary file directory $p_tmpfile_dir" | tee -a $p_log
    fi


    # Temporary file
    p_tmp=$(mktemp -t tmp_XXXX)

    # For building the fscale string.
    p_fscale_str=""
    
    # Search for navigated file in $p_navdir, excluding its subdirectories
    # because I store old/bad versions of the navigated files in those.
    p_nav_file_arr=( $(find ../$p_navdir -maxdepth 1 -name "nav*-$filter-*.fit*" | xargs) )
    p_num_imgs=${#p_nav_file_arr[@]}
    p_num_arr+=( $p_num_imgs )
    p_texp=0 # Total exposure
    echo "  Found $p_num_imgs navigated fits files." | tee -a $p_log
    
    if [ $p_num_imgs -eq 0 ]; then
        echo "  Moving on to next filter..." | tee -a $p_log
        p_fname_arr+=( "None" )
        p_fstat_arr+=( 0 )
        p_texp_arr+=( 0 )
        p_time_arr+=( 0 )
        cd $p_odir
        continue
    else
        p_tmpfile_fullpath=$(pwd)/$p_tmpfile_dir
    fi
    
    # Skip processing if p_clean==0 and the resampled file exists
    if [ $p_clean -eq 0 ]; then
        if [ -e $p_res_img ]; then
            echo "  Skipping processing, as resampled co-added image $p_res_img exists" | tee -a $p_log
            echo "  and noclean specified on the command line" | tee -a $p_log
            p_fstat_arr+=( 0 )
            p_time_arr+=( 0 )
            p_fname_arr+=( $p_res_img )
            cd $p_odir
            continue
        fi
    fi
    
    for p_nav_file in ${p_nav_file_arr[@]}; do
        # Set up names for expected inputs and outputs ---------------------
    
        # Basename of input calibrated file.
        p_base_name=$(basename ${p_nav_file})
        
        # astropy's fitsheader tool has some peculiarities...
        # EXPOSURE or EXPTIME or ONTIME or LIVETIME
        if [ -e $p_tmp ]; then
            rm $p_tmp
        fi
        $p_fitsheader $p_nav_file --extension=0 --table ascii >& $p_tmp
        if [ $? -ne 0 ]; then
            # fitsheader does not use normal unix return status codes, but
            # an untreated exception might generate a non-zero return code.
            echo "Error, fitsheader threw an exception with $p_nav_file." | tee -a $p_log
            exit 20
        fi
        p_nexpkw=$(grep -c "EXPOSURE\|EXPTIME\|ONTIME\|LIVETIME" $p_tmp)
        if [ $p_nexpkw -eq 0 ]; then
            echo "Error, could not find EXPOSURE, EXPTIME, ONTIME or LIVETIME keywords in $p_nav_file." | tee -a $p_log
            exit 20
        fi
        # Ascii table format is "filenamepath extensionnum keyword value"
        p_exp=$(grep "EXPOSURE\|EXPTIME\|ONTIME\|LIVETIME" $p_tmp | sort -k 3,3 | head -1 | awk '{print $4}')
        p_fscale=$(echo $p_exp | awk '{printf"%.6f\n",1/$1}')
        p_texp=$(echo $p_texp $p_exp | awk '{print $1+$2}')
        
        printf "  File %-40s EXPOSURE %8.3f FSCALE %8.6f \n" $p_base_name \
            $p_exp $p_fscale | tee -a $p_log
        
        # Build fscale string
        p_fscale_str=$(echo $p_fscale_str","$p_fscale)
        
        # Done for each navigated file of that filter
    done
    
    # If sum mode, don't flux scale. Otherwise, trim off first "," in fscale_str
    if [ $p_addmode -eq 2 ]; then
        p_fscale_str=1.0
    else
        p_fscale_str=$(echo "$p_fscale_str" | sed -e 's/,//')
    fi
    echo "  Flux scalings: $p_fscale_str" | tee -a $p_log
    
    # Convert total exposure from seconds to minutes
    p_texp=$(echo $p_texp | awk '{print $1/60.0}')
    echo "  Total exposure in $filter filter:" $p_texp "minutes." | tee -a $p_log
    echo ""  | tee -a $p_log
    p_texp_arr+=( $p_texp )
    
    # Name of output files: target filter size (resampled or weights) co-add_type
    p_res_img=$p_targ'_'$filter'_'$(echo $p_image_size | sed -e 's/,/x/')'_resamp_'$p_add_str'.fits'
    p_wgt_img=$p_targ'_'$filter'_'$(echo $p_image_size | sed -e 's/,/x/')'_weights_'$p_add_str'.fits'
    p_fname_arr+=( $p_res_img )
    
    echo "  Beginning swarp for filter $filter at" $(date) | tee -a $p_log
    p_res_start_time=$(dks_time)
    if [ $p_verbose -eq 1 ]; then
        echo "Running the following command: "  $p_swarp ${p_nav_file_arr[@]} \
        -SUBTRACT_BACK N -COMBINE_TYPE $p_combine_type \
        -FSCALASTRO_TYPE $p_fscalastro_type -VERBOSE_TYPE FULL \
        -NOPENFILES_MAX p_nopenmax -GAIN_KEYWORD $p_gain_kw \
        -GAIN_DEFAULT 1.0 -FSCALE_DEFAULT $p_fscale_str \
        -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE $p_pixscale \
        -CENTER_TYPE MANUAL -CENTER "$p_ra,$p_dec" \
        -RESAMPLING_TYPE $p_res_type -OVERSAMPLING $p_oversampling \
        -IMAGE_SIZE $p_image_size -PROJECTION_TYPE $p_projection \
        -IMAGEOUT_NAME $p_res_img \
        -WRITE_FILEINFO $p_fileinfo -WRITE_XML N \
        -DELETE_TMPFILES $p_delete_tmp -RESAMPLE_DIR $p_tmpfile_dir \
        -WEIGHTOUT_NAME $p_wgt_img  | tee -a $p_log
    fi
    $p_swarp ${p_nav_file_arr[@]} \
        -SUBTRACT_BACK N -COMBINE_TYPE $p_combine_type \
        -FSCALASTRO_TYPE $p_fscalastro_type -VERBOSE_TYPE FULL \
        -NOPENFILES_MAX p_nopenmax -GAIN_KEYWORD $p_gain_kw \
        -FSCALE_DEFAULT $p_fscale_str \
        -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE $p_pixscale \
        -CENTER_TYPE MANUAL -CENTER "$p_ra,$p_dec" \
        -RESAMPLING_TYPE $p_res_type -OVERSAMPLING $p_oversampling \
        -IMAGE_SIZE $p_image_size -PROJECTION_TYPE $p_projection \
        -IMAGEOUT_NAME $p_res_img \
        -WRITE_FILEINFO $p_fileinfo -WRITE_XML N \
        -DELETE_TMPFILES $p_delete_tmp -RESAMPLE_DIR $p_tmpfile_dir \
        -WEIGHTOUT_NAME $p_wgt_img >& $p_tmp
    p_status=$?
    cat $p_tmp >> $p_log
    if [ $p_status -ne 0 ]; then
        echo "Error, $p_swarp failed with status $p_status" | tee -a $p_log
        
    fi
    p_fstat_arr+=( $p_status )
    
    # Check that resampled image and weights exist.
    if [ -e $p_res_img ]; then
        echo "  Generated resampled co-added image $p_res_img" | tee -a $p_log
    else
        echo "  Warning, could not find resampled co-added image $p_res_img" | tee -a $p_log
    fi
    if [ -e $p_wgt_img ]; then
        echo "  Generated co-added weights image $p_wgt_img" | tee -a $p_log
    else
        echo "  Warning, could not find co-added weights image $p_wgt_img" | tee -a $p_log
    fi
    
    # Timing info. Might be useful if we switch parameters or even the
    # resampling method used.
    p_res_end_time=$(dks_time)
    p_res_time=$(dks_time $p_res_start_time $p_res_end_time)
    p_time_arr+=( $p_res_time )
    echo "  Swarp for filter $filter completed at" $(date)", took $p_res_time seconds." | tee -a $p_log

    # Remove temp file
    if [ -e $p_tmp ]; then
        rm $p_tmp
    fi

    cd $p_odir
    echo ""  | tee -a $p_log
    # Done for each filter
done


# Print summary of number of files successfully processed, skipped or failed.
echo "" | tee -a $p_log
echo "--------------------------------------------------------------------------" | tee -a $p_log
echo "Run summary:" | tee -a $p_log
p_num=${#p_fname_arr[@]}
idx=0
printf "%-10s  %-65s  %-6s  %8s  %8s  %-6s\n" "Filter" \
    "Resampled Output" "Ninput" "TexpMin" "SwrpTSec" "Status"| tee -a $p_log
while [ $idx -lt $p_num ]; do
    printf "%-10s  %-65s  %6d  %8.2f  %8.2f  %6d\n" ${p_filtname_arr[$idx]} \
        ${p_fname_arr[$idx]} ${p_num_arr[$idx]} ${p_texp_arr[$idx]} \
        ${p_time_arr[$idx]} ${p_fstat_arr[$idx]} | tee -a $p_log
    ((idx++))
done
echo "Master log file: $p_log" | tee -a $p_log

#-----------------------------------------------------------------------
# Clean up and exit.
##echo "Deactivating virtual environment." | tee -a $p_log
##echo "" | tee -a $p_log
##deactivate

# All done
p_scriptname=$(basename $0)
p_tmp_mb=$(du -sm $p_tmpfile_fullpath | cut -f 1)
p_end=$(date)
p_end_time=$(dks_time)
p_run_time=$(dks_time $p_start_time $p_end_time)
echo "" | tee -a $p_log
echo "--------------------------------------------------------------------------" | tee -a $p_log
echo "  Temporary files occupy ${p_tmp_mb} MB in $p_tmpfile_fullpath" | tee -a $p_log
echo "  Commands were logged to $p_log" | tee -a $p_log
echo "  Resampled images were written to $p_resdir" | tee -a $p_log
echo "$p_scriptname finished at $p_end, run time $p_run_time seconds." | tee -a $p_log
exit 0
