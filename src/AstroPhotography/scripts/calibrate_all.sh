#!/usr/bin/env bash
#-----------------------------------------------------------------------
# 
# calibrate_all.sh [target] [itelescope] [skybg|noskybg] [noclean|clean]
# 
# Perform basic calibration on a set of iTelescope files using
# ap_calibrate.py
#
# - If 'clean' is specified existing files will be removed and regenerated.
#   By default 'noclean' is set so existing files are not rebuilt.
#
# Limitations:
# - Needs to be hardwired with the bias, dark, badpix and flat
#   calibration files to use.
# - Can not separately add to metadata without redoing entire calibration
# - clean option does not properly work with sky BG subtraction if some
#   files have already been calibrated.
#
# History:
# 2020-12-09 dks : Initial version using ap_calibrate.py
# 2020-12-20 dks : Added metadata update to support ap_find_stars.py
# 2021-01-16 dks : Add --fixcosmic to ap_calibrate call.
# 2021-02-21 dks : Make target a command line argument.
# 2021-02-27 dks : Add 2021-02-14 calibration, --dark_still_biased CLI flag.
# 2021-06-07 dks : Add option to calculate and subtract the sky background.
# 2021-07-22 dks : Work on better supporting other iTelescope telescopes.
# 2021-08-19 dks : Actually apply sky background correction.
# 2021-08-25 dks : Simpler calibration setup for Cygnus Loop mosaic.
#
#-----------------------------------------------------------------------
# Initialization

# iTelescope I have calibration data for.
p_supported=("T05" "T09" "T14" "T16" "T20" "T32" )

p_usage="$0 [target_name] [T05|T20] [skybg|noskybg] [noclean|clean]"

if [ $# -lt 3 ]; then
    echo "Error: Expecting at least 2 command line arguments, got $#"
    echo "  usage: $p_usage"
    exit 1
fi

# Target name. This is name used in the iTelescope file and directory
# names.
# Examples: ngc_6888, M82
p_targ="$1"
echo "Processing data for target: $p_targ"

# Name of the iTelescope we're processing. This is used to find the original
# data AND the calibration files.
# Note that we've switched to iTelescope's recent mode of having the telescope
# number in uppercase for both data and calibration.
p_telescope="$2"
# Check the telescope is in the list we expect.
if [[ ! " ${p_supported[@]} " =~ " ${p_telescope} " ]]; then
    echo "Error, telescope $p_telescope is not one of the expect/supported telescopes."
    echo "  Supported telescopes: ${p_supported[@]}"
    exit 1
fi

# Whether to calculate and subtract an estimate of the "sky" background.
# This may also be necessary when the bias/dark/flat correction leaves
# artifacts in the images.
if [ -z $3 ]; then
    echo "Error, specify one of skybg or noskybg."
    echo "  usage: $p_usage"
    exit 1
else
    if [[ $3 == "skybg" ]]; then
        echo "Sky background will be calculated and subtracted from the images."
        p_dosky=1
    else
        echo "No sky background estimation will be calculated."
        p_dosky=0
    fi
fi
    
# Clean run, remove existing outputs before running again.
if [ -z $4 ]; then
    p_clean=0
else
    if [[ "$4" == "clean" ]]; then
        p_clean=1
    else
        p_clean=0
    fi
fi

# Log files for entire runs
p_log='tmp_cal_all.log'
if [ -e $p_log ]; then
    rm $p_log
fi

p_start=$(date)
echo "$0 started at $p_start" | tee -a $p_log
echo "Running from" $(pwd) | tee -a $p_log
echo "Logging to $p_log" | tee -a $p_log
echo "Activating virtual environment." | tee -a $p_log
source ~/venv/astro38/bin/activate

# Scripts that do the work.
p_apcal="$HOME/git/AstroPhotography/src/AstroPhotography/scripts/ap_calibrate.py"
p_apmeta="$HOME/git/AstroPhotography/src/AstroPhotography/scripts/ap_add_metadata.py"
p_apskybg="$HOME/git/AstroPhotography/src/AstroPhotography/scripts/ap_measure_background.py"
p_apimarith="$HOME/git/AstroPhotography/src/AstroPhotography/scripts/ap_imarith.py"
echo "Checking that scripts exist..." | tee -a $p_log
p_err=0
for p_script in $p_apcal $p_apmeta $p_apskybg $p_apimarith; do
    if [ -e $p_script ]; then
        echo "  Found $p_script" | tee -a $p_log
    else
        echo "  Error: $p_script does not exist." | tee -a $p_log
        p_err=1
    fi
done
if [ $p_err -gt 0 ]; then
    exit 1
fi

# Common settings for the scripts.
p_loglevel='DEBUG'

# Filters to process:
p_filter_arr=("Red" "Green" "Blue" "Ha" "OIII" "SII" "Clear" "Luminance" "B" "V" "I" )
##p_filter_arr=("Red") # Testing purposes only.

# Calibration files to use
p_itel=/Users/dks/Tmp/iTelescopeScratch
p_data_dir=$p_itel/$p_telescope/$p_targ/
p_cal_dir=$p_itel/calibration-library/$p_telescope/Masters
declare -A p_flat_arr

echo "Checking that base-level data and calibration directories exist:" | tee -a $p_log
for p_dir in $p_data_dir $p_cal_dir; do
    if [ -d $p_dir ]; then
        echo "  Directory $p_dir found." | tee -a $p_log
    else
        echo "  Error, directory $p_dir not found." | tee -a $p_log
        p_err=1
    fi
done
if [ $p_err -gt 0 ]; then
    exit 1
fi

# The following is observation specific and doesn't work very well
# in this form in a shell script.

p_cygloop21=("CygnusLoop_x1_y1" "CygnusLoop_x1_y2" "CygnusLoop_x1_y3" "CygnusLoop_x1_y4")


p_cal_to_use="undefined"
if [[ $p_targ == "ngc_6888" ]]; then
    p_cal_to_use="2020-03"
    
    p_cal_date=2020-03
    
    # Masters
    p_mdark=$p_cal_dir/Darks/$p_cal_date/Master_Dark_1_2184x1472_Bin1x1_Temp-10C_ExpTime300s.fit
    p_mbadp=$p_cal_dir/Darks/$p_cal_date/Master_Badpix_1_2184x1472_Bin1x1_Temp-10C_ExpTime300s.fit
    p_mbias=$p_cal_dir/Bias/$p_cal_date/Master_Bias_1_2184x1472_Bin1x1_Temp-10C_ExpTime0ms.fit
    p_dark_still_biased="--dark_still_biased"
    
    
    # The flats are filter-specific. The associative array contains the file names,
    # but not the path. These names are true for 2020-03, but may not work for
    # other dates.                                                                   
    p_flat_dir=$p_cal_dir/Flats/$p_cal_date
    
    p_flat_arr["B"]="Master_Flat_B_1_B_2184x1472_Bin1x1_Temp-10C_ExpTime15s.fit"
    p_flat_arr["I"]="Master_Flat_I_1_I_2184x1472_Bin1x1_Temp-10C_ExpTime13s.fit"
    p_flat_arr["Blue"]="Master_Flat_Blue_1_Blue_2184x1472_Bin1x1_Temp-10C_ExpTime5s.fit"
    p_flat_arr["OIII"]="Master_Flat_OIII_1_OIII_2184x1472_Bin1x1_Temp-10C_ExpTime12s.fit"
    p_flat_arr["Clear"]="Master_Flat_Clear_1_Clear_2184x1472_Bin1x1_Temp-10C_ExpTime23s.fit"
    p_flat_arr["Red"]="Master_Flat_Red_1_Red_2184x1472_Bin1x1_Temp-10C_ExpTime4s.fit"
    p_flat_arr["Green"]="Master_Flat_Green_1_Green_2184x1472_Bin1x1_Temp-10C_ExpTime5s.fit"
    p_flat_arr["SII"]="Master_Flat_SII_1_SII_2184x1472_Bin1x1_Temp-10C_ExpTime67s.fit"
    p_flat_arr["Ha"]="Master_Flat_Ha_1_Ha_2184x1472_Bin1x1_Temp-10C_ExpTime46s.fit"
    p_flat_arr["V"]="Master_Flat_V_1_V_2184x1472_Bin1x1_Temp-10C_ExpTime5s.fit"


    
elif [[ $p_targ == "M82" ]]; then
    p_cal_to_use="2021-02-14"
    
    p_cal_date=2021-02-14
    
    # Masters
    p_mdark=$p_cal_dir/Darks/$p_cal_date/Master_Dark_1_2184x1472_Bin1x1_Temp-10C_ExpTime900s.fit
    p_mbadp=$p_cal_dir/Darks/$p_cal_date/Master_Badpix_1_2184x1472_Bin1x1_Temp-10C_ExpTime900s.fit
    p_mbias=$p_cal_dir/Bias/$p_cal_date/Master_Bias_1_2184x1472_Bin1x1_Temp-10C_ExpTime0ms.fit
    p_dark_still_biased=""
    
    # The flats are filter-specific. The associative array contains the file names,
    # but not the path. These names are true for 2020-03, but may not work for
    # other dates.                                                                   
    p_flat_dir=$p_cal_dir/Flats/$p_cal_date
    
    p_flat_arr["B"]="Master_Flat_B_1_B_2184x1472_Bin1x1_Temp-10C_ExpTime14s.fit"
    p_flat_arr["I"]="Master_Flat_I_1_I_2184x1472_Bin1x1_Temp-10C_ExpTime22s.fit"
    p_flat_arr["Blue"]="Master_Flat_Blue_1_Blue_2184x1472_Bin1x1_Temp-10C_ExpTime6s.fit"
    p_flat_arr["OIII"]="Master_Flat_OIII_1_OIII_2184x1472_Bin1x1_Temp-10C_ExpTime3s.fit"
    p_flat_arr["Clear"]="Master_Flat_Clear_1_Clear_2184x1472_Bin1x1_Temp-10C_ExpTime4s.fit"
    p_flat_arr["Red"]="Master_Flat_Red_1_Red_2184x1472_Bin1x1_Temp-10C_ExpTime4s.fit"
    p_flat_arr["Green"]="Master_Flat_Green_1_Green_2184x1472_Bin1x1_Temp-10C_ExpTime5s.fit"
    p_flat_arr["SII"]="Master_Flat_SII_1_SII_2184x1472_Bin1x1_Temp-10C_ExpTime3s.fit"
    p_flat_arr["Ha"]="Master_Flat_Ha_1_Ha_2184x1472_Bin1x1_Temp-10C_ExpTime3s.fit"
    p_flat_arr["V"]="Master_Flat_V_1_V_2184x1472_Bin1x1_Temp-10C_ExpTime6s.fit"
elif [[ $p_targ == "M81" ]]; then
    p_cal_to_use="2021-02-14"
    
    p_cal_date=2021-02-14
    
    # Masters
    p_mdark=$p_cal_dir/Darks/$p_cal_date/Master_Dark_1_2184x1472_Bin1x1_Temp-10C_ExpTime900s.fit
    p_mbadp=$p_cal_dir/Darks/$p_cal_date/Master_Badpix_1_2184x1472_Bin1x1_Temp-10C_ExpTime900s.fit
    p_mbias=$p_cal_dir/Bias/$p_cal_date/Master_Bias_1_2184x1472_Bin1x1_Temp-10C_ExpTime0ms.fit
    p_dark_still_biased=""
    
    # The flats are filter-specific. The associative array contains the file names,
    # but not the path. These names are true for 2020-03, but may not work for
    # other dates.                                                                   
    p_flat_dir=$p_cal_dir/Flats/$p_cal_date
    
    p_flat_arr["B"]="Master_Flat_B_1_B_2184x1472_Bin1x1_Temp-10C_ExpTime14s.fit"
    p_flat_arr["I"]="Master_Flat_I_1_I_2184x1472_Bin1x1_Temp-10C_ExpTime22s.fit"
    p_flat_arr["Blue"]="Master_Flat_Blue_1_Blue_2184x1472_Bin1x1_Temp-10C_ExpTime6s.fit"
    p_flat_arr["OIII"]="Master_Flat_OIII_1_OIII_2184x1472_Bin1x1_Temp-10C_ExpTime3s.fit"
    p_flat_arr["Clear"]="Master_Flat_Clear_1_Clear_2184x1472_Bin1x1_Temp-10C_ExpTime4s.fit"
    p_flat_arr["Red"]="Master_Flat_Red_1_Red_2184x1472_Bin1x1_Temp-10C_ExpTime4s.fit"
    p_flat_arr["Green"]="Master_Flat_Green_1_Green_2184x1472_Bin1x1_Temp-10C_ExpTime5s.fit"
    p_flat_arr["SII"]="Master_Flat_SII_1_SII_2184x1472_Bin1x1_Temp-10C_ExpTime3s.fit"
    p_flat_arr["Ha"]="Master_Flat_Ha_1_Ha_2184x1472_Bin1x1_Temp-10C_ExpTime3s.fit"
    p_flat_arr["V"]="Master_Flat_V_1_V_2184x1472_Bin1x1_Temp-10C_ExpTime6s.fit"
elif [[ $p_targ == "gaia" ]]; then
    p_cal_to_use="2021-02-14"
    p_true_target_name="HE 2300-0630"    
    p_cal_date=2021-02-14
    
    # Masters
    p_mdark=$p_cal_dir/Darks/$p_cal_date/Master_Dark_1_2184x1472_Bin1x1_Temp-10C_ExpTime900s.fit
    p_mbadp=$p_cal_dir/Darks/$p_cal_date/Master_Badpix_1_2184x1472_Bin1x1_Temp-10C_ExpTime900s.fit
    p_mbias=$p_cal_dir/Bias/$p_cal_date/Master_Bias_1_2184x1472_Bin1x1_Temp-10C_ExpTime0ms.fit
    p_dark_still_biased=""
    
    # The flats are filter-specific. The associative array contains the file names,
    # but not the path. These names are true for 2020-03, but may not work for
    # other dates.                                                                   
    p_flat_dir=$p_cal_dir/Flats/$p_cal_date
    
    p_flat_arr["B"]="Master_Flat_B_1_B_2184x1472_Bin1x1_Temp-10C_ExpTime14s.fit"
    p_flat_arr["I"]="Master_Flat_I_1_I_2184x1472_Bin1x1_Temp-10C_ExpTime22s.fit"
    p_flat_arr["Blue"]="Master_Flat_Blue_1_Blue_2184x1472_Bin1x1_Temp-10C_ExpTime6s.fit"
    p_flat_arr["OIII"]="Master_Flat_OIII_1_OIII_2184x1472_Bin1x1_Temp-10C_ExpTime3s.fit"
    p_flat_arr["Clear"]="Master_Flat_Clear_1_Clear_2184x1472_Bin1x1_Temp-10C_ExpTime4s.fit"
    p_flat_arr["Red"]="Master_Flat_Red_1_Red_2184x1472_Bin1x1_Temp-10C_ExpTime4s.fit"
    p_flat_arr["Green"]="Master_Flat_Green_1_Green_2184x1472_Bin1x1_Temp-10C_ExpTime5s.fit"
    p_flat_arr["SII"]="Master_Flat_SII_1_SII_2184x1472_Bin1x1_Temp-10C_ExpTime3s.fit"
    p_flat_arr["Ha"]="Master_Flat_Ha_1_Ha_2184x1472_Bin1x1_Temp-10C_ExpTime3s.fit"
    p_flat_arr["V"]="Master_Flat_V_1_V_2184x1472_Bin1x1_Temp-10C_ExpTime6s.fit"
elif [[ " ${p_cygloop21[@]} " =~ " ${p_targ} " ]]; then
    p_cal_to_use="2020-04"
    p_true_target_name="Cygnus Loop"
    p_cal_date=2020-04
    
    # Masters
    p_mdark=$p_cal_dir/Darks/$p_cal_date/Master_Dark_1_4008x2672_Bin1x1_Temp-15C_ExpTime900s.fit
    p_mbadp=$p_cal_dir/Darks/$p_cal_date/Master_Badpix_1_4008x2672_Bin1x1_Temp-15C_ExpTime900s.fit
    p_mbias=$p_cal_dir/Bias/$p_cal_date/Master_Bias_1_4008x2672_Bin1x1_Temp-15C_ExpTime0ms.fit
    p_dark_still_biased="--dark_still_biased"
    
    # The flats are filter-specific. The associative array contains the file names,
    # but not the path.                                                              
    p_flat_dir=$p_cal_dir/Flats/$p_cal_date
    
    p_flat_arr["V"]="Master_Flat_V_1_V_4008x2672_Bin1x1_Temp-15C_ExpTime22s.fit"
    p_flat_arr["Blue"]="Master_Flat_Blue_1_Blue_4008x2672_Bin1x1_Temp-15C_ExpTime25s.fit"
    p_flat_arr["OIII"]="Master_Flat_OIII_1_OIII_4008x2672_Bin1x1_Temp-15C_ExpTime4s.fit"
    p_flat_arr["Red"]="Master_Flat_Red_1_Red_4008x2672_Bin1x1_Temp-15C_ExpTime9s.fit"
    p_flat_arr["Green"]="Master_Flat_Green_1_Green_4008x2672_Bin1x1_Temp-15C_ExpTime81s.fit"
    p_flat_arr["SII"]="Master_Flat_SII_1_SII_4008x2672_Bin1x1_Temp-15C_ExpTime98s.fit"
    p_flat_arr["Ha"]="Master_Flat_Ha_1_Ha_4008x2672_Bin1x1_Temp-15C_ExpTime46s.fit"
    p_flat_arr["Luminance"]="Master_Flat_Luminance_1_Luminance_4008x2672_Bin1x1_Temp-15C_ExpTime9s.fit"

else
    echo "Error, calibration for target $p_targ is not defined."
    exit 1
fi
    
# Argh. This is fugly code that isn't going to work well with 
# different telescopes... We really need some auto-discovery style
# code, e.g. using CCDPROC.
echo "Using $p_cal_to_use calibration with $p_targ"
echo "Calibration files defined..."

#-----------------------------------------------------------------------
# Main work area

# Check filter-independent masters exist
echo "Checking master calibration files exist..." | tee -a $p_log
p_nmissing=0
for file in $p_mbias $p_mbadp $p_mdark; do
    if [ -f $file ]; then
        echo "  Found $file, OK" | tee -a $p_log
    else
        echo "  Warning, $file not found." | tee -a $p_log
        ((p_nmissing++))
    fi
done
for filter in $p_filter_arr; do
    file=$p_flat_dir/${p_flat_arr[$filter]}
    if [ -f $file ]; then
        echo "  Found $file, OK" | tee -a $p_log
    else
        echo "  Warning, $file not found." | tee -a $p_log
        ((p_nmissing++))
    fi    
done
if [ $p_nmissing -ne 0 ]; then
    echo "Error: $p_nmissing calibration files were not found." | tee -a $p_log
    echo "  Check the listed paths." | tee -a $p_log
    exit 2
fi

# Counters to store the number files successfully calibrated, those that were
# skipped, and those that failed to process.
p_nproc=0
p_nskip=0
p_nfail=0

# Array for file names and processing statuses.
p_fname_arr=()
p_fstat_arr=()

# Begin main processing loop.
echo "-----------------------------------------------------------------"  | tee -a $p_log
echo "About to loop through the following filters: ${p_filter_arr[@]}" | tee -a $p_log
for filter in ${p_filter_arr[@]}; do
    echo "Processing $filter filter images:" | tee -a $p_log
    
    # Name of master flat
    p_mflat=$p_flat_dir/${p_flat_arr[$filter]}

    p_raw_file_arr=( $(find . -name "raw-*-$filter-*.fit*" | xargs) )
    echo "  Found ${#p_raw_file_arr[@]} raw fits files." | tee -a $p_log
    for p_raw_file in ${p_raw_file_arr[@]}; do
        # Set up names for expected inputs and outputs ---------------------
    
        # Basename of input calibrated file with suffix removed
        p_base_name=$(basename ${p_raw_file%.*})
        
        # Name for the calibrated file. Note change to .fits from .fit
        p_cal_file=$(echo ${p_raw_file%.*} | sed -e 's/raw-/cal-/')'.fits'

        # If sky background subtraction is to be performed, the following
        # files are used: a background estimate, and a copy of the original
        # calibrated image (without sky background subtraction).
        p_skybg_file=$(echo ${p_cal_file%.*} | sed -e 's/cal-/skybg-/')'.fits'
        p_orig_cal_file=$(echo ${p_cal_file%.*} | sed -e 's/cal-/withskybg_skybg-/')'.fits'

        # Name for log file
        p_log_file=$(echo ${p_raw_file%.*} | sed -e 's/raw-/logcal-/')'.log'

        echo "  Raw file $p_raw_file:" | tee -a $p_log
        p_fname_arr+=( $p_raw_file )
                    
        # Does calibrated file already exist?
        if [ -e $p_cal_file ] && [ $p_clean -eq 0 ]; then
            # Skip calibration
            echo "    Calibrated file already exists: $p_cal_file" | tee -a $p_log
            ((p_nskip++))
            p_fstat_arr+=( "SKIP" )
            continue
        elif [ -e $p_cal_file ] && [ $p_clean -eq 1 ]; then
            # Wipe existing file and regenerate.
            echo "    Removing existing calibrated files and regenerating them..." | tee -a $p_log
            for tfile in $p_cal_file $p_log_file $p_orig_cal_file $p_skybg_file; do
                if [ -e $tfile ]; then
                    rm $tfile
                    echo "      Removed $tfile" | tee -a $p_log
                fi
            done

        else
            echo "    Generating calibrated file: $p_cal_file" | tee -a $p_log
        fi
        
        # Run ap_calibrate.py, performing bias subtraction, dark 
        # subtraction, flat fielding, bad pixel removal, cosmic ray
        # removal.
        $p_apcal $p_raw_file \
            $p_mbias $p_mdark $p_cal_file \
            --master_flat=$p_mflat \
            --master_badpix=$p_mbadp \
            $p_dark_still_biased --fixcosmic \
            --loglevel=DEBUG &> $p_log_file
        p_status=$?

        # Run ap_add_metadata.py
        if [ $p_status -eq 0 ]; then
            # Check whether p_true_target_name is set, in which case
            # we know that file name resolution will fail.
            if [ -n "$p_true_target_name" ]; then
                echo "" &>> $p_log_file
                $p_apmeta $p_cal_file --loglevel=DEBUG \
                    --target="$p_true_target_name" &>> $p_log_file
                p_add_status=$?
            else
                # Use file name based target name
                $p_apmeta $p_cal_file --loglevel=DEBUG &>> $p_log_file
                p_add_status=$?
            fi
            
            if [ $p_add_status -ne 0 ]; then
                p_status=$p_add_status
            fi
        fi 

        # Sky background estimation and subtraction.
        if [ $p_dosky -eq 1 ]; then
            # If clean is set, rebuild skybg and subtract.
            # Else check that the skybg and pre-skybg-subtracted file
            # exist, and only if so can we skip this step.
            
            if [ $p_clean -eq 1 ] || [ ! -e $p_skybg_file ]; then
                echo "      Performing sky background estimation on $p_cal_file" | tee -a $p_log
        
                # Copy the original calibrated file to $p_orig_cal_file,
                # then calculate sky background on that
                echo "        Creating copy of original calibrated file: $p_orig_cal_file" | tee -a $p_log
                echo "" &>> $p_log_file
                cp -v $p_cal_file $p_orig_cal_file &>> $p_log_file
                p_status=$?
            
                echo "        Generating sky background $p_skybg_file" | tee -a $p_log
                echo "" &>> $p_log_file
                $p_apskybg $p_orig_cal_file $p_skybg_file -l DEBUG \
                    --nbg_cols=16 --nbg_rows=32 &>> $p_log_file
                p_status=$?
            
                echo "        Subtracting sky background to (re)generate $p_cal_file" | tee -a $p_log
                echo "" &>> $p_log_file
                $p_apimarith $p_orig_cal_file SUB $p_skybg_file \
                    $p_cal_file -l DEBUG &>> $p_log_file
                p_status=$?
            else
                echo "      Skipping sky background estimation for $p_cal_file"
            fi
        fi

        if [ $p_status -eq 0 ]; then
            ((p_nproc++))
            p_fstat_arr+=( "OK" )
            echo "      Calibration succesful, logged to $p_log_file" | tee -a $p_log
        else
            ((p_nfail++))
            p_fstat_arr+=( "ERROR" )
            echo "      Calibration FAILED, see $p_log_file" | tee -a $p_log
        fi
    
        # Done for each raw file of that filter
    done

    # Done for each filter
done


# Print summary of number of files successfully processed, skipped or failed.
echo "" | tee -a $p_log
echo "-----------------------------------------------------------------" | tee -a $p_log
echo "Run summary:" | tee -a $p_log
p_num=${#p_fname_arr[@]}
idx=0
printf "%3s  %-90s  %-6s\n" "Idx" "FileName" "Status" | tee -a $p_log
while [ $idx -lt $p_num ]; do
    printf "%03d  %-90s  %-6s\n" $idx ${p_fname_arr[$idx]} ${p_fstat_arr[$idx]} | tee -a $p_log
    ((idx++))
done
echo "Succesfully calibrated $p_nproc files, skipped $p_nskip files, $p_nfail failed." | tee -a $p_log
echo "Master log file: $p_log" | tee -a $p_log

#-----------------------------------------------------------------------
# Clean up and exit.
echo "Deactivating virtual environment." | tee -a $p_log
deactivate

# All done
p_end=$(date)
echo "$0 started at $p_start, finished at $p_end" | tee -a $p_log
exit 0


#-----------------------------------------------------------------------
# Clean up and shut down

# Deactivate virtual env.
exit 0
