#!/usr/bin/env bash
#-----------------------------------------------------------------------
# 
# composite_all.sh [prefix] [suffix] [one or more color selectors from "sho", "rgb"]
# 
# Use stiff to form color composite images from separate fits files.
# For example, given a directory with the following images...
#   n6888_Blue_2560x1920_resamp.fits
#   n6888_OIII_2560x1920_resamp.fits
#   n6888_Green_2560x1920_resamp.fits
#   n6888_Red_2560x1920_resamp.fits
#   n6888_Ha_2560x1920_resamp.fits
#   n6888_SII_2560x1920_resamp.fits
# ...running this script with "n6888" as prefix, "2560x1920_resamp.fits"
# as suffix, and "sho" as the band selector will create a variety of
# 3-color composite images. The different outputs will share the same
# input images and color selections, but have different sets of values
# for optional stiff command line flags.
# 
# Color Selectors:
# - sho:     SII as red, Ha as green, OIII as blue.
# - rgb:     Red as red, Green as green, Blue as blue.
# - hgb:     Ha as red, Green as green, Blue as blue
#
# Assumes:
# That stiff in installed.
# All the images have the same size and alignment, i.e. were navigated
# and resampled to the same WCS.
#
# Limitations:
# - Rather clunky way of determining which images to use, because I
#   expect to have multiple different resampled versions of a given
#   filter present in a given directory.
# - As of 2024 an RPM version of stiff is not available for Fedora 39,
#   and the github source fails to compile and would require a lot of
#   modification to get compiling.
#
# History:
# 2021-01-26 dks : Initial version begun.
# 2021-08-28 dks : Added hgb color selections.
# 2024-01-30 dks : Note stiff no longer available via RPM and that the
#                  github version fails to compile.  
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
p_usage="$0 [prefix] [suffix] [one or more color selectors from 'sho', 'rgb']"

if [ $# -lt 3 ]; then
    echo "Error: Expecting at least 3 command line arguments, got $#"
    echo "  usage: $p_usage"
    exit 1
fi

# File prefix, the part before the filter specification
p_pref=$1

# File suffix, the part after the filter specification
p_suff=$2

# Color selection array
p_arg_arr=( "$@" ) # This is args with $0 stripped off, so indices go down by 1
p_colsel_arr=( )
i=2
while [ $i -lt $# ]; do
    p_colsel_arr+=( ${p_arg_arr[$i]} )
    (( i++ ))
done

# Log files for entire runs
p_log='tmp_composite_all.log'
if [ -e $p_log ]; then
    rm $p_log
fi

p_start=$(date)
p_start_time=$(dks_time)
echo "$0 started at $p_start" | tee -a $p_log
echo "Running from" $(pwd) | tee -a $p_log
echo "Logging file:     $p_log" | tee -a $p_log
echo "File prefix:      $p_pref" | tee -a $p_log
echo "Fill suffix:      $p_suff" | tee -a $p_log
echo "Color selections:" ${p_colsel_arr[@]} | tee -a $p_log

# Name of stiff executable
p_stiff=stiff

# Check that the programs we need are present.
for p_exe in $p_stiff; do
    which $p_exe >& /dev/null
    if [ $? -ne 0 ]; then
        echo "Error, $p_exe not found in your path."
        exit 1
    else
        echo "Using $p_exe from" $(which $p_exe)
    fi
done

# Allowed 3-color combination for color selection
p_allowed_colsel=( 'sho' 'rgb' 'hgb' )

#-----------------------------------------------------------------------
# Main work area

# Counters to store the number files successfully calibrated, those that were
# skipped, and those that failed to process.
p_nproc=0
p_nskip=0
p_nfail=0

# Array for file names and processing statuses.
p_fname_arr=()
p_time_arr=()
p_fstat_arr=()

# Suffix without file extension
p_short_suff=${p_suff%.*}

# Original directory
p_odir=$(pwd)

# Temporary file
p_tmp=$(mktemp -t tmp_XXXX)

# Gamma, saturation, other options to loop over...
# The labels are the strings to be used in the output file names.
p_color_sat_arr=( 1.0 1.5 2.0 )
p_color_sat_labels=( "cs10" "cs15" "cs20")

p_gamma_fac_arr=( 1.0 1.2 1.4 )
p_gamma_fac_labels=( "gf10" "gf12" "gf14" )

##p_bitpix_arr=( 8 16 )
p_bitpix_arr=( 8 )
p_bitpix_labels=( "b8" "b16" )

# Hardwired stiff options.
p_gamma=2.2
p_gamma_type=POWER-LAW
p_description=$p_pref
p_copyright=$(id -F)
p_verbose=FULL # Writes some annoying control characters...
p_write_xml=N

# Important scaling parameters:
# The use of quantile for the min is non-default, stiff normally uses GRAYEVEL
p_max_level="0.995,0.995,0.995"
p_max_type=="QUANTILE,QUANTILE,QUANTILE"
p_min_level="0.510,0.80,0.510"
p_min_type=="QUANTILE,QUANTILE,QUANTILE"
if [[ $p_pref == "m82" ]]; then
    p_max_level="5.0,4.5,4.0"
    p_max_type=="MANUAL,MANUAL,MANUAL"
    p_min_level="0.8,0.94,0.65"
    p_min_type=="MANUAL,MANUAL,MANUAL"
fi
##p_min_level="0.005,0.005,0.005"
##p_min_type=="MANUAL,MANUAL,MANUAL"
p_max_level="0.999,0.999,0.999"
p_max_type=="QUANTILE,QUANTILE,QUANTILE"
p_min_level="0.60,0.60,0.60"
p_min_type=="QUANTILE,QUANTILE,QUANTILE"



# Loop over 3-color combinations to process
for colsel in ${p_colsel_arr[@]}; do
    echo "-----------------------------------------------------------------"  | tee -a $p_log
    echo "Processing 3-color combination $colsel" | tee -a $p_log

    # Filters to process in red, green, blue order
    if [[ $colsel == 'sho' ]]; then
        p_filter_arr=("SII" "Ha" "OIII" )
    elif [[ $colsel == 'rgb' ]]; then
        p_filter_arr=("Red" "Green" "Blue")
    elif [[ $colsel == 'hgb' ]]; then
        p_filter_arr=("Ha" "Green" "Blue")
    else
        echo "Error, unexpected 3-color combination $colsel" | tee -a $p_log
        echo "  Allowed values are: ${p_allowed_colsel[@]}" | tee -a $p_log
        exit 4
    fi
    echo "  Filters in red, green, blue order: ${p_filter_arr[@]}" | tee -a $p_log
    
    echo "  Checking that input FITS files exist:" | tee -a $p_log
    p_files_arr=()
    p_nmissing=0
    for filter in ${p_filter_arr[@]}; do
        p_fname=$p_pref'_'$filter'_'$p_suff
        p_files_arr+=( $p_fname )
        if [ -e $p_fname ]; then
            echo "    Found $p_fname" | tee -a $p_log
        else
            echo "    Error, cannot find $p_fname" | tee -a $p_log
            (( p_nmissing++ ))
        fi
    done # end for filter
    if [ $p_nmissing -ne 0 ]; then
        echo "Error, missing $p_nmissing required files." | tee -a $p_log
        echo "  Current directory: $p_odir" | tee -a $p_log
        exit 8
    fi

    # Collapse filter name array into single string
    p_colstr=$(echo ${p_filter_arr[@]} | sed -e 's/ //g')

    # Loop over options for bits per pix, gamma_fac, and color_sat
    idx_bp=0
    for bitpix in ${p_bitpix_arr[@]}; do
        p_bitpix_str=${p_bitpix_labels[$idx_bp]}

        # Loop over luminance gamma factor
        idx_gf=0
        for p_gamma_fac in ${p_gamma_fac_arr[@]}; do
            p_gamma_fac_str=${p_gamma_fac_labels[$idx_gf]}
        
            # Loop over color saturation
            idx_cs=0
            for p_color_sat in ${p_color_sat_arr[@]}; do
                p_color_sat_str=${p_color_sat_labels[$idx_cs]}
        
                echo "    Bits-per-pixel=$bitpix, gamma-fac=$p_gamma_fac, color_sat=$p_color_sat"  | tee -a $p_log
                p_opt_str=$p_gamma_fac_str'_'$p_color_sat_str'_'$p_bitpix_str
        
                # Name of output TIFF file.
                p_colimg=$p_pref'_'$p_colstr'_'$p_short_suff'_'$p_opt_str'.tiff'
                p_fname_arr+=( $p_colimg )
                if [ -e $p_colimg ]; then
                    rm $p_colimg
                fi
                        
                echo "    Running stiff, generating $p_colimg at" $(date) | tee -a $p_log
                echo "    stiff command: " $p_stiff ${p_files_arr[@]} -OUTFILE_NAME $p_colimg \
                    -GAMMA_TYPE  $p_gamma_type -GAMMA $p_gamma \
                    -GAMMA_FAC $p_gamma_fac -COLOUR_SAT $p_color_sat \
                    -MAX_TYPE $p_max_type -MAX_LEVEL $p_max_level \
                    -MIN_TYPE $p_min_type -MIN_LEVEL $p_min_level \
                    -BITS_PER_CHANNEL $bitpix -VERBOSE_TYPE "$p_verbose" \
                    -DESCRIPTION $p_description -COPYRIGHT "$p_copyright" \
                    -WRITE_XML $p_write_xml >> $p_log
                p_tiff_start_time=$(dks_time)
                $p_stiff ${p_files_arr[@]} -OUTFILE_NAME $p_colimg \
                    -GAMMA_TYPE  $p_gamma_type -GAMMA $p_gamma \
                    -GAMMA_FAC $p_gamma_fac -COLOUR_SAT $p_color_sat \
                    -MAX_TYPE $p_max_type -MAX_LEVEL $p_max_level \
                    -MIN_TYPE $p_min_type -MIN_LEVEL $p_min_level \
                    -BITS_PER_CHANNEL $bitpix -VERBOSE_TYPE "$p_verbose" \
                    -DESCRIPTION $p_description -COPYRIGHT "$p_copyright" \
                    -WRITE_XML $p_write_xml >& $p_tmp
                p_status=$?
                cat $p_tmp >> $p_log
                if [ $p_status -ne 0 ]; then
                    echo "Error, $p_stiff failed with status $p_status" | tee -a $p_log
                fi
                p_fstat_arr+=( $p_status )
                
                # Check that output image exists.
                if [ -e $p_colimg ]; then
                    echo "      Generated $p_colimg" | tee -a $p_log
                else
                    echo "      Warning, cannot find $p_colimg" | tee -a $p_log
                fi
                
                # Timing info. Might be useful if we switch parameters or even the
                # resampling method used.
                p_tiff_end_time=$(dks_time)
                p_tiff_time=$(dks_time $p_tiff_start_time $p_tiff_end_time)
                p_time_arr+=( $p_tiff_time )
                echo "      Stiff took $p_tiff_time seconds." | tee -a $p_log
                echo ""  | tee -a $p_log

                # Increment color sat counter
                (( idx_cs++ ))
            done # end for color_sat

            # Increment gamma factor counter
            (( idx_gf++ ))
        done # end for gamma_fac

        # Increment bitpix counter.
        (( idx_bp++ ))
    done # for each bitpix

done # end for colsel

# Print summary of number of files successfully processed, skipped or failed.
echo "" | tee -a $p_log
echo "-----------------------------------------------------------------" | tee -a $p_log
echo "Run summary:" | tee -a $p_log
p_num=${#p_fname_arr[@]}
idx=0
printf "%-3s  %-60s  %8s  %6s\n" "Img" \
    "Resampled Output" "Time" "Status"| tee -a $p_log
while [ $idx -lt $p_num ]; do
    printf "%-3d  %-60s  %8.3f  %6d\n" $idx \
        ${p_fname_arr[$idx]} ${p_time_arr[$idx]} \
        ${p_fstat_arr[$idx]} | tee -a $p_log
    ((idx++))
done
echo "Master log file: $p_log" | tee -a $p_log

#-----------------------------------------------------------------------
# Clean up and exit.

# Remove temp file
if [ -e $p_tmp ]; then
    rm $p_tmp
fi

# All done
p_end=$(date)
p_end_time=$(dks_time)
p_run_time=$(dks_time $p_start_time $p_end_time)
echo "$0 finished at $p_end, run time $p_run_time seconds." | tee -a $p_log
exit 0

