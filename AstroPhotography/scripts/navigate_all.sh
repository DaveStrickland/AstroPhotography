#!/usr/bin/env bash
#-----------------------------------------------------------------------
# 
# navigate_all.sh <dosrc|nosrc> <donav|nonav> [noclean|clean]
# 
# Perform basic processing on a set of iTelescope files:
# 1. Run ap_find_stars on calibrated images to generate a list of 
#    detected star positions, approximate magnitudes, image quality 
#    metrics, and image-format ds9 region files.
# 2. Run ap_astrometry.py using the star positions to generate an
#    astrometric solution for the image, writing a new navigated
#    image file with a WCS solution AND generate RA, Dec for the
#    detected stars.
#
# - If 'dosrc' is supplied on the command line source searching will
#   be performed.
# - if 'donav' is supplied on the command line astrometry will be
#   performed.
# - If 'clean' is specified existing files will be removed and regenerated.
#   By default 'noclean' is set so existing files are not rebuilt.
#
# Limitations:
# - Only processes files matching the hard-wired pattern.
#
# History:
# 2020-07-?? dks : Initial version using iTelescope calibrated files.
# 2020-12-09 dks : Modified for use with ApCalibrate calibrated files.
#                  Add quality summary at end.
# 2021-02-21 dks : Make sleep between calls a variable.
# 2021-05-26 dks : Write log file locations to screen.
# 2021-08-24 dks : Improve logging.
#
#-----------------------------------------------------------------------
# Initialization

p_usage="$0 <dosrc|nosrc> <donav|nonav> [noclean|clean]"


# Log files for entire runs
p_log='tmp_nav_all.log'
if [ -e $p_log ]; then
    rm $p_log
fi

if [ $# -lt 2 ]; then
    echo "Error: Expecting at least 2 command line arguments, got $#" | tee -a $p_log
    echo "  usage: $p_usage" | tee -a $p_log
    exit 1
fi

# Whether to do source search (ap_star_find), navigation (ap_astrometry)
# an a quality summary (ap_quality_summary).
# These are defaults for src and nav, they're updated from the command line.
p_do_src=1
p_do_nav=0
p_do_qual=1

# Do source searching?
if [[ $1 == "dosrc" ]]; then
    p_do_src=1
    echo "Source searching will be performed." | tee -a $p_log
else
    p_do_src=0
    echo "Source searching disabled." | tee -a $p_log
fi

# Do astrometry?
if [[ $2 == "donav" ]]; then
    p_do_nav=1
    echo "Astrometry enabled. Assumes source searching completed." | tee -a $p_log
else
    p_do_nav=0
    echo "Astrometry disabled." | tee -a $p_log
fi

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

echo "$0 started at" $(date) | tee -a $p_log
echo "Running from" $(pwd) | tee -a $p_log
echo "Activating virtual environment." | tee -a $p_log
source ~/venv/astro38/bin/activate

# Log files for python runs
p_srclog='tmp_find.log'
p_navlog='tmp_astr.log'

# Scripts that do the work.
p_find_stars="$HOME/git/AstroPhotography/src/AstroPhotography/scripts/ap_find_stars.py"
p_astrometry="$HOME/git/AstroPhotography/src/AstroPhotography/scripts/ap_astrometry.py"
p_quality="$HOME/git/AstroPhotography/src/AstroPhotography/scripts/ap_quality_summary.py"
echo "Checking that scripts exist..."
p_err=0
for p_script in $p_find_stars $p_astrometry $p_quality; do
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
p_maxstars=200

# Amount of time to sleep between calls to Astrometry.net to avoid
# spamming the service (once it saturates all processing fails, so
# slow processing is much better than no processing.)
p_sleep=30

p_src_dir=SourceLists
p_nav_dir=NavigatedImages
p_meta_dir=MetaData
p_qualcsv=nav_all_summary.csv

for p_dir in $p_src_dir $p_nav_dir $p_meta_dir; do
    if [ -d $p_dir ]; then
        echo "Reusing existing $p_dir directory." | tee -a $p_log
    else
        echo "Creating $p_dir" | tee -a $p_log
        mkdir -p $p_dir
    fi
done

#-----------------------------------------------------------------------
# Main work area

# Process each calibrated file
for p_cal_file in $(find . -name "cal-*fits"); do
    # Set up names for expected inputs and outputs ---------------------

    # Basename of input calibrated file with suffix removed
    p_base_name=$(basename ${p_cal_file%.*})
    
    # Name of file for sourcelist
    p_src_file=$(basename $p_cal_file | sed -e 's/cal-/srclist_/')
    p_src_path=$p_src_dir/$p_src_file
    
    # Name and path of file for bitmapped image with sources marked.
    p_implot_file=$(echo $p_base_name | sed -e 's/cal-/implot_/')
    p_implot_path=$p_src_dir/$p_implot_file'.png'
    
    # Name and path of file for bitmapped star FWHM fitting
    p_fwhmplot_file=$(echo $p_base_name | sed -e 's/cal-/fwhmplot_/')
    p_fwhmplot_path=$p_src_dir/$p_fwhmplot_file'.png'
    
    # Name and path of file for ds9 region files
    p_ds9_file=$(echo $p_base_name | sed -e 's/cal-/ds9_/')
    p_ds9_path=$p_src_dir/$p_ds9_file'.reg'
    
    # Name and path of file for quality reports
    p_qual_file=$(echo $p_base_name | sed -e 's/cal-/qual_/')
    p_qual_path=$p_meta_dir/$p_qual_file'.yaml'
    
    # Name and path of file for source searching log
    p_srclog_file=$(echo $p_base_name | sed -e 's/cal-/srclog_/')
    p_srclog_path=$p_meta_dir/$p_srclog_file'.log'

    # Name and path of file for source searching log
    p_navlog_file=$(echo $p_base_name | sed -e 's/cal-/navlog_/')
    p_navlog_path=$p_meta_dir/$p_navlog_file'.log'
    
    # Name of file for navigated image
    p_nav_file=$(basename $p_cal_file | sed -e 's/cal-/nav_/')
    p_nav_path=$p_nav_dir/$p_nav_file

    echo "Processing $p_cal_file at" $(date)

    p_allfiles=( $p_src_path $p_implot_path $p_ds9_path $p_qual_path $p_srclog_path $p_nav_path $p_navlog_path )
    
    # Clean? This forces a full re-run of all processing.
    if [ $p_clean -ne 0 ]; then
        echo "  Cleaning all files..."  | tee -a $p_log
        for file in ${p_allfiles[@]}; do
            if [ -e $file ]; then
                rm -v $file  | tee -a $p_log
            fi
        done
    fi

    # Do star finding --------------------------------------------------
    if [ $p_do_src -gt 0 ]; then
        # Check if a source list exists
        if [ -e $p_src_path ]; then
            echo "  Skipping star detection as sourcelist $p_src_path exists." | tee -a $p_log
        else
            echo "  Performing star detection on $p_cal_file" | tee -a $p_log
            python3 $p_find_stars -l $p_loglevel -m $p_maxstars \
                $p_cal_file $p_src_path --retain_saturated \
                --plotfile=$p_implot_path --fwhm_plot=$p_fwhmplot_path \
                --quality_report=$p_qual_path --ds9=$p_ds9_path >& $p_srclog
            p_status=$?
            if [ $p_status -ne 0 ]; then
                echo "Error: ap_find_stars.py exited with error code $p_status" | tee -a $p_log
                echo "  Log file is $p_srclog" | tee -a $p_log
                exit 8
            elif [ ! -e $p_src_path ]; then
                echo "Error, could not finded expected output source list." | tee -a $p_log
                echo "  Failed to find $p_src_path" | tee -a $p_log
                echo "  Current dir:" $(pwd) | tee -a $p_log
                exit 16
            else
                echo "    Star detection completed successfully at" $(date) | tee -a $p_log
            fi
            mv $p_srclog $p_srclog_path
            echo "      Log file written to $p_srclog_path" | tee -a $p_log
            
            # end if whether srclist already exists
        fi
        # end if p_do_src
    fi
    
    # Do astrometry ----------------------------------------------------
    if [ $p_do_nav -gt 0 ]; then
        # Check if a navigated image exists
        if [ -e $p_nav_path ]; then
            echo "  Skipping astrometry as navigated image $p_nav_path exists." | tee -a $p_log
        else
            # Check sourcelist exists
            if [ ! -e $p_src_path ]; then
                echo "Error, cannot find source list $p_src_path" | tee -a $p_log
                echo "  Current path:" $(pwd) | tee -a $p_log
                exit 16
            fi
            
            echo "  Performing astrometry on $p_cal_file" | tee -a $p_log
            python3 $p_astrometry -l $p_loglevel \
                $p_cal_file $p_src_path $p_nav_path >& $p_navlog
            p_status=$?
            if [ $p_status -ne 0 ]; then
                echo "Error: ap_astrometry.py exited with error code $p_status" | tee -a $p_log
                echo "  Log file is $p_navlog" | tee -a $p_log
                exit 32
            else
                echo "    Astrometry completed successfully at" $(date) | tee -a $p_log
                echo "    Sleeping for $p_sleep seconds..." | tee -a $p_log
                sleep $p_sleep
            fi
            mv $p_navlog $
            echo "      Log file written to $p_navlog_path" | tee -a $p_log
                
            # end if whether navigated image exists
        fi
        # end if p_do_nav
    fi
    # for each calibrated file.
done

if [ $p_do_qual -eq 1 ]; then
    echo "About to generate a CSV summary from the find stars quality files." | tee -a $p_log
    if [ -e $p_qualcsv ]; then
        rm $p_qualcsv
    fi
    python3 $p_quality $p_meta_dir $p_qualcsv --loglevel=DEBUG
    if [ $? -eq 0 ]; then
        echo "  Succesfully generated $p_qualcsv" | tee -a $p_log
    else
        echo "Error, failed to generate $p_qualcsv" | tee -a $p_log
    fi
fi

#-----------------------------------------------------------------------
# Clean up and shut down

# Deactivate virtual env.
echo "Deactivating virtual environment." | tee -a $p_log
deactivate

# All done
echo "$0 finished at" $(date) | tee -a $p_log
echo "Run log file written to" $p_log
exit 0
