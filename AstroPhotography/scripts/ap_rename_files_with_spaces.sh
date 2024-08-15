#!/bin/bash
#
#-----------------------------------------------------------------------
# bash ap_rename_files_with_space.sh
#
# Attempts to find all files (not directories) within the current
# directory tree that have spaces in their name, replacing " " with "_".
#
# @history 2020-08-03 dks : Clean up old script.
# @history 2024-08-13 dks : Give chance to quit if running in home dir
#-----------------------------------------------------------------------

# Running this script in your home directory can mess up hidden application
# settings, e.g. in .local, and so on. Give chance to quit out....

if [[ "$HOME" == "$(pwd)" ]]; then
    echo "WARNING: You are running this script in your HOME directory $HOME"
    echo "  This can cause problems by renaming hidden application settings."
    echo "  You have 10 seconds to Control-C to cancel out of this script."
    for i in {1..9}; do
        echo -n "."
        sleep 1
    done
    echo "."
    sleep 1
    echo "Continuing with rename..."
fi

find . -type f -name "* *" | while read filename; do
    p_out=$(echo "$filename" | sed -e 's/\ /_/g')
    #p_raw=$(echo "$filename" | sed -e 's/\ /\\\ /g')
    echo "$p_out"
    #echo "$p_raw"
    mv "$filename" "$p_out"
done

# All done.
exit 0

