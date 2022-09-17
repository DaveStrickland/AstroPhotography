#!/bin/bash
#
#-----------------------------------------------------------------------
# bash ap_rename_files_with_space.sh
#
# Attempts to find all files (not directories) within the current
# directory tree that have spaces in their name, replacing " " with "_".
#
# @history 2020-08-03 dks : Clean up old script.
#-----------------------------------------------------------------------
find . -type f -name "* *" | while read filename; do
    p_out=$(echo "$filename" | sed -e 's/\ /_/g')
    #p_raw=$(echo "$filename" | sed -e 's/\ /\\\ /g')
    echo "$p_out"
    #echo "$p_raw"
    mv "$filename" "$p_out"
done

# All done.
exit 0

