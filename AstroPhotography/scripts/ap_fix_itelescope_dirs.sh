#!/bin//bash
#
#-----------------------------------------------------------------------
# ap_fix_itelescope_dirs.sh 
#
# Recursively 
# (a) Change the permissions of all the current directory's 
#     sub-directories so that the :
#     - user permisions are rwx, 
#     - guest and other can rx but not w.
# (b) Replace spaces in directory names with underscores
#
# Verified to work on directories with spaces in their names.
#
# @history 2020-07-30 dks : Coded to fix iTelescope data directories permissions.
#-----------------------------------------------------------------------

# This must be done recursively
for p_depth in 1 2 3 4 5; do
    echo "Traversing directory tree depth $p_depth"
    find . -maxdepth $p_depth -mindepth $p_depth -type d | while read p_dir; do 
        # Fix permissions
        echo "Changing permissions on $p_dir"
        chmod u+rwx,go+rx,go-x "$p_dir"
        
        # Remove any spaces in the file name
        re="[[:space:]]+"
        if [[ "$p_dir" =~ $re ]]; then
            p_out=$(echo "$p_dir" | sed -e 's/\ /_/g')
            echo "Renaming $p_dir to remove whitespace: $p_out"
            mv "$p_dir" "$p_out"
        fi
    done
done

exit 0
