#!/bin/bash
#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

# Play AVI file of the true density field.
PLAY_TRUE_FIELD=false

# Play AVI file of both - the true density field and data assimilation solution.
PLAY_BOTH_FIELDS=false

# Plot the relative difference between the true density field and data assimilation solution.
PLOT_REL_DIFF=false

# Parsing the command-line arguments.
while getopts 'tbd' OPTION
do
    case ${OPTION} in
        t) PLAY_TRUE_FIELD=true;;
        b) PLAY_BOTH_FIELDS=true;;
        d) PLOT_REL_DIFF=true;;
        ?) echo "unknown option: ${OPTION}";;
    esac
done

if ${PLAY_TRUE_FIELD}; then
    echo "playing the true density field AVI"
    mplayer -nosound -x 512 -y 512 ./output/true_field.avi
fi

if ${PLAY_BOTH_FIELDS}; then
    echo "playing both density fields AVI"
    mplayer -nosound -x 1024 -y 512 ./output/both_fields.avi
fi

if ${PLOT_REL_DIFF}; then
    echo "plotting the relative difference between both fields"
    python3 ./scripts/python/plot_rel_diff.py
fi


