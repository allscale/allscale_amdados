#!/bin/bash
#------------------------------------------------------------------------------
# Author   : Albert Akhriev (albert_akhriev@ie.ibm.com)
# Copyright: IBM Research Ireland, 2017-2018
#------------------------------------------------------------------------------

# Clear any previous build.
/bin/rm -fr build
sync

# Create the "build" and "output" directories.
# Note, the output directory (parameter "output_dir") in your configuration
# file might be different.
OUTDIR=$(grep 'output_dir' amdados.conf | awk '{print $2}')
mkdir -p build
mkdir -p ${OUTDIR}
sync

# Optionally load a fresh Allscale API and Armadillo library for unit tests.
# bash ./scripts/download.sh

# Number of CPU utilized.
NCPU=2

# Build and run unit tests.
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../code
make -j ${NCPU}
ctest -j ${NCPU}
cd ../

