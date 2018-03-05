#!/bin/bash
#------------------------------------------------------------------------------
# Author   : Albert Akhriev (albert_akhriev@ie.ibm.com)
# Copyright: IBM Research Ireland, 2017-2018
#------------------------------------------------------------------------------

# Clear any previous build.
/bin/rm -fr build
sync
mkdir -p build
mkdir -p output
sync

# Optionally load a fresh Allscale API
# and Armadillo library for unit tests.
# bash ./scripts/download.sh

# Number of CPU utilized.
NCPU=2

# Build and run unit tests.
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../code
make -j ${NCPU}
ctest -j ${NCPU}
cd ../

