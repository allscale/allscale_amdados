#!/bin/bash
#------------------------------------------------------------------------------
# Author   : Albert Akhriev (albert_akhriev@ie.ibm.com)
# Copyright: IBM Research Ireland, 2017-2018
#------------------------------------------------------------------------------

# Check we are in the root directory.
rootDir="allscale_amdados"
currDir=$(basename "$PWD")
if [ "${rootDir}" != "${currDir}" ]
then
    echo "ERROR: the script must be executed from the project root directory: ${rootDir}"
    exit 1
fi

# Remove (clear) and then make the destination folder.
/bin/rm -fr api
sync
mkdir -p api
sync

# Download and build Armadillo library (for some tests).
bash ./scripts/armadillo.sh

cd api

# Download the google-test suit.
curl -L -o gtest-1.8.0.tar.gz http://insieme-compiler.org/ext_libs/gtest-1.8.0.tar.gz

# Clone the AllScale API project.
/bin/rm -fr allscale_api
git clone https://github.com/allscale/allscale_api

cd ../

