#!/bin/bash
#--------------------------------------------------------------------------------------------------
# Author   : Albert Akhriev (albert_akhriev@ie.ibm.com)
# Copyright: IBM Research Ireland, 2014 - 2017
#--------------------------------------------------------------------------------------------------

# Check we are in the root directory.
rootDir="allscale_amdados"
currDir=$(basename "$PWD")
if [ "${rootDir}" != "${currDir}" ]
then
    echo "ERROR: the script must be executed from the project root directory: ${rootDir}"
    exit 1
fi

# Make the destination folder, if needed.
mkdir -p api

# Download and build Armadillo library (for some tests).
bash ./scripts/armadillo.sh

cd api

# Download the google-test suit.
curl -L -o gtest-1.8.0.tar.gz http://insieme-compiler.org/ext_libs/gtest-1.8.0.tar.gz

# Clone the AllScale API project.
git clone https://github.com/allscale/allscale_api

cd ../

