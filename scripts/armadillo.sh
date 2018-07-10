#!/bin/bash
#------------------------------------------------------------------------------
# Author   : Albert Akhriev (albert_akhriev@ie.ibm.com)
# Copyright: IBM Research Ireland, 2017-2018
#------------------------------------------------------------------------------

VER=8.600.0

FNAME="armadillo-${VER}.tar.xz"
DNAME="armadillo-${VER}"
ADDR=http://sourceforge.net/projects/arma/files/${FNAME}

############################################################
# Utilities.
############################################################
export MY_ARCHIVE="./archives"

# Error handling.
THIS=`basename $0`
Verify() {
    if [ ${1} -ne 0 ]; then
        msg="ERROR at: ${THIS}:${2}"
        echo; echo ${msg}; echo ${msg} >> error.log; echo; echo; exit 1
    fi
}

Warning() {
    if [ ${1} -ne 0 ]; then
        echo; echo "WARNING at ${THIS}: ${2}" >> warning.log; echo
    fi
}

# Function loads an archive file. Using: LoadArchive "file name" "url address".
LoadArchive() {
    mkdir -p ${MY_ARCHIVE}                     ; Verify $? $LINENO
    if [ ! -f "${MY_ARCHIVE}/${1}" ]; then
        curl -L -o "${MY_ARCHIVE}/${1}" "${2}" ; Verify $? $LINENO
    fi
}

# Function removes specified directory. TODO: check the argument
RemoveDir() {
    /bin/rm -fr ${1}
}
############################################################

# Set pristine environmental variables, otherwise the installer searches
# my local library folders.
export CPATH=
export C_INCLUDE_PATH=
export CPLUS_INCLUDE_PATH=
export LIBRARY_PATH=
export LD_LIBRARY_PATH=
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
# My local libraries live in "$SOFT"; use the system-wide libraries only.
export SOFT="."

# Create destination directory.
DEST="./api/armadillo"
mkdir -p "${DEST}"

# Load the library.
LoadArchive ${FNAME} ${ADDR}
if [[ ${DOWNLOAD_ONLY} == true ]] ; then exit 0; fi

# Unpack, build and delete distribution.
RemoveDir ${DNAME}                              ;  Verify $? $LINENO
echo "${MY_ARCHIVE}/${FNAME}"
tar xJf "${MY_ARCHIVE}/${FNAME}"                ;  Verify $? $LINENO
CWD=$(pwd)                                      ;  Verify $? $LINENO
cd ${DNAME}                                     ;  Verify $? $LINENO
./configure -DCMAKE_INSTALL_PREFIX="../${DEST}" ;  Verify $? $LINENO
make                                            ;  Verify $? $LINENO
make install                                    ;  Verify $? $LINENO
cd "${CWD}"                                     ;  Verify $? $LINENO
RemoveDir ${DNAME}                              ;  Verify $? $LINENO
RemoveDir ${MY_ARCHIVE}                         ;  Verify $? $LINENO

# -DOpenBLAS_LIBRARY="${SOFT}/superlu/"
# -DSuperLU_LIBRARY="${SOFT}/superlu/"


