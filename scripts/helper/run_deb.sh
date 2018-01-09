#!/bin/bash

echo "-----------------------------------------"
echo "Running Amdados application in Debug mode"
echo "-----------------------------------------"

./mybuild -d

DIR=$(cat amdados.conf | grep 'output_dir' | awk '{print $2}')
mkdir -p ${DIR}
/bin/rm -f ${DIR}/*.png
/bin/rm -f ${DIR}/*.pgm
/bin/rm -f ${DIR}/*.jpg
/bin/rm -f ${DIR}/*.avi
/bin/rm -f ${DIR}/*.log
/bin/rm -f ${DIR}/*.out
/bin/rm -f ${DIR}/*profile*.txt
sync

lldb ./build/app/amdados
