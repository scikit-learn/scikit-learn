#!/bin/sh
# Script to do a local install of threadpoolctl
set +x
export LC_ALL=C
INSTALL_FOLDER=threadpoolctl_install
rm -rf _threadpoolctl.py $INSTALL_FOLDER 2> /dev/null
if [ -z "$1" ]
then
    # Grab the latest stable release from PyPI
    THREADPOOLCTL=threadpoolctl
else
    THREADPOOLCTL=$1
fi
pip install --no-cache $THREADPOOLCTL --target $INSTALL_FOLDER
cp $INSTALL_FOLDER/threadpoolctl.py _threadpoolctl.py
rm -rf $INSTALL_FOLDER