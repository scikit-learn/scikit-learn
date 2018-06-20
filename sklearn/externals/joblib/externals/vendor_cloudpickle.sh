#!/bin/sh
# Script to do a local install of cloudpickle
export LC_ALL=C
INSTALL_FOLDER=tmp/cloudpickle_install
rm -rf cloudpickle $INSTALL_FOLDER 2> /dev/null
if [ -z "$1" ]
then
    # Grab the latest stable release from PyPI
    CLOUDPICKLE=cloudpickle
else
    CLOUDPICKLE=$1
fi
pip install $CLOUDPICKLE --target $INSTALL_FOLDER
cp -r $INSTALL_FOLDER/cloudpickle .
rm -rf $INSTALL_FOLDER

# Needed to rewrite the doctests
# Note: BSD sed -i needs an argument unders OSX
# so first renaming to .bak and then deleting backup files
sed -i.bak "s/from cloudpickle.cloudpickle/from .cloudpickle/" cloudpickle/__init__.py
find cloudpickle -name "*.bak" | xargs rm

rm -r tmp
