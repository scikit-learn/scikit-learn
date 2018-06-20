#!/bin/sh
# Script to do a local install of loky
set +x
export LC_ALL=C
INSTALL_FOLDER=tmp/loky_install
rm -rf loky $INSTALL_FOLDER 2> /dev/null
if [ -z "$1" ]
then
    # Grab the latest stable release from PyPI
    LOKY=loky
else
    LOKY=$1
fi
pip install $LOKY --target $INSTALL_FOLDER
cp -r $INSTALL_FOLDER/loky .
rm -rf $INSTALL_FOLDER

# Needed to rewrite the doctests
# Note: BSD sed -i needs an argument unders OSX
# so first renaming to .bak and then deleting backup files
find loky -name "*.py" | xargs sed -i.bak "s/from loky/from joblib.externals.loky/"
find loky -name "*.bak" | xargs rm

sed -i "s/loky.backend.popen_loky/joblib.externals.loky.backend.popen_loky/" loky/backend/popen_loky_posix.py
