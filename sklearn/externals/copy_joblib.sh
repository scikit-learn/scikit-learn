#!/bin/sh
# Script to do a local install of joblib
export LC_ALL=C
INSTALL_FOLDER=tmp/joblib_install
rm -rf joblib $INSTALL_FOLDER
pip install joblib --target $INSTALL_FOLDER
cp -r $INSTALL_FOLDER/joblib .
rm -rf $INSTALL_FOLDER

# Needed to rewrite the doctests
# Note: BSD sed -i needs an argument unders OSX
# so first renaming to .bak and then deleting backup files
find joblib -name "*.py" | xargs sed -i.bak "s/from joblib/from sklearn.externals.joblib/"
find joblib -name "*.bak" | xargs rm

# Remove the tests folders to speed-up test time for scikit-learn.
# joblib is already tested on its own CI infrastructure upstream.
rm -r joblib/test

# Remove joblib/testing.py which is only used in tests and has a
# pytest dependency (needed until we drop nose)
rm joblib/testing.py
