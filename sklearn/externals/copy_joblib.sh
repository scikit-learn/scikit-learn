#!/bin/sh
# Script to do a local install of joblib
rm -rf tmp joblib
mkdir -p tmp/lib/python2.7/site-packages
ln -s tmp/lib/python2.7 tmp/lib/python2.6
mkdir -p tmp/bin
export PYTHONPATH=$(pwd)/tmp/lib/python2.7/site-packages:$(pwd)/tmp/lib/python2.6/site-packages
easy_install -Zeab tmp joblib
old_pwd=$(pwd)
#cd /home/varoquau/dev/joblib/
cd tmp/joblib/
python setup.py install --prefix $old_pwd/tmp
cd $old_pwd
cp -r tmp/lib/python2.7/site-packages/joblib-*.egg/joblib .
rm -rf tmp
# Needed to rewrite the doctests
find joblib -name "*.py" | xargs sed -i.bak "s/from joblib/from sklearn.externals.joblib/"
find joblib -name "*.bak" | xargs rm

# Remove the tests folders to speed-up test time for scikit-learn.
# joblib is already tested on its own CI infrastructure upstream.
rm -r joblib/test

chmod -x joblib/*.py
