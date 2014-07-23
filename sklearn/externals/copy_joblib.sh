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

# Make it possible to use the system joblib based on an environment variable
mv joblib/__init__.py joblib/embedded_joblib_init.py
NEW_INIT="
import os
import warnings
if int(os.environ.get('SKLEARN_SYSTEM_JOBLIB', '0')):
    # Try to import joblib from the system for development purpose
    from joblib import *
    from joblib import __version__
    from .embedded_joblib_init import __version__ as orig_version
    warnings.warn('Using system joblib %s instead of embedded joblib %s'
                  % (__version__, orig_version))
else:
    from .embedded_joblib_init import *
"
echo "$NEW_INIT" > joblib/__init__.py

chmod -x joblib/*.py
