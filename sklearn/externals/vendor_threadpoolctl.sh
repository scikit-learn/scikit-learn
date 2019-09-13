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

# Needed to rewrite the doctests
# Note: BSD sed -i needs an argument unders OSX
# so first renaming to .bak and then deleting backup files
#find loky -name "*.py" | xargs sed -i.bak "s/from loky/from joblib.externals.loky/"
#find loky -name "*.bak" | xargs rm

#for f in $(git grep -l "cloudpickle" loky); do
#     echo $f;
#     sed -i 's/import cloudpickle/from joblib.externals import cloudpickle/' $f
#     sed -i 's/from cloudpickle import/from joblib.externals.cloudpickle import/' $f
# done

# sed -i "s/loky.backend.popen_loky/joblib.externals.loky.backend.popen_loky/" loky/backend/popen_loky_posix.py
