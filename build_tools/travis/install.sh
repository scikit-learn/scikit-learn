#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See https://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD

# Travis clone scikit-learn/scikit-learn repository in to a local repository.
# We use a cached directory with three scikit-learn repositories (one for each
# matrix entry) from which we pull from local Travis repository. This allows
# us to keep build artefact for gcc + cython, and gain time

set -e

# Fail fast
build_tools/travis/travis_fastfail.sh

echo "List files from cached directories"
echo "pip:"
ls $HOME/.cache/pip

export CC=/usr/lib/ccache/gcc
export CXX=/usr/lib/ccache/g++
# Useful for debugging how ccache is used
# export CCACHE_LOGFILE=/tmp/ccache.log
# ~60M is used by .ccache when compiling from scratch at the time of writing
ccache --max-size 100M --show-stats

# Deactivate the travis-provided virtual environment and setup a
# conda-based environment instead
# If Travvis has language=generic, deactivate does not exist. `|| :` will pass.
deactivate || :

# Install miniconda
fname=Miniconda3-latest-Linux-x86_64.sh
wget https://repo.continuum.io/miniconda/$fname -O miniconda.sh
MINICONDA_PATH=$HOME/miniconda
chmod +x miniconda.sh && ./miniconda.sh -b -p $MINICONDA_PATH
export PATH=$MINICONDA_PATH/bin:$PATH
conda update --yes conda

# Create environment and install dependencies
conda create -n testenv --yes python=3.7
source activate testenv

pip install --upgrade pip setuptools
echo "Installing numpy and scipy master wheels"
dev_url=https://7933911d6844c6c53a7d-47bd50c35cd79bd838daf386af554a83.ssl.cf2.rackcdn.com
pip install --pre --upgrade --timeout=60 -f $dev_url numpy scipy pandas cython
echo "Installing joblib master"
pip install https://github.com/joblib/joblib/archive/master.zip
echo "Installing pillow master"
pip install https://github.com/python-pillow/Pillow/archive/master.zip
pip install pytest==4.6.4 pytest-cov

# Build scikit-learn in the install.sh script to collapse the verbose
# build output in the travis output when it succeeds.
python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"

python setup.py develop

ccache --show-stats
# Useful for debugging how ccache is used
# cat $CCACHE_LOGFILE

# fast fail
build_tools/travis/travis_fastfail.sh
