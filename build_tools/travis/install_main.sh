#!/bin/bash

# Travis clone "scikit-learn/scikit-learn" repository into
# a local repository. We use a cached directory with three
# scikit-learn repositories (one for each matrix entry for
# non continuous deployment jobs) from which we pull local
# Travis repository. This allows us to keep build artifact
# for GCC + Cython, and gain time.

set -e

echo "CPU Arch: $TRAVIS_CPU_ARCH."

# Import "get_dep"
source build_tools/shared.sh

echo "List files from cached directories."
echo "pip:"
ls $HOME/.cache/pip

export CC=/usr/lib/ccache/gcc
export CXX=/usr/lib/ccache/g++

# Useful for debugging how ccache is used
# export CCACHE_LOGFILE=/tmp/ccache.log

# 60MB are (more or less) used by .ccache, when
# compiling from scratch at the time of writing
ccache --max-size 100M --show-stats

# Deactivate the default virtual environment
# to setup a conda-based environment instead
deactivate

MINICONDA_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh"

# Install Miniconda
wget $MINICONDA_URL -O miniconda.sh
MINICONDA_PATH=$HOME/miniconda
chmod +x miniconda.sh && ./miniconda.sh -b -p $MINICONDA_PATH
export PATH=$MINICONDA_PATH/bin:$PATH
conda update --yes conda

# Create environment and install dependencies
conda create -n testenv --yes python=3.7

source activate testenv
conda install -y scipy numpy pandas cython
pip install joblib threadpoolctl

pip install $(get_dep pytest $PYTEST_VERSION) pytest-xdist

# Build scikit-learn in this script to collapse the
# verbose build output in the Travis output when it
# succeeds
python --version
python -c "import numpy; print(f'numpy {numpy.__version__}')"
python -c "import scipy; print(f'scipy {scipy.__version__}')"

pip install -e .
python setup.py develop

ccache --show-stats

# Useful for debugging how ccache is used
# cat $CCACHE_LOGFILE
