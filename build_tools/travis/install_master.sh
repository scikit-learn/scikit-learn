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

if [[ $TRAVIS_CPU_ARCH == arm64 ]]; then
    # Different Miniconda URL for ARM64 architectures
    MINICONDA_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh"
else
    MINICONDA_URL="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
fi

# Install Miniconda
wget $MINICONDA_URL -O miniconda.sh
MINICONDA_PATH=$HOME/miniconda
chmod +x miniconda.sh && ./miniconda.sh -b -p $MINICONDA_PATH
export PATH=$MINICONDA_PATH/bin:$PATH
conda update --yes conda

# Create environment and install dependencies
conda create -n testenv --yes python=3.7

source activate testenv

if [[ $TRAVIS_CPU_ARCH == amd64 ]]; then
    echo "Upgrading pip and setuptools."
    pip install --upgrade pip setuptools
    echo "Installing numpy, scipy and pandas master wheels."
    dev_anaconda_url=https://pypi.anaconda.org/scipy-wheels-nightly/simple
    pip install --pre --upgrade --timeout=60 --extra-index $dev_anaconda_url numpy scipy pandas
    echo "Installing cython pre-release wheels."
    pip install --pre cython
    echo "Installing joblib master."
    pip install https://github.com/joblib/joblib/archive/master.zip
    echo "Installing pillow master."
    pip install https://github.com/python-pillow/Pillow/archive/master.zip
else
    conda install -y scipy numpy pandas cython
    pip install joblib threadpoolctl
fi

pip install $(get_dep pytest $PYTEST_VERSION) pytest-xdist

# Build scikit-learn in this script to collapse the 
# verbose build output in the Travis output when it
# succeeds
python --version
python -c "import numpy; print(f'numpy {numpy.__version__}')"
python -c "import scipy; print(f'scipy {scipy.__version__}')"

if [[ $BUILD_WITH_ICC == true ]]; then
    wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
    sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
    rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
    sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
    sudo apt-get update
    sudo apt-get install intel-oneapi-icc
    source /opt/intel/oneapi/setvars.sh

    # The "build_clib" command is implicitly used to build "libsvm-skl".
    # To compile with a different compiler, we also need to specify the
    # compiler for this command
    python setup.py build_ext --compiler=intelem -i build_clib --compiler=intelem
else
    pip install -e .
fi

python setup.py develop

ccache --show-stats

# Useful for debugging how ccache is used
# cat $CCACHE_LOGFILE
