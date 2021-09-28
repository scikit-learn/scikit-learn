#!/bin/bash

set -e
set -x

UNAMESTR=`uname`

setup_ccache() {
    echo "Setting up ccache"
    mkdir /tmp/ccache/
    which ccache
    for name in gcc g++ cc c++ x86_64-linux-gnu-gcc x86_64-linux-gnu-c++; do
      ln -s $(which ccache) "/tmp/ccache/${name}"
    done
    export PATH="/tmp/ccache/:${PATH}"
    ccache -M 256M
}

# imports get_dep
source build_tools/shared.sh

sudo add-apt-repository --remove ppa:ubuntu-toolchain-r/test
sudo apt-get update

# Setup conda environment
MINICONDA_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-aarch64.sh"

# Install Mambaforge
wget $MINICONDA_URL -O mambaforge.sh
MINICONDA_PATH=$HOME/miniconda
chmod +x mambaforge.sh && ./mambaforge.sh -b -p $MINICONDA_PATH
export PATH=$MINICONDA_PATH/bin:$PATH
mamba update --yes conda

# Create environment and install dependencies
mamba create -n testenv --yes python=3.7
source activate testenv

# Use the latest by default
mamba install --verbose -y  ccache \
                            pip \
                            $(get_dep numpy $NUMPY_VERSION) \
                            $(get_dep scipy $SCIPY_VERSION) \
                            $(get_dep cython $CYTHON_VERSION) \
                            $(get_dep joblib $JOBLIB_VERSION) \
                            $(get_dep threadpoolctl $THREADPOOLCTL_VERSION) \
                            $(get_dep pytest $PYTEST_VERSION) \
                            $(get_dep pytest-xdist $PYTEST_XDIST_VERSION)
setup_ccache

if [[ "$COVERAGE" == "true" ]]; then
    mamba install --verbose -y codecov pytest-cov
fi

if [[ "$TEST_DOCSTRINGS" == "true" ]]; then
    # numpydoc requires sphinx
    mamba install --verbose -y sphinx
    mamba install --verbose -y numpydoc
fi

python --version
conda list

# Set parallelism to 3 to overlap IO bound tasks with CPU bound tasks on CI
# workers with 2 cores when building the compiled extensions of scikit-learn.
export SKLEARN_BUILD_PARALLEL=3

pip install --verbose --editable . --no-build-isolation
ccache -s
python -c "import sklearn; sklearn.show_versions()"
python -m threadpoolctl --import sklearn
python -m pytest sklearn
