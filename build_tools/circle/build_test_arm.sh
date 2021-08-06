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

MINICONDA_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-arm64.sh"

# Install Miniconda
wget $MINICONDA_URL -O miniconda.sh
MINICONDA_PATH=$HOME/miniconda
chmod +x miniconda.sh && ./miniconda.sh -b -p $MINICONDA_PATH
export PATH=$MINICONDA_PATH/bin:$PATH
conda update --yes conda

# Create environment and install dependencies
conda create -n testenv --yes python=3.7
source activate testenv

# Use the latest by default
conda install --verbose -c conda-forge -y ccache numpy scipy cython pip
pip install joblib threadpoolctl
pip install $(get_dep pytest $PYTEST_VERSION) pytest-xdist

setup_ccache
python -m pip install $(get_dep cython $CYTHON_VERSION) \
                      $(get_dep joblib $JOBLIB_VERSION)
python -m pip install $(get_dep threadpoolctl $THREADPOOLCTL_VERSION) \
                      $(get_dep pytest $PYTEST_VERSION) \
                      $(get_dep pytest-xdist $PYTEST_XDIST_VERSION)

if [[ "$COVERAGE" == "true" ]]; then
    python -m pip install codecov pytest-cov
fi

if [[ "$PYTEST_XDIST_VERSION" != "none" ]]; then
    python -m pip install pytest-xdist
fi

if [[ "$TEST_DOCSTRINGS" == "true" ]]; then
    # numpydoc requires sphinx
    python -m pip install sphinx
    python -m pip install numpydoc
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
