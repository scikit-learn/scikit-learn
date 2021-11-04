#!/bin/bash

set -e
set -x

UNAMESTR=`uname`
N_CORES=`nproc --all`


setup_ccache() {
    echo "Setting up ccache"
    mkdir /tmp/ccache/
    which ccache
    for name in gcc g++ cc c++ x86_64-linux-gnu-gcc x86_64-linux-gnu-c++; do
      ln -s $(which ccache) "/tmp/ccache/${name}"
    done
    export PATH="/tmp/ccache:${PATH}"
    # Unset ccache limits
    ccache -F 0
    ccache -M 0
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
mamba init --all --verbose
mamba update --yes conda

# Create environment and install dependencies
mamba create -n testenv --yes $(get_dep python $PYTHON_VERSION)
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
                            $(get_dep pytest-xdist $PYTEST_XDIST_VERSION) \
                            "setuptools<58.5" # TODO: Remove this line once setuptools#2849 is resolved.

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

# Set parallelism to $N_CORES + 1 to overlap IO bound tasks with CPU bound tasks on CI
# workers with $N_CORES cores when building the compiled extensions of scikit-learn.
export SKLEARN_BUILD_PARALLEL=$(($N_CORES + 1))

# Disable the build isolation and build in the tree so that the same folder can be
# cached between CI runs.
# TODO: remove the '--use-feature' flag when made obsolete in pip 21.3.
pip install --verbose --no-build-isolation --use-feature=in-tree-build .

# Report cache usage
ccache -s --verbose

mamba list

# Changing directory not to have module resolution use scikit-learn source
# directory but to the installed package.
cd /tmp
python -c "import sklearn; sklearn.show_versions()"
python -m threadpoolctl --import sklearn
# Test using as many workers as available cores
pytest --pyargs -n $N_CORES sklearn
