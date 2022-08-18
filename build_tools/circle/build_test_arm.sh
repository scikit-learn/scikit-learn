#!/bin/bash

set -e
set -x

UNAMESTR=`uname`
N_CORES=`nproc --all`

# defines the get_dep and show_installed_libraries functions
source build_tools/shared.sh

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

MINICONDA_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-aarch64.sh"

# Install Mambaforge
wget $MINICONDA_URL -O mambaforge.sh
MINICONDA_PATH=$HOME/miniconda
chmod +x mambaforge.sh && ./mambaforge.sh -b -p $MINICONDA_PATH
export PATH=$MINICONDA_PATH/bin:$PATH
mamba init --all --verbose
mamba update --yes mamba
mamba update --yes conda
mamba install "$(get_dep conda-lock min)" -y
conda-lock install --name $CONDA_ENV_NAME $LOCK_FILE
source activate $CONDA_ENV_NAME

setup_ccache

python --version

# Set parallelism to $N_CORES + 1 to overlap IO bound tasks with CPU bound tasks on CI
# workers with $N_CORES cores when building the compiled extensions of scikit-learn.
export SKLEARN_BUILD_PARALLEL=$(($N_CORES + 1))

# Disable the build isolation and build in the tree so that the same folder can be
# cached between CI runs.
pip install --verbose --no-build-isolation .

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
