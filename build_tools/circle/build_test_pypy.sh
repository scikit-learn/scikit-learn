#!/usr/bin/env bash
set -x
set -e

# System build tools
apt-get -yq update
apt-get -yq install wget bzip2 build-essential ccache

# Install pypy and all the scikit-learn dependencies from conda-forge. In
# particular, we want to install pypy compatible binary packages for numpy and
# scipy as it would be to costly to build those from source.
conda install -y mamba
mamba create -n pypy -y \
    pypy numpy scipy cython \
    joblib threadpoolctl pillow pytest \
    sphinx numpydoc docutils

eval "$(conda shell.bash hook)"
conda activate pypy

# Check that we are running PyPy instead of CPython in this environment.
python --version
which python
python -c "import platform; assert platform.python_implementation() == 'PyPy'"

# Build and install scikit-learn in dev mode
ccache -M 512M
export CCACHE_COMPRESS=1
export PATH=/usr/lib/ccache:$PATH
export LOKY_MAX_CPU_COUNT="2"
export OMP_NUM_THREADS="1"
# Set parallelism to 3 to overlap IO bound tasks with CPU bound tasks on CI
# workers with 2 cores when building the compiled extensions of scikit-learn.
export SKLEARN_BUILD_PARALLEL=3
pip install --no-build-isolation -e .

python -m pytest sklearn
