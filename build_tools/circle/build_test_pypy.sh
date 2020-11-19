#!/usr/bin/env bash
set -x
set -e

# System build tools
sudo apt-get -yq update
sudo apt-get -yq install wget bzip2 build-essential ccache

# Init micromamba, a fast and lightweight alternative to conda for CI.
rm -rf ~/micromamba
mkdir -p $HOME/bin
pushd $HOME
wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
popd
$HOME/bin/micromamba shell init -s bash -p $HOME/micromamba
source ~/.bashrc
micromamba activate

# Install pypy and all the scikit-learn dependencies from conda-forge. In
# particular, we want to install pypy compatible binary packages for numpy and
# scipy as it would be to costly to build those from source.
micromamba install -y -c conda-forge \
    pypy numpy scipy cython \
    joblib threadpoolctl pillow pytest \
    sphinx numpydoc docutils

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
python -m pytest doc/sphinxext/
python -m pytest $(find doc -name '*.rst' | sort)
