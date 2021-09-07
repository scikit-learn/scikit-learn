#!/bin/bash

set -e
set -x

sudo apt-get update
sudo apt-get install -y wget

# Install Python 3.8 with conda from miniforge because this is
# the simplest way to get a working anaconda-client for the wheel
# upload step to the staging area on the anaconda.org service.
export CONDA=$HOME/mambaforge
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh
bash Mambaforge-$(uname)-$(uname -m).sh -b -f -p $CONDA
echo 'export CONDA=$CONDA' >> $BASH_ENV  # to be reused in a later step
export PATH=$CONDA/bin:$PATH
mamba update --yes conda
pip install --upgrade pip
for pyversion in 37 38 39; do
    pywheel=`echo "cp"$pyversion"-manylinux_aarch64"`
    export CIBW_BUILD=$pywheel
    python -m pip install cibuildwheel
    python -m cibuildwheel --output-dir wheelhouse
done
