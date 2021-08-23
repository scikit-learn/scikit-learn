#!/bin/bash

set -e
set -x

sudo add-apt-repository --remove ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install python3-virtualenv ccache
python3 -m virtualenv --system-site-packages --python=python3 testenv
source testenv/bin/activate
pip install --upgrade pip
for pyversion in 37 38 39; do
    pywheel=`echo "cp"$pyversion"-manylinux_aarch64"`
    export CIBW_BUILD=$pywheel
    python -m pip install cibuildwheel
        python -m cibuildwheel --output-dir wheelhouse
done

