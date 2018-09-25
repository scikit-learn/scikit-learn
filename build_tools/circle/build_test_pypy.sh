#!/usr/bin/env bash
set -x
set -e

apt-get -yq update
apt-get -yq install libatlas-dev libatlas-base-dev liblapack-dev gfortran ccache

pip install virtualenv

if command -v pypy3; then
    virtualenv -p $(command -v pypy3) pypy-env
elif command -v pypy; then
    virtualenv -p $(command -v pypy) pypy-env
fi

source pypy-env/bin/activate

python --version
which python

pip install --extra-index https://antocuni.github.io/pypy-wheels/ubuntu numpy Cython pytest
pip install "scipy>=1.1.0" sphinx numpydoc docutils

ccache -M 512M
export CCACHE_COMPRESS=1
export PATH=/usr/lib/ccache:$PATH
export LOKY_MAX_CPU_COUNT="2"

pip install -vv -e . 

python -m pytest sklearn/
python -m pytest doc/sphinxext/
python -m pytest $(find doc -name '*.rst' | sort)
