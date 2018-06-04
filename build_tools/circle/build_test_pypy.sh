#!/usr/bin/env bash
set -x
set -e

apt-get -yq update
apt-get -yq install libatlas-dev libatlas-base-dev liblapack-dev gfortran

pip install virtualenv

if command -v pypy3; then
    virtualenv -p $(command -v pypy3) pypy-env
elif command -v pypy; then
    virtualenv -p $(command -v pypy) pypy-env
fi

source pypy-env/bin/activate

python --version
which python

pip install --extra-index https://antocuni.github.io/pypy-wheels/ubuntu numpy Cython Tempita pytest
pip install "scipy>=1.1.0"
pip install -e .

pytest sklearn/
