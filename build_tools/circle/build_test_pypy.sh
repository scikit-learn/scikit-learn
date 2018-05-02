#!/usr/bin/env bash
set -x
set -e


apt-get -yq update
apt-get -yq install libatlas-dev libatlas-base-dev liblapack-dev gfortran

pip install virtualenv

# python executable in the PyPy container is incorrect
# https://github.com/docker-library/pypy/issues/21
if command -v pypy3; then
    virtualenv -p $(command -v pypy3) pypy-env
elif command -v pypy; then
    virtualenv -p $(command -v pypy) pypy-env
fi

source activate pypy-env/bin/activate

python --version
which python

# run via python -m pip to be sure; creating a virtualenv fails with the PyPy2 container
pip install --upgrade pip
pip install --extra-index https://antocuni.github.io/pypy-wheels/ubuntu numpy Cython Tempita pytest
pip install scipy==1.1.0rc1
pip install -e .


pytest sklearn/
