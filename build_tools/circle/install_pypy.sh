#!/usr/bin/env bash
set -x
set -e


apt-get -yq update
apt-get -yq install libatlas-dev libatlas-base-dev liblapack-dev gfortran

pip install virtualenv

# python executable in the PyPy container is incorrect
# https://github.com/docker-library/pypy/issues/21
if command -v pypy3; then
    ln -s $(command -v pypy3) /usr/local/bin/python
elif command -v pypy; then
    ln -s $(command -v pypy) /usr/local/bin/python
fi

python --version
which python

# run via python -m pip to be sure; creating a virtualenv fails with the PyPy2 container
python -m pip install --upgrade pip
python -m pip install --extra-index https://antocuni.github.io/pypy-wheels/ubuntu numpy Cython Tempita pytest
python -m pip install scipy==1.1.0rc1
python -m pip install -e .


python -m pytest sklearn/
