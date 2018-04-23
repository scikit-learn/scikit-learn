#!/usr/bin/env bash
set -x
set -e


apt-get -yq update
apt-get -yq install libatlas-dev libatlas-base-dev liblapack-dev gfortran

# python executable is incorrect by default
if command -v pypy3; then ln -s $(command -v pypy3) /usr/local/bin/python; fi

python --version
which python

python -m pip install --upgrade pip

which pip


pip install --extra-index https://antocuni.github.io/pypy-wheels/ubuntu numpy Cython Tempita
pip install scipy==1.1.0rc1
pip install -e .


pytest sklearn/
