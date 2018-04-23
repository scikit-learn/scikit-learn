#!/usr/bin/env bash
set -x
set -e

python --version
which python
which pip

#if command -v pypy3; then ln -s $(command -v pypy3) /usr/local/bin/python; fi

apt-get -yq update
apt-get -y install libatlas-dev libatlas-base-dev liblapack-dev gfortran

pip install --extra-index https://antocuni.github.io/pypy-wheels/ubuntu numpy Cython Tempita
pip install scipy==1.1.0rc1
pip install -e .
