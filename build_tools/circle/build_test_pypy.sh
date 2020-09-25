#!/usr/bin/env bash
set -x
set -e

apt-get -yq update
apt-get -yq install libatlas-base-dev liblapack-dev gfortran ccache libopenblas-dev

pip install virtualenv

if command -v pypy3; then
    virtualenv -p $(command -v pypy3) pypy-env
elif command -v pypy; then
    virtualenv -p $(command -v pypy) pypy-env
fi

source pypy-env/bin/activate

python --version
which python

pip install -U pip

# pins versions to install wheel from https://antocuni.github.io/pypy-wheels/manylinux2010
pip install --extra-index-url https://antocuni.github.io/pypy-wheels/manylinux2010 numpy==1.18.0 scipy==1.3.2

# Install Cython directly
pip install https://antocuni.github.io/pypy-wheels/ubuntu/Cython/Cython-0.29.14-py3-none-any.whl
pip install sphinx numpydoc docutils joblib pillow pytest

ccache -M 512M
export CCACHE_COMPRESS=1
export PATH=/usr/lib/ccache:$PATH
export LOKY_MAX_CPU_COUNT="2"
export OMP_NUM_THREADS="1"

python setup.py build_ext --inplace -j 3
pip install --no-build-isolation -e .

# Check that Python implementation is PyPy
python - << EOL
import platform
from sklearn.utils import IS_PYPY
assert IS_PYPY is True, "platform={}!=PyPy".format(platform.python_implementation())
EOL

python -m pytest sklearn/
python -m pytest doc/sphinxext/
python -m pytest $(find doc -name '*.rst' | sort)
