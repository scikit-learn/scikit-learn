#!/bin/bash

set -e

UNAMESTR=`uname`

make_conda() {
    # Install Miniconda
    if [[ "$UNAMESTR" == 'Linux' ]]; then
        wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    elif [[ "$UNAMESTR" == 'Darwin' ]]; then
        wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
    else
        exit 1
    fi
    chmod +x miniconda.sh
    ./miniconda.sh -b
}

echo "$VIRTUALENV_DIR"

if [[ "$DISTRIB" == "ubuntu" ]]; then
    pip install virtualenv numpy==$NUMPY_VERSION scipy==$SCIPY_VERSION
    virtualenv --system-site-packages --python=python3 $VIRTUALENV_DIR
    source $VIRTUALENV_DIR/bin/activate
    pip install pytest  cython joblib==$JOBLIB_VERSION
fI

# python --version
# python -c "import numpy; print('numpy %s' % numpy.__version__)"
# python -c "import scipy; print('scipy %s' % scipy.__version__)"
# python setup.py develop
