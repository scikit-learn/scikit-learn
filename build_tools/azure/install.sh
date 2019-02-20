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

if [[ "$DISTRIB" == "ubuntu" ]]; then
    sudo apt-get install python3-scipy libatlas3-base libatlas-base-dev libatlas-dev
    pip install virtualenv
    virtualenv --system-site-packages --python=python3 $VIRTUALENV_DIR
    source $VIRTUALENV_DIR/bin/activate
    pip install pytest pytest-cov cython joblib==$JOBLIB_VERSION
fi

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python -c "\
try:
    import pandas
    print('pandas %s' % pandas.__version__)
except ImportError:
    pass
"
python setup.py develop
