#!/bin/bash

set -e

make_conda() {
    TO_INSTALL="$@"
    # Install Miniconda
    UNAMESTR=`uname`
    if [[ "$UNAMESTR" == 'Linux' ]]; then
        wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    elif [[ "$UNAMESTR" == 'Darwin' ]]; then
        wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
    else
        exit 1
    fi
    chmod +x miniconda.sh
    ./miniconda.sh -b

    export PATH=$HOME/miniconda3/bin:$PATH
    conda update --yes conda
    conda create -n $VIRTUALENV --yes $TO_INSTALL
    source activate $VIRTUALENV
}

if [[ "$DISTRIB" == "conda" ]]; then
    TO_INSTALL="python=$PYTHON_VERSION pip pytest pytest-cov \
                numpy=$NUMPY_VERSION scipy=$SCIPY_VERSION \
                cython=$CYTHON_VERSION"

    if [[ "$INSTALL_MKL" == "true" ]]; then
        TO_INSTALL="$TO_INSTALL mkl"
    else
        TO_INSTALL="$TO_INSTALL nomkl"
    fi

    if [[ -n "$PANDAS_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL pandas=$PANDAS_VERSION"
    fi

    if [[ -n "$PYAMG_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL pyamg=$PYAMG_VERSION"
    fi

    if [[ -n "$PILLOW_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL pillow=$PILLOW_VERSION"
    fi

    if [[ -n "$JOBLIB_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL joblib=$JOBLIB_VERSION"
    fi

	make_conda $TO_INSTALL

elif [[ "$DISTRIB" == "ubuntu" ]]; then
    set -x
    sudo apt-get install python3-scipy libatlas3-base libatlas-base-dev libatlas-dev python3-virtualenv
    python3 -c "import numpy; print('numpy %s' % numpy.__version__)"
    virtualenv --system-site-packages --python=python3 $VIRTUALENV
    source $VIRTUALENV/bin/activate
    pip install pytest pytest-cov cython joblib==$JOBLIB_VERSION
    set +x
fi

pip list
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
