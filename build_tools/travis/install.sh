#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD

# Travis clone scikit-learn/scikit-learn repository in to a local repository.
# We use a cached directory with three scikit-learn repositories (one for each
# matrix entry) from which we pull from local Travis repository. This allows
# us to keep build artefact for gcc + cython, and gain time

set -e

# Fix the compilers to workaround avoid having the Python 3.4 build
# lookup for g++44 unexpectedly.
export CC=gcc
export CXX=g++

echo 'List files from cached directories'
echo 'pip:'
ls $HOME/.cache/pip


if [[ "$DISTRIB" == "conda" ]]; then
    # Deactivate the travis-provided virtual environment and setup a
    # conda-based environment instead
    deactivate

    # Use the miniconda installer for faster download / install of conda
    # itself
    pushd .
    cd
    mkdir -p download
    cd download
    echo "Cached in $HOME/download :"
    ls -l
    echo
    if [[ ! -f miniconda.sh ]]
        then
        wget http://repo.continuum.io/miniconda/Miniconda-3.6.0-Linux-x86_64.sh \
            -O miniconda.sh
        fi
    chmod +x miniconda.sh && ./miniconda.sh -b
    cd ..
    export PATH=/home/travis/miniconda/bin:$PATH
    conda update --yes conda
    popd

    # Configure the conda environment and put it in the path using the
    # provided versions
    if [[ "$INSTALL_MKL" == "true" ]]; then
        conda create -n testenv --yes python=$PYTHON_VERSION pip nose pytest \
            numpy=$NUMPY_VERSION scipy=$SCIPY_VERSION \
            mkl cython=$CYTHON_VERSION \
            ${PANDAS_VERSION+pandas=$PANDAS_VERSION}
            
    else
        conda create -n testenv --yes python=$PYTHON_VERSION pip nose pytest \
            numpy=$NUMPY_VERSION scipy=$SCIPY_VERSION \
            nomkl cython=$CYTHON_VERSION \
            ${PANDAS_VERSION+pandas=$PANDAS_VERSION}
    fi
    source activate testenv

    # Install nose-timer via pip
    pip install nose-timer

elif [[ "$DISTRIB" == "ubuntu" ]]; then
    # At the time of writing numpy 1.9.1 is included in the travis
    # virtualenv but we want to use the numpy installed through apt-get
    # install.
    deactivate
    # Create a new virtualenv using system site packages for python, numpy
    # and scipy
    virtualenv --system-site-packages testvenv
    source testvenv/bin/activate
    pip install nose nose-timer cython

elif [[ "$DISTRIB" == "scipy-dev-wheels" ]]; then
    # Set up our own virtualenv environment to avoid travis' numpy.
    # This venv points to the python interpreter of the travis build
    # matrix.
    virtualenv --python=python ~/testvenv
    source ~/testvenv/bin/activate
    pip install --upgrade pip setuptools

    echo "Installing numpy and scipy master wheels"
    dev_url=https://7933911d6844c6c53a7d-47bd50c35cd79bd838daf386af554a83.ssl.cf2.rackcdn.com
    pip install --pre --upgrade --timeout=60 -f $dev_url numpy scipy
    pip install nose nose-timer cython
fi

if [[ "$COVERAGE" == "true" ]]; then
    pip install coverage codecov
fi

if [[ "$SKIP_TESTS" == "true" ]]; then
    echo "No need to build scikit-learn when not running the tests"
else
    # Build scikit-learn in the install.sh script to collapse the verbose
    # build output in the travis output when it succeeds.
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
fi

if [[ "$RUN_FLAKE8" == "true" ]]; then
    # flake8 version is temporarily set to 2.5.1 because the next
    # version available on conda (3.3.0) has a bug that checks non
    # python files and cause non meaningful flake8 errors
    conda install --yes flake8=2.5.1
fi
