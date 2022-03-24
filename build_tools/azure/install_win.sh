#!/bin/bash

set -e
set -x

if [[ "$PYTHON_ARCH" == "64" ]]; then
    conda create -n $VIRTUALENV -q -y python=$PYTHON_VERSION numpy scipy cython matplotlib wheel pillow joblib

    source activate $VIRTUALENV

    pip install threadpoolctl

    if [[ "$PYTEST_VERSION" == "*" ]]; then
        pip install pytest
    else
        pip install pytest==$PYTEST_VERSION
    fi
else
    pip install numpy scipy cython pytest wheel pillow joblib threadpoolctl
fi

if [[ "$PYTEST_XDIST_VERSION" != "none" ]]; then
    pip install pytest-xdist
fi

if [[ "$COVERAGE" == "true" ]]; then
    # XXX: coverage is temporary pinned to 6.2 because 6.3 is not fork-safe
    # cf. https://github.com/nedbat/coveragepy/issues/1310
    pip install coverage codecov pytest-cov coverage==6.2
fi

python --version
pip --version

# Build scikit-learn
python setup.py bdist_wheel

# Install the generated wheel package to test it
pip install --pre --no-index --find-links dist scikit-learn
