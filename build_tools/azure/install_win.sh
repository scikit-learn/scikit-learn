#!/bin/bash

set -e
set -x

# defines the get_dep and show_installed_libraries functions
source build_tools/shared.sh

if [[ -n "$LOCK_FILE" ]]; then
    if [[ "$DISTRIB" == "conda" ]]; then
        conda update -n base conda -y
        conda install pip -y
        conda list
        # FIXME install conda-lock dev version with a fixed commit while waiting
        # for the release
        python -m pip install git+https://github.com/conda-incubator/conda-lock@4203aef
        conda-lock install --name $VIRTUALENV $LOCK_FILE
        source activate $VIRTUALENV
    fi
else
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
fi

show_installed_libraries

# Build scikit-learn
python setup.py bdist_wheel

# Install the generated wheel package to test it
pip install --pre --no-index --find-links dist scikit-learn
