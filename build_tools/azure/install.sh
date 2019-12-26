#!/bin/bash

set -e

UNAMESTR=`uname`

make_conda() {
    TO_INSTALL="$@"
    conda create -n $VIRTUALENV --yes $TO_INSTALL
    source activate $VIRTUALENV
}

version_ge() {
    # The two version numbers are separated with a new line is piped to sort
    # -rV. The -V activates for version number sorting and -r sorts in
    # descending order. If the first argument is the top element of the sort, it
    # is greater than or equal to the second argument.
    test "$(printf "${1}\n${2}" | sort -rV | head -n 1)" == "$1"
}

if [[ "$DISTRIB" == "conda" ]]; then

    TO_INSTALL="python=$PYTHON_VERSION pip \
                numpy=$NUMPY_VERSION scipy=$SCIPY_VERSION \
                cython=$CYTHON_VERSION joblib=$JOBLIB_VERSION\
                blas[build=$BLAS]"

    if [[ -n "$PANDAS_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL pandas=$PANDAS_VERSION"
    fi

    if [[ -n "$PYAMG_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL pyamg=$PYAMG_VERSION"
    fi

    if [[ -n "$PILLOW_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL pillow=$PILLOW_VERSION"
    fi

    if [[ -n "$MATPLOTLIB_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL matplotlib=$MATPLOTLIB_VERSION"
    fi

    if [[ "$UNAMESTR" == "Darwin" ]]; then
        if [[ "$SKLEARN_TEST_NO_OPENMP" != "true" ]]; then
            # on macOS, install an OpenMP-enabled clang/llvm from conda-forge.
            TO_INSTALL="$TO_INSTALL conda-forge::compilers \
                        conda-forge::llvm-openmp"
        fi
    fi

    # Old packages coming from the 'free' conda channel have been removed but
    # we are using them for testing Python 3.5. See
    # https://www.anaconda.com/why-we-removed-the-free-channel-in-conda-4-7/
    # for more details. restore_free_channel is defined starting from conda 4.7
    conda_version=$(conda -V | awk '{print $2}')
    if version_ge "$conda_version" "4.7.0" && [[ "$PYTHON_VERSION" == "3.5" ]]; then
        conda config --set restore_free_channel true
    fi

	make_conda $TO_INSTALL

    if [[ "$PYTEST_VERSION" == "*" ]]; then
        python -m pip install pytest
    else
        python -m pip install pytest=="$PYTEST_VERSION"
    fi

    if [[ "$PYTHON_VERSION" == "*" ]]; then
        python -m pip install pytest-xdist
    fi

elif [[ "$DISTRIB" == "ubuntu" ]]; then
    sudo add-apt-repository --remove ppa:ubuntu-toolchain-r/test
    sudo apt-get update
    sudo apt-get install python3-scipy python3-matplotlib libatlas3-base libatlas-base-dev libatlas-dev python3-virtualenv
    python3 -m virtualenv --system-site-packages --python=python3 $VIRTUALENV
    source $VIRTUALENV/bin/activate
    python -m pip install pytest==$PYTEST_VERSION pytest-cov cython joblib==$JOBLIB_VERSION
elif [[ "$DISTRIB" == "ubuntu-32" ]]; then
    apt-get update
    apt-get install -y python3-dev python3-scipy python3-matplotlib libatlas3-base libatlas-base-dev libatlas-dev python3-virtualenv
    python3 -m virtualenv --system-site-packages --python=python3 $VIRTUALENV
    source $VIRTUALENV/bin/activate
    python -m pip install pytest==$PYTEST_VERSION pytest-cov cython joblib==$JOBLIB_VERSION
elif [[ "$DISTRIB" == "conda-pip-latest" ]]; then
    # Since conda main channel usually lacks behind on the latest releases,
    # we use pypi to test against the latest releases of the dependencies.
    # conda is still used as a convenient way to install Python and pip.
    make_conda "python=$PYTHON_VERSION"
    python -m pip install -U pip
    python -m pip install numpy scipy cython joblib
    python -m pip install pytest==$PYTEST_VERSION pytest-cov pytest-xdist
    python -m pip install pandas matplotlib pyamg
fi

if [[ "$COVERAGE" == "true" ]]; then
    python -m pip install coverage codecov pytest-cov
fi

if [[ "$TEST_DOCSTRINGS" == "true" ]]; then
    # numpydoc requires sphinx
    # FIXME: until jinja2 2.10.2 is released with a fix the import station for
    # collections.abc so as to not raise a spurious deprecation warning
    python -m pip install sphinx==2.1.2
    python -m pip install numpydoc
fi

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python -c "\
try:
    import pandas
    print('pandas %s' % pandas.__version__)
except ImportError:
    print('pandas not installed')
"
python -m pip list

# Use setup.py instead of `pip install -e .` to be able to pass the -j flag
# to speed-up the building multicore CI machines.
python setup.py build_ext --inplace -j 3
python setup.py develop
