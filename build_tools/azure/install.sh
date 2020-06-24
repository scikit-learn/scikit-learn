#!/bin/bash

set -e
set -x

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

    if [[ -n "$SCIKIT_IMAGE_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL scikit-image=$SCIKIT_IMAGE_VERSION"
    fi

    if [[ -n "$MATPLOTLIB_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL matplotlib=$MATPLOTLIB_VERSION"
    fi

    if [[ "$UNAMESTR" == "Darwin" ]]; then
        if [[ "$SKLEARN_TEST_NO_OPENMP" != "true" ]]; then
            # on macOS, install an OpenMP-enabled clang/llvm from conda-forge.
            TO_INSTALL="$TO_INSTALL conda-forge::compilers>=1.0.4 \
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

    pip install threadpoolctl==$THREADPOOLCTL_VERSION

    if [[ "$PYTEST_VERSION" == "*" ]]; then
        python -m pip install pytest
    else
        python -m pip install pytest=="$PYTEST_VERSION"
    fi

elif [[ "$DISTRIB" == "ubuntu" ]]; then
    sudo add-apt-repository --remove ppa:ubuntu-toolchain-r/test
    sudo apt-get update
    sudo apt-get install python3-scipy python3-matplotlib libatlas3-base libatlas-base-dev python3-virtualenv
    python3 -m virtualenv --system-site-packages --python=python3 $VIRTUALENV
    source $VIRTUALENV/bin/activate
    python -m pip install pytest==$PYTEST_VERSION pytest-cov cython joblib==$JOBLIB_VERSION threadpoolctl==$THREADPOOLCTL_VERSION
elif [[ "$DISTRIB" == "ubuntu-32" ]]; then
    apt-get update
    apt-get install -y python3-dev python3-scipy python3-matplotlib libatlas3-base libatlas-base-dev python3-virtualenv
    python3 -m virtualenv --system-site-packages --python=python3 $VIRTUALENV
    source $VIRTUALENV/bin/activate
    python -m pip install pytest==$PYTEST_VERSION pytest-cov cython joblib==$JOBLIB_VERSION threadpoolctl==$THREADPOOLCTL_VERSION
elif [[ "$DISTRIB" == "conda-pip-latest" ]]; then
    # Since conda main channel usually lacks behind on the latest releases,
    # we use pypi to test against the latest releases of the dependencies.
    # conda is still used as a convenient way to install Python and pip.
    make_conda "python=$PYTHON_VERSION"
    python -m pip install -U pip
    python -m pip install pytest==$PYTEST_VERSION pytest-cov

    python -m pip install pandas matplotlib pyamg scikit-image
    # do not install dependencies for lightgbm since it requires scikit-learn
    python -m pip install lightgbm --no-deps
elif [[ "$DISTRIB" == "conda-pip-scipy-dev" ]]; then
    make_conda "python=$PYTHON_VERSION"
    python -m pip install -U pip
    python -m pip install pytest==$PYTEST_VERSION pytest-cov
    echo "Installing numpy and scipy master wheels"
    dev_anaconda_url=https://pypi.anaconda.org/scipy-wheels-nightly/simple
    pip install --pre --upgrade --timeout=60 --extra-index $dev_anaconda_url numpy scipy pandas
    # Cython nightly build should be still fetched from the Rackspace container
    dev_rackspace_url=https://7933911d6844c6c53a7d-47bd50c35cd79bd838daf386af554a83.ssl.cf2.rackcdn.com
    pip install --pre --upgrade --timeout=60 -f $dev_rackspace_url cython
    echo "Installing joblib master"
    pip install https://github.com/joblib/joblib/archive/master.zip
    echo "Installing pillow master"
    pip install https://github.com/python-pillow/Pillow/archive/master.zip
fi

if [[ "$COVERAGE" == "true" ]]; then
    python -m pip install coverage codecov pytest-cov
fi

if [[ "$PYTEST_XDIST" == "true" ]]; then
    python -m pip install pytest-xdist
fi

if [[ "$TEST_DOCSTRINGS" == "true" ]]; then
    # numpydoc requires sphinx
    python -m pip install sphinx
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

if [[ "$DISTRIB" == "conda-pip-latest" ]]; then
    # Check that pip can automatically install the build dependencies from
    # pyproject.toml using an isolated build environment:
    pip install --verbose --editable .
else
    # Use the pre-installed build dependencies and build directly in the
    # current environment.
    # Use setup.py instead of `pip install -e .` to be able to pass the -j flag
    # to speed-up the building multicore CI machines.
    python setup.py build_ext --inplace -j 3
    python setup.py develop
fi
