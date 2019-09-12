#!/bin/bash

set -e
set -x

UNAMESTR=`uname`

if [[ "$UNAMESTR" == "Darwin" ]]; then
    # install OpenMP not present by default on osx
    HOMEBREW_NO_AUTO_UPDATE=1 brew install libomp

    # enable OpenMP support for Apple-clang
    export CC=/usr/bin/clang
    export CXX=/usr/bin/clang++
    export CPPFLAGS="$CPPFLAGS -Xpreprocessor -fopenmp"
    export CFLAGS="$CFLAGS -I/usr/local/opt/libomp/include"
    export CXXFLAGS="$CXXFLAGS -I/usr/local/opt/libomp/include"
    export LDFLAGS="$LDFLAGS -L/usr/local/opt/libomp/lib -lomp"
    export DYLD_LIBRARY_PATH=/usr/local/opt/libomp/lib
fi

make_conda() {
    TO_INSTALL="$@"
    conda create -n $VIRTUALENV --yes $TO_INSTALL
    source activate $VIRTUALENV
}

if [[ "$DISTRIB" == "conda" ]]; then

    TO_INSTALL="python=$PYTHON_VERSION pip pytest=$PYTEST_VERSION \
                pytest-cov numpy=$NUMPY_VERSION scipy=$SCIPY_VERSION \
                cython=$CYTHON_VERSION joblib=$JOBLIB_VERSION"

    # TODO: Remove when openssl gets fixed. This is a temporary fix to get
    # `fetch_*` to work on the CI.
    TO_INSTALL="$TO_INSTALL openssl<=1.1.1c"

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

    if [[ -n "$MATPLOTLIB_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL matplotlib=$MATPLOTLIB_VERSION"
    fi

    # Old packages coming from the 'free' conda channel have been removed but
    # we are using them for testing Python 3.5. See
    # https://www.anaconda.com/why-we-removed-the-free-channel-in-conda-4-7/
    # for more details. For Python 3.5 we use the conda-forge channel
    # as a workaround.

    make_conda $TO_INSTALL
    if [[ "$PYTHON_VERSION" == "*" ]]; then
        pip install pytest-xdist
    fi

elif [[ "$DISTRIB" == "ubuntu" ]]; then
    sudo add-apt-repository --remove ppa:ubuntu-toolchain-r/test
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
elif [[ "$DISTRIB" == "conda-latest" ]]; then
    # since conda main channel usually lacks behind on the latest releases,
    # we use pypi to test against the latest releases of the dependencies.
    make_conda "python=$PYTHON_VERSION" "openssl<=1.1.1c"
    # TODO: Remove openssl ssl pinning above. This is a temporary fix to get
    # `fetch_*` to work on the CI.
    python -m pip install numpy scipy joblib cython
    python -m pip install pytest==$PYTEST_VERSION pytest-cov pytest-xdist
    python -m pip install pandas matplotlib pyamg pillow
fi

if [[ "$COVERAGE" == "true" ]]; then
    python -m pip install coverage codecov
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
pip list
python setup.py build_ext --inplace -j 3
python setup.py develop
