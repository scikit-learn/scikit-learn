#!/bin/bash

set -e

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
    sudo apt-get install python3-scipy libatlas3-base libatlas-base-dev libatlas-dev python3-virtualenv
    python3 -m virtualenv --system-site-packages --python=python3 $VIRTUALENV
    source $VIRTUALENV/bin/activate
    python -m pip install pytest pytest-cov cython joblib==$JOBLIB_VERSION

elif [[ "$DISTRIB" == "32bit" ]]; then
    # TODO: Also choose python version?
    TO_INSTALL="pytest pytest-cov \
                numpy scipy \
                cython"

    if [[ "$INSTALL_MKL" == "true" ]]; then
        TO_INSTALL="$TO_INSTALL mkl"
    else
        TO_INSTALL="$TO_INSTALL nomkl"
    fi

    # Don't specify the versions for now for 32-bit
    if [[ -n "$PANDAS_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL pandas"
    fi

    if [[ -n "$PYAMG_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL pyamg"
    fi

    if [[ -n "$PILLOW_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL pillow"
    fi

    if [[ -n "$JOBLIB_VERSION" ]]; then
        TO_INSTALL="$TO_INSTALL joblib"
    fi

    echo "TO_INSTALL: $TO_INSTALL"

    # update ubuntu
    apt-get update && apt-get upgrade

    # TODO: make python version a param
    conda create -n python python=3.5 -y
    echo "source activate python" >> ~/.bashrc

    python -m pip install --upgrade pip
    python -m pip install $TO_INSTALL
fi


if [[ "$COVERAGE" == "true" ]]; then
    python -m pip install coverage codecov
fi

if [[ "$TEST_DOCSTRINGS" == "true" ]]; then
    python -m pip install sphinx numpydoc  # numpydoc requires sphinx
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
export BITS=`python -c 'import struct; print(8 * struct.calcsize("P"))'`
echo "Architecture: $BITS bits"

pip list
python setup.py develop
