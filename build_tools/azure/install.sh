#!/bin/bash

set -e
set -x

UNAMESTR=`uname`

make_conda() {
    TO_INSTALL="$@"
    conda create -n $VIRTUALENV --yes $TO_INSTALL
    source activate $VIRTUALENV
}

setup_ccache() {
    echo "Setting up ccache"
    mkdir /tmp/ccache/
    which ccache
    for name in gcc g++ cc c++ x86_64-linux-gnu-gcc x86_64-linux-gnu-c++; do
      ln -s $(which ccache) "/tmp/ccache/${name}"
    done
    export PATH="/tmp/ccache/:${PATH}"
    ccache -M 256M
}

# imports get_dep
source build_tools/shared.sh

if [[ "$DISTRIB" == "conda" ]]; then

    if [[ "$CONDA_CHANNEL" != "" ]]; then
        TO_INSTALL="-c $CONDA_CHANNEL"
    else
        TO_INSTALL=""
    fi

    TO_INSTALL="$TO_INSTALL python=$PYTHON_VERSION ccache pip blas[build=$BLAS]"

    TO_INSTALL="$TO_INSTALL $(get_dep numpy $NUMPY_VERSION)"
    TO_INSTALL="$TO_INSTALL $(get_dep scipy $SCIPY_VERSION)"
    TO_INSTALL="$TO_INSTALL $(get_dep cython $CYTHON_VERSION)"
    TO_INSTALL="$TO_INSTALL $(get_dep joblib $JOBLIB_VERSION)"
    TO_INSTALL="$TO_INSTALL $(get_dep pandas $PANDAS_VERSION)"
    TO_INSTALL="$TO_INSTALL $(get_dep pyamg $PYAMG_VERSION)"
    TO_INSTALL="$TO_INSTALL $(get_dep Pillow $PILLOW_VERSION)"
    TO_INSTALL="$TO_INSTALL $(get_dep matplotlib $MATPLOTLIB_VERSION)"

    if [[ "$UNAMESTR" == "Darwin" ]]; then
        if [[ "$SKLEARN_TEST_NO_OPENMP" != "true" ]]; then
            # on macOS, install an OpenMP-enabled clang/llvm from conda-forge.
            # TODO: Remove !=1.1.0 when the following is fixed:
            # sklearn/svm/_libsvm.cpython-38-darwin.so,
            # 2): Symbol not found: _svm_check_parameter error
            TO_INSTALL="$TO_INSTALL compilers>=1.0.4,!=1.1.0 llvm-openmp"
        fi
    fi
	make_conda $TO_INSTALL
    setup_ccache

elif [[ "$DISTRIB" == "ubuntu" ]]; then
    sudo add-apt-repository --remove ppa:ubuntu-toolchain-r/test
    sudo apt-get update
    sudo apt-get install python3-scipy python3-matplotlib libatlas3-base libatlas-base-dev python3-virtualenv ccache
    python3 -m virtualenv --system-site-packages --python=python3 $VIRTUALENV
    source $VIRTUALENV/bin/activate
    setup_ccache
    python -m pip install $(get_dep cython $CYTHON_VERSION) \
                          $(get_dep joblib $JOBLIB_VERSION)

elif [[ "$DISTRIB" == "ubuntu-32" ]]; then
    apt-get update
    apt-get install -y python3-dev python3-scipy python3-matplotlib libatlas3-base libatlas-base-dev python3-virtualenv python3-pandas ccache

    python3 -m virtualenv --system-site-packages --python=python3 $VIRTUALENV
    source $VIRTUALENV/bin/activate
    setup_ccache
    python -m pip install $(get_dep cython $CYTHON_VERSION) \
                          $(get_dep joblib $JOBLIB_VERSION)

elif [[ "$DISTRIB" == "conda-pip-latest" ]]; then
    # Since conda main channel usually lacks behind on the latest releases,
    # we use pypi to test against the latest releases of the dependencies.
    # conda is still used as a convenient way to install Python and pip.
    make_conda "ccache python=$PYTHON_VERSION"
    setup_ccache
    python -m pip install -U pip

    # Do not build scikit-image from source because it is an optional dependency
    python -m pip install --only-binary :all: scikit-image || true

    python -m pip install pandas matplotlib pyamg
    # do not install dependencies for lightgbm since it requires scikit-learn
    # and install a version less than 3.0.0 until the issue #18316 is solved.
    python -m pip install "lightgbm<3.0.0" --no-deps
elif [[ "$DISTRIB" == "conda-pip-scipy-dev" ]]; then
    make_conda "ccache python=$PYTHON_VERSION"
    python -m pip install -U pip
    echo "Installing numpy and scipy master wheels"
    dev_anaconda_url=https://pypi.anaconda.org/scipy-wheels-nightly/simple
    pip install --pre --upgrade --timeout=60 --extra-index $dev_anaconda_url numpy pandas

    # issue with metadata in scipy dev builds https://github.com/scipy/scipy/issues/13196
    # --use-deprecated=legacy-resolver needs to be included
    pip install --pre --upgrade --timeout=60 --extra-index $dev_anaconda_url scipy --use-deprecated=legacy-resolver
    pip install --pre cython
    setup_ccache
    echo "Installing joblib master"
    pip install https://github.com/joblib/joblib/archive/master.zip
    echo "Installing pillow master"
    pip install https://github.com/python-pillow/Pillow/archive/master.zip
fi

python -m pip install $(get_dep threadpoolctl $THREADPOOLCTL_VERSION) \
                      $(get_dep pytest $PYTEST_VERSION) \
                      $(get_dep pytest-xdist $PYTEST_XDIST_VERSION)

if [[ "$COVERAGE" == "true" ]]; then
    python -m pip install codecov pytest-cov
fi

if [[ "$PYTEST_XDIST_VERSION" != "none" ]]; then
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
# Set parallelism to 3 to overlap IO bound tasks with CPU bound tasks on CI
# workers with 2 cores when building the compiled extensions of scikit-learn.
export SKLEARN_BUILD_PARALLEL=3

python -m pip list
if [[ "$DISTRIB" == "conda-pip-latest" ]]; then
    # Check that pip can automatically build scikit-learn with the build
    # dependencies specified in pyproject.toml using an isolated build
    # environment:
    pip install --verbose --editable .
else
    if [[ "$BUILD_WITH_ICC" == "true" ]]; then
        wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
        sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
        rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
        sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
        sudo apt-get update
        sudo apt-get install intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic
        source /opt/intel/oneapi/setvars.sh

        # The "build_clib" command is implicitly used to build "libsvm-skl".
        # To compile with a different compiler, we also need to specify the
        # compiler for this command
        python setup.py build_ext --compiler=intelem -i build_clib --compiler=intelem
    fi
    # Use the pre-installed build dependencies and build directly in the
    # current environment.
    python setup.py develop
fi
ccache -s
