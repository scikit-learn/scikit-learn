#!/bin/bash

set -e
set -x

UNAMESTR=`uname`

if [[ "$DISTRIB" == "conda-mamba-pypy3" ]]; then
    # condaforge/mambaforge-pypy3 needs compilers
    apt-get -yq update
    apt-get -yq install build-essential
fi

make_conda() {
    TO_INSTALL="$@"
    if [[ "$DISTRIB" == *"mamba"* ]]; then
        mamba create -n $VIRTUALENV --yes $TO_INSTALL
    else
        conda config --show
        conda create -n $VIRTUALENV --yes $TO_INSTALL
    fi
    source activate $VIRTUALENV
}

setup_ccache() {
    echo "Setting up ccache with CCACHE_DIR=${CCACHE_DIR}"
    mkdir /tmp/ccache/
    which ccache
    for name in gcc g++ cc c++ clang clang++ i686-linux-gnu-gcc i686-linux-gnu-c++ x86_64-linux-gnu-gcc x86_64-linux-gnu-c++ x86_64-apple-darwin13.4.0-clang x86_64-apple-darwin13.4.0-clang++; do
      ln -s $(which ccache) "/tmp/ccache/${name}"
    done
    export PATH="/tmp/ccache/:${PATH}"
    ccache -M 256M
}

# defines the get_dep and show_installed_libraries functions
source build_tools/shared.sh

if [[ "$DISTRIB" == "conda" || "$DISTRIB" == *"mamba"* ]]; then

    if [[ "$CONDA_CHANNEL" != "" ]]; then
        TO_INSTALL="--override-channels -c $CONDA_CHANNEL"
    else
        TO_INSTALL=""
    fi

    if [[ "$DISTRIB" == *"pypy"* ]]; then
        TO_INSTALL="$TO_INSTALL pypy"
    else
        TO_INSTALL="$TO_INSTALL python=$PYTHON_VERSION"
    fi

    TO_INSTALL="$TO_INSTALL ccache pip blas[build=$BLAS]"

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
            TO_INSTALL="$TO_INSTALL compilers llvm-openmp"
        else
            # Without openmp, we use the system clang. Here we use /usr/bin/ar
            # instead because llvm-ar errors
            export AR=/usr/bin/ar
        fi
    else
        # FIXME: temporary fix to link against system libraries on linux
        export LDFLAGS="$LDFLAGS -Wl,--sysroot=/"
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

elif [[ "$DISTRIB" == "debian-32" ]]; then
    apt-get update
    apt-get install -y python3-dev python3-numpy python3-scipy python3-matplotlib libatlas3-base libatlas-base-dev python3-virtualenv python3-pandas ccache

    python3 -m virtualenv --system-site-packages --python=python3 $VIRTUALENV
    source $VIRTUALENV/bin/activate
    setup_ccache
    python -m pip install $(get_dep cython $CYTHON_VERSION) \
                          $(get_dep joblib $JOBLIB_VERSION)

elif [[ "$DISTRIB" == "conda-pip-latest" ]]; then
    # FIXME: temporary fix to link against system libraries on linux
    export LDFLAGS="$LDFLAGS -Wl,--sysroot=/"
    # Since conda main channel usually lacks behind on the latest releases,
    # we use pypi to test against the latest releases of the dependencies.
    # conda is still used as a convenient way to install Python and pip.
    make_conda "ccache python=$PYTHON_VERSION"
    setup_ccache
    python -m pip install -U pip

    # Do not build scikit-image from source because it is an optional dependency
    python -m pip install --only-binary :all: scikit-image || true

    python -m pip install pandas matplotlib pyamg
    # do not install dependencies for lightgbm since it requires scikit-learn.
    python -m pip install "lightgbm>=3.0.0" --no-deps
elif [[ "$DISTRIB" == "conda-pip-scipy-dev" ]]; then
    # FIXME: temporary fix to link against system libraries on linux
    export LDFLAGS="$LDFLAGS -Wl,--sysroot=/"
    make_conda "ccache python=$PYTHON_VERSION"
    python -m pip install -U pip
    echo "Installing numpy and scipy master wheels"
    dev_anaconda_url=https://pypi.anaconda.org/scipy-wheels-nightly/simple
    pip install --pre --upgrade --timeout=60 --extra-index $dev_anaconda_url numpy pandas scipy
    pip install --pre cython
    setup_ccache
    echo "Installing joblib master"
    pip install https://github.com/joblib/joblib/archive/master.zip
    echo "Installing pillow master"
    pip install https://github.com/python-pillow/Pillow/archive/main.zip
fi

python -m pip install $(get_dep threadpoolctl $THREADPOOLCTL_VERSION) \
                      $(get_dep pytest $PYTEST_VERSION) \
                      $(get_dep pytest-xdist $PYTEST_XDIST_VERSION)

if [[ "$COVERAGE" == "true" ]]; then
    # XXX: coverage is temporary pinned to 6.2 because 6.3 is not fork-safe
    # cf. https://github.com/nedbat/coveragepy/issues/1310
    python -m pip install codecov pytest-cov coverage==6.2
fi

if [[ "$TEST_DOCSTRINGS" == "true" ]]; then
    # numpydoc requires sphinx
    python -m pip install sphinx
    # TODO: update the docstring checks to be compatible with new
    # numpydoc versions
    python -m pip install "numpydoc<1.2"
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

show_installed_libraries

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
