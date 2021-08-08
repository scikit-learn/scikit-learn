#!/bin/bash

set -e
set -x

UNAMESTR=`uname`

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

get_build_type() {
    if [ -z "$CIRCLE_SHA1" ]
    then
        echo SKIP: undefined CIRCLE_SHA1
        return
    fi
    commit_msg=$(git log --format=%B -n 1 $CIRCLE_SHA1)
    if [[ "$commit_msg" =~ \[cd\ build\] ]]
    then
        echo BUILD: [cd build] marker found
        BUILD_WHEEL=true
        return
    fi
}

get_build_type

# imports get_dep
source build_tools/shared.sh

sudo add-apt-repository --remove ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install python3-virtualenv ccache
python3 -m virtualenv --system-site-packages --python=python3 testenv
source testenv/bin/activate
pip install --upgrade pip
if [[ "$BUILD_WHEEL" == "true" ]]
then
    for pyversion in 37 38 39; do
        pywheel=`echo "cp"$pyversion"-manylinux_aarch64"`
        export CIBW_BUILD=$pywheel
        python -m pip install cibuildwheel
        python -m cibuildwheel --output-dir wheelhouse
    done
else
    setup_ccache
    python -m pip install $(get_dep cython $CYTHON_VERSION) \
                          $(get_dep joblib $JOBLIB_VERSION)
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

    # Set parallelism to 3 to overlap IO bound tasks with CPU bound tasks on CI
    # workers with 2 cores when building the compiled extensions of scikit-learn.
    export SKLEARN_BUILD_PARALLEL=3

    python -m pip list
    pip install --verbose --editable .
    ccache -s
    python -c "import sklearn; sklearn.show_versions()"
    python -m threadpoolctl --import sklearn
    python -m pytest sklearn
fi
