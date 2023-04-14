#!/bin/bash

set -e
set -x

# OpenMP is not present on macOS by default
if [[ $(uname) == "Darwin" ]]; then
    # Make sure to use a libomp version binary compatible with the oldest
    # supported version of the macos SDK as libomp will be vendored into the
    # scikit-learn wheels for macos.

    if [[ "$CIBW_BUILD" == *-macosx_arm64 ]]; then
        if [[ $(uname -m) == "x86_64" ]]; then
            # arm64 builds must cross compile because the CI instance is x86
            # This turns off the computation of the test program in
            # sklearn/_build_utils/pre_build_helpers.py
            export PYTHON_CROSSENV=1
        fi
        # SciPy requires 12.0 on arm to prevent kernel panics
        # https://github.com/scipy/scipy/issues/14688
        # We use the same deployment target to match SciPy.
        export MACOSX_DEPLOYMENT_TARGET=12.0
        OPENMP_URL="https://anaconda.org/conda-forge/llvm-openmp/11.1.0/download/osx-arm64/llvm-openmp-11.1.0-hf3c4609_1.tar.bz2"
    else
        export MACOSX_DEPLOYMENT_TARGET=10.9
        OPENMP_URL="https://anaconda.org/conda-forge/llvm-openmp/11.1.0/download/osx-64/llvm-openmp-11.1.0-hda6cdc1_1.tar.bz2"
    fi

    sudo conda create -n build $OPENMP_URL
    PREFIX="$CONDA_HOME/envs/build"

    export CC=/usr/bin/clang
    export CXX=/usr/bin/clang++
    export CPPFLAGS="$CPPFLAGS -Xpreprocessor -fopenmp"
    export CFLAGS="$CFLAGS -I$PREFIX/include"
    export CXXFLAGS="$CXXFLAGS -I$PREFIX/include"
    export LDFLAGS="$LDFLAGS -Wl,-rpath,$PREFIX/lib -L$PREFIX/lib -lomp"

    if [[ $(uname -m) == "arm64" && "$CIBW_BUILD" == "cp38-macosx_arm64" ]]; then
        # Enables native building and testing for macosx arm on Python 3.8. For details see:
        # https://cibuildwheel.readthedocs.io/en/stable/faq/#macos-building-cpython-38-wheels-on-arm64
        curl -o /tmp/Python38.pkg https://www.python.org/ftp/python/3.8.10/python-3.8.10-macos11.pkg
        sudo installer -pkg /tmp/Python38.pkg -target /
        sh "/Applications/Python 3.8/Install Certificates.command"
    fi
fi

# The version of the built dependencies are specified
# in the pyproject.toml file, while the tests are run
# against the most recent version of the dependencies

python -m pip install cibuildwheel
python -m cibuildwheel --output-dir wheelhouse
