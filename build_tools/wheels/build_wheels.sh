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
fi


if [[ "$GITHUB_EVENT_NAME" == "schedule" || "$CIRRUS_CRON" == "nightly" ]]; then
    # Nightly build:  See also `../github/upload_anaconda.sh` (same branching).
    # To help with NumPy 2.0 transition, ensure that we use the NumPy 2.0
    # nightlies.  This lives on the edge and opts-in to all pre-releases.
    # That could be an issue, in which case no-build-isolation and a targeted
    # NumPy install may be necessary, instead.
    export CIBW_BUILD_FRONTEND='pip; args: --pre --extra-index-url "https://pypi.anaconda.org/scientific-python-nightly-wheels/simple"'
fi

# The version of the built dependencies are specified
# in the pyproject.toml file, while the tests are run
# against the most recent version of the dependencies

python -m pip install cibuildwheel
python -m cibuildwheel --output-dir wheelhouse
