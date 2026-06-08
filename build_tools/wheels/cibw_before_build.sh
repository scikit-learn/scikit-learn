#!/bin/bash

set -exo pipefail

PROJECT_DIR="$1"
LICENSE_FILE="$PROJECT_DIR/COPYING"

echo "" >>"$LICENSE_FILE"
echo "----" >>"$LICENSE_FILE"
echo "" >>"$LICENSE_FILE"

if [[ $RUNNER_OS == "Linux" ]]; then
    cat $PROJECT_DIR/build_tools/wheels/LICENSE_linux.txt >>"$LICENSE_FILE"
elif [[ $RUNNER_OS == "macOS" ]]; then
    cat $PROJECT_DIR/build_tools/wheels/LICENSE_macos.txt >>"$LICENSE_FILE"
elif [[ $RUNNER_OS == "Windows" ]]; then
    cat $PROJECT_DIR/build_tools/wheels/LICENSE_windows.txt >>"$LICENSE_FILE"
fi

# Set environment variables to make our wheel build easier to reproduce byte
# for byte from source. See https://reproducible-builds.org/. The long term
# motivation would be to be able to detect supply chain attacks.
#
# In particular we set SOURCE_DATE_EPOCH to the commit date of the last commit.
#
# XXX: setting those environment variables is not enough. See the following
# issue for more details on what remains to do:
# https://github.com/scikit-learn/scikit-learn/issues/28151
echo SOURCE_DATE_EPOCH=$(git log -1 --pretty=%ct) >> "$GITHUB_ENV"
echo PYTHONHASHSEED=0 >> "$GITHUB_ENV"

# OpenMP is not present on macOS by default
if [[ $(uname) == "Darwin" ]]; then
    # Make sure to use a libomp version binary compatible with the oldest
    # supported version of the macos SDK as libomp will be vendored into the
    # scikit-learn wheels for macos.

    if [[ "$CIBW_BUILD" == *-macosx_arm64 ]]; then
        # SciPy requires 12.0 on arm to prevent kernel panics
        # https://github.com/scipy/scipy/issues/14688
        # We use the same deployment target to match SciPy.
        export MACOSX_DEPLOYMENT_TARGET=12.0
        OPENMP_URL="https://anaconda.org/conda-forge/llvm-openmp/11.1.0/download/osx-arm64/llvm-openmp-11.1.0-hf3c4609_1.tar.bz2"
    else
        export MACOSX_DEPLOYMENT_TARGET=10.9
        OPENMP_URL="https://anaconda.org/conda-forge/llvm-openmp/11.1.0/download/osx-64/llvm-openmp-11.1.0-hda6cdc1_1.tar.bz2"
    fi

    conda create -n build $OPENMP_URL
    PREFIX="$HOME/miniconda3/envs/build"

    echo C=/usr/bin/clang >> "$GITHUB_ENV"
    echo CXX=/usr/bin/clang++ >> "$GITHUB_ENV"
    echo CPPFLAGS="$CPPFLAGS -Xpreprocessor -fopenmp" >> "$GITHUB_ENV"
    echo CFLAGS="$CFLAGS -I$PREFIX/include" >> "$GITHUB_ENV"
    echo CXXFLAGS="$CXXFLAGS -I$PREFIX/include" >> "$GITHUB_ENV"
    echo LDFLAGS="$LDFLAGS -Wl,-rpath,$PREFIX/lib -L$PREFIX/lib -lomp" >> "$GITHUB_ENV"
fi
