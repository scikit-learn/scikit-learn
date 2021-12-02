#!/bin/bash

set -e
set -x

# OpenMP is not present on macOS by default
if [[ "$RUNNER_OS" == "macOS" ]]; then
    # Make sure to use a libomp version binary compatible with the oldest
    # supported version of the macos SDK as libomp will be vendored into the
    # scikit-learn wheels for macos. The list of binaries are in
    # https://packages.macports.org/libomp/.
    if [[ "$CIBW_BUILD" == *-macosx_arm64 ]]; then
        # arm64 builds must cross compile because CI is on x64
        export PYTHON_CROSSENV=1
        # SciPy requires 12.0 on arm to prevent kernel panics
        # https://github.com/scipy/scipy/issues/14688
        # We use the same deployment target to match SciPy.
        export MACOSX_DEPLOYMENT_TARGET=12.0
        wget https://packages.macports.org/libomp/libomp-11.0.1_0.darwin_20.arm64.tbz2 -O libomp.tbz2
    else
        # Currently, the oldest supported macos version is: High Sierra / 10.13.
        # Note that Darwin_17 == High Sierra / 10.13.
        export MACOSX_DEPLOYMENT_TARGET=10.13
        wget https://packages.macports.org/libomp/libomp-11.0.1_0+universal.darwin_17.i386-x86_64.tbz2 -O libomp.tbz2
    fi
    sudo tar -C / -xvjf libomp.tbz2 opt

    export CC=/usr/bin/clang
    export CXX=/usr/bin/clang++
    export CPPFLAGS="$CPPFLAGS -Xpreprocessor -fopenmp"
    export CFLAGS="$CFLAGS -I/opt/local/include/libomp"
    export CXXFLAGS="$CXXFLAGS -I/opt/local/include/libomp"
    export LDFLAGS="$LDFLAGS -Wl,-rpath,/opt/local/lib/libomp -L/opt/local/lib/libomp -lomp"
fi

# The version of the built dependencies are specified
# in the pyproject.toml file, while the tests are run
# against the most recent version of the dependencies

python -m pip install cibuildwheel
python -m cibuildwheel --output-dir wheelhouse
