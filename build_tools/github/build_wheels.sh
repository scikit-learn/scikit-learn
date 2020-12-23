#!/bin/bash

set -e
set -x

# OpenMP is not present on macOS by default
if [[ "$RUNNER_OS" == "macOS" ]]; then
    # Use the oldest support version of macos SDK as libomp will be vendored
    # into the scikit-learn wheels on macos.
    wget https://homebrew.bintray.com/bottles/libomp-11.0.0.high_sierra.bottle.tar.gz
    brew install libomp-11.0.0.high_sierra.bottle.tar.gz
    export CC=/usr/bin/clang
    export CXX=/usr/bin/clang++
    export CPPFLAGS="$CPPFLAGS -Xpreprocessor -fopenmp"
    export CFLAGS="$CFLAGS -I/usr/local/opt/libomp/include"
    export CXXFLAGS="$CXXFLAGS -I/usr/local/opt/libomp/include"
    export LDFLAGS="$LDFLAGS -Wl,-rpath,/usr/local/opt/libomp/lib -L/usr/local/opt/libomp/lib -lomp"
fi

# The version of the built dependencies are specified
# in the pyproject.toml file, while the tests are run
# against the most recent version of the dependencies

python -m pip install cibuildwheel
python -m cibuildwheel --output-dir wheelhouse
