#!/bin/bash

set -e
set -x

# OpenMP is not present on macOS by default
if [ "$RUNNER_OS" == "macOS" ]; then
    brew install libomp
    echo "CC=/usr/bin/clang" >> $GITHUB_ENV
    echo "CXX=/usr/bin/clang++" >> $GITHUB_ENV
    echo "CPPFLAGS=$CPPFLAGS -Xpreprocessor -fopenmp" >> $GITHUB_ENV
    echo "CFLAGS=$CFLAGS -I/usr/local/opt/libomp/include" >> $GITHUB_ENV
    echo "CXXFLAGS=$CXXFLAGS -I/usr/local/opt/libomp/include" >> $GITHUB_ENV
    echo "LDFLAGS=$LDFLAGS -Wl,-rpath,/usr/local/opt/libomp/lib -L/usr/local/opt/libomp/lib -lomp" >> $GITHUB_ENV
fi

# The version of the built dependencies are specified
# in the pyproject.toml file, while the tests are run
# against the most recent version of the dependencies

python -m pip install cibuildwheel
python -m cibuildwheel --output-dir wheelhouse
