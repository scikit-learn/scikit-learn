#!/bin/bash

pip install --upgrade pip || travis_terminate $?
pip install pytest pytest-xdist || travis_terminate $?

# Test that there are no links to system libraries in the threadpoolctl
# section of the show_versions output.
python -c "import sklearn; sklearn.show_versions()" || travis_terminate $?
python -m pytest -n $CPU_COUNT --pyargs sklearn || travis_terminate $?
