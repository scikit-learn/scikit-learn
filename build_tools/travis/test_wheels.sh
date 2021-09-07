#!/bin/bash

pip install --upgrade pip || travis_terminate $?
pip install pytest pytest-xdist || travis_terminate $?

# XXX: limit to sklearn.cluster instead of full sklearn to see if overall
# test duration is causing the problem.
python -m pytest -n $CPU_COUNT --pyargs sklearn.cluster || travis_terminate $?

# Test that there are no links to system libraries
python -m threadpoolctl -i sklearn || travis_terminate $?
