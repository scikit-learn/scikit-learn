#!/bin/bash

pip install --upgrade pip || travis_terminate $?
pip install pytest pytest-xdist || travis_terminate $?

python -m pytest -n $CPU_COUNT --pyargs sklearn || travis_terminate $?

# Test that there are no links to system libraries
python -m threadpoolctl -i sklearn || travis_terminate $?
