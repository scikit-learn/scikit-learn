#!/bin/bash

pip install --upgrade pip || travis_terminate $?
pip install pytest pytest-xdist || travis_terminate $?

XDIST_ARGS=""
if [[ "$CPU_COUNT" != "1" ]]; then
    XDIST_ARGS="-n $CPU_COUNT"
fi

python -m pytest $XDIST_ARGS --pyargs sklearn || travis_terminate $?

# Test that there are no links to system libraries
python -m threadpoolctl -i sklearn || travis_terminate $?
