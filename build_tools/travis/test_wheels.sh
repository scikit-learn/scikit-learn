#!/bin/bash

set -e
set -x

TEST_CMD="pytest --pyargs"

if [[ $TRAVIS_CPU_ARCH == arm64 ]]; then
    TEST_CMD="$TEST_CMD -n $SKLEARN_BUILD_PARALLEL"
fi

$TEST_CMD sklearn

# Test that there are no links to system libraries
python -m threadpoolctl -i sklearn
