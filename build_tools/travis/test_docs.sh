#!/bin/bash

set -e
set -x

if [[ $TRAVIS_CPU_ARCH != arm64 ]]; then
    # Faster run of the documentation tests
    PYTEST="pytest -n $CPU_COUNT" make test-doc
fi
