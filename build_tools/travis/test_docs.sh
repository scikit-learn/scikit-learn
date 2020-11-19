#!/bin/bash

set -e
set -x

if [[ $TRAVIS_CPU_ARCH != arm64 ]]; then
    # Faster tests for the documentation
    PYTEST="pytest -n $SKLEARN_BUILD_PARALLEL" make test-doc
fi
